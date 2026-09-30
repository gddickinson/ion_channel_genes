"""s5_rescue.py — tblastn over a genome for the cells miniprot left empty.

miniprot is a spliced aligner tuned for homologues it can chain; a distant
or fragmented copy can leave no alignment at all. tblastn is the second bar:
it runs, per genome, on the **informative** zero cells miniprot found no
locus for (a family present at high confidence in another species of the
same group — or, for a genome-only species, any family present in its group),
with that family's baits minus the genome's own species.

A hit is a **trace**, never a call: an HSP is kept only if it lies outside
every locus the sweep recorded (any family, any call), so a Kv voltage
sensor inside an Hv1 locus or a cNMP domain inside a CNG locus cannot
become a trace of the other family. Traces are counted and located; the
cell is `trace`, and S5b / S16 decide what they are.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta, write_fasta  # noqa: E402
from s5_lib import BAITS_FAA, s5_dir  # noqa: E402

BIN = "/opt/anaconda3/envs/piezo1/bin"      # D18: BLAST+ lives in piezo1
RESCUE_E = 1e-5
HSP_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore"
TRACE_MERGE_BP = 20_000       # HSPs of one family this close are one trace


def blast_db(acc: str, fna: Path) -> Path:
    db = s5_dir(acc, "blastdb") / acc
    if not Path(str(db) + ".njs").exists() and not Path(str(db) + ".nal").exists():
        r = subprocess.run([f"{BIN}/makeblastdb", "-dbtype", "nucl", "-in", str(fna),
                            "-out", str(db), "-max_file_sz", "4GB"],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError(f"makeblastdb: {r.stderr[:300]}")
    return db


def run_tblastn(acc: str, fna: Path, baits: list[str], tag: str,
                threads: int = 8, faa: Path = BAITS_FAA) -> Path:
    """One tblastn over the genome for the listed baits; cached by tag."""
    out = s5_dir(acc, "rescue") / f"{tag}.tsv"
    if out.exists():
        return out
    want = set(baits)
    q = s5_dir(acc, "rescue") / f"{tag}.faa"
    write_fasta(q, [(h.split()[0], s) for h, s in iter_fasta(faa)
                    if h.split()[0] in want])
    db = blast_db(acc, fna)
    tmp = out.with_suffix(".partial")
    r = subprocess.run([f"{BIN}/tblastn", "-query", str(q), "-db", str(db),
                        "-evalue", str(RESCUE_E), "-outfmt", f"6 {HSP_FIELDS}",
                        "-num_threads", str(threads), "-max_target_seqs", "50",
                        "-out", str(tmp)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"tblastn: {r.stderr[:300]}")
    tmp.replace(out)
    return out


def read_hsps(path: Path) -> list[dict]:
    keys = HSP_FIELDS.split()
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) != len(keys):
            continue
        r = dict(zip(keys, f))
        s, e = int(r["sstart"]), int(r["send"])
        r.update(family=r["qseqid"].split("|")[0], contig=r["sseqid"],
                 start=min(s, e), end=max(s, e),
                 evalue=float(r["evalue"]), bitscore=float(r["bitscore"]),
                 pident=float(r["pident"]))
        rows.append(r)
    return rows


def outside_loci(hsps: list[dict], loci: list[dict]) -> list[dict]:
    """HSPs overlapping no recorded locus (either strand)."""
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for L in loci:
        by_contig.setdefault(L["contig"], []).append((int(L["start"]), int(L["end"])))
    return [h for h in hsps
            if not any(s <= h["end"] and h["start"] <= e
                       for s, e in by_contig.get(h["contig"], []))]


def traces(hsps: list[dict]) -> dict[str, list[dict]]:
    """family → merged trace regions (best HSP of each)."""
    out: dict[str, list[dict]] = {}
    for fam in sorted({h["family"] for h in hsps}):
        hs = sorted((h for h in hsps if h["family"] == fam),
                    key=lambda h: (h["contig"], h["start"]))
        cur = None
        for h in hs:
            if cur and h["contig"] == cur["contig"] and h["start"] <= cur["end"] + TRACE_MERGE_BP:
                cur["end"] = max(cur["end"], h["end"])
                cur["n_hsps"] += 1
                if h["bitscore"] > cur["bitscore"]:
                    cur.update(bitscore=h["bitscore"], evalue=h["evalue"],
                               pident=h["pident"], bait=h["qseqid"])
            else:
                cur = {"family": fam, "contig": h["contig"], "start": h["start"],
                       "end": h["end"], "n_hsps": 1, "bitscore": h["bitscore"],
                       "evalue": h["evalue"], "pident": h["pident"],
                       "bait": h["qseqid"]}
                out.setdefault(fam, []).append(cur)
    return out
