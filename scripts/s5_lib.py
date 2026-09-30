"""s5_lib.py — paths, the miniprot run, alignment parsing and loci for S5.

The alignment layer of the genomic sweep, ported from
`../ip3r_genes/scripts/s5_sweep_lib.py` and changed where one family became
ninety-one:

* **A locus is a set of alignments that overlap on the genome**, not
  alignments within 10 kb. With one family the parent could merge
  neighbours; here the K2P, Kv and iGluR tandem arrays put *different*
  families' genes within 10 kb of each other, and merging them would hand
  one locus to whichever family scored best.
* **The locus is not called from its baits.** Each locus's best
  translations are scored against the S3a profiles and assigned by D32
  (`s5_classify.py`), so the genome and the proteome are judged by one
  instrument. The bait family is recorded beside it as a second witness.
* **Self-exclusion**: alignments from baits drawn from the genome's own
  species are dropped at parse time (the genome-scale leave-one-out).
"""

from __future__ import annotations

import csv
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s5_genome_io import build_fai, genomes_root, read_fai  # noqa: E402

OUT_DIR = ROOT / "results" / "genome_sweep"
MANIFEST = ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv"
BAITS_FAA = OUT_DIR / "baits.faa"
BAITS_TSV = OUT_DIR / "baits.tsv"
#: Census revision r4's families get their own bait panel and their own
#: miniprot run per genome (`s5r4_sweep.py`), so S5b's panel, runs and loci
#: stay byte-identical (D43).
BAITS_R4_FAA = OUT_DIR / "baits_r4.faa"
BAITS_R4_TSV = OUT_DIR / "baits_r4.tsv"
LIVE = ROOT / "results" / "session_live.json"

#: miniprot's max intron (-G). A too-small value does not lose a gene, it
#: *splits* one into apparent fragments, so this is a threshold on the call.
#: S5a runs the pilot generously and measures the widest annotated intron at
#: family-called loci (`s5_calibrate.py`); S5b takes its value from there.
MAX_INTRON_PROK = 2_000          # no spliceosomal introns; miniprot needs a value
MAX_INTRON_EUK = 200_000         # miniprot's default
#: Measured before S5b (D38): 1.5 Mb on mouse lost 9 high-confidence loci
#: (Hvcn1 among them) to chaining across genes and gained no call needing an
#: intron > 1 Mb; Asic2's 996 kb intron is already recovered intact at 1 Mb.
#: A wider intron splits a gene into partial loci — `introns.tsv` `over_G`.
MAX_INTRON_LARGE = 1_000_000     # genomes >= LARGE_GENOME_BP
LARGE_GENOME_BP = 1_000_000_000

#: miniprot memory (measured by the parent project) and the chunk size that
#: keeps one chunk under the 32 GB machine.
MINIPROT_GB_RSS_PER_GBP = 6.32
CHUNK_BP = 2_500_000_000
THREADS = 8

#: Locus floor. Recorded, not tuned: every cluster is written with its best
#: identity, and the floor applies downstream so it can be measured (the
#: parent's S23b lesson — a floor measured on already-filtered loci is
#: circular).
RECORD_MIN_IDENTITY = 0.15

EDGE_BP = 10_000                 # locus within this of a contig end -> edge
PROK_GROUPS = ("prokaryote", "virus")


def s5_dir(*parts: str) -> Path:
    p = genomes_root() / "s5"
    for x in parts:
        p = p / x
    p.mkdir(parents=True, exist_ok=True)
    return p


def manifest() -> dict[str, dict]:
    with open(MANIFEST) as fh:
        return {r["species"]: r for r in csv.DictReader(fh, delimiter="\t")}


def sweep_assembly(row: dict) -> str:
    """The assembly S5 searches: the current version of S4's assembly.

    The proteome is compared by sequence, not by coordinates, so the newest
    version of the same assembly is the better search target.
    """
    return row.get("current_assembly") or row["assembly_id"]


def max_intron_for(group: str, genome_bp: int) -> int:
    if group in PROK_GROUPS:
        return MAX_INTRON_PROK
    return MAX_INTRON_LARGE if genome_bp >= LARGE_GENOME_BP else MAX_INTRON_EUK


def load_bait_meta() -> dict[str, dict]:
    out = {}
    for p in (BAITS_TSV, BAITS_R4_TSV):
        if p.exists():
            with open(p) as fh:
                out.update({r["bait"]: r for r in csv.DictReader(fh, delimiter="\t")})
    return out


def r4_families() -> frozenset[str]:
    """The families census revision r4 added (their loci live in `r4/`)."""
    from s3r4_sweep import NEW_PROFILES
    return frozenset(NEW_PROFILES)


# ---------------------------------------------------------------- miniprot

def run_miniprot(fna: Path, baits: Path, out_gff: Path, max_intron: int,
                 threads: int = THREADS) -> None:
    """One reference, the whole panel; the GFF is written **atomically**.

    A run killed mid-write leaves a GFF cut at a line boundary that parses
    perfectly and is short of loci; the `.partial` rename is the only way
    the next run can tell.
    """
    tmp = out_gff.with_suffix(out_gff.suffix + ".partial")
    cmd = ["miniprot", "-t", str(threads), "--gff", "--trans", "--outs=0.3",
           "-N", "60", "-G", str(max_intron), str(fna), str(baits)]
    try:
        with open(tmp, "w") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                  text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"miniprot failed: {proc.stderr.strip()[:300]}")
        tmp.replace(out_gff)
    finally:
        tmp.unlink(missing_ok=True)


def split_by_contig(fna: Path, out_dir: Path, chunk_bp: int = CHUNK_BP) -> list[Path]:
    """Whole-contig chunks, planned from the index first so none overshoots."""
    out_dir.mkdir(parents=True, exist_ok=True)
    existing = sorted(out_dir.glob("chunk_*.fna"))
    if existing and (out_dir / ".planned").exists():
        return existing
    for stale in out_dir.glob("chunk_*"):
        stale.unlink()
    idx = read_fai(build_fai(fna))
    plan, cur, cur_bp = [], [], 0
    for name, (length, *_r) in idx.items():
        if cur and cur_bp + length > chunk_bp:
            plan.append(cur)
            cur, cur_bp = [], 0
        cur.append(name)
        cur_bp += length
    if cur:
        plan.append(cur)
    where = {n: i for i, names in enumerate(plan) for n in names}
    handles = [open(out_dir / f"chunk_{i + 1:03d}.fna", "w") for i in range(len(plan))]
    try:
        fh = None
        with open(fna) as src:
            for line in src:
                if line.startswith(">"):
                    fh = handles[where[line[1:].split()[0]]]
                fh.write(line)
    finally:
        for h in handles:
            h.close()
    (out_dir / ".planned").write_text(f"{len(plan)} chunks of <= {chunk_bp} bp\n")
    return sorted(out_dir.glob("chunk_*.fna"))


def run_miniprot_any(fna: Path, baits: Path, out_gff: Path, max_intron: int,
                     genome_bp: int, work: Path) -> int:
    """Whole genome if it fits in memory, else contig chunks. Returns chunks."""
    if genome_bp <= CHUNK_BP:
        run_miniprot(fna, baits, out_gff, max_intron)
        return 1
    parts = []
    for i, chunk in enumerate(split_by_contig(fna, work), 1):
        part = work / f"chunk_{i:03d}.gff"
        if not part.exists() or part.stat().st_size == 0:
            run_miniprot(chunk, baits, part, max_intron)
        parts.append(part)
    tmp = out_gff.with_suffix(out_gff.suffix + ".partial")
    with open(tmp, "w") as out:
        out.write("##gff-version 3\n")
        for part in parts:
            with open(part) as fh:
                out.writelines(l for l in fh if not l.startswith("##gff-version"))
    tmp.replace(out_gff)
    for c in work.glob("chunk_*.fna"):
        c.unlink()
    return len(parts)


# ---------------------------------------------------------------- alignments

@dataclass
class Aln:
    contig: str
    start: int
    end: int
    strand: str
    score: float
    identity: float
    bait: str
    family: str
    q_start: int
    q_end: int
    bait_len: int
    frameshifts: int = 0
    stop_codons: int = 0
    mp_id: str = ""
    translation: str = ""
    aligned_aa: int = 0
    cds_blocks: list = field(default_factory=list)

    @property
    def coverage(self) -> float:
        if not self.bait_len:
            return 0.0
        n = self.aligned_aa or (self.q_end - self.q_start + 1)
        return min(1.0, n / self.bait_len)

    @property
    def max_intron(self) -> int:
        b = sorted(self.cds_blocks)
        return max((b[i + 1][0] - b[i][1] - 1 for i in range(len(b) - 1)),
                   default=0)


def _attr(attrs: str, key: str) -> str:
    for part in attrs.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


def _union_len(spans: list[tuple[int, int]]) -> int:
    if not spans:
        return 0
    spans = sorted(spans)
    total, (cs, ce) = 0, spans[0]
    for s, e in spans[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    return total + ce - cs + 1


def parse_miniprot_gff(path: Path, bait_meta: dict[str, dict],
                       exclude_species: str = "") -> tuple[list[Aln], int]:
    """mRNA + CDS + ##STA → alignments; returns (alignments, n_self_dropped).

    ##STA precedes its mRNA. Chunked GFFs repeat miniprot's IDs, so CDS lines
    bind to the mRNA just seen. Alignments from baits of `exclude_species`
    are dropped (the genome-scale leave-one-out).
    """
    alns, by_id, q_spans = [], {}, {}
    pending, cur_raw, cur_key, dropped, skip = "", "", "", 0, False
    for line in open(path):
        if line.startswith("##STA"):
            pending = line.rstrip("\n").split("\t", 1)[1] if "\t" in line else ""
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        parts = _attr(f[8], "Target").split()
        if f[2] == "mRNA":
            bait = parts[0] if parts else ""
            meta = bait_meta.get(bait, {})
            if exclude_species and meta.get("species") == exclude_species:
                dropped += 1
                skip, pending = True, ""
                continue
            skip = False
            a = Aln(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                    score=float(f[5]) if f[5] not in (".", "") else 0.0,
                    identity=float(_attr(f[8], "Identity") or 0.0),
                    bait=bait, family=meta.get("family", "unknown"),
                    q_start=int(parts[1]) if len(parts) > 2 else 0,
                    q_end=int(parts[2]) if len(parts) > 2 else 0,
                    bait_len=int(meta.get("length") or 0),
                    frameshifts=int(_attr(f[8], "Frameshift") or 0),
                    stop_codons=int(_attr(f[8], "StopCodon") or 0),
                    mp_id=_attr(f[8], "ID"), translation=pending)
            pending = ""
            cur_raw = cur_key = a.mp_id
            if cur_key in by_id:
                cur_key = f"{cur_raw}#{len(alns)}"
            a.mp_id = cur_key
            alns.append(a)
            by_id[cur_key] = a
        elif f[2] == "CDS" and not skip:
            parent = _attr(f[8], "Parent")
            key = cur_key if parent == cur_raw else parent
            if len(parts) > 2:
                q_spans.setdefault(key, []).append((int(parts[1]), int(parts[2])))
            if key in by_id:
                by_id[key].cds_blocks.append((int(f[3]), int(f[4])))
    for k, spans in q_spans.items():
        if k in by_id:
            by_id[k].aligned_aa = _union_len(spans)
    return alns, dropped


# ---------------------------------------------------------------- loci

@dataclass
class Locus:
    contig: str
    strand: str
    start: int
    end: int
    alns: list[Aln] = field(default_factory=list)

    @property
    def best(self) -> Aln:
        return max(self.alns, key=lambda a: a.score)

    def best_per_family(self, n: int = 2) -> list[Aln]:
        """The top `n` alignments of every bait family present here — the
        translations the profiles are asked to judge."""
        by_f: dict[str, list[Aln]] = {}
        for a in self.alns:
            by_f.setdefault(a.family, []).append(a)
        out = []
        for lst in by_f.values():
            out += sorted(lst, key=lambda a: -a.score)[:n]
        return out

    def bait_family_margin(self) -> tuple[str, str, float]:
        """(winning bait family, runner-up, D7 relative margin) on miniprot
        score — the second witness beside the profile call."""
        best: dict[str, float] = {}
        for a in self.alns:
            best[a.family] = max(best.get(a.family, 0.0), a.score)
        ranked = sorted(best.items(), key=lambda kv: -kv[1])
        if len(ranked) == 1:
            return ranked[0][0], "", 1.0
        (w, ws), (r, rs) = ranked[0], ranked[1]
        return w, r, round((ws - rs) / ws, 4) if ws > 0 else 0.0


def cluster_loci(alns: list[Aln]) -> list[Locus]:
    """Alignments that overlap on the same contig and strand are one locus."""
    by_key: dict[tuple, list[Aln]] = {}
    for a in alns:
        if a.identity >= RECORD_MIN_IDENTITY:
            by_key.setdefault((a.contig, a.strand), []).append(a)
    loci = []
    for (contig, strand), group in by_key.items():
        group.sort(key=lambda a: a.start)
        cur = None
        for a in group:
            if cur and a.start <= cur.end:
                cur.end = max(cur.end, a.end)
                cur.alns.append(a)
            else:
                cur = Locus(contig, strand, a.start, a.end, [a])
                loci.append(cur)
    return loci
