"""S10b shared pieces: archived fetches, the S4 panel as a BLAST reference,
and the region/lineage logic of D51 (1).

Bulk under `<data root>/raw_api/s10b/` (fetches) and
`<data root>/blast_db/s4_panel/` (the panel BLAST DB).
"""
from __future__ import annotations

import json
import subprocess
import time
import urllib.request
from pathlib import Path

from s3_hmm_lib import read_tsv

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv"
ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")
METAZOA_GROUPS = {"basal_metazoan", "cnidarian", "invertebrate", "deuterostome",
                  "vertebrate"}
#: D51 (1): a region's lineage is its best hit's only past D7's margin
LINEAGE_MARGIN = 0.10
EVALUE = 1e-10
LOCUS_PAD = 500
WINDOW = 100_000          # D51 (1): contigs > 200 kb are searched at locus ± 100 kb
WINDOW_ABOVE = 200_000


def data_root() -> Path:
    from src.utils.data_root import require_data_root
    return Path(require_data_root())


def raw_dir(*parts: str) -> Path:
    d = data_root() / "raw_api" / "s10b"
    for p in parts:
        d = d / p
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch(url: str, dest: Path, attempts: int = 4) -> str:
    """GET once, archive, and read from the archive thereafter."""
    if dest.exists() and dest.stat().st_size:
        return dest.read_text()
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ion-channel-census/S10b"})
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode()
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(text)
            return text
        except Exception as exc:       # retried, then raised — never swallowed
            last = exc
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"fetch failed: {url}: {last}")


def fetch_json(url: str, dest: Path) -> dict:
    return json.loads(fetch(url, dest))


# ------------------------------------------------------------- lineages
def lineage_of_group(group: str) -> str:
    if group == "prokaryote":
        return "prokaryote"
    if group == "virus":
        return "virus"
    return "Metazoa" if group in METAZOA_GROUPS else "other_eukaryote"


def panel_taxa() -> dict[str, dict]:
    """proteome OX taxid → {species, lineage}."""
    return {r["proteome_taxid"]: {"species": r["species"],
                                  "lineage": lineage_of_group(r["group"])}
            for r in read_tsv(MANIFEST) if r["status"] == "proteome"}


def panel_fasta() -> Path:
    return data_root() / "proteomes" / "s4" / "panel_refprot.fasta"


def panel_db() -> Path:
    """makeblastdb of S4's panel FASTA (built once) + id → OX map."""
    d = data_root() / "blast_db" / "s4_panel"
    d.mkdir(parents=True, exist_ok=True)
    db = d / "panel_refprot"
    if not (d / "panel_refprot.pdb").exists() and not (d / "panel_refprot.00.pin").exists() \
            and not (d / "panel_refprot.pin").exists():
        subprocess.run([str(ENV_BIN / "makeblastdb"), "-in", str(panel_fasta()),
                        "-dbtype", "prot", "-out", str(db)], check=True,
                       stdout=subprocess.DEVNULL)
    return db


def panel_ox() -> dict[str, str]:
    """first header token → OX taxid, cached beside the DB."""
    cache = data_root() / "blast_db" / "s4_panel" / "ox.json"
    if cache.exists():
        return json.loads(cache.read_text())
    out = {}
    with open(panel_fasta()) as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].split()[0]
                ox = line.split(" OX=")[1].split()[0] if " OX=" in line else ""
                out[tok] = ox
    cache.write_text(json.dumps(out))
    return out


# ------------------------------------------------------------- BLAST
def run_blast(prog: str, query: Path, out: Path, *, db: Path | None = None,
              subject: Path | None = None, threads: int = 9,
              extra: list[str] | None = None) -> list[dict]:
    fmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
    if not (out.exists() and out.stat().st_size):
        cmd = [str(ENV_BIN / prog), "-query", str(query), "-outfmt", fmt,
               "-evalue", "1e-5", "-out", str(out) + ".part"]
        cmd += ["-db", str(db), "-num_threads", str(threads)] if db else \
               ["-subject", str(subject)]
        cmd += extra or []
        subprocess.run(cmd, check=True)
        Path(str(out) + ".part").rename(out)
    rows = []
    for line in out.read_text().splitlines():
        f = line.split("\t")
        rows.append({"q": f[0], "s": f[1], "pident": float(f[2]), "len": int(f[3]),
                     "qs": int(f[4]), "qe": int(f[5]), "ss": int(f[6]), "se": int(f[7]),
                     "evalue": float(f[8]), "bits": float(f[9])})
    return rows


def merge_regions(hsps: list[dict]) -> list[list[dict]]:
    """HSPs merged by overlap of their query (contig) interval."""
    iv = sorted(hsps, key=lambda h: min(h["qs"], h["qe"]))
    regions, end = [], -1
    for h in iv:
        lo, hi = min(h["qs"], h["qe"]), max(h["qs"], h["qe"])
        if regions and lo <= end:
            regions[-1].append(h)
            end = max(end, hi)
        else:
            regions.append([h])
            end = hi
    return regions


def region_lineage(hsps: list[dict], ox: dict, taxa: dict) -> dict:
    """Best hit per lineage; the region's lineage past LINEAGE_MARGIN."""
    best: dict[str, dict] = {}
    for h in hsps:
        t = taxa.get(ox.get(h["s"], ""), {})
        lin = t.get("lineage", "unknown")
        if lin not in best or h["bits"] > best[lin]["bits"]:
            best[lin] = {**h, "species": t.get("species", "")}
    ranked = sorted(best.items(), key=lambda kv: -kv[1]["bits"])
    top_lin, top = ranked[0]
    second = ranked[1][1]["bits"] if len(ranked) > 1 else 0.0
    margin = (top["bits"] - second) / top["bits"]
    lo = min(min(h["qs"], h["qe"]) for h in hsps)
    hi = max(max(h["qs"], h["qe"]) for h in hsps)
    return {"start": lo, "end": hi, "n_hsps": len(hsps),
            "best_lineage": top_lin, "best_hit": top["s"], "best_species": top["species"],
            "best_bits": top["bits"], "best_pident": top["pident"],
            "runner_lineage": ranked[1][0] if len(ranked) > 1 else "",
            "runner_bits": second, "margin": round(margin, 3),
            "lineage": top_lin if margin >= LINEAGE_MARGIN else "unclear"}
