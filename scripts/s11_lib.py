"""s11_lib.py — S11a shared pieces (D52): species, species tree, human gene
symbols, OHNOLOGS v2. Every fetch is archived under `<data root>/raw_api/s11/`
and an offline rerun reads the archive.
"""

from __future__ import annotations

import csv
import gzip
import html
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s10_lib import build_tree, fetch_lineages, raw          # noqa: E402
from s11_recon import stree_from                               # noqa: E402

MANIFEST = ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv"
TIER1 = ROOT / "results" / "phylogeny" / "tier1_trees.tsv"
TREES = ROOT / "results" / "phylogeny" / "tier1"
MEMBERS = ROOT / "results" / "alignments" / "members.tsv"
OUT = ROOT / "results" / "duplication"
LIVE = ROOT / "results" / "session_live.json"
S7D_UNDEFINED = {"cng", "k2p", "kv_kcnq", "kv_shaker", "kca_slo", "kv_eag"}
OG_PREFIX = "OG_"
HGNC_URL = ("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
            "hgnc_complete_set.txt")
OHNO_URL = "http://ohnologs.curie.fr/cgi-bin/Browse.cgi"
OHNO_CRITERIA = {"strict": "A", "relaxed": "C"}       # D52 (6)


def s11raw(*parts: str) -> Path:
    return raw(*parts, task="s11")


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------- species
def code(species: str) -> str:
    """The species suffix `s6_lib.label_for` puts on every tip."""
    return "".join(w[:3] for w in species.split()[:2])


def panel() -> list[dict]:
    rows = read_tsv(MANIFEST)
    codes = [code(r["species"]) for r in rows]
    if len(set(codes)) != len(codes):
        raise RuntimeError("two panel species share a tip code")
    return rows


def species_of_tip(tip: str, by_code: dict[str, str]) -> str:
    return by_code[tip.rsplit("__", 1)[1]]


def species_tree():
    """D52 (2): NCBI Taxonomy for the 52 panel taxids, polytomies kept."""
    rows = panel()
    lin = fetch_lineages([r["panel_taxid"] for r in rows], name="s4_panel", task="s11")
    tips = {}
    for r in rows:
        tips[r["species"]] = [("1", "root", "no rank")] + lin[r["panel_taxid"]]
    t = build_tree(tips, root="root")
    return t, stree_from(t)


# ---------------------------------------------------------------- archive
def fetch_archived(url: str, path: Path, params: dict | None = None,
                   refresh: bool = False, timeout: int = 300) -> str:
    if path.exists() and not refresh:
        return path.read_text(encoding="latin-1")
    for attempt in range(5):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code == 200 and r.content:
                break
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    else:
        raise RuntimeError(f"fetch failed: {url} {params}")
    path.write_bytes(r.content)
    (path.parent / (path.name + ".date")).write_text(time.strftime("%Y-%m-%d") + "\n")
    return r.content.decode("latin-1")


# ---------------------------------------------------------------- human genes
def hgnc_table(refresh: bool = False) -> list[dict]:
    text = fetch_archived(HGNC_URL, s11raw("hgnc", "hgnc_complete_set.txt"), refresh=refresh)
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def hgnc_maps(refresh: bool = False) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    """→ (UniProt acc → HGNC id, Ensembl gene → HGNC id, HGNC id → symbol)."""
    acc, ens, sym = {}, {}, {}
    for r in hgnc_table(refresh):
        hid = r["hgnc_id"]
        sym[hid] = r["symbol"]
        for a in filter(None, r.get("uniprot_ids", "").split("|")):
            acc.setdefault(a, hid)
        if r.get("ensembl_gene_id"):
            ens[r["ensembl_gene_id"]] = hid
    return acc, ens, sym


def human_gn() -> dict[str, str]:
    """Human proteome accession → `GN=` symbol (S4's file, offline)."""
    from src.utils.data_root import require_data_root
    d = require_data_root() / "proteomes" / "s4"
    path = next(d.glob("UP000005640_9606.fasta.gz"))
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line[1:].split("|")[1]
                m = re.search(r" GN=(\S+)", line)
                out[acc] = m.group(1) if m else ""
    return out


# ---------------------------------------------------------------- OHNOLOGS v2
def ohnolog_pairs(criterion: str, refresh: bool = False) -> list[dict]:
    """Human 2R ohnologue pairs, OHNOLOGS v2 (Singh & Isambert 2020)."""
    c = OHNO_CRITERIA[criterion]
    path = s11raw("ohnologs", f"hsapiens_2R_pairs_{criterion}.html")
    text = fetch_archived(f"{OHNO_URL}?crit=[{c}]&org=hsapiens&opt=pairs&wgd=2R", path,
                          refresh=refresh)
    out = []
    for row in re.findall(r"<tr>(.*?)</tr>", text, re.S):
        cells = [html.unescape(re.sub(r"<[^>]+>", "", x)).strip()
                 for x in re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)]
        if len(cells) < 8:
            continue
        e1, e2 = (re.search(r"ENSG\d+", x).group(0) for x in cells[:2])
        out.append({"ens1": e1, "ens2": e2, "sym1": cells[2], "sym2": cells[3],
                    "time": cells[7]})
    if not out:
        raise RuntimeError(f"no OHNOLOGS pairs parsed from {path}")
    return out
