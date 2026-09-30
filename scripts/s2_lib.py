"""S2 shared pieces — paths, the taxonomic shards, polite fetching, record parsing.

The census space is every UniProtKB protein carrying at least one of the
catalogue's pore signatures (`registry.pore_signatures()`, never `PF00520`
alone — H7). Measured on UniProt 2026_03 that is ~1.25 million proteins,
~2,500 pages of 500. One cursor walk would take ~3.5 h, so the query is cut
into **taxonomic shards that partition it exactly** — every protein falls in
one and only one — and the shards walk in parallel. The partition is checked,
not assumed: `s2_enumerate.py` refuses to proceed unless the shard counts sum
to the union count.

Why UniProt and not InterPro's per-signature listing (D31): the architecture
tier needs each protein's *complete* Pfam list with copy numbers, including
the 19 accessions the rules consult that are not pore signatures (forbid
lists, subfamily markers). InterPro's per-signature listing returns only the
queried entry. UniProt's JSON carries every Pfam cross-reference with its
`MatchStatus` copy count, the transmembrane features and the sequence, in
one record. InterPro's count for each signature is still fetched, as the
independent check.
"""

from __future__ import annotations

import gzip
import json
import random
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.catalogue.registry import pore_signatures          # noqa: E402
from src.utils.data_root import require_data_root           # noqa: E402

OUT_DIR = ROOT / "results" / "census_v2"
LIVE = ROOT / "results" / "session_live.json"

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"
FIELDS = ("accession,reviewed,organism_id,organism_name,length,fragment,"
          "xref_pfam,ft_transmem,gene_primary,sequence")
PAGE = 500

#: Taxonomic shards. Each is (key, clause); the clauses partition UniProtKB.
_B, _A, _E, _V = "taxonomy_id:2", "taxonomy_id:2157", "taxonomy_id:2759", "taxonomy_id:10239"
SHARDS: list[tuple[str, str]] = [
    ("bact_pseudomonadota", f"({_B}) AND (taxonomy_id:1224)"),
    ("bact_other",          f"({_B}) NOT (taxonomy_id:1224)"),
    ("archaea",             f"({_A})"),
    ("virus",               f"({_V})"),
    ("mammalia",            "(taxonomy_id:40674)"),
    ("actinopterygii",      "(taxonomy_id:7898)"),
    ("chordata_other",      "(taxonomy_id:7711) NOT (taxonomy_id:40674) NOT "
                            "(taxonomy_id:7898)"),
    ("metazoa_other",       "(taxonomy_id:33208) NOT (taxonomy_id:7711)"),
    ("viridiplantae",       "(taxonomy_id:33090)"),
    ("fungi",               "(taxonomy_id:4751)"),
    ("euk_other",           f"({_E}) NOT (taxonomy_id:33208) NOT "
                            "(taxonomy_id:33090) NOT (taxonomy_id:4751)"),
    ("unplaced",            f"NOT ({_B}) NOT ({_A}) NOT ({_E}) NOT ({_V})"),
]


#: Delta revisions of census v2 (D43): r4 the families added after S20, r5 the
#: viroporin signatures (S4's row). Each walked by `s2r4_delta.py --rev`, per
#: taxonomic shard: records carrying a new signature and none enumerated before.
REVISIONS = ("r4", "r5")
R4_SHARDS: list[tuple[str, str]] = [(f"r4_{k}", c) for k, c in SHARDS]


def _walked(key: str) -> bool:
    try:
        from src.utils.data_root import get_data_root
        return (get_data_root() / "raw_api" / "s2" / "pages" / key / "state.json").exists()
    except Exception:                 # no drive: only what is known to exist
        return False


#: Every shard the census is assembled from: the 12, and each revision's
#: shards once walked (r4's always — census v2 r4 is built).
CENSUS_SHARDS: list[tuple[str, str]] = SHARDS + [
    (f"{rev}_{k}", c) for rev in REVISIONS for k, c in SHARDS
    if rev == "r4" or _walked(f"{rev}_{k}")]


def signatures() -> list[str]:
    return sorted(s.accession for s in pore_signatures())


def union_clause() -> str:
    return " OR ".join(f"(xref:pfam-{a})" for a in signatures())


def shard_query(clause: str) -> str:
    return f"({union_clause()}) AND ({clause})" if not clause.startswith("NOT") \
        else f"({union_clause()}) {clause}"


def raw_dir() -> Path:
    """Bulk output lives on the external drive (D1) — never in the repo."""
    d = require_data_root() / "raw_api" / "s2"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------- fetching
def get(url: str, params: dict | None = None, timeout_s: int = 120,
        tries: int = 10) -> requests.Response:
    """GET with jittered exponential backoff; raises after `tries` (strict)."""
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout_s)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r
        except (requests.RequestException, ValueError) as exc:
            last = exc
            time.sleep(min(120, 2 ** i) * (0.5 + random.random()))
    raise RuntimeError(f"GET failed after {tries} tries: {url} ({last})")


def uniprot_count(query: str) -> tuple[int, str]:
    r = get(UNIPROT, {"query": query, "size": 0, "format": "json"})
    return int(r.headers["x-total-results"]), r.headers.get("x-uniprot-release", "")


def interpro_count(pfam: str) -> int | None:
    r = get(f"{INTERPRO}/protein/uniprot/entry/pfam/{pfam}/", {"page_size": 1})
    return r.json().get("count")


def next_link(r: requests.Response) -> str | None:
    link = r.headers.get("link", "")
    if 'rel="next"' not in link:
        return None
    return link.split(";")[0].strip().strip("<>")


# ------------------------------------------------------------- parsing
def parse_record(x: dict) -> dict:
    """One UniProt JSON result → one flat census row (sequence included)."""
    org = x.get("organism", {})
    lin = org.get("lineage", [])
    pf = {}
    for ref in x.get("uniProtKBCrossReferences", []):
        if ref.get("database") != "Pfam":
            continue
        props = {p["key"]: p["value"] for p in ref.get("properties", [])}
        try:
            pf[ref["id"]] = int(props.get("MatchStatus", "1"))
        except ValueError:
            pf[ref["id"]] = 1
    tm = sum(1 for f in x.get("features", []) if f.get("type") == "Transmembrane")
    genes = x.get("genes") or [{}]
    seq = x.get("sequence", "")
    if isinstance(seq, dict):
        seq = seq.get("value", "")
    return {
        "accession": x.get("primaryAccession", ""),
        "reviewed": "reviewed" if "Swiss-Prot" in x.get("entryType", "") else "unreviewed",
        "taxon_id": org.get("taxonId", ""),
        "organism": org.get("scientificName", ""),
        "domain": lin[0] if lin else "",
        "group": lin[1] if len(lin) > 1 else "",
        "phylum": lin[2] if len(lin) > 2 else "",
        "length": len(seq),
        "fragment": (x.get("proteinDescription") or {}).get("flag", ""),
        "gene": (genes[0].get("geneName") or {}).get("value", ""),
        "pfam": ";".join(f"{k}:{v}" for k, v in sorted(pf.items())),
        "tm_count": tm,
        "sequence": seq,
    }


def pfam_dict(s: str) -> dict[str, int]:
    out = {}
    for tok in filter(None, s.split(";")):
        k, _, v = tok.partition(":")
        out[k] = int(v or 1)
    return out


def iter_pages(shard: str):
    """Yield parsed records from one shard's archived pages, in order."""
    for p in sorted((raw_dir() / "pages" / shard).glob("page_*.json.gz")):
        with gzip.open(p, "rt") as fh:
            for x in json.load(fh).get("results", []):
                yield parse_record(x)


def live(steps: list[tuple[str, bool]], workers: int = 1) -> None:
    LIVE.parent.mkdir(parents=True, exist_ok=True)
    LIVE.write_text(json.dumps({"task": "S2", "workers": workers, "steps": [
        {"label": l, "done": d} for l, d in steps]}, indent=2))
