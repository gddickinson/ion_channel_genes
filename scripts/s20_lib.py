"""s20_lib.py — S20 shared pieces: the three database channelomes and the
catalogue, joined on HGNC ids.

**What a "published total" is here.** Three curated lists of human ion
channels that are maintained by databases and can be re-fetched, each
archived byte-for-byte under `<data root>/raw_api/s20/` with its release:

* **GtoPdb** — IUPHAR/BPS Guide to Pharmacology `targets_and_families.csv`,
  rows of type `vgic`, `lgic` and `other_ic`.
* **HGNC** — the gene group "Ion channels" (id 177) and every group under it.
* **UniProt** — reviewed human entries carrying keyword KW-0407 "Ion channel".

Every list and every catalogue gene is reduced to an **HGNC id**, so a
renamed symbol (AQP0 → MIP) cannot make a gene look missing. A catalogue
symbol HGNC cannot resolve is reported, never dropped.
"""

from __future__ import annotations

import csv
import io
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import Fetcher  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "auxiliary"
GTOP_URL = "https://www.guidetopharmacology.org/DATA/targets_and_families.csv"
HGNC_GROUP_URL = "https://www.genenames.org/cgi-bin/genegroup/download?id=177&type=branch"
UNIPROT_URL = ("https://rest.uniprot.org/uniprotkb/stream?query=keyword:KW-0407"
               "+AND+organism_id:9606+AND+reviewed:true&format=tsv"
               "&fields=accession,gene_primary,xref_hgnc,protein_name")
HGNC_REST = "https://rest.genenames.org/fetch"
GTOP_TYPES = ("vgic", "lgic", "other_ic")
LISTS = ("gtopdb", "hgnc", "uniprot")


def raw_dir(*parts: str) -> Path:
    p = require_data_root() / "raw_api" / "s20"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def fetch_archived(url: str, name: str, refresh: bool = False) -> tuple[str, dict]:
    """Body + response headers, archived; offline rerun reads the archive."""
    path, meta = raw_dir() / name, raw_dir() / f"{name}.headers.json"
    if path.exists() and meta.exists() and not refresh:
        return path.read_text(), json.loads(meta.read_text())
    req = urllib.request.Request(url, headers={"User-Agent": "ion-channel-census/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read().decode()
        headers = {k: v for k, v in r.headers.items()
                   if k.lower() in ("x-uniprot-release", "x-uniprot-release-date",
                                    "last-modified", "date", "x-total-results")}
    path.write_text(body)
    meta.write_text(json.dumps({"url": url, **headers}, indent=1))
    return body, headers


def gtopdb_list(refresh: bool = False) -> tuple[dict[str, dict], str]:
    body, _ = fetch_archived(GTOP_URL, "gtopdb_targets_and_families.csv", refresh)
    lines = body.splitlines()
    version = lines[0].strip('"# ')
    out: dict[str, dict] = {}
    for r in csv.DictReader(io.StringIO("\n".join(lines[1:]))):
        if r["Type"] in GTOP_TYPES and r["HGNC id"]:
            out.setdefault(f"HGNC:{r['HGNC id']}", {
                "symbol": r["HGNC symbol"], "group": _strip_html(r["Family name"]),
                "type": r["Type"]})
    return out, version


def hgnc_list(refresh: bool = False) -> tuple[dict[str, dict], str]:
    body, h = fetch_archived(HGNC_GROUP_URL, "hgnc_group177.tsv", refresh)
    out: dict[str, dict] = {}
    for r in csv.DictReader(io.StringIO(body), delimiter="\t"):
        if r["Status"] == "Approved":
            out.setdefault(r["HGNC ID"], {"symbol": r["Approved symbol"],
                                          "group": r["Group name"],
                                          "type": r["Locus type"]})
    return out, f"HGNC group 177 branch, fetched {h.get('date', h.get('Date', ''))}"


def uniprot_list(refresh: bool = False) -> tuple[dict[str, dict], str, list[str]]:
    body, h = fetch_archived(UNIPROT_URL, "uniprot_kw0407_human.tsv", refresh)
    out: dict[str, dict] = {}
    no_hgnc = []
    for r in csv.DictReader(io.StringIO(body), delimiter="\t"):
        ids = [x for x in r["HGNC"].split(";") if x]
        if not ids:
            no_hgnc.append(f"{r['Entry']}:{r['Gene Names (primary)']}")
            continue
        out.setdefault(ids[0], {"symbol": r["Gene Names (primary)"],
                                "group": r["Protein names"][:60], "type": r["Entry"]})
    rel = h.get("X-UniProt-Release", h.get("x-uniprot-release", ""))
    return out, f"UniProtKB {rel}", no_hgnc


def _strip_html(s: str) -> str:
    import re
    return re.sub(r"<[^>]+>", "", s).replace("  ", " ").strip()


def hgnc_resolve(symbols: list[str], refresh: bool = False) -> dict[str, dict]:
    """symbol → {hgnc_id, approved, how}: approved, else previous, else alias.

    One archived JSON per symbol; an ambiguous previous/alias symbol (more
    than one HGNC record) is left unresolved rather than guessed.
    """
    f = Fetcher(sleep_s=0.12)
    d = raw_dir("hgnc")
    out = {}
    for sym in symbols:
        path = d / f"{sym}.json"
        if path.exists() and not refresh:
            rec = json.loads(path.read_text())
        else:
            rec = {}
            for field in ("symbol", "prev_symbol", "alias_symbol"):
                j = f.json(f"{HGNC_REST}/{field}/{sym}")
                docs = (j or {}).get("response", {}).get("docs", [])
                if docs:
                    rec = {"field": field, "docs": [{"hgnc_id": x["hgnc_id"],
                                                     "symbol": x["symbol"]} for x in docs]}
                    break
            if f.failures:
                raise RuntimeError(f"HGNC lookup failed: {f.failures[-1]}")
            path.write_text(json.dumps(rec))
        docs = rec.get("docs", [])
        if len(docs) == 1:
            out[sym] = {"hgnc_id": docs[0]["hgnc_id"], "approved": docs[0]["symbol"],
                        "how": rec["field"]}
        else:
            out[sym] = {"hgnc_id": "", "approved": "",
                        "how": "ambiguous" if docs else "unresolved"}
    return out


def catalogue_genes() -> list[dict]:
    """Every human gene the catalogue names, with its family and status."""
    rows = []
    for fam in CATALOGUE.values():
        for g in fam.human_genes:
            rows.append({"symbol": g, "family": fam.key,
                         "superfamily": fam.superfamily, "status": fam.status.name,
                         "census": fam.census_member()})
    return rows


def category(status: str, census: bool) -> str:
    """The S20 reading of a catalogue status."""
    if census:
        return "pore_census"
    return {"CHANNEL_ASSOCIATED": "auxiliary", "TRANSPORTER": "transporter",
            "NON_CHANNEL_HOMOLOG": "non_channel_homolog",
            "OUT_OF_SCOPE": "out_of_scope"}.get(status, status.lower())
