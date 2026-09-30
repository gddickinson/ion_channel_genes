"""S0 helpers — the live lookups the catalogue verification runs on.

Stdlib only, deliberately: the verification is the gate every later task
depends on, and it has to run in whatever interpreter the session starts in.
Each function returns a plain dict/list so the driver can write TSVs without
knowing anything about HTTP.

Politeness: one short sleep per request, backed-off retries, and a `Fetcher`
that records every failure rather than swallowing it. A verification run that
quietly dropped half its requests would report a clean catalogue, which is
the one outcome worse than a dirty one.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org"
BACKOFFS = (2, 5, 15, 30)


@dataclass
class Fetcher:
    """A polite JSON client that keeps a ledger of what failed."""
    sleep_s: float = 0.12
    timeout_s: int = 40
    max_retries: int = 3
    n_requests: int = 0
    failures: list[tuple[str, str]] = field(default_factory=list)

    def json(self, url: str, allow_404: bool = True):
        for attempt in range(self.max_retries + 1):
            try:
                req = urllib.request.Request(
                    url, headers={"Accept": "application/json",
                                  "User-Agent": "ion-channel-census/0.1"})
                with urllib.request.urlopen(req, timeout=self.timeout_s) as r:
                    self.n_requests += 1
                    # InterPro answers "no matches" with 204 and an empty
                    # body, which urlopen returns as success (CCDC51, the
                    # first catalogued protein with no Pfam domain).
                    body = r.read()
                    payload = json.loads(body) if (body.strip() and r.status != 204) else None
                time.sleep(self.sleep_s)
                return payload
            except urllib.error.HTTPError as e:
                if e.code in (204, 404) and allow_404:
                    self.n_requests += 1
                    time.sleep(self.sleep_s)
                    return None
                if attempt >= self.max_retries:
                    self.failures.append((url, f"HTTP {e.code}"))
                    return None
            except Exception as e:                     # network, JSON, timeout
                if attempt >= self.max_retries:
                    self.failures.append((url, type(e).__name__))
                    return None
            time.sleep(BACKOFFS[min(attempt, len(BACKOFFS) - 1)])
        return None

    def text(self, url: str) -> str:
        for attempt in range(self.max_retries + 1):
            try:
                with urllib.request.urlopen(url, timeout=self.timeout_s) as r:
                    self.n_requests += 1
                    out = r.read().decode()
                time.sleep(self.sleep_s)
                return out
            except Exception as e:
                if attempt >= self.max_retries:
                    self.failures.append((url, type(e).__name__))
                    return ""
            time.sleep(BACKOFFS[min(attempt, len(BACKOFFS) - 1)])
        return ""


# ------------------------------------------------------------ InterPro
def pfam_entry(f: Fetcher, accession: str) -> dict | None:
    """`{accession, short, name, type}` for a Pfam entry, or None if absent."""
    d = f.json(f"{INTERPRO}/entry/pfam/{accession}/")
    if not d:
        return None
    m = d.get("metadata", {})
    n = m.get("name", {})
    if isinstance(n, str):
        n = {"name": n, "short": ""}
    return {"accession": m.get("accession", accession),
            "short": n.get("short", ""), "name": n.get("name", ""),
            "type": m.get("type", "")}


def protein_pfams(f: Fetcher, accession: str) -> dict[str, int]:
    """Observed `{pfam: copies}` for a UniProt entry."""
    d = f.json(f"{INTERPRO}/entry/pfam/protein/uniprot/{accession}/?page_size=100")
    if not d:
        return {}
    out: dict[str, int] = {}
    for e in d.get("results", []):
        acc = e["metadata"]["accession"]
        n = sum(len(p.get("entry_protein_locations", []) or [])
                for p in e.get("proteins", []))
        out[acc] = max(1, n)
    return out


def pfam_protein_count(f: Fetcher, accession: str) -> int | None:
    """How many UniProt proteins carry this Pfam — the census denominator."""
    d = f.json(f"{INTERPRO}/protein/UniProt/entry/pfam/{accession}/?page_size=1")
    return None if d is None else d.get("count")


# ------------------------------------------------------------- UniProt
def resolve_gene(f: Fetcher, gene: str, taxon_id: int,
                 reviewed_only: bool = True) -> dict | None:
    """Gene symbol + taxon → the reviewed entry, longest wins on a tie."""
    # `taxonomy_id` includes strain-level descendants; `organism_id` is the
    # exact node. Measured 2026-08-19: GLIC (`glvI`, Gloeobacter violaceus)
    # is filed under strain PCC 7421 (251221), so `organism_id:33072` returns
    # nothing and `taxonomy_id:33072` finds it. Every prokaryotic exemplar
    # was invisible to the resolver until this changed.
    parts = [f"gene_exact:{gene}", f"taxonomy_id:{taxon_id}"]
    if reviewed_only:
        parts.append("reviewed:true")
    q = urllib.parse.quote(" AND ".join(parts))
    d = f.json(f"{UNIPROT}/uniprotkb/search?query={q}"
               f"&fields=accession,id,protein_name,length,gene_primary,"
               f"organism_name,reviewed&size=5")
    rows = (d or {}).get("results", [])
    if not rows:
        return None

    def primary(r: dict) -> str:
        return ((r.get("genes") or [{}])[0].get("geneName") or {}).get("value", "")

    # `gene_exact` matches synonyms too, and the longest hit used to win:
    # TRPC7 → TRPM2, GRIK2 → GRIK5, KCNG3 → KCNG4 (S3a). An entry whose
    # *primary* name is the query now wins; a synonym-only match is kept as
    # the fallback and flagged, never silently preferred.
    own = [r for r in rows if primary(r).upper() == gene.upper()]
    pool = own or rows
    pool.sort(key=lambda r: -(r.get("sequence", {}).get("length") or 0))
    r = pool[0]
    genes = r.get("genes") or [{}]
    return {
        "accession": r["primaryAccession"],
        "entry": r.get("uniProtkbId", ""),
        "gene": (genes[0].get("geneName") or {}).get("value", ""),
        "length": (r.get("sequence") or {}).get("length"),
        "organism": (r.get("organism") or {}).get("scientificName", ""),
        "taxon_id": (r.get("organism") or {}).get("taxonId"),
        "protein": ((r.get("proteinDescription", {}) or {})
                    .get("recommendedName", {}).get("fullName", {})
                    .get("value", "")),
        "reviewed": r.get("entryType", "").startswith("UniProtKB reviewed"),
        "synonym_match": not own,
        "n_candidates": len(rows),
    }


def entry_by_accession(f: Fetcher, accession: str) -> dict | None:
    d = f.json(f"{UNIPROT}/uniprotkb/{accession}.json?fields=accession,"
               f"gene_primary,protein_name,length,organism_name,reviewed")
    if not d:
        return None
    genes = d.get("genes") or [{}]
    return {
        "accession": d.get("primaryAccession", accession),
        "gene": (genes[0].get("geneName") or {}).get("value", ""),
        "length": (d.get("sequence") or {}).get("length"),
        "organism": (d.get("organism") or {}).get("scientificName", ""),
        "taxon_id": (d.get("organism") or {}).get("taxonId"),
        "protein": ((d.get("proteinDescription", {}) or {})
                    .get("recommendedName", {}).get("fullName", {})
                    .get("value", "")),
        "reviewed": d.get("entryType", "").startswith("UniProtKB reviewed"),
    }


def sequence(f: Fetcher, accession: str) -> str:
    text = f.text(f"{UNIPROT}/uniprotkb/{accession}.fasta")
    return "".join(l.strip() for l in text.splitlines()
                   if l and not l.startswith(">"))


def taxon(f: Fetcher, taxon_id: int) -> dict | None:
    d = f.json(f"{UNIPROT}/taxonomy/{taxon_id}.json")
    if not d:
        return None
    return {"taxon_id": d.get("taxonId"),
            "scientific": d.get("scientificName", ""),
            "common": d.get("commonName", ""),
            "rank": d.get("rank", ""),
            "lineage": ";".join(x.get("scientificName", "")
                                for x in (d.get("lineage") or [])[:4])}


# ------------------------------------------------------------- helpers
def write_tsv(path, header, rows) -> None:
    from pathlib import Path
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if x is None else str(x) for x in r) + "\n")


def read_tsv(path) -> list[dict]:
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        return []
    lines = p.read_text().splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[1:] if l.strip()]


def live_progress(path, task: str, steps: list[tuple[str, bool]],
                  workers: int = 1) -> None:
    """Write `results/session_live.json` so the dashboard shows real progress."""
    from pathlib import Path
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        "task": task, "workers": workers,
        "steps": [{"label": l, "done": d} for l, d in steps]}, indent=2))
