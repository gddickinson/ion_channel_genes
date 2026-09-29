"""S4 shared pieces: the proteome scope's rules, fetchers and file checks.

The denominator of every absence claim this project makes is the set of
proteomes it declares here, one row per species in `src/utils/species.py`.
Three rules decide that set, and each is code rather than prose:

* **Candidates** are UniProt *reference* proteomes under the panel taxon
  (strains included — `taxonomy_id:` is a subtree query) **that the pinned
  release actually ships** (`README` of the reference-proteome FTP). The live
  API runs ahead of the release, and a proteome the release lacks would give
  metadata with no file behind it.
* **Selection** (`select_proteome`): most Swiss-Prot entries, then BUSCO
  score, then gene count, then UPID. The first key picks the strain the
  literature was done on (S288C, 3D7, Nipponbare, HB8, PR8) without anyone
  naming it; the UPID makes it total.
* **A species with no reference proteome** is kept as a row, not dropped:
  `genome_only` with its best NCBI assembly (`assembly_rank`, ported from
  `../ip3r_genes/scripts/s4_manifest_lib.py:rank_key`), or `none`. Its
  families can then only be found by S5's genomic sweep, and the manifest
  says so.

Every API response is archived under `<data root>/raw_api/s4/`, so a rerun
is offline unless `--refresh`. Proteome files go to
`<data root>/proteomes/s4/`. No gene symbol is read for any purpose.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import requests  # noqa: E402

from src.utils.data_root import require_data_root  # noqa: E402

RELEASE = "2026_03"            # census v2's release (D31); the README must agree
RELEASE_DATE = "2026-09-03"    # the README's own Last-Modified; a file newer
                               # than this is a mid-release reissue (D35)
OUT_DIR = ROOT / "results" / "proteome_scope"
LIVE = ROOT / "results" / "session_live.json"

PROTEOMES_API = "https://rest.uniprot.org/proteomes/search"
FTP_ROOT = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
            "knowledgebase/reference_proteomes")
NCBI_GENOME = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"

FTP_DOMAIN = {"eukaryota": "Eukaryota", "bacteria": "Bacteria",
              "archaea": "Archaea", "viruses": "Viruses"}
LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}


# ------------------------------------------------------------------ paths
def raw_dir() -> Path:
    p = require_data_root() / "raw_api" / "s4"
    p.mkdir(parents=True, exist_ok=True)
    return p


def proteome_dir() -> Path:
    p = require_data_root() / "proteomes" / "s4"
    p.mkdir(parents=True, exist_ok=True)
    return p


# ------------------------------------------------------------------ network
def get(url: str, params: dict | None = None, timeout_s: int = 120,
        tries: int = 8, stream: bool = False) -> requests.Response:
    """GET with jittered backoff; a 404 is returned, anything else raises."""
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout_s, stream=stream)
            if r.status_code == 404:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r
        except requests.RequestException as exc:
            last = exc
            time.sleep(min(60, 2 ** i) * (0.5 + random.random()))
    raise RuntimeError(f"GET failed after {tries} tries: {url} ({last})")


def archived_json(name: str, fetch, refresh: bool = False):
    """Return the archived JSON `name`, fetching (and archiving) if absent."""
    path = raw_dir() / name
    if path.exists() and not refresh:
        return json.loads(path.read_text())
    data = fetch()
    path.write_text(json.dumps(data))
    return data


def fetch_reference_proteomes(taxon_id: int, refresh: bool = False) -> list[dict]:
    """Every UniProt reference proteome in the taxon's subtree (raw JSON)."""
    def fetch():
        r = get(PROTEOMES_API, {"query": f"(taxonomy_id:{taxon_id}) AND (reference:true)",
                                "format": "json", "size": 500})
        return r.json().get("results", [])
    return archived_json(f"proteomes_taxon_{taxon_id}.json", fetch, refresh)


def fetch_release_readme(refresh: bool = False) -> str:
    path = raw_dir() / f"reference_proteomes_README_{RELEASE}.txt"
    if path.exists() and not refresh:
        return path.read_text()
    text = get(f"{FTP_ROOT}/README").text
    path.write_text(text)
    return text


def fetch_assembly(accession: str, refresh: bool = False) -> dict | None:
    def fetch():
        # a proteome is often annotated on a superseded assembly version,
        # which the default (current-only) query silently omits
        r = get(f"{NCBI_GENOME}/accession/{accession}/dataset_report",
                {"filters.assembly_version": "all_assemblies"})
        return {} if r.status_code == 404 else r.json()
    reps = archived_json(f"assembly_{accession}.json", fetch, refresh).get("reports", [])
    return reps[0] if reps else None


def fetch_taxon_assemblies(taxon_id: int, refresh: bool = False) -> list[dict]:
    def fetch():
        r = get(f"{NCBI_GENOME}/taxon/{taxon_id}/dataset_report",
                {"page_size": 100, "filters.assembly_version": "current"})
        return {} if r.status_code == 404 else r.json()
    return archived_json(f"assemblies_taxon_{taxon_id}.json", fetch,
                         refresh).get("reports", [])


def fetch_metalink(upid: str, superregnum: str, refresh: bool = False) -> dict[str, str]:
    """{file name: md5} from the proteome directory's RELEASE.metalink."""
    def fetch():
        url = f"{FTP_ROOT}/{FTP_DOMAIN[superregnum]}/{upid}/RELEASE.metalink"
        text = get(url).text
        return dict(re.findall(r'<file name="([^"]+)">.*?<hash type="md5">([0-9a-f]+)</hash>',
                               text, flags=re.S))
    return archived_json(f"metalink_{upid}.json", fetch, refresh)


def fetch_last_modified(url: str, name: str, refresh: bool = False) -> str:
    """ISO date of a served file (HTTP Last-Modified), archived."""
    def fetch():
        for i in range(5):
            try:
                r = requests.head(url, timeout=60)
                r.raise_for_status()
                return {"last_modified": r.headers.get("Last-Modified", "")}
            except requests.RequestException:
                time.sleep(2 ** i)
        raise RuntimeError(f"HEAD failed: {url}")
    raw = archived_json(f"head_{name}.json", fetch, refresh)["last_modified"]
    return time.strftime("%Y-%m-%d", time.strptime(raw, "%a, %d %b %Y %H:%M:%S %Z")) if raw else ""


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ parsing
def parse_readme(text: str) -> tuple[str, dict[str, dict]]:
    """(release line, {upid: {taxid, superregnum, n_canonical, n_additional,
    n_gene2acc, name}}) from the release README's statistics table."""
    release = ""
    table: dict[str, dict] = {}
    in_table = False
    for line in text.splitlines():
        if line.startswith("Release ") and not release:
            release = line.strip()
        if line.startswith("Proteome_ID\t"):
            in_table = True
            continue
        if in_table and line.startswith("UP"):
            f = line.split("\t")
            table[f[0]] = {"taxid": int(f[1]), "superregnum": f[3],
                           "n_canonical": int(f[4]), "n_additional": int(f[5]),
                           "n_gene2acc": int(f[6]), "name": f[7]}
    return release, table


def flatten_proteome(rec: dict) -> dict:
    """The fields of one UniProt proteome record the manifest keeps (D9)."""
    tax = rec.get("taxonomy", {})
    busco = rec.get("proteomeCompletenessReport", {}).get("buscoReport", {})
    cpd = rec.get("proteomeCompletenessReport", {}).get("cpdReport", {})
    asm = rec.get("genomeAssembly", {})
    stats = rec.get("proteomeStatistics", {})
    busco_pct = (round(100 * busco["complete"] / busco["total"], 1)
                 if busco.get("total") else None)
    return {
        "upid": rec.get("id", ""),
        "proteome_taxid": tax.get("taxonId"),
        "proteome_organism": tax.get("scientificName", ""),
        "superregnum": rec.get("superkingdom", ""),
        "modified": rec.get("modified", ""),
        "gene_count": rec.get("geneCount") or 0,
        "protein_count": rec.get("proteinCount") or 0,
        "reviewed": stats.get("reviewedProteinCount") or 0,
        "busco_complete_pct": busco_pct,
        "busco_fragmented": busco.get("fragmented"),
        "busco_missing": busco.get("missing"),
        "busco_lineage": busco.get("lineageDb", ""),
        "cpd_status": cpd.get("status", ""),
        "annotation_source": (rec.get("genomeAnnotation") or {}).get("source", ""),
        "assembly_id": asm.get("assemblyId", ""),
        "assembly_source": asm.get("source", ""),
        "assembly_level_uniprot": asm.get("level", ""),
    }


def select_proteome(cands: list[dict]) -> dict | None:
    """The one declared rule: most reviewed → BUSCO → genes → UPID."""
    if not cands:
        return None
    return sorted(cands, key=lambda c: (-c["reviewed"],
                                        -(c["busco_complete_pct"] or 0),
                                        -c["gene_count"], c["upid"]))[0]


def flatten_assembly(rec: dict) -> dict:
    info = rec.get("assembly_info", {})
    stats = rec.get("assembly_stats", {})
    return {
        "accession": rec.get("accession", ""),
        "assembly_name": info.get("assembly_name", ""),
        "level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "release_date": info.get("release_date", ""),
        "total_length": int(stats.get("total_sequence_length", 0) or 0),
        "scaffold_n50": int(stats.get("scaffold_n50", 0) or 0),
        "contig_n50": int(stats.get("contig_n50", 0) or 0),
        "annotated": "annotation_info" in rec,
        "status": info.get("assembly_status", ""),
        "current_accession": rec.get("current_accession", ""),
    }


def assembly_rank(asm: dict) -> tuple:
    """annotated > RefSeq > level > scaffold N50 (the parent's rank_key)."""
    return (1 if asm["annotated"] else 0,
            1 if asm["accession"].startswith("GCF_") else 0,
            LEVEL_RANK.get(asm["level"], 0), asm["scaffold_n50"])


# ------------------------------------------------------------------ files
def proteome_files(upid: str, taxid: int, superregnum: str) -> dict[str, tuple[str, Path]]:
    """{kind: (url, local path)} for the canonical FASTA and gene2acc map."""
    base = f"{FTP_ROOT}/{FTP_DOMAIN[superregnum]}/{upid}/{upid}_{taxid}"
    local = proteome_dir() / f"{upid}_{taxid}"
    return {"fasta": (f"{base}.fasta.gz", Path(f"{local}.fasta.gz")),
            "gene2acc": (f"{base}.gene2acc.gz", Path(f"{local}.gene2acc.gz"))}


def gzip_ok(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass
        return True
    except (OSError, EOFError):
        return False


def download(url: str, dest: Path) -> str:
    """ok / cached / failed:<why>. Written via .part, validated as gzip."""
    if gzip_ok(dest):
        return "cached"
    r = get(url, stream=True, timeout_s=300)
    if r.status_code == 404:
        return "failed:404"
    tmp = dest.with_name(dest.name + ".part")
    with tmp.open("wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    tmp.rename(dest)
    if gzip_ok(dest):
        return "ok"
    dest.unlink(missing_ok=True)
    return "failed:gzip-invalid"


def fasta_accessions(path: Path) -> list[str]:
    """Accessions of a UniProt FASTA (`>sp|ACC|…`), in file order."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                out.append(line[1:].split("|", 2)[1])
    return out


def gene_groups(path: Path) -> dict:
    """Count gene2acc rows and gene groups. A group is column 3; a row whose
    group is '-' (gene unknown) is its own group — counted, not merged."""
    rows = unknown = 0
    groups: set[str] = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            rows += 1
            if f[2] == "-":
                unknown += 1
            else:
                groups.add(f[2])
    return {"rows": rows, "named_groups": len(groups), "unknown_gene_rows": unknown,
            "gene_groups": len(groups) + unknown}


def file_verdict(md5_ok: bool, readme_agrees: bool, reissued: bool,
                 entries: int, gene_count: int) -> str:
    """exact / reissued / FAIL:<why>. Integrity (metalink MD5) is never
    waived. The README's counts are waived only for a file served *after* the
    release date whose entry count equals the proteome's declared gene count
    — a documented mid-release reissue (D35), not a download we cannot explain."""
    if not md5_ok:
        return "FAIL:md5"
    if readme_agrees:
        return "exact"
    if reissued and entries == gene_count:
        return "reissued"
    return "FAIL:counts"
