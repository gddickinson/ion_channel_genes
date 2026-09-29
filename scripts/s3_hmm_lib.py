"""S3 shared pieces: paths, FASTA IO, HMMER table parsing, census loaders.

Ported from `../ip3r_genes/scripts/s3_hmm_lib.py` (method, not results) and
widened from two profiles to one per catalogue family. Nothing here reads a
gene symbol for any purpose other than *scoring* an assignment it did not
make (H15).

Bulk files (profiles, alignments, domtblouts, the non-redundant sweep FASTA)
live under `<data root>/hmmer/s3/`; committed tables under
`results/census_v3/`.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "census_v3"
LIVE = ROOT / "results" / "session_live.json"
REFERENCE_PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
S1_CALLS = ROOT / "results" / "benchmark_controls" / "calls.tsv"

csv.field_size_limit(1 << 30)


def s3_dir(*parts: str) -> Path:
    """`<data root>/hmmer/s3/<parts>` — raises if the drive is absent (D1)."""
    p = require_data_root() / "hmmer" / "s3"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def s2_raw() -> Path:
    return require_data_root() / "raw_api" / "s2"


# ------------------------------------------------------------------ FASTA
def _open(path: Path, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def iter_fasta(path: Path):
    """Yield (header-without->, sequence) from a plain or gzipped FASTA."""
    head, chunks = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if head is not None:
                    yield head, "".join(chunks)
                head, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if head is not None:
        yield head, "".join(chunks)


def read_fasta(path: Path) -> dict[str, str]:
    return {h.split()[0]: s for h, s in iter_fasta(path)}


def write_fasta(path: Path, items, width: int = 0) -> None:
    with _open(path, "wt") as fh:
        for name, seq in items:
            fh.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def seq_id(seq: str) -> str:
    """Stable id for an exact sequence (the non-redundant sweep key)."""
    return "s" + hashlib.md5(seq.encode()).hexdigest()[:16]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------ domtblout IO
DOMTBL_FIELDS = [
    "target_name", "target_acc", "tlen", "query_name", "query_acc", "qlen",
    "full_evalue", "full_score", "full_bias", "dom_num", "dom_of",
    "dom_cvalue", "dom_ivalue", "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
]
_INT = {"tlen", "qlen", "dom_num", "dom_of", "hmm_from", "hmm_to",
        "ali_from", "ali_to", "env_from", "env_to"}
_FLOAT = {"full_evalue", "full_score", "full_bias", "dom_cvalue",
          "dom_ivalue", "dom_score", "dom_bias", "acc"}

# A domain row counts towards profile coverage only if it is significant on
# its own. Low-scoring spurious domains would otherwise inflate coverage.
DOM_IVALUE_MAX = 1e-3


def iter_domtblout(path: Path):
    """One dict per domain row (plain or gzipped --domtblout)."""
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(maxsplit=len(DOMTBL_FIELDS))
            row = dict(zip(DOMTBL_FIELDS, parts))
            for k in _INT:
                row[k] = int(row[k])
            for k in _FLOAT:
                row[k] = float(row[k])
            yield row


def _merged(spans: list[tuple[int, int]]) -> int:
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(spans):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def hits_by_target(rows) -> dict[tuple[str, str], dict]:
    """(target, profile) → full-sequence score + merged profile coverage.

    The margin is measured on the **full-sequence** bit score; coverage is
    the fraction of the profile's match states spanned by the significant
    domains, merged so overlapping domains are not counted twice.
    """
    acc: dict[tuple[str, str], dict] = {}
    for r in rows:
        key = (r["target_name"], r["query_name"])
        h = acc.get(key)
        if h is None:
            h = acc[key] = {"target": r["target_name"],
                            "profile": r["query_name"],
                            "score": r["full_score"],
                            "evalue": r["full_evalue"],
                            "qlen": r["qlen"], "tlen": r["tlen"],
                            "hmm_spans": [], "ali_spans": []}
        if r["dom_ivalue"] <= DOM_IVALUE_MAX:
            h["hmm_spans"].append((r["hmm_from"], r["hmm_to"]))
            h["ali_spans"].append((r["ali_from"], r["ali_to"]))
    for h in acc.values():
        h["hmm_coverage"] = round(_merged(h.pop("hmm_spans")) / h["qlen"], 4)
        h["target_coverage"] = round(_merged(h.pop("ali_spans")) / h["tlen"], 4)
    return acc


# ------------------------------------------------------------------ tables
def read_tsv(path: Path) -> list[dict]:
    with _open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with _open(path, "wt") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
            n += 1
    return n


def iter_census_v2(columns: tuple[str, ...] | None = None):
    """Stream census v2 rows (28 MB gz, 1.25 M rows) from the data root."""
    with gzip.open(s2_raw() / "census_v2.tsv.gz", "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield {k: row[k] for k in columns} if columns else row


def census_v2_fasta() -> Path:
    return s2_raw() / "census_v2.fasta.gz"
