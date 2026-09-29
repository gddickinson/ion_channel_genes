"""S3b shared pieces: paths, the panel sweep DB and its header parsing.

S3b sweeps S3a's frozen profile library over S4's declared denominator — the
concatenated reference proteomes of the 52-species panel (D35) — and runs
jackhmmer from one derived seed per census family (D10). Bulk output lives
under `<data root>/hmmer/s3b/`; committed tables under `results/panel_sweep/`.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import read_tsv  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "panel_sweep"
MANIFEST = ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv"

_OX = re.compile(r"\bOX=(\d+)")


def s3b_dir(*parts: str) -> Path:
    """`<data root>/hmmer/s3b/<parts>` — raises if the drive is absent (D1)."""
    p = require_data_root() / "hmmer" / "s3b"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def panel_db() -> Path:
    """S4's concatenated sweep DB (822,499 canonical entries, 50 proteomes)."""
    return require_data_root() / "proteomes" / "s4" / "panel_refprot.fasta"


def accession(target: str) -> str:
    """`sp|P00350|6PGD_ECOLI` → `P00350` (the HMMER target name is the
    first word of the FASTA header)."""
    parts = target.split("|")
    return parts[1] if len(parts) >= 3 else target


def header_taxid(header: str) -> str:
    m = _OX.search(header)
    return m.group(1) if m else ""


def species_by_taxid() -> dict[str, dict]:
    """proteome taxon id → manifest row (species, group, upid), for the
    rows S4 declared as proteomes. The proteome's taxon (e.g. E. coli K-12,
    83333) is what the FASTA's OX= carries, not the panel taxon."""
    out = {}
    for r in read_tsv(MANIFEST):
        if r["status"] == "proteome":
            out[r["proteome_taxid"]] = r
    return out


# ------------------------------------------------------------- jackhmmer log
# Ported from ../ip3r_genes/scripts/s3_hmm_lib.py. jackhmmer's log has no
# "@@ Round: 1" marker: each round ENDS with "@@ New targets included: X",
# so those lines, in order, are rounds 1..k. Each round's score table lists
# the included targets above a "------ inclusion threshold ------" divider;
# those per-round lists are what D10 is evaluated on and exist only here.
_NEWT_RE = re.compile(r"^@@ New targets included:\s+(\d+)")
_TOTT_RE = re.compile(r"^@@ Included in MSA:.*from (\d+) targets")
_CONV_RE = re.compile(r"CONVERGED", re.IGNORECASE)
_INCL_DIVIDER = "inclusion threshold"
_SCORES_HDR = re.compile(r"^Scores for complete sequences")
_TABLE_END = ("Domain annotation", "Internal pipeline", "Query:", "//")
_NAME_COL = 8


def _looks_numeric(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _parse_score_row(stripped: str) -> str | None:
    toks = stripped.split()
    if toks and toks[0] in ("+", "-"):
        toks = toks[1:]
    if len(toks) > _NAME_COL and _looks_numeric(toks[0]):
        return toks[_NAME_COL]
    return None


def parse_jackhmmer_log(path: Path) -> dict:
    """{"rounds": [{"round", "new_targets", "targets_in_msa", "included"}],
    "converged", "trailing_incomplete_round"}. A score table with no closing
    "New targets" line is an unfinished round and is reported separately,
    never appended as a round with 0 new targets — that is what a converged
    final round looks like."""
    new_targets: list[int] = []
    cumulative: list[int] = []
    per_round: list[list[str]] = []
    converged = in_table = below = False
    current: list[str] = []
    with path.open(errors="replace") as f:
        for line in f:
            if _SCORES_HDR.match(line):
                in_table, below = True, False
                continue
            if in_table:
                s = line.strip()
                if not s:
                    continue
                if _INCL_DIVIDER in s:
                    below = True
                    continue
                if s.startswith(_TABLE_END):
                    in_table = False
                elif not below:
                    name = _parse_score_row(s)
                    if name:
                        current.append(name)
                    continue
                else:
                    continue
            m = _NEWT_RE.match(line)
            if m:
                new_targets.append(int(m.group(1)))
                per_round.append(current)
                current = []
                continue
            m = _TOTT_RE.match(line)
            if m:
                cumulative.append(int(m.group(1)))
                continue
            if _CONV_RE.search(line):
                converged = True
    rounds = [{"round": i, "new_targets": nt,
               "targets_in_msa": cumulative[i - 1] if i - 1 < len(cumulative) else "",
               "included": per_round[i - 1]}
              for i, nt in enumerate(new_targets, 1)]
    return {"rounds": rounds, "converged": converged,
            "trailing_incomplete_round": len(current)}
