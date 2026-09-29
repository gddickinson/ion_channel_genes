"""s7_lib.py — S7 shared pieces: where tier-1 files live, alignment metrics,
trimAl column masks.

Tier-1 trees are built on S6's family L-INS-i alignments
(`<data root>/alignments/s6/family/<fam>.linsi.fasta`). Outgroup rows are
added with `mafft --add --keeplength`, so the family's columns never move,
and the trimming is a **column mask computed on the family alignment
alone** and then applied to the outgroup rows — the ingroup's trimmed
alignment is the same whether or not an outgroup is attached.

Bulk (alignments with outgroups, IQ-TREE output): `<data root>/trees/s7/`.
Committed: `results/phylogeny/tier1_*.tsv` and one treefile per family
under `results/phylogeny/tier1/`.
"""

from __future__ import annotations

import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta  # noqa: E402
from s6_lib import s6_dir  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "phylogeny"
TREE_DIR = OUT_DIR / "tier1"
LIVE = ROOT / "results" / "session_live.json"

#: The trimming candidates compared before any tree is built (S6 emergent row).
TRIM_METHODS = {
    "none": [],
    "automated1": ["-automated1"],
    "gappyout": ["-gappyout"],
    "gt0.5": ["-gt", "0.5"],
}


def s7_dir(*parts: str) -> Path:
    """`<data root>/trees/s7/...` — raises without the drive (D1)."""
    p = require_data_root() / "trees" / "s7"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def family_alignment(fam: str) -> Path:
    return s6_dir("family") / f"{fam}.linsi.fasta"


def read_aln(path: Path) -> list[tuple[str, str]]:
    return [(h.split()[0], s) for h, s in iter_fasta(path)]


def trim_columns(aln: Path, method: str) -> list[int]:
    """0-based columns trimAl keeps under `method` (all of them for "none")."""
    if method == "none":
        return list(range(len(read_aln(aln)[0][1])))
    p = subprocess.run(["trimal", "-in", str(aln), "-colnumbering",
                        *TRIM_METHODS[method]],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"trimal {method} failed on {aln.name}: "
                           f"{p.stderr.strip()[:300]}")
    m = re.search(r"#ColumnsMap\s*(.*)", p.stdout)
    if not m:
        raise RuntimeError(f"trimal {method}: no #ColumnsMap for {aln.name}")
    return [int(x) for x in m.group(1).replace(",", " ").split()]


def apply_mask(rows: list[tuple[str, str]], cols: list[int]) -> list[tuple[str, str]]:
    return [(label, "".join(seq[c] for c in cols)) for label, seq in rows]


def is_gap(ch: str) -> bool:
    return ch in "-.Xx?*"


def metrics(rows: list[tuple[str, str]], cols: list[int]) -> dict:
    """Columns, parsimony-informative sites, gap fraction, residue retention."""
    n = len(rows)
    inf = 0
    gaps = 0
    for c in cols:
        states = Counter(s[c].upper() for _, s in rows if not is_gap(s[c]))
        gaps += n - sum(states.values())
        if sum(1 for v in states.values() if v >= 2) >= 2:
            inf += 1
    kept = set(cols)
    ret = []
    for _, s in rows:
        total = sum(1 for ch in s if not is_gap(ch))
        k = sum(1 for i, ch in enumerate(s) if i in kept and not is_gap(ch))
        ret.append(k / total if total else 0.0)
    ret.sort()
    return {"cols": len(cols), "informative": inf,
            "gap_frac": round(gaps / (n * len(cols)), 4) if cols else 1.0,
            "retention_median": round(ret[len(ret) // 2], 4),
            "retention_min": round(ret[0], 4)}


def write_aln(path: Path, rows: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for label, seq in rows:
            fh.write(f">{label}\n{seq}\n")
