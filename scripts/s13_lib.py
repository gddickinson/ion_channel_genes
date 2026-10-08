"""s13_lib.py — S13 shared pieces (D55): anchors, paths, CDS validation, statistics.

The four anchors are the only gene names S13 declares, by UniProt accession
(D55 (1)); every set is derived from S7b's trees (`s13_sets.py`). Ported
from the IP3R project's S9 (`../ip3r_genes/scripts/s9_cds_lib.py`,
`s9_codeml_lib.py`, `s9_relax.py`): CDS validation by translation with every
disagreeing codon masked (IP3R D36), the χ² survival function, BH, and the
HyPhy JSON loader that repairs the bare `inf`/`nan` HyPhy writes.

Bulk (raw fetches, codeml/HyPhy work directories) under
`<data root>/selection/s13/`; committed tables under `results/selection/`.
"""

from __future__ import annotations

import json
import math
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "selection"
LIVE = ROOT / "results" / "session_live.json"
MEMBERS = ROOT / "results" / "alignments" / "members.tsv"
TREES = ROOT / "results" / "phylogeny" / "tier1"

#: D55 (1) — the clinical anchors: family → (UniProt accession, display name).
#: The display name is a label for tables and figures, never an input.
ANCHORS = {
    "cftr": ("P13569", "CFTR"),
    "nav": ("P35498", "SCN1A"),
    "kv_kcnq": ("P51787", "KCNQ1"),
    "ryr": ("P21817", "RYR1"),
}
HUMAN = "Homsap"
VERTEBRATE = "vertebrate"
CYCLOSTOME = {"Petmar"}          # D55 (8): the saturation sensitivity drop
DS_SATURATED = 3.0               # D55 (8)
OMEGA_BOUND = 999.0              # codeml's upper bound on ω (IP3R D38)
BS_OMEGA_STARTS = (0.5, 1.5, 2.5, 4.0)   # D55 (7), IP3R D37
ENV_BIN = Path("/opt/anaconda3/envs/piezo1/bin")

#: The standard code (NCBI table 1), built here so the selftest needs no
#: Biopython (`src.analysis` imports it).
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {a + b + c: _AA[16 * i + 4 * j + k]
                for i, a in enumerate(_BASES) for j, b in enumerate(_BASES)
                for k, c in enumerate(_BASES)}


def s13_dir(*parts: str) -> Path:
    """`<data root>/selection/s13/...` — raises without the drive (D1)."""
    p = require_data_root() / "selection" / "s13"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def tool_bin(name: str) -> str:
    """codeml, hyphy and pal2nal.pl live in the piezo1 env (D18); a missing
    one stops the run (D28)."""
    found = shutil.which(name)
    if found:
        return found
    env = ENV_BIN / name
    if env.exists():
        return str(env)
    raise SystemExit(f"{name} not found on PATH or in {ENV_BIN} (D28)")


def species_code(label: str) -> str:
    return label.rsplit("__", 1)[1]


def anchor_label(fam: str) -> str:
    return f"{ANCHORS[fam][0]}__{HUMAN}"


# ------------------------------------------------------------------ CDS
def translate(cds: str) -> str:
    """Full translation, stops as `*`, never truncated at one."""
    cds = cds.upper().replace("U", "T")
    return "".join(GENETIC_CODE.get(cds[i:i + 3], "X")
                   for i in range(0, len(cds) - len(cds) % 3, 3))


def validate_and_mask(cds: str, prot: str) -> tuple[str | None, dict]:
    """IP3R D36: accept iff the translation has the protein's length (after
    stripping a terminal stop) and ≥ 99 % identity at comparable positions;
    mask **every** disagreeing codon to NNN."""
    cds = cds.upper().replace("U", "T")
    cds = cds[: len(cds) - len(cds) % 3]
    tr = translate(cds)
    if tr.endswith("*") and len(tr) == len(prot) + 1:
        tr, cds = tr[:-1], cds[:-3]
    st = {"cds_len": len(cds), "prot_len": len(prot), "n_mismatch": "",
          "n_masked": "", "n_comparable": "", "reason": ""}
    if len(tr) != len(prot):
        st["reason"] = f"length {len(tr)} vs {len(prot)}"
        return None, st
    out, mism, masked, comp = list(cds), 0, 0, 0
    for i, (a, b) in enumerate(zip(tr, prot)):
        if b != "X" and a not in "*X":
            comp += 1
        if a == b:
            continue
        out[3 * i:3 * i + 3] = "NNN"
        masked += 1
        if a not in "*X" and b != "X":
            mism += 1
    st.update(n_mismatch=mism, n_masked=masked, n_comparable=comp)
    if not comp or mism / comp > 0.01:
        st["reason"] = f"identity {1 - mism / max(comp, 1):.4f} < 0.99"
        return None, st
    return "".join(out), st


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# ------------------------------------------------------------------ statistics
def chi2_sf(x: float, k: int) -> float:
    """Upper tail of χ²_k (regularised incomplete gamma Q), stdlib."""
    if x <= 0:
        return 1.0
    a, xx = k / 2.0, x / 2.0
    if xx < a + 1.0:
        term = total = 1.0 / a
        for n in range(1, 10_000):
            term *= xx / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, min(1.0, 1.0 - total * math.exp(-xx + a * math.log(xx) - math.lgamma(a))))
    tiny = 1e-300
    b, c, d = xx + 1.0 - a, 1.0 / tiny, 1.0 / (xx + 1.0 - a)
    h = d
    for i in range(1, 10_000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        h *= d * c
        if abs(d * c - 1.0) < 1e-15:
            break
    return max(0.0, min(1.0, h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))))


def lrt_p(lnl0: float, lnl1: float, df: int, mixture: bool = False) -> tuple[float, float]:
    """(2ΔlnL, p). `mixture`: the 50:50 χ²₀/χ²₁ boundary null (IP3R D37)."""
    stat = max(0.0, 2.0 * (lnl1 - lnl0))
    p = 0.5 * chi2_sf(stat, 1) if mixture else chi2_sf(stat, df)
    return stat, (1.0 if mixture and stat == 0 else p)


def bh(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg q-values, input order kept."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q, prev = [0.0] * n, 1.0
    for rank, idx in enumerate(reversed(order)):
        prev = min(prev, pvals[idx] * n / (n - rank))
        q[idx] = prev
    return q


_NONFINITE = re.compile(r"(?<=[:\[,])(\s*)(-?)inf\b|(?<=[:\[,])(\s*)nan\b")


def load_hyphy_json(path: Path) -> dict:
    """HyPhy JSON with its bare `inf` / `nan` literals made legal."""
    def sub(m: re.Match) -> str:
        if m.group(3) is not None:
            return f"{m.group(3)}NaN"
        return f"{m.group(1)}{m.group(2)}Infinity"
    return json.loads(_NONFINITE.sub(sub, path.read_text()))
