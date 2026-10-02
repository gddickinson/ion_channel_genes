"""s8_lib.py — S8 shared pieces: tier-2 units, per-cell identity, the
representative rule (D8, D48).

**Units.** A tier-2 tree is built for every *alignable* superfamily with two
or more census families (`units()`); a single-family superfamily's tier-2
tree is its tier-1 tree (S7) and is not rebuilt. Units with a declared
`module_rule` (P-loop, iGluR, Ca²⁺-release, innexin clan) are built on S6's
extracted pore modules — one tip per module, so a four-repeat chain gives
four tips — and every output calls them *pore-module trees* (D27). Units
without one (Cys-loop, P2X, DEG/ENaC) are full-length, as the catalogue
declares them alignable end to end.

**Representatives (D8).** A cell is family × panel group × module index.
Within a cell the members are ordered by *centrality* (mean identity to the
rest of the cell) and clustered greedily: each member joins the first
representative it matches at ≥ the threshold, else becomes one. Identity is
measured inside the family, where it is meaningful: over the module's match
states of S6's `hmmalign` to the family profile (module units) or over the
family's L-INS-i alignment columns (full-length units), counting only
columns where both sequences have a residue and requiring ≥ 50 % of the
shorter to be covered (else the pair is "not similar"). Central-first
ordering keeps a long-branch outlier from becoming the face of a cluster.

Bulk: `<data root>/trees/s8/`. Committed: `results/phylogeny/tier2_*.tsv`.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta, read_tsv  # noqa: E402
from s6_lib import s6_dir  # noqa: E402
from src.catalogue import CATALOGUE, SUPERFAMILIES, census_families  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "phylogeny"
ALN_DIR = ROOT / "results" / "alignments"
LIVE = ROOT / "results" / "session_live.json"
MIN_PAIR_COVER = 0.5


def s8_dir(*parts: str) -> Path:
    """`<data root>/trees/s8/...` — raises without the drive (D1)."""
    p = require_data_root() / "trees" / "s8"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def _census_keys() -> set[str]:
    fams = census_families()
    return {f if isinstance(f, str) else f.key for f in fams}


def unit_families(sf: str) -> list[str]:
    keys = _census_keys()
    return sorted(k for k, f in CATALOGUE.items()
                  if f.superfamily == sf and k in keys)


def units() -> list[str]:
    """Alignable superfamilies with ≥ 2 census families — the tier-2 units."""
    return sorted(k for k, sf in SUPERFAMILIES.items()
                  if sf.alignable and len(unit_families(k)) >= 2)


def single_family_units() -> list[str]:
    return sorted(k for k, sf in SUPERFAMILIES.items()
                  if sf.alignable and len(unit_families(k)) == 1)


def is_module_unit(sf: str) -> bool:
    return bool(SUPERFAMILIES[sf].module_rule)


def members() -> dict[tuple[str, str], dict]:
    return {(r["family"], r["label"]): r
            for r in read_tsv(ALN_DIR / "members.tsv") if r["verdict"] == "include"}


def base_label(label: str) -> str:
    return re.sub(r"_m\d+$", "", label)


# ---------------------------------------------------------------- identity

def _a2m_states(a2m: str) -> str:
    """Match-state string of an A2M row (inserts dropped, '-' kept)."""
    return "".join(c for c in a2m if c.isupper() or c == "-")


def module_strings(fam: str, spans: dict[str, tuple[int, int]]) -> dict[str, str]:
    """label_mN → its module's match states, from S6's hmmalign of the family."""
    rows = {h.split()[0]: _a2m_states(s) for h, s in
            iter_fasta(s6_dir("project") / f"{fam}.members.a2m")}
    out = {}
    for lab, (k0, k1) in spans.items():
        st = rows.get(base_label(lab))
        if st is not None:
            out[lab] = st[k0 - 1:k1]
    return out


def family_alignment_rows(fam: str) -> dict[str, str] | None:
    """The family's L-INS-i rows; None for a family S6 did not align (< 4)."""
    p = s6_dir("family") / f"{fam}.linsi.fasta"
    if not p.exists():
        return None
    return {h.split()[0]: s.upper() for h, s in iter_fasta(p)}


def identity(a: str, b: str) -> float | None:
    """Identity over mutually covered columns; None if coverage is too low."""
    both = same = 0
    na = nb = 0
    for x, y in zip(a, b):
        gx, gy = x in "-.", y in "-."
        na += not gx
        nb += not gy
        if gx or gy:
            continue
        both += 1
        same += x == y
    short = min(na, nb)
    if not short or both / short < MIN_PAIR_COVER:
        return None
    return same / both


def matrix(seqs: dict[str, str]) -> dict[tuple[str, str], float]:
    labs = sorted(seqs)
    m = {}
    for i, a in enumerate(labs):
        for b in labs[i + 1:]:
            v = identity(seqs[a], seqs[b])
            m[(a, b)] = m[(b, a)] = 0.0 if v is None else v
    return m


def centrality(labs: list[str], m: dict) -> dict[str, float]:
    if len(labs) == 1:
        return {labs[0]: 1.0}
    return {a: sum(m[(a, b)] for b in labs if b != a) / (len(labs) - 1)
            for a in labs}


def greedy(labs: list[str], m: dict, threshold: float) -> list[tuple[str, list[str]]]:
    """Central-first greedy clustering → [(representative, members)]."""
    cen = centrality(labs, m)
    order = sorted(labs, key=lambda x: (-cen[x], x))
    reps: list[tuple[str, list[str]]] = []
    for lab in order:
        for rep, mem in reps:
            if m[(lab, rep)] >= threshold:
                mem.append(lab)
                break
        else:
            reps.append((lab, [lab]))
    return reps


def cells(sf: str) -> dict[tuple[str, str, str], dict[str, str]]:
    """(family, group, module) → {tip label: comparable string} for a unit."""
    mem = members()
    out: dict[tuple, dict[str, str]] = defaultdict(dict)
    if is_module_unit(sf):
        rows = [r for r in read_tsv(ALN_DIR / "modules.tsv")
                if r["unit"] == sf and r["quality"] == "full"]
        by_fam: dict[str, dict] = defaultdict(dict)
        for r in rows:
            lab = r["label"]
            by_fam[r["family"]][(lab, r["module"])] = (int(r["k_start"]), int(r["k_end"]))
        for fam, d in by_fam.items():
            strings = module_strings(fam, {lab: ks for (lab, _), ks in d.items()})
            for (lab, mod), _ in d.items():
                g = mem[(fam, base_label(lab))]["group"]
                out[(fam, g, mod)][lab] = strings[lab]
    else:
        for fam in unit_families(sf):
            rows = family_alignment_rows(fam)
            for (f, lab), r in mem.items():
                if f != fam:
                    continue
                if rows is None:      # unaligned family: every member its own cell
                    out[(fam, r["group"], f"0:{lab}")][lab] = ""
                elif lab in rows:
                    out[(fam, r["group"], "0")][lab] = rows[lab]
    return dict(out)
