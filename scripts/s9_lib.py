"""s9_lib.py — S9 shared pieces: the filter coordinate system and the
congruence statistics (D49).

**One coordinate system.** Every P-loop pore module S6 extracted is added to
S8's untrimmed P-loop tier-2 L-INS-i alignment with `--keeplength`, so a
column means the same thing for a Kv module, a Nav repeat and a TRP module.
Two anchors fix where the filter is:

* the **K window** — the five columns KcsA's T75–G79 (`TVGYG`) occupy;
* the **repeat locus** — the column each of Nav1.5's D372 / E898 / K1419 /
  A1711 occupies (one per repeat), read in a module at its own repeat's
  column.

**Congruence.** A character (one filter string per tip; `None` = missing) is
laid on a tree and scored by Fitch parsimony — exact on a binary tree, and an
IQ-TREE tree is binary once its drawn root trifurcation is read as an
unrooted node, which leaves the length unchanged — against tip-label
permutations, with the retention index beside it. Missing tips are pruned
before any clade is read.

Bulk: `<data root>/alignments/s9/`. Committed: `results/filter_atlas/`.
"""

from __future__ import annotations

import random
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta, read_tsv  # noqa: E402
from s7_newick import Node  # noqa: E402
from src.classify.motifs import FOUR_REPEAT_ANCHOR  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "filter_atlas"
ALN_DIR = ROOT / "results" / "alignments"
PHY_DIR = ROOT / "results" / "phylogeny"
LIVE = ROOT / "results" / "session_live.json"

UNIT = "ploop"
KCSA = ("P0A334", 75, "TVGYG")          # accession, first position, motif
KCSA_ROW = "kcsa_prok__EX_Sl_KcsA"      # its module row in the tier-2 alignment
NAV15_LABEL = "Q14524__Homsap"           # D39 label of human Nav1.5
REPEAT_FAMILIES = FOUR_REPEAT_ANCHOR.applies_to
N_PERM = 1000
SEED = 1
MISSING = re.compile(r"[-X?]")


def s9_dir(*parts: str) -> Path:
    """`<data root>/alignments/s9/...` — raises without the drive (D1)."""
    p = require_data_root() / "alignments" / "s9"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_aln(path: Path) -> list[tuple[str, str]]:
    return [(h.split()[0], s.upper()) for h, s in iter_fasta(path)]


def col_of_residue(row: str, pos: int) -> int | None:
    """0-based alignment column of the row's `pos`-th residue (1-based)."""
    n = 0
    for i, ch in enumerate(row):
        if ch not in "-.":
            n += 1
            if n == pos:
                return i
    return None


def residue_index(row: str, col: int) -> int | None:
    """1-based residue number at a column (None if the row has a gap there)."""
    if row[col] in "-.":
        return None
    return sum(1 for ch in row[:col + 1] if ch not in "-.")


def module_rows() -> list[dict]:
    """S6's P-loop module rows (one per member × module)."""
    return [r for r in read_tsv(ALN_DIR / "modules.tsv") if r["unit"] == UNIT]


def module_key(family: str, label: str, module: str, n_modules: int) -> str:
    """The key S6's `<unit>.modules.fasta` uses: `fam|label` or `fam|label_mN`."""
    return f"{family}|{label}" if n_modules == 1 else f"{family}|{label}_m{module}"


def tip_to_key(tip: str) -> str:
    """Tier-2 tip `fam__ACC__Spe[_mN]` → module key `fam|ACC__Spe[_mN]`."""
    fam, rest = tip.split("__", 1)
    return f"{fam}|{rest}"


def is_missing(state: str | None) -> bool:
    return state is None or state == "" or bool(MISSING.search(state))


# ----------------------------------------------------------------- trees

def prune(tree: Node, keep: set[str]) -> Node | None:
    """The subtree on `keep` leaves; unary nodes collapse into their child."""
    if not tree.children:
        return Node(tree.name, tree.length) if tree.name in keep else None
    kids = [k for k in (prune(c, keep) for c in tree.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        return kids[0]
    return Node(tree.name, tree.length, kids)


def _binary_post(tree: Node) -> list[tuple[int, list[int]]]:
    """Post-order (node id, child ids) with every multifurcation resolved
    left-to-right — exact for Fitch only where the multifurcation is a drawn
    unrooted root, which is the only one IQ-TREE writes."""
    order: list[tuple[int, list[int]]] = []
    leaves: dict[int, str] = {}
    counter = [0]

    def walk(n: Node) -> int:
        if not n.children:
            i = counter[0]
            counter[0] += 1
            leaves[i] = n.name
            return i
        ids = [walk(c) for c in n.children]
        while len(ids) > 2:
            i = counter[0]
            counter[0] += 1
            order.append((i, ids[:2]))
            ids = [i] + ids[2:]
        i = counter[0]
        counter[0] += 1
        order.append((i, ids))
        return i

    walk(tree)
    return order, leaves


def multifurcations(tree: Node, top: bool = True) -> int:
    """Internal nodes below the drawn root with more than two children."""
    own = 0 if top else int(len(tree.children) > 2)
    return own + sum(multifurcations(c, False) for c in tree.children)


class Fitch:
    """Fitch parsimony on one tree, re-scorable under permuted tip states."""

    def __init__(self, tree: Node):
        self.order, self.leaves = _binary_post(tree)

    def length(self, states: dict[str, str | None]) -> int:
        idx = {s: 1 << i for i, s in enumerate(sorted({v for v in states.values()
                                                         if not is_missing(v)}))}
        full = (1 << len(idx)) - 1
        sets = {i: (idx[states[n]] if not is_missing(states.get(n)) else full)
                for i, n in self.leaves.items()}
        cost = 0
        for i, kids in self.order:
            a, b = sets[kids[0]], sets[kids[1]]
            if a & b:
                sets[i] = a & b
            else:
                sets[i] = a | b
                cost += 1
        return cost


def congruence(tree: Node, states: dict[str, str | None]) -> dict:
    """Fitch length, minimum and star-tree lengths, RI, permutation p."""
    known = {t: s for t, s in states.items() if not is_missing(s)}
    tips = set(tree.leaves())
    known = {t: s for t, s in known.items() if t in tips}
    if len(known) < 3:
        return {"n_tips": len(known), "n_states": len(set(known.values()))}
    sub = prune(tree, set(known))
    fitch = Fitch(sub)
    obs = fitch.length(known)
    counts = Counter(known.values())
    m = len(counts) - 1
    g = len(known) - max(counts.values())
    out = {"n_tips": len(known), "n_states": len(counts), "changes": obs,
           "min_changes": m, "star_changes": g,
           "ri": round((g - obs) / (g - m), 4) if g > m else None}
    if g == m:           # one state, or one tip per extra state: nothing to test
        return out
    rng = random.Random(SEED)
    names, vals = list(known), list(known.values())
    null = []
    for _ in range(N_PERM):
        rng.shuffle(vals)
        null.append(fitch.length(dict(zip(names, vals))))
    out["null_median"] = sorted(null)[len(null) // 2]
    out["null_min"] = min(null)
    out["p_perm"] = round((1 + sum(1 for x in null if x <= obs)) / (N_PERM + 1), 4)
    return out


def origins(tree: Node, carriers: set[str], exclude: frozenset) -> int:
    """Maximal carrier-only clades in the rooted sense (tree drawn so that
    `exclude` — the outgroup — is outside every clade read)."""
    cl = [frozenset([t]) for t in carriers]
    from s8_parse import clades  # noqa: PLC0415 — s8_parse imports s8_lib at load
    cl += [s for s, _ in clades(tree, exclude) if s <= carriers]
    maximal = {s for s in cl if not any(s < o for o in cl)}
    return len(maximal)
