"""s8_parse.py — S8b: read the tier-2 trees with the tests fixed in D48.

Called as `python3 scripts/s8_tier2.py parse`. Writes

* `results/phylogeny/tier2_trees.tsv` — one row per unit: model, outgroup
  one clade or not, root UFBoot, the ingroup's basal split, support summary;
* `results/phylogeny/tier2_families.tsv` — one row per family (and, in pore-
  module units, per family × module): its tips one clade or not, the clade's
  UFBoot, and the fewest other tips that would have to join it to make it
  one (`intruders`, with their families) — a descriptive near-miss measure;
* `results/phylogeny/tier2_placement.tsv` — where each tip the kingdom rule
  took out of the outgroup falls (the four animal `plgic_prok` members): the
  nested clades above it, up to the first that holds another family, with
  their composition and UFBoot;
* treefiles under `results/phylogeny/tier2/`.

"One clade" is read in the rooted sense when the outgroup is one clade (a
side of an edge that holds no outgroup tip) and in the unrooted sense
otherwise (either side of an edge) — every row records which. Nothing here
changes a D48 rule after a tree was read.
"""

from __future__ import annotations

import json
import re
import shutil
from collections import Counter

from s3_hmm_lib import read_tsv, write_tsv
from s7_newick import Node, ingroup_support, parse, root_partition, split_support, splits
from s8_lib import OUT_DIR, s8_dir

TREE_DIR = OUT_DIR / "tier2"
TREE_FIELDS = ["unit", "kind", "n_tips", "n_families", "cols", "informative", "model",
               "root_family", "n_outgroup", "outgroup_monophyletic", "root_ufboot",
               "rooted", "root_clade_sizes", "root_clade_ufboot", "root_small_clade",
               "internal_edges", "ufboot_median", "frac_ge95", "frac_lt70", "seconds"]
FAM_FIELDS = ["unit", "family", "module", "n_tips", "sense", "one_clade", "ufboot",
              "n_intruders", "intruder_families"]
PLACE_FIELDS = ["unit", "tip", "level", "clade_size", "ufboot", "composition"]
MAX_LEVELS = 3     # clades reported above a placed tip, from the first mixed one


def family_of(tip: str) -> str:
    return tip.split("__", 1)[0]


def module_of(tip: str) -> str:
    m = re.search(r"_m(\d+)$", tip)
    return f"m{m.group(1)}" if m else ""


def _model(path) -> str:
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", path.read_text())
    return m.group(1) if m else ""


def clades(tree: Node, exclude: frozenset) -> list[tuple[frozenset, float | None]]:
    """Every edge side holding no `exclude` tip (both sides when exclude is empty)."""
    leaves = frozenset(tree.leaves())
    out = []
    for side, sup in splits(tree):
        for s in (side, leaves - side):
            if len(s) > 1 and not (s & exclude):
                out.append((s, sup))
    return out


def clade_test(tree: Node, group: set[str], exclude: frozenset,
               cl: list) -> tuple[bool, float | None, frozenset]:
    """(one clade, its UFBoot, fewest other tips sharing a clade with the group)."""
    g = frozenset(group)
    if len(g) == 1:
        return True, None, frozenset()
    if exclude:
        present = [(s, u) for s, u in cl if s == g]
        one, sup = (True, present[0][1]) if present else (False, None)
    else:
        one, sup = split_support(tree, g)
    if one:
        return True, sup, frozenset()
    best = frozenset(tree.leaves()) - g - exclude
    for s, _ in cl:
        if g <= s and len(s - g) < len(best):
            best = s - g
    return False, None, best


def _comp(tips) -> str:
    c = Counter(family_of(t) for t in tips)
    return ",".join(f"{f}:{n}" for f, n in sorted(c.items(), key=lambda x: (-x[1], x[0])))


def placement(unit: str, tree: Node, tips: list[str], exclude: frozenset,
              cl: list) -> list[dict]:
    """The nested clades above each placed tip, smallest first, until
    MAX_LEVELS of them hold a tip outside the placed set."""
    placed = set(tips)
    rows = []
    for tip in tips:
        above = sorted((x for x in cl if tip in x[0]), key=lambda x: len(x[0]))
        level, mixed = 0, 0
        for s, u in above:
            if mixed >= MAX_LEVELS:
                break
            level += 1
            if s - placed:
                mixed += 1
            rows.append({"unit": unit, "tip": tip, "level": level, "clade_size": len(s),
                         "ufboot": "" if u is None else f"{u:g}",
                         "composition": _comp(s - {tip})})
    return rows


def unit_rows(row: dict) -> tuple[dict, list[dict], list[dict]]:
    sf = row["unit"]
    d = s8_dir("iqtree", sf)
    run = json.loads((d / "run.json").read_text())
    if run["input_sha256"] != row["input_sha256"]:
        raise RuntimeError(f"{sf}: tree built on a different input")
    tree = parse((d / f"{sf}.treefile").read_text())
    TREE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy(d / f"{sf}.treefile", TREE_DIR / f"{sf}.treefile")
    og_rec = json.loads((s8_dir("unit") / f"{sf}.outgroup.json").read_text())
    og = set(og_rec["outgroup"])
    leaves = set(tree.leaves())
    if not og <= leaves:
        raise RuntimeError(f"{sf}: outgroup tips missing from the tree")
    mono, sup = split_support(tree, og) if og else (False, None)
    part = root_partition(tree, og) if mono else None
    rooted = bool(part)
    exclude = frozenset(og) if rooted else frozenset()
    cl = clades(tree, exclude)
    s = ingroup_support(tree, og)
    small = min(part, key=lambda q: len(q[0]))[0] if part else frozenset()
    trow = {"unit": sf, "kind": row["kind"], "n_tips": len(leaves),
            "n_families": len({family_of(t) for t in leaves}),
            "cols": row["cols"], "informative": row["informative"],
            "model": _model(d / f"{sf}.iqtree"), "root_family": og_rec["root_family"],
            "n_outgroup": len(og), "outgroup_monophyletic": mono if og else "",
            "root_ufboot": "" if sup is None else f"{sup:g}", "rooted": rooted,
            "root_clade_sizes": ",".join(str(len(p)) for p, _ in part) if part else "",
            "root_clade_ufboot": ",".join("" if u is None else f"{u:g}"
                                          for _, u in part) if part else "",
            "root_small_clade": _comp(small) if part else "",
            **{k: ("" if v is None else v) for k, v in s.items()},
            "seconds": run["seconds"]}
    groups: dict[tuple[str, str], set] = {}
    for t in leaves:
        groups.setdefault((family_of(t), ""), set()).add(t)
        if module_of(t):
            groups.setdefault((family_of(t), module_of(t)), set()).add(t)
    frows = []
    for (fam, mod), g in sorted(groups.items()):
        if fam == og_rec["root_family"] and rooted:
            continue      # the root family is the outgroup; its test is the root row
        one, u, intr = clade_test(tree, g, exclude, cl)
        frows.append({"unit": sf, "family": fam, "module": mod, "n_tips": len(g),
                      "sense": "rooted" if rooted else "unrooted", "one_clade": one,
                      "ufboot": "" if u is None else f"{u:g}", "n_intruders": len(intr),
                      "intruder_families": _comp(intr)})
    prows = placement(sf, tree, sorted(og_rec["excluded"]), exclude, cl)
    return trow, frows, prows


def cmd_parse(_a) -> None:
    trees, fams, place = [], [], []
    for row in read_tsv(OUT_DIR / "tier2_inputs.tsv"):
        if not (s8_dir("iqtree", row["unit"]) / f"{row['unit']}.treefile").exists():
            print(f"{row['unit']}: no tree yet", flush=True)
            continue
        t, f, p = unit_rows(row)
        trees.append(t)
        fams += f
        place += p
        print(f"{t['unit']:13s} {t['model']:16s} og one clade={t['outgroup_monophyletic']!s:5} "
              f"root UFBoot={t['root_ufboot'] or '-':>4}  split={t['root_clade_sizes'] or '-':9}"
              f" families one clade {sum(r['one_clade'] for r in f if not r['module'])}"
              f"/{sum(not r['module'] for r in f)}", flush=True)
    write_tsv(OUT_DIR / "tier2_trees.tsv", TREE_FIELDS, trees)
    write_tsv(OUT_DIR / "tier2_families.tsv", FAM_FIELDS, fams)
    write_tsv(OUT_DIR / "tier2_placement.tsv", PLACE_FIELDS, place)
    print(f"{len(trees)} trees → {OUT_DIR / 'tier2_trees.tsv'}")
