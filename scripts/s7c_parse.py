"""s7c_parse.py — S7c/S7d: read the re-rooted trees against the root criterion.

Called as `python3 scripts/s7c_reroot.py parse`. Writes

* `results/phylogeny/tier1_reroot_trees.tsv` — one row per tree (family ×
  root set): model, outgroup clade, root UFBoot, the ingroup's basal split
  and its support, ingroup support summary;
* `results/phylogeny/tier1_reroot.tsv` — one row per family: the three
  criteria (D47), the verdict (`resolved` / `unresolved` and why), and the
  S7b KcsA-rooted result alongside;
* treefiles under `results/phylogeny/tier1/reroot/`.

"Same root" is exact: the ingroup's basal split, as a set of leaf sets, is
identical in the two trees. A near miss is reported (the smaller basal
clade's Jaccard between trees) but never accepted.
"""

from __future__ import annotations

import json
import re
import shutil

from s3_hmm_lib import read_tsv, write_tsv
from s7_lib import OUT_DIR, TREE_DIR
from s7_newick import (ingroup_support, outgroup_intruders, parse, root_partition,
                       split_support)
from s7c_basal import s7c_dir
from s7c_reroot import MIN_ROOT_UFBOOT, tree_id

TREE_FIELDS = ["family", "root_set", "n_ingroup", "n_outgroup", "cols",
               "informative", "model", "outgroup_monophyletic", "root_ufboot",
               "root_clades", "root_clade_sizes", "root_clade_ufboot",
               "root_small_clade_basal", "og_intruders", "og_intruders_basal",
               "og_intruder_labels", "internal_edges", "ufboot_median", "frac_ge95", "seconds"]
FAM_FIELDS = ["family", "root_sets", "n_ingroup", "n_basal", "a_one_clade",
              "b_root_ufboot_ge95", "c_same_root", "smaller_clade_jaccard",
              "verdict", "reason", "root_split_sizes",
              "s7b_outgroup", "s7b_one_clade", "s7b_root_ufboot",
              "s7b_split_on_shared"]


def _facts(path) -> dict:
    txt = path.read_text()

    def grab(pat):
        m = re.search(pat, txt)
        return m.group(1).strip() if m else ""
    return {"model": grab(r"Best-fit model according to BIC:\s*(\S+)")}


def _label(side: frozenset) -> str:
    """A short readable name for a basal clade: its size and two members."""
    return f"{len(side)}:{','.join(sorted(side)[:2])}"


def _restrict(part, keep: set) -> set[frozenset]:
    return {frozenset(s & keep) for s, _ in part if s & keep}


def basal_picks(fam: str) -> set[str]:
    """Labels of the S7c basal sequences added to a family's set (D47 (2))."""
    return {r["label"] for r in read_tsv(OUT_DIR / "reroot_basal.tsv")
            if r["family"] == fam and r["verdict"] == "add"}


def tree_row(row: dict) -> tuple[dict, list | None]:
    tid = tree_id(row)
    d = s7c_dir("iqtree", tid)
    run = json.loads((d / "run.json").read_text())
    if run["input_sha256"] != row["input_sha256"]:
        raise RuntimeError(f"{tid}: tree built on a different input")
    tree = parse((d / f"{tid}.treefile").read_text())
    og = set(row["outgroup"].split(","))
    mono, sup = split_support(tree, og)
    part = root_partition(tree, og) if mono else None
    dest = TREE_DIR / "reroot"
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy(d / f"{tid}.treefile", dest / f"{tid}.treefile")
    s = ingroup_support(tree, og)
    basal = basal_picks(row["family"])
    intr = outgroup_intruders(tree, og)
    small = min(part, key=lambda q: len(q[0]))[0] if part else frozenset()
    return ({"family": row["family"], "root_set": row["root_set"],
             "n_ingroup": row["n_ingroup"], "n_outgroup": len(og),
             "cols": row["cols"], "informative": row["informative"],
             **_facts(d / f"{tid}.iqtree"),
             "outgroup_monophyletic": mono,
             "root_ufboot": "" if sup is None else sup,
             "root_clades": "" if part is None else ";".join(_label(p) for p, _ in part),
             "root_clade_sizes": "" if part is None else ",".join(str(len(p)) for p, _ in part),
             "root_clade_ufboot": "" if part is None else ",".join(
                 "" if u is None else f"{u:g}" for _, u in part),
             "root_small_clade_basal": f"{len(small & basal)}/{len(small)}" if part else "",
             "og_intruders": len(intr),
             "og_intruders_basal": len(intr & basal),
             "og_intruder_labels": ",".join(sorted(intr)[:5]) if 0 < len(intr) <= 5 else "",
             **{k: ("" if v is None else v) for k, v in s.items()},
             "seconds": run["seconds"]}, part)


def s7b_root(fam: str, shared: set) -> dict:
    """The S7b (KcsA/MthK/NaK) result, its root split restricted to shared leaves."""
    old = {r["family"]: r for r in read_tsv(OUT_DIR / "tier1_trees.tsv")}.get(fam)
    inp = {r["family"]: r for r in read_tsv(OUT_DIR / "tier1_inputs.tsv")}.get(fam)
    if not old or not inp:
        return {}
    og = {l for l in inp["outgroup"].split(",") if l}
    part = root_partition(parse((TREE_DIR / f"{fam}.treefile").read_text()), og)
    return {"s7b_outgroup": inp["outgroup_family"],
            "s7b_one_clade": old["outgroup_monophyletic"],
            "s7b_root_ufboot": old["root_ufboot"],
            "s7b_split": None if part is None else _restrict(part, shared)}


def criteria(trees: list[tuple[dict, list | None]]) -> tuple:
    """D47's three tests over a family's two trees → (a, b, c, jaccard, reason)."""
    a = all(t["outgroup_monophyletic"] for t, _ in trees)
    b = a and all(t["root_ufboot"] != "" and float(t["root_ufboot"]) >= MIN_ROOT_UFBOOT
                  for t, _ in trees)
    parts = [p for _, p in trees]
    c = a and len(parts) == 2 and {s for s, _ in parts[0]} == {s for s, _ in parts[1]}
    jac = ""
    if a and len(parts) == 2:
        x, y = (min(p, key=lambda q: len(q[0]))[0] for p in parts)
        jac = round(len(x & y) / len(x | y), 3)
    reason = ("" if a and b and c else
              "outgroup not one clade" if not a else
              "root edge UFBoot < 95" if not b else "root moves with the outgroup")
    return a, b, c, jac, reason


def family_row(fam: str, trees: list[tuple[dict, list | None]]) -> dict:
    a, b, c, jac, reason = criteria(trees)
    parts = [p for _, p in trees]
    n_in = int(trees[0][0]["n_ingroup"])
    inputs = {r["root_set"]: r for r in read_tsv(OUT_DIR / "tier1_reroot_inputs.tsv")
              if r["family"] == fam}
    d39 = set(_labels(fam))
    old = s7b_root(fam, d39)
    new_split = _restrict(parts[0], d39) if c else None
    return {"family": fam, "root_sets": ",".join(t["root_set"] for t, _ in trees),
            "n_ingroup": n_in,
            "n_basal": next(iter(inputs.values()))["n_basal"],
            "a_one_clade": a, "b_root_ufboot_ge95": b, "c_same_root": c,
            "smaller_clade_jaccard": jac,
            "verdict": "resolved" if a and b and c else "unresolved",
            "reason": reason,
            "root_split_sizes": trees[0][0]["root_clade_sizes"] if c else "",
            "s7b_outgroup": old.get("s7b_outgroup", ""),
            "s7b_one_clade": old.get("s7b_one_clade", ""),
            "s7b_root_ufboot": old.get("s7b_root_ufboot", ""),
            "s7b_split_on_shared": ("" if not c or old.get("s7b_split") is None
                                    else "same" if old["s7b_split"] == new_split
                                    else "different")}


def _labels(fam: str) -> list[str]:
    from s6_lib import s6_dir
    from s3_hmm_lib import iter_fasta
    return [h.split()[0] for h, _ in iter_fasta(s6_dir("family") / f"{fam}.fasta")]


def cmd_parse(_a) -> None:
    by_fam: dict[str, list] = {}
    for row in read_tsv(OUT_DIR / "tier1_reroot_inputs.tsv"):
        d = s7c_dir("iqtree", tree_id(row))
        if not (d / f"{tree_id(row)}.treefile").exists():
            print(f"{tree_id(row)}: no tree yet", flush=True)
            continue
        by_fam.setdefault(row["family"], []).append(tree_row(row))
    trees = [t for v in by_fam.values() for t, _ in v]
    write_tsv(OUT_DIR / "tier1_reroot_trees.tsv", TREE_FIELDS, trees)
    fams = [family_row(f, v) for f, v in sorted(by_fam.items()) if len(v) == 2]
    write_tsv(OUT_DIR / "tier1_reroot.tsv", FAM_FIELDS, fams)
    for r in fams:
        print(f"{r['family']:10s} {r['verdict']:10s} {r['reason']}")
    print(f"{len(trees)} trees, {len(fams)} families → {OUT_DIR / 'tier1_reroot.tsv'}")
