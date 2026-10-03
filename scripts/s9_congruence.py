"""s9_congruence.py — S9: do selectivity filters track phylogeny? (Q5, D49)

    python3 scripts/s9_congruence.py

Lays S9's filter reads (`results/filter_atlas/filter_chains.tsv`,
`filter_modules.tsv`) on

* every P-loop **tier-1** tree (S7b): one tip per chain, the chain's string
  (four-repeat families: the locus residues of repeats I…n; the rest: the K
  windows of its modules) — rooted on the outgroup where S7b found it one
  clade, unrooted otherwise and for the six families whose root S7d left
  undefined (D47);
* the P-loop **tier-2** tree (S8b, unrooted): one tip per module, two
  characters — the residue at the shared repeat-locus column, and the K
  window.

Per tree × character: Fitch changes, the minimum (states − 1) and star-tree
length, retention index, permutation p (`s9_lib.congruence`). Per string
with ≥ 2 carriers: one clade or not, UFBoot, intruders (S8b's
`clade_test`), and origins (maximal carrier-only clades, rooted trees only).

Writes `results/filter_atlas/congruence.tsv`, `string_clades.tsv`.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict

from s3_hmm_lib import read_tsv, write_tsv
from s7_newick import parse
from s7c_basal import reroot_families
from s8_parse import clade_test, clades
from s9_lib import OUT_DIR, PHY_DIR, UNIT, congruence, is_missing, origins, prune

CONG_FIELDS = ["tree", "family", "character", "rooted", "n_tips", "n_missing", "n_states",
               "changes", "min_changes", "star_changes", "ri", "null_median", "null_min",
               "p_perm"]
CLADE_FIELDS = ["tree", "family", "character", "string", "n_carriers", "sense", "one_clade",
                "ufboot", "n_intruders", "intruder_strings", "origins", "carriers"]
MAX_LISTED = 12      # carriers listed per row (all are counted)


def tier1_trees() -> list[dict]:
    rows = {r["family"]: r for r in read_tsv(PHY_DIR / "tier1_trees.tsv")
            if r["superfamily"] == UNIT}
    inputs = {r["family"]: r for r in read_tsv(PHY_DIR / "tier1_inputs.tsv")}
    reroot = set(reroot_families())
    out = []
    for fam, r in sorted(rows.items()):
        og = [x for x in inputs[fam]["outgroup"].split(",") if x]
        # S7d left the re-rooted families' roots undefined (D47): read unrooted
        rooted = (r["root_rule"] == "rooted" and r["outgroup_monophyletic"] == "True"
                  and fam not in reroot)
        out.append({"family": fam, "path": PHY_DIR / "tier1" / f"{fam}.treefile",
                    "outgroup": frozenset(og), "rooted": rooted})
    return out


def string_rows(tree_name, fam, char, tree, states, outgroup, rooted) -> list[dict]:
    known = {t: s for t, s in states.items() if not is_missing(s)}
    keep = set(known) | (set(outgroup) if rooted else set())
    sub = prune(tree, keep)
    exclude = frozenset(outgroup) if rooted else frozenset()
    cl = clades(sub, exclude)
    by = defaultdict(set)
    for t, s in known.items():
        by[s].add(t)
    out = []
    for s, carriers in sorted(by.items(), key=lambda x: (-len(x[1]), x[0])):
        if len(carriers) < 2:
            continue
        one, sup, intr = clade_test(sub, carriers, exclude, cl)
        icount = Counter(known.get(t, "?") for t in intr)
        out.append({"tree": tree_name, "family": fam, "character": char, "string": s,
                    "n_carriers": len(carriers), "sense": "rooted" if rooted else "unrooted",
                    "one_clade": one, "ufboot": "" if sup is None else f"{sup:g}",
                    "n_intruders": len(intr),
                    "intruder_strings": ",".join(f"{k}:{v}" for k, v in icount.most_common(6)),
                    "origins": origins(sub, carriers, exclude) if rooted else "",
                    "carriers": ",".join(sorted(carriers)[:MAX_LISTED])})
    return out


def cong_row(tree_name, fam, char, tree, states, rooted, n_all) -> dict:
    c = congruence(tree, states)
    row = {"tree": tree_name, "family": fam, "character": char, "rooted": rooted,
           "n_missing": n_all - c["n_tips"]}
    row.update({k: ("" if c.get(k) is None else c[k]) for k in CONG_FIELDS if k in c})
    return row


def tier1(chains) -> tuple[list[dict], list[dict]]:
    cong, strings = [], []
    for t in tier1_trees():
        fam = t["family"]
        tree = parse(t["path"].read_text())
        ingroup = set(tree.leaves()) - t["outgroup"]
        states = {lab: chains.get((fam, lab)) for lab in ingroup}
        char = next((r["character"] for (f, _), r in chains_rows.items() if f == fam), "")
        ing = prune(tree, ingroup)
        cong.append(cong_row("tier1", fam, char, ing, states, t["rooted"], len(ingroup)))
        strings += string_rows("tier1", fam, char, tree, states, t["outgroup"], t["rooted"])
        print(f"tier1 {fam}: {cong[-1].get('changes')} changes, RI {cong[-1].get('ri')}, "
              f"p {cong[-1].get('p_perm')}", flush=True)
    return cong, strings


def tier2(modules) -> tuple[list[dict], list[dict]]:
    tree = parse((PHY_DIR / "tier2" / f"{UNIT}.treefile").read_text())
    tips = tree.leaves()
    cong, strings = [], []
    for char in ("locus", "window"):
        states = {}
        for tip in tips:
            fam, rest = tip.split("__", 1)
            r = modules.get((fam, rest))
            states[tip] = r[char] if r else None
        cong.append(cong_row("tier2", UNIT, char, tree, states, False, len(tips)))
        strings += string_rows("tier2", UNIT, char, tree, states, frozenset(), False)
        print(f"tier2 {char}: {cong[-1]}", flush=True)
    return cong, strings


chains_rows: dict = {}


def main() -> None:
    global chains_rows
    chains_rows = {(r["family"], r["chain"]): r for r in read_tsv(OUT_DIR / "filter_chains.tsv")}
    chains = {k: r["string"] for k, r in chains_rows.items()}
    modules = {(r["family"], r["label"]): r for r in read_tsv(OUT_DIR / "filter_modules.tsv")}
    c1, s1 = tier1(chains)
    c2, s2 = tier2(modules)
    write_tsv(OUT_DIR / "congruence.tsv", CONG_FIELDS, c1 + c2)
    write_tsv(OUT_DIR / "string_clades.tsv", CLADE_FIELDS, s1 + s2)


if __name__ == "__main__":
    sys.exit(main())
