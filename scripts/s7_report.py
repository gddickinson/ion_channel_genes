"""s7_report.py — render results/phylogeny/tier1_report.md from S7's tables (D13).

    python3 scripts/s7_report.py
"""

from __future__ import annotations

import json
import statistics as st
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv  # noqa: E402
from s7_lib import OUT_DIR, ROOT, TRIM_METHODS  # noqa: E402

A = ROOT / "results" / "alignments" / "alignments.tsv"


def _f(x) -> float | None:
    return float(x) if x not in ("", None) else None


def trim_section(L: list[str]) -> dict:
    rows = read_tsv(OUT_DIR / "tier1_trim_compare.tsv")
    L += ["## 1. Trimming, decided before any tree (D41)", "",
          "S6 trimmed with trimAl `-automated1`, which on low-identity families "
          "switches to its strict mode and removes the variable-but-aligned "
          "columns that carry the within-family signal. Four candidates were "
          "measured on all 63 family alignments on alignment properties alone "
          "(`tier1_trim_compare.tsv`); **no tree was built or read to choose.**", "",
          "| method | median columns | median informative sites | median gap fraction "
          "| median share of a member's residues kept | median of the per-family minimum |",
          "|---|---|---|---|---|---|"]
    med = {}
    for m in TRIM_METHODS:
        r = [x for x in rows if x["method"] == m]
        g = {k: st.median(float(x[k]) for x in r)
             for k in ("cols", "informative", "gap_frac", "retention_median",
                       "retention_min")}
        med[m] = g
        L.append(f"| {m} | {g['cols']:.0f} | {g['informative']:.0f} | "
                 f"{g['gap_frac']:.3f} | {g['retention_median']:.2f} | "
                 f"{g['retention_min']:.2f} |")
    a = {x["family"]: x for x in rows if x["method"] == "automated1"}
    g = {x["family"]: x for x in rows if x["method"] == "gt0.5"}
    ge = sum(int(g[k]["informative"]) >= int(a[k]["informative"]) for k in a)
    worst = min(a, key=lambda k: int(g[k]["informative"]) / max(1, int(a[k]["informative"])))
    L += ["", f"**Chosen: `-gt 0.5`** — keep every column in which at least half the "
          f"family has a residue, for every family. It keeps a median "
          f"{med['gt0.5']['informative']:.0f} informative sites against "
          f"{med['automated1']['informative']:.0f} for automated1, at least as many in "
          f"**{ge}/{len(a)}** families (the exception, {worst}: "
          f"{g[worst]['informative']} vs {a[worst]['informative']}), and keeps a median "
          f"{med['gt0.5']['retention_median']:.0%} of each member's own residues against "
          f"{med['automated1']['retention_median']:.0%}, at a median gap fraction of "
          f"{med['gt0.5']['gap_frac']:.1%}. The three families S6 flagged: "
          + ", ".join(f"{k} {a[k]['cols']} → {g[k]['cols']} columns"
                      for k in ("k2p", "nachr", "cng")) + ".", ""]
    return med


def rooting_section(L: list[str], inputs: list[dict]) -> None:
    c = Counter(r["root_rule"] for r in inputs)
    L += ["## 2. Rooting — from the catalogue only (D41)", "",
          "A family is rooted on the catalogue exemplars of its superfamily's declared "
          "outgroup family (`Superfamily.root_with`), added to its S6 alignment with "
          "`mafft --localpair --maxiterate 1000 --add --keeplength` (every ingroup row "
          "checked unchanged) and trimmed by the same D41 column mask. Nothing is "
          "midpoint-rooted.", "",
          "| rooting | families |", "|---|---|"]
    for k, v in sorted(c.items()):
        fams = ", ".join(r["family"] for r in inputs if r["root_rule"] == k)
        L.append(f"| {k} ({v}) | {fams} |")
    too_few = [r["family"] for r in read_tsv(A) if r["status"] != "aligned"]
    L += ["", f"No tree at all (< 4 sequences in the D39 set): {', '.join(too_few)}.", "",
          "**Outgroup occupancy** is the share of the family's trimmed columns an "
          "outgroup residue fills. The P-loop families are rooted on KcsA, MthK and "
          "NaK — two-helix pores — so in a four-repeat or TRP alignment the outgroup "
          "is mostly gap:", "",
          "| family | outgroup | occupancy |", "|---|---|---|"]
    for r in sorted(inputs, key=lambda r: r["family"]):
        if r["outgroup"]:
            L.append(f"| {r['family']} | {r['outgroup'].replace('OG_', '')} | "
                     f"{r['og_occupancy']} |")
    L.append("")


def trees_section(L: list[str], trees: list[dict], inputs: list[dict]) -> dict:
    L += ["## 3. Trees", "",
          "IQ-TREE 2 (`-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -bnni -seed 1`) on every "
          "family; the model is ModelFinder's BIC choice among four general empirical "
          "matrices with their rate and frequency variants. Support is summarised over the family's own "
          "internal edges (the outgroup excluded).", ""]
    if not trees:
        L += ["*No trees yet.*", ""]
        return {}
    rooted = [t for t in trees if t["outgroup_family"]]
    mono = [t for t in rooted if t["outgroup_monophyletic"] == "True"]
    strong = [t for t in mono if (_f(t["root_ufboot"]) or 0) >= 95]
    fr = [_f(t["frac_ge95"]) for t in trees if _f(t["frac_ge95"]) is not None]
    secs = sum(float(t["seconds"]) for t in trees)
    models = Counter(t["model"].split("+")[0] for t in trees)
    s = {"trees": len(trees), "rooted": len(rooted), "outgroup_clade": len(mono),
         "root_ge95": len(strong), "frac_ge95_median": st.median(fr) if fr else None,
         "cpu_hours": secs / 3600}
    L += [f"**{len(trees)} trees.** Of the {len(rooted)} with an outgroup, the outgroup "
          f"forms one clade — the family is monophyletic with respect to it and the root "
          f"is defined — in **{len(mono)}**, with UFBoot ≥ 95 on the root edge in "
          f"**{len(strong)}**. Median share of a family's internal edges at UFBoot ≥ 95: "
          f"**{s['frac_ge95_median']:.2f}**. Substitution matrices chosen: "
          + ", ".join(f"{k} {v}" for k, v in models.most_common()) +
          f". Wall time summed over families: {s['cpu_hours']:.1f} h.", ""]
    notclade = [t for t in rooted if t["outgroup_monophyletic"] != "True"]
    if notclade:
        L += ["**Outgroup not one clade** (root undefined — reported, not repaired): "
              + ", ".join(f"{t['family']} (outgroup {t['outgroup_family']})"
                          for t in notclade) + ".", ""]
    L += ["| family | seqs | cols | model | rooting | root UFBoot | internal edges | "
          "median UFBoot | ≥ 95 | < 70 |", "|---|---|---|---|---|---|---|---|---|---|"]
    for t in sorted(trees, key=lambda t: -int(t["n_ingroup"])):
        if t["outgroup_family"]:
            root = (f"{t['root_ufboot'] or '—'}" if t["outgroup_monophyletic"] == "True"
                    else "not a clade")
            rule = f"on {t['outgroup_family']}"
        else:
            root, rule = "—", t["root_rule"].split(":")[1].replace("_", " ")
        L.append(f"| {t['family']} | {t['n_ingroup']} | {t['cols']} | {t['model']} | "
                 f"{rule} | {root} | {t['internal_edges']} | {t['ufboot_median']} | "
                 f"{t['frac_ge95']} | {t['frac_lt70']} |")
    missing = sorted({r["family"] for r in inputs} - {t["family"] for t in trees})
    if missing:
        L += ["", f"Not yet built: {', '.join(missing)}."]
    L.append("")
    return s


def main() -> int:
    inputs = read_tsv(OUT_DIR / "tier1_inputs.tsv")
    tp = OUT_DIR / "tier1_trees.tsv"
    trees = read_tsv(tp) if tp.exists() else []
    L = ["# S7 — tier-1 phylogenies", "",
         "Rendered by `scripts/s7_report.py` from `tier1_trim_compare.tsv`, "
         "`tier1_inputs.tsv` and `tier1_trees.tsv` (D13). Treefiles: `tier1/`. Bulk "
         "(alignments with outgroups, IQ-TREE output): `<data root>/trees/s7/`.", "",
         f"![S7](figures/{'tier1_trees' if trees else 'tier1_trim'}.png)", ""]
    med = trim_section(L)
    rooting_section(L, inputs)
    s = trees_section(L, trees, inputs)
    L += ["## 4. What these trees are and are not", "",
          "* A **protein tree** of one family (tier 1, full length after D41 trimming). "
          "Cross-family structure is S8's pore-module trees; cross-superfamily "
          "relationships are the fold network, never a tree (D27).",
          "* Sets are S6's D39 sets: high-confidence profile calls and intact genome "
          "loci only. A family's tree therefore holds what the panel's 52 species carry, "
          "not every sequence in UniProt.",
          "* Where the outgroup is a two-helix prokaryotic pore and the family is a "
          "24-helix chain, the root rests on the pore columns alone; the occupancy table "
          "above says how few. Root support is reported per tree, never assumed.", ""]
    (OUT_DIR / "tier1_report.md").write_text("\n".join(L))
    (OUT_DIR / "tier1_summary.json").write_text(json.dumps(
        {"trim": med, "rooting": dict(Counter(r["root_rule"] for r in inputs)),
         "trees": s}, indent=1))
    print(f"→ {OUT_DIR / 'tier1_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
