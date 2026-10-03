"""s8_report_trees.py — § 6 of `results/phylogeny/tier2_report.md` (S8b),
rendered from `tier2_trees.tsv`, `tier2_families.tsv` and
`tier2_placement.tsv` only (D13). Imported by `s8_report.py`.
"""

from __future__ import annotations

from s3_hmm_lib import read_tsv


def _table(rows: list[dict], cols: list[str]) -> str:
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def tested_groups(fams: list[dict]) -> list[dict]:
    """The rows a clade test is read on: per module for families with several
    modules (a four-repeat chain's repeats are not expected to group), the whole
    family otherwise."""
    multi = {(r["unit"], r["family"]) for r in fams if r["module"]}
    return [r for r in fams if bool(r["module"]) == ((r["unit"], r["family"]) in multi)]


def section(d) -> list[str]:
    trees = read_tsv(d / "tier2_trees.tsv")
    fams = read_tsv(d / "tier2_families.tsv")
    place = read_tsv(d / "tier2_placement.tsv")
    tg = tested_groups(fams)
    rooted = [t for t in trees if t["rooted"] == "True"]
    split_og = [t for t in trees if t["outgroup_monophyletic"] == "False"]
    one = sum(r["one_clade"] == "True" for r in tg)
    L = ["## 6. Tier-2 trees (S8b)", "",
         f"**{len(trees)} trees, 0 failures** (IQ-TREE 2, S7's `IQ_ARGS`; "
         f"{sum(float(t['seconds']) for t in trees) / 3600:.1f} run-hours summed, "
         f"Cys-loop the longest at {max(float(t['seconds']) for t in trees) / 3600:.1f} h). "
         f"Read with `scripts/s8_tier2.py parse` (`s8_parse.py`) under the D48 rules, "
         f"none revisited after a tree was seen. Treefiles: `results/phylogeny/tier2/`.", "",
         "### 6.1 Roots", "",
         f"**{len(rooted)} of {len(trees)} trees are rooted** on their declared outgroup "
         f"({', '.join(t['unit'] for t in rooted)}); "
         f"**{len(split_og)} have an outgroup the ingroup splits** "
         f"({', '.join(t['unit'] for t in split_og)}) and are read unrooted, never "
         f"re-rooted on a substitute; the innexin clan declares none. A single-tip "
         f"outgroup (iGluR on GluR0) is a terminal edge with no support value.", "",
         _table(trees, ["unit", "n_tips", "informative", "model", "root_family",
                        "n_outgroup", "outgroup_monophyletic", "root_ufboot",
                        "root_clade_sizes", "root_clade_ufboot", "root_small_clade"]), "",
         "Where a root exists, the P2X and Ca²⁺-release trees' basal split isolates "
         "a single tip (one metazoan P2X against 24; one RyR module against four), "
         "the pattern S7d recorded at tier 1; the iGluR tree's separates seven of "
         "the eight non-vertebrate iGluR modules from the rest (UFBoot 96 / 82).", "",
         "### 6.2 Support", "",
         "Share of ingroup internal edges at UFBoot ≥ 95 and < 70:", "",
         _table(trees, ["unit", "internal_edges", "ufboot_median", "frac_ge95",
                        "frac_lt70"]), "",
         "The P-loop pore-module tree is the weakest by construction — 319 tips on 90 "
         "informative columns (the D48 cap of 4 tips per site, met at the 0.3 "
         "threshold): its deep branching order is not resolved, and nothing in "
         "§ 6.3 should be read as a statement about the order of the P-loop "
         "families' divergence.", "",
         "### 6.3 Families as clades", "",
         f"Each family's tips tested as one clade — per module for families with "
         f"several modules (four-repeat chains, K2P, TPC), whose repeats are not "
         f"expected to group. Rooted sense where the tree is rooted, unrooted (one "
         f"side of any edge) otherwise. **{one} of {len(tg)} groups are one clade.** "
         f"`n_intruders` is the fewest other tips that would have to join a group to "
         f"make it one — a descriptive near-miss measure, never a pass.", "",
         _table(tg, ["unit", "family", "module", "n_tips", "sense", "one_clade",
                     "ufboot", "n_intruders", "intruder_families"]), "",
         "Reading, from the table: every four-repeat chain's repeats group by repeat "
         "within Nav and NALCN (each repeat one clade at UFBoot ≥ 95), while Cav's "
         "repeats are scattered across the four-repeat and TPC tips — at 90 columns "
         "that is a resolution limit, not a finding (S11 owns the repeat order). "
         "Near misses are mostly a related family nested inside: the silent Kv "
         "modifier inside Shaker, ENaC inside the invertebrate degenerins, HCN "
         "with one CNG tip, 5-HT3 and ZAC inside nAChR, and in the iGluR module "
         "tree the vertebrate non-NMDA families beside NMDA tips.", "",
         "### 6.4 The four animal `plgic_prok` tips (S7a row)", "",
         "D48's kingdom rule kept the four animal members the prokaryotic-pLGIC "
         "profile calls (*Branchiostoma*, *Aplysia*, *Lottia* ×2) out of the "
         "Cys-loop outgroup as ingroup tips. Each placed tip's nested clades, "
         "smallest first:", "",
         _table(place, ["tip", "level", "clade_size", "ufboot", "composition"]), "",
         ]
    pp = next((r for r in fams if r["unit"] == "cysloop" and r["family"] == "plgic_prok"),
              None)
    if pp:
        L += [f"**They are nested in no Cys-loop family.** The four form one clade, "
              f"which joins ELIC (UFBoot 69) and then GLIC — the whole "
              f"`plgic_prok` set ({pp['n_tips']} tips) is one side of an edge at "
              f"UFBoot {pp['ufboot']}. Because they sit between GLIC and ELIC, the "
              f"declared two-tip prokaryotic outgroup is split and the Cys-loop tree "
              f"stays unrooted under D48 (not re-rooted on the six-tip clade, which "
              f"would be choosing a root after reading the tree). Either a "
              f"bacterial-type pLGIC lineage in these animals (horizontal transfer, "
              f"or retained from an ancient lineage) or divergent sequences drawn "
              f"to the long prokaryotic branches; the tree cannot tell the two "
              f"apart — S21 (genomic context, contamination check) carries the "
              f"question. (The 32-tip row is the side of an edge in an unrooted "
              f"tree — it also holds the anion-selective receptors — not a sister "
              f"relation.)", ""]
    L += ["Figure: `results/phylogeny/figures/tier2_trees.png`.", ""]
    return L
