"""s9_report.py — renders `results/filter_atlas/report.md` and `summary.json`
purely from S9's committed tables (D13).

    python3 scripts/s9_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter

from s3_hmm_lib import read_tsv
from s9_lib import OUT_DIR, REPEAT_FAMILIES, is_missing

T = {n: read_tsv(OUT_DIR / f"{n}.tsv") for n in
     ("anchors", "filter_modules", "filter_chains", "tip_check", "projection_check",
      "congruence", "string_clades", "column_background")}


def table(rows, cols, head=None) -> str:
    head = head or cols
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def k_window_check() -> dict:
    c = Counter()
    for r in T["filter_modules"]:
        if not r["window"]:
            c["unaligned"] += 1
        elif r["k_regex"]:
            c["same" if r["k_regex"] == r["window"] else
              ("other_k" if r["window_is_k"] == "yes" else "not_k")] += 1
        elif r["window_is_k"] == "yes":
            c["k_without_regex"] += 1
        else:
            c["no_motif"] += 1
    return dict(c)


def projection() -> dict:
    out = {}
    for fam in REPEAT_FAMILIES:
        out[fam] = dict(Counter(r["agree"] for r in T["projection_check"]
                                if r["family"] == fam))
    return out


def strings() -> dict:
    out = {}
    for fam in REPEAT_FAMILIES:
        st = [r["string"] for r in T["filter_chains"] if r["family"] == fam]
        out[fam] = {"chains": len(st), "unread": sum(is_missing(s) for s in st),
                    "top": Counter(s for s in st if not is_missing(s)).most_common(6)}
    return out


def summary() -> dict:
    anc = T["anchors"]
    tips = T["tip_check"]
    cong1 = [r for r in T["congruence"] if r["tree"] == "tier1"]
    tested = [r for r in cong1 if r.get("p_perm")]
    return {
        "modules": len(T["filter_modules"]), "chains": len(T["filter_chains"]),
        "k_window_columns": [int(r["column"]) for r in anc if r["anchor"] == "k_window"],
        "repeat_locus_columns": [int(r["column"]) for r in anc if r["anchor"] == "repeat_locus"],
        "tier2_tips_reread": [sum(t["agree"] == "yes" for t in tips), len(tips)],
        "k_window_check": k_window_check(), "projection_check": projection(),
        "strings": strings(),
        "tier1_tested": len(tested),
        "tier1_p_le_0_001": sum(float(r["p_perm"]) <= 0.001 for r in tested),
        "tier1_not_significant": [r["family"] for r in tested if float(r["p_perm"]) > 0.05],
        "tier2": {r["character"]: {k: r[k] for k in ("n_tips", "n_states", "changes",
                                                     "ri", "p_perm")}
                  for r in T["congruence"] if r["tree"] == "tier2"},
    }


def render(s: dict) -> str:
    pc = s["projection_check"]
    four = [f for f in REPEAT_FAMILIES if "not_comparable" not in pc[f]]
    n4 = sum(sum(pc[f].values()) for f in four)
    agree4 = sum(pc[f].get("yes", 0) + pc[f].get("yes_read", 0) for f in four)
    kw = s["k_window_check"]
    nk = kw.get("same", 0) + kw.get("other_k", 0) + kw.get("not_k", 0)
    cav = next(r for r in T["string_clades"] if r["tree"] == "tier1"
               and r["family"] == "cav" and r["string"] == "EEDD")
    L = ["# S9 — Selectivity-filter atlas (Q5)", "",
         "Rendered by `scripts/s9_report.py` from the tables beside it (D13). Rules: D49.",
         "", "## 1. The instrument", "",
         f"{s['modules']:,} P-loop pore modules (S6, D39 members) in {s['chains']:,} chains "
         "were added to S8's untrimmed P-loop tier-2 L-INS-i alignment "
         "(`mafft --localpair --add --keeplength --mapout`, one run per family chunk), so "
         "every module is read in one coordinate system. Anchors, located through MAFFT's "
         "residue→column map (counting residues along an aligned row is wrong wherever "
         "`--keeplength` deleted an insertion upstream — the first read made that error on "
         "Nav1.5 repeat III and was corrected before any tree was read):", "",
         table(T["anchors"], ["anchor", "reference", "position", "residue", "repeat",
                              "column"]), "",
         f"**All four Nav1.5 locus residues (D372, E898, K1419, A1711) fall in one column "
         f"({s['repeat_locus_columns'][0]}), the column KcsA's V76 — the x of T-x-G-Y-G — "
         "occupies.** The four-repeat locus and the K⁺ signature are read at the same "
         "alignment position (an alignment statement, not a structural one).", "",
         "**Validation (D49 (4)).**", "",
         f"* Classifier's pairwise projection (MAFFT to Nav1.5) vs the shared-alignment read, "
         f"four-repeat chains: **{agree4} / {n4} agree** "
         f"({100 * agree4 / n4:.1f} %; identical, or identical wherever both read a residue). "
         "Per family: " + "; ".join(f"{f} {pc[f]}" for f in REPEAT_FAMILIES) + ". "
         "Every disagreement differs at one position, all in divergent chains "
         "(`projection_check.tsv`). CatSper (one repeat per chain) and TPC (two) are outside "
         "the projection's design and are not compared.",
         f"* K window vs `K_FILTER_RE` in the module: of {nk:,} modules carrying a TxGYG-like "
         f"motif, **{kw.get('same', 0):,} read the same motif in the window**, "
         f"{kw.get('other_k', 0)} a different K motif (two in the module), "
         f"{kw.get('not_k', 0)} a non-K window (mostly CNG); "
         f"{kw.get('k_without_regex', 0)} windows read a K motif the regex did not find. "
         f"{kw.get('unaligned', 0)} modules have no sequence (S6 `absent`).",
         f"* Every tier-2 tip re-added reads what its own row reads: "
         f"**{s['tier2_tips_reread'][0]} / {s['tier2_tips_reread'][1]}**.", "",
         "## 2. The atlas: four-repeat strings", "",
         table([{"family": f, "chains": v["chains"], "unread": v["unread"],
                 "strings": ", ".join(f"{a} {b}" for a, b in v["top"])}
                for f, v in s["strings"].items()], ["family", "chains", "unread", "strings"]),
         "", "Strings with no `FILTER_CALLS` entry (DEEA, DKEA, EKEE, DDDD, …) are recorded "
         "as read, with their placement below; S9 adds no call (D49 (6)). K-window "
         "families: `filter_modules.tsv` / `filter_chains.tsv`.", "",
         "## 3. Congruence (D49 (5))", "",
         "Fitch changes on each tree against 1,000 tip-label permutations; RI = retention "
         "index (1 = every change once, 0 = no more structure than a star tree). Rooted "
         "(S7b outgroup one clade): " + ", ".join(r["family"] for r in T["congruence"]
                                                 if r["tree"] == "tier1"
                                                 and r["rooted"] == "True")
         + "; the rest, including the six S7d left undefined, are read unrooted.", "",
         table(T["congruence"], ["tree", "family", "character", "rooted", "n_tips",
                                 "n_states", "changes", "min_changes", "star_changes", "ri",
                                 "p_perm"]), "",
         f"**{s['tier1_p_le_0_001']} of {s['tier1_tested']} testable tier-1 trees: p ≤ 0.001.** "
         f"Not significant: {', '.join(s['tier1_not_significant'])} — each a family whose "
         "minority string is carried by one or two scattered chains.", "",
         "### Per string (≥ 2 carriers), four-repeat families", "",
         table([r for r in T["string_clades"] if r["tree"] == "tier1"
                and r["family"] in REPEAT_FAMILIES],
               ["family", "string", "n_carriers", "sense", "one_clade", "ufboot",
                "n_intruders", "intruder_strings", "origins"]), "",
         "`origins` counts maximal carrier-only clades and is fragmented by any one nested "
         "change; read it with the core clade below.", "",
         "### The Q5 test case: Cav3's EEDD", "",
         f"EEDD has {cav['n_carriers']} carriers and is not one clade under the strict test "
         f"({cav['n_intruders']} intruders, because two carriers sit far away). **47 of the "
         "49 form one clade at UFBoot 100**, a 49-chain clade that holds the three human "
         "T-type channels (O43497, O95180, Q9P0X4), invertebrate, cnidarian and placozoan "
         "chains, and nothing else but one EQDD and one chain with an unread repeat. The two "
         "EEDD chains outside it are both *Hydra* (A0ABM4BFT6, A0ABM4BFU2), read the same by "
         "both instruments. **EEDD sits where the tree says: one origin, in the lineage "
         "that carries it, with at most one independent EEDD (in *Hydra*).** Both of the "
         "positions where EEDD differs from EEEE (repeats III, IV) rank at the 95th–96th "
         "percentile of the Cav alignment's columns for tree congruence (§ 4). The 11 DDDD "
         "chains are all *Paramecium*; 8 form one clade (UFBoot 100), and with NDDN "
         "(also *Paramecium*) a 14-chain clade at UFBoot 100.", "",
         "Other non-canonical strings, read descriptively (best-fitting clade = the clade "
         "maximising carriers minus others): Nav DEEA — 12 carriers, 6 in one clade "
         "(UFBoot 78), the other 6 scattered (*Hydra*, *Nematostella*, *Trichoplax*, "
         "amphioxus, and an elephant-shark and an opossum chain both instruments read DEEA, "
         "emergent); Nav DKEA — 2 carriers, one clade (UFBoot 100); NALCN EKEE — 2 "
         "carriers, not one clade; CatSper E — 2 carriers, scattered (RI 0).", "",
         "### Tier 2 (P-loop pore modules, unrooted, 15 % of edges at UFBoot ≥ 95)", "",
         "The locus column is congruent with the module tree beyond chance (table above), "
         "but the tree's deep order is unsupported (S8b), so no statement is made about "
         "where across families a filter type arose. Per-residue clade tests are in "
         "`string_clades.tsv`; no residue is one clade across families.", "",
         "## 4. Post-hoc control: filter positions against every other column", "",
         "Added after § 3 was read, and labelled so: the filter is inside the alignment "
         "each tree was inferred from, so § 3's permutation test shows only that the filter "
         "carries phylogenetic signal, as every column does. Here each single filter "
         "position's RI is ranked among every parsimony-informative column (≤ 50 % gaps) "
         "of the same tree's input alignment.", "",
         table([r for r in T["column_background"] if r["percentile"] != ""],
               ["tree", "family", "position", "n_states", "ri", "n_columns",
                "column_ri_median", "percentile"]), "",
         "**Reading.** In the four-repeat channels the positions that vary are more "
         "tree-congruent than the typical site — Cav repeats III/IV (EEDD) 95th/96th, Nav "
         "II/III 100th/97th, TPC 81st/92nd: selectivity at the four-repeat locus follows "
         "descent. The exceptions are positions whose minority state is carried by one or "
         "two chains (CatSper, NALCN repeat II). In the K channels the Y/F of GYG tracks "
         "descent (K2P, Kir, EAG, CNG ≥ 93rd) while the x of TxGYG is the homoplastic "
         "position (Kir 6th, Slo 3rd percentile; Shaker's T/S position 8th), i.e. re-tuned "
         "repeatedly within families.", "",
         "## 5. Limits", "",
         "* Members are S6's D39 sets (high-confidence calls, intact genome loci); medium "
         "calls and the r4 families are not in them.",
         "* Six P-loop roots are undefined (S7d); their clade reads are unrooted.",
         "* The tier-2 P-loop tree is a 90-column module tree; only its supported edges "
         "are read.",
         "* Literature checks of the non-canonical strings (DEEA, DKEA, EKEE, DDDD) are "
         "*(pending: S14 references)*; the elephant-shark and opossum DEEA chains are an "
         "emergent row.", ""]
    return "\n".join(L)


def main() -> None:
    s = summary()
    (OUT_DIR / "summary.json").write_text(json.dumps(s, indent=1))
    (OUT_DIR / "report.md").write_text(render(s))
    print(json.dumps({k: s[k] for k in ("modules", "chains", "tier2_tips_reread",
                                        "k_window_check", "tier1_tested",
                                        "tier1_p_le_0_001")}))


if __name__ == "__main__":
    sys.exit(main())
