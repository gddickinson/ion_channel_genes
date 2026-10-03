"""s11_report.py — renders `results/duplication/report.md` and `summary.json`
from S11a's committed tables only (D13).

    python3 scripts/s11_report.py
"""
from __future__ import annotations

import json
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s11_lib import OUT, S7D_UNDEFINED, read_tsv     # noqa: E402

T = lambda name: read_tsv(OUT / name)                 # noqa: E731


def summary() -> dict:
    trees, dups = T("recon_trees.tsv"), T("duplications.tsv")
    pairs, genes = T("human_pairs.tsv"), T("human_genes.tsv")
    sup = [d for d in dups if d["supported"] == "True"]
    deep = [d for d in sup if d["species_specific"] == "False"]
    ss = [d for d in sup if d["species_specific"] == "True"]
    dec = [t for t in trees if t["root_mode"] == "declared"]
    rec = [t for t in trees if t["root_mode"] == "reconciliation"]
    strict = [p for p in pairs if p["ohno_strict"] == "True"]
    sw = Counter(p["window"] for p in strict)
    score = lambda w: statistics.median(                       # noqa: E731
        float(d["score"]) for p in strict if p["window"] == w
        for d in dups if d["family"] == p["family"] and d["dup"] == p["node"])
    gs = Counter((g["has_2R"], g["ohno_strict"]) for g in genes)
    return {
        "trees": len(trees), "tips": sum(int(t["n_tips"]) for t in trees),
        "declared_roots": len(dec),
        "declared_is_optimum": sum(t["declared_optimal"] == "True" for t in dec),
        "reconciliation_roots": len(rec),
        "reconciliation_unique": sum(t["n_optimal_roots"] == "1" for t in rec),
        "dups_stated": len(dups), "dups_supported": len(sup),
        "dups_unstated_root_dependent": sum(int(t["dups_unstated"]) for t in trees),
        "supported_species_specific": len(ss), "supported_deeper": len(deep),
        "ss_with_genome_tip": sum(int(d["n_genome_tips"]) > 0 for d in ss),
        "deep_by_taxon": dict(Counter(d["taxon"] for d in deep).most_common()),
        "ss_by_species": dict(Counter(d["taxon"] for d in ss).most_common()),
        "three_r_by_family": dict(Counter(d["family"] for d in deep
                                          if d["taxon"] == "Clupeocephala").most_common()),
        "two_r_by_family": dict(Counter(d["family"] for d in deep
                                        if d["taxon"] in ("Vertebrata", "Gnathostomata")).most_common()),
        "human_pairs": len(pairs), "pairs_by_window": dict(Counter(p["window"] for p in pairs)),
        "ohno_strict_pairs_in_families": len(strict),
        "ohno_strict_by_window": dict(sw),
        "ohno_strict_older_by_family": dict(Counter(p["family"] for p in strict
                                                    if p["window"] == "older").most_common()),
        "score_median_2R": score("2R"), "score_median_older": score("older"),
        "older_strict_by_taxon": dict(Counter(
            d["taxon"] for p in strict if p["window"] == "older" for d in dups
            if d["family"] == p["family"] and d["dup"] == p["node"]).most_common()),
        "genes": len(genes),
        "genes_ohno_strict": sum(g["ohno_strict"] == "True" for g in genes),
        "genes_ohno_strict_with_2R": gs[("True", "True")],
        "genes_2R": sum(g["has_2R"] == "True" for g in genes),
        "ohno_cross_family": {r["criterion"]: int(r["pairs"]) for r in T("ohnolog_compare.tsv")
                              if r["window"] == "across_families"},
        "s7d": {t["family"]: {"n_optimal_roots": int(t["n_optimal_roots"]),
                              "root_split": t["root_split"], "small_side": t["root_small_side"]}
                for t in trees if t["family"] in S7D_UNDEFINED},
    }


def table(rows, cols) -> str:
    out = ["| " + " | ".join(c for c, _ in cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(f(r)) for _, f in cols) + " |" for r in rows]
    return "\n".join(out)


def render(s: dict) -> str:
    pct = lambda a, b: f"{100 * a / b:.0f} %" if b else "–"     # noqa: E731
    sw = s["ohno_strict_by_window"]
    L = [
        "# S11a — Duplication history by reconciliation",
        "",
        "Rendered by `scripts/s11_report.py` from the tables beside it (D13). Rules: **D52**,"
        " fixed before any tree was reconciled. Driver: `scripts/s11_duplications.py`"
        " (`reconcile`, `pairs`); logic `scripts/s11_recon.py`; self-test `selftest_s11.py`.",
        "",
        "## 1. What was reconciled",
        "",
        f"S7b's **{s['trees']} tier-1 gene trees** ({s['tips']:,} ingroup tips, outgroups pruned)"
        " against the NCBI Taxonomy tree of the 52 S4 panel species (D15; polytomies kept —"
        " Chordata is a three-way star of amphioxus, *Ciona* and vertebrates, because NCBI has"
        " no Olfactores). Duplications by species overlap; losses (tie-break only) by LCA"
        " mapping. Families with no tier-1 tree (the five < 4-sequence families and the r4"
        " families, which S6 predates) are outside this task.",
        "",
        "## 2. Roots",
        "",
        f"* **{s['declared_roots']}** trees keep their S7b root (outgroup one clade, UFBoot ≥ 95"
        f" or one sequence, not an S7d family). The declared root is itself a reconciliation"
        f" optimum in **{s['declared_is_optimum']} / {s['declared_roots']}** — the outgroup root"
        " and the fewest-duplications root mostly disagree (reported, never replaced).",
        f"* **{s['reconciliation_roots']}** trees are rooted by reconciliation (no declared"
        f" outgroup, a split or weak outgroup, or S7d's six); the optimum is a single edge in"
        f" **{s['reconciliation_unique']}**. A duplication that holds under only some tied"
        f" optima is not stated: {s['dups_unstated_root_dependent']} such nodes in all.",
        "",
        "The six roots S7d left undefined, as reconciliation places them (for S11b to compare"
        " with an outgroup-free root):",
        "",
        table(sorted(s["s7d"].items()),
              [("family", lambda kv: kv[0]), ("optimal roots", lambda kv: kv[1]["n_optimal_roots"]),
               ("basal split", lambda kv: kv[1]["root_split"] or "–"),
               ("smaller side (species)", lambda kv: kv[1]["small_side"] or "tied")]),
        "",
        "## 3. Duplications",
        "",
        f"**{s['dups_stated']:,} stated duplications, {s['dups_supported']:,} supported** (each"
        f" child a leaf or UFBoot ≥ 95). Of the supported, **{s['supported_species_specific']:,}"
        f" are species-specific** ({s['ss_with_genome_tip']} involve a genome-locus tip) and"
        f" **{s['supported_deeper']} are shared by more than one species**.",
        "",
        "Supported duplications shared by more than one species, by mapped taxon:",
        "",
        table(list(s["deep_by_taxon"].items())[:15], [("taxon", lambda kv: kv[0]),
                                                      ("supported", lambda kv: kv[1])]),
        "",
        "Species-specific supported duplications, top species (a species-specific node can"
        " also be an allelic or redundant gene model; the source of every tip is in"
        " `duplications.tsv`):",
        "",
        table(list(s["ss_by_species"].items())[:12], [("species", lambda kv: f"*{kv[0]}*"),
                                                      ("supported", lambda kv: kv[1])]),
        "",
        "Teleost 3R window (mapped to Clupeocephala, *Danio* + *Takifugu*): "
        + ", ".join(f"{f} {n}" for f, n in s["three_r_by_family"].items()) + ".",
        "",
        "2R window (Vertebrata or Gnathostomata): "
        + ", ".join(f"{f} {n}" for f, n in s["two_r_by_family"].items()) + ".",
        "",
        "## 4. Human paralogue pairs and OHNOLOGS v2",
        "",
        f"{s['human_pairs']:,} human paralogue pairs inside one family, by the window of the"
        " stated duplication separating them: "
        + ", ".join(f"{w} {n}" for w, n in sorted(s["pairs_by_window"].items())) + ".",
        "",
        f"**OHNOLOGS v2 strict 2R pairs inside our families: {s['ohno_strict_pairs_in_families']}."
        f" The gene trees date {sw.get('2R', 0)} ({pct(sw.get('2R', 0), s['ohno_strict_pairs_in_families'])})"
        f" to the 2R window and {sw.get('older', 0)} ({pct(sw.get('older', 0), s['ohno_strict_pairs_in_families'])})"
        " older than vertebrates**; the rest: "
        + ", ".join(f"{w} {n}" for w, n in sw.items() if w not in ("2R", "older")) + ". "
        f"No OHNOLOGS pair joins two of our families (strict {s['ohno_cross_family'].get('strict')},"
        f" relaxed {s['ohno_cross_family'].get('relaxed')}).",
        "",
        "The too-old pairs concentrate in a few families ("
        + ", ".join(f"{f} {n}" for f, n in s["ohno_strict_older_by_family"].items())
        + ") and map to " + ", ".join(f"{t} {n}" for t, n in s["older_strict_by_taxon"].items())
        + f". Their duplications have a median species-overlap score of {s['score_median_older']}"
        f" against {s['score_median_2R']} for the pairs dated to 2R: a non-vertebrate sequence"
        " nested among vertebrate paralogues (gene-tree incongruence, or a real pre-vertebrate"
        " duplication) drags the node to a deeper taxon. D52 is not changed for it; a"
        " synteny-free tree reconciliation dates 2R reliably only where the vertebrate"
        " paralogue clades are clean.",
        "",
        f"Per gene ({s['genes']} human genes on tier-1 trees): OHNOLOGS strict lists"
        f" {s['genes_ohno_strict']} as 2R ohnologues; our trees give"
        f" {s['genes_ohno_strict_with_2R']} of them ({pct(s['genes_ohno_strict_with_2R'], s['genes_ohno_strict'])})"
        f" a 2R-window paralogue. Our trees give {s['genes_2R']} genes a 2R-window paralogue,"
        f" {s['genes_ohno_strict_with_2R']} ({pct(s['genes_ohno_strict_with_2R'], s['genes_2R'])})"
        " of them on the strict list (OHNOLOGS demands synteny with outgroup genomes).",
        "",
        "## 5. Limits",
        "",
        "* Reconciliation trusts the gene tree. Supported (UFBoot ≥ 95) child edges do not"
        " protect against a misplaced long branch: § 4 measures how often that dates a"
        " synteny-supported 2R pair too deep.",
        "* The species tree is NCBI's, with its polytomies (D15); a duplication mapped to a"
        " polytomy node is dated to that node, not resolved.",
        "* Copy number is the D39 set: proteome genes (one per gene) and intact genome loci."
        " Species-specific duplications include redundant models and genome-locus splits;"
        " they are counted apart from the headline.",
        "* The repeat-duplication order of the 4×6TM channels and outgroup-free roots are"
        " S11b's.",
        "",
        "## Files",
        "",
        "`recon_trees.tsv` (per tree: root mode, optima, declared-root agreement, counts,"
        " reconciliation basal split) · `duplications.tsv` (every stated duplication) ·"
        " `expansions.tsv` (supported per family × taxon) · `human_pairs.tsv` ·"
        " `human_genes.tsv` · `ohnolog_compare.tsv` · `summary.json` · figure"
        " `figures/duplications.png`. Archives under `<data root>/raw_api/s11/` (NCBI"
        " Taxonomy, HGNC complete set, OHNOLOGS v2 human 2R pairs).",
        "",
    ]
    return "\n".join(L)


def main() -> int:
    s = summary()
    (OUT / "summary.json").write_text(json.dumps(s, indent=1) + "\n")
    (OUT / "report.md").write_text(render(s))
    print(OUT / "report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
