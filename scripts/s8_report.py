"""s8_report.py — renders `results/phylogeny/tier2_report.md` from the
committed S8 tables only (D13): tier-2 units, representatives, inputs, the
ITPR span check, the tier-3 fold network, and (S8b, `s8_report_trees.py`)
the tier-2 trees once `tier2_trees.tsv` exists.

    python3 scripts/s8_report.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv  # noqa: E402

D = Path(__file__).resolve().parents[1] / "results" / "phylogeny"
NET = D / "fold_network"


def table(rows: list[dict], cols: list[str]) -> str:
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def main() -> int:
    units = read_tsv(D / "tier2_units.tsv")
    thr = read_tsv(D / "tier2_threshold.tsv")
    inp = read_tsv(D / "tier2_inputs.tsv")
    itpr = read_tsv(D / "tier2_itpr_span.tsv")
    struct = read_tsv(NET / "structures.tsv")
    nodes = read_tsv(NET / "units.tsv")
    edges = read_tsv(NET / "literature_edges.tsv")
    med = read_tsv(NET / "superfamily_medians.tsv")
    pairs = read_tsv(NET / "pairs.tsv")
    trees = D / "tier2_trees.tsv"

    treed = [u for u in units if u["tier2"] == "tree"]
    same = [u for u in units if u["tier2"] == "= tier 1"]
    refused = [u for u in units if u["tier2"].startswith("refused")]
    lit = {tuple(sorted((e["a"], e["b"]))) for e in edges}
    unasserted = sorted((m for m in med if m["a"] != m["b"]
                         and (m["a"], m["b"]) not in lit),
                        key=lambda m: -float(m["median_tm"]))[:8]
    missing = [s for s in struct if not s["model"]]
    L = ["# S8 — tier-2 pore-module trees and the tier-3 fold network", "",
         "Rendered by `scripts/s8_report.py` from the committed tables (D13). "
         "Rules: **D48** (fixed before any tier-2 tree or structure comparison).", "",
         "## 1. What gets a tier-2 tree", "",
         f"**{len(treed)} units get a tier-2 tree** (alignable, ≥ 2 census families); "
         f"**{len(same)} single-family superfamilies** are their tier-1 tree "
         f"({sum('no tier-1' in u['note'] for u in same)} of them have none: < 4 "
         f"sequences or added after S6); **{len(refused)} are refused by D27** "
         f"(`build_tier2` raises `NotAlignable`) and go to the fold network — a result, "
         f"not an error.", "",
         table(units, ["superfamily", "tier2", "families", "kind", "note"]), "",
         "## 2. Representatives (D8, D48)", "",
         "Cell = family × panel group × module (repeat). Within a cell, members are "
         "ordered by centrality (mean within-family identity) and clustered greedily; "
         "a member joins the first representative at ≥ the threshold. One threshold "
         "(0.5) for every unit, except where it leaves more than 4 tips per "
         "informative site: the P-loop pore module is ~105 residues and trims to "
         "85 columns, so 840 tips at 0.5 became 317 at 0.3 — fixed on alignment "
         "properties before any tree. Tip counts at every candidate threshold:", "",
         table(thr, list(thr[0])), "",
         "## 3. Inputs", "",
         "MAFFT L-INS-i over each unit (rows, raggedness and residues checked), "
         "trimAl `-gt 0.5` (D41), IQ-TREE 2 with S7's `IQ_ARGS`. **Root**: the "
         "superfamily's `root_with` family; its tips are the outgroup only within "
         "the kingdom of its catalogue exemplars, and the exemplars are added as tips "
         "(module units: cut by S6's method) unless a tip already carries the "
         "sequence. The four animal `plgic_prok` members are therefore *ingroup* "
         "tips in the Cys-loop tree — where they fall is the answer to the S7a row. "
         "innexin clan: no declared outgroup, unrooted.", "",
         table(inp, ["unit", "kind", "n_tips", "n_families", "root_family", "n_outgroup",
                     "cols", "informative", "gap_frac"]), "",
         "Outgroup tips excluded by the kingdom rule: "
         + "; ".join(f"{r['unit']}: {r['outgroup_excluded']}" for r in inp
                     if r["outgroup_excluded"]) + ".", "",
         "## 4. ITPR's module span, checked against structure (S6 row)", "",
         "S6 fixed ITPR's span by a RyR-seeded vote plus the helix snap. On the three "
         "human ITPRs the extracted module starts at UniProt's TM5 and ends 2–3 "
         "residues past TM6 — the convention of every annotated family — and contains "
         "the GVGD filter; on cryo-EM 6DQJ (human ITPR3) it starts inside the TM5 "
         "helix, spans the pore helix and filter, and ends inside TM6.", "",
         table(itpr, list(itpr[0])), "",
         "## 5. The tier-3 fold network (D27, D48)", "",
         f"One AlphaFold DB model per census family (first S0 exemplar with a model of "
         f"exactly that accession), cut to its comparison unit (module 1 for module "
         f"superfamilies, whole model otherwise; pLDDT ≥ 70 only), plus a Kv voltage-"
         f"sensor node. **{len([n for n in nodes if n['residues']])} nodes, "
         f"{len(pairs)} pairs**, TM-align (average-length TM-score) and Foldseek "
         f"(exhaustive). Unmeasured: "
         + "; ".join(f"{s['node']} ({s['note']})" for s in missing) + ".", "",
         "**Verdict rule, fixed first**: an edge is `supported` if the median family-"
         "pair TM-score is ≥ 0.5 *and* each side's best other superfamily is the other "
         "side; `not_distinguished` otherwise. An edge is structural similarity, never "
         "phylogeny — `network.json` carries no branch lengths or support values.", "",
         table(edges, ["a", "b", "median_tm", "n_pairs", "rank_b_for_a", "rank_a_for_b",
                       "best_other", "best_other_tm", "verdict"]), "",
         "Reading: small helical pore modules score 0.4–0.6 against almost any helical "
         "unit of their size, which is why the rank test, not the bar alone, carries "
         "the verdict. TMEM16/OSCA/TMC misses the bar (0.484) while sitting far above "
         "its best outside partner — `not_distinguished` as written. ITPR's module "
         "(65 confident residues) is as close to the iGluR pore as to the P-loop "
         "pore: it is a P-loop-like pore, but this measurement cannot place it nearer "
         "the P-loop than the iGluR.", "",
         "**Strongest pairs no literature edge asserts** (for S18, not claims):", "",
         table(unasserted, ["a", "b", "median_tm"]), "",
         "Figure: `results/phylogeny/figures/fold_network.png`.", ""]
    if trees.exists():
        from s8_report_trees import section
        L += section(D)
    else:
        L += ["## 6. Tier-2 trees", "",
              "Running detached (`scripts/s8_tier2.py run`, log "
              "`<data root>/trees/s8/run.log`); parsed in S8b.", ""]
    (D / "tier2_report.md").write_text("\n".join(L))
    print(D / "tier2_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
