"""s10_report.py — S10a's report, rendered from the committed tables (D13).

    python3 scripts/s10_report.py   → results/repertoire/report.md + summary.json
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict

from s3_hmm_lib import read_tsv
from s10_lib import BUSCO_FLOOR, GAINS, OUT, PRIMARY

NODES = ["Eukaryota", "Opisthokonta", "Fungi", "Metazoa", "Bilateria", "Ecdysozoa",
         "Lophotrochozoa", "Deuterostomia", "Vertebrata", "Mammalia", "Viridiplantae",
         "Embryophyta", "Sar", "Amoebozoa", "Discoba"]


def main() -> int:
    tree = json.loads((OUT / "species_tree.json").read_text())
    ch = read_tsv(OUT / "characters.tsv")
    fam = {(r["family"], r["cost"]): r for r in read_tsv(OUT / "families.tsv")}
    fams = sorted({f for f, _ in fam})
    nodes = defaultdict(dict)
    for r in read_tsv(OUT / "nodes.tsv"):
        nodes[(r["family"], r["node"])][r["cost"]] = r["state"]
    ev = read_tsv(OUT / "events.tsv")
    losses = read_tsv(OUT / "losses.tsv")
    place = read_tsv(OUT / "s4_species_orders.tsv")

    on_tree = {r["order"] for r in ch} - set(o for o, *_ in tree["orders_unplaced"])
    states = Counter(r["state"] for r in ch if r["order"] in on_tree)
    why = Counter("absence in a proteome below the BUSCO floor" if r["why"].startswith("no call")
                  else r["why"] for r in ch if r["order"] in on_tree and r["state"] == "?")
    leca = {k: Counter(fam[(f, k)]["leca"] for f in fams) for k in GAINS}
    leca_robust = [f for f in fams if all(fam[(f, k)]["leca"] == "present" for k in GAINS)]
    leca_cost = [f for f in fams if fam[(f, PRIMARY)]["leca"] == "present" and f not in leca_robust]
    absent_all = [f for f in fams if all(fam[(f, k)]["leca"] == "absent" for k in GAINS)]
    tot = {k: (sum(int(fam[(f, k)]["gains"]) for f in fams), sum(int(fam[(f, k)]["losses"]) for f in fams),
               sum(int(fam[(f, k)]["ambiguous_edges"]) for f in fams)) for k in GAINS}
    prim = [e for e in ev if e["cost"] == PRIMARY and e["stated"] == "1"]
    robust = Counter(e["event"] for e in prim if e["robust"] == "1")
    strength = Counter(r["strength"] for r in losses)
    node_counts = {n: Counter(nodes[(f, n)].get(PRIMARY, "?") for f in fams) for n in NODES}
    node_stable = {n: sum(1 for f in fams if nodes[(f, n)].get(PRIMARY) == "1"
                          and len(set(nodes[(f, n)].values())) == 1) for n in NODES}
    top_loss = sorted(fams, key=lambda f: -int(fam[(f, PRIMARY)]["losses"]))[:8]

    L = ["# S10a — Repertoire reconstruction", "",
         "*Rendered by `scripts/s10_report.py` from the tables beside it (D13). Rules: D50, "
         "fixed before any family was reconstructed; addenda (a)–(c) labelled.*", "",
         "## Headline", "",
         f"Every census family's presence was reconstructed on **{tree['tips']} eukaryotic orders** "
         f"(S4b, one reference proteome each) over the NCBI Taxonomy tree "
         f"({tree['internal_nodes']} internal nodes, {tree['polytomies']} polytomies, "
         f"a {len(tree['root_children'])}-way star at the eukaryotic root). "
         f"Under the primary cost (gain = 2 losses) **{leca['g2']['present']} of {len(fams)} families "
         f"are reconstructed present in the last eukaryotic common ancestor (LECA)**, "
         f"{len(leca_robust)} of them under all three costs; "
         f"{tot['g2'][0]} gains and {tot['g2'][1]} losses are stated, "
         f"{robust['gain']} gains and {robust['loss']} losses under all three costs. "
         f"**Of the {len(losses)} stated losses, {strength['controlled']} are controlled by S5's genome "
         f"sweep, {strength['contradicted']} contradicted by it (all S5 proteome misses), and "
         f"{strength['proteome-only']} are proteome-only** — the reconstruction's absences are, "
         "with few exceptions, annotation-level (D46).", "",
         "## Inputs", "",
         f"- **Species tree** (D15): {tree['source']}, retrieved {tree['retrieved']}, archived under "
         "`<data root>/raw_api/s10/taxonomy/`; `species_tree.nwk`, `species_tree.tsv`. Root children: "
         + ", ".join(tree["root_children"]) + ".",
         f"- **Off the tree** (addendum (c)): " + "; ".join(f"{o} (taxid {t})" for o, t, _ in tree["orders_unplaced"])
         + " — NCBI places each in an order already on the tree.",
         f"- **Characters** (`characters.tsv`, D50 (2)): present {states['1']:,}, absent {states['0']:,}, "
         f"missing {states['?']:,} (" + ", ".join(f"{k}: {v:,}" for k, v in why.most_common())
         + f"; BUSCO floor {BUSCO_FLOOR:.0f} %).",
         f"- **S4 species on the tree** for the overlay: {tree['s4_species_on_tree']} of "
         f"{tree['s4_species']} eukaryotic panel species (off: "
         + ", ".join(r["species"] for r in place if r["on_tree"] == "0") + ").", "",
         "## Costs compared", "",
         "| cost | LECA present | absent | ambiguous | stated gains | stated losses | ambiguous edges |",
         "|---|---:|---:|---:|---:|---:|---:|"]
    for k in GAINS:
        L.append(f"| {k} | {leca[k]['present']} | {leca[k]['absent']} | {leca[k]['ambiguous']} | "
                 f"{tot[k][0]} | {tot[k][1]} | {tot[k][2]} |")
    L += ["", f"**Present at LECA under all three costs ({len(leca_robust)})**: " + ", ".join(leca_robust) + ".",
          f"**Present at LECA under the primary cost only ({len(leca_cost)})**: " + ", ".join(leca_cost) + ".",
          f"**Absent at LECA under all three ({len(absent_all)})** — lineage-specific origins.", "",
          "*Reading*: a LECA presence here is a profile call on proteins outside animals. Where a "
          "deep-lineage member is a family member by descent or an architecture-matched relative "
          "(plant/fungal K2P-like two-pore channels, plant CNBD channels called CNG) is a tree "
          "question this reconstruction does not ask (emergent row).", "",
          "## Ancestral nodes (primary cost)", "",
          "| node | present | absent | ambiguous | present under all costs |", "|---|---:|---:|---:|---:|"]
    for n in NODES:
        c = node_counts[n]
        L.append(f"| {n} | {c['1']} | {c['0']} | {c['?']} | {node_stable[n]} |")
    L += ["", "## Families with the most stated losses (primary)", "",
          "| family | present orders | stated losses | robust | gains | LECA (g1 / g2 / Dollo) |",
          "|---|---:|---:|---:|---:|---|"]
    rob_f = Counter(e["family"] for e in prim if e["event"] == "loss" and e["robust"] == "1")
    for f in top_loss:
        r = fam[(f, PRIMARY)]
        L.append(f"| {f} | {r['present']} | {r['losses']} | {rob_f[f]} | {r['gains']} | "
                 + " / ".join(fam[(f, k)]["leca"] for k in GAINS) + " |")
    L += ["", "## Losses against S5's genome cells (D50 (4))", "",
          "| family | loss on | orders below | strength | S5 verdicts in the absent region |",
          "|---|---|---:|---|---|"]
    for r in sorted(losses, key=lambda r: (r["strength"], r["family"])):
        if r["strength"] != "proteome-only":
            L.append(f"| {r['family']} | {r['child']} | {r['child_tips']} | {r['strength']} | {r['s5_verdicts']} |")
    L += ["", "Every `contradicted` loss is a cell S5 already reported as a proteome miss "
          "(`genome_found`) or a genome-only presence (*Cornu*): the proteome-level reconstruction "
          "states a loss the genome refutes, which is the reason D46 keeps the two strengths apart.", "",
          "## Not done here", "",
          "- Q4 (animal MscS) and the emergent absence checks — **S10b**, under D50 (6).",
          "- Family-specific readings of gains (horizontal transfer, contamination: e.g. the one "
          "sponge order with a KcsA-like call) — S10b / S21.", "",
          "## Files", "", "`species_tree.{nwk,tsv,json}`, `s4_species_orders.tsv`, `characters.tsv`, "
          "`families.tsv` (per family × cost), `nodes.tsv` (every internal node × family × cost), "
          "`events.tsv` (every edge where a change is possible; `stated`, `robust`), `losses.tsv` "
          "(primary stated losses with strength), `figures/repertoire.png`."]
    (OUT / "report.md").write_text("\n".join(L) + "\n")
    summ = {"orders": tree["tips"], "families": len(fams), "leca": leca, "leca_robust": leca_robust,
            "totals": tot, "robust_events": dict(robust), "loss_strength": dict(strength),
            "characters": dict(states)}
    (OUT / "summary.json").write_text(json.dumps(summ, indent=2) + "\n")
    print("\n".join(L[:12]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
