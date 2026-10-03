"""selftest_s11b.py — S11b's offline invariants (D53: the non-reversible model
mapping, the rootstrap / root-test readers, the repeat constraint and the AU
verdict), called by `selftest.py`."""

from __future__ import annotations

import tempfile
from pathlib import Path


def run(check) -> None:
    from s11b_repeats import au_table, constraint, tip_class, verdict
    from s11b_roots import canon, nonrev_model, rootstrap_tree

    check("S11b: NQ.pfam keeps S7b's +I/+R terms and drops +F",
          [nonrev_model(m) for m in ("Q.pfam+F+R10", "JTT+I+R5", "LG+F+I+G4")],
          ["NQ.pfam+R10", "NQ.pfam+I+R5", "NQ.pfam+I+G4"])
    check("S11b: repeat classes — I–IV for 4×6TM, tI/tII for TPC",
          (tip_class("nav", 3), tip_class("tpc", 2)), ("III", "tII"))
    tips = ["I__nav__a", "II__nav__a", "III__nav__a", "IV__nav__a",
            "I__cav__b", "III__cav__b"]
    check("S11b: H13's constraint is the bipartition {I,III}|{II,IV} over 4×6TM tips",
          constraint(tips, "H13"),
          "((I__nav__a,III__nav__a,I__cav__b,III__cav__b),(II__nav__a,IV__nav__a));\n")

    from s11b_repeats import free_tip, satisfies
    check("S11b: a tip is freed only when the constraint would hold exactly 64 taxa",
          (free_tip([f"I__x{i}" for i in range(64)]) != "", free_tip(tips)), (True, ""))
    check("S11b: the freed constraint omits the free tip",
          "I__cav__b" in constraint(tips, "H13", free="I__cav__b")
          .strip("();\n").replace("(", "").replace(")", "").split(","), False)
    check("S11b: a constrained tree is checked for the full bipartition afterwards",
          (satisfies("((I__a,III__a),(II__a,IV__a),tI__b);", "H13"),
           satisfies("((I__a,II__a),(III__a,IV__a),tI__b);", "H13")), (True, False))
    check("S11b: a free TPC tip nested inside one side does not break the bipartition",
          satisfies("(((I__a,tI__b),III__a),(II__a,IV__a),tII__b);", "H13"), True)

    nex = ('#NEXUS\nbegin trees;\n  tree tree_1 = ((a:1[&id="3",rootstrap="10"],'
           'b:1[&id="4",rootstrap="5"]):1[&id="2",rootstrap="97"],(c:1[&id="6",'
           'rootstrap="0"],d:1[&id="7",rootstrap="0"]):1[&id="5",rootstrap="97"])'
           ':0[&id="0",rootstrap="97"];\nend;\n')
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "t.nex"
        p.write_text(nex)
        tree, below = rootstrap_tree(p)
    leaves = frozenset("abcd")
    check("S11b: the rootstrap reader maps a branch id to the leaves below it",
          (sorted(below["2"][0]), below["2"][1], sorted(below["4"][0])),
          (["a", "b"], 97.0, ["b"]))
    check("S11b: the root's two half-edges are one split (canonical orientation)",
          canon(below["2"][0], leaves) == canon(below["5"][0], leaves), True)

    iq = ("USER TREES\n----------\n\nTree      logL    deltaL  bp-RELL    p-KH     p-SH"
          "       c-ELW       p-AU\n----\n  1 -100.0   0  0.6 +  0.7 +  1 +  0.6 +  0.8 + \n"
          "  2 -110.5  10.5  0.01 -  0.02 -  0.03 -  0.01 -  0.004 - \n")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "au.iqtree"
        p.write_text(iq)
        rows = au_table(p)
    check("S11b: the AU table reader takes p-AU from the eighth number",
          [(r["tree"], r["p_AU"]) for r in rows], [(1, 0.8), (2, 0.004)])
    au = [{"unit": "u", "tree": "ML", "rejected": "False"},
          {"unit": "u", "tree": "H13", "rejected": "False"},
          {"unit": "u", "tree": "H12", "rejected": "True"},
          {"unit": "u", "tree": "H14", "rejected": "False"}]
    check("S11b: two unrejected pairings → unresolved, never a pick (D53 (4))",
          verdict(au, "u"), "unresolved:H13,H14")
    au[3]["rejected"] = "True"
    check("S11b: one unrejected pairing names the order", verdict(au, "u"), "H13")
