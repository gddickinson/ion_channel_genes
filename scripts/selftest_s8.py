"""selftest_s8.py — S8b's offline invariants, called by `selftest.py`
(which is already over the file budget, so new checks live here)."""

from __future__ import annotations


def run(check) -> None:
    from s7_newick import parse
    from s8_parse import clade_test, clades, module_of
    from s8_report_trees import tested_groups
    t = parse("((a1:1,a2:1)90:1,((b1:1,x:1)80:1,b2:1)70:1,(o1:1,o2:1)99:1);")
    og = frozenset({"o1", "o2"})
    cl = clades(t, og)
    check("S8b: a family that is one clade (rooted sense) carries its UFBoot",
          clade_test(t, {"a1", "a2"}, og, cl)[:2], (True, 90.0))
    one, _, intr = clade_test(t, {"b1", "b2"}, og, cl)
    check("S8b: a split family names the fewest tips in the way", (one, sorted(intr)),
          (False, ["x"]))
    check("S8b: a rooted clade test never counts a side holding the outgroup",
          any(og & s for s, _ in cl), False)
    check("S8b: module suffix read from the tip name",
          (module_of("nav__Q14524__Homsap_m3"), module_of("kir__P1__Homsap")), ("m3", ""))
    rows = [{"unit": "u", "family": "nav", "module": ""}, {"unit": "u", "family": "nav",
            "module": "m1"}, {"unit": "u", "family": "kir", "module": ""}]
    check("S8b: multi-module families are tested per module, others whole",
          [(r["family"], r["module"]) for r in tested_groups(rows)],
          [("nav", "m1"), ("kir", "")])
