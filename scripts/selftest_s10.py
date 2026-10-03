"""selftest_s10.py — S10a's offline invariants (the D50 reconstruction), called
by `selftest.py`."""

from __future__ import annotations


def run(check) -> None:
    from s10_lib import TNode, build_tree, events, sankoff

    def T(name, *ch):
        return TNode(name, children=list(ch))

    def L(n):
        return TNode(n, tip=n)

    t = T("R", T("X", L("a"), L("b")), T("Y", L("c"), L("d")))
    st = {"a": 1, "b": 1, "c": 1, "d": 0}
    r2 = sankoff(t, st, 2)
    check("S10: g = 2 — root present, one stated loss on d",
          (r2["score"], r2["root"], [(e["child"], e["event"], e["stated"]) for e in events(r2)]),
          (3, 1, [("d", "loss", 1)]))
    r1 = sankoff(t, st, 1)
    check("S10: g = 1 — tied optima leave the root and the events unstated",
          (r1["score"], r1["root"], sorted({e["stated"] for e in events(r1)})), (2, None, [0]))
    star = T("R", L("a"), L("b"), L("c"), L("d"))
    check("S10: a polytomy is scored exactly (star, two present)",
          sankoff(star, {"a": 1, "b": 0, "c": 1, "d": None}, 2)["score"], 3)
    check("S10: a missing tip costs nothing in either state",
          sankoff(t, {**st, "d": None}, 2)["score"], 2)
    check("S10: Dollo-like cost allows one origin only",
          sankoff(t, {"a": 1, "b": 0, "c": 1, "d": 0}, 10 ** 6)["score"] // 10 ** 6, 1)
    lin = {"o1": [("1", "Eukaryota", "superkingdom"), ("2", "A", "clade"), ("3", "o1", "order")],
           "o2": [("1", "Eukaryota", "superkingdom"), ("2", "A", "clade"), ("4", "B", "clade"),
                  ("5", "o2", "order")]}
    tr = build_tree(lin)
    check("S10: species tree collapses unary chains, keeps the deeper name",
          (tr.name, sorted(c.name for c in tr.children)), ("A", ["o1", "o2"]))
