"""selftest_s9.py — S9's offline invariants (filter atlas statistics), called
by `selftest.py`."""

from __future__ import annotations


def run(check) -> None:
    from s7_newick import parse
    from s9_atlas import compare
    from s9_lib import Fitch, congruence, is_missing, prune
    t = parse("(((a:1,b:1)90:1,(c:1,d:1)80:1)70:1,(e:1,f:1)99:1,g:1);")
    st = {"a": "K", "b": "K", "c": "Y", "d": "Y", "e": "Y", "f": "Y", "g": "Y"}
    check("S9: Fitch counts one change for a clade-bound state",
          Fitch(t).length(st), 1)
    st2 = {"a": "K", "b": "Y", "c": "K", "d": "Y", "e": "K", "f": "Y", "g": "Y"}
    check("S9: Fitch counts each scattered carrier", Fitch(t).length(st2), 3)
    check("S9: a missing tip costs nothing",
          Fitch(t).length({**st, "c": None}), 1)
    check("S9: RI is 1 when every state arises once", congruence(t, st)["ri"], 1.0)
    check("S9: pruning keeps the named leaves and nothing else",
          sorted(prune(t, {"a", "c", "g"}).leaves()), ["a", "c", "g"])
    check("S9: gap, unread and X strings are missing data",
          [is_missing(s) for s in ("DEKA", "DE-A", "?EKA", "DXKA", "")],
          [False, True, True, True, True])
    check("S9: projection comparison — identical, identical where both read, differ",
          [compare("DEKA", "DEKA"), compare("?EKA", "-EKA"), compare("EEDD", "ETDD"),
           compare("AN", "--PN")],
          ["yes", "yes_read", "no", "not_comparable"])
