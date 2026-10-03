"""selftest_s11.py — S11a's offline invariants (D52: reconciliation, roots,
stated duplications, ohnologue windows), called by `selftest.py`."""

from __future__ import annotations


def _st(spec):
    """Toy species tree from nested tuples (name, [children]) / tip strings."""
    from s10_lib import TNode
    from s11_recon import stree_from

    def mk(x):
        if isinstance(x, str):
            return TNode(x, tip=x)
        return TNode(x[0], children=[mk(c) for c in x[1]])
    return stree_from(mk(spec))


def run(check) -> None:
    from s7_newick import parse
    from s11_recon import Reconciler, utree
    from s11_duplications import accepted, window

    st = _st(("root", [("X", ["A", "B"]), "C"]))
    sp = {"a1": "A", "a2": "A", "b1": "B", "b2": "B", "c1": "C"}
    t = utree(parse("((a1,b1)100,(a2,b2)80,c1);"))
    rc = Reconciler(t, st, sp)
    opt, best = rc.optimal()
    check("S11a: one duplication, no loss, a unique root on the outgroup species' edge",
          (len(opt), best), (1, (1, 0)))
    ids = rc.dup_ids(opt[0])
    (ident, s), = ids.items()
    check("S11a: the duplication maps to the LCA of its species",
          st.names[s.m], "X")
    check("S11a: a duplication with a child edge under UFBoot 95 is not supported",
          rc.supported(ident, s), False)
    t2 = utree(parse("((a1,b1)100,(a2,b2)99,c1);"))
    rc2 = Reconciler(t2, st, sp)
    (i2, s2), = rc2.dup_ids(rc2.optimal()[0][0]).items()
    check("S11a: both child edges ≥ 95 → supported; species-overlap score 1.0",
          (rc2.supported(i2, s2), rc2.score(s2)), (True, 1.0))

    star = _st(("root", ["A", "B", "C"]))
    rc3 = Reconciler(utree(parse("(a1,b1,c1);")), star, sp)
    opt3, best3 = rc3.optimal()
    check("S11a: a polytomy never yields a negative loss; every root optimal when tied",
          (len(opt3), best3), (3, (0, 0)))

    t4 = utree(parse("((a1,OG_x)90,b1,c1);"), drop={"OG_x"})
    check("S11a: a pruned outgroup leaves no degree-2 node",
          sorted(len(v) for k, v in t4.adj.items() if k not in t4.leaf), [3])

    # stated = present under every optimal root: three paralogues of one species
    # tie on every root (two duplications each) but share no duplication node
    rc5 = Reconciler(utree(parse("(a1,a2,a3);")), st, {"a1": "A", "a2": "A", "a3": "A"})
    opt5, best5 = rc5.optimal()
    per = [set(rc5.dup_ids(e)) for e in opt5]
    check("S11a: a duplication not present under every optimal root is not stated",
          (len(opt5), best5, len(set.intersection(*per)), len(set.union(*per)) > 0),
          (3, (2, 0), 0, True))

    vt = _st(("root", [("Chordata", ["Branchiostoma", ("Vertebrata", [
        "Petromyzon", ("Gnathostomata", ["Callorhinchus", ("Euteleostomi", [
            ("Mammalia", ["Homo", "Mus"]), ("Clupeocephala", ["Danio", "Takifugu"])])])])])]))
    names = {n: i for i, n in enumerate(vt.names)}
    check("S11a: ohnologue windows by mapped taxon (D52 (6))",
          [window(vt, names[n]) for n in
           ("Vertebrata", "Gnathostomata", "Euteleostomi", "Clupeocephala", "Chordata",
            "Mammalia")],
          ["2R", "2R", "bony_vertebrate", "3R", "older", "younger"])
    row = {"family": "x", "root_rule": "rooted", "outgroup_monophyletic": "True",
           "n_outgroup": "3", "root_ufboot": "89.0"}
    check("S11a: an S7b root is used only if one clade with UFBoot ≥ 95, or one sequence",
          (accepted(row), accepted({**row, "root_ufboot": "100"}),
           accepted({**row, "n_outgroup": "1", "root_ufboot": ""}),
           accepted({**row, "family": "kv_eag", "root_ufboot": "100"})),
          (False, True, True, False))
