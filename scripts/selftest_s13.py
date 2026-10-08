"""selftest_s13.py — S13's offline invariants (D55), called by `selftest.py`.

Checks on refusal as much as on success: a CDS that encodes another protein
is refused, a disagreeing codon is masked rather than kept, a labelled set
that is not a clade raises (codeml would otherwise mark a larger one), and
the orthogroup rule's support floor and pruning bound behave as written.
"""

from __future__ import annotations


def run(check) -> None:
    from s13_lib import bh, lrt_p, validate_and_mask
    from s13_tree import UTree, drawn_for, orthogroup, paml_newick, prune

    # -- CDS validation (IP3R D36)
    cds = "ATGAAATTTGGG"                                   # M K F G
    check("S13: an exact CDS validates with nothing masked",
          validate_and_mask(cds + "TAA", "MKFG")[0], cds)
    check("S13: a CDS with one disagreeing codon in 4 is refused (< 99 %)",
          validate_and_mask(cds, "MKFA")[0], None)
    long_p = "M" + "K" * 199
    long_c = "ATG" + "AAA" * 198 + "TTT"                   # last residue F, not K
    masked, st = validate_and_mask(long_c, long_p)
    check("S13: a single disagreement in 200 is masked NNN, not kept",
          (masked[-3:], st["n_masked"]), ("NNN", 1))
    check("S13: an internal stop where the protein has X is masked",
          validate_and_mask("ATGTAAGGG", "MXG")[0], "ATGNNNGGG")
    check("S13: a CDS of a different length is refused",
          validate_and_mask(cds, "MKF")[0], None)

    # -- statistics
    check("S13: the 50:50 boundary mixture halves the χ²₁ p-value",
          round(lrt_p(-10.0, -8.0, 1, mixture=True)[1] * 2, 6),
          round(lrt_p(-10.0, -8.0, 1)[1], 6))
    check("S13: Benjamini–Hochberg on a known vector",
          [round(q, 4) for q in bh([0.01, 0.04, 0.03, 0.5])], [0.04, 0.0533, 0.0533, 0.5])

    # -- trees: pruning keeps the larger UFBoot on a merged edge
    t = UTree.from_newick("((A:1,(B:1,(C:1,x:1)60:1)99:1)80:1,D:1,E:1);")
    p = prune(t, {"A", "B", "C", "D", "E"})
    sup = {frozenset(s): v for s, v, _ in __import__("s13_tree").sides(p) if len(s) == 2}
    check("S13: pruning a tip merges edges and keeps the larger UFBoot",
          sup.get(frozenset({"B", "C"})), 99.0)

    # -- orthogroup rule: largest one-anchor side with UFBoot ≥ 95
    t = UTree.from_newick("(((h1:1,a:1)100:1,b:1)90:1,((h2:1,c:1)100:1,d:1)100:1,(h3:1,e:1)100:1);")
    g1, _ = orthogroup(t, "h1", {"h2", "h3"})
    g2, _ = orthogroup(t, "h2", {"h1", "h3"})
    check("S13: an orthogroup stops below a bounding edge under UFBoot 95",
          sorted(g1), ["a", "h1"])
    check("S13: an orthogroup takes the largest supported one-anchor side",
          sorted(g2), ["c", "d", "h2"])
    check("S13: a lone anchor's orthogroup is every tip",
          len(orthogroup(t, "h1", set())[0]), 8)
    try:
        orthogroup(t, "h1", {"h2"})
        two = False
    except ValueError:
        two = True
    check("S13: the rule refuses exactly two anchor tips (degenerate side)", two, True)

    # -- labels: within-branches exclude the stem; stems label one branch
    drawn = drawn_for(t, [frozenset(g2), frozenset(g1)])
    nwk = paml_newick(drawn, {frozenset(g2): "#1"})
    check("S13: within-labels mark every branch inside the set, never its stem",
          (nwk.count("#1"), "(h2#1,c#1)#1" in nwk, ")#1)" in nwk or nwk.endswith(")#1;")),
          (4, True, False))
    stem = paml_newick(drawn, {}, stem={frozenset(g1): "#1"})
    check("S13: a stem label marks exactly one branch, the set's own",
          (stem.count("#1"), "(h1,a)#1" in stem), (1, True))
    try:
        paml_newick(drawn, {frozenset({"h1", "c"}): "#1"})
        refused = False
    except ValueError:
        refused = True
    check("S13: a labelled set that is not a clade is refused", refused, True)
    check("S13: HyPhy labels use {Name}",
          "h2{Test}" in paml_newick(drawn, {frozenset(g2): "Test"}, style="hyphy"), True)

    # -- the codon map places one codon per aligned residue
    from s13_codon import codon_map
    check("S13: codon map — gaps become ---, residues take codons in order",
          codon_map({"s": "M-K"}, {"s": "ATGAAA"}), {"s": "ATG---AAA"})
