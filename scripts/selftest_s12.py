"""selftest_s12.py — S12's offline invariants (D54: the D48 edge rule reused on
re-read pair tables, the isoform map's canonical identity), called by
`selftest.py`."""

from __future__ import annotations


def run(check) -> None:
    from s8_fold_network import edge_rows
    from s12_lib import isoform_to_canonical

    def pair(a, b, sa, sb, tm):
        return {"a": a, "b": b, "sf_a": sa, "sf_b": sb, "tm_avg": tm}
    # iGluR–P-loop at 0.6, each the other's best → supported; P-loop's best other
    # is iGluR, so Ca-release–P-loop at 0.55 ranks 2nd for P-loop → not distinguished.
    pairs = [pair("ampa", "nav", "iglur", "ploop", 0.6),
             pair("itpr", "nav", "ca_release", "ploop", 0.55),
             pair("itpr", "ampa", "ca_release", "iglur", 0.3)]
    v = {(r["a"], r["b"]): r["verdict"] for r in edge_rows(pairs)}
    check("S12: edge_rows — mutual best at ≥ 0.5 is supported",
          v[("iglur", "ploop")], "supported")
    check("S12: edge_rows — ≥ 0.5 but not mutual best is not distinguished",
          v[("ca_release", "ploop")], "not_distinguished")
    check("S12: edge_rows — an edge with no measured pair is unmeasured",
          v[("innexin_like", "connexin")], "unmeasured")
    check("S12: canonical and '-1' isoform numbering map to themselves",
          (isoform_to_canonical("P1", "P1").get(17), isoform_to_canonical("P1-1", "P1").get(5)),
          (17, 5))
