"""selftest_s10b.py — S10b's offline invariants (D51 (1): regions and their
lineage), called by `selftest.py`."""

from __future__ import annotations


def run(check) -> None:
    from s10b_lib import lineage_of_group, merge_regions, region_lineage

    def H(qs, qe, s, bits):
        return {"qs": qs, "qe": qe, "s": s, "bits": bits, "pident": 50.0}

    hs = [H(100, 400, "a", 90), H(350, 600, "b", 80), H(900, 700, "c", 50)]
    check("S10b: HSPs merge by contig overlap, either strand",
          [len(r) for r in merge_regions(hs)], [2, 1])
    ox = {"a": "1", "b": "2", "c": "3"}
    taxa = {"1": {"species": "x", "lineage": "Metazoa"},
            "2": {"species": "y", "lineage": "prokaryote"},
            "3": {"species": "z", "lineage": "prokaryote"}}
    check("S10b: a lineage inside D7's 10 % margin is unclear",
          region_lineage([H(1, 9, "a", 100), H(1, 9, "b", 95)], ox, taxa)["lineage"],
          "unclear")
    check("S10b: a lineage past the margin is the region's",
          region_lineage([H(1, 9, "a", 100), H(1, 9, "b", 80)], ox, taxa)["lineage"],
          "Metazoa")
    check("S10b: holozoan protists are not animals",
          (lineage_of_group("holozoa"), lineage_of_group("basal_metazoan")),
          ("other_eukaryote", "Metazoa"))
