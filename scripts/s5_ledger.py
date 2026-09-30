"""s5_ledger.py — the per-cell ledger: species × census family, genome verdict.

Every cell of S3b's presence matrix gets a genome status from its loci:

  found       a locus the S3a profiles call to this family (D32, high/medium)
  partial     a locus this family's baits win that the profiles do not call
              (module / superfamily-only / low score) — evidence, not a call
  gap         a partial locus at a contig edge or beside an N-run (D4's local bar)
  no_locus    neither

**Present cells are the positive control.** A cell S3b called at high
confidence in the proteome should be `found` in the genome; the fraction that
is, per genome, is that genome's recall, with the genome's own species'
baits excluded (D29). A zero cell's absence is only as good as that number.

Zero cells then get a verdict:

  genome_found     the genome carries a called locus the proteome has no call for
  partial / gap    as above
  trace            no locus; a tblastn HSP outside every recorded locus
  absent           no locus, no trace, contiguity bar met, control passed
  absent_below_bar no locus, no trace, but the assembly's N50 is below the
                   family's measured gene span (D4)
  uncontrolled     no locus, no trace, the genome's control recall is below
                   CONTROL_FLOOR — no absence can be read
  no_locus_unrescued  a non-informative zero cell (the family is absent from
                   the whole group); miniprot only, tblastn not run
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from src.catalogue import registry  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402
from s3b_lib import s3b_dir  # noqa: E402
from s5_lib import BAITS_TSV, OUT_DIR, s5_dir  # noqa: E402

NGAP_RUN = 100
CALLED = {"high", "medium"}
PARTIAL_CALLS = {"module", "superfamily_only", "low_score", "ambiguous", "no_hit"}

CELL_FIELDS = ["species", "group", "family", "superfamily", "proteome_records",
               "proteome_high", "proteome_band", "informative", "matched", "control", "genome",
               "n_found", "best_confidence", "n_strong", "n_partial", "n_gap", "n_traces",
               "trace_best_bits", "verdict", "found_loci"]


def proteome_cells() -> dict:
    """(species, family) → count of all / high-confidence S3b profile calls."""
    calls = defaultdict(Counter)
    for r in read_tsv(s3b_dir() / "census_v3.tsv.gz"):
        if r["p_call"] == "family":
            calls[(r["species"], r["p_family"])]["all"] += 1
            if r["p_confidence"] == "high":
                calls[(r["species"], r["p_family"])]["high"] += 1
    return calls


def band_by_species() -> Counter:
    """(species, family) → entries of that species with a band call naming it."""
    sp_of = {r["target"]: r["species"]
             for r in read_tsv(s3b_dir() / "panel_universe.tsv.gz")}
    out = Counter()
    for r in csv.DictReader(gzip.open(s3b_dir() / "panel_profile_calls.tsv.gz", "rt"),
                            delimiter="\t"):
        if r["p_call"] in ("superfamily_only", "ambiguous") and r["band"]:
            for f in r["band"].split(","):
                out[(sp_of.get(r["target"], ""), f)] += 1
    return out


def load_loci(assembly: str) -> list[dict]:
    """S5b's loci, plus census revision r4's (`r4/loci.tsv.gz`) where they
    bear on an r4 family: called to one, or a partial on one of its baits.

    An r4 bait that lands on an old family's gene (a KChIP bait on
    calmodulin) is dropped here, so a new bait panel can never change an old
    family's cell (D43). `r4_dropped_loci()` counts what was dropped.
    """
    p = s5_dir(assembly) / "loci.tsv.gz"
    old = read_tsv(p) if p.exists() else []
    return old + [L for L in _r4_loci(assembly) if _r4_keep(L)]


def _r4_loci(assembly: str) -> list[dict]:
    p = s5_dir(assembly) / "r4" / "loci.tsv.gz"
    return read_tsv(p) if p.exists() else []


def _r4_keep(L: dict) -> bool:
    from s5_lib import r4_families
    r4 = r4_families()
    return (L["p_family"] in r4 if L["p_call"] == "family" else L["bait_family"] in r4)


def r4_dropped_loci(assembly: str) -> list[dict]:
    return [L for L in _r4_loci(assembly) if not _r4_keep(L)]


def genome_status(fam: str, loci: list[dict]) -> dict:
    found = [L for L in loci if L["p_call"] == "family" and L["p_family"] == fam
             and L["p_confidence"] in CALLED]
    partial = [L for L in loci if L["bait_family"] == fam
               and L["p_call"] in PARTIAL_CALLS]
    gap = [L for L in partial if L["at_edge"] == "1"
           or int(L["longest_n_run"] or 0) >= NGAP_RUN]
    conf = ("high" if any(L["p_confidence"] == "high" for L in found)
            else "medium" if found else "")
    strong = [L for L in found if L["p_confidence"] == "high"
              and L.get("call_intact") == "1"]
    status = ("found" if found else "gap" if gap else "partial" if partial
              else "no_locus")
    return {"genome": status, "n_found": len(found), "best_confidence": conf,
            "n_strong": len(strong), "n_partial": len(partial), "n_gap": len(gap),
            "found_loci": ";".join(L["locus"] for L in found[:5])}


def informative_cells(calls: dict, group_of: dict) -> set[tuple[str, str]]:
    """(group, family) pairs where some species has a high-confidence call."""
    return {(group_of[s], f) for (s, f), c in calls.items()
            if c["high"] and s in group_of}


def controls(cells: list[dict]) -> list[dict]:
    """Per genome: high-confidence proteome cells recovered as `found`."""
    out = []
    for sp in sorted({c["species"] for c in cells}):
        mine = [c for c in cells if c["species"] == sp and c["control"]]
        rec = [c for c in mine if c["genome"] == "found"]
        part = [c for c in mine if c["genome"] in ("partial", "gap")]
        miss = [c for c in mine if c["genome"] == "no_locus"]
        out.append({"species": sp, "group": cells[[c["species"] for c in cells].index(sp)]["group"],
                    "control_cells": len(mine), "found": len(rec),
                    "partial": len(part), "no_locus": len(miss),
                    "recall": round(len(rec) / len(mine), 4) if mine else "",
                    "missed": ",".join(sorted(c["family"] for c in miss)),
                    "partial_families": ",".join(sorted(c["family"] for c in part))})
    return out


def matched_baits() -> set[tuple[str, str, str]]:
    """(group, family, bait species) for every bait — a cell is *matched*
    when a bait of its family comes from another species of its group.

    Measured in the pilot: neither miniprot nor tblastn reaches a relative
    from another group at 20–30 % identity (7 animal/ciliate TPC baits, 0
    HSPs on Arabidopsis TPC1), so detection is a property of the nearest
    bait. An absence is judged only where the nearest bait is in-group, and
    the control that bounds it is measured on cells in the same condition.
    """
    from s5_lib import load_bait_meta
    return {(r["group"], r["family"], r["species"])
            for r in load_bait_meta().values()}


def build_cells(species_rows: list[dict], man: dict) -> list[dict]:
    calls = proteome_cells()
    band = band_by_species()
    group_of = {s: r["group"] for s, r in man.items()}
    inf = informative_cells(calls, group_of)
    mb = matched_baits()
    # a genome-only species has no proteome; its control cells are the
    # group's core families — high-confidence in every proteome species of it
    prot_sp = defaultdict(set)
    for s, r in man.items():
        if r["status"] == "proteome":
            prot_sp[r["group"]].add(s)
    core = {(g, f) for g, sps in prot_sp.items()
            for f in {k[1] for k in calls}
            if all(calls.get((s, f), Counter())["high"] for s in sps)}
    fams = registry.census_families()
    cells = []
    for run in species_rows:
        if run.get("note", "").startswith("FAILED"):
            continue
        sp, loci = run["species"], load_loci(run["assembly"])
        for f in fams:
            c = calls.get((sp, f.key), Counter())
            row = {"species": sp, "group": run["group"], "family": f.key,
                   "superfamily": f.superfamily, "proteome_records": c["all"],
                   "proteome_high": c["high"], "proteome_band": band[(sp, f.key)],
                   "informative": int((run["group"], f.key) in inf),
                   "matched": int(any(g == run["group"] and fam == f.key and bs != sp
                                      for g, fam, bs in mb))}
            genome_only = man[sp]["status"] != "proteome"
            row["control"] = int(c["high"] > 0 if not genome_only
                                 else (run["group"], f.key) in core)
            row.update(genome_status(f.key, loci))
            cells.append(row)
    return cells
