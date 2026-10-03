"""S10b — *Daphnia pulex* on its current assembly, under S5's instrument (D51 (3)).

    bin/envpy scripts/s10b_daphnia.py pick     # the assembly, by S4's assembly_rank
    bin/envpy scripts/s10b_daphnia.py run      # fetch → miniprot (both bait panels) → calls
    bin/envpy scripts/s10b_daphnia.py cells    # S5's cells/rescue/control/D4 → compare

Nothing here is a new method: every step calls S5's own function, with the
assembly swapped. The 2011 assembly's cells stay in `results/genome_sweep/`;
this writes `results/repertoire/daphnia_*.tsv` beside them.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import s4_proteome_lib as s4  # noqa: E402
import s5_sweep  # noqa: E402
import s5r4_sweep  # noqa: E402
from s3_hmm_lib import read_tsv, sha256, write_tsv  # noqa: E402
from s5_genome_io import write_json  # noqa: E402
from s5_ledger import build_cells, controls  # noqa: E402
from s5_lib import (BAITS_FAA, BAITS_R4_FAA, OUT_DIR, load_bait_meta,  # noqa: E402
                    manifest)
from s5_verdict import d4_bar, rescue, spans, verdict  # noqa: E402
from src.catalogue import registry  # noqa: E402

SPECIES = "Daphnia pulex"
TAXID = 6669
OUT = ROOT / "results" / "repertoire"
PICK = OUT / "daphnia_assembly.json"
CELLS = OUT / "daphnia_cells.tsv"
COMPARE = OUT / "daphnia_compare.tsv"


def cmd_pick(_a) -> dict:
    asms = [s4.flatten_assembly(r) for r in s4.fetch_taxon_assemblies(TAXID)]
    asms.sort(key=s4.assembly_rank, reverse=True)
    old = s5_sweep.sweep_assembly(manifest()[SPECIES])
    best = asms[0]
    rec = {"species": SPECIES, "old_assembly": old, "new_assembly": best["accession"],
           "rank": list(s4.assembly_rank(best)), "level": best["level"],
           "scaffold_n50": best["scaffold_n50"], "candidates": len(asms),
           "rule": "S4 assembly_rank: annotated > RefSeq > level > scaffold N50"}
    write_json(PICK, rec)
    print(json.dumps(rec, indent=1))
    return rec


def _row() -> dict:
    pick = json.loads(PICK.read_text())
    return {**manifest()[SPECIES], "current_assembly": pick["new_assembly"]}


def cmd_run(_a) -> None:
    row, meta = _row(), load_bait_meta()
    rec = s5_sweep.run_one(SPECIES, row, sha256(BAITS_FAA), meta)
    print("S5b baits:", {k: rec[k] for k in ("loci", "family_called_loci", "n50")})
    r4 = s5r4_sweep.run_one(SPECIES, row, sha256(BAITS_R4_FAA), meta)
    print("r4 baits:", {k: r4[k] for k in ("loci", "r4_family_loci")})
    write_json(OUT / "daphnia_run.json", {"s5": rec, "r4": r4})


def cells_for(run: dict) -> tuple[list[dict], dict]:
    """S5's verdict pipeline, verbatim, on one genome."""
    man, meta = manifest(), load_bait_meta()
    cells = build_cells([run], man)
    traces = rescue(run, cells, meta)
    per = defaultdict(list)
    for t in traces:
        per[t["family"]].append(t)
    for c in cells:
        ts = per.get(c["family"], [])
        c["n_traces"] = len(ts)
        c["trace_best_bits"] = max((t["bitscore"] for t in ts), default="")
        c["detected"] = int(c["genome"] != "no_locus" or bool(ts))
    ctl = controls(cells)[0]
    mine = [c for c in cells if c["control"]]
    m = [c for c in mine if c["matched"]]
    ctl["matched_detection"] = (round(sum(c["detected"] for c in m) / len(m), 4)
                                if m else "")
    ctl["matched_cells"] = len(m)
    ctl["matched_undetected"] = ",".join(sorted(c["family"] for c in m
                                                if not c["detected"]))
    span = spans()
    band = {f.key: f.length_band_aa[1] for f in registry.CATALOGUE.values()}
    for c in cells:
        bp, src = d4_bar(c["family"], c["group"], span, band.get(c["family"], 0))
        c["bar_bp"], c["bar_source"] = ("" if bp is None else bp), src
        c["verdict"] = verdict(c, run, {SPECIES: ctl})
    return cells, ctl


def cmd_cells(_a) -> None:
    side = json.loads((OUT / "daphnia_run.json").read_text())["s5"]
    run = {"species": SPECIES, "group": side["group"], "assembly": side["assembly"],
           "n50": side["n50"], "status": side["status"]}
    cells, ctl = cells_for(run)
    keep = ["family", "superfamily", "proteome_records", "proteome_high", "control",
            "matched", "genome", "n_found", "best_confidence", "n_partial", "n_gap",
            "n_traces", "verdict", "found_loci", "bar_bp", "bar_source"]
    write_tsv(CELLS, keep, cells)
    old = {c["family"]: c for c in read_tsv(OUT_DIR / "cells.tsv")
           if c["species"] == SPECIES}
    old_ctl = next(r for r in read_tsv(OUT_DIR / "controls.tsv")
                   if r["species"] == SPECIES)
    rows = []
    for c in cells:
        o = old.get(c["family"], {})
        rows.append({"family": c["family"], "proteome_high": c["proteome_high"],
                     "old_genome": o.get("genome", ""), "old_verdict": o.get("verdict", ""),
                     "new_genome": c["genome"], "new_verdict": c["verdict"],
                     "new_best_confidence": c["best_confidence"],
                     "new_found_loci": c["found_loci"],
                     "changed": int(o.get("verdict", "") != c["verdict"])})
    write_tsv(COMPARE, list(rows[0]), rows)
    summary = {"old_assembly": json.loads(PICK.read_text())["old_assembly"],
               "new_assembly": run["assembly"], "new_n50": run["n50"],
               "old_matched_detection": old_ctl["matched_detection"],
               "new_matched_detection": ctl["matched_detection"],
               "new_matched_cells": ctl["matched_cells"],
               "new_matched_undetected": ctl["matched_undetected"],
               "replaces_old": (ctl["matched_detection"] != "" and
                                float(ctl["matched_detection"]) >=
                                float(old_ctl["matched_detection"] or 0)),
               "changed_cells": sum(r["changed"] for r in rows)}
    write_json(OUT / "daphnia_summary.json", summary)
    print(json.dumps(summary, indent=1))


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("pick", cmd_pick), ("run", cmd_run), ("cells", cmd_cells)):
        sub.add_parser(name).set_defaults(fn=fn)
    a = ap.parse_args()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
