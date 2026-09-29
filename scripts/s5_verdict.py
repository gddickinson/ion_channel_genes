"""s5_verdict.py — rescue, controls, D4, and the verdict on every cell.

    python3 scripts/s5_verdict.py            # over every genome in genome_runs.tsv

1. **Rescue** (`s5_rescue.py`): one tblastn per genome, over every cell the
   sweep left at `no_locus` that is either a control cell (proteome-present
   at high confidence) or an informative zero cell. Control cells are
   rescued too, because an absence verdict requires the *whole* instrument
   to return nothing, so the control must measure the whole instrument.
2. **Controls**: per genome, the share of control cells the instrument
   *detects* (found / partial / gap / trace) — one minus the false-absence
   rate an absent zero cell inherits — and the share it *finds* (a profile
   call).
3. **D4**: a family's contiguity bar is its measured gene span
   (`spans.tsv`, from annotation — `s5_calibrate.py`); a genome whose
   searched-file N50 is below it cannot carry an absence.
4. **Verdict** per zero cell, as listed in `s5_ledger.py`.
"""

from __future__ import annotations

import hashlib
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s5_genome_io import fna_path  # noqa: E402
from s5_ledger import CELL_FIELDS, build_cells, controls, load_loci  # noqa: E402
from s5_lib import OUT_DIR, load_bait_meta, manifest  # noqa: E402
from s5_rescue import outside_loci, read_hsps, run_tblastn, traces  # noqa: E402

#: A genome whose instrument detects fewer than this share of its own
#: proteome-present families cannot carry an absence. Set from the pilot's
#: measured detection (see report); not tuned on the zero cells it judges.
CONTROL_FLOOR = 0.90
TRACE_FIELDS = ["species", "family", "contig", "start", "end", "n_hsps",
                "bitscore", "evalue", "pident", "bait", "cell_kind"]
CONTROL_FIELDS = ["species", "group", "control_cells", "found", "partial",
                  "no_locus", "recall", "missed", "partial_families",
                  "detected", "detection", "undetected", "matched_cells",
                  "matched_detected", "matched_detection", "matched_undetected",
                  "unmatched_cells", "unmatched_detected"]


def rescue(run: dict, cells: list[dict], meta: dict) -> list[dict]:
    """tblastn for this genome's rescuable no-locus cells → trace rows."""
    sp = run["species"]
    want = [c for c in cells if c["species"] == sp and c["genome"] == "no_locus"
            and (c["control"] or
                 (int(c["proteome_records"]) == 0 and c["informative"] == 1))]
    if not want:
        return []
    fams = {c["family"] for c in want}
    baits = sorted(b for b, m in meta.items()
                   if m["family"] in fams and m["species"] != sp)
    fna = fna_path(run["assembly"])
    tag = "rescue_" + hashlib.sha256("\n".join(baits).encode()).hexdigest()[:12]
    hsps = read_hsps(run_tblastn(run["assembly"], fna, baits, tag))
    hsps = [h for h in hsps if h["family"] in fams]
    kept = outside_loci(hsps, load_loci(run["assembly"]))
    kind = {c["family"]: ("control" if c["control"] else "zero")
            for c in want}
    rows = []
    for fam, regions in traces(kept).items():
        for t in regions:
            rows.append({**t, "species": sp, "cell_kind": kind[fam]})
    return rows


def spans() -> dict[str, int]:
    p = OUT_DIR / "spans.tsv"
    if not p.exists():
        return {}
    return {r["family"]: int(float(r["bar_bp"])) for r in read_tsv(p) if r["bar_bp"]}


def verdict(c: dict, run: dict, ctl: dict, bar: dict) -> str:
    if int(c["proteome_records"]) > 0:
        return "present"
    g = c["genome"]
    if g == "found":
        if run.get("status", "proteome") != "proteome":
            return "genome_present"       # no proteome to have missed it
        # a claim that the proteome missed a gene needs a high-confidence
        # call on an intact reading frame (no stop, no frameshift) — else a
        # retrocopy or a chained low-identity model reads as a missed gene
        return "genome_found" if int(c.get("n_strong", 0)) else "genome_weak"
    if g in ("partial", "gap"):
        return g
    if int(c["n_traces"]):
        return "trace"
    if c["informative"] != 1:
        return "no_locus_unrescued"
    if c["matched"] != 1:
        return "unmatched"
    det = ctl.get(c["species"], {}).get("matched_detection")
    if det == "" or det is None or float(det) < CONTROL_FLOOR:
        return "uncontrolled"
    need = bar.get(c["family"])
    if need is None:
        return "absent_bar_unmeasured"
    return "absent" if int(run["n50"]) >= need else "absent_below_bar"


def main() -> int:
    man = manifest()
    runs = [r for r in read_tsv(OUT_DIR / "genome_runs.tsv")
            if not r.get("note", "").startswith("FAILED")]
    meta = load_bait_meta()
    cells = build_cells(runs, man)
    trace_rows = []
    for run in runs:
        print(f"rescue {run['species']}", flush=True)
        trace_rows += rescue(run, cells, meta)
    per_cell = defaultdict(list)
    for t in trace_rows:
        per_cell[(t["species"], t["family"])].append(t)
    for c in cells:
        ts = per_cell.get((c["species"], c["family"]), [])
        c["n_traces"] = len(ts)
        c["trace_best_bits"] = max((t["bitscore"] for t in ts), default="")
        c["detected"] = int(c["genome"] != "no_locus" or bool(ts))
    ctl_rows = controls(cells)
    for r in ctl_rows:
        mine = [c for c in cells if c["species"] == r["species"] and c["control"]]
        det = [c for c in mine if c["detected"]]
        r["detected"] = len(det)
        r["detection"] = round(len(det) / len(mine), 4) if mine else ""
        r["undetected"] = ",".join(sorted(c["family"] for c in mine
                                          if not c["detected"]))
        m = [c for c in mine if c["matched"]]
        u = [c for c in mine if not c["matched"]]
        r.update(matched_cells=len(m), matched_detected=sum(c["detected"] for c in m),
                 matched_detection=(round(sum(c["detected"] for c in m) / len(m), 4)
                                    if m else ""),
                 matched_undetected=",".join(sorted(c["family"] for c in m
                                                    if not c["detected"])),
                 unmatched_cells=len(u), unmatched_detected=sum(c["detected"] for c in u))
    ctl = {r["species"]: r for r in ctl_rows}
    bar = spans()
    run_of = {r["species"]: r for r in runs}
    for c in cells:
        c["verdict"] = verdict(c, run_of[c["species"]], ctl, bar)
    write_tsv(OUT_DIR / "cells.tsv", CELL_FIELDS + ["detected"], cells)
    write_tsv(OUT_DIR / "controls.tsv", CONTROL_FIELDS, ctl_rows)
    write_tsv(OUT_DIR / "traces.tsv", TRACE_FIELDS, trace_rows)
    for r in ctl_rows:
        print(r["species"], r["recall"], r["detection"], r["undetected"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
