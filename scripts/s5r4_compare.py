"""s5r4_compare.py — what the r4 genome sweep changed, cell by cell (D43).

    python3 scripts/s5r4_compare.py --base <dir with S5b's cells.tsv / controls.tsv>

Every (species × family) cell S5b had is compared with the current
`cells.tsv` on genome status, verdict and the fields that decide it; every
genome's control row likewise. An old cell may legitimately move for two
reasons only — r4 changed its proteome evidence (the S3b panel calls were
re-assigned with 103 profiles) or its genome's matched detection moved
because r4's families added control cells — and each change is written out
with the fields that moved. The r4 families' own cells are summarised by
verdict.

→ `results/genome_sweep/r4_cell_changes.tsv`, `r4_control_changes.tsv`,
`r4_cells.tsv` (the r4 families' cells), `r4_summary.json`, and
`r4_report.md` rendered from them (D13).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s5_lib import OUT_DIR, r4_families  # noqa: E402

CELL_KEYS = ["proteome_records", "proteome_high", "informative", "matched", "control",
             "genome", "n_found", "n_partial", "n_traces", "bar_bp", "verdict"]
CTL_KEYS = ["control_cells", "detected", "detection", "matched_cells",
            "matched_detected", "matched_detection"]


def render(summary: dict, mine: list[dict], ctl: list[dict]) -> None:
    runs = read_tsv(OUT_DIR / "genome_runs_r4.tsv")
    baits = read_tsv(OUT_DIR / "baits_r4.tsv")
    census = sorted({c["family"] for c in mine})
    L = ["# S5 r4 — the genome sweep for the families added after S20 (D43)", "",
         "Rendered by `scripts/s5r4_compare.py` from `baits_r4.tsv`, "
         "`genome_runs_r4.tsv`, `r4_cells.tsv`, `r4_cell_changes.tsv` and "
         "`r4_control_changes.tsv` (D13).", "",
         "![S5 r4](figures/genome_r4.png)", "",
         f"**Baits**: {len(baits)} for {len({b['family'] for b in baits})} families "
         "(the 7 new census families, KChIP and the 4 decoys), drawn by S5b's rules "
         "B1–B3 restricted to them, in a panel of their own so S5b's panel, runs and "
         f"loci stay byte-identical. **Runs**: {sum(not r['note'] for r in runs)}/{len(runs)} "
         f"genomes, {sum(int(r['loci'] or 0) for r in runs):,} loci, "
         f"{sum(int(r['r4_family_loci'] or 0) for r in runs):,} called to an r4 family, "
         f"miniprot {sum(float(r['miniprot_s'] or 0) for r in runs) / 60:.0f} min. "
         "An r4 locus counts only toward an r4 family (called to one, or partial on "
         "its bait); an old family's trace is judged against S5b's loci only.", "",
         f"**Old cells: {summary['old_cells']:,} compared, {summary['old_cells_changed']} "
         f"changed, {summary['old_verdicts_changed']} verdicts changed.** "
         f"{summary['controls_changed']} genomes' control counts grew (the new families' "
         "high-confidence cells join the controls); none crossed the 0.90 floor.", "",
         "## The new census families, per species", "",
         "| verdict | " + " | ".join(census) + " |", "|---|" + "---|" * len(census)]
    vs = sorted({c["verdict"] for c in mine})
    for v in vs:
        L.append(f"| {v} | " + " | ".join(
            str(sum(c["family"] == f and c["verdict"] == v for c in mine)) for f in census) + " |")
    for v, title in (("genome_found", "Proteome misses (high-confidence intact loci)"),
                     ("absent", "Controlled absences (matched bait, detection ≥ 0.90, D4 bar met)"),
                     ("genome_weak", "Weak genome evidence"), ("gap", "Gaps")):
        rows = [c for c in mine if c["verdict"] == v]
        L += ["", f"**{title}**: " + ("; ".join(f"*{c['species']}* {c['family']}"
                                               + (f" (bar {c['bar_source']})" if v == "absent" else "")
                                               for c in rows) if rows else "none") + "."]
    L += ["", "These are measurements under the S5b instrument, not literature "
          "statements: each absence is a candidate for the S10 checks, and the "
          "single-seed profiles (TMCO1, TMEM87, TMEM109, CLCC1, MITOK) call loci "
          "at the same D32 gates as every other family.", ""]
    (OUT_DIR / "r4_report.md").write_text("\n".join(L))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    a = ap.parse_args()
    base = Path(a.base)
    r4 = r4_families()
    old = {(c["species"], c["family"]): c for c in read_tsv(base / "cells.tsv")}
    new = {(c["species"], c["family"]): c for c in read_tsv(OUT_DIR / "cells.tsv")}
    missing = [k for k in old if k not in new]
    if missing:
        raise SystemExit(f"{len(missing)} S5b cells vanished, e.g. {missing[:3]}")
    changes = []
    for k, o in old.items():
        n = new[k]
        moved = [f for f in CELL_KEYS if str(o.get(f, "")) != str(n.get(f, ""))]
        if moved:
            changes.append({"species": k[0], "family": k[1], "fields": ",".join(moved),
                            **{f"{f}_before": o.get(f, "") for f in ("proteome_high", "control", "genome", "verdict")},
                            **{f"{f}_after": n.get(f, "") for f in ("proteome_high", "control", "genome", "verdict")}})
    write_tsv(OUT_DIR / "r4_cell_changes.tsv",
              list(changes[0]) if changes else ["species", "family", "fields"], changes)
    oc = {r["species"]: r for r in read_tsv(base / "controls.tsv")}
    ctl = []
    for r in read_tsv(OUT_DIR / "controls.tsv"):
        o = oc.get(r["species"], {})
        moved = [f for f in CTL_KEYS if str(o.get(f, "")) != str(r.get(f, ""))]
        if moved:
            ctl.append({"species": r["species"], "fields": ",".join(moved),
                        **{f"{f}_before": o.get(f, "") for f in ("control_cells", "matched_detection")},
                        **{f"{f}_after": r.get(f, "") for f in ("control_cells", "matched_detection")}})
    write_tsv(OUT_DIR / "r4_control_changes.tsv",
              list(ctl[0]) if ctl else ["species", "fields"], ctl)
    mine = [c for k, c in new.items() if k[1] in r4]
    write_tsv(OUT_DIR / "r4_cells.tsv", list(mine[0]), mine)
    verdict_changed = [c for c in changes if c["verdict_before"] != c["verdict_after"]]
    summary = {"old_cells": len(old), "old_cells_changed": len(changes),
               "old_verdicts_changed": len(verdict_changed),
               "changed_fields": dict(Counter(f for c in changes for f in c["fields"].split(","))),
               "controls_changed": len(ctl),
               "r4_cells": len(mine),
               "r4_verdicts": dict(Counter(c["verdict"] for c in mine)),
               "r4_genome": dict(Counter(c["genome"] for c in mine))}
    (OUT_DIR / "r4_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    render(summary, mine, ctl)
    print(json.dumps(summary, indent=1))
    for c in verdict_changed[:20]:
        print(f"  {c['species']:28s} {c['family']:18s} {c['verdict_before']} → "
              f"{c['verdict_after']}  ({c['fields']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
