"""s5_report.py — render results/genome_sweep/report.md from its tables (D13).

    python3 scripts/s5_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv  # noqa: E402
from s5_lib import OUT_DIR  # noqa: E402
from s5_verdict import CONTROL_FLOOR  # noqa: E402

#: Absences the literature reports for pilot genomes — an external check on
#: the verdicts, never an input to them. References are pending verification
#: through the review's reference pipeline before any is cited.
EXPECTED_ABSENT = [
    ("Mus musculus", "zac", "ZAC has no rodent orthologue"),
    ("Caenorhabditis elegans", "nav", "nematodes lack voltage-gated Na+ channels"),
    ("Caenorhabditis elegans", "p2x", "no P2X receptor in C. elegans"),
    ("Drosophila melanogaster", "p2x", "no P2X receptor in Drosophila"),
]


def _t(rows: list[dict], cols: list[str], heads: list[str] | None = None) -> str:
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(heads)]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def main() -> int:
    baits = read_tsv(OUT_DIR / "baits.tsv")
    runs = read_tsv(OUT_DIR / "genome_runs.tsv")
    ctl = read_tsv(OUT_DIR / "controls.tsv")
    cells = read_tsv(OUT_DIR / "cells.tsv")
    traces = read_tsv(OUT_DIR / "traces.tsv")
    spans = read_tsv(OUT_DIR / "spans.tsv")
    introns = read_tsv(OUT_DIR / "introns.tsv")
    ann = read_tsv(OUT_DIR / "annotated_loci.tsv")

    zero = [c for c in cells if c["proteome_records"] == "0"]
    verdicts = Counter(c["verdict"] for c in zero)
    inf_zero = [c for c in zero if c["informative"] == "1"]
    inf_v = Counter(c["verdict"] for c in inf_zero)
    m_cells = sum(int(r["matched_cells"] or 0) for r in ctl)
    m_det = sum(int(r["matched_detected"] or 0) for r in ctl)
    u_cells = sum(int(r["unmatched_cells"] or 0) for r in ctl)
    u_det = sum(int(r["unmatched_detected"] or 0) for r in ctl)
    c_cells = sum(int(r["control_cells"] or 0) for r in ctl)
    c_found = sum(int(r["found"] or 0) for r in ctl)
    summary = {
        "baits": len(baits), "bait_families": len({b["family"] for b in baits}),
        "bait_species": len({b["species"] for b in baits}),
        "genomes": len(runs), "genome_bp": sum(int(r["total_bp"] or 0) for r in runs),
        "loci": sum(int(r["loci"] or 0) for r in runs),
        "family_called_loci": sum(int(r["family_called_loci"] or 0) for r in runs),
        "control_cells": c_cells, "control_found": c_found,
        "matched_control_cells": m_cells, "matched_control_detected": m_det,
        "unmatched_control_cells": u_cells, "unmatched_control_detected": u_det,
        "zero_cells": len(zero), "zero_verdicts": dict(verdicts),
        "informative_zero_cells": len(inf_zero), "informative_verdicts": dict(inf_v),
        "traces": len(traces), "families_with_span": len(spans),
        "control_floor": CONTROL_FLOOR}
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    by_key = {(c["species"], c["family"]): c for c in cells}
    exp_rows = []
    for sp, fam, why in EXPECTED_ABSENT:
        c = by_key.get((sp, fam))
        if c:
            exp_rows.append({"species": sp, "family": fam, "why": why,
                             "proteome": c["proteome_records"],
                             "genome": c["genome"], "traces": c["n_traces"],
                             "verdict": c["verdict"]})
    gf = [c for c in zero if c["verdict"] in ("genome_found", "genome_weak",
                                               "partial", "gap", "trace")]
    gp = [c for c in zero if c["verdict"] == "genome_present"]
    absent = [c for c in zero if c["verdict"].startswith("absent") or c["verdict"] == "uncontrolled"]

    md = [
        "# S5a — the genomic sweep's instrument, and a seven-genome pilot",
        "",
        "Rendered by `scripts/s5_report.py` from the tables in this directory "
        "(D13). Bulk output (genomes, miniprot GFFs, per-genome loci, tblastn) "
        "is under `<data root>/genomes/`.",
        "",
        "## What was measured",
        "",
        f"- **Bait panel**: {summary['baits']} baits for {summary['bait_families']} "
        f"catalogue families (controls included) from {summary['bait_species']} "
        "panel species — rules B1–B3 in `s5_baits.py`, one bait per species (D37).",
        f"- **Pilot**: {summary['genomes']} genomes, "
        f"{summary['genome_bp'] / 1e9:.2f} Gbp; {summary['loci']} loci, "
        f"{summary['family_called_loci']} called to a family by the S3a profiles (D32).",
        f"- **Positive control, the whole instrument**: of {c_cells} control cells "
        f"(families the genome's own proteome carries at high confidence, or a "
        f"genome-only species' group-core families), {c_found} are *found* "
        "(a profile-called locus) with the genome's own species' baits excluded.",
        f"- **Matched detection** — the number an absence inherits: **{m_det} / "
        f"{m_cells}** control cells with a bait from another species of the same "
        f"group are detected (found, partial, gap or trace); **{u_det} / {u_cells}** "
        "without one. Detection is a property of the nearest bait.",
        "",
        "## Per-genome control",
        "",
        _t(ctl, ["species", "group", "control_cells", "found", "detected",
                 "matched_cells", "matched_detected", "matched_detection",
                 "matched_undetected", "unmatched_cells", "unmatched_detected"],
           ["species", "group", "control", "found", "detected", "matched",
            "m. detected", "m. rate", "m. undetected", "unmatched", "u. detected"]),
        "",
        f"The control floor is {CONTROL_FLOOR:.0%} matched detection; a genome below "
        "it has every absence read as `uncontrolled`.",
        "",
        "## Zero cells: what the genomes say about S3b's absences",
        "",
        f"{len(zero)} zero cells in the pilot genomes; {len(inf_zero)} informative "
        "(the family is present at high confidence in another species of the group).",
        "",
        _t([{"verdict": k, "all": v, "informative": inf_v.get(k, 0)}
            for k, v in verdicts.most_common()], ["verdict", "all", "informative"]),
        "",
        "### Zero cells with genomic evidence",
        "",
        _t([c for c in gf if c["species"] not in {x["species"] for x in gp}],
           ["species", "family", "informative", "proteome_band", "genome",
                "n_found", "best_confidence", "n_strong", "n_partial", "n_traces",
                "trace_best_bits", "verdict"]),
        "",
        f"Genome-only species (no proteome to compare): {len(gp)} families present "
        "by a profile-called locus — "
        + ", ".join(sorted(f"{c['species'].split()[0]} {c['family']}" for c in gp)) + ".",
        "",
        "### Zero cells read as absent (or blocked)",
        "",
        _t(absent, ["species", "family", "matched", "verdict"]),
        "",
        "### Literature-expected absences (external check, references pending)",
        "",
        _t(exp_rows, ["species", "family", "why", "proteome", "genome",
                      "traces", "verdict"]),
        "",
        "## Calibration from annotation (not from miniprot)",
        "",
        "Widest annotated intron at a high-confidence family-called locus, against "
        "the `-G` used:",
        "",
        _t(introns, ["species", "annotation_accession", "called_loci",
                     "annotated_loci", "widest_annotated_intron", "widest_gene",
                     "max_intron_used", "over_G"]),
        "",
        f"Gene spans ({len(ann)} annotated loci) → D4's bar per family "
        f"(median span; {len(spans)} families measured):",
        "",
        _t(spans, ["family", "n_genes", "n_species", "median_span", "min_span",
                   "max_span", "by_group"]),
        "",
        "## Runs",
        "",
        _t(runs, ["species", "assembly", "total_bp", "n50", "max_intron", "chunks",
                  "alignments", "self_dropped", "loci", "family_called_loci",
                  "miniprot_s", "call_s", "note"]),
        "",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(md))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
