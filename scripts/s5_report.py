"""s5_report.py — render results/genome_sweep/report.md from its tables (D13).

    python3 scripts/s5_report.py

At panel scale (52 genomes, 3,536 cells) the report carries summaries and
the cells a reader must see — proteome misses, weak loci, absences — and
leaves every other cell to `cells.tsv`.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv  # noqa: E402
from s5_lib import OUT_DIR  # noqa: E402
from s5_verdict import CONTROL_FLOOR, GROUP_BAR_MIN_GENES  # noqa: E402

#: Absences the literature reports — an external check on the verdicts,
#: never an input to them. References pending verification through the
#: review's reference pipeline before any is cited.
EXPECTED_ABSENT = [
    ("Mus musculus", "zac", "ZAC has no rodent orthologue"),
    ("Rattus norvegicus", "zac", "ZAC has no rodent orthologue"),
    ("Caenorhabditis elegans", "nav", "nematodes lack voltage-gated Na+ channels"),
    ("Caenorhabditis elegans", "p2x", "no P2X receptor in C. elegans"),
    ("Drosophila melanogaster", "p2x", "no P2X receptor in Drosophila"),
]
ABSENT_LIKE = ("absent", "absent_below_bar", "absent_bar_unmeasured", "uncontrolled")


def _t(rows: list[dict], cols: list[str], heads: list[str] | None = None) -> str:
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(heads)]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def _n(rows: list[dict], key: str) -> int:
    return sum(int(r.get(key) or 0) for r in rows)


def _short(species: str) -> str:
    g, _, s = species.partition(" ")
    return f"{g[0]}. {s}" if s else g


def main() -> int:
    baits = read_tsv(OUT_DIR / "baits.tsv")
    runs = read_tsv(OUT_DIR / "genome_runs.tsv")
    ctl = read_tsv(OUT_DIR / "controls.tsv")
    cells = read_tsv(OUT_DIR / "cells.tsv")
    traces = read_tsv(OUT_DIR / "traces.tsv")
    spans = read_tsv(OUT_DIR / "spans.tsv")
    introns = read_tsv(OUT_DIR / "introns.tsv")
    ann = read_tsv(OUT_DIR / "annotated_loci.tsv")
    check = read_tsv(OUT_DIR / "genome_found_check.tsv")
    v4 = json.loads((OUT_DIR / "census_v4.json").read_text())

    zero = [c for c in cells if c["proteome_records"] == "0"]
    inf_zero = [c for c in zero if c["informative"] == "1"]
    verdicts, inf_v = Counter(c["verdict"] for c in zero), Counter(c["verdict"] for c in inf_zero)
    failed = [r for r in runs if r.get("note", "").startswith("FAILED")]
    summary = {
        "baits": len(baits), "bait_families": len({b["family"] for b in baits}),
        "bait_species": len({b["species"] for b in baits}),
        "genomes": len(runs), "genomes_failed": len(failed),
        "genome_bp": _n(runs, "total_bp"), "loci": _n(runs, "loci"),
        "family_called_loci": _n(runs, "family_called_loci"),
        "miniprot_cpu_h": round(sum(float(r["miniprot_s"] or 0) for r in runs) / 3600, 2),
        "cells": len(cells),
        "control_cells": _n(ctl, "control_cells"), "control_found": _n(ctl, "found"),
        "control_detected": _n(ctl, "detected"),
        "matched_control_cells": _n(ctl, "matched_cells"),
        "matched_control_detected": _n(ctl, "matched_detected"),
        "unmatched_control_cells": _n(ctl, "unmatched_cells"),
        "unmatched_control_detected": _n(ctl, "unmatched_detected"),
        "genomes_below_floor": sorted(r["species"] for r in ctl if r["matched_detection"]
                                      and float(r["matched_detection"]) < CONTROL_FLOOR),
        "genomes_no_matched_control": sorted(r["species"] for r in ctl
                                             if not r["matched_detection"]),
        "zero_cells": len(zero), "zero_verdicts": dict(verdicts),
        "informative_zero_cells": len(inf_zero), "informative_verdicts": dict(inf_v),
        "traces": len(traces), "families_with_span": len(spans),
        "annotated_loci": len(ann), "control_floor": CONTROL_FLOOR, "census_v4": v4}
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    g_of = {r["species"]: r["group"] for r in runs}
    groups = list(dict.fromkeys(r["group"] for r in runs))
    vcols = [v for v, _ in inf_v.most_common()]
    by_group = []
    for g in groups:
        cg = Counter(c["verdict"] for c in inf_zero if c["group"] == g)
        if cg:
            by_group.append({"group": g, **{v: cg.get(v, "") for v in vcols},
                             "total": sum(cg.values())})
    by_fam = defaultdict(list)
    for c in zero:
        if c["verdict"] in ABSENT_LIKE:
            by_fam[(c["family"], c["verdict"])].append(c)
    absent_rows = [{"family": f, "verdict": v, "n": len(cs),
                    "species": ", ".join(_short(c["species"]) for c in cs),
                    "bar": ", ".join(sorted({c["bar_source"] for c in cs}))}
                   for (f, v), cs in sorted(by_fam.items())]
    key = {(c["species"], c["family"]): c for c in cells}
    exp_rows = [{"species": sp, "family": f, "why": why, **{k: key[(sp, f)][k] for k in
                 ("proteome_records", "genome", "n_traces", "verdict")}}
                for sp, f, why in EXPECTED_ABSENT if (sp, f) in key]
    found = [c for c in zero if c["verdict"] == "genome_found"]
    weak = [c for c in zero if c["verdict"] == "genome_weak"]
    gp = defaultdict(list)
    for c in zero:
        if c["verdict"] == "genome_present":
            gp[c["species"]].append(c["family"])
    g_used = {r["species"]: int(r["max_intron"]) for r in runs if r.get("max_intron")}
    over = [{**a, "max_intron_used": g_used[a["species"]]} for a in ann
            if a["species"] in g_used and int(a["annotated_max_intron"]) > g_used[a["species"]]]

    md = [
        "# S5 — the genomic sweep: what 52 genomes say about the proteome census",
        "",
        "Rendered by `scripts/s5_report.py` from the tables in this directory (D13). "
        "S5a built and measured the instrument on seven genomes; S5b ran it over the "
        "whole S4 panel. Bulk output (genomes, miniprot GFFs, per-genome loci, tblastn, "
        "census v4) is under `<data root>/genomes/`.",
        "",
        "## What was measured",
        "",
        f"- **Bait panel**: {summary['baits']} baits for {summary['bait_families']} catalogue "
        f"families (controls included) from {summary['bait_species']} panel species, one per "
        "species (B1–B3, D37).",
        f"- **Sweep**: {summary['genomes']} genomes, {summary['genome_bp'] / 1e9:.1f} Gbp, "
        f"{len(failed)} failures; {summary['loci']:,} loci, **{summary['family_called_loci']:,} "
        f"called to a family by the S3a profiles** (D32); miniprot "
        f"{summary['miniprot_cpu_h']} h. `-G` 1 Mb for genomes ≥ 1 Gbp, measured not assumed (D38).",
        f"- **Positive control**: {summary['control_found']:,} of {summary['control_cells']:,} "
        f"control cells *found* (a profile-called locus), {summary['control_detected']:,} "
        "detected by the whole instrument, with the genome's own species' baits excluded.",
        f"- **Matched detection — the number every absence inherits: "
        f"{summary['matched_control_detected']:,} / {summary['matched_control_cells']:,}** "
        f"control cells with an in-group non-self bait; unmatched "
        f"{summary['unmatched_control_detected']} / {summary['unmatched_control_cells']}. "
        f"Below the {CONTROL_FLOOR:.0%} floor: {', '.join(summary['genomes_below_floor']) or 'none'}. "
        "No matched control at all (single-species groups — no absence readable): "
        f"{', '.join(summary['genomes_no_matched_control']) or 'none'}.",
        f"- **Census v4**: census v3's {v4['census_v3_rows']:,} proteome rows unchanged + "
        f"**{v4['genome_rows']} genome loci** ({v4['genome_only_rows']} from the genome-only "
        f"species, {v4['genome_found_rows']} proteome misses) in {v4['genome_cells_added']} cells "
        "(`s5_census_v4.py`; SHA-256 in `census_v4.json`). A genome row is a locus, never "
        "merged into a proteome call.",
        "",
        "## Per-genome control",
        "",
        _t(ctl, ["species", "group", "control_cells", "found", "detected", "matched_cells",
                 "matched_detected", "matched_detection", "matched_undetected",
                 "unmatched_cells", "unmatched_detected"],
           ["species", "group", "control", "found", "detected", "matched", "m. detected",
            "m. rate", "m. undetected", "unmatched", "u. detected"]),
        "",
        "## Zero cells: what the genomes say about S3b's absences",
        "",
        f"{len(zero):,} zero cells (census family × species with no proteome call); "
        f"{len(inf_zero)} informative (the family is present at high confidence in another "
        "species of the group). Every cell is in `cells.tsv`.",
        "",
        _t([{"verdict": k, "all": v, "informative": inf_v.get(k, 0)}
            for k, v in verdicts.most_common()], ["verdict", "all", "informative"]),
        "",
        "Informative zero cells by group:",
        "",
        _t(by_group, ["group"] + vcols + ["total"]),
        "",
        "### Proteome misses: a high-confidence intact gene the proteome has no call for",
        "",
        "`genome_found` (D37 (4)). The check beside each is reported, not applied "
        "(`genome_found_check.tsv`, added after the cells were read): the genome's loci of the "
        "family, and proteome entries whose runner-up profile is that family — a gene the "
        "proteome might hold under a sister family's call.",
        "",
        _t(check, ["species", "family", "genome_loci", "genome_loci_called",
                   "proteome_near", "near_min_margin"]),
        "",
        "### Weak loci (not claims)",
        "",
        "A profile-called locus without a high-confidence intact frame: retrocopies, "
        "frameshifted or chained models. Never read as a proteome miss.",
        "",
        _t(weak, ["species", "family", "informative", "n_found", "best_confidence",
                  "n_strong", "found_loci"]),
        "",
        "### Absences (and blocked absences), by family",
        "",
        f"`absent` = no locus, no trace, a matched bait, a controlled genome and N50 ≥ the D4 "
        f"bar; the bar is `cds` (prokaryote/virus: 3 × band), `group` (the group's own median "
        f"span, ≥ {GROUP_BAR_MIN_GENES} genes annotated) or `pooled` (D38).",
        "",
        _t(absent_rows, ["family", "verdict", "n", "species", "bar"]),
        "",
        "### Literature-expected absences (external check, references pending)",
        "",
        _t(exp_rows, ["species", "family", "why", "proteome_records", "genome", "n_traces",
                      "verdict"]),
        "",
        "### Genome-only species",
        "",
    ]
    for sp, fams in sorted(gp.items()):
        md.append(f"- *{sp}*: **{len(fams)} census families present** by a profile-called "
                  f"locus — {', '.join(sorted(fams))}.")
    md += [
        "",
        "## Calibration from annotation (not from miniprot)",
        "",
        "Widest annotated intron at a high-confidence family-called locus, against the `-G` "
        "used (`over_G` = annotated genes with an intron wider than it):",
        "",
        _t(introns, ["species", "annotation_accession", "called_loci", "annotated_loci",
                     "widest_annotated_intron", "widest_gene", "max_intron_used", "over_G"]),
        "",
        f"Annotated genes with an intron beyond `-G` ({len(over)}) — each still called, "
        "which is what D38 predicts (a split gene, not a lost one):",
        "",
        _t(over, ["species", "family", "gene_name", "gene_span", "annotated_max_intron",
                  "max_intron_used", "locus"]),
        "",
        f"Gene spans ({len(ann):,} annotated loci) → D4's bar (D38; {len(spans)} families "
        "measured):",
        "",
        _t(spans, ["family", "n_genes", "n_species", "median_span", "min_span", "max_span",
                   "by_group", "by_group_n"]),
        "",
        "## Runs",
        "",
        _t(runs, ["species", "assembly", "total_bp", "n50", "max_intron", "chunks",
                  "alignments", "self_dropped", "loci", "family_called_loci", "miniprot_s",
                  "call_s", "note"]),
        "",
    ]
    (OUT_DIR / "report.md").write_text("\n".join(md))
    print(json.dumps({k: v for k, v in summary.items() if k != "census_v4"}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
