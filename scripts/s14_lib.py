"""S14 shared helpers — manuscript assembly, figures, deposit manifest, claims.

Used by `s14_figures.py`, `s14_deposit.py`, `s14_claims.py` and the driver
`s14_assemble.py`. Stdlib only: the manuscript package must build on a bare
Python so a reviewer can regenerate it without the analysis toolchain.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

# ---------------------------------------------------------------- page size
#: Printed width of a full-width figure, in inches. A4 at 2.0 cm margins gives
#: a 17.0 cm text block. `figstyle` imports these so figures are drawn at the
#: size they are placed at.
W_FULL = 6.7
W_HALF = 3.25
#: Tallest a figure may be before the page height limits it.
H_MAX = 8.6

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
MS = ROOT / "manuscript"
MS_FIGS = MS / "figures"

# ---------------------------------------------------------------- manuscript

#: Section files, in the order they are stitched into `manuscript.md`.
#: The manuscript does not exist until S14 writes it; `s14_assemble.py`
#: exits non-zero while any of these is missing, which is the intended
#: behaviour for a package that is not yet written.
SECTION_ORDER = [
    "00_frontmatter.md",
    "01_main.md",
    "02_results_scope.md",       # the declared scope, and what it excludes
    "03_results_classifier.md",  # the benchmark, per family and per hazard
    "04_results_census.md",      # the census across the declared space
    "05_results_forest.md",      # tier-1 and tier-2 phylogenies
    "06_results_network.md",     # the fold network, and what is convergent
    "07_results_repertoire.md",  # repertoire evolution across the tree of life
    "08_discussion.md",
    "09_methods_search.md",
    "10_methods_analysis.md",
    "11_figure_legends.md",
    "12_extended_data.md",
    "13_supplementary.md",
    "14_data_availability.md",
    "15_references.md",
]

SECTION_FRONTMATTER = "00_frontmatter.md"
SECTION_FIGURE_LEGENDS = "11_figure_legends.md"
SECTION_EXTENDED_DATA = "12_extended_data.md"
SECTION_SUPPLEMENTARY = "13_supplementary.md"

# ------------------------------------------------------------------ figures

#: (number, slug, source path relative to results/, one-line caption stub).
#: This is the *plan*: the source paths are where the roadmap's tasks write
#: their figures, and `s14_figures.py` fails loudly on any that is missing
#: rather than quietly shipping an incomplete set. Renumber freely — the
#: number here is the publication number, and stale files in
#: `manuscript/figures/` are deleted on every build.
MAIN_FIGURES = [
    (1, "scope", "s0_baseline/figures/catalogue_scope",
     "What counts as an ion channel: 90 families, 25 superfamilies, and "
     "the 22 catalogued so they can be excluded"),
    (2, "classifier", "benchmark_controls/figures/benchmark",
     "Classification without gene symbols: recall per family, specificity "
     "per hazard, and which tier made each call"),
    (3, "census", "census/figures/census_by_superfamily",
     "The census across the declared search space"),
    (4, "forest", "phylogeny/figures/forest_overview",
     "A forest, not a tree: tier-1 families inside tier-2 pore-module "
     "superfamilies"),
    (5, "filter", "filter_atlas/figures/filter_vs_tree",
     "The selectivity filter against the pore-module phylogeny"),
    (6, "network", "fold_network/figures/fold_network",
     "What relates the superfamilies: structural similarity, with no "
     "branch lengths and no common ancestor implied"),
    (7, "repertoire", "repertoire/figures/repertoire_heatmap",
     "Channel repertoires across the tree of life"),
]

#: Extended Data figures. Several bundle more than one source panel file.
EXTENDED_FIGURES = [
    (1, "signature_sharing", ["s0_baseline/figures/signature_sharing"],
     "Which domain signatures are shared, and which cross the "
     "channel / non-channel boundary"),
    (2, "hazards", ["benchmark_controls/figures/hazard_matrix"],
     "The sixteen hazards, the test that closes each, and the ones no "
     "panel member exercises"),
    (3, "method_contribution", ["methods/figures/recall_curve"],
     "What each search method was worth, per superfamily"),
    (4, "tier1_trees", ["phylogeny/figures/tier1_panel"],
     "Tier-1 phylogenies, one per census family"),
    (5, "tier2_trees", ["phylogeny/figures/tier2_panel"],
     "Tier-2 pore-module phylogenies, one per alignable superfamily"),
    (6, "refusals", ["phylogeny/figures/refused_superfamilies"],
     "The superfamilies no tree can span, and why"),
    (7, "annotation_audit", ["annotation_audit/figures/annotation_audit"],
     "How the channels are recorded across four databases"),
    (8, "auxiliary", ["auxiliary/figures/auxiliary_inflation"],
     "The excluded fifteen per cent: auxiliary subunits by family"),
    (9, "structures", ["structures/figures/afdb_coverage"],
     "Structural coverage per family"),
    (10, "mechanism", ["mechanism/figures/clc_ano_trees"],
     "Mechanism against clade in the CLC and anoctamin families"),
    (11, "clinical", ["clinical/figures/variant_map"],
     "Channelopathy variants on the pore module"),
]

#: Supplementary figures — the alignments and structures the main figures
#: rest on. Nothing in them is re-aligned or re-rendered from new data, so a
#: supplementary panel that disagreed with its main figure would be a bug.
SUPPLEMENTARY_FIGURES = [
    (1, "alignment_tier1", ["alignments/figures/supp_aln_tier1"],
     "The within-family alignments behind the tier-1 trees"),
    (2, "alignment_pore", ["alignments/figures/supp_aln_pore_modules"],
     "The pore modules, residue by residue, with the extraction method "
     "marked per sequence"),
    (3, "filter_alignment", ["alignments/figures/supp_aln_filter"],
     "The four-repeat filter locus, projected from the reference"),
    (4, "reference_panel", ["s0_baseline/figures/reference_panel"],
     "The reference panel: what it covers and what it does not"),
    (5, "structure_folds", ["structures/figures/supp_fold_gallery"],
     "One representative structure per pore fold"),
    (6, "network_evidence", ["fold_network/figures/supp_network_evidence"],
     "Every fold-network edge with its measurement, and every pair not "
     "measured"),
]

# ------------------------------------------------------------------ deposit

#: Directories deposited wholesale (relative to results/).
DEPOSIT_DIRS = [
    "alignments", "annotation_audit", "annotation_validation", "architecture",
    "benchmark_controls", "census_v2", "census_v3", "census_v4", "census_v5",
    "constraint", "duplication", "expression", "figures", "hmm_sweep",
    "loss_dynamics", "methods", "msa_v2", "phylogeny", "reconciliation",
    "s23_baits", "s5_baits", "selection", "structures", "synteny",
]

#: Individual files at the results/ root that are deposited.
DEPOSIT_ROOT_GLOBS = ["*.tsv", "*.csv", "*.md", "*.txt", "*.faa", "*.json"]

#: Extensions never deposited (working files, caches, binaries).
DEPOSIT_SKIP_SUFFIXES = {".pyc", ".log"}

#: Exact filenames never deposited.
DEPOSIT_SKIP_NAMES = {".DS_Store", "Thumbs.db",
                      ".pdf_build.md", ".pdf_preamble.tex"}

#: Bulk data classes deliberately excluded, with how to regenerate each.
#: Filled in by the tasks that create each class — an entry here without a
#: working regeneration command is worse than no entry, so add one only when
#: the command has been run.
BULK_EXCLUSIONS: list[tuple[str, str, str, str, str]] = [
    # (what, size, source, manifest committed in the repo, command)
]


def deposit_group(rel_parts: tuple[str, ...]) -> str:
    """Group label for the deposit manifest: one row group per source tree."""
    if rel_parts[0] == "results" and len(rel_parts) > 1:
        return f"results/{rel_parts[1]}" if len(rel_parts) > 2 else "results"
    return rel_parts[0] if len(rel_parts) > 1 else "root"


def png_size_inches(path: Path) -> tuple[float, float]:
    """(width, height) in inches from a PNG's IHDR and pHYs chunks.

    Stdlib only, so the manuscript package still builds without matplotlib.
    Falls back to 400 dpi (the project's savefig default) when a file carries
    no physical-size chunk.
    """
    import struct

    with path.open("rb") as fh:
        if fh.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"{path} is not a PNG")
        px_w = px_h = 0
        dpi = 400.0
        while True:
            head = fh.read(8)
            if len(head) < 8:
                break
            length, kind = struct.unpack(">I4s", head)
            data = fh.read(length)
            fh.read(4)                                   # CRC
            if kind == b"IHDR":
                px_w, px_h = struct.unpack(">II", data[:8])
            elif kind == b"pHYs" and len(data) >= 9 and data[8] == 1:
                ppm_x = struct.unpack(">I", data[:4])[0]  # pixels per metre
                if ppm_x:
                    dpi = ppm_x * 0.0254
            elif kind in (b"IDAT", b"IEND"):
                break
    if not px_w:
        raise ValueError(f"{path}: no IHDR")
    return px_w / dpi, px_h / dpi


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_json(path: Path):
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def size_str(n_bytes: int) -> str:
    """Compact human size (1 decimal above KB)."""
    x = float(n_bytes)
    for unit in ("B", "KB", "MB", "GB"):
        if x < 1024 or unit == "GB":
            return f"{x:.0f} {unit}" if unit == "B" else f"{x:.1f} {unit}"
        x /= 1024.0
    return f"{x:.1f} GB"
