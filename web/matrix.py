"""The S5 presence matrix: every census family in every panel species, with its genome verdict.

Read from results/genome_sweep/family_by_species_v4.tsv (census v4). Verdict labels and their best-evidence-first
order are the repo's own, from scripts/s5r4_figures.py:VERDICT (copied here because that module imports matplotlib;
tests/test_web.py checks the copy still matches). 'unmatched' appears in the table but not in that scale.
"""
from __future__ import annotations

from . import provenance, stages, tables

SOURCE = "results/genome_sweep/family_by_species_v4.tsv"

VERDICTS: list[tuple[str, str]] = [
    ("present", "in the proteome"),
    ("genome_found", "genome only: proteome miss"),
    ("genome_present", "genome-only species"),
    ("genome_weak", "genome, weak"),
    ("partial", "partial / gap"),
    ("gap", "partial / gap"),
    ("trace", "trace"),
    ("absent", "controlled absence"),
    ("no_locus_unrescued", "not informative"),
]


def payload() -> dict:
    recs = tables.read_records(SOURCE)
    fams: dict[str, str] = {}
    species: dict[str, str] = {}
    for r in recs:
        fams.setdefault(r["family"], r["superfamily"])
        species.setdefault(r["species"], r["group"])
    known = [v for v, _ in VERDICTS]
    extra = sorted({r["verdict"] for r in recs} - set(known))
    order = known + extra
    vidx = {v: i for i, v in enumerate(order)}
    fi = {f: i for i, f in enumerate(fams)}
    si = {s: i for i, s in enumerate(species)}
    cells = [[None] * len(fams) for _ in species]
    for r in recs:
        cells[si[r["species"]]][fi[r["family"]]] = [vidx[r["verdict"]], int(r["proteome_records"]),
                                                    int(r["proteome_high"]), int(r["genome_loci_added"])]
    labels = dict(VERDICTS)
    return {"families": [{"id": f, "superfamily": sf} for f, sf in fams.items()],
            "species": [{"name": s, "group": g} for s, g in species.items()],
            "verdicts": [{"id": v, "label": labels.get(v, v), "in_repo_scale": v in labels} for v in order],
            "cell_fields": ["verdict", "proteome_records", "proteome_high", "genome_loci_added"],
            "cells": cells, "provenance": provenance.of(SOURCE), "stage": stages.stage_of(SOURCE)}
