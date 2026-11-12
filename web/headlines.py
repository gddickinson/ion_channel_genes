"""Headline numbers for the home page, each read at build time from one committed file.

A headline is declared, never typed in: `file` + `key` (a path into a JSON summary) or `rows` (the data-row count of
a TSV), optionally `of` (a second key, shown as "x of y"). The value shown is whatever the committed file holds now,
with its provenance and the ledger status of the task that wrote it. Tasks the ledger lists as not finished are shown
as pending tiles with no number.
"""
from __future__ import annotations

import json

from . import provenance, stages, tables

HEADLINES: list[dict] = [
    {"id": "families", "file": "results/s0_baseline/verification_summary.json",
     "key": ["catalogue", "families_census"], "of": ["catalogue", "families_total"],
     "label": "census families (of all catalogued families)", "task": "S0"},
    {"id": "human_genes", "file": "results/s0_baseline/verification_summary.json",
     "key": ["catalogue", "human_genes_census"], "label": "human pore-forming census genes", "task": "S0"},
    {"id": "superfamilies", "file": "results/s0_baseline/verification_summary.json",
     "key": ["catalogue", "superfamilies"], "label": "superfamilies in the catalogue", "task": "S0"},
    {"id": "benchmark", "file": "results/benchmark_controls/summary.json",
     "key": ["positives_correct"], "of": ["positives"], "label": "S1 benchmark: positives named correctly",
     "task": "S1"},
    {"id": "decoys", "file": "results/benchmark_controls/summary.json",
     "key": ["decoys_rejected"], "of": ["decoys"], "label": "S1 benchmark: decoys rejected", "task": "S1"},
    {"id": "records", "file": "results/census_v2/release.json", "key": ["r5", "union"],
     "label": "UniProtKB records carrying a pore signature (release 2026_03, revision r5)", "task": "S2e"},
    {"id": "species", "file": "results/proteome_scope/summary.json", "key": ["species"],
     "label": "panel species (the declared search space)", "task": "S4"},
    {"id": "entries", "file": "results/proteome_scope/summary.json", "key": ["entries"],
     "label": "genes in the panel's reference proteomes", "task": "S4"},
    {"id": "genomes", "file": "results/genome_sweep/summary.json", "key": ["genomes"],
     "label": "genomes swept", "task": "S5b"},
    {"id": "loci", "file": "results/genome_sweep/summary.json", "key": ["family_called_loci"], "of": ["loci"],
     "label": "genome loci called to a family (of all loci)", "task": "S5b"},
    {"id": "orders", "file": "results/panel_density/order_panel.tsv", "rows": True,
     "label": "eukaryotic orders in the dense panel", "task": "S4b"},
    {"id": "members", "file": "results/method_contribution/summary.json", "key": ["proteome_members"],
     "label": "final-census members in the panel proteomes", "task": "S15"},
    {"id": "trees", "file": "results/phylogeny/tier1_summary.json", "key": ["trees", "trees"],
     "label": "tier-1 family trees", "task": "S7b"},
    {"id": "dups", "file": "results/duplication/summary.json", "key": ["dups_supported"], "of": ["dups_stated"],
     "label": "supported gene duplications (of all stated)", "task": "S11a"},
    {"id": "afdb", "file": "results/structures/summary.json", "key": ["coverage", "model"],
     "of": ["coverage", "members"], "label": "census members with an AlphaFold DB model", "task": "S12"},
    {"id": "selection", "pending": True, "label": "selection results (codeml / HyPhy)", "task": "S13b"},
]


def _get(obj, key: list[str]):
    for k in key:
        obj = obj[k]
    return obj


def resolve(spec: dict) -> dict:
    """One headline with its value read from the committed file now."""
    out = {k: spec[k] for k in ("id", "label", "task") if k in spec}
    out["stage"] = stages.badge(spec["task"])
    if spec.get("pending"):
        out.update(value=None, pending=True)
        return out
    path = spec["file"]
    if not provenance.is_servable(path):
        raise FileNotFoundError(f"headline {spec['id']}: {path} is not a committed file")
    if spec.get("rows"):
        _cols, rows = tables.read_tsv(path)
        out.update(value=len(rows), source="data rows in the table")
    else:
        data = json.loads((provenance.ROOT / path).read_text())
        out.update(value=_get(data, spec["key"]), source="key " + ".".join(spec["key"]))
        if "of" in spec:
            out["of"] = _get(data, spec["of"])
            out["source"] += " of " + ".".join(spec["of"])
    out["provenance"] = provenance.of(path)
    return out


def payload() -> dict:
    return {"headlines": [resolve(s) for s in HEADLINES], "head": provenance.head()}


def files() -> list[str]:
    return sorted({s["file"] for s in HEADLINES if "file" in s})
