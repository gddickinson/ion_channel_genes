"""Committed Newick treefiles (IQ-TREE output) with the statistics row each stage wrote for them.

Groups: S7b tier-1 family trees (tier1_trees.tsv), S7c/S7d re-rooted trees (tier1_reroot_trees.tsv), S8b tier-2
superfamily trees (tier2_trees.tsv), and two single trees (S2f GOST, S11b repeat units). Internal-node labels are
UFBoot support, as written by IQ-TREE.
"""
from __future__ import annotations

from functools import lru_cache

from . import provenance, stages, tables

P = "results/phylogeny/"
GROUPS = [
    {"id": "tier1", "title": "Tier-1 family trees", "table": P + "tier1_trees.tsv",
     "file": lambda r: f"{P}tier1/{r['family']}.treefile", "name": lambda r: r["family"]},
    {"id": "reroot", "title": "Re-rooted tier-1 trees (two declared outgroups each)",
     "table": P + "tier1_reroot_trees.tsv",
     "file": lambda r: f"{P}tier1/reroot/{r['family']}.{r['root_set']}.treefile",
     "name": lambda r: f"{r['family']} · {r['root_set']}"},
    {"id": "tier2", "title": "Tier-2 superfamily trees", "table": P + "tier2_trees.tsv",
     "file": lambda r: f"{P}tier2/{r['unit']}.treefile", "name": lambda r: r["unit"]},
]
SINGLES = [
    {"group": "other", "name": "GOST superfamily (TMEM87 lineage)", "path": P + "gost/gost.treefile"},
    {"group": "other", "name": "four-repeat channels: whole repeats", "path": "results/duplication/repeats/repeat.treefile"},
    {"group": "other", "name": "four-repeat channels: pore modules", "path": "results/duplication/repeats/pore.treefile"},
]


def _tips(path: str) -> int:
    return (provenance.ROOT / path).read_text().count(",") + 1


@lru_cache(maxsize=1)
def catalog() -> list[dict]:
    out = []
    for g in GROUPS:
        for r in tables.read_records(g["table"]):
            path = g["file"](r)
            if provenance.is_servable(path):
                out.append({"path": path, "group": g["id"], "name": g["name"](r), "tips": _tips(path),
                            "stats": r, "stats_table": g["table"], "stage": stages.stage_of(path)})
    for s in SINGLES:
        if provenance.is_servable(s["path"]):
            out.append({**s, "tips": _tips(s["path"]), "stats": None, "stats_table": None,
                        "stage": stages.stage_of(s["path"])})
    return out


def paths() -> list[str]:
    return [t["path"] for t in catalog()]


def index_payload() -> dict:
    titles = {g["id"]: g["title"] for g in GROUPS} | {"other": "Other committed trees"}
    return {"groups": titles, "trees": [{k: v for k, v in t.items()} for t in catalog()]}


def load(path: str) -> dict | None:
    meta = next((t for t in catalog() if t["path"] == path), None)
    if meta is None:
        return None
    out = {**meta, "newick": (provenance.ROOT / path).read_text().strip(), "provenance": provenance.of(path),
           }
    if meta["stats_table"]:
        out["stats_provenance"] = provenance.of(meta["stats_table"])
    return out
