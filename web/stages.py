"""The roadmap ledger and the ledger task(s) each results file belongs to.

Task states are the repo's own status markers from PUBLICATION_ROADMAP.md (`completed <date>`, `in_progress`,
`pending`), read with the dashboard's own reader (`scripts/dashboard.py:parse_ledger`), so the site and the session
dashboard can never disagree about what is finished.
"""
from __future__ import annotations

import re
import sys
from functools import lru_cache

from . import provenance

sys.path.insert(0, str(provenance.ROOT / "scripts"))
import dashboard  # noqa: E402  (stdlib-only; its FIGURES list comes from figure_notes.py)

LEDGER = "PUBLICATION_ROADMAP.md"

#: (regex on the path below results/, ledger tasks). First match wins; the mapping is taken from INTERFACE.md's
#: description of which stage writes which file.
FILE_TASKS: list[tuple[str, list[str]]] = [
    (r"^s0_baseline/outgroups", ["S7c"]),
    (r"^s0_baseline/", ["S0"]),
    (r"^toolchain_manifest", ["S1"]),
    (r"^benchmark_controls/", ["S1"]),
    (r"^census_v2/s2b_", ["S2b"]),
    (r"^census_v2/s2c_", ["S2c"]),
    (r"^census_v2/r4_", ["S2d"]),
    (r"^census_v2/r5_", ["S2e"]),
    (r"^census_v2/", ["S2"]),
    (r"^census_v3/(r4_|figures/census_r4)", ["S2d"]),
    (r"^census_v3/r5_", ["S2e"]),
    (r"^census_v3/r6_", ["S2f"]),
    (r"^census_v3/", ["S3a"]),
    (r"^proteome_scope/", ["S4"]),
    (r"^panel_sweep/", ["S3b"]),
    (r"^panel_density/(options|instrument|proteomes|summary)", ["S4b"]),
    (r"^panel_density/", ["S4b"]),
    (r"^genome_sweep/(r4_|baits_r4|genome_runs_r4|figures/genome_r4)", ["S5c"]),
    (r"^genome_sweep/", ["S5a", "S5b"]),
    (r"^alignments/", ["S6"]),
    (r"^phylogeny/gost/", ["S2f"]),
    (r"^phylogeny/(tier1_trim|figures/tier1_trim)", ["S7a"]),
    (r"^phylogeny/(reroot_basal|tier1_reroot|tier1/reroot/)", ["S7c", "S7d"]),
    (r"^phylogeny/(tier1|figures/tier1_trees)", ["S7b"]),
    (r"^phylogeny/(fold_network|figures/fold_network)", ["S8a"]),
    (r"^phylogeny/(tier2_(trees|families|placement|report)|tier2/|figures/tier2_trees)", ["S8b"]),
    (r"^phylogeny/", ["S8a"]),
    (r"^filter_atlas/", ["S9"]),
    (r"^repertoire/(contig|daphnia|k2p|ortholog|absent|q4|s10b|figures/s10b)", ["S10b"]),
    (r"^repertoire/", ["S10a"]),
    (r"^duplication/roots_", ["S11c"]),
    (r"^duplication/(repeat|s11b|figures/s11b)", ["S11b"]),
    (r"^duplication/", ["S11a"]),
    (r"^structures/", ["S12"]),
    (r"^selection/", ["S13a", "S13b"]),
    (r"^method_contribution/", ["S15"]),
    (r"^auxiliary/", ["S20"]),
]


@lru_cache(maxsize=1)
def ledger() -> list[dict]:
    return dashboard.parse_ledger()


def task_index() -> dict[str, dict]:
    return {t["id"]: t for t in ledger()}


def badge(task_id: str) -> dict:
    """{id, state, date, title} for one ledger task; state 'unknown' if the ledger has no such row."""
    t = task_index().get(task_id)
    if t is None:
        return {"id": task_id, "state": "unknown", "date": "", "title": ""}
    return {"id": t["id"], "state": t["state"], "date": t["date"], "title": t["title"]}


def tasks_for(path: str) -> list[str]:
    rel = re.sub(r"^results/", "", path)
    for pat, tasks in FILE_TASKS:
        if re.search(pat, rel):
            return tasks
    return []


def stage_of(path: str) -> list[dict]:
    """Ledger badges for the task(s) that wrote a results file."""
    return [badge(t) for t in tasks_for(path)]


def ledger_payload() -> dict:
    tasks = ledger()
    counts = {s: sum(t["state"] == s for t in tasks) for s in ("completed", "in_progress", "pending")}
    return {"tasks": tasks, "counts": counts, "total": len(tasks), "provenance": provenance.of(LEDGER),
            "note": "Task states and headline text are quoted from the roadmap ledger (PUBLICATION_ROADMAP.md). "
                    "The headline column is the project's own prose; the numbers on the other pages are read "
                    "from the committed tables."}
