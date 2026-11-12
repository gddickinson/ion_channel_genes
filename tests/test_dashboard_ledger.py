"""Regression tests for scripts/dashboard.py:parse_ledger (the roadmap ledger reader).

A results cell may contain a literal '|' (bipartitions such as 'CNG 142|171'); the reader must not let that shift the
columns and misread a completed task as pending.

Run: python -m pytest -q tests/test_dashboard_ledger.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import dashboard  # noqa: E402

TEXT = """
| ID | Task (one session each) | Depends | Status | Results (headline) |
|----|-------------------------|---------|--------|--------------------|
| S11a | Reconcile | S7b | completed 2026-10-03 | CNG 142|171 (smaller side), K2P 98|282 |
| S13b | Selection results | S13a | pending | |

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | Method contribution | S5b | high | completed 2026-09-29 | split {I,II}|{III,IV} p = 0.1 |
| S16 | Annotation audit | S5b | medium | pending | |
"""


def by_id(tasks):
    return {t["id"]: t for t in tasks}


def test_pipe_inside_results_cell_keeps_status():
    t = by_id(dashboard.parse_ledger(TEXT))
    assert t["S11a"]["state"] == "completed" and t["S11a"]["date"] == "2026-10-03"
    assert t["S11a"]["results"] == "CNG 142|171 (smaller side), K2P 98|282"
    assert t["S15"]["state"] == "completed"
    assert t["S15"]["results"] == "split {I,II}|{III,IV} p = 0.1"


def test_plain_rows_unchanged():
    t = by_id(dashboard.parse_ledger(TEXT))
    assert t["S13b"]["state"] == "pending" and t["S13b"]["results"] == ""
    assert t["S16"]["state"] == "pending" and t["S16"]["depends"] == "S5b"


def test_real_roadmap_s11_rows_completed():
    t = by_id(dashboard.parse_ledger())
    assert t["S11a"]["state"] == "completed"
    assert t["S11b"]["state"] == "completed"
