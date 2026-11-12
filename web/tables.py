"""The committed results tables, read at request time (and so at static-build time).

Only files git tracks under results/ are listed. TSVs are served as columns + rows of strings; small JSON summaries
are flattened to key/value rows. Large tables are capped at MAX_ROWS rows in file order, and the page says so and
links the full file at its commit.
"""
from __future__ import annotations

import csv
import json
from functools import lru_cache

from . import provenance, stages

MAX_ROWS = 1500
#: JSON files larger than this are bulk records (sweep databases, matrices, sequences), not summaries.
MAX_JSON_BYTES = 20_000


def _kind(path: str) -> str | None:
    if path.endswith(".tsv"):
        return "tsv"
    if path.endswith(".json") and (provenance.ROOT / path).stat().st_size <= MAX_JSON_BYTES:
        return "json"
    return None


@lru_cache(maxsize=1)
def catalog() -> list[dict]:
    out = []
    for path in sorted(provenance.tracked()):
        if not path.startswith("results/") or "/figures/" in path or not provenance.is_servable(path):
            continue
        kind = _kind(path)
        if kind is None:
            continue
        p = provenance.ROOT / path
        if kind == "tsv":
            with p.open(newline="") as fh:
                n = max(sum(1 for _ in fh) - 1, 0)
        else:
            n = len(flatten(json.loads(p.read_text())))
        rel = path[len("results/"):]
        out.append({"path": path, "dir": rel.rsplit("/", 1)[0] if "/" in rel else "", "name": rel.rsplit("/", 1)[-1],
                    "kind": kind, "rows": n, "bytes": p.stat().st_size, "tasks": stages.tasks_for(path)})
    return out


def paths() -> list[str]:
    return [t["path"] for t in catalog()]


def flatten(obj, prefix: str = "") -> list[list[str]]:
    """Nested JSON -> [[dotted.key, value], ...]; lists of scalars joined, lists of objects indexed."""
    rows: list[list[str]] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            rows += flatten(v, f"{prefix}{k}.")
    elif isinstance(obj, list) and any(isinstance(v, (dict, list)) for v in obj):
        for i, v in enumerate(obj):
            rows += flatten(v, f"{prefix}{i}.")
    elif isinstance(obj, list):
        rows.append([prefix[:-1], ", ".join(map(str, obj))])
    else:
        rows.append([prefix[:-1], "" if obj is None else str(obj)])
    return rows


def read_tsv(path: str) -> tuple[list[str], list[list[str]]]:
    with (provenance.ROOT / path).open(newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return (rows[0] if rows else []), rows[1:]


def read_records(path: str) -> list[dict]:
    cols, rows = read_tsv(path)
    return [dict(zip(cols, r)) for r in rows]


def load(path: str, max_rows: int = MAX_ROWS) -> dict | None:
    """One table for the explorer, or None if it is not a committed results table."""
    meta = next((t for t in catalog() if t["path"] == path), None)
    if meta is None:
        return None
    if meta["kind"] == "tsv":
        cols, rows = read_tsv(path)
    else:
        cols, rows = ["key", "value"], flatten(json.loads((provenance.ROOT / path).read_text()))
    return {**meta, "columns": cols, "rows": rows[:max_rows], "total_rows": len(rows),
            "truncated": len(rows) > max_rows, "provenance": provenance.of(path), "stage": stages.stage_of(path)}


def index_payload() -> dict:
    return {"tables": catalog(), "max_rows": MAX_ROWS, "max_json_bytes": MAX_JSON_BYTES}
