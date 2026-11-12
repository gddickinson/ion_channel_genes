"""Prose and figures: FINDINGS.md, the stage reports and the headline figures with their captions.

The stage reports (results/**/report*.md) are rendered by each stage purely from its committed tables (D13).
FINDINGS.md is the project's hand-written narrative; its "(pending: task)" markers are counted and highlighted on
the site. Figure captions come from scripts/figure_notes.py, the single source the README and dashboard also use.
"""
from __future__ import annotations

import re
from functools import lru_cache

from . import provenance, stages

FINDINGS = "FINDINGS.md"
PENDING_RE = re.compile(r"\(pending:[^)]*\)")


def findings_payload() -> dict:
    text = (provenance.ROOT / FINDINGS).read_text()
    return {"markdown": text, "pending_markers": PENDING_RE.findall(text), "provenance": provenance.of(FINDINGS),
            "note": "The project's plain-language narrative, written at the end of each task. Claims that depend "
                    "on unfinished work carry the repo's own '(pending: task)' marker, highlighted here."}


@lru_cache(maxsize=1)
def report_paths() -> list[str]:
    return sorted(p for p in provenance.tracked()
                  if p.startswith("results/") and re.search(r"(^|/)(\w+_)?report\.md$", p.rsplit("/", 1)[-1])
                  and provenance.is_servable(p))


def _title(path: str) -> str:
    for line in (provenance.ROOT / path).read_text().splitlines():
        if line.startswith("#"):
            return line.lstrip("#").strip()
    return path


def reports_payload() -> dict:
    return {"reports": [{"path": p, "title": _title(p), "stage": stages.stage_of(p)} for p in report_paths()]}


def report_payload(path: str) -> dict | None:
    if path not in report_paths():
        return None
    return {"path": path, "title": _title(path), "markdown": (provenance.ROOT / path).read_text(),
            "provenance": provenance.of(path), "stage": stages.stage_of(path)}


@lru_cache(maxsize=1)
def figures() -> list[dict]:
    out = []
    for f in stages.dashboard._NOTES:
        if provenance.is_servable(f["path"]):
            out.append({k: f.get(k, "") for k in ("path", "task", "script", "title", "shows", "how", "read")})
    return out


def figure_paths() -> list[str]:
    return [f["path"] for f in figures()]


def figures_payload() -> dict:
    notes = "scripts/figure_notes_late.py" if provenance.is_servable("scripts/figure_notes_late.py") else None
    return {"figures": [{**f, "provenance": provenance.of(f["path"]),
                         "stage": [stages.badge(t) for t in re.findall(r"S\d+[a-z]?", f["task"])]}
                        for f in figures()],
            "captions": [provenance.of("scripts/figure_notes.py")] + ([provenance.of(notes)] if notes else [])}
