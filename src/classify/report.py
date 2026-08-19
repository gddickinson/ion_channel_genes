"""Benchmark tables for the classifier — built here, rendered by the S1 scripts.

Decision **D13**, inherited: a report is rendered from committed tables and
never hand-written alongside them, so the prose and the data cannot drift.
This module owns the tables; `scripts/s1_report.py` owns the prose that
reads them back.

Four tables:

`calls_table`      one row per classified protein, with the decisive tier
                   and the full evidence string.
`confusion_table`  expected family × called family, for a panel whose truth
                   is known (the catalogue's own gene lists).
`recall_table`     per-family recall, plus which tier made each call — the
                   number that says whether the architecture tier is doing
                   the work or the reference tier is carrying it.
`hazard_table`     per hazard: how many decoys it rejected and whether it
                   rejected anything real. A hazard whose rule never fires
                   is not solved, it is untested.
"""

from __future__ import annotations

from collections import Counter, defaultdict

from ..catalogue import CATALOGUE, HAZARDS, census_families
from .classifier import ChannelCall

CALLS_HEADER = ("accession", "gene_symbol", "expected_family", "called_family",
                "superfamily", "status", "confidence", "margin",
                "filter_string", "decisive_tier", "hazards", "conflicts",
                "evidence")


def _decisive(call: ChannelCall) -> str:
    return next((e.tier for e in call.evidence if e.decisive), "")


def calls_table(calls: list[ChannelCall],
                expected: dict[str, str] | None = None) -> list[tuple]:
    expected = expected or {}
    rows = []
    for c in calls:
        rows.append((
            c.query, c.gene_symbol, expected.get(c.query, ""), c.family,
            c.superfamily, c.status, c.confidence, f"{c.margin:.3f}",
            c.filter_string, _decisive(c), ";".join(c.hazards),
            ";".join(c.conflicts),
            " | ".join(f"{e.tier}: {e.detail}" for e in c.evidence),
        ))
    return rows


def confusion_table(calls: list[ChannelCall],
                    expected: dict[str, str]) -> list[tuple]:
    """`(expected, called, n)` — every cell that is not empty."""
    cells: Counter = Counter()
    for c in calls:
        exp = expected.get(c.query, "")
        if not exp:
            continue
        cells[(exp, c.family or "unassigned")] += 1
    return [(e, g, n) for (e, g), n in sorted(cells.items())]


def recall_table(calls: list[ChannelCall],
                 expected: dict[str, str]) -> list[tuple]:
    """Per expected family: n, correct, recall, and the tier that called it."""
    per: dict[str, list[ChannelCall]] = defaultdict(list)
    for c in calls:
        exp = expected.get(c.query, "")
        if exp:
            per[exp].append(c)
    rows = []
    for fam in sorted(per):
        got = per[fam]
        correct = [c for c in got if c.family == fam]
        tiers = Counter(_decisive(c) for c in correct)
        conf = Counter(c.confidence for c in correct)
        rows.append((
            fam,
            CATALOGUE[fam].name if fam in CATALOGUE else "",
            len(got), len(correct), f"{len(correct) / len(got):.3f}",
            ";".join(f"{t or 'none'}={n}" for t, n in tiers.most_common()),
            ";".join(f"{k}={v}" for k, v in conf.most_common()),
            ";".join(sorted({c.family or "unassigned" for c in got
                             if c.family != fam})),
        ))
    return rows


def hazard_table(calls: list[ChannelCall],
                 expected: dict[str, str]) -> list[tuple]:
    """Per hazard: did its rule fire, on what, and did it cost anything?"""
    rows = []
    for h in HAZARDS:
        involved = [c for c in calls if h.hid in c.hazards]
        correct = sum(1 for c in involved
                      if expected.get(c.query, "") == c.family)
        wrong = [c for c in involved
                 if expected.get(c.query, "") not in ("", c.family)]
        rows.append((
            h.hid, h.severity, ";".join(h.families), len(involved), correct,
            len(wrong), ";".join(c.query for c in wrong[:6]),
            h.test_owner, h.discriminator.replace("\n", " ")[:160],
        ))
    return rows


def coverage_table(calls: list[ChannelCall]) -> list[tuple]:
    """Which census families the panel actually exercised — the honest denominator."""
    seen = Counter(c.family for c in calls if c.family)
    rows = []
    for fam in census_families():
        rows.append((fam.key, fam.superfamily, len(fam.human_genes),
                     seen.get(fam.key, 0)))
    return rows


def summary_counts(calls: list[ChannelCall],
                   expected: dict[str, str] | None = None) -> dict[str, int]:
    expected = expected or {}
    conf = Counter(c.confidence for c in calls)
    out = {
        "n": len(calls),
        "called": sum(1 for c in calls if c.family),
        "superfamily_only": sum(1 for c in calls if not c.family and c.superfamily),
        "unassigned": sum(1 for c in calls if not c.family and not c.superfamily),
        "with_conflict": sum(1 for c in calls if c.conflicts),
        "census_members": sum(1 for c in calls if c.census_member()),
    }
    out.update({f"conf_{k}": v for k, v in conf.items()})
    if expected:
        judged = [c for c in calls if expected.get(c.query)]
        out["judged"] = len(judged)
        out["correct"] = sum(1 for c in judged if c.family == expected[c.query])
    return out


def write_tsv(path, header: tuple, rows: list[tuple]):
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path
