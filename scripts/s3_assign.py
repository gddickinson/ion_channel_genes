"""S3a — best-profile assignment across every family profile (D7, D14).

Every target is scored against all family profiles. The best-scoring
profile is the call only when it passes three gates, fixed before any
census result was seen:

1. **Score** — the winner's full-sequence bit score ≥ `MIN_SCORE` (30 bits,
   the parent project's floor).
2. **Coverage** — the winner's significant domains span ≥ `MIN_COVERAGE`
   (0.30) of the profile's match states. This is D30's number applied to
   profiles: a protein sharing one module with a family (the cNMP domain of
   HCN/CNG/EAG, the T1 domain of Kv/KCTD, SPRY in RyR) matches a small
   fraction of that family's profile and is *module* evidence, not a call.
   Relative, not absolute, because the profiles run from ~60 to ~5,000
   match states.
3. **Margin** — the winner beats the runner-up (the best *other* profile,
   whatever its coverage) by ≥ `REL_MARGIN` (0.10, D7) of its own score.
   Relative because bit scores scale with alignable length.

Inside the margin the target is not a family call. If every profile within
the margin belongs to one superfamily the call is `superfamily_only` at
that superfamily — the same honest partial answer S2 gives — otherwise
`ambiguous`. Nothing here reads a gene symbol, a name or a length.
"""

from __future__ import annotations

import heapq
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import hits_by_target, iter_domtblout  # noqa: E402
from src.catalogue import registry  # noqa: E402

MIN_SCORE = 30.0
MIN_COVERAGE = 0.30
REL_MARGIN = 0.10
HIGH_MARGIN = 0.30          # confidence: high = this margin and ≥ half the profile
KEEP = 6                    # hits kept per target (winner, runner-up, band)

ASSIGN_FIELDS = ["target", "p_call", "p_family", "p_superfamily",
                 "p_confidence", "win_score", "win_coverage",
                 "win_target_coverage", "runner", "runner_score",
                 "rel_margin", "band", "p_reason"]


def fam_superfamily() -> dict[str, str]:
    return {f.key: f.superfamily for f in registry.families()}


def collect_hits(domtbl_paths, keep: int = KEEP,
                 targets: set[str] | None = None) -> dict[str, list[tuple]]:
    """target → top-`keep` hits as (score, profile, hmm_cov, tcov, evalue).

    Streams one profile's domtblout at a time so memory holds only the
    per-target top lists, not every (target, profile) row of the sweep.
    """
    top: dict[str, list[tuple]] = {}
    for p in domtbl_paths:
        rows = iter_domtblout(p)
        if targets is not None:
            rows = (r for r in rows if r["target_name"] in targets)
        for (t, prof), h in hits_by_target(rows).items():
            item = (h["score"], prof, h["hmm_coverage"],
                    h["target_coverage"], h["evalue"])
            lst = top.setdefault(t, [])
            if len(lst) < keep:
                heapq.heappush(lst, item)
            elif item > lst[0]:
                heapq.heapreplace(lst, item)
    return {t: sorted(v, reverse=True) for t, v in top.items()}


def assign_one(target: str, hits: list[tuple], sf_of: dict[str, str],
               exclude: set[str] = frozenset()) -> dict:
    """One target's call from its sorted hit list (best first)."""
    hits = [h for h in hits if h[1] not in exclude]
    row = {k: "" for k in ASSIGN_FIELDS}
    row["target"] = target
    if not hits:
        row.update(p_call="no_hit", p_confidence="none",
                   p_reason="no profile reported a hit")
        return row
    w = hits[0]
    r = hits[1] if len(hits) > 1 else None
    win, lose = w[0], (r[0] if r else 0.0)
    rel = (win - lose) / win if win > 0 else 0.0
    band = [h[1] for h in hits if win > 0 and h[0] >= win * (1 - REL_MARGIN)]
    row.update(win_score=round(win, 1), win_coverage=w[2],
               win_target_coverage=w[3],
               runner=r[1] if r else "", runner_score=round(lose, 1) if r else "",
               rel_margin=round(rel, 4), band=",".join(band))
    if win < MIN_SCORE:
        row.update(p_call="low_score", p_confidence="none",
                   p_reason=f"best {w[1]} {win:.1f} bits < {MIN_SCORE:.0f}")
    elif w[2] < MIN_COVERAGE:
        row.update(p_call="module", p_confidence="none",
                   p_reason=f"best {w[1]} spans {w[2]:.0%} of the profile "
                            f"(< {MIN_COVERAGE:.0%}) — shared module")
    elif rel < REL_MARGIN:
        sfs = {sf_of.get(p, "") for p in band}
        if len(sfs) == 1 and "" not in sfs:
            sf = sfs.pop()
            row.update(p_call="superfamily_only", p_superfamily=sf,
                       p_confidence="low",
                       p_reason=f"{', '.join(band)} within {rel:.1%} "
                                f"(< {REL_MARGIN:.0%}, D7); one superfamily")
        else:
            row.update(p_call="ambiguous", p_confidence="none",
                       p_reason=f"{', '.join(band)} within {rel:.1%} "
                                "across superfamilies")
    else:
        high = rel >= HIGH_MARGIN and w[2] >= 0.5
        row.update(p_call="family", p_family=w[1],
                   p_superfamily=sf_of.get(w[1], ""),
                   p_confidence="high" if high else "medium",
                   p_reason=f"{w[1]} wins by {win - lose:.1f} bits "
                            f"({rel:.1%} of {win:.1f}) over "
                            f"{w[2]:.0%} of the profile")
    return row


def assign_all(top: dict[str, list[tuple]]) -> dict[str, dict]:
    sf_of = fam_superfamily()
    return {t: assign_one(t, h, sf_of) for t, h in top.items()}
