"""S3b — D10's coded kill criterion for the per-family jackhmmer runs.

Ported from `../ip3r_genes/scripts/s3_kill.py` (method and constants, not
results), where it guarded two sister families. Here a run is seeded in one
of 68 census families and can drift into any of the other 90 the profile
library knows, so K1 is widened from "the sister family" to "any other
family" — measured with the S3b profile calls, which exist for every panel
entry before any jackhmmer run starts, so the rule is evaluated on evidence
fixed in advance.

  K1  Other-family *drift*. The share of the round's included targets that
      the profile library calls to a family other than the seed's rises more
      than MAX_OTHER_RISE above its round-1 value. Round 1 is the baseline
      and only the rise counts: the parent project measured that a single
      ITPR1 sequence already returns 32 % RyRs at round 1 — real homology,
      not drift — and a flat ceiling fired on every run.
  K2  Explosive growth: the included set grows past MAX_GROWTH × the
      previous round's, once past a floor of MIN_PREV_FOR_K2.
  K3  Ceiling: MAX_ITER rounds *without* converging. A run converging on
      its last allowed round has converged.

The first rule to fire ends the run's usable rounds; later rounds are
reported and excluded. Targets no profile calls to any family (module-only,
no hit, ambiguous) are reported per round and never a kill rule — that
bucket also holds real members too fragmentary for D32's coverage gate,
which is exactly what an iterative search exists to find.
"""

from __future__ import annotations

MAX_OTHER_RISE = 0.10
MIN_ROUND_FOR_K1 = 2
MAX_GROWTH = 10.0
MIN_PREV_FOR_K2 = 20
MAX_ITER = 10


def evaluate(rounds: list[dict], family_of: dict[str, str],
             own_family: str, converged: bool = False) -> dict:
    """Apply K1–K3 to one run's parsed rounds.

    `rounds` are `s3b_lib.parse_jackhmmer_log()` dicts, each with an
    `included` list of target names; `family_of` maps a target to its S3b
    profile family call ("" when the profile makes none).
    """
    per_round: list[dict] = []
    killed_at = None
    rule = reason = ""
    prev = None
    baseline: float | None = None
    for rd in rounds:
        inc = rd.get("included") or []
        n = len(inc)
        calls = [family_of.get(t, "") for t in inc]
        n_own = sum(1 for c in calls if c == own_family)
        n_other = sum(1 for c in calls if c and c != own_family)
        frac = round(n_other / n, 4) if n else 0.0
        if baseline is None and n:
            baseline = frac
        rise = round(frac - baseline, 4) if baseline is not None else 0.0
        growth = round(n / prev, 2) if prev else None
        per_round.append({
            "round": rd["round"], "new_targets": rd["new_targets"],
            "n_included": n, "n_own": n_own, "n_other_family": n_other,
            "n_uncalled": n - n_own - n_other, "other_frac": frac,
            "other_rise": rise, "growth": growth if growth is not None else ""})
        if killed_at is None:
            if rd["round"] >= MIN_ROUND_FOR_K1 and rise > MAX_OTHER_RISE:
                killed_at, rule = rd["round"], "K1"
                reason = (f"round {rd['round']}: other-family share {frac:.1%} "
                          f"({n_other}/{n}), {rise:+.1%} on the round-1 "
                          f"baseline {baseline:.1%} (limit "
                          f"{MAX_OTHER_RISE:.0%} points)")
            elif prev and prev >= MIN_PREV_FOR_K2 and growth \
                    and growth > MAX_GROWTH:
                killed_at, rule = rd["round"], "K2"
                reason = (f"round {rd['round']} grew the included set to {n}, "
                          f"{growth}× the previous {prev} (limit {MAX_GROWTH}×)")
        prev = n or prev
    if killed_at is None and not converged and len(rounds) >= MAX_ITER:
        killed_at, rule = MAX_ITER, "K3"
        reason = (f"{MAX_ITER}-round ceiling without converging; no "
                  "completeness claim may rest on this run")
    return {"verdict": "killed" if killed_at else "clean",
            "accepted_rounds": (killed_at - 1) if killed_at else len(rounds),
            "killed_at": killed_at or "", "rule": rule,
            "reason": reason or "no kill rule fired", "per_round": per_round}


def accepted_targets(rounds: list[dict], accepted_rounds: int) -> set[str]:
    """Targets included in the last accepted round (empty if none)."""
    for rd in reversed(rounds):
        if rd["round"] <= accepted_rounds:
            return set(rd.get("included") or [])
    return set()
