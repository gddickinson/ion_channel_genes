"""S3b — render `results/panel_sweep/report.md` purely from the committed
tables (D13).

    python3 scripts/s3b_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s3_assign import MIN_COVERAGE, MIN_SCORE, REL_MARGIN  # noqa: E402
from scripts.s3_hmm_lib import read_tsv                            # noqa: E402
from scripts.s3_report import pct, table                           # noqa: E402
from scripts.s3b_jackhmmer import EVALUE                           # noqa: E402
from scripts.s3b_kill import (                                     # noqa: E402
    MAX_GROWTH, MAX_ITER, MAX_OTHER_RISE, MIN_PREV_FOR_K2,
)
from scripts.s3b_lib import OUT_DIR                                # noqa: E402

CENSUS = ("channel", "channel_contested")
TABLES = ("sweep_runs", "jackhmmer_seeds", "jackhmmer_runs", "jackhmmer_rounds",
          "jackhmmer_completeness", "jackhmmer_candidates", "instrument_check",
          "census_status", "domain_search_missed", "missed_by_group",
          "missed_records", "family_by_species", "human_recall")


def load():
    T = {n: read_tsv(OUT_DIR / f"{n}.tsv") for n in TABLES}
    J = {n: json.loads((OUT_DIR / f"{n}.json").read_text())
         for n in ("summary", "sweep_db")}
    return T, J


def section_sweep(T, J) -> list[str]:
    runs = T["sweep_runs"]
    secs = sum(float(r["seconds"]) for r in runs)
    inst = T["instrument_check"]
    same = sum(int(r["records"]) for r in inst if r["family"] == "same")
    tot = sum(int(r["records"]) for r in inst)
    diff = [r for r in inst if r["family"] == "differs"]
    return [
        "## 1. The sweep",
        "",
        f"S3a's **{len(runs)} profiles**, unchanged (SHA-256 per profile in "
        f"`sweep_runs.tsv`), searched against S4's declared denominator: "
        f"**{J['sweep_db']['sequences']:,} canonical entries from "
        f"{J['sweep_db']['species']} reference proteomes** (release 2026_03, "
        f"D35; DB SHA-256 `{J['sweep_db']['db_sha256'][:16]}…`). `-Z` fixed "
        f"to the DB size; {secs / 3600:.1f} profile-hours. Assignment is "
        f"S3a's D32 rule with its constants imported, not restated: "
        f"≥ {MIN_SCORE:.0f} bits, ≥ {MIN_COVERAGE:.0%} of the profile's match "
        f"states, ≥ {REL_MARGIN:.0%} relative margin over the best other "
        "profile.",
        "",
        f"**Instrument check.** {tot:,} panel entries are also census v2 "
        f"records, so S3a already scored the same sequence with the same "
        f"profiles; only `-Z` differs, which moves domain i-E-values and so "
        f"coverage. **Same verdict and family on {same:,} / {tot:,} "
        f"({pct(same, tot)}).** The differences:",
        "",
        table(diff, ["s3a_p_call", "s3b_p_call", "records"]) if diff else "_none_",
        "",
    ]


def section_census(T, J) -> list[str]:
    S = J["summary"]
    st = T["census_status"]
    by_basis = defaultdict(int)
    for r in st:
        by_basis[r["v3_basis"].split(":")[0] if r["v3_basis"].startswith("v3a")
                 else r["v3_basis"]] += int(r["records"])
    fams = T["domain_search_missed"]
    census_f = [r for r in fams if r["status"] in CENSUS]
    aux = [r for r in fams if r["status"] not in CENSUS]
    aux_miss = sum(int(r["missed_by_domain_search"]) for r in aux)
    called = sum(int(r["called"]) for r in census_f)
    missed = sum(int(r["missed_by_domain_search"]) for r in census_f)
    m_high = sum(int(r["missed_high"]) for r in census_f)
    top = sorted(census_f, key=lambda r: -int(r["missed_by_domain_search"]))
    grp = defaultdict(Counter)
    for r in T["missed_by_group"]:
        if r["family"] not in {f["family"] for f in census_f}:
            continue
        grp[r["group"]]["v2"] += int(r["in_census_v2"])
        grp[r["group"]]["miss"] += int(r["missed_by_domain_search"])
        grp[r["group"]]["high"] += int(r["missed_high"])
    grows = [{"group": g, "in_census_v2": c["v2"],
              "missed_by_domain_search": c["miss"], "missed_high": c["high"],
              "missed": pct(c["miss"], c["v2"] + c["miss"])}
             for g, c in sorted(grp.items(), key=lambda x: -x[1]["miss"])]
    return [
        "## 2. Census v3 on the panel — what domain search missed",
        "",
        "The merge rule (fixed before the tables were read): an entry in census "
        "v2 keeps its v3a call; an entry outside it takes the S3b profile call "
        "(`panel_profile`); an entry only a jackhmmer run reaches is a "
        "`candidate`, never a family call (D14, D33).",
        "",
        f"**{S['with_evidence']:,} of {S['panel_entries']:,} panel entries** "
        f"carry evidence; {S['in_census_v2']:,} are census v2 records.",
        "",
        table(st, ["v3_basis", "v3_status", "records"]),
        "",
        f"**Channel-family calls: {called:,}, of which {missed:,} "
        f"({pct(missed, called)}) are outside census v2** — members of a "
        "census family that carry none of the 67 enumerated pore signatures, "
        "so no domain search over those signatures could have found them. "
        f"**Only {m_high} of them are high-confidence profile calls** (margin "
        "≥ 30 % and ≥ half the profile). The other "
        f"{missed - m_high} are medium, and the largest medium blocks are "
        "**not channels**: the TRPN profile is mostly its ankyrin-repeat "
        "array and the LRRC8 profile half leucine-rich repeat, so an "
        "ankyrin- or LRR-repeat protein (ANKRD52, ankyrin, IκB, titin…) can "
        "cover 30 % of the profile with the repeat module alone and pass "
        "D32's coverage gate. D32 assumed a shared module is a small part of "
        "the profile; for repeat-dominated profiles it is not (emergent). "
        "Read the high column as the finding and the rest as an upper bound. "
        "Per family, largest first:",
        "",
        table(top, ["family", "superfamily", "called", "in_census_v2",
                    "missed_by_domain_search", "missed_high", "missed_frac"], 25),
        "",
        f"Control families are reported apart: {aux_miss:,} of their "
        f"{sum(int(r['called']) for r in aux):,} panel calls are outside census "
        "v2 — expected, since auxiliary subunits and non-channel homologues "
        "carry no pore signature by definition (D23).",
        "",
        "By panel group (census families only):",
        "",
        table(grows, ["group", "in_census_v2", "missed_by_domain_search",
                      "missed_high", "missed"]),
        "",
        "Every missed record is listed in `missed_records.tsv`; the family × "
        "species counts S10 starts from are `family_by_species.tsv`.",
        "",
    ]


def section_jackhmmer(T) -> list[str]:
    runs = T["jackhmmer_runs"]
    v = Counter(r["verdict"] for r in runs)
    rules = Counter(r["rule"] for r in runs if r["rule"])
    conv = sum(r["converged"] == "True" for r in runs)
    killed = [r for r in runs if r["verdict"] == "killed"]
    hrs = sum(float(r["seconds"] or 0) for r in runs) / 3600
    return [
        "## 3. jackhmmer to convergence (D10)",
        "",
        f"One run per census family from a **derived seed**: the panel entry "
        f"with the highest-scoring high-confidence profile call to that family "
        f"(`jackhmmer_seeds.tsv`, written before any run). `-N {MAX_ITER} -E "
        f"{EVALUE} --incE {EVALUE}` over the panel DB; {hrs:.1f} run-hours.",
        "",
        f"The kill criterion (`s3b_kill.py`, ported from the IP3R project with "
        f"its constants): **K1** the share of included targets the profiles "
        f"call to *another* family rises > {MAX_OTHER_RISE:.0%} points above "
        f"round 1; **K2** the included set grows > {MAX_GROWTH:.0f}× (from ≥ "
        f"{MIN_PREV_FOR_K2}); **K3** {MAX_ITER} rounds without converging. "
        "Rounds from the first firing on are excluded.",
        "",
        f"**{len(runs)} runs: {v['clean']} clean, {v['killed']} killed "
        f"({', '.join(f'{k} {n}' for k, n in sorted(rules.items())) or 'none'})"
        f", {v['no_run']} not run; {conv} converged.**",
        "",
        table(killed, ["family", "rule", "killed_at", "accepted_rounds",
                       "n_accepted", "reason"]) if killed else "",
        "",
    ]


def section_completeness(T) -> list[str]:
    rows = [r for r in T["jackhmmer_completeness"] if r["verdict"] != "no_run"]
    prof = sum(int(r["profile_calls"]) for r in rows)
    rec = sum(int(r["recovered"] or 0) for r in rows)
    low = sorted([r for r in rows if r["recall"] and float(r["recall"]) < 0.9],
                 key=lambda r: float(r["recall"]))
    cand, cand_f = Counter(), Counter()
    for r in T["jackhmmer_candidates"]:
        cand[r["profile_verdict"]] += int(r["records"])
        cand_f[r["family"]] += int(r["records"])
    clean = [r for r in rows if r["verdict"] == "clean"]
    cp = sum(int(r["profile_calls"]) for r in clean)
    cr = sum(int(r["recovered"] or 0) for r in clean)
    cu = sum(int(r["jh_accepted"]) - int(r["recovered"] or 0)
             - int(r["jh_other_family"]) for r in clean)
    cols = ["family", "verdict", "profile_calls", "recovered", "recall",
            "jh_accepted", "jh_other_family", "jh_own_superfamily_only",
            "jh_module", "jh_low_score", "jh_no_hit", "top_other_families"]
    return [
        "## 4. The completeness argument",
        "",
        "Two directions. **Does a single-sequence iterated search reach what "
        "the profile reaches?** Recall of each family's panel profile calls by "
        f"its own accepted jackhmmer set: **{rec:,} / {prof:,} "
        f"({pct(rec, prof)})**. **Does it reach what the profile cannot call?** "
        "Accepted targets the profile library calls to no family are the upper "
        "bound on what the sweep misses; those outside census v2 are the "
        "`candidate` rows, by profile verdict: "
        + (", ".join(f"{k} {n:,}" for k, n in cand.most_common()) or "none") + ".",
        "",
        f"**The claim rests on the {len(clean)} clean runs**, which recover "
        f"**{cr:,} / {cp:,} ({pct(cr, cp)})** of their families' profile calls "
        f"and include only **{cu:,}** targets the profile library calls to no "
        "family (at most, what the sweep misses in those families). Killed "
        "runs are reported with their pre-kill rounds; a K3 run contributes "
        "no candidates (D36).",
        "",
        "Candidate rows by the run that found them (a target in several runs "
        "counts once per run): "
        + ", ".join(f"{k} {n:,}" for k, n in cand_f.most_common(6)) + ". "
        "The TRPA and TRPN runs dominate: their pre-K1 rounds already hold "
        "thousands of ankyrin-repeat proteins the profiles call `module` or not "
        "at all. **K1 cannot see this drift** — uncalled targets are exempt by "
        "design, as in the parent project, whose S19 measured the same blind "
        "spot. The candidates are counted, never called.",
        "",
        f"Families whose jackhmmer run recovers < 90 % of their profile calls "
        f"({len(low)}):",
        "",
        table(low, cols) if low else "_none_",
        "",
        "All families: `jackhmmer_completeness.tsv`.",
        "",
    ]


def section_human(T) -> list[str]:
    hr = T["human_recall"]
    v = Counter(r["verdict"] for r in hr)
    bad = [r for r in hr if r["verdict"] != "right_family"]
    return [
        "## 5. Human positive control",
        "",
        f"The {len(hr)} human census genes against the human reference proteome "
        f"(gene symbol used to *score* the call, never to make it — H15): "
        + ", ".join(f"**{k} {n}**" for k, n in v.most_common()) + ".",
        "",
        table(bad, ["gene", "family", "target", "v3_family", "v3_status",
                    "v3_basis", "p_call", "p_family", "verdict"]) if bad else "",
        "",
    ]


def main() -> int:
    T, J = load()
    out = ["# S3b — profile sweep of the panel proteomes + jackhmmer → census v3",
           "",
           "_Rendered by `scripts/s3b_report.py` from the tables in this "
           "directory (D13). Do not edit by hand._", ""]
    for part in (section_sweep(T, J), section_census(T, J),
                 section_jackhmmer(T), section_completeness(T), section_human(T)):
        out += part
    (OUT_DIR / "report.md").write_text("\n".join(out) + "\n")
    print(f"wrote {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
