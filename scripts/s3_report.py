"""S3a — render `results/census_v3/report.md` purely from the committed tables (D13).

    python3 scripts/s3_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s3_assign import MIN_COVERAGE, MIN_SCORE, REL_MARGIN  # noqa: E402
from scripts.s3_hmm_lib import OUT_DIR, read_tsv                   # noqa: E402


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f} %" if b else "—"


def table(rows: list[dict], cols: list[str], limit: int | None = None) -> str:
    rows = rows[:limit] if limit else rows
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def load():
    T = {n: read_tsv(OUT_DIR / f"{n}.tsv") for n in (
        "profile_build", "seed_manifest", "benchmark_calls", "benchmark_recall",
        "v3_basis", "v3_status", "v3_families", "calibration",
        "calibration_by_domain", "conflicts", "superfamily_only_resolved",
        "unassigned_fate", "human_recall_v3", "v2_status_x_profile",
        "sweep_runs", "manifest", "call_support")}
    J = {n: json.loads((OUT_DIR / f"{n}.json").read_text())
         for n in ("summary", "benchmark_summary", "sweep_db")}
    return T, J


def section_profiles(T) -> list[str]:
    pb = T["profile_build"]
    thin = [r for r in pb if int(r["n_seeds"]) <= 2]
    rules = Counter(r["rule"] for r in T["seed_manifest"])
    return [
        "## 1. The profile library",
        "",
        f"**{len(pb)} profiles**, one per catalogue family (census and control "
        f"families alike, so a decoy has a profile to win), from "
        f"**{len(T['seed_manifest'])} seeds**: R1 curated human genes "
        f"{rules['R1']}, R2 catalogue exemplars {rules['R2']}, R3 reviewed S2 "
        f"family calls (one per species, round-robin over groups) {rules['R3']}. "
        f"Seeds are disjoint across families by construction (the build "
        f"aborts otherwise). MAFFT L-INS-i single-threaded, hmmbuild; "
        f"SHA-256 of every seed set, alignment and profile in "
        f"`profile_build.tsv`. Match states total "
        f"{sum(int(r['match_states']) for r in pb):,}.",
        "",
        f"**{len(thin)} profiles rest on one or two sequences** — the families "
        "with no human genes and no S2 family calls. They are reported as thin, "
        "not padded:",
        "",
        table(sorted(thin, key=lambda r: r["family"]),
              ["family", "superfamily", "n_seeds", "match_states"]),
        "",
    ]


def section_benchmark(T, J) -> list[str]:
    B = J["benchmark_summary"]
    A, O = B["A_s1_panel"], B["B_orthologues"]
    calls = T["benchmark_calls"]
    miss = [r for r in calls if r["outcome"] not in ("correct", "decoy_ok")]
    dec = [r for r in calls if r["set"] == "A_s1_panel" and r["outcome"].startswith("decoy")]
    own = sum(r["p_family"] == r["expected_family"] for r in dec)
    nohit = sum(r["p_call"] == "no_hit" for r in dec)
    return [
        "## 2. The instrument, measured before use",
        "",
        f"Gates fixed before any census result: best profile ≥ {MIN_SCORE:.0f} "
        f"bits, spanning ≥ {MIN_COVERAGE:.0%} of its match states (D30's "
        f"number, applied to profiles), and beating the runner-up by ≥ "
        f"{REL_MARGIN:.0%} of its own score (D7).",
        "",
        f"**A. S1 control panel, leave-one-out** ({A['loo']} rebuilt profiles; "
        f"{A['sole_seed']} queries are their family's only seed and have no "
        f"LOO profile): **{A['correct']}/{A['positives']} positives correct**, "
        f"{A['decoy_ok']}/{A['decoy_ok'] + A['decoy_called_channel']} decoys "
        f"not called to a channel family. For comparison, S1's three-tier "
        f"classifier scored 50/72 and S2's reference-free calls 45 on the same "
        f"proteins. (S1 counted 72 positives and 25 decoys; here the CLC "
        f"transporter, catalogue status `transporter`, is scored as a decoy.) "
        f"The stronger decoy test — called to its **own** control family, a "
        f"positive call rather than a mere non-channel one — holds for "
        f"**{own}/{len(dec)}**; the other {len(dec) - own} "
        f"({nohit} with no profile hit at all) come from auxiliary-subunit "
        f"families that pool *unrelated* proteins under one catalogue key — "
        f"`assoc_k_beta` holds an aldo-keto reductase (Kvβ), the BK β "
        f"subunits, the LRRC γ subunits and DPP6/10 — so once the query "
        f"is left out, no remaining seed is its homologue. One profile cannot "
        f"represent a family that is not a family; the census is unaffected "
        f"(none of these proteins carries a pore signature).",
        "",
        f"**B. Held-out orthologues** — reviewed, non-human, non-seed census "
        f"records whose gene symbol is a catalogue human gene (the symbol "
        f"scores the call, it never makes it — H15): **{O['correct']}/"
        f"{O['positives']} correct**, {O['decoy_ok']} control-family "
        f"orthologues, {O['decoy_called_channel']} called to a channel. "
        f"**This set is almost entirely vertebrate** (see "
        f"`benchmark_calls.tsv`), so it tests the profiles on the lineage "
        f"they were seeded from; section 3 is the non-vertebrate test.",
        "",
        "Every benchmark miss:",
        "",
        table(miss, ["set", "accession", "expected_family", "loo", "outcome",
                     "p_call", "p_family", "win_score", "win_coverage",
                     "runner", "rel_margin"]),
        "",
    ]


def section_census(T, J) -> list[str]:
    S = J["summary"]
    n = S["records"]
    basis = {r["basis"]: int(r["records"]) for r in T["v3_basis"]}
    fam_calls = basis.get("both", 0) + basis.get("s2_only", 0) + basis.get("profile_only", 0)
    cal = T["calibration"]
    ag = sum(int(r["agree"]) for r in cal)
    dis = sum(int(r["disagree"]) for r in cal)
    tot = sum(int(r["v2_calls"]) for r in cal)
    sfo = defaultdict(Counter)
    for r in T["superfamily_only_resolved"]:
        sfo[r["v2_superfamily"]][r["v3_outcome"]] += int(r["records"])
    sfo_rows = []
    for sf, c in sorted(sfo.items(), key=lambda x: -sum(x[1].values())):
        tot_sf = sum(c.values())
        resolved = sum(v for k, v in c.items() if k not in ("superfamily_only", "unassigned", "conflict"))
        top = ", ".join(f"{k} {v:,}" for k, v in c.most_common(4))
        sfo_rows.append({"superfamily": sf, "v2_superfamily_only": f"{tot_sf:,}",
                         "resolved_to_family": f"{resolved:,} ({pct(resolved, tot_sf)})",
                         "outcomes": top})
    un = T["unassigned_fate"]
    un_tot = sum(int(r["records"]) for r in un)
    un_by = Counter()
    for r in un:
        un_by[r["p_call"]] += int(r["records"])
    return [
        "## 3. Calibration against S2 (seeds excluded)",
        "",
        f"On the {tot:,} records S2 called to a family (seed accessions "
        f"removed), the profile makes the **same call on {ag:,}** and a "
        f"**different family call on {dis:,}** — agreement "
        f"{pct(ag, ag + dis)} where both instruments speak. The two read "
        f"different evidence (Pfam architecture and filter motif vs a "
        f"full-length profile margin), so agreement is a measurement, not a "
        f"tautology.",
        "",
        table(T["calibration_by_domain"], ["domain", "v2_calls", "agree",
              "disagree", "profile_superfamily_only", "profile_module",
              "profile_low_score", "profile_no_hit", "agreement"]),
        "",
        "Per family (largest 30; full table `calibration.tsv`):",
        "",
        table(cal, ["family", "v2_calls", "agree", "disagree",
                    "profile_superfamily_only", "profile_ambiguous",
                    "profile_module", "agreement"], 30),
        "",
        "## 4. Census v3a",
        "",
        f"**{n:,} records**; **{fam_calls:,} ({pct(fam_calls, n)}) now carry "
        f"a family call** — both instruments {basis.get('both', 0):,}, S2 "
        f"only {basis.get('s2_only', 0):,}, profile only "
        f"{basis.get('profile_only', 0):,} — against 346,627 (27.8 %) in S2. "
        f"Superfamily only {basis.get('superfamily_only', 0):,}; unassigned "
        f"{basis.get('unassigned', 0):,}; **conflict "
        f"{basis.get('conflict', 0):,}**, kept and counted.",
        "",
        "**What the profile-only calls have not been tested on.** The "
        "superfamily splits below (Cys-loop, iGluR, DEG/ENaC, P2X, CLC, "
        "TMEM16) are families S2 never called, so section 3 cannot check "
        "them, and section 2's orthologue test is vertebrate. Where those "
        "calls fall in invertebrates, plants and protists they rest on "
        "profiles seeded mostly from human genes and have no independent "
        "check yet.",
        "",
        table(T["v3_status"], ["status", "records"]),
        "",
        "### What happened to S2's superfamily-only calls",
        "",
        table(sfo_rows, ["superfamily", "v2_superfamily_only",
                         "resolved_to_family", "outcomes"]),
        "",
        "### What happened to S2's unassigned records",
        "",
        f"{un_tot:,} records S2 left unassigned: "
        + "; ".join(f"`{k}` {v:,} ({pct(v, un_tot)})" for k, v in un_by.most_common())
        + ". A `module` verdict is a profile match over less than "
        f"{MIN_COVERAGE:.0%} of the family profile and is deliberately not a "
        "call; `no_hit` means no profile reported the record at all. S2 "
        "traced most of its unassigned records to three co-domain signatures "
        "(cNMP, SBP_bac_3, PAS); which of these two verdicts those records "
        "received is not broken down here.",
        "",
        table(un, ["p_call", "v3_family", "records"], 15),
        "",
        "### Calls that rest on S2 alone",
        "",
        "`s2_only` is a family call S2 made that the profile did not "
        "confirm — usually because the profile abstained (`module`: the "
        "record matches under 30 % of the family profile). These calls are "
        "only as good as S2's rule, and section 6 shows one rule that is "
        "not good enough. Largest groups (`call_support.tsv`):",
        "",
        table(sorted((r for r in T["call_support"] if r["v3_basis"] == "s2_only"),
                     key=lambda r: -int(r["records"])),
              ["v3_family", "s2_tier", "p_call", "records"], 12),
        "",
        "### Conflicts",
        "",
        "Where S2 and the profile name different families. S2's deciding "
        "tier is shown because a filter-motif call and a profile call "
        "disagreeing is a different finding from an architecture call and a "
        "profile call disagreeing.",
        "",
        table(T["conflicts"], ["s2_call", "profile_family", "s2_tier", "records"], 25),
        "",
    ]


def section_human(T, J) -> list[str]:
    hr = T["human_recall_v3"]
    S = J["summary"]
    wrong = [r for r in hr if r["correct_v3"] == "no"]
    return [
        "## 5. Human census genes",
        "",
        f"**{S['human_correct_v3']}/{S['human_genes']} called to the right "
        f"family in v3a**, against {S['human_correct_s2']}/{S['human_genes']} "
        f"in S2. Most human genes are R1 seeds of their own profile, so this "
        f"row measures the merge, not the profiles' generalisation — "
        f"section 2 is the generalisation test. Every human gene not "
        f"called right:",
        "",
        table(wrong, ["gene", "expected_family", "s2_family", "p_call",
                      "p_family", "v3_family", "v3_basis"]),
        "",
    ]


def section_external() -> list[str]:
    """The parent projects' censuses as an external check (never an input)."""
    p = OUT_DIR / "external_check.tsv"
    if not p.exists():
        return []
    ext = read_tsv(p)
    rev = read_tsv(OUT_DIR / "external_reverse.tsv")
    src = read_tsv(OUT_DIR / "external_sources.tsv")
    dis = read_tsv(OUT_DIR / "external_disagreements.tsv")
    ip = [r for r in ext if r["project"].startswith("ip3r") and r["parent_call"] in ("ITPR", "RYR")]
    comparable = [r for r in ip if r["verdict"] not in ("not_a_uniprot_accession", "not_in_census_v2")]
    agree = sum(int(r["records"]) for r in comparable if r["verdict"] == "agree")
    n_cmp = sum(int(r["records"]) for r in comparable)
    swapped = sum(int(r["records"]) for r in comparable if r["verdict"] == "swapped_itpr_ryr")
    other = sum(int(r["records"]) for r in comparable if r["verdict"] == "called_other_family")
    sw = [r for r in dis if r["verdict"] == "swapped_itpr_ryr"]
    sw_h4 = sum(r["v2_family"] == "itpr" and r["v2_tier"] == "hazard"
                and r["p_call"] in ("module", "no_hit", "low_score") for r in sw)
    sw_arch = sum(r["parent_reason"].startswith("architecture only") for r in sw)
    pz_un = sum(int(r["records"]) for r in ext
                if r["project"].startswith("piezo") and r["verdict"] == "v3_unassigned")
    sw_len = sorted(int(r["length"]) for r in sw)[len(sw) // 2] if sw else 0
    return [
        "## 6. External check — the parent projects' censuses",
        "",
        "The IP3R and PIEZO projects censused three of this catalogue's "
        "families with their own seeds, profiles and search spaces "
        "(vertebrate reference proteomes and genome sweeps). Their results "
        "are **compared, never imported** (CLAUDE.md: port the method, never "
        "a result). Only UniProt accessions can be matched; genome-derived "
        "models and Ensembl proteins are counted as unmatched.",
        "",
        f"**IP3R census v6, call against call:** of {n_cmp:,} parent ITPR/RYR "
        f"calls on accessions census v2 also holds, census v3a makes the "
        f"**same call on {agree:,} ({pct(agree, n_cmp)})**; ITPR and RYR "
        f"swapped on {swapped:,}; another family on {other:,}; the rest "
        "carry no family call here (verdicts beginning `v3_`).",
        "",
        table(ext, ["project", "parent_call", "verdict", "records"]),
        "",
        f"**The {swapped:,} swapped records are all one error, and it is "
        f"S2's.** {sw_h4:,} of them were called ITPR by S2's hazard rule "
        f"`H4-itpr` — `PF08709` *without* `PF02026`/`PF06459` — while the "
        f"profile abstained (median length {sw_len} aa). `PF08709` is a "
        "domain ITPR and RyR *share* (H4: they share every diagnostic "
        "domain), so the rule is an absence test, which this project's "
        "conventions forbid, and it calls N-terminal RyR fragments IP3 "
        "receptors. On evidence both families carry, the right answer is "
        "no family call — not ITPR, and not necessarily RYR either. The "
        "parent project's basis for each RYR call is in `parent_reason`; "
        f"{sw_arch:,} of the {swapped:,} are its architecture-only calls. "
        "The profile does not repeat the error; the merge inherits it "
        "because an abstaining profile leaves S2's call standing.",
        "",
        f"PIEZO census v5 is a membership list that includes fragments and "
        f"short-motif hits; {pz_un:,} of its records are unassigned here. "
        "Whether those are fragments below this census's gates or real "
        "PIEZOs missed has not been checked.",        "",
        "Reverse direction — of this census's own itpr / ryr / piezo calls, "
        "how many the parent census holds at all. A `no` is expected where "
        "this census reaches beyond the parent's search space "
        "(non-vertebrates for IP3R before its S20 sweep, UniProt entries "
        "outside reference proteomes, a newer UniProt release):",
        "",
        table(rev, ["our_family", "v3_basis", "in_parent_census", "records"]),
        "",
        f"{len(dis):,} individual disagreements are listed in "
        "`external_disagreements.tsv` with both sides' evidence. Parent "
        "files and their SHA-256:",
        "",
        table(src, ["file", "rows", "sha256"]),
        "",
    ]


def main() -> int:
    T, J = load()
    db = J["sweep_db"]
    L = ["# S3a — profile library and best-profile assignment → census v3a",
         "",
         "*Rendered from the tables in this directory by `scripts/s3_report.py`. "
         "Census v2 = UniProtKB 2026_03.*",
         "",
         f"Sweep database: census v2 collapsed to **{db['unique_sequences']:,} "
         f"unique sequences** ({db['records']:,} records, "
         f"{db['unique_residues']:,} residues); every profile searched with "
         f"`-Z` fixed to that size. Wall-clock per profile in `sweep_runs.tsv` "
         f"(total {sum(float(r['seconds']) for r in T['sweep_runs']) / 3600:.1f} "
         f"profile-hours).",
         ""]
    L += section_profiles(T) + section_benchmark(T, J) + section_census(T, J) + section_human(T, J) + section_external()
    L += ["## Files", "", table(T["manifest"], ["file", "bytes", "sha256", "records"]), ""]
    (OUT_DIR / "report.md").write_text("\n".join(L))
    print(f"wrote {OUT_DIR / 'report.md'} ({len(L)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
