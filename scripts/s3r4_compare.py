"""s3r4_compare.py — what census v3a r4 changed against r3, record by record.

    python3 scripts/s3r4_compare.py

Reads the archived r3 census v3a (`<data root>/hmmer/s3/census_v3.r3.tsv.gz`)
and the current one. For every r3 record: its v3 call before and after.
The new profiles can move an r3 record only if they score it, so every
change must involve one of the new families — as the new call, the profile
winner, or the runner-up that pulled the winner inside the D7 margin; a
change that does not is a hard failure. For the
26,783 r4 records: their v3 calls, by family.

Also every **reviewed** r4 record with its UniProt gene symbol, seed flag and
profile call — the H17–H19 look-alike test beyond S3a's set B, which skips
seeds and "Precursor"-flagged records (the symbol only *scores* a call made
without it, H15).

→ `results/census_v3/r4_transitions.tsv`, `r4_new_records.tsv`,
`r4_reviewed_calls.tsv`, `r4_summary.json`, and `r4_report.md` rendered from
them (D13).
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import OUT_DIR, s3_dir, write_tsv  # noqa: E402
from scripts.s3r4_sweep import NEW_PROFILES                 # noqa: E402


def rows(p: Path):
    with gzip.open(p, "rt", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def label(r: dict) -> str:
    return (r["v3_family"] or (f"[{r['v3_superfamily']}]" if r.get("v3_superfamily")
                               else r["v3_status"]))


def reviewed_calls() -> None:
    from scripts.s3_hmm_lib import read_tsv, s2_raw
    seeds = {r["accession"] for r in read_tsv(OUT_DIR / "seed_manifest.tsv")}
    meta = {r["accession"]: r for r in rows(s2_raw() / "census_v2.tsv.gz")
            if r["reviewed"] == "reviewed" and r["shard"].startswith("r4_")}
    out = []
    for r in rows(s3_dir() / "census_v3.tsv.gz"):
        m = meta.get(r["accession"])
        if m:
            out.append({"accession": r["accession"], "gene": m["gene"],
                        "organism": m["organism"], "seed": int(r["accession"] in seeds),
                        "p_family": r["p_family"], "p_confidence": r["p_confidence"],
                        "rel_margin": r["rel_margin"], "runner": r["runner"],
                        "v3_family": r["v3_family"], "v3_status": r["v3_status"]})
    write_tsv(OUT_DIR / "r4_reviewed_calls.tsv", list(out[0]), out)


def render(summary: dict) -> None:
    from scripts.s3_hmm_lib import read_tsv
    tr = read_tsv(OUT_DIR / "r4_transitions.tsv")
    nr = read_tsv(OUT_DIR / "r4_new_records.tsv")
    rv = read_tsv(OUT_DIR / "r4_reviewed_calls.tsv")
    b = json.loads((OUT_DIR / "benchmark_summary.json").read_text())["B_orthologues"]
    fams = Counter()
    for r in nr:
        fams[r["v3_call"]] += int(r["records"])
    L = ["# Census v2/v3a revision r4 — the families added after S20 (D43)", "",
         "Rendered by `scripts/s3r4_compare.py` from `r4_transitions.tsv`, "
         "`r4_new_records.tsv`, `r4_reviewed_calls.tsv` and `benchmark_summary.json` (D13).",
         "", f"**{summary['r4_records']:,} records added** to census v2 (delta walk, "
         "`s2r4_delta.py`), swept by all profiles; **12 profiles added** (7 new census "
         "families; KChIP; the EMC3, GOST, BRI3BP and neuronal-calcium-sensor decoys), "
         "each swept over r3's NR database with `-Z` held at r3's size, so no r3 hit "
         f"moved. **{summary['r3_calls_changed']} of {summary['r3_records']:,} r3 calls "
         "changed**, every one with a new family as the new call, the winner or the "
         "runner-up that pulled the winner inside the D7 margin (checked, hard failure "
         "otherwise):", "", "| before | after | basis | profile winner | records |",
         "|---|---|---|---|---|"]
    L += [f"| {r['before']} | {r['after']} | {r['v3_basis']} | {r['profile_family'] or '—'} "
          f"| {r['records']} |" for r in tr]
    L += ["", "## The new records' calls", "", "| v3 call | records |", "|---|---|"]
    L += [f"| {k} | {v:,} |" for k, v in fams.most_common()]
    L += ["", "## The look-alike tests (H17–H20)", "",
          f"S3a's held-out orthologue benchmark (set B) now scores {b['positives']} "
          f"channel orthologues ({b['correct']} correct) and {b['decoy_ok'] + b['decoy_called_channel']} "
          f"decoys ({b['decoy_called_channel']} called a channel). Set B skips seeds and "
          "'Precursor'-flagged records, so every reviewed r4 record is listed here with "
          "its call (the gene symbol scores, never makes, the call — H15):", "",
          "| gene | organism | seed | profile call | confidence | margin | runner-up |",
          "|---|---|---|---|---|---|---|"]
    L += [f"| {r['gene'] or '—'} | {r['organism'][:28]} | {'yes' if r['seed'] == '1' else ''} "
          f"| {r['p_family'] or '—'} | {r['p_confidence']} | {r['rel_margin']} | "
          f"{r['runner'] or '—'} |" for r in sorted(rv, key=lambda r: (r['p_family'], r['gene']))]
    (OUT_DIR / "r4_report.md").write_text("\n".join(L) + "\n")


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--rev", default="r4")
    ap.add_argument("--involved", default="",
                    help="families whose profiles this revision changed (r6: "
                         "tmem87,nonchannel_gost,nonchannel_tmem87b)")
    a = ap.parse_args()
    rev = a.rev
    prev = {"r4": "r3", "r5": "r4", "r6": "r5"}[rev]
    # r4 added profiles, so an old call may move only with one involved; r5
    # changed none, so no earlier call may move at all; r6 names its own
    new_profiles = (set(NEW_PROFILES) if rev == "r4" else
                    set(filter(None, a.involved.split(","))))
    old = {r["accession"]: r for r in rows(s3_dir() / f"census_v3.{prev}.tsv.gz")}
    trans, new_fam, outside = Counter(), Counter(), []
    new_status = Counter()
    for r in rows(s3_dir() / "census_v3.tsv.gz"):
        o = old.get(r["accession"])
        if o is None:
            new_fam[(label(r), r["v3_basis"], r["p_confidence"])] += 1
            new_status[r["v3_status"]] += 1
            continue
        x, y = label(o), label(r)
        if x != y:
            trans[(x, y, r["v3_basis"], r["p_family"])] += 1
            # a new profile can move a call by winning or by becoming the
            # runner-up inside the D7 margin (→ ambiguous / superfamily_only)
            if not ({r["p_family"], r["v3_family"], o["p_family"], r["runner"],
                     o["runner"]} & new_profiles):     # a retired profile, too
                outside.append(r["accession"])
    if outside:
        raise SystemExit(f"{len(outside)} r3 calls changed with no new family "
                         f"involved, e.g. {outside[:5]}")
    write_tsv(OUT_DIR / f"{rev}_transitions.tsv",
              ["before", "after", "v3_basis", "profile_family", "records"],
              [{"before": a, "after": b, "v3_basis": c, "profile_family": d,
                "records": n} for (a, b, c, d), n in trans.most_common()])
    write_tsv(OUT_DIR / f"{rev}_new_records.tsv",
              ["v3_call", "v3_basis", "p_confidence", "records"],
              [{"v3_call": a, "v3_basis": b, "p_confidence": c, "records": n}
               for (a, b, c), n in new_fam.most_common()])
    summary = {f"{prev}_records": len(old), f"{prev}_calls_changed": sum(trans.values()),
               f"{rev}_records": sum(new_fam.values()), f"{rev}_status": dict(new_status)}
    (OUT_DIR / f"{rev}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if rev == "r4":
        reviewed_calls()
        render({"r3_records": len(old), "r3_calls_changed": sum(trans.values()),
                "r4_records": sum(new_fam.values())})
    print(json.dumps(summary))
    for (a, b, c, d), n in trans.most_common(12):
        print(f"  {a} → {b} ({c}, profile {d}): {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
