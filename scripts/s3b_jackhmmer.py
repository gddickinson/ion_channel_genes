"""S3b step 2 — jackhmmer to convergence, one derived seed per census family.

    python3 scripts/s3b_jackhmmer.py seeds
    python3 scripts/s3b_jackhmmer.py run [--jobs 4 --cpu 2] [--only nav,cav]
    python3 scripts/s3b_jackhmmer.py parse

**The seed is derived, not chosen** (the parent project's S20 rule): for
each census family, the panel entry the S3b profile sweep calls to that
family at `high` confidence with the highest full-sequence score, ties
broken by target name. `seeds` writes `jackhmmer_seeds.tsv` before anything
runs, so a changed sweep shows up as a changed seed. A family with no
high-confidence call on the panel has no seed and no run — a missing result
with its reason (D28), never a substitute seed.

`run` is one `jackhmmer -N 10 -E 1e-5 --incE 1e-5` per seed over the panel
DB, logs archived under `<data root>/hmmer/s3b/jackhmmer/`; a run is reused
only if its `.done` marker records the same seed and DB SHA-256. `parse`
re-derives every table from the archived logs (no search), applies D10
(`s3b_kill.py`) and writes the per-run, per-round and completeness tables
plus the accepted target sets the census merge reads — from clean runs
and from the pre-kill rounds of K1/K2 runs, never from a K3 run (D36).
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_hmm_lib import LIVE, iter_fasta, read_tsv, write_tsv  # noqa: E402
from scripts.s3b_kill import MAX_ITER, accepted_targets, evaluate  # noqa: E402
from scripts.s3b_lib import (  # noqa: E402
    OUT_DIR, panel_db, parse_jackhmmer_log, s3b_dir,
)
from scripts.s3b_sweep import calls_path, universe_path  # noqa: E402
from src.catalogue import registry  # noqa: E402

EVALUE = "1e-5"
SEED_FIELDS = ["family", "superfamily", "seed", "species", "length",
               "win_score", "win_coverage", "rel_margin", "status"]
RUN_FIELDS = ["family", "seed", "rounds", "converged", "verdict", "rule",
              "killed_at", "accepted_rounds", "n_accepted", "seconds", "reason"]
ROUND_FIELDS = ["family", "round", "new_targets", "n_included", "n_own",
                "n_other_family", "n_uncalled", "other_frac", "other_rise",
                "growth", "n_other_superfamily", "other_sf_frac",
                "uncalled_frac"]
COMPLETE_FIELDS = ["family", "superfamily", "verdict", "profile_calls",
                   "recovered", "recall", "sf_profile_calls", "sf_recovered",
                   "sf_recall", "uncalled_frac", "jh_accepted", "jh_other_family",
                   "jh_own_superfamily_only", "jh_module", "jh_ambiguous",
                   "jh_low_score", "jh_no_hit", "top_other_families"]


def jh_dir() -> Path:
    return s3b_dir("jackhmmer")


def load_calls() -> dict[str, dict]:
    return {r["target"]: r for r in read_tsv(calls_path())}


def seeds() -> list[dict]:
    calls = load_calls()
    uni = {r["target"]: r for r in read_tsv(universe_path())}
    best: dict[str, dict] = {}
    for t, c in calls.items():
        if c["p_call"] != "family" or c["p_confidence"] != "high":
            continue
        f = c["p_family"]
        key = (float(c["win_score"]), t)
        if f not in best or key > (float(best[f]["win_score"]), best[f]["target"]):
            best[f] = c
    rows = []
    for fam in registry.census_families():
        c = best.get(fam.key)
        if c is None:
            rows.append({"family": fam.key, "superfamily": fam.superfamily,
                         "status": "no_seed: no high-confidence profile call "
                                   "on the panel"})
            continue
        u = uni[c["target"]]
        rows.append({"family": fam.key, "superfamily": fam.superfamily,
                     "seed": c["target"], "species": u["species"],
                     "length": u["length"], "win_score": c["win_score"],
                     "win_coverage": c["win_coverage"],
                     "rel_margin": c["rel_margin"], "status": "seeded"})
    write_tsv(OUT_DIR / "jackhmmer_seeds.tsv", SEED_FIELDS, rows)
    n = sum(r["status"] == "seeded" for r in rows)
    print(f"[seeds] {n}/{len(rows)} census families seeded")
    want = {r["seed"] for r in rows if r["status"] == "seeded"}
    got = {}
    for head, seq in iter_fasta(panel_db()):
        t = head.split(maxsplit=1)[0]
        if t in want:
            got[t] = (head, seq)
    for r in rows:
        if r["status"] == "seeded":
            head, seq = got[r["seed"]]
            (jh_dir() / f"{r['family']}.seed.faa").write_text(
                f">{head}\n{seq}\n")
    return rows


def _run_one(fam: str, seed: str, cpu: int, db_sha: str) -> dict:
    d = jh_dir()
    log, done = d / f"{fam}.log", d / f"{fam}.done"
    stamp = {"seed": seed, "db_sha256": db_sha, "evalue": EVALUE,
             "max_iter": MAX_ITER}
    if done.exists():
        rec = json.loads(done.read_text())
        if all(rec.get(k) == v for k, v in stamp.items()):
            return {**rec, "reused": True}
    cmd = ["jackhmmer", "--cpu", str(cpu), "-N", str(MAX_ITER), "-E", EVALUE,
           "--incE", EVALUE, "--noali", "-o", str(log),
           str(d / f"{fam}.seed.faa"), str(panel_db())]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"{fam}: jackhmmer exit {r.returncode}: {r.stderr[:300]}")
    rec = {**stamp, "family": fam, "seconds": round(time.time() - t0, 1),
           "command": " ".join(cmd[:-2] + ["<seed>", "<db>"])}
    done.write_text(json.dumps(rec, indent=2))
    return rec


def run(jobs: int, cpu: int, only: set[str]) -> None:
    db_sha = json.loads((OUT_DIR / "sweep_db.json").read_text())["db_sha256"]
    todo = [r for r in read_tsv(OUT_DIR / "jackhmmer_seeds.tsv")
            if r["status"] == "seeded" and (not only or r["family"] in only)]
    t0 = time.time()
    with ThreadPoolExecutor(jobs) as ex:
        futs = {ex.submit(_run_one, r["family"], r["seed"], cpu, db_sha): r
                for r in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            print(f"[run] {i}/{len(todo)} {rec['family']:22s} {rec['seconds']:7.1f} s"
                  f"{' (reused)' if rec.get('reused') else ''}"
                  f"  [{time.time() - t0:.0f} s]", flush=True)
            s0_lib.live_progress(LIVE, "S3b", [
                ("panel sweep", True), ("assign", True),
                (f"jackhmmer {i}/{len(todo)} families", i == len(todo)),
                ("census v3", False)], workers=jobs)


def _verdict_class(c: dict | None, own: str, own_sf: str) -> str:
    if c is None:
        return "no_hit"
    if c["p_call"] == "family":
        return "own" if c["p_family"] == own else "other_family"
    if c["p_call"] == "superfamily_only":
        return "own_superfamily_only" if c["p_superfamily"] == own_sf \
            else "ambiguous"
    return {"module": "module", "ambiguous": "ambiguous",
            "low_score": "low_score"}.get(c["p_call"], "no_hit")


def parse() -> None:
    calls = load_calls()
    family_of = {t: (c["p_family"] if c["p_call"] == "family" else "")
                 for t, c in calls.items()}
    sf_call = {t: (c["p_superfamily"] if c["p_call"] in
                   ("family", "superfamily_only") else "")
               for t, c in calls.items()}
    by_family: dict[str, set[str]] = {}
    by_sf: dict[str, set[str]] = {}
    for t, f in family_of.items():
        if f:
            by_family.setdefault(f, set()).add(t)
    for t, sfk in sf_call.items():
        if sfk:
            by_sf.setdefault(sfk, set()).add(t)
    runs, rounds_out, complete = [], [], []
    accepted_rows = []
    for s in read_tsv(OUT_DIR / "jackhmmer_seeds.tsv"):
        fam, sf = s["family"], s["superfamily"]
        prof = by_family.get(fam, set())
        base = {"family": fam, "superfamily": sf,
                "profile_calls": len(prof)}
        log, done = jh_dir() / f"{fam}.log", jh_dir() / f"{fam}.done"
        if s["status"] != "seeded" or not done.exists():
            why = s["status"] if s["status"] != "seeded" else "not run"
            runs.append({"family": fam, "seed": s.get("seed", ""),
                         "verdict": "no_run", "reason": why})
            complete.append({**base, "verdict": "no_run"})
            continue
        conv = parse_jackhmmer_log(log)
        ev = evaluate(conv["rounds"], family_of, fam, conv["converged"])
        acc = accepted_targets(conv["rounds"], ev["accepted_rounds"])
        runs.append({"family": fam, "seed": s["seed"],
                     "rounds": len(conv["rounds"]),
                     "converged": conv["converged"], "verdict": ev["verdict"],
                     "rule": ev["rule"], "killed_at": ev["killed_at"],
                     "accepted_rounds": ev["accepted_rounds"],
                     "n_accepted": len(acc),
                     "seconds": json.loads(done.read_text())["seconds"],
                     "reason": ev["reason"]})
        # Reported, never applied (D10's rules are fixed): did the run
        # leave its superfamily, and how much of it no profile calls at all?
        for pr, rd in zip(ev["per_round"], conv["rounds"]):
            inc = rd["included"]
            o = sum(1 for t in inc if sf_call.get(t) and sf_call[t] != sf)
            pr["n_other_superfamily"] = o
            pr["other_sf_frac"] = round(o / len(inc), 4) if inc else 0.0
            pr["uncalled_frac"] = round(pr["n_uncalled"] / len(inc), 4) \
                if inc else 0.0
        rounds_out += [{"family": fam, **pr} for pr in ev["per_round"]]
        sfp = by_sf.get(sf, set())
        cls = Counter(_verdict_class(calls.get(t), fam, sf) for t in acc)
        others = Counter(family_of[t] for t in acc
                         if family_of.get(t) and family_of[t] != fam)
        rec = len(prof & acc)
        complete.append({
            **base, "verdict": ev["verdict"], "recovered": rec,
            "recall": round(rec / len(prof), 4) if prof else "",
            "sf_profile_calls": len(sfp), "sf_recovered": len(sfp & acc),
            "sf_recall": round(len(sfp & acc) / len(sfp), 4) if sfp else "",
            "uncalled_frac": round((cls["module"] + cls["low_score"]
                                    + cls["no_hit"] + cls["ambiguous"])
                                   / len(acc), 4) if acc else "",
            "jh_accepted": len(acc), "jh_other_family": cls["other_family"],
            "jh_own_superfamily_only": cls["own_superfamily_only"],
            "jh_module": cls["module"], "jh_ambiguous": cls["ambiguous"],
            "jh_low_score": cls["low_score"], "jh_no_hit": cls["no_hit"],
            "top_other_families": ",".join(f"{k}:{v}"
                                           for k, v in others.most_common(3))})
        # A K3 run never converged: its rounds are reported, and it
        # contributes no candidates (D36). K1/K2 runs contribute the rounds
        # before the rule fired; clean runs everything.
        if ev["rule"] != "K3":
            accepted_rows += [{"family": fam, "target": t,
                               "profile_class": _verdict_class(calls.get(t), fam, sf)}
                              for t in sorted(acc)]
    write_tsv(OUT_DIR / "jackhmmer_runs.tsv", RUN_FIELDS, runs)
    write_tsv(OUT_DIR / "jackhmmer_rounds.tsv", ROUND_FIELDS, rounds_out)
    write_tsv(OUT_DIR / "jackhmmer_completeness.tsv", COMPLETE_FIELDS, complete)
    write_tsv(s3b_dir() / "jackhmmer_accepted.tsv.gz",
              ["family", "target", "profile_class"], accepted_rows)
    v = Counter(r["verdict"] for r in runs)
    print(f"[parse] {dict(v)}; {len(accepted_rows):,} accepted (family, target)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["seeds", "run", "parse"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    if a.step == "seeds":
        seeds()
    elif a.step == "run":
        run(a.jobs, a.cpu, set(filter(None, a.only.split(","))))
    else:
        parse()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
