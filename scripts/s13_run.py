"""s13_run.py — run S13's job list, resumably, and say how far it has got.

    nohup caffeinate -i python3 -u scripts/s13_run.py --workers 7 \
        > <data root>/selection/s13/run.log 2>&1 &
    python3 scripts/s13_run.py --status

Jobs come from `s13_jobs.all_jobs()`, largest first. A finished job (its
output parses) is skipped and a claimed one left to its owner, so the same
command resumes an interrupted suite. `--only PREFIX,…` restricts by job
name prefix (e.g. to send the RyR branch-site runs to another machine).
Progress goes to `results/session_live.json` for the dashboard.
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s0_lib import live_progress  # noqa: E402
from s13_codeml import run  # noqa: E402
from s13_jobs import all_jobs  # noqa: E402
from s13_lib import LIVE  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--hyphy-threads", type=int, default=2)
    ap.add_argument("--only", default="")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    jobs = all_jobs()
    if a.only:
        jobs = [j for j in jobs if any(j.name.startswith(p) for p in a.only.split(","))]
    if a.status:
        for j in jobs:
            state = "done" if j.done() else ("running" if (j.workdir / "RUNNING").exists() else "-")
            print(f"{state:8s} {j.name}")
        print(f"{sum(j.done() for j in jobs)}/{len(jobs)} done")
        return
    state = {j.name: j.done() for j in jobs}

    def progress() -> None:
        live_progress(LIVE, "S13a", [(n, d) for n, d in state.items()], workers=a.workers)

    progress()
    t0 = time.time()
    todo = [j for j in jobs if not state[j.name]]
    print(f"{len(jobs)} jobs, {len(todo)} to run, {a.workers} workers", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(run, j, a.hyphy_threads): j for j in todo}
        for f in as_completed(futs):
            j = futs[f]
            try:
                verdict = f.result()
            except Exception as exc:  # noqa: BLE001
                verdict = f"error: {exc}"
            state[j.name] = j.done()
            print(f"[{(time.time() - t0) / 3600:6.2f} h] {verdict:8s} {j.name}", flush=True)
            progress()
    failed = [n for n, d in state.items() if not d]
    print(f"finished: {len(jobs) - len(failed)}/{len(jobs)} done; failed {failed}", flush=True)


if __name__ == "__main__":
    main()
