"""s4b_sweep.py — the S3a profile library over S4b's dense panel (one proteome per order).

    python3 scripts/s4b_sweep.py search [--jobs 4 --cpu 2]
    python3 scripts/s4b_sweep.py assign
    python3 scripts/s4b_sweep.py matrix          # → results/panel_density/order_matrix.tsv

Same instrument as S3b (`s3_sweep._search`, `s3_assign.assign_all`, D32),
pointed at `<data root>/proteomes/s4b/dense_panel.fasta` with `-Z` = its size,
resumable on each profile's SHA-256. The matrix is census family × order:
high-confidence and any-confidence profile calls per proteome. Every cell is
a **proteome** call — S5's genome control exists only for the 52 S4 species
(D37/D38); an empty cell here is an annotation-level absence, never a
controlled one.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_assign import ASSIGN_FIELDS, assign_all, collect_hits  # noqa: E402
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s3_sweep import _search  # noqa: E402
from s4b_dense_panel import OUT, bulk  # noqa: E402
from src.catalogue import registry  # noqa: E402


def cmd_search(a) -> None:
    n = json.loads((OUT / "dense_db.json").read_text())["sequences"]
    db = bulk() / "dense_panel.fasta"
    profs = [f.key for f in registry.families()]
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(_search, p, db, n, a.cpu, bulk("domtbl")): p for p in profs}
        for i, f in enumerate(as_completed(futs), 1):
            r = f.result()
            print(f"[search] {i}/{len(profs)} {r['profile']:26s} {r['targets_reported']:8,d} "
                  f"{r['seconds']:7.1f} s{' (reused)' if r.get('reused') else ''}", flush=True)


def cmd_assign(_a) -> None:
    paths = sorted(bulk("domtbl").glob("*.domtbl.gz"))
    if len(paths) != len(registry.families()):
        raise SystemExit(f"{len(paths)} domtbls for {len(registry.families())} profiles")
    calls = assign_all(collect_hits(paths))
    n = write_tsv(bulk() / "calls.tsv.gz", ASSIGN_FIELDS, (calls[t] for t in sorted(calls)))
    print(f"[assign] {n:,} entries with ≥ 1 profile hit")


def cmd_matrix(_a) -> None:
    uni = {}
    with gzip.open(bulk() / "universe.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            t, upid, order, cls, phylum, kingdom = line.rstrip("\n").split("\t")
            uni[t] = (upid, order)
    census = {f.key for f in registry.census_families()}
    cell = defaultdict(Counter)
    for r in read_tsv(bulk() / "calls.tsv.gz"):
        if r["p_call"] == "family" and r["p_family"] in census and r["target"] in uni:
            upid, order = uni[r["target"]]
            cell[(order, r["p_family"])]["any"] += 1
            if r["p_confidence"] == "high":
                cell[(order, r["p_family"])]["high"] += 1
    panel = read_tsv(OUT / "order_panel.tsv")
    rows = []
    for p in panel:
        for fam in sorted(census):
            c = cell.get((p["order"], fam), Counter())
            rows.append({"order": p["order"], "class": p["class"], "phylum": p["phylum"],
                         "kingdom": p["kingdom"], "upid": p["upid"], "organism": p["organism"],
                         "family": fam, "high": c["high"], "any": c["any"],
                         "present": int(c["high"] > 0)})
    write_tsv(OUT / "order_matrix.tsv", list(rows[0]), rows)
    fam_orders = Counter(r["family"] for r in rows if r["present"])
    summary = {"orders": len(panel), "families": len(census),
               "present_cells": sum(r["present"] for r in rows),
               "families_by_orders_present": dict(fam_orders.most_common())}
    (OUT / "order_matrix.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "families_by_orders_present"}))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["search", "assign", "matrix"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpu", type=int, default=2)
    a = ap.parse_args()
    {"search": cmd_search, "assign": cmd_assign, "matrix": cmd_matrix}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
