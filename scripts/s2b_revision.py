"""S2b / S2c — what a hazard-rule rewrite changed, call by call.

    python3 scripts/s2b_revision.py --head r2            # S2b: r1 → r2
    python3 scripts/s2b_revision.py --tag s2c --base r2 \
        --accessions PF02214,PF00664                      # S2c: r2 → r3

Reads the archived call files of the base revision
(`<data root>/raw_api/s2/calls_<base>/`) and the current ones, and counts
every transition of (family, superfamily). Only records carrying one of the
re-checked accessions were re-classified; everything else is identical by
construction, and this script checks that rather than assuming it.

Writes `results/census_v2/<tag>_transitions.tsv` and `<tag>_revision.json`.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s2_lib import OUT_DIR, SHARDS, raw_dir   # noqa: E402
from scripts.s3_hmm_lib import write_tsv              # noqa: E402

S2B_RECHECKED = "PF02931,PF08709,PF08016,PF20519"


def rows(p: Path):
    with gzip.open(p, "rt", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def label(r: dict) -> str:
    return r["family"] or (f"[{r['superfamily']}]" if r["superfamily"] else "unassigned")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="s2b")
    ap.add_argument("--base", default="r1")
    ap.add_argument("--accessions", default=S2B_RECHECKED)
    ap.add_argument("--head", default="",
                    help="archived revision to compare to (default: current calls)")
    a = ap.parse_args()
    RECHECKED = set(a.accessions.split(","))
    trans, n, outside_changed, rechecked = Counter(), 0, 0, 0
    for key, _ in SHARDS:
        old = {r["accession"]: r for r in rows(raw_dir() / f"calls_{a.base}" / f"{key}.tsv.gz")}
        head = raw_dir() / (f"calls_{a.head}" if a.head else "calls")
        for r in rows(head / f"{key}.tsv.gz"):
            n += 1
            o = old[r["accession"]]
            carries = bool({p.split(":")[0] for p in r["pfam"].split(";")} & RECHECKED)
            rechecked += carries
            x, y = label(o), label(r)
            if x != y:
                if not carries:
                    outside_changed += 1
                trans[(x, y, r["hazards"] or "-")] += 1
    if outside_changed:
        raise SystemExit(f"{outside_changed} calls changed outside the re-checked "
                         "set — the recheck was not the only change")
    write_tsv(OUT_DIR / f"{a.tag}_transitions.tsv",
              ["before", "after", "hazards_after", "records"],
              [{"before": x, "after": y, "hazards_after": h, "records": k}
               for (x, y, h), k in trans.most_common()])
    summary = {"base": a.base, "head": a.head or "current",
               "accessions": sorted(RECHECKED),
               "records": n, "rechecked": rechecked,
               "changed": sum(trans.values()), "changed_outside_recheck": 0}
    (OUT_DIR / f"{a.tag}_revision.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
