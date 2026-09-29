"""S2b — what the H2 / H4 / H13 rewrite changed, call by call (r1 → r2).

    python3 scripts/s2b_revision.py

Reads the archived r1 call files (`<data root>/raw_api/s2/calls_r1/`) and
the current ones, and counts every transition of (family, superfamily)
between the two revisions. Only records carrying `PF02931`, `PF08709`,
`PF08016` or `PF20519` were re-classified; everything else is identical by
construction, and this script checks that rather than assuming it.

Writes `results/census_v2/s2b_transitions.tsv` and `s2b_revision.json`.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s2_lib import OUT_DIR, SHARDS, raw_dir   # noqa: E402
from scripts.s3_hmm_lib import write_tsv              # noqa: E402

RECHECKED = {"PF02931", "PF08709", "PF08016", "PF20519"}


def rows(p: Path):
    with gzip.open(p, "rt", newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def label(r: dict) -> str:
    return r["family"] or (f"[{r['superfamily']}]" if r["superfamily"] else "unassigned")


def main() -> int:
    trans, n, outside_changed, rechecked = Counter(), 0, 0, 0
    for key, _ in SHARDS:
        old = {r["accession"]: r for r in rows(raw_dir() / "calls_r1" / f"{key}.tsv.gz")}
        for r in rows(raw_dir() / "calls" / f"{key}.tsv.gz"):
            n += 1
            o = old[r["accession"]]
            carries = bool({p.split(":")[0] for p in r["pfam"].split(";")} & RECHECKED)
            rechecked += carries
            a, b = label(o), label(r)
            if a != b:
                if not carries:
                    outside_changed += 1
                trans[(a, b, r["hazards"] or "-")] += 1
    if outside_changed:
        raise SystemExit(f"{outside_changed} calls changed outside the re-checked "
                         "set — the recheck was not the only change")
    write_tsv(OUT_DIR / "s2b_transitions.tsv",
              ["r1_call", "r2_call", "r2_hazards", "records"],
              [{"r1_call": a, "r2_call": b, "r2_hazards": h, "records": k}
               for (a, b, h), k in trans.most_common()])
    summary = {"records": n, "rechecked": rechecked,
               "changed": sum(trans.values()), "changed_outside_recheck": 0}
    (OUT_DIR / "s2b_revision.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
