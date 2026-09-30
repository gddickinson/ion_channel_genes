"""s3_piezo_unassigned.py — why census v3a leaves 2,258 PIEZO-census records unassigned.

    python3 scripts/s3_piezo_unassigned.py

S3a's external check found 2,258 records the parent PIEZO project's census
(v5) counts as Piezo-like that census v3a calls nothing. Each is looked up in
census v3a (profile outcome, winning score and coverage, length, fragment
flag) and binned by *why* it has no call: below D32's score gate, below its
coverage gate (a partial Piezo), inside the margin, no hit at all. Read-only
comparison (the parent census is never an input, D0).

→ `results/census_v3/piezo_unassigned.tsv` (bins) and
`piezo_unassigned_records.tsv` (every record).
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import OUT_DIR, s3_dir, write_tsv  # noqa: E402

PARENT = ROOT.parent / "piezo_genes" / "results" / "census_v5" / "piezo_like_census.csv"
PIEZO_LEN = 2000        # a complete Piezo is ~2,100–2,900 aa


def main() -> int:
    parent = {r["accession"]: r for r in csv.DictReader(open(PARENT))}
    rows = []
    with gzip.open(s3_dir() / "census_v3.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] in parent and r["v3_status"] == "unassigned":
                length = int(r["length"] or 0)
                frag = "Fragment" in (r["fragment"] or "")
                why = {"low_score": "below the 30-bit score gate",
                       "module": "below the 30 % coverage gate (partial)",
                       "no_hit": "no profile hit",
                       "ambiguous": "inside the margin, two superfamilies"}.get(
                    r["p_call"], r["p_call"] or "no profile hit")
                rows.append({"accession": r["accession"], "length": length,
                             "fragment": int(frag), "p_call": r["p_call"],
                             "p_family": r["p_family"], "win_score": r["win_score"],
                             "win_coverage": r["win_coverage"], "reason": why,
                             "size": ("fragment-flagged" if frag else
                                      "< 500 aa" if length < 500 else
                                      "500–2,000 aa" if length < PIEZO_LEN else
                                      "≥ 2,000 aa"),
                             "parent_source": parent[r["accession"]]["source"]})
    write_tsv(OUT_DIR / "piezo_unassigned_records.tsv", list(rows[0]), rows)
    c = Counter((r["reason"], r["size"]) for r in rows)
    write_tsv(OUT_DIR / "piezo_unassigned.tsv", ["reason", "size", "records"],
              [{"reason": a, "size": b, "records": n} for (a, b), n in c.most_common()])
    print(len(rows))
    for (a, b), n in c.most_common():
        print(f"  {n:5d}  {a}  |  {b}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
