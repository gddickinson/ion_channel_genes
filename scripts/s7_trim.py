"""s7_trim.py — S7 step 1: measure the trimming candidates before any tree.

    python3 scripts/s7_trim.py            # → results/phylogeny/tier1_trim_compare.tsv

S6 trimmed every family with trimAl `-automated1` and found it keeps 3–4 %
of columns on the largest divergent families (K2P 145 / 5,149). The trimming
used for the census trees is chosen here, on alignment properties only —
columns, parsimony-informative sites, gap fraction and the share of each
member's own residues retained — and fixed (D41) before any tier-1 tree
exists, so no tree is read to choose it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s7_lib import (OUT_DIR, ROOT, TRIM_METHODS, family_alignment,  # noqa: E402
                    metrics, read_aln, trim_columns)

FIELDS = ["family", "n_seqs", "method", "cols", "informative", "gap_frac",
          "retention_median", "retention_min"]


def main() -> int:
    fams = [r["family"] for r in read_tsv(ROOT / "results" / "alignments" /
                                          "alignments.tsv")
            if r["status"] == "aligned"]
    out = []
    for fam in fams:
        aln = family_alignment(fam)
        rows = read_aln(aln)
        for method in TRIM_METHODS:
            m = metrics(rows, trim_columns(aln, method))
            out.append({"family": fam, "n_seqs": len(rows), "method": method, **m})
        print(fam, len(rows), " ".join(f"{r['method']}={r['cols']}/{r['informative']}"
                                      for r in out[-len(TRIM_METHODS):]), flush=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT_DIR / "tier1_trim_compare.tsv", FIELDS, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
