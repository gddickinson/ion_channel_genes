"""S3a — census v3a against the parent projects' censuses (an external check).

    python3 scripts/s3_external_check.py

The IP3R project (`../ip3r_genes`) and the PIEZO project (`../piezo_genes`)
each ran a census of families this catalogue also holds, with their own
seeds, their own profiles and their own search spaces (vertebrate reference
proteomes plus genome sweeps; UniProt releases of their own dates). CLAUDE.md:
"Port the method; never port a result" — so nothing here feeds a call. The
parent censuses are read only to **compare**, on the UniProt accessions both
contain, and to count what each holds that the other does not.

What is compared, and why it is not symmetrical:

* **IP3R census v6** carries a positive per-record call (`ITPR` / `RYR` /
  `unassigned` / `conflict`) from two instruments, so it is compared call
  against call: parent ITPR → our `itpr`, parent RYR → our `ryr`.
* **PIEZO census v5** is a *membership* list ("piezo-like": Pfam carriers
  plus profile and jackhmmer hits, fragments included) with no per-record
  family call. It is compared as membership: how many of its UniProt
  accessions are in census v2 at all, and what census v3a calls them.

Rows whose accession is not a UniProt accession (Ensembl/Compara proteins,
genome-derived gene models) cannot be matched and are counted as such.

Writes `external_check.tsv`, `external_reverse.tsv`,
`external_disagreements.tsv` and `external_sources.tsv` to
`results/census_v3/`.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import OUT_DIR, read_tsv, s3_dir, sha256, write_tsv  # noqa: E402

IP3R_V6 = ROOT.parent / "ip3r_genes" / "results" / "census_v6" / "census_v6.tsv"
PIEZO_V5 = ROOT.parent / "piezo_genes" / "results" / "census_v5" / "piezo_like_census.csv"
# UniProt accession grammar (uniprot.org/help/accession_numbers).
UNIPROT_ACC = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
IP3R_MAP = {"ITPR": "itpr", "RYR": "ryr"}
OURS = ("itpr", "ryr", "piezo")


def parent_rows() -> list[dict]:
    rows = []
    for r in read_tsv(IP3R_V6):
        rows.append({"project": "ip3r_genes/census_v6", "accession": r["accession"],
                     "parent_call": r["call"], "parent_source": r["source"],
                     "parent_reason": r["reason"]})
    with PIEZO_V5.open(newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append({"project": "piezo_genes/census_v5", "accession": r["accession"],
                         "parent_call": "piezo_like", "parent_source": r["source"],
                         "parent_reason": r["description"][:120]})
    return rows


def ours_index(wanted: set[str]) -> tuple[dict[str, dict], dict[str, tuple]]:
    """Census v3a rows for the wanted accessions + our own itpr/ryr/piezo calls."""
    idx, ours = {}, {}
    for r in read_tsv(s3_dir() / "census_v3.tsv.gz"):
        a = r["accession"]
        if a in wanted:
            idx[a] = r
        if r["v3_family"] in OURS:
            ours[a] = (r["v3_family"], r["v3_basis"])
    return idx, ours


def verdict(parent: str, ours: dict | None) -> str:
    if ours is None:
        return "not_in_census_v2"
    fam = ours["v3_family"]
    if parent == "piezo_like":
        return ("called_piezo" if fam == "piezo"
                else f"called_other_family" if fam
                else f"v3_{ours['v3_basis']}")
    want = IP3R_MAP.get(parent)
    if want is None:                                  # parent unassigned / conflict
        return f"parent_{parent}__ours_{fam or ours['v3_basis']}"
    if fam == want:
        return "agree"
    if fam in IP3R_MAP.values():
        return "swapped_itpr_ryr"
    return "called_other_family" if fam else f"v3_{ours['v3_basis']}"


def main() -> int:
    prow = parent_rows()
    uni = [r for r in prow if UNIPROT_ACC.match(r["accession"])]
    idx, ours = ours_index({r["accession"] for r in uni})

    cnt, dis = Counter(), []
    for r in prow:
        if not UNIPROT_ACC.match(r["accession"]):
            cnt[(r["project"], r["parent_call"], "not_a_uniprot_accession")] += 1
            continue
        o = idx.get(r["accession"])
        v = verdict(r["parent_call"], o)
        cnt[(r["project"], r["parent_call"], v)] += 1
        if o is not None and v not in ("agree", "called_piezo") \
                and not v.startswith("parent_"):
            dis.append({**r, "verdict": v, "v2_family": o["v2_family"],
                        "v2_tier": o["v2_tier"], "p_call": o["p_call"],
                        "p_family": o["p_family"], "win_score": o["win_score"],
                        "win_coverage": o["win_coverage"], "runner": o["runner"],
                        "rel_margin": o["rel_margin"], "v3_family": o["v3_family"],
                        "v3_basis": o["v3_basis"], "length": o["length"],
                        "fragment": o["fragment"]})

    write_tsv(OUT_DIR / "external_check.tsv",
              ["project", "parent_call", "verdict", "records"],
              [{"project": p, "parent_call": c, "verdict": v, "records": n}
               for (p, c, v), n in sorted(cnt.items())])
    # Reverse direction: of our own itpr / ryr / piezo calls, how many does
    # the parent census hold at all (IP3R for itpr/ryr, PIEZO for piezo)?
    held = {"itpr": {r["accession"] for r in uni if r["project"].startswith("ip3r")},
            "piezo": {r["accession"] for r in uni if r["project"].startswith("piezo")}}
    held["ryr"] = held["itpr"]
    rev = Counter((fam, basis, a in held[fam]) for a, (fam, basis) in ours.items())
    write_tsv(OUT_DIR / "external_reverse.tsv",
              ["our_family", "v3_basis", "in_parent_census", "records"],
              [{"our_family": f, "v3_basis": b, "in_parent_census": "yes" if h else "no",
                "records": n} for (f, b, h), n in sorted(rev.items())])
    write_tsv(OUT_DIR / "external_disagreements.tsv",
              ["project", "accession", "parent_call", "parent_source", "verdict",
               "v2_family", "v2_tier", "p_call", "p_family", "win_score",
               "win_coverage", "runner", "rel_margin", "v3_family", "v3_basis",
               "length", "fragment", "parent_reason"],
              sorted(dis, key=lambda r: (r["project"], r["verdict"], r["accession"])))
    write_tsv(OUT_DIR / "external_sources.tsv", ["file", "rows", "sha256"],
              [{"file": str(p.relative_to(ROOT.parent)), "rows": n, "sha256": sha256(p)}
               for p, n in ((IP3R_V6, sum(r["project"].startswith("ip3r") for r in prow)),
                            (PIEZO_V5, sum(r["project"].startswith("piezo") for r in prow)))])
    print(f"[external] {len(prow):,} parent rows, {len(uni):,} UniProt accessions, "
          f"{len(idx):,} in census v2; {len(dis):,} disagreements")
    for k, n in sorted(cnt.items()):
        print(f"  {k[0]:24s} {k[1]:12s} {k[2]:40s} {n:7,d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
