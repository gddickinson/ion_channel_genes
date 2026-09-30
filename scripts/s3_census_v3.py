"""S3a step 4 — merge S2's calls with the profile calls → census v3a.

    python3 scripts/s3_census_v3.py

Every census v2 record gains a profile verdict (via its sequence id), and
the merged call is a stated rule over the two instruments, which read
different evidence (Pfam architecture + filter motif vs a full-length
profile margin):

* both make the same family call               → that family, basis `both`
* S2 makes a family call, the profile does not → S2's call, `s2_only`
* the profile makes one, S2 does not           → the profile's, `profile_only`
  (only if it lies in the superfamily S2 named, when S2 named one)
* they name different families, or the profile's family lies outside the
  superfamily S2 named                         → `conflict`, kept and counted
* neither makes a family call → the narrower superfamily either offers,
  `superfamily_only`, else `unassigned`

Bulk: `<data root>/hmmer/s3/census_v3.tsv.gz` + `profile_calls.tsv.gz`.
Committed tables → `results/census_v3/` (rendered by `s3_report.py`, D13).
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_assign import ASSIGN_FIELDS, assign_all, collect_hits  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    LIVE, OUT_DIR, iter_census_v2, read_tsv, s3_dir, sha256, write_tsv,
)
from src.catalogue import registry  # noqa: E402

V2_COLS = ("accession", "reviewed", "taxon_id", "domain", "length", "fragment",
           "family", "superfamily", "status", "decisive_tier", "gene")
V3_FIELDS = ["accession", "reviewed", "taxon_id", "domain", "length",
             "fragment", "v2_family", "v2_superfamily", "v2_status",
             "v2_tier", "p_call", "p_family", "p_superfamily", "p_confidence",
             "win_score", "win_coverage", "runner", "rel_margin",
             "v3_family", "v3_superfamily", "v3_status", "v3_basis"]


def merge(v2: dict, p: dict, fam_status: dict[str, str]) -> dict:
    f2, sf2 = v2["family"], v2["superfamily"]
    fp = p.get("p_family", "") if p.get("p_call") == "family" else ""
    sfp = p.get("p_superfamily", "")
    if f2 and fp:
        if f2 == fp:
            return {"v3_family": f2, "v3_superfamily": sf2 or sfp, "v3_basis": "both"}
        return {"v3_family": "", "v3_superfamily": "", "v3_basis": "conflict"}
    if f2:
        return {"v3_family": f2, "v3_superfamily": sf2, "v3_basis": "s2_only"}
    if fp:
        if sf2 and sfp != sf2:
            return {"v3_family": "", "v3_superfamily": "", "v3_basis": "conflict"}
        return {"v3_family": fp, "v3_superfamily": sfp, "v3_basis": "profile_only"}
    sf = sf2 or (sfp if p.get("p_call") == "superfamily_only" else "")
    return {"v3_family": "", "v3_superfamily": sf,
            "v3_basis": "superfamily_only" if sf else "unassigned"}


def v3_status(row: dict, fam_status: dict[str, str]) -> str:
    if row["v3_basis"] == "conflict":
        return "conflict"
    if row["v3_family"]:
        return fam_status[row["v3_family"]]
    return "superfamily_only" if row["v3_superfamily"] else "unassigned"


def main() -> int:
    s0_lib.live_progress(LIVE, "S3a", [("profiles built", True),
                                       ("census sweep", True),
                                       ("assigning", False)])
    paths = sorted(s3_dir("domtbl").glob("*.domtbl.gz"))
    n_prof = len(read_tsv(OUT_DIR / "profile_build.tsv"))
    if len(paths) != n_prof:
        raise SystemExit(f"{len(paths)} domtblouts for {n_prof} profiles — "
                         "run s3_sweep.py search to completion first")
    # Census v2 r4 (D43): the delta's own sequences, searched by every
    # profile with -Z held at r3's size (s3r4_sweep.py). A target sequence
    # is in exactly one of the two databases, so the hit lists never overlap.
    deltas = []
    for d in sorted(s3_dir().glob("domtbl_r*")):          # r4, r5 … (D43)
        got = sorted(d.glob("*.domtbl.gz"))
        if len(got) != n_prof:
            raise SystemExit(f"{len(got)} {d.name} domtblouts for {n_prof} profiles — "
                             "run s3r4_sweep.py --rev … search to completion first")
        deltas += got
    top = collect_hits(paths + deltas)
    calls = assign_all(top)
    print(f"[assign] {len(calls):,} sequences with ≥ 1 profile hit")
    pc = s3_dir() / "profile_calls.tsv.gz"
    write_tsv(pc, ASSIGN_FIELDS, (calls[t] for t in sorted(calls)))

    sid = {}
    for mp in ["census_v2.nr_map.tsv.gz"] + sorted(
            q.name for q in s3_dir().glob("census_v2.r*_delta.nr_map.tsv.gz")):
        if not (s3_dir() / mp).exists():
            continue
        with gzip.open(s3_dir() / mp, "rt") as fh:
            next(fh)
            for line in fh:
                a, s, _ = line.rstrip("\n").split("\t")
                sid[a] = s
    fam_status = {f.key: f.status.value for f in registry.families()}
    seeds = {r["accession"] for r in read_tsv(OUT_DIR / "seed_manifest.tsv")}
    no_hit = {"p_call": "no_hit", "p_confidence": "none"}

    basis, status, pcall = Counter(), Counter(), Counter()
    fam = defaultdict(Counter)
    calib = defaultdict(Counter)            # v2 family → agreement, seeds excluded
    calib_dom = defaultdict(Counter)        # the same, by taxonomic domain
    conflicts, sfo, unas = Counter(), Counter(), Counter()
    support = Counter()                     # what each v3 family call rests on
    human = {}

    def rows():
        for v2 in iter_census_v2(V2_COLS):
            p = calls.get(sid[v2["accession"]], no_hit)
            m = merge(v2, p, fam_status)
            row = {**{k: v2[k] for k in ("accession", "reviewed", "taxon_id",
                                          "domain", "length", "fragment")},
                   "v2_family": v2["family"], "v2_superfamily": v2["superfamily"],
                   "v2_status": v2["status"], "v2_tier": v2["decisive_tier"],
                   **{k: p.get(k, "") for k in ("p_call", "p_family",
                      "p_superfamily", "p_confidence", "win_score",
                      "win_coverage", "runner", "rel_margin")}, **m}
            row["v3_status"] = v3_status(row, fam_status)
            basis[row["v3_basis"]] += 1
            status[row["v3_status"]] += 1
            pcall[(v2["status"], row["p_call"])] += 1
            if row["v3_family"]:
                support[(row["v3_family"], row["v3_basis"], v2["decisive_tier"] or "-",
                         row["p_call"])] += 1
                fam[row["v3_family"]][row["v3_basis"]] += 1
                fam[row["v3_family"]][f"dom_{v2['domain']}"] += 1
            if v2["family"] and v2["accession"] not in seeds:
                k = ("agree" if row["p_family"] == v2["family"] and row["p_call"] == "family"
                     else "disagree" if row["p_call"] == "family"
                     else f"profile_{row['p_call']}")
                calib[v2["family"]][k] += 1
                calib_dom[v2["domain"] or "unplaced"][k] += 1
            if row["v3_basis"] == "conflict":
                conflicts[(v2["family"] or f"[{v2['superfamily'] or '-'}]",
                           row["p_family"], v2["decisive_tier"] or "-")] += 1
            if v2["status"] == "superfamily_only":
                sfo[(v2["superfamily"], row["v3_family"] or row["v3_basis"])] += 1
            if v2["status"] == "unassigned":
                unas[(row["p_call"], row["v3_family"] or "-")] += 1
            if v2["taxon_id"] == "9606" and v2["reviewed"] == "reviewed":
                human[v2["accession"]] = row
            yield row

    out = s3_dir() / "census_v3.tsv.gz"
    n = write_tsv(out, V3_FIELDS, rows())
    write_tables(n, basis, status, pcall, fam, calib, conflicts, sfo, unas,
                 human, fam_status, calib_dom)
    write_tsv(OUT_DIR / "call_support.tsv",
              ["v3_family", "v3_basis", "s2_tier", "p_call", "records"],
              [{"v3_family": f, "v3_basis": b, "s2_tier": t, "p_call": c, "records": v}
               for (f, b, t, c), v in sorted(support.items())])
    write_tsv(OUT_DIR / "manifest.tsv", ["file", "bytes", "sha256", "records"], [
        {"file": str(p.relative_to(s3_dir().parents[1])), "bytes": p.stat().st_size,
         "sha256": sha256(p), "records": r}
        for p, r in ((out, n), (pc, len(calls)))])
    s0_lib.live_progress(LIVE, "S3a", [("profiles built", True),
                                       ("census sweep", True),
                                       ("assigning", True)])
    print(f"[v3] {n:,} records; basis {dict(basis)}")
    return 0


def write_tables(n, basis, status, pcall, fam, calib, conflicts, sfo, unas,
                 human, fam_status, calib_dom) -> None:
    O = OUT_DIR
    write_tsv(O / "v3_basis.tsv", ["basis", "records"],
              [{"basis": k, "records": v} for k, v in basis.most_common()])
    write_tsv(O / "v3_status.tsv", ["status", "records"],
              [{"status": k, "records": v} for k, v in status.most_common()])
    write_tsv(O / "v2_status_x_profile.tsv", ["v2_status", "p_call", "records"],
              [{"v2_status": a, "p_call": b, "records": v}
               for (a, b), v in sorted(pcall.items())])
    doms = ("Eukaryota", "Bacteria", "Archaea", "Viruses")
    write_tsv(O / "v3_families.tsv",
              ["family", "superfamily", "catalogue_status", "records", "both",
               "s2_only", "profile_only"] + [d.lower() for d in doms],
              sorted(({"family": f, "superfamily": registry.family(f).superfamily,
                       "catalogue_status": fam_status[f],
                       "records": sum(c[b] for b in ("both", "s2_only", "profile_only")),
                       "both": c["both"], "s2_only": c["s2_only"],
                       "profile_only": c["profile_only"],
                       **{d.lower(): c[f"dom_{d}"] for d in doms}}
                      for f, c in fam.items()), key=lambda r: -r["records"]))
    kinds = ["agree", "disagree", "profile_superfamily_only", "profile_ambiguous",
             "profile_module", "profile_low_score", "profile_no_hit"]
    write_tsv(O / "calibration.tsv", ["family", "v2_calls"] + kinds + ["agreement"],
              sorted(({"family": f, "v2_calls": sum(c.values()),
                       **{k: c[k] for k in kinds},
                       "agreement": round(c["agree"] / (c["agree"] + c["disagree"]), 4)
                       if c["agree"] + c["disagree"] else ""}
                      for f, c in calib.items()), key=lambda r: -r["v2_calls"]))
    write_tsv(O / "calibration_by_domain.tsv", ["domain", "v2_calls"] + kinds + ["agreement"],
              [{"domain": d, "v2_calls": sum(c.values()), **{k: c[k] for k in kinds},
                "agreement": round(c["agree"] / (c["agree"] + c["disagree"]), 4)
                if c["agree"] + c["disagree"] else ""}
               for d, c in sorted(calib_dom.items())])
    write_tsv(O / "conflicts.tsv", ["s2_call", "profile_family", "s2_tier", "records"],
              [{"s2_call": a, "profile_family": b, "s2_tier": t, "records": v}
               for (a, b, t), v in conflicts.most_common()])
    write_tsv(O / "superfamily_only_resolved.tsv", ["v2_superfamily", "v3_outcome", "records"],
              [{"v2_superfamily": a, "v3_outcome": b, "records": v}
               for (a, b), v in sorted(sfo.items(), key=lambda x: (x[0][0], -x[1]))])
    write_tsv(O / "unassigned_fate.tsv", ["p_call", "v3_family", "records"],
              [{"p_call": a, "v3_family": b, "records": v}
               for (a, b), v in unas.most_common()])
    hr = []
    for r in read_tsv(ROOT / "results" / "census_v2" / "human_recall.tsv"):
        row = human.get(r["accession"], {})
        v3f = row.get("v3_family", "")
        hr.append({"gene": r["gene"], "expected_family": r["expected_family"],
                   "accession": r["accession"], "s2_family": r["s2_family"],
                   "p_call": row.get("p_call", ""), "p_family": row.get("p_family", ""),
                   "v3_family": v3f, "v3_basis": row.get("v3_basis", "not_enumerated"),
                   "correct_s2": r["correct"],
                   "correct_v3": "yes" if v3f == r["expected_family"] else "no"})
    write_tsv(O / "human_recall_v3.tsv", list(hr[0]), hr)
    (O / "summary.json").write_text(json.dumps({
        "records": n, "basis": dict(basis), "status": dict(status),
        "human_correct_s2": sum(r["correct_s2"] == "yes" for r in hr),
        "human_correct_v3": sum(r["correct_v3"] == "yes" for r in hr),
        "human_genes": len(hr)}, indent=2) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
