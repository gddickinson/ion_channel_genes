"""S3b step 3 — census v3 over the declared denominator (the S4 panel).

    python3 scripts/s3b_census_v3.py

Every panel entry carrying any evidence gets one row. The merge is a stated
rule, fixed before the tables were read:

* the entry is in census v2 (it carries one of the 67 enumerated pore
  signatures)            → its census v3a call, unchanged (`v3a:<basis>`)
* it is not, and the S3b profile sweep makes a family call
                         → that family, basis `panel_profile` — **what
                           domain search missed**
* not in v2, profile `superfamily_only` → that superfamily, `panel_profile`
* not in v2, no profile call, but inside a *clean or accepted* jackhmmer
  round of some family's run → `candidate`, basis `jackhmmer_only`. Never a
  family call: iterated homology to a seed is not a positive test (D14,
  D33); it is the upper bound on what the profile library cannot call.

Two checks run on the way: the S3a and S3b profile verdicts on the very
same sequences (in census v2 *and* the panel) must agree — the same
profiles, only `-Z` differs — and the 320 human census genes are scored
against the human proteome by gene symbol (scoring only, H15).

Bulk: `<data root>/hmmer/s3b/census_v3.tsv.gz`. Committed tables →
`results/panel_sweep/` (rendered by `s3b_report.py`, D13).
"""

from __future__ import annotations

import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    LIVE, iter_fasta, read_tsv, s3_dir, seq_id, write_tsv,
)
from scripts.s3b_lib import OUT_DIR, s3b_dir  # noqa: E402
from scripts.s3b_sweep import calls_path, universe_path  # noqa: E402
from src.catalogue import registry  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

FIELDS = ["target", "accession", "species", "group", "length", "in_v2",
          "v3a_family", "v3a_superfamily", "v3a_basis",
          "p_call", "p_family", "p_superfamily", "p_confidence", "win_score",
          "win_coverage", "runner", "rel_margin", "jh_families",
          "v3_family", "v3_superfamily", "v3_status", "v3_basis"]
PCALL_KEYS = ("p_call", "p_family", "p_superfamily", "p_confidence",
              "win_score", "win_coverage", "runner", "rel_margin")
_GN = re.compile(r"\bGN=(\S+)")


def merge(in_v2: bool, v3a: dict, p: dict, jh: str) -> dict:
    if in_v2:
        return {"v3_family": v3a.get("v3_family", ""),
                "v3_superfamily": v3a.get("v3_superfamily", ""),
                "v3_basis": "v3a:" + v3a.get("v3_basis", "missing")}
    if p.get("p_call") == "family":
        return {"v3_family": p["p_family"], "v3_superfamily": p["p_superfamily"],
                "v3_basis": "panel_profile"}
    if p.get("p_call") == "superfamily_only":
        return {"v3_family": "", "v3_superfamily": p["p_superfamily"],
                "v3_basis": "panel_profile"}
    if jh:
        return {"v3_family": "", "v3_superfamily": "",
                "v3_basis": "jackhmmer_only"}
    return {}


def status_of(row: dict, fam_status: dict[str, str]) -> str:
    b = row["v3_basis"]
    if b.endswith("conflict"):
        return "conflict"
    if b == "jackhmmer_only":
        return "candidate"
    if row["v3_family"]:
        return fam_status[row["v3_family"]]
    return "superfamily_only" if row["v3_superfamily"] else "unassigned"


def load_v3a(accs: set[str]) -> dict[str, dict]:
    out = {}
    with gzip.open(s3_dir() / "census_v3.tsv.gz", "rt") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            a = line.split("\t", 1)[0]
            if a in accs:
                out[a] = dict(zip(hdr, line.rstrip("\n").split("\t")))
    return out


def load_s3a_calls(seqs: dict[str, str]) -> dict[str, dict]:
    """S3a profile verdicts for the given accession → seq_id."""
    want = set(seqs.values())
    return {r["target"]: r for r in read_tsv(s3_dir() / "profile_calls.tsv.gz")
            if r["target"] in want}


def human_recall(calls: dict, v3_by_target: dict) -> list[dict]:
    fa = require_data_root() / "proteomes" / "s4" / "UP000005640_9606.fasta.gz"
    by_gene: dict[str, list[str]] = defaultdict(list)
    for head, _ in iter_fasta(fa):
        m = _GN.search(head)
        if m:
            by_gene[m.group(1)].append(head.split(maxsplit=1)[0])
    rows = []
    for gene, fam in sorted(registry.human_genes().items()):
        tg = by_gene.get(gene, [])
        if not tg:
            rows.append({"gene": gene, "family": fam, "verdict": "not_in_proteome"})
            continue
        t = tg[0]
        v = v3_by_target.get(t, {})
        got = v.get("v3_family", "")
        verdict = ("right_family" if got == fam else
                   "wrong_family" if got else
                   v.get("v3_status") or "no_evidence")
        rows.append({"gene": gene, "family": fam, "target": t,
                     "v3_family": got, "v3_status": v.get("v3_status", ""),
                     "v3_basis": v.get("v3_basis", ""),
                     "p_call": calls.get(t, {}).get("p_call", "no_hit"),
                     "p_family": calls.get(t, {}).get("p_family", ""),
                     "verdict": verdict})
    return rows


def main() -> int:
    s0_lib.live_progress(LIVE, "S3b", [("panel sweep", True), ("assign", True),
                                       ("jackhmmer", True), ("census v3", False)])
    fam_status = {f.key: f.status.value for f in registry.families()}
    census = {f.key for f in registry.census_families()}
    uni = {r["target"]: r for r in read_tsv(universe_path())}
    calls = {r["target"]: r for r in read_tsv(calls_path())}
    jh: dict[str, set[str]] = defaultdict(set)
    for r in read_tsv(s3b_dir() / "jackhmmer_accepted.tsv.gz"):
        jh[r["target"]].add(r["family"])

    acc_to_t = {u["accession"]: t for t, u in uni.items()}
    v2_seq: dict[str, str] = {}
    with gzip.open(s3_dir() / "census_v2.nr_map.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            a, s, _ = line.rstrip("\n").split("\t")
            if a in acc_to_t:
                v2_seq[a] = s
    v3a = load_v3a(set(v2_seq))
    s3a = load_s3a_calls(v2_seq)
    print(f"[v3] {len(v2_seq):,} panel entries are census v2 records")

    rows, by_t = [], {}
    for t in sorted(set(calls) | set(jh) | {acc_to_t[a] for a in v2_seq}):
        u = uni[t]
        a = u["accession"]
        p = calls.get(t, {"p_call": "no_hit"})
        m = merge(a in v2_seq, v3a.get(a, {}), p, ",".join(sorted(jh.get(t, ()))))
        if not m:
            continue
        row = {"target": t, "accession": a, "species": u["species"],
               "group": u["group"], "length": u["length"],
               "in_v2": int(a in v2_seq),
               "v3a_family": v3a.get(a, {}).get("v3_family", ""),
               "v3a_superfamily": v3a.get(a, {}).get("v3_superfamily", ""),
               "v3a_basis": v3a.get(a, {}).get("v3_basis", ""),
               **{k: p.get(k, "") for k in PCALL_KEYS},
               "jh_families": ",".join(sorted(jh.get(t, ()))), **m}
        row["v3_status"] = status_of(row, fam_status)
        rows.append(row)
        by_t[t] = row
    write_tsv(s3b_dir() / "census_v3.tsv.gz", FIELDS, rows)
    print(f"[v3] {len(rows):,} panel entries with evidence")

    # ---- instrument check: S3a vs S3b profile verdict on the same sequence
    inst = Counter()
    for a, s in v2_seq.items():
        p3a = s3a.get(s, {"p_call": "no_hit", "p_family": ""})
        p3b = calls.get(acc_to_t[a], {"p_call": "no_hit", "p_family": ""})
        same = (p3a["p_call"], p3a.get("p_family", "")) == \
               (p3b["p_call"], p3b.get("p_family", ""))
        inst[(p3a["p_call"], p3b["p_call"], "same" if same else "differs")] += 1
    write_tsv(OUT_DIR / "instrument_check.tsv",
              ["s3a_p_call", "s3b_p_call", "family", "records"],
              [{"s3a_p_call": k[0], "s3b_p_call": k[1], "family": k[2],
                "records": v} for k, v in sorted(inst.items(), key=lambda x: -x[1])])

    # ---- status / basis
    sb = Counter((r["v3_basis"], r["v3_status"]) for r in rows)
    write_tsv(OUT_DIR / "census_status.tsv", ["v3_basis", "v3_status", "records"],
              [{"v3_basis": b, "v3_status": s, "records": n}
               for (b, s), n in sorted(sb.items(), key=lambda x: -x[1])])

    # ---- what domain search missed, per family and per group
    fam = defaultdict(Counter)
    grp = defaultdict(Counter)
    cell = defaultdict(Counter)
    for r in rows:
        f = r["v3_family"]
        if not f:
            continue
        k = "in_v2" if r["in_v2"] else "panel_profile"
        fam[f][k] += 1
        grp[(f, r["group"])][k] += 1
        if k == "panel_profile" and r["p_confidence"] == "high":
            fam[f]["high"] += 1
            grp[(f, r["group"])]["high"] += 1
        cell[(f, r["species"])][k] += 1
    sf_of = {f.key: f.superfamily for f in registry.families()}
    frows = []
    for f in sorted(fam, key=lambda f: (sf_of[f], f)):
        c = fam[f]
        tot = c["in_v2"] + c["panel_profile"]
        frows.append({"family": f, "superfamily": sf_of[f],
                      "status": fam_status[f], "called": tot,
                      "in_census_v2": c["in_v2"],
                      "missed_by_domain_search": c["panel_profile"],
                      "missed_high": c["high"],
                      "missed_frac": round(c["panel_profile"] / tot, 4)})
    write_tsv(OUT_DIR / "domain_search_missed.tsv",
              ["family", "superfamily", "status", "called", "in_census_v2",
               "missed_by_domain_search", "missed_high", "missed_frac"], frows)
    write_tsv(OUT_DIR / "missed_by_group.tsv",
              ["family", "group", "in_census_v2", "missed_by_domain_search",
               "missed_high"],
              [{"family": f, "group": g, "in_census_v2": c["in_v2"],
                "missed_by_domain_search": c["panel_profile"],
                "missed_high": c["high"]}
               for (f, g), c in sorted(grp.items())])
    write_tsv(OUT_DIR / "missed_records.tsv",
              ["target", "species", "group", "length", "v3_family", "p_confidence",
               "win_score", "win_coverage", "runner", "rel_margin"],
              (r for r in rows if r["v3_basis"] == "panel_profile" and r["v3_family"]))

    # ---- family × species presence (census families only), for S10
    species = sorted({u["species"] for u in uni.values()})
    write_tsv(OUT_DIR / "family_by_species.tsv",
              ["family", "superfamily", "species", "records",
               "missed_by_domain_search"],
              [{"family": f, "superfamily": sf_of[f], "species": s,
                "records": cell[(f, s)]["in_v2"] + cell[(f, s)]["panel_profile"],
                "missed_by_domain_search": cell[(f, s)]["panel_profile"]}
               for f in sorted(census) for s in species])

    # ---- jackhmmer-only candidates, by the family whose run found them
    cand = Counter()
    for r in rows:
        if r["v3_basis"] == "jackhmmer_only":
            for f in r["jh_families"].split(","):
                cand[(f, r["p_call"])] += 1
    write_tsv(OUT_DIR / "jackhmmer_candidates.tsv",
              ["family", "profile_verdict", "records"],
              [{"family": f, "profile_verdict": p, "records": n}
               for (f, p), n in sorted(cand.items())])

    hr = human_recall(calls, by_t)
    write_tsv(OUT_DIR / "human_recall.tsv",
              ["gene", "family", "target", "v3_family", "v3_status", "v3_basis",
               "p_call", "p_family", "verdict"], hr)
    summ = {"panel_entries": len(uni), "with_evidence": len(rows),
            "in_census_v2": len(v2_seq),
            "status": dict(Counter(r["v3_status"] for r in rows)),
            "basis": dict(Counter(r["v3_basis"] for r in rows)),
            "human_recall": dict(Counter(r["verdict"] for r in hr)),
            "instrument_same": sum(v for k, v in inst.items() if k[2] == "same"),
            "instrument_total": sum(inst.values())}
    (OUT_DIR / "summary.json").write_text(json.dumps(summ, indent=2) + "\n")
    print(json.dumps(summ, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
