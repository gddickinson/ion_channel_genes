"""S3a step 3 — measure the profile instrument before trusting it.

    python3 scripts/s3_benchmark.py

Two test sets, scored with the same `s3_assign.assign_one()` the census uses:

**A. The S1 control panel** (97 proteins: 72 channel positives, 25 decoys
from the control families). Most are seeds (R1/R2), so every panel protein
that seeded a profile is scored **leave-one-out** (D29): its family profile
is rebuilt from the same alignment with its row removed and all-gap columns
dropped, and the query is scored against that profile instead. A family
whose only seed is the query has no LOO profile; that is recorded as
`sole_seed` and the query is scored with the family's profile absent.

**B. Held-out orthologues.** Reviewed, non-human census v2 records that are
not seeds and whose UniProt gene symbol equals a catalogue human gene
(case-insensitive) → that gene's family. The symbol is used here only to
*score* a call made without it (H15), exactly as S2's `human_recall.tsv`
does. This is the non-human test the S1 panel lacks.

A decoy (control-family) query passes if it is not called to a census
(channel) family. Outputs: `benchmark_calls.tsv`, `benchmark_recall.tsv`,
`benchmark_summary.json` in `results/census_v3/`.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.s3_assign import ASSIGN_FIELDS, assign_one, collect_hits, fam_superfamily  # noqa: E402
from scripts.s3_build_profiles import gather_sequences  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    OUT_DIR, S1_CALLS, census_v2_fasta, iter_census_v2, iter_fasta, read_fasta,
    read_tsv, s3_dir, write_fasta, write_tsv,
)
from src.catalogue import registry  # noqa: E402
from src.catalogue.schema import CENSUS_STATUSES  # noqa: E402

CALL_FIELDS = ["set", "accession", "expected_family", "expected_superfamily",
               "loo", "outcome"] + ASSIGN_FIELDS[1:]


def _z() -> str:
    return str(json.loads((OUT_DIR / "sweep_db.json").read_text())["unique_sequences"])


def hmmsearch(hmm: Path, fasta: Path, out: Path) -> Path:
    cmd = ["hmmsearch", "--cpu", "8", "--noali", "-E", "1e-3", "--domE",
           "1e-3", "-Z", _z(), "--domZ", _z(), "-o", "/dev/null",
           "--domtblout", str(out), str(hmm), str(fasta)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"hmmsearch failed on {hmm.name}: {r.stderr[:300]}")
    return out


def loo_profile(fam: str, acc: str, work: Path) -> Path | None:
    """The family profile rebuilt without `acc` (same alignment, row dropped)."""
    aln = read_fasta(s3_dir("aln") / f"{fam}.afa")
    rest = {k: v for k, v in aln.items() if k != acc}
    if not rest:
        return None
    width = len(next(iter(rest.values())))
    keep = [i for i in range(width) if any(s[i] != "-" for s in rest.values())]
    afa = work / f"{fam}__{acc}.afa"
    write_fasta(afa, [(k, "".join(s[i] for i in keep)) for k, s in rest.items()])
    hmm = work / f"{fam}__{acc}.hmm"
    r = subprocess.run(["hmmbuild", "--cpu", "1", "-n", fam, "--informat",
                        "afa", str(hmm), str(afa)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"LOO hmmbuild {fam}/{acc}: {r.stderr[:300]}")
    return hmm


def outcome(row: dict, exp: str, exp_sf: str, census: set[str]) -> str:
    if exp not in census:                           # a decoy
        return ("decoy_called_channel" if row["p_family"] in census
                else "decoy_ok")
    if row["p_call"] == "family":
        return "correct" if row["p_family"] == exp else (
            "wrong_family_same_sf" if row["p_superfamily"] == exp_sf
            else "wrong_superfamily")
    if row["p_call"] == "superfamily_only":
        return "sf_only_correct" if row["p_superfamily"] == exp_sf else "wrong_superfamily"
    return "no_call"


def test_sets() -> tuple[list[dict], dict[str, str]]:
    genes = {g.upper(): f for g, f in registry.human_genes(census_only=False).items()}
    seeds = {r["accession"]: r["family"] for r in read_tsv(OUT_DIR / "seed_manifest.tsv")}
    rows = [{"set": "A_s1_panel", "accession": r["accession"],
             "expected_family": r["expected_family"]} for r in read_tsv(S1_CALLS)]
    cols = ("accession", "reviewed", "taxon_id", "fragment", "gene")
    for r in iter_census_v2(cols):
        if (r["reviewed"] == "reviewed" and r["taxon_id"] != "9606"
                and "Fragment" not in (r["fragment"] or "")   # "Precursor" is complete
                and r["accession"] not in seeds
                and r["gene"].upper() in genes):
            rows.append({"set": "B_orthologues", "accession": r["accession"],
                         "expected_family": genes[r["gene"].upper()]})
    return rows, seeds


def main() -> int:
    work = s3_dir("bench")
    rows, seeds = test_sets()
    accs = {r["accession"] for r in rows}
    seqs = {h.split("|", 1)[0]: s for h, s in iter_fasta(census_v2_fasta())
            if h.split("|", 1)[0] in accs}
    seqs.update(gather_sequences(accs - set(seqs)))
    fasta = work / "bench.fasta"
    write_fasta(fasta, sorted((a, seqs[a]) for a in accs))
    print(f"[bench] {sum(r['set'] == 'A_s1_panel' for r in rows)} panel + "
          f"{sum(r['set'] == 'B_orthologues' for r in rows)} orthologues")

    dt = hmmsearch(s3_dir("profiles") / "all.hmm", fasta, work / "all.domtbl")
    top = collect_hits([dt], keep=12)

    sf_of = fam_superfamily()
    census = {f.key for f in registry.families() if f.status in CENSUS_STATUSES}
    out = []
    for r in rows:
        acc, exp = r["accession"], r["expected_family"]
        hits = list(top.get(acc, []))
        loo = ""
        if acc in seeds:
            fam = seeds[acc]
            hits = [h for h in hits if h[1] != fam]
            one = work / f"q_{acc}.fa"
            write_fasta(one, [(acc, seqs[acc])])
            hmm = loo_profile(fam, acc, work)
            if hmm is None:
                loo = "sole_seed"
            else:
                loo = "loo"
                q = collect_hits([hmmsearch(hmm, one, work / f"q_{acc}.domtbl")])
                hits.extend(q.get(acc, []))
                hits.sort(reverse=True)
        call = assign_one(acc, hits, sf_of)
        exp_sf = sf_of.get(exp, "")
        out.append({**r, "expected_superfamily": exp_sf, "loo": loo,
                    "outcome": outcome(call, exp, exp_sf, census),
                    **{k: v for k, v in call.items() if k != "target"}})

    write_tsv(OUT_DIR / "benchmark_calls.tsv", CALL_FIELDS, out)
    per = defaultdict(Counter)
    for r in out:
        per[(r["set"], r["expected_family"])][r["outcome"]] += 1
    kinds = ["correct", "sf_only_correct", "wrong_family_same_sf",
             "wrong_superfamily", "no_call", "decoy_ok", "decoy_called_channel"]
    rec = [{"set": s, "family": f, "n": sum(c.values()), **{k: c[k] for k in kinds}}
           for (s, f), c in sorted(per.items())]
    write_tsv(OUT_DIR / "benchmark_recall.tsv", ["set", "family", "n"] + kinds, rec)
    summary = {}
    for s in ("A_s1_panel", "B_orthologues"):
        c = Counter(r["outcome"] for r in out if r["set"] == s)
        pos = sum(c[k] for k in kinds[:5])
        summary[s] = {"positives": pos, **{k: c[k] for k in kinds},
                      "loo": sum(r["loo"] == "loo" for r in out if r["set"] == s),
                      "sole_seed": sum(r["loo"] == "sole_seed" for r in out
                                       if r["set"] == s)}
    (OUT_DIR / "benchmark_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
