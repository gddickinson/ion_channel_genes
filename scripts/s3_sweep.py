"""S3a step 2 — every family profile against every census v2 sequence.

    python3 scripts/s3_sweep.py prep                 # non-redundant sweep DB
    python3 scripts/s3_sweep.py search [--jobs 5 --cpu 2] [--only nav,cav]

`prep` collapses census v2 to its **unique sequences** (UniProt carries
many identical sequences under separate accessions; scoring them once is
the same result for less work) and writes the accession → sequence-id map
beside it, so every accession recovers its hits.

`search` runs `hmmsearch` per profile, resumably: a profile's domtblout is
reused only if the profile's SHA-256 matches the one recorded with it, so a
rebuilt profile is re-searched and an unchanged one is not. `-Z` is fixed
to the database size so E-values are comparable across profiles and runs.
A failed hmmsearch is a hard failure (D28).

Everything lands under `<data root>/hmmer/s3/`; the committed record is
`results/census_v3/sweep_runs.tsv`.
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    LIVE, OUT_DIR, census_v2_fasta, iter_fasta, read_tsv, s3_dir, seq_id,
    sha256, write_tsv,
)
from src.catalogue import registry  # noqa: E402

# Reporting thresholds. A target is reported for a profile at full-sequence
# E ≤ 1e-3; the assignment's own gates (score, coverage, margin) are applied
# afterwards in s3_assign.py, so these only bound the output size.
REPORT_E = "1e-3"
REPORT_DOME = "1e-3"
RUN_FIELDS = ["profile", "targets_reported", "seconds", "hmm_sha256",
              "db_sequences", "command"]


def nr_fasta() -> Path:
    return s3_dir() / "census_v2.nr.fasta"


def nr_map() -> Path:
    return s3_dir() / "census_v2.nr_map.tsv.gz"


def prep() -> dict:
    seen: set[str] = set()
    n_rec = n_res = 0
    with nr_fasta().open("w") as fa, gzip.open(nr_map(), "wt") as mp:
        mp.write("accession\tseq_id\tlength\n")
        for head, seq in iter_fasta(census_v2_fasta()):
            acc = head.split("|", 1)[0]
            sid = seq_id(seq)
            mp.write(f"{acc}\t{sid}\t{len(seq)}\n")
            n_rec += 1
            if sid not in seen:
                seen.add(sid)
                n_res += len(seq)
                fa.write(f">{sid}\n{seq}\n")
    stats = {"records": n_rec, "unique_sequences": len(seen),
             "unique_residues": n_res, "nr_sha256": sha256(nr_fasta())}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "sweep_db.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(f"[prep] {n_rec:,} records → {len(seen):,} unique sequences "
          f"({n_res:,} residues)")
    return stats


def _search(prof: str, db: Path, n_db: int, cpu: int,
            out_dir: Path | None = None) -> dict:
    hmm = s3_dir("profiles") / f"{prof}.hmm"
    out = (out_dir or s3_dir("domtbl")) / f"{prof}.domtbl"
    side = out.with_suffix(".json")
    h = sha256(hmm)
    if side.exists() and json.loads(side.read_text()).get("hmm_sha256") == h \
            and Path(str(out) + ".gz").exists():
        return {**json.loads(side.read_text()), "reused": True}
    cmd = ["hmmsearch", "--cpu", str(cpu), "--noali", "-E", REPORT_E,
           "--domE", REPORT_DOME, "-Z", str(n_db), "--domZ", str(n_db),
           "-o", "/dev/null", "--domtblout", str(out), str(hmm), str(db)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"{prof}: hmmsearch exit {r.returncode}: {r.stderr[:300]}")
    targets = {l.split(maxsplit=1)[0] for l in out.open() if not l.startswith("#")}
    subprocess.run(["gzip", "-f", str(out)], check=True)
    rec = {"profile": prof, "targets_reported": len(targets),
           "seconds": round(time.time() - t0, 1), "hmm_sha256": h,
           "db_sequences": n_db, "command": " ".join(cmd[:-2] + ["<hmm>", "<db>"])}
    side.write_text(json.dumps(rec, indent=2))
    return rec


def search(jobs: int, cpu: int, only: set[str], db: Path | None = None,
           n_db: int | None = None) -> list[dict]:
    db = db or nr_fasta()
    if n_db is None:
        n_db = json.loads((OUT_DIR / "sweep_db.json").read_text())["unique_sequences"]
    build = {r["family"]: int(r["match_states"])
             for r in read_tsv(OUT_DIR / "profile_build.tsv")}
    profs = [f.key for f in registry.families()
             if f.key in build and (not only or f.key in only)]
    profs.sort(key=lambda p: -build[p])          # longest first
    rows, t0 = [], time.time()
    with ThreadPoolExecutor(jobs) as ex:
        futs = {ex.submit(_search, p, db, n_db, cpu): p for p in profs}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            rows.append(rec)
            print(f"[search] {i}/{len(profs)} {rec['profile']:26s} "
                  f"{rec['targets_reported']:8,d} targets "
                  f"{rec['seconds']:7.1f} s{' (reused)' if rec.get('reused') else ''}"
                  f"  [{time.time() - t0:.0f} s]", flush=True)
            s0_lib.live_progress(LIVE, "S3a", [
                ("profiles built", True),
                (f"census sweep {i}/{len(profs)} profiles", i == len(profs))],
                workers=jobs)
    if only:        # merge into the record of the full sweep
        rows += [r for r in read_tsv(OUT_DIR / "sweep_runs.tsv")
                 if r["profile"] not in only]
    rows.sort(key=lambda r: r["profile"])
    write_tsv(OUT_DIR / "sweep_runs.tsv", RUN_FIELDS, rows)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prep", "search"])
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--cpu", type=int, default=2)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    if a.step == "prep":
        prep()
    else:
        search(a.jobs, a.cpu, set(filter(None, a.only.split(","))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
