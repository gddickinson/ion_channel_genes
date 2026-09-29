"""S3b step 1 — the frozen S3a profile library over the S4 panel proteomes.

    python3 scripts/s3b_sweep.py db                  # index the sweep DB
    python3 scripts/s3b_sweep.py search [--jobs 4 --cpu 2]
    python3 scripts/s3b_sweep.py assign

`db` streams S4's `panel_refprot.fasta` once and writes one row per entry
(accession, taxon, species, length, SHA-256 of the file) — the universe
every later completeness statement is made over.

`search` is S3a's `hmmsearch` step (`s3_sweep._search`) pointed at the panel
DB: every profile, `-Z` fixed to the panel size, resumable on each profile's
SHA-256. The profiles are S3a's, unchanged (91, incl. the S3a2 ABC decoy).

`assign` applies S3a's best-profile rule (D32) — the same three gates, the
same constants, imported rather than restated — to every panel entry any
profile reported.

Bulk → `<data root>/hmmer/s3b/`; committed → `results/panel_sweep/`.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_assign import ASSIGN_FIELDS, assign_all, collect_hits  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    LIVE, OUT_DIR as S3A_DIR, iter_fasta, read_tsv, sha256, write_tsv,
)
from scripts.s3_sweep import RUN_FIELDS, _search  # noqa: E402
from scripts.s3b_lib import (  # noqa: E402
    OUT_DIR, accession, header_taxid, panel_db, s3b_dir, species_by_taxid,
)
from src.catalogue import registry  # noqa: E402


def universe_path() -> Path:
    return s3b_dir() / "panel_universe.tsv.gz"


def calls_path() -> Path:
    return s3b_dir() / "panel_profile_calls.tsv.gz"


def build_db() -> dict:
    sp = species_by_taxid()
    per_sp: dict[str, int] = {}
    n = unmapped = 0
    with gzip.open(universe_path(), "wt") as out:
        out.write("target\taccession\ttaxon_id\tspecies\tgroup\tlength\n")
        for head, seq in iter_fasta(panel_db()):
            tgt = head.split(maxsplit=1)[0]
            tax = header_taxid(head)
            row = sp.get(tax)
            if row is None:
                unmapped += 1
            name = row["species"] if row else ""
            per_sp[name] = per_sp.get(name, 0) + 1
            out.write(f"{tgt}\t{accession(tgt)}\t{tax}\t{name}\t"
                      f"{row['group'] if row else ''}\t{len(seq)}\n")
            n += 1
    stats = {"db": str(panel_db()), "db_sha256": sha256(panel_db()),
             "sequences": n, "unmapped_taxon": unmapped,
             "species": len([k for k in per_sp if k])}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "sweep_db.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(f"[db] {n:,} entries, {stats['species']} species, "
          f"{unmapped} with an unmapped OX=")
    if unmapped:
        raise SystemExit("every panel entry must map to a manifest proteome")
    return stats


def search(jobs: int, cpu: int) -> list[dict]:
    n_db = json.loads((OUT_DIR / "sweep_db.json").read_text())["sequences"]
    build = {r["family"]: int(r["match_states"])
             for r in read_tsv(S3A_DIR / "profile_build.tsv")}
    profs = sorted((f.key for f in registry.families() if f.key in build),
                   key=lambda p: -build[p])
    if len(profs) != len(registry.families()):
        raise SystemExit("every catalogue family needs its S3a profile")
    out_dir, db = s3b_dir("domtbl"), panel_db()
    rows, t0 = [], time.time()
    with ThreadPoolExecutor(jobs) as ex:
        futs = {ex.submit(_search, p, db, n_db, cpu, out_dir): p for p in profs}
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            rows.append(rec)
            print(f"[search] {i}/{len(profs)} {rec['profile']:26s} "
                  f"{rec['targets_reported']:8,d} targets "
                  f"{rec['seconds']:7.1f} s{' (reused)' if rec.get('reused') else ''}"
                  f"  [{time.time() - t0:.0f} s]", flush=True)
            s0_lib.live_progress(LIVE, "S3b", [
                (f"panel sweep {i}/{len(profs)} profiles", i == len(profs)),
                ("assign", False), ("jackhmmer", False), ("census v3", False)],
                workers=jobs)
    rows.sort(key=lambda r: r["profile"])
    write_tsv(OUT_DIR / "sweep_runs.tsv", RUN_FIELDS, rows)
    return rows


def assign() -> int:
    paths = sorted(s3b_dir("domtbl").glob("*.domtbl.gz"))
    if len(paths) != len(registry.families()):
        raise SystemExit(f"{len(paths)} domtbls for "
                         f"{len(registry.families())} profiles — run search")
    calls = assign_all(collect_hits(paths))
    n = write_tsv(calls_path(), ASSIGN_FIELDS,
                  (calls[t] for t in sorted(calls)))
    print(f"[assign] {n:,} panel entries with ≥ 1 profile hit")
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["db", "search", "assign"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpu", type=int, default=2)
    a = ap.parse_args()
    {"db": build_db, "assign": assign,
     "search": lambda: search(a.jobs, a.cpu)}[a.step]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
