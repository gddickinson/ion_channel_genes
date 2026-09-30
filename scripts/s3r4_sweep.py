"""s3r4_sweep.py — census v2 r4 through the profile library (D43).

    python3 scripts/s3r4_sweep.py [--rev r4] prep     # delta NR + map
    python3 scripts/s3r4_sweep.py [--rev r4] search [--jobs 4 --cpu 2]

`--rev r5` (the viroporin signatures) adds no profile: every profile is
searched over r5's delta alone, into `domtbl_r5/`.

r4 added 26,783 records to census v2 and 11 profiles to the library. Rather
than re-sweep all 102 profiles over 1.2 M sequences, and so move every
existing E-value with a new `-Z`, r4 searches only what is new:

* the **11 new profiles over r3's NR database** (into `hmmer/s3/domtbl/`,
  beside the 91 frozen domtblouts, which are not touched);
* **all 102 profiles over the delta NR** — r4 sequences whose sequence id
  is not already in r3's NR (into `hmmer/s3/domtbl_r4/`).

`-Z` stays at r3's NR size for every search, so each (sequence, profile)
pair has exactly one E-value, comparable to every other, and no r3 hit
changes. Every target sequence lives in exactly one of the two databases,
so `s3_census_v3.py` reads both directories without double counting.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s3_hmm_lib import (OUT_DIR, census_v2_fasta, iter_fasta,  # noqa: E402
                                s3_dir, seq_id, sha256, write_tsv)
from scripts.s3_sweep import _search, nr_fasta, nr_map                 # noqa: E402
from src.catalogue import registry                                      # noqa: E402

NEW_PROFILES = ("pacc", "tmco1", "tmem87", "tmem109", "clcc1", "mitok", "gphr",
                "assoc_kchip", "nonchannel_emc3", "nonchannel_gost",
                "nonchannel_bri3bp", "nonchannel_ncs")


REV = "r4"          # set from --rev


def delta_fasta(rev: str | None = None) -> Path:
    return s3_dir() / f"census_v2.{rev or REV}_delta.nr.fasta"


def delta_map(rev: str | None = None) -> Path:
    return s3_dir() / f"census_v2.{rev or REV}_delta.nr_map.tsv.gz"


def earlier_maps() -> list[Path]:
    """r3's NR map and every earlier revision's delta map."""
    out = [nr_map()]
    for rev in ("r4", "r5", "r6"):
        if rev >= REV:
            break
        out.append(delta_map(rev))
    return out


def r3_z() -> int:
    return json.loads((OUT_DIR / "sweep_db.json").read_text())["unique_sequences"]


def cmd_prep(_a) -> None:
    old_ids, old_acc = set(), set()
    for mp in earlier_maps():
        with gzip.open(mp, "rt") as fh:
            next(fh)
            for line in fh:
                acc, sid, _ = line.rstrip("\n").split("\t")
                old_acc.add(acc)
                old_ids.add(sid)
    new_ids: set[str] = set()
    n_rec = 0
    with delta_fasta().open("w") as fa, gzip.open(delta_map(), "wt") as mp:
        mp.write("accession\tseq_id\tlength\n")
        for head, seq in iter_fasta(census_v2_fasta()):
            acc = head.split("|", 1)[0]
            if acc in old_acc:
                continue
            sid = seq_id(seq)
            mp.write(f"{acc}\t{sid}\t{len(seq)}\n")
            n_rec += 1
            if sid not in old_ids and sid not in new_ids:
                new_ids.add(sid)
                fa.write(f">{sid}\n{seq}\n")
    stats = {"records": n_rec, "new_unique_sequences": len(new_ids),
             "delta_sha256": sha256(delta_fasta()), "z_held_at_r3": r3_z()}
    db = json.loads((OUT_DIR / "sweep_db.json").read_text())
    db[REV] = stats
    (OUT_DIR / "sweep_db.json").write_text(json.dumps(db, indent=2) + "\n")
    print(f"[prep] {n_rec:,} {REV} records → {len(new_ids):,} sequences not searched before")


def cmd_search(a) -> None:
    z = r3_z()
    fams = [f.key for f in registry.families()]
    missing = [p for p in NEW_PROFILES if p not in fams]
    if missing:
        raise SystemExit(f"not catalogue families: {missing}")
    # r4 added profiles: they go over r3's NR too (and over r4's delta below).
    # A later revision adds none, but any profile rebuilt since must also be
    # re-searched over earlier deltas — `_search` reuses a domtbl only when
    # the profile's SHA-256 matches, so re-running r4's search does that.
    jobs = ([(p, nr_fasta(), s3_dir("domtbl")) for p in NEW_PROFILES]
            if REV == "r4" else [])
    jobs += [(p, delta_fasta(), s3_dir(f"domtbl_{REV}")) for p in fams
             if (s3_dir("profiles") / f"{p}.hmm").exists()]
    rows = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(_search, p, db, z, a.cpu, out): (p, db.name)
                for p, db, out in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            p, dbn = futs[fut]
            rec = fut.result()
            rows.append({**rec, "db": dbn})
            print(f"[search] {i}/{len(jobs)} {p:26s} {dbn:32s} "
                  f"{rec['targets_reported']:7,d} targets {rec['seconds']:6.1f} s"
                  f"{' (reused)' if rec.get('reused') else ''}", flush=True)
    rows.sort(key=lambda r: (r["db"], r["profile"]))
    write_tsv(OUT_DIR / f"sweep_runs_{REV}.tsv",
              ["profile", "db", "targets_reported", "seconds", "hmm_sha256",
               "db_sequences", "command"], rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rev", default="r4")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prep")
    p = sub.add_parser("search")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--cpu", type=int, default=2)
    a = ap.parse_args()
    global REV
    REV = a.rev
    {"prep": cmd_prep, "search": cmd_search}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
