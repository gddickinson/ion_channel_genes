"""S2 step 1 — enumerate the census space, completely, and prove it.

    python3 scripts/s2_enumerate.py counts      # per-signature + per-shard counts
    python3 scripts/s2_enumerate.py walk [-j 4] # cursor-walk every shard (resumable)
    python3 scripts/s2_enumerate.py verify      # fetched == count, everywhere

Three completeness checks, each a table in `results/census_v2/`:

1. **Partition** (`shard_counts.tsv`): the taxonomic shards' counts sum to
   the union count, so walking the shards is walking the union.
2. **Pagination** (`shard_counts.tsv`): records fetched per shard == UniProt's
   own `x-total-results` for that shard.
3. **Per signature** (`signature_counts.tsv`): for each of the 67 pore
   signatures, the number of fetched records carrying it == UniProt's count
   for that signature alone — and InterPro's count alongside, with the
   difference, because the two services are on different release cycles.

Raw pages are archived gzipped under `<data root>/raw_api/s2/pages/<shard>/`,
with a `state.json` per shard holding the resume cursor. A walk interrupted
at page 400 resumes at page 401.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s0_lib import read_tsv, write_tsv                  # noqa: E402
from scripts.s2_lib import (                                    # noqa: E402
    FIELDS, OUT_DIR, PAGE, SHARDS, UNIPROT, get, interpro_count, iter_pages,
    live, next_link, pfam_dict, raw_dir, shard_query, signatures,
    union_clause, uniprot_count,
)


# ------------------------------------------------------------- counts
def cmd_counts() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    union, release = uniprot_count(union_clause())
    rows, total = [], 0
    for key, clause in SHARDS:
        n, _ = uniprot_count(shard_query(clause))
        total += n
        rows.append([key, clause, n, "", ""])
        print(f"  shard {key:22s} {n:>9,}")
    rows.append(["UNION", "any pore signature", union, "", ""])
    write_tsv(OUT_DIR / "shard_counts.tsv",
              ["shard", "clause", "uniprot_count", "fetched", "status"], rows)
    print(f"union {union:,}  shard sum {total:,}  release {release}")
    sig_rows = []
    for a in signatures():
        u, _ = uniprot_count(f"(xref:pfam-{a})")
        ip = interpro_count(a)
        sig_rows.append([a, u, ip if ip is not None else "",
                         (u - ip) if ip is not None else "", "", ""])
        print(f"  {a}  uniprot {u:>8,}  interpro {ip}")
    write_tsv(OUT_DIR / "signature_counts.tsv",
              ["pfam", "uniprot_count", "interpro_count", "uniprot_minus_interpro",
               "fetched", "status"], sig_rows)
    (OUT_DIR / "release.json").write_text(json.dumps(
        {"uniprot_release": release, "counted": time.strftime("%Y-%m-%d")}, indent=2))
    if total != union:
        print(f"PARTITION FAILS: shards sum to {total:,}, union is {union:,}")
        return 1
    return 0


# ------------------------------------------------------------- walk
def _state_path(shard: str) -> Path:
    return raw_dir() / "pages" / shard / "state.json"


def walk_shard(key: str, clause: str, query: str | None = None) -> dict:
    """Walk one shard to exhaustion, resuming from its saved cursor.

    `query` replaces the union-scoped shard query outright (r4's delta
    shards, `s2r4_delta.py`)."""
    d = raw_dir() / "pages" / key
    d.mkdir(parents=True, exist_ok=True)
    sp = _state_path(key)
    st = json.loads(sp.read_text()) if sp.exists() else {
        "next": None, "pages": 0, "fetched": 0, "done": False, "count": None}
    if st["done"]:
        return st
    url, params = (st["next"], None) if st["next"] else (
        UNIPROT, {"query": query or shard_query(clause), "fields": FIELDS,
                  "format": "json", "size": PAGE})
    while True:
        r = get(url, params)
        if st["count"] is None:
            st["count"] = int(r.headers["x-total-results"])
        data = r.json()
        st["pages"] += 1
        with gzip.open(d / f"page_{st['pages']:05d}.json.gz", "wt") as fh:
            json.dump(data, fh)
        st["fetched"] += len(data.get("results", []))
        st["next"] = next_link(r)
        st["done"] = st["next"] is None
        sp.write_text(json.dumps(st))
        if st["done"]:
            return st
        url, params = st["next"], None


def cmd_walk(jobs: int) -> int:
    states: dict[str, dict] = {}

    def progress() -> None:
        steps = []
        for k, _ in SHARDS:
            sp = _state_path(k)
            s = json.loads(sp.read_text()) if sp.exists() else {}
            lab = f"{k}: {s.get('fetched', 0):,}/{s.get('count') or '?'}"
            steps.append((lab, bool(s.get("done"))))
        live(steps, workers=jobs)

    progress()
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(walk_shard, k, c): k for k, c in SHARDS}
        pending = set(futs)
        while pending:
            done = {f for f in pending if f.done()}
            for f in done:
                k = futs[f]
                states[k] = f.result()          # re-raises a strict failure
                print(f"  {k}: {states[k]['fetched']:,}/{states[k]['count']:,} "
                      f"in {states[k]['pages']} pages", flush=True)
            pending -= done
            progress()
            time.sleep(15)
    progress()
    return 0


# ------------------------------------------------------------- verify
def cmd_verify() -> int:
    shard_rows = read_tsv(OUT_DIR / "shard_counts.tsv")
    sig_rows = read_tsv(OUT_DIR / "signature_counts.tsv")
    carried: Counter = Counter()
    seen: set[str] = set()
    dupes = 0
    bad = 0
    for row in shard_rows:
        if row["shard"] == "UNION":
            continue
        n = 0
        for rec in iter_pages(row["shard"]):
            n += 1
            if rec["accession"] in seen:
                dupes += 1
                continue
            seen.add(rec["accession"])
            for a in pfam_dict(rec["pfam"]):
                carried[a] += 1
        row["fetched"] = n
        row["status"] = "ok" if n == int(row["uniprot_count"]) else "SHORT"
        bad += row["status"] != "ok"
        print(f"  {row['shard']:22s} {n:>9,} / {int(row['uniprot_count']):>9,}  {row['status']}")
    for row in shard_rows:
        if row["shard"] == "UNION":
            row["fetched"] = len(seen)
            row["status"] = "ok" if len(seen) == int(row["uniprot_count"]) else "SHORT"
            bad += row["status"] != "ok"
    for row in sig_rows:
        row["fetched"] = carried[row["pfam"]]
        row["status"] = "ok" if carried[row["pfam"]] == int(row["uniprot_count"]) else "MISMATCH"
        bad += row["status"] != "ok"
    cols = ["shard", "clause", "uniprot_count", "fetched", "status"]
    write_tsv(OUT_DIR / "shard_counts.tsv", cols, [[r[c] for c in cols] for r in shard_rows])
    cols = ["pfam", "uniprot_count", "interpro_count", "uniprot_minus_interpro",
            "fetched", "status"]
    write_tsv(OUT_DIR / "signature_counts.tsv", cols, [[r[c] for c in cols] for r in sig_rows])
    print(f"unique {len(seen):,}; cross-shard duplicates {dupes}; problems {bad}")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["counts", "walk", "verify"])
    ap.add_argument("-j", "--jobs", type=int, default=4)
    a = ap.parse_args()
    return {"counts": cmd_counts, "verify": cmd_verify,
            "walk": lambda: cmd_walk(a.jobs)}[a.cmd]()


if __name__ == "__main__":
    sys.exit(main())
