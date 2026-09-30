"""s2r4_delta.py — census v2 revision r4: enumerate the signatures added after S20.

    python3 scripts/s2r4_delta.py counts     # delta per shard; partition + union checks
    python3 scripts/s2r4_delta.py walk [-j 4]
    python3 scripts/s2r4_delta.py record     # append r4 rows to shard/signature counts

**What r4 adds (D43).** The catalogue additions of 2026-09-29 (PACC1, TMCO1,
TMEM87A, TMEM109, CLCC1, GPHRA/B; `src/catalogue/proposed.py`) bring eight
Pfam signatures into the census search space. Census v2 r3 enumerated 67.
Re-walking all 1.25 M records would re-fetch what is already archived, so r4
walks the **delta** — records carrying a new signature and none of the old
67 — in the same taxonomic shards, as shards `r4_<shard>`:

    (new₁ OR … OR new₈) NOT (old₁ OR … OR old₆₇) AND <shard clause>

Checked, not assumed: (1) the delta shards partition the delta; (2) r3's
union + the delta == UniProt's count for the 75-signature union; (3) after
the walk, `s2_enumerate.py verify` re-counts every signature over all
shards, old and new, against UniProt's own per-signature count. Records
already in r3 that carry a new signature keep their row and are
re-classified (`s2_classify.py --recheck … --archive r3`).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s0_lib import read_tsv, write_tsv                        # noqa: E402
from scripts.s2_enumerate import walk_shard                           # noqa: E402
from scripts.s2_lib import (OUT_DIR, SHARDS, interpro_count, raw_dir,  # noqa: E402
                            signatures, uniprot_count, union_clause)

R3_RELEASE = "2026_03"


def old_signatures() -> list[str]:
    """The 67 census v2 r3 enumerated — from its own committed record."""
    return sorted(r["pfam"] for r in read_tsv(OUT_DIR / "signature_counts.tsv")
                  if not r.get("revision"))


def new_signatures() -> list[str]:
    return sorted(set(signatures()) - set(old_signatures()))


def _or(accs: list[str]) -> str:
    return " OR ".join(f"(xref:pfam-{a})" for a in accs)


def delta_shards() -> list[tuple[str, str]]:
    """(key, full query) per shard: new signature, no old one, in the shard."""
    base = f"({_or(new_signatures())}) NOT ({_or(old_signatures())})"
    out = []
    for key, clause in SHARDS:
        q = f"{base} AND ({clause})" if not clause.startswith("NOT") else f"{base} {clause}"
        out.append((f"r4_{key}", q))
    return out


def state_file() -> Path:
    return raw_dir() / "r4_counts.json"


def cmd_counts() -> int:
    new, old = new_signatures(), old_signatures()
    if len(old) != 67 or len(new) != 8:
        raise SystemExit(f"expected 67 old + 8 new signatures, got {len(old)} + {len(new)}")
    delta, rel = uniprot_count(f"({_or(new)}) NOT ({_or(old)})")
    if rel != R3_RELEASE:
        raise SystemExit(f"UniProt is on {rel}; census v2 is {R3_RELEASE} — a "
                         "delta across releases is not a revision (D35)")
    union75, _ = uniprot_count(union_clause())
    union67, _ = uniprot_count(_or(old))
    shards = {}
    for key, q in delta_shards():
        shards[key] = uniprot_count(q)[0]
        print(f"  {key:26s} {shards[key]:>7,}")
    st = {"release": rel, "delta": delta, "union_75": union75, "union_67": union67,
          "shards": shards, "shard_sum": sum(shards.values()),
          "new_signatures": new}
    state_file().write_text(json.dumps(st, indent=1))
    ok = (st["shard_sum"] == delta) and (union67 + delta == union75)
    print(f"delta {delta:,}  shards sum {st['shard_sum']:,}  "
          f"union67 {union67:,} + delta = {union67 + delta:,} vs union75 {union75:,}  "
          f"→ {'ok' if ok else 'FAIL'}")
    return 0 if ok else 1


def cmd_walk(jobs: int) -> int:
    from concurrent.futures import ThreadPoolExecutor
    st = json.loads(state_file().read_text())
    shards = delta_shards()
    with ThreadPoolExecutor(jobs) as ex:
        results = list(ex.map(lambda kq: walk_shard(kq[0], "", query=kq[1]), shards))
    bad = 0
    for (key, _), res in zip(shards, results):
        ok = res["fetched"] == res["count"] == st["shards"][key]
        bad += not ok
        print(f"  {key:26s} {res['fetched']:>7,} / {res['count']:>7,} "
              f"(counted {st['shards'][key]:,}) {'ok' if ok else 'MISMATCH'}")
    return 1 if bad else 0


def cmd_record() -> int:
    st = json.loads(state_file().read_text())
    shard_rows = read_tsv(OUT_DIR / "shard_counts.tsv")
    if any(r["shard"].startswith("r4_") for r in shard_rows):
        raise SystemExit("r4 rows already recorded")
    out = [r for r in shard_rows if r["shard"] != "UNION"]
    for key, q in delta_shards():
        out.append({"shard": key, "clause": "r4 delta: new signature, no r3 signature",
                    "uniprot_count": st["shards"][key], "fetched": "", "status": ""})
    out.append({"shard": "UNION", "clause": "any pore signature (75, r4)",
                "uniprot_count": st["union_75"], "fetched": "", "status": ""})
    write_tsv(OUT_DIR / "shard_counts.tsv", list(out[0]), [list(r.values()) for r in out])
    sig_rows = read_tsv(OUT_DIR / "signature_counts.tsv")
    cols = list(sig_rows[0]) + (["revision"] if "revision" not in sig_rows[0] else [])
    for r in sig_rows:
        r.setdefault("revision", "")
    for a in st["new_signatures"]:
        u, _ = uniprot_count(f"(xref:pfam-{a})")
        ip = interpro_count(a)
        sig_rows.append({"pfam": a, "uniprot_count": u,
                         "interpro_count": ip if ip is not None else "",
                         "uniprot_minus_interpro": (u - ip) if ip is not None else "",
                         "fetched": "", "status": "", "revision": "r4"})
        print(f"  {a} uniprot {u:,} interpro {ip}")
    write_tsv(OUT_DIR / "signature_counts.tsv", cols,
              [[r.get(c, "") for c in cols] for r in sig_rows])
    rel = json.loads((OUT_DIR / "release.json").read_text())
    rel["r4"] = {"delta_records": st["delta"], "union_75": st["union_75"],
                 "new_signatures": st["new_signatures"]}
    (OUT_DIR / "release.json").write_text(json.dumps(rel, indent=2))
    print("recorded; now: s2_enumerate.py verify")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["counts", "walk", "record"])
    ap.add_argument("-j", "--jobs", type=int, default=4)
    a = ap.parse_args()
    return {"counts": cmd_counts, "walk": lambda: cmd_walk(a.jobs),
            "record": cmd_record}[a.cmd]()


if __name__ == "__main__":
    sys.exit(main())
