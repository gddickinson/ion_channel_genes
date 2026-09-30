"""s2r4_delta.py — census v2 revision r4: enumerate the signatures added after S20.

    python3 scripts/s2r4_delta.py [--rev r4] counts   # delta per shard; partition + union checks
    python3 scripts/s2r4_delta.py [--rev r4] walk [-j 4]
    python3 scripts/s2r4_delta.py [--rev r4] record   # append rows to shard/signature counts

`--rev` names the revision (r4: the families added after S20; r5: the
viroporin signatures, S4's emergent row). "Old" is every signature an
earlier revision recorded in `signature_counts.tsv`; "new" is what the
catalogue now enumerates beyond them.

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
REV = "r4"          # set from --rev


def old_signatures() -> list[str]:
    """Every signature an earlier revision enumerated — from the committed
    record (r3's rows carry no revision tag)."""
    return sorted(r["pfam"] for r in read_tsv(OUT_DIR / "signature_counts.tsv")
                  if (r.get("revision") or "r3") < REV)


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
        out.append((f"{REV}_{key}", q))
    return out


def state_file() -> Path:
    return raw_dir() / f"{REV}_counts.json"


def cmd_counts() -> int:
    new, old = new_signatures(), old_signatures()
    if not new:
        raise SystemExit(f"{REV}: no signature beyond the {len(old)} already enumerated")
    print(f"{REV}: {len(old)} old + {len(new)} new signatures: {', '.join(new)}")
    delta, rel = uniprot_count(f"({_or(new)}) NOT ({_or(old)})")
    if rel != R3_RELEASE:
        raise SystemExit(f"UniProt is on {rel}; census v2 is {R3_RELEASE} — a "
                         "delta across releases is not a revision (D35)")
    union_new, _ = uniprot_count(union_clause())
    union_old, _ = uniprot_count(_or(old))
    shards = {}
    for key, q in delta_shards():
        shards[key] = uniprot_count(q)[0]
        print(f"  {key:26s} {shards[key]:>7,}")
    st = {"release": rel, "delta": delta, "union_new": union_new, "union_old": union_old,
          "shards": shards, "shard_sum": sum(shards.values()),
          "new_signatures": new}
    state_file().write_text(json.dumps(st, indent=1))
    ok = (st["shard_sum"] == delta) and (union_old + delta == union_new)
    print(f"delta {delta:,}  shards sum {st['shard_sum']:,}  "
          f"old union {union_old:,} + delta = {union_old + delta:,} vs new union {union_new:,}  "
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
    if any(r["shard"].startswith(f"{REV}_") for r in shard_rows):
        raise SystemExit(f"{REV} rows already recorded")
    out = [r for r in shard_rows if r["shard"] != "UNION"]
    for key, q in delta_shards():
        out.append({"shard": key, "clause": f"{REV} delta: new signature, none earlier",
                    "uniprot_count": st["shards"][key], "fetched": "", "status": ""})
    out.append({"shard": "UNION",
                "clause": f"any pore signature ({len(signatures())}, {REV})",
                "uniprot_count": st["union_new"], "fetched": "", "status": ""})
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
                         "fetched": "", "status": "", "revision": REV})
        print(f"  {a} uniprot {u:,} interpro {ip}")
    write_tsv(OUT_DIR / "signature_counts.tsv", cols,
              [[r.get(c, "") for c in cols] for r in sig_rows])
    rel = json.loads((OUT_DIR / "release.json").read_text())
    rel[REV] = {"delta_records": st["delta"], "union": st["union_new"],
                 "new_signatures": st["new_signatures"]}
    (OUT_DIR / "release.json").write_text(json.dumps(rel, indent=2))
    print("recorded; now: s2_enumerate.py verify")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["counts", "walk", "record"])
    ap.add_argument("-j", "--jobs", type=int, default=4)
    ap.add_argument("--rev", default="r4")
    a = ap.parse_args()
    global REV
    REV = a.rev
    return {"counts": cmd_counts, "walk": lambda: cmd_walk(a.jobs),
            "record": cmd_record}[a.cmd]()


if __name__ == "__main__":
    sys.exit(main())
