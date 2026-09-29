"""S2 step 2 — a positive call on every census record.

    python3 scripts/s2_classify.py [-j 8] [--shard KEY] [--limit N]

Every record the enumeration fetched is classified by the same `classify()`
the S1 benchmark measured, fed from the UniProt record itself: the complete
Pfam architecture with copy numbers, the transmembrane-feature count and the
sequence. Nothing is looked up per protein.

**Which tiers run, and why (D31).**

* *Architecture + hazard rules* — every record. Dictionary lookups.
* *K+ filter regex* — every record. A regex.
* *Four-repeat filter projection* (MAFFT against Nav1.5, ~2 s) — only on
  records carrying `PF00520` × ≥ 4, which is the trigger of hazard H1 and
  the only architecture on which DEKA / EEEE / EEKE can change the call.
  Everything else is scanned with the regex alone, and the gate is a column
  in the output so the reader can count what was not projected.
* *Reference identity* — **not run**. Measured in S1 at ~30 s per protein,
  it would take ~10,000 CPU-hours here. Per D28 that is a missing result
  with a note on every row (`reference_tier = not_run`), never a weaker
  substitute. The replacement is S3's profile assignment.

Output (data root, gzipped, one file per shard so a run resumes by shard):
`<data root>/raw_api/s2/calls/<shard>.tsv.gz`. `s2_census.py` assembles
them into the census and the committed tables.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s2_lib import ROOT, SHARDS, iter_pages, live, pfam_dict, raw_dir  # noqa: E402
from src.classify.classifier import ChannelQuery, classify                     # noqa: E402
from src.classify.motifs import FOUR_REPEAT_ANCHOR, verify_anchor              # noqa: E402
from src.utils.mafft import parse_fasta                                        # noqa: E402

FOUR_REPEAT_MIN = 4
NAV_LABEL = FOUR_REPEAT_ANCHOR.reference_label

COLS = ["accession", "reviewed", "taxon_id", "organism", "domain", "group",
        "phylum", "length", "fragment", "tm_count", "pfam", "seed_signatures",
        "family", "superfamily", "status", "confidence", "decisive_tier",
        "filter_string", "k_filter", "filter_projected", "reference_tier",
        "hazards", "conflicts", "notes", "evidence"]

_NAV = ""
_K_HITS = re.compile(r"\b[1-9]\d* K\+ filter match")
_SEEDS: frozenset = frozenset()


def nav_reference() -> str:
    panel = dict(parse_fasta((ROOT / "results" / "s0_baseline" /
                              "reference_panel.fasta").read_text()))
    for k, v in panel.items():
        if k.split()[0].split("|")[0] == NAV_LABEL or NAV_LABEL in k:
            return v
    raise SystemExit(f"{NAV_LABEL} not in reference_panel.fasta")


def _init(nav: str, seeds: frozenset) -> None:
    global _NAV, _SEEDS
    _NAV, _SEEDS = nav, seeds


def classify_row(rec: dict) -> dict:
    counts = pfam_dict(rec["pfam"])
    project = counts.get("PF00520", 0) >= FOUR_REPEAT_MIN
    q = ChannelQuery(accession=rec["accession"], sequence=rec["sequence"],
                     gene_symbol=rec["gene"], pfam_counts=counts,
                     # 0 TM helices is a measurement, not a missing value
                     # (S2b: a truthiness test here made it None).
                     tm_count=rec["tm_count"] if rec["tm_count"] not in ("", None) else None,
                     length_aa=rec["length"],
                     # UniProt's flag: "Fragment(s)" is a fragment; "Precursor" is not
                     fragment="Fragment" in (rec["fragment"] or ""))
    c = classify(q, refs=None, nav_reference=_NAV if project else "")
    decisive = next((e.tier for e in c.evidence if e.decisive), "")
    motif = next((e.detail for e in c.evidence if e.tier == "motif"), "")
    out = {k: rec.get(k, "") for k in COLS}
    out.update({
        "seed_signatures": ";".join(sorted(a for a in counts if a in _SEEDS)),
        "family": c.family, "superfamily": c.superfamily,
        "status": c.status or ("superfamily_only" if c.superfamily and not c.family
                               else "unassigned"),
        "confidence": c.confidence, "decisive_tier": decisive,
        "filter_string": c.filter_string,
        "k_filter": "yes" if _K_HITS.search(motif) else "",
        "filter_projected": "yes" if project else "",
        "reference_tier": "not_run",
        "hazards": ",".join(c.hazards), "conflicts": ";".join(c.conflicts),
        "notes": c.notes,
        "evidence": " | ".join(f"{e.tier}: {e.detail}" for e in c.evidence),
    })
    return out


def run_shard(key: str, pool: Pool, limit: int | None) -> int:
    dest = raw_dir() / "calls" / f"{key}.tsv.gz"
    if dest.exists() and limit is None:
        print(f"  {key}: done already ({dest.name})")
        return 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    recs = []
    for r in iter_pages(key):
        recs.append(r)
        if limit and len(recs) >= limit:
            break
    t0 = time.time()
    with gzip.open(tmp, "wt", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for i, row in enumerate(pool.imap(classify_row, recs, chunksize=200), 1):
            w.writerow(row)
            if i % 20000 == 0:
                print(f"    {key}: {i:,}/{len(recs):,}  {time.time()-t0:.0f}s", flush=True)
    if limit is None:
        tmp.rename(dest)
    print(f"  {key}: {len(recs):,} classified in {time.time()-t0:.0f}s", flush=True)
    return len(recs)


def recheck_shard(key: str, pool: Pool, accs: frozenset, archive: str) -> tuple[int, int]:
    """Re-classify only the records carrying one of `accs`; copy the rest.

    For a rule change that can only affect records carrying one of a known
    set of accessions (every rewritten rule requires one of them), the other
    calls are pure functions of unchanged evidence and are copied verbatim.
    The previous call file is moved to `calls_<archive>/`, never overwritten.
    """
    src = raw_dir() / "calls" / f"{key}.tsv.gz"
    if (raw_dir() / f"calls_{archive}" / src.name).exists():
        raise SystemExit(f"{key}: calls_{archive}/{src.name} exists — this shard "
                         "was already re-checked; refusing to overwrite the archive")
    old = {}
    with gzip.open(src, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            old[r["accession"]] = r
    recs, order = [], []
    for r in iter_pages(key):
        order.append(r["accession"])
        if set(pfam_dict(r["pfam"])) & accs:
            recs.append(r)
    if len(order) != len(old):
        raise SystemExit(f"{key}: {len(order)} page records vs {len(old)} calls")
    new = {row["accession"]: row for row in pool.imap(classify_row, recs, chunksize=200)}
    arch = raw_dir() / f"calls_{archive}"
    arch.mkdir(exist_ok=True)
    tmp = src.with_suffix(".part")
    with gzip.open(tmp, "wt", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for a in order:
            w.writerow(new.get(a) or old[a])
    src.rename(arch / src.name)
    tmp.rename(src)
    changed = sum(new[a]["family"] != old[a]["family"] or
                  new[a]["superfamily"] != old[a]["superfamily"] for a in new)
    print(f"  {key}: {len(recs):,} re-classified, {changed:,} calls changed", flush=True)
    return len(recs), changed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-j", "--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument("--shard", action="append")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--recheck", default="",
                    help="comma-separated Pfam accessions: re-classify only the "
                         "records carrying one, copy every other call")
    ap.add_argument("--archive", default="r1",
                    help="suffix for the archived previous call files (--recheck)")
    a = ap.parse_args()
    nav = nav_reference()
    if not verify_anchor(nav, FOUR_REPEAT_ANCHOR):
        raise SystemExit("four-repeat anchor does not verify — refusing to run (S1 rule)")
    print(f"anchor {FOUR_REPEAT_ANCHOR.key} verified against {NAV_LABEL}")
    from scripts.s2_lib import signatures
    seeds = frozenset(signatures())
    keys = a.shard or [k for k, _ in SHARDS]
    with Pool(a.jobs, initializer=_init, initargs=(nav, seeds)) as pool:
        for n, key in enumerate(keys):
            live([(f"classify {k}", i < n) for i, k in enumerate(keys)], a.jobs)
            if a.recheck:
                recheck_shard(key, pool, frozenset(a.recheck.split(",")), a.archive)
            else:
                run_shard(key, pool, a.limit)
    live([(f"classify {k}", True) for k in keys], a.jobs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
