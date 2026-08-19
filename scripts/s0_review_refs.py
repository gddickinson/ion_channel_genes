"""Resolve the review's source titles against Europe PMC → references.tsv.

Every reference in `docs/channel_review_2026.md` passes through here first.
The rule the whole review rests on: **a citation exists only if a live
literature service returned a record whose title matches the one we asked
for.** A misremembered title fails to resolve and is reported; it never
becomes a plausible-looking entry in a bibliography.

Matching is deliberately strict. Europe PMC will happily return *something*
for almost any query, so a result is accepted only when its title, stripped
of punctuation and case, is at least `MIN_TITLE_MATCH` similar to the query
title. A near-miss is reported as `TITLE_MISMATCH` with what came back, so a
person can decide whether the source list has a typo or the paper is not the
one intended.

    python3 scripts/s0_review_refs.py [--out results/s0_baseline]
                                      [--limit N] [--recheck]

Output: `references.tsv` (key, pmid, doi, year, journal, authors, title,
status) and a summary line. Exits non-zero if anything failed to resolve, so
it can gate the review build.
"""

from __future__ import annotations

import argparse
import difflib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from review_sources import SOURCES
from scripts.s0_lib import read_tsv, write_tsv

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

#: How close the returned title must be to the requested one.
MIN_TITLE_MATCH = 0.90

HEADER = ("key", "pmid", "doi", "year", "journal", "authors", "title",
          "status", "match")


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9 ]+", " ", (s or "").lower()).strip()


def similarity(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, norm(a), norm(b)).ratio()


def search(title: str, extra: str = "", page_size: int = 5,
           timeout_s: int = 40) -> list[dict]:
    q = f'TITLE:"{title}"' + (f" AND ({extra})" if extra else "")
    url = (f"{EPMC}?query=" + urllib.parse.quote(q)
           + f"&format=json&pageSize={page_size}")
    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as r:
            return json.load(r)["resultList"]["result"]
    except Exception:
        return []


def search_loose(title: str, timeout_s: int = 40) -> list[dict]:
    """Fallback: unquoted query, for titles with punctuation EPMC dislikes."""
    url = (f"{EPMC}?query=" + urllib.parse.quote(title)
           + "&format=json&pageSize=5")
    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as r:
            return json.load(r)["resultList"]["result"]
    except Exception:
        return []


def best(results: list[dict], title: str) -> tuple[dict | None, float]:
    scored = [(similarity(title, r.get("title", "")), r) for r in results]
    if not scored:
        return None, 0.0
    scored.sort(key=lambda t: -t[0])
    return scored[0][1], scored[0][0]


def resolve(key: str, title: str, extra: str = "",
            sleep_s: float = 0.25) -> tuple:
    rec, score = best(search(title, extra), title)
    if rec is None or score < MIN_TITLE_MATCH:
        time.sleep(sleep_s)
        rec2, score2 = best(search_loose(title), title)
        if score2 > score:
            rec, score = rec2, score2
    time.sleep(sleep_s)

    if rec is None:
        return (key, "", "", "", "", "", title, "NOT_FOUND", "0.00")
    status = "ok" if score >= MIN_TITLE_MATCH else "TITLE_MISMATCH"
    authors = rec.get("authorString", "")
    # Nature-style short author form: "Doyle, D. A. et al."
    if authors.count(",") > 2:
        authors = authors.split(",")[0].strip() + " et al."
    return (key, rec.get("pmid", ""), rec.get("doi", ""), rec.get("pubYear", ""),
            rec.get("journalTitle", ""), authors, rec.get("title", "").rstrip("."),
            status, f"{score:.2f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "s0_baseline")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--recheck", action="store_true",
                    help="re-resolve keys already marked ok (default: keep them)")
    args = ap.parse_args()

    path = args.out / "references.tsv"
    existing = {r["key"]: r for r in read_tsv(path)}
    sources = SOURCES[:args.limit] if args.limit else SOURCES

    rows, failed = [], []
    for i, src in enumerate(sources, 1):
        key, title = src[0], src[1]
        extra = src[2] if len(src) > 2 else ""
        prev = existing.get(key)
        if prev and prev.get("status") == "ok" and not args.recheck:
            rows.append(tuple(prev.get(h, "") for h in HEADER))
            continue
        row = resolve(key, title, extra)
        rows.append(row)
        if row[7] != "ok":
            failed.append((key, row[7], row[6][:60]))
        if i % 10 == 0:
            print(f"[refs] {i}/{len(sources)}  ({len(failed)} unresolved)",
                  flush=True)

    write_tsv(path, HEADER, rows)
    ok = sum(1 for r in rows if r[7] == "ok")
    print(f"[refs] {ok}/{len(rows)} resolved → {path}")
    if failed:
        print(f"[refs] {len(failed)} did NOT resolve — fix the title in "
              f"scripts/review_sources.py or drop the source:")
        for key, status, got in failed:
            print(f"   {key:22s} {status:16s} got: {got}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
