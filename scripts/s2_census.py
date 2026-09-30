"""S2 step 3 — assemble census v2 and the tables its report is rendered from.

    python3 scripts/s2_census.py

Streams the per-shard call files once. Bulk output goes to the data root
(too large for the repo — ~1.25 M rows); the repo keeps the tables and a
manifest with the SHA-256 of each bulk file, so the committed numbers can be
tied to the exact files they came from.

Data root, `<root>/raw_api/s2/`:
    census_v2.tsv.gz     one row per record, every call column + provenance
    census_v2.fasta.gz   every record's sequence, header `acc|family|status|taxon`

`results/census_v2/`:
    census_status.tsv          status × confidence
    census_families.tsv        per family: records, reviewed, by domain, human, confidence
    census_superfamily_only.tsv  records reaching a superfamily and no family
    tier_attribution.tsv       which tier decided each family-level call
    unassigned_signatures.tsv  what brought the unassigned records in
    four_repeat_filters.tsv    the projected filter strings and what they called
    human_recall.tsv           the 320 human census genes: found? called right?
    s1_panel.tsv               every S1 panel protein: S1 call vs S2 call
    partial_architectures.tsv  unassigned records carrying part of exactly one
                               family's declared architecture
    manifest.tsv, summary.json
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s0_lib import read_tsv, write_tsv                        # noqa: E402
from scripts.s2_classify import COLS                                  # noqa: E402
from scripts.s2_lib import CENSUS_SHARDS as SHARDS, OUT_DIR, ROOT, iter_pages, live, raw_dir  # noqa: E402
from src.catalogue import CATALOGUE, Level                            # noqa: E402
from src.catalogue.registry import human_genes                        # noqa: E402

DOMAINS = ("Eukaryota", "Bacteria", "Archaea", "Viruses")
CONF = ("gold", "silver", "bronze", "unassigned")
PROV = ["source", "release", "shard"]
# r1 = S2 as run 2026-09-28. r2 = S2b: hazard rules H2, H4, H13 rewritten as
# positive tests; the 108,407 records carrying PF02931/PF08709/PF08016/PF20519
# re-classified (`s2_classify.py --recheck`), r1 calls and bulk archived.
# r3 = S2c: H11, H12 likewise; PF02214/PF00664 carriers re-classified, r2
# calls archived under calls_r2/.
REVISION = "r3 (S2c, 2026-09-28: H2/H4/H11/H12/H13 positive tests)"

#: {family: its FAMILY-level accessions} — what a derived rule requires in full.
FAMILY_SIGS = {k: {s.accession for s in f.signatures if s.level is Level.FAMILY}
               for k, f in CATALOGUE.items()}


def partial_candidates(seeds: str) -> list[str]:
    """Families whose declared architecture this record carries *part* of."""
    have = set(filter(None, seeds.split(";")))
    return sorted(k for k, sig in FAMILY_SIGS.items() if sig & have)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def iter_calls():
    for key, _ in SHARDS:
        p = raw_dir() / "calls" / f"{key}.tsv.gz"
        if not p.exists():
            raise SystemExit(f"missing {p} — run s2_classify.py first")
        with gzip.open(p, "rt", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["shard"] = key
                yield row


def main() -> int:
    release = json.loads((OUT_DIR / "release.json").read_text())["uniprot_release"]
    rd = raw_dir()
    status_conf: Counter = Counter()
    fam = defaultdict(Counter)
    sfo = defaultdict(Counter)
    tiers: Counter = Counter()
    unassigned_combo: Counter = Counter()
    unassigned_sig: Counter = Counter()
    partial: Counter = Counter()
    filters: Counter = Counter()
    projected = 0
    human: dict[str, list[dict]] = defaultdict(list)
    s1 = {r["accession"]: r for r in read_tsv(ROOT / "results" /
                                               "benchmark_controls" / "calls.tsv")}
    s1_hits: dict[str, dict] = {}
    seen: set[str] = set()
    n = 0
    live([("assemble census_v2", False), ("write FASTA", False)])
    # The call files do not carry the gene symbol (the classifier never reads
    # it — H15); it is joined back here from the archived records so a reader
    # can see it and the human-recall table can score by it.
    genes = {rec["accession"]: rec["gene"] for key, _ in SHARDS
             for rec in iter_pages(key) if rec["gene"]}
    out_tsv = rd / "census_v2.tsv.gz"
    with gzip.open(out_tsv, "wt", newline="") as fh:
        w = csv.DictWriter(fh, COLS + ["gene"] + PROV, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in iter_calls():
            if r["accession"] in seen:
                continue
            seen.add(r["accession"])
            n += 1
            r["source"], r["release"] = "UniProtKB", release
            r["gene"] = genes.get(r["accession"], "")
            w.writerow(r)
            st, cf = r["status"], r["confidence"]
            status_conf[(st, cf)] += 1
            dom = r["domain"] if r["domain"] in DOMAINS else "other"
            if r["family"]:
                c = fam[r["family"]]
                c["records"] += 1
                c[dom] += 1
                c[cf] += 1
                c["reviewed"] += r["reviewed"] == "reviewed"
                c["fragment"] += bool(r["fragment"])
                c["human"] += r["taxon_id"] == "9606"
                tiers[(r["decisive_tier"], st)] += 1
            elif st == "superfamily_only":
                c = sfo[r["superfamily"]]
                c["records"] += 1
                c["ambiguous"] += r["notes"].startswith("ambiguous between")
                c[dom] += 1
                c["reviewed"] += r["reviewed"] == "reviewed"
            else:
                unassigned_combo[r["seed_signatures"]] += 1
                for a in r["seed_signatures"].split(";"):
                    unassigned_sig[(a, dom)] += 1
                cands = partial_candidates(r["seed_signatures"])
                partial[cands[0] if len(cands) == 1 else
                        ("(several)" if cands else "(none)")] += 1
            if r["filter_projected"]:
                projected += 1
                filters[(r["filter_string"], r["family"])] += 1
            if r["taxon_id"] == "9606" and r["reviewed"] == "reviewed" and r["gene"]:
                human[r["gene"].upper()].append(r)
            if r["accession"] in s1:
                s1_hits[r["accession"]] = r
    live([("assemble census_v2", True), ("write FASTA", False)])

    out_fa = rd / "census_v2.fasta.gz"
    calls = {}
    with gzip.open(out_tsv, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            calls[r["accession"]] = (r["family"] or r["superfamily"] or "-", r["status"])
    written: set[str] = set()
    with gzip.open(out_fa, "wt") as fh:
        for key, _ in SHARDS:
            for rec in iter_pages(key):
                a = rec["accession"]
                if a in written:
                    continue
                written.add(a)
                f, st = calls.get(a, ("-", "-"))
                fh.write(f">{a}|{f}|{st}|{rec['taxon_id']}\n{rec['sequence']}\n")
    live([("assemble census_v2", True), ("write FASTA", True)])

    write_tables(n, status_conf, fam, sfo, tiers, unassigned_combo, unassigned_sig,
                 filters, projected, human, s1, s1_hits, release)
    write_tsv(OUT_DIR / "partial_architectures.tsv",
              ["family", "declared_family_signatures", "unassigned_records"],
              [[k, ",".join(sorted(FAMILY_SIGS.get(k, ()))), v]
               for k, v in sorted(partial.items(), key=lambda x: -x[1])])
    write_tsv(OUT_DIR / "manifest.tsv", ["file", "bytes", "sha256", "records"], [
        [str(p.relative_to(rd.parents[1])), p.stat().st_size, sha256(p), n]
        for p in (out_tsv, out_fa)])
    print(f"census_v2: {n:,} records → {out_tsv}")
    return 0


def write_tables(n, status_conf, fam, sfo, tiers, uc, us, filters, projected,
                 human, s1, s1_hits, release) -> None:
    write_tsv(OUT_DIR / "census_status.tsv", ["status", "confidence", "records"],
              [[s, c, k] for (s, c), k in sorted(status_conf.items(),
                                                   key=lambda x: -x[1])])
    rows = []
    for key, c in sorted(fam.items(), key=lambda x: -x[1]["records"]):
        f = CATALOGUE.get(key)
        rows.append([key, f.superfamily if f else "", f.status.value if f else "",
                     c["records"], c["reviewed"], c["fragment"], c["human"],
                     *[c[d] for d in DOMAINS], c["other"], *[c[x] for x in CONF[:3]]])
    write_tsv(OUT_DIR / "census_families.tsv",
              ["family", "superfamily", "catalogue_status", "records", "reviewed",
               "fragment", "human", *[d.lower() for d in DOMAINS], "other",
               *CONF[:3]], rows)
    write_tsv(OUT_DIR / "census_superfamily_only.tsv",
              ["superfamily", "records", "ambiguous", "reviewed",
               *[d.lower() for d in DOMAINS], "other"],
              [[k, c["records"], c["ambiguous"], c["reviewed"],
                *[c[d] for d in DOMAINS], c["other"]]
               for k, c in sorted(sfo.items(), key=lambda x: -x[1]["records"])])
    write_tsv(OUT_DIR / "tier_attribution.tsv", ["decisive_tier", "status", "records"],
              [[t, s, k] for (t, s), k in sorted(tiers.items(), key=lambda x: -x[1])])
    by_sig = defaultdict(Counter)
    for (a, d), k in us.items():
        by_sig[a][d] += k
    write_tsv(OUT_DIR / "unassigned_signatures.tsv",
              ["pfam", "records", *[d.lower() for d in DOMAINS], "other", "solo_records"],
              [[a, sum(c.values()), *[c[d] for d in DOMAINS], c["other"], uc.get(a, 0)]
               for a, c in sorted(by_sig.items(), key=lambda x: -sum(x[1].values()))])
    write_tsv(OUT_DIR / "four_repeat_filters.tsv", ["filter_string", "family", "records"],
              [[s or "(none)", f or "-", k] for (s, f), k in
               sorted(filters.items(), key=lambda x: -x[1])])
    hrows, found, right = [], 0, 0
    for g, want in sorted(human_genes().items(), key=lambda x: (x[1], x[0])):
        recs = human.get(g.upper(), [])
        best = next((r for r in recs if r["family"] == want), recs[0] if recs else None)
        ok = bool(best and best["family"] == want)
        found += bool(best)
        right += ok
        hrows.append([g, want, best["accession"] if best else "",
                      best["family"] if best else "", best["status"] if best else "absent",
                      best["confidence"] if best else "", best["decisive_tier"] if best else "",
                      "yes" if ok else "no"])
    write_tsv(OUT_DIR / "human_recall.tsv",
              ["gene", "expected_family", "accession", "s2_family", "s2_status",
               "confidence", "decisive_tier", "correct"], hrows)
    prows, s2_right, s1_right = [], 0, 0
    for a, r in s1.items():
        h = s1_hits.get(a)
        exp = r["expected_family"]
        s1_ok, s2_ok = r["called_family"] == exp, bool(h and h["family"] == exp)
        s1_right += s1_ok
        s2_right += s2_ok
        prows.append([a, r["gene_symbol"], exp, r["called_family"],
                      h["family"] if h else "", h["status"] if h else "not_enumerated",
                      h["decisive_tier"] if h else "", "yes" if s1_ok else "no",
                      "yes" if s2_ok else "no"])
    write_tsv(OUT_DIR / "s1_panel.tsv",
              ["accession", "gene", "expected_family", "s1_family", "s2_family",
               "s2_status", "s2_decisive_tier", "s1_correct", "s2_correct"], prows)
    fam_calls = sum(c["records"] for c in fam.values())
    no_ref = sum(k for (t, _), k in tiers.items() if t in ("architecture", "hazard"))
    (OUT_DIR / "summary.json").write_text(json.dumps({
        "uniprot_release": release, "records": n, "family_calls": fam_calls,
        "revision": REVISION,
        "family_calls_channel": sum(k for (t, s), k in tiers.items() if s == "channel"),
        "superfamily_only": sum(c["records"] for c in sfo.values()),
        "unassigned": sum(k for (s, _), k in status_conf.items() if s == "unassigned"),
        "family_calls_by_architecture_or_hazard": no_ref,
        "family_calls_by_motif": sum(k for (t, _), k in tiers.items() if t == "motif"),
        "four_repeat_projected": projected,
        "human_genes": len(human_genes()), "human_found": found, "human_correct": right,
        "s1_panel": len(s1), "s1_enumerated": len(s1_hits),
        "s1_correct_in_s1": s1_right, "s1_correct_in_s2": s2_right,
    }, indent=2))


if __name__ == "__main__":
    sys.exit(main())
