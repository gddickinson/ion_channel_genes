"""S1 — measure the classifier against positives and the decoys that look like them.

The parent project's equivalent benchmark is the reason this project has a
hazard registry at all: on its first run, every one of six ryanodine-receptor
decoys was promoted as a candidate IP3 receptor, because RyRs carry all four
IP3R-diagnostic Pfam domains. Specificity was 25/31 and every failure was the
same mistake.

Here there are sixteen such mistakes available, so the benchmark is built
around them. `presets/controls_benchmark.json` names a positive from every
census family with a resolvable human gene, and a decoy for every hazard —
KCTD1 for the Kv T1 domain, GRM1 for the iGluR clamshell, POMT1 for MIR,
TPTE for the pore-module annotation on a phosphatase, ABCC8 for CFTR's
architecture. The tables answer four questions:

    recall      does the classifier find each family's own member?
    precision   does it reject each decoy?
    per hazard  did the rule that closes each hazard actually fire?
    per tier    which tier made each correct call — is the architecture tier
                doing the work, or is the reference tier carrying it?

The last one matters most. A classifier whose calls all come from
reference identity is a nearest-neighbour lookup with extra steps, and it
will fail on exactly the sequences a census cares about: the ones with no
close relative in the panel.

    python3 scripts/s1_benchmark.py [--preset controls_benchmark]
                                    [--limit N] [--no-reference]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import live_progress
from src.catalogue import HAZARDS
from src.classify import ReferenceSet, classify_all
from src.classify.motifs import FOUR_REPEAT_ANCHOR, verify_anchor
from src.classify.reference import fetch_uniprot_sequence
from src.classify.report import (CALLS_HEADER, calls_table, confusion_table,
                                 coverage_table, hazard_table, recall_table,
                                 summary_counts, write_tsv)
from src.cli_channel import build_query, resolve_symbol
from src.utils.mafft import mafft_available, mafft_version


def log(m: str) -> None:
    print(f"[s1] {m}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="controls_benchmark")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "benchmark_controls")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-reference", action="store_true")
    args = ap.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    live = ROOT / "results" / "session_live.json"
    t0 = time.time()

    preset = json.loads((ROOT / "presets" / f"{args.preset}.json").read_text())
    targets = list(preset["classify"])
    expected_by_symbol = dict(preset["expected"])
    positives = set(preset.get("positives", []))
    if args.limit:
        targets = targets[:args.limit]

    # -- gate 1: does the four-repeat anchor still validate? --------------
    log("verifying the four-repeat filter anchor "
        f"({FOUR_REPEAT_ANCHOR.reference_label} "
        f"{FOUR_REPEAT_ANCHOR.reference_uniprot} "
        f"{FOUR_REPEAT_ANCHOR.positions} == {FOUR_REPEAT_ANCHOR.expected})")
    nav_ref = fetch_uniprot_sequence(FOUR_REPEAT_ANCHOR.reference_uniprot)
    anchor_ok = bool(nav_ref) and verify_anchor(nav_ref)
    log(f"  anchor {'validates' if anchor_ok else 'DOES NOT VALIDATE'} "
        f"(reference {len(nav_ref)} aa)")
    if not anchor_ok:
        log("  → the motif tier is disabled for this run; the filter columns "
            "in the tables will be empty rather than wrong")
        nav_ref = ""

    # -- gate 2: is MAFFT there? -----------------------------------------
    log(f"mafft: {mafft_version() or 'NOT AVAILABLE'}")
    if not mafft_available():
        log("  → both the motif and reference tiers are unavailable; this "
            "run measures the architecture tier alone (D28)")

    # -- reference panel --------------------------------------------------
    refs = None
    if not args.no_reference:
        cache = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
        if cache.exists():
            refs = ReferenceSet.from_fasta(cache)
            log(f"reference panel: {len(refs.sequences)} exemplars from "
                f"{cache.relative_to(ROOT)}")
        else:
            log(f"no reference panel at {cache.relative_to(ROOT)} — run "
                f"scripts/s0_catalogue_verify.py first")

    # -- fetch evidence ---------------------------------------------------
    steps = [(t, False) for t in targets]
    live_progress(live, "S1", steps)
    queries, expected, alias = [], {}, {}
    for i, sym in enumerate(targets):
        acc = sym if (len(sym) in (6, 10) and sym[1:2].isdigit()) \
            else resolve_symbol(sym)
        if not acc:
            log(f"  [{i+1}/{len(targets)}] {sym}: UNRESOLVED")
            continue
        q = build_query(acc)
        alias[acc] = sym
        if sym in expected_by_symbol:
            expected[acc] = expected_by_symbol[sym]
        queries.append(q)
        steps[i] = (sym, True)
        if (i + 1) % 5 == 0:
            live_progress(live, "S1", steps)
        log(f"  [{i+1}/{len(targets)}] {sym:10s} {acc:10s} "
            f"{q.length_aa or '?':>5} aa  {len(q.pfam_counts)} pfam  "
            f"tm={q.tm_count}")
    live_progress(live, "S1", steps)

    # -- classify ---------------------------------------------------------
    n_self = 0
    if refs:
        panel_accs = {q.accession for q in queries}
        n_self = len(panel_accs & set(refs.accession_of.values()))
        log(f"{n_self}/{len(queries)} panel proteins are themselves catalogue "
            f"exemplars — classifying leave-one-out so the reference tier "
            f"cannot score a protein against itself")
    log(f"classifying {len(queries)} protein(s)… "
        f"(each is up to {12} MAFFT alignments against the panel)")
    calls = classify_all(queries, refs, nav_ref, leave_one_out=True,
                         on_progress=lambda m: log(f"  {m}"))
    for c in calls:
        c.gene_symbol = alias.get(c.query, c.gene_symbol)

    # -- tables -----------------------------------------------------------
    write_tsv(out / "calls.tsv", CALLS_HEADER, calls_table(calls, expected))
    write_tsv(out / "confusion.tsv", ("expected", "called", "n"),
              confusion_table(calls, expected))
    write_tsv(out / "recall.tsv",
              ("family", "name", "n", "correct", "recall", "decisive_tiers",
               "confidence", "mis_calls"), recall_table(calls, expected))
    write_tsv(out / "hazards.tsv",
              ("hazard", "severity", "families", "n_touched", "n_correct",
               "n_wrong", "wrong_examples", "test_owner", "discriminator"),
              hazard_table(calls, expected))
    write_tsv(out / "coverage.tsv",
              ("family", "superfamily", "human_genes", "n_in_panel"),
              coverage_table(calls))

    pos_calls = [c for c in calls if alias.get(c.query) in positives]
    dec_calls = [c for c in calls if alias.get(c.query) not in positives]
    correct_pos = sum(1 for c in pos_calls if c.family == expected.get(c.query))
    rejected = sum(1 for c in dec_calls if not c.census_member())

    summary = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "preset": args.preset,
        "anchor_validates": anchor_ok,
        "mafft": mafft_version(),
        "reference_panel": len(refs.sequences) if refs else 0,
        "panel_proteins_that_are_exemplars": n_self,
        "leave_one_out": True,
        "n_targets": len(targets),
        "n_classified": len(calls),
        "positives": len(pos_calls),
        "positives_correct": correct_pos,
        "recall": round(correct_pos / len(pos_calls), 4) if pos_calls else None,
        "decoys": len(dec_calls),
        "decoys_rejected": rejected,
        "specificity": round(rejected / len(dec_calls), 4) if dec_calls else None,
        "hazards_defined": len(HAZARDS),
        "hazards_exercised": sum(
            1 for h in HAZARDS if any(h.hid in c.hazards for c in calls)),
        "elapsed_s": round(time.time() - t0, 1),
    }
    summary.update({f"counts_{k}": v
                    for k, v in summary_counts(calls, expected).items()})
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    log("")
    log(f"recall      {correct_pos}/{len(pos_calls)} "
        f"({summary['recall'] or 0:.1%}) on positives")
    log(f"specificity {rejected}/{len(dec_calls)} "
        f"({summary['specificity'] or 0:.1%}) — decoys not called as channels")
    log(f"hazards     {summary['hazards_exercised']}/{len(HAZARDS)} exercised "
        f"by this panel")
    log(f"done in {summary['elapsed_s']}s → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
