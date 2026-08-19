"""S1 — render `results/benchmark_controls/report.md` from the committed tables.

D13 again: the prose here reads the TSVs back and states what they say. If a
re-run changes recall, this report changes with it; if it does not, the
report is stale and the mismatch is visible in git.

    python3 scripts/s1_report.py [--dir results/benchmark_controls]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import read_tsv
from src.catalogue import HAZARD_BY_ID


def table(header, rows) -> str:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def build(d: Path) -> str:
    s = json.loads((d / "summary.json").read_text())
    calls = read_tsv(d / "calls.tsv")
    recall = read_tsv(d / "recall.tsv")
    hazards = read_tsv(d / "hazards.tsv")
    coverage = read_tsv(d / "coverage.tsv")

    L: list[str] = []
    A = L.append
    A("# S1 — classifier control benchmark")
    A("")
    A(f"*Generated {s['generated']} by `scripts/s1_report.py` from the tables "
      f"in `{d.relative_to(ROOT)}` (D13). Panel: `{s['preset']}`. "
      f"{s['n_classified']}/{s['n_targets']} targets classified in "
      f"{s['elapsed_s']} s.*")
    A("")

    A("## Headline")
    A("")
    A(table(["measure", "value"], [
        ["recall (positives called correctly)",
         f"**{s['positives_correct']}/{s['positives']}"
         f"** ({(s['recall'] or 0):.1%})"],
        ["specificity (decoys not called as channels)",
         f"**{s['decoys_rejected']}/{s['decoys']}"
         f"** ({(s['specificity'] or 0):.1%})"],
        ["hazards exercised",
         f"{s['hazards_exercised']}/{s['hazards_defined']}"],
        ["four-repeat anchor validates",
         "yes" if s["anchor_validates"] else "**no — motif tier disabled**"],
        ["MAFFT", s["mafft"] or "**not available**"],
        ["reference panel", f"{s['reference_panel']} exemplar sequences"],
        ["panel proteins that are themselves exemplars",
         f"{s.get('panel_proteins_that_are_exemplars', 0)} — classified "
         f"leave-one-out" if s.get("leave_one_out") else "not excluded"],
    ]))
    A("")
    n_self = s.get("panel_proteins_that_are_exemplars", 0)
    A(f"**Leave-one-out matters here.** {n_self} of the {s['n_classified']} "
      f"panel proteins are themselves catalogue exemplars — both are drawn "
      f"from the same curated gene lists — so without excluding the query's "
      f"own accession the reference tier would score almost the whole panel "
      f"against itself at 100 % identity and the benchmark would measure "
      f"nothing but that overlap (**D29**).")
    A("")

    A("## Which tier made the call")
    A("")
    tiers = Counter(r["decisive_tier"] or "none" for r in calls)
    A("This is the number that says whether the classifier is doing anything "
      "a nearest-neighbour lookup could not. A call made by the architecture "
      "or hazard tier rests on domain composition and works on a sequence "
      "with no close relative in the panel; a call made only by the reference "
      "tier does not.")
    A("")
    A(table(["decisive tier", "n calls"],
            [[k, v] for k, v in tiers.most_common()]))
    A("")
    conf = Counter(r["confidence"] for r in calls)
    A(table(["confidence", "n"], [[k, v] for k, v in conf.most_common()]))
    A("")

    A("## Per-family recall")
    A("")
    pos = [r for r in recall if r.get("panel_role", "positive") == "positive"]
    dec = [r for r in recall if r.get("panel_role") == "decoy"]
    perfect = [r for r in pos if r["recall"] == "1.000"]
    missed = [r for r in pos if r["recall"] != "1.000"]
    A(f"{len(perfect)}/{len(pos)} **positive** families were called correctly "
      f"for every member. The {len(dec)} decoy families are scored by "
      f"specificity, not recall — a decoy left `unassigned` is a success, and "
      f"listing it as a recall failure would invert the result.")
    A("")
    if missed:
        A("**Positive families the classifier did not call correctly:**")
        A("")
        A(table(["family", "n", "correct", "recall", "called instead"],
                [[r["family"], r["n"], r["correct"], r["recall"],
                  r["mis_calls"] or "—"] for r in missed]))
        A("")
    if dec:
        wrong_dec = [r for r in dec if r["mis_calls"] and
                     "unassigned" not in r["mis_calls"]]
        A(f"**Decoys:** {len(dec)} families; "
          f"{len(dec) - len(wrong_dec)} were left unassigned or called as "
          f"themselves, {len(wrong_dec)} were called as something else.")
        if wrong_dec:
            A("")
            A(table(["decoy family", "called instead"],
                    [[r["family"], r["mis_calls"]] for r in wrong_dec]))
        A("")
    A("Full table: `recall.tsv`.")
    A("")

    A("## Per-hazard")
    A("")
    A("A hazard is only closed by a test. `n touched` counts panel members "
      "whose classification involved that hazard's rule; a hazard with zero "
      "is not solved, it is **untested** by this panel.")
    A("")
    rows = []
    for h in hazards:
        hz = HAZARD_BY_ID.get(h["hazard"])
        rows.append([
            h["hazard"], h["severity"], h["n_touched"], h["n_correct"],
            h["n_wrong"], h["wrong_examples"] or "—",
            (hz.title if hz else "")[:56]])
    A(table(["hazard", "severity", "n touched", "correct", "wrong",
             "examples", "what it is"], rows))
    A("")
    untested = [h for h in hazards if h["n_touched"] == "0"]
    if untested:
        A(f"**{len(untested)} hazard(s) carry no measurement from this panel:** "
          + ", ".join(h["hazard"] for h in untested)
          + ". Extending the panel to exercise them is a task, not a footnote.")
        A("")

    A("## Panel coverage")
    A("")
    covered = [r for r in coverage if r["n_in_panel"] != "0"]
    A(f"{len(covered)}/{len(coverage)} census families have at least one "
      f"member in this panel. The families with none are mostly the "
      f"non-human ones — prokaryotic channels, invertebrate degenerins, "
      f"plant OSCA — and a benchmark that only measures human proteins "
      f"cannot claim anything about a census that is not human-only.")
    A("")
    uncovered = [r["family"] for r in coverage if r["n_in_panel"] == "0"]
    if uncovered:
        A("Families with no panel member: " + ", ".join(f"`{u}`" for u in uncovered))
        A("")

    A("## Every call")
    A("")
    A("`calls.tsv` carries one row per protein with the full evidence string "
      "from every tier, including the tiers that did not decide. The point of "
      "keeping the losing evidence is that a wrong call is diagnosable "
      "without re-running anything.")
    A("")
    interesting = [r for r in calls if r["conflicts"]]
    if interesting:
        A(f"**{len(interesting)} call(s) had disagreeing tiers:**")
        A("")
        A(table(["protein", "gene", "expected", "called", "confidence",
                 "conflict"],
                [[r["accession"], r["gene_symbol"], r["expected_family"],
                  r["called_family"] or "—", r["confidence"],
                  r["conflicts"][:70]] for r in interesting]))
        A("")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", type=Path,
                    default=ROOT / "results" / "benchmark_controls")
    args = ap.parse_args()
    if not (args.dir / "summary.json").exists():
        print(f"[s1_report] no summary.json in {args.dir} — run "
              f"scripts/s1_benchmark.py first")
        return 1
    out = args.dir / "report.md"
    out.write_text(build(args.dir))
    print(f"[s1_report] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
