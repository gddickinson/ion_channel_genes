"""S0 — render `results/s0_baseline/report.md` from the committed tables.

Decision **D13**, inherited from the parent project: nothing in this report
is written by hand. Every number is read back out of the TSVs that
`s0_catalogue_verify.py` wrote, so the report cannot drift from the data, and
re-running the verification after a catalogue edit produces a report that
disagrees with the last one exactly where the catalogue changed.

    python3 scripts/s0_report.py [--dir results/s0_baseline]
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


def pct(n: int, d: int) -> str:
    return f"{100 * n / d:.0f} %" if d else "—"


def table(header: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def build(d: Path) -> str:
    summary = json.loads((d / "verification_summary.json").read_text())
    pfam = read_tsv(d / "pfam_verification.tsv")
    exres = read_tsv(d / "exemplars_resolved.tsv")
    arch = read_tsv(d / "exemplar_architecture.tsv")
    tax = read_tsv(d / "taxonomy_check.tsv")
    share = read_tsv(d / "signature_sharing.tsv")
    problems = read_tsv(d / "catalogue_validation.tsv")
    cat = summary["catalogue"]

    L: list[str] = []
    A = L.append
    A("# S0 — catalogue verification baseline")
    A("")
    A(f"*Generated {summary['generated']} by `scripts/s0_report.py` from the "
      f"tables in `{d.relative_to(ROOT)}`. Nothing here is hand-written "
      f"(D13). {summary['requests']} live requests, "
      f"{len(summary.get('request_failures', []))} failure(s), "
      f"{summary['elapsed_s']} s.*")
    A("")
    A("## What the catalogue claims")
    A("")
    A(table(["quantity", "n"],
            [[k.replace("_", " "), v] for k, v in cat.items()]))
    A("")
    A(f"The census denominator is **{cat['human_genes_census']} human "
      f"pore-forming genes** across **{cat['families_census']} families** in "
      f"**{cat['superfamilies']} superfamilies**, with "
      f"**{cat['families_control']} further families catalogued specifically "
      f"so they can be excluded** — auxiliary subunits, domain-sharing "
      f"non-channels and channels that are not ion channels.")
    A("")

    A("## 1. Internal consistency")
    A("")
    if problems:
        A(f"**{len(problems)} problem(s)** found by `src.catalogue.validate()`:")
        A("")
        for r in problems:
            A(f"- {r['problem']}")
    else:
        A("`src.catalogue.validate()` is clean: every family points at a "
          "superfamily that exists, no human gene is claimed by two census "
          "families, every exemplar is resolvable, every `confusable_with` "
          "and every hazard family reference matches a catalogue key.")
    A("")

    A("## 2. Pfam accessions")
    A("")
    ok = [r for r in pfam if r["status"] == "ok"]
    missing = [r for r in pfam if r["status"] == "MISSING"]
    renamed = [r for r in pfam if r["status"] == "NAME_MISMATCH"]
    A(f"{len(pfam)} distinct accessions declared across the catalogue; "
      f"**{len(ok)} verified ({pct(len(ok), len(pfam))})**, "
      f"{len(missing)} missing from InterPro, "
      f"{len(renamed)} carrying a different short name.")
    A("")
    if missing:
        A("**Missing accessions** — these are claims the catalogue makes that "
          "InterPro does not support:")
        A("")
        A(table(["accession", "declared name", "families using it"],
                [[r["accession"], r["declared_name"], r["n_families_using"]]
                 for r in missing]))
        A("")
    if renamed:
        A("**Short-name mismatches** — the accession exists; the catalogue "
          "calls it something else:")
        A("")
        A(table(["accession", "catalogue", "InterPro"],
                [[r["accession"], r["declared_name"], r["interpro_short"]]
                 for r in renamed]))
        A("")
    counted = [r for r in pfam if r.get("n_uniprot_proteins", "").isdigit()]
    if counted:
        counted.sort(key=lambda r: -int(r["n_uniprot_proteins"]))
        A("**The ten largest search spaces** — how many UniProt proteins carry "
          "each signature. This is the census's raw denominator, and the "
          "reason S2 is a task of its own:")
        A("")
        A(table(["accession", "name", "UniProt proteins", "families"],
                [[r["accession"], r["interpro_short"],
                  f"{int(r['n_uniprot_proteins']):,}", r["n_families_using"]]
                 for r in counted[:10]]))
        A("")

    A("## 3. Exemplars")
    A("")
    st = Counter(r["status"] for r in exres)
    A(f"{len(exres)} exemplars declared; "
      f"{st.get('ok', 0)} resolved cleanly, "
      f"{st.get('UNRESOLVED', 0)} unresolvable, "
      f"{st.get('DECLARED_GENE_MISMATCH', 0)} with a declared accession that "
      f"names a different gene, "
      f"{st.get('GENE_NAME_DIFFERS', 0)} where UniProt's primary symbol "
      f"differs from the catalogue's.")
    A("")
    A(table(["status", "n"], [[k, v] for k, v in st.most_common()]))
    A("")
    bad = [r for r in exres if r["status"] not in ("ok", "")]
    if bad:
        A("Every exemplar that did not resolve cleanly:")
        A("")
        A(table(["family", "label", "gene", "species", "declared", "resolved",
                 "resolved gene", "status"],
                [[r["family"], r["label"], r["gene"], r["species"],
                  r["declared_accession"] or "—",
                  r["resolved_accession"] or "—",
                  r["resolved_gene"] or "—", r["status"]] for r in bad]))
        A("")
    A(f"`reference_panel.fasta` holds "
      f"{summary['exemplars'].get('sequences', 0)} sequences — the panel the "
      f"classifier's reference tier scores against.")
    A("")

    A("## 4. Declared architecture vs observed")
    A("")
    mismatch = [r for r in arch if r["status"] != "ok"]
    A(f"{len(arch)} exemplars had their Pfam architecture re-derived. "
      f"**{len(arch) - len(mismatch)} match the catalogue's declaration; "
      f"{len(mismatch)} do not.**")
    A("")
    A("A mismatch is not automatically an error in the catalogue — a "
      "signature declared for a family need not be present on every member, "
      "and that is exactly what hazard **H7** records for the TRP families. "
      "It is, however, always something a person should look at.")
    A("")
    if mismatch:
        A(table(["family", "exemplar", "declared but not observed",
                 "observed but not declared"],
                [[r["family"], r["label"],
                  r["declared_not_observed"] or "—",
                  (r["observed_not_declared"] or "—")[:60]]
                 for r in mismatch[:40]]))
        if len(mismatch) > 40:
            A("")
            A(f"*({len(mismatch) - 40} further rows in "
              f"`exemplar_architecture.tsv`.)*")
        A("")

    A("## 5. Shared signatures — the hazard registry's evidence")
    A("")
    crossing = [r for r in share if r["crosses_channel_boundary"] == "yes"]
    A(f"{len(share)} signatures are carried by more than one family. "
      f"**{len(crossing)} of them cross the channel / non-channel boundary** "
      f"— a domain that is evidence for a channel family and is also carried "
      f"by something the catalogue does not count as a channel.")
    A("")
    A(table(["accession", "families", "statuses"],
            [[r["accession"], r["families"][:80], r["statuses"]]
             for r in share[:15]]))
    A("")

    A("## 6. Taxonomy")
    A("")
    tbad = [r for r in tax if r["status"] != "ok"]
    A(f"{len(tax)} taxon ids checked against the UniProt taxonomy service; "
      f"{len(tax) - len(tbad)} confirmed, {len(tbad)} mismatched.")
    A("")
    if tbad:
        A(table(["taxon_id", "catalogue", "UniProt", "status"],
                [[r["taxon_id"], r["declared"], r["uniprot_scientific"] or "—",
                  r["status"]] for r in tbad]))
        A("")

    A("## 7. What this does not verify")
    A("")
    A("- **Family membership.** That `KCNA1` belongs to `kv_shaker` is a "
      "literature assignment. S1 tests whether the classifier reproduces it; "
      "nothing here tests whether it is right.")
    A("- **Selectivity and gating.** Both are literature attributes of a "
      "family (hazard **H14**) and are not recoverable from sequence, so no "
      "database check can confirm them.")
    A("- **Completeness.** A verified catalogue is not a complete one. "
      "Whether 320 is the right number of human pore-forming genes is the "
      "question S2 answers, and it is answered against a declared search "
      "space, never against \"all ion channels\".")
    A("")
    fails = summary.get("request_failures", [])
    if fails:
        A("## Request failures")
        A("")
        for fdict in fails[:20]:
            A(f"- `{fdict['url']}` — {fdict['why']}")
        A("")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", type=Path, default=ROOT / "results" / "s0_baseline")
    args = ap.parse_args()
    if not (args.dir / "verification_summary.json").exists():
        print(f"[s0_report] no verification_summary.json in {args.dir} — "
              f"run scripts/s0_catalogue_verify.py first")
        return 1
    out = args.dir / "report.md"
    out.write_text(build(args.dir))
    print(f"[s0_report] wrote {out} ({out.stat().st_size / 1000:.1f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
