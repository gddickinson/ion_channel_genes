"""s6_report.py — renders results/alignments/report.md from the S6 tables (D13).

    python3 scripts/s6_report.py
"""

from __future__ import annotations

import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv  # noqa: E402

D = ROOT / "results" / "alignments"


def _t(header: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in rows]
    return "\n".join(out)


def sets_section() -> list[str]:
    m = read_tsv(D / "members.tsv")
    v = Counter(r["verdict"].split(":")[0] for r in m)
    inc = [r for r in m if r["verdict"] == "include"]
    fams = Counter(r["family"] for r in inc)
    genome = sum(1 for r in inc if r["source"] == "genome")
    return [
        "## 1. The alignment sets (D39)", "",
        f"**{len(m):,} census v4 rows are called to a census family; "
        f"{v['include']:,} enter an alignment** ({genome} of them genome loci) "
        f"across {len(fams)} families. Excluded, each counted in `members.tsv`:", "",
        _t(["reason", "rows"], [[k, f"{n:,}"] for k, n in v.most_common() if k != "include"]),
        "",
        "The rule was fixed before any alignment ran: a sequence enters its "
        "family's alignment only on a **high-confidence** S3a profile call "
        "(D32 margin ≥ 0.30, ≥ half the profile covered), and a genome locus "
        "only on an intact reading frame. Medium calls are where S3b found "
        "the ankyrin/LRR upper bound and where sister families sit inside the "
        "margin; a tier-1 tree is the wrong place to carry either.", ""]


def align_section() -> list[str]:
    rows = read_tsv(D / "alignments.tsv")
    st = Counter(r["status"].split(" ")[0] for r in rows)
    al = [r for r in rows if r["status"] == "aligned"]
    out = ["## 2. Family alignments — MAFFT L-INS-i + trimAl", "",
           f"**{st['aligned']} families aligned**, {st['too_few']} with fewer than "
           f"4 sequences (no alignment), {st['not_run']} not run. One method for "
           "every family — `mafft --localpair --maxiterate 1000`, then "
           "`trimal -automated1`; each run checked for row count, raggedness and "
           "unchanged residues, and recorded by SHA-256 (input, alignment, "
           "trimmed) so a rerun resumes by family.", ""]
    if al:
        secs = sum(float(r["mafft_seconds"]) for r in al)
        keep = [int(r["trim_cols"]) / int(r["aln_cols"]) for r in al]
        out += [f"L-INS-i wall time summed over families: {secs / 3600:.1f} h. "
                f"trimAl keeps a median {statistics.median(keep):.0%} of columns "
                f"(range {min(keep):.0%}–{max(keep):.0%}).", "",
                _t(["family", "sequences", "L-INS-i cols", "gap frac", "trimmed cols",
                    "trim gap", "L-INS-i s"],
                   [[r["family"], r["included"], r["aln_cols"], r["aln_gap"],
                     r["trim_cols"], r["trim_gap"], r["mafft_seconds"]]
                    for r in sorted(al, key=lambda r: -int(r["included"]))]), ""]
    few = [r["family"] for r in rows if r["status"].startswith("too_few")]
    if few:
        out += [f"Too few sequences for an alignment: {', '.join(few)}.", ""]
    return out


def module_section() -> list[str]:
    refs = read_tsv(D / "module_refs.tsv")
    spans = read_tsv(D / "module_spans.tsv")
    loo = read_tsv(D / "module_loo.tsv")
    mods = read_tsv(D / "modules.tsv")
    val = [r for r in read_tsv(D / "module_validation.tsv") if r["status"] == "ok"]
    models = read_tsv(D / "module_models.tsv")
    fam_basis = {r["family"]: r["basis"] for r in spans}
    ann = sorted(f for f, b in fam_basis.items() if b == "annotated")
    vote = sorted(f for f, b in fam_basis.items() if b.startswith("unit_hmm_vote"))
    q = Counter(r["quality"] for r in mods)
    unit_full = Counter(r["unit"] for r in mods if r["quality"] == "full")
    kf = Counter(r["k_filter"] for r in mods if r["k_filter"] and r["quality"] == "full")
    j = [float(r["jaccard"]) for r in val]
    fm = [r for r in mods if r["k_filter"] == "no" and r["quality"] == "full"]
    miss = Counter(r["k_filter_in_chain"] for r in fm)
    miss_k2p = sum(1 for r in fm if r["k_filter_in_chain"] == "yes"
                   and r["family"] == "k2p")
    errs = [int(r["vote_error"]) for r in spans if r["vote_error"]
            and fam_basis[r["family"]] == "annotated"
            and r["unit"] in ("ploop", "iglur")]
    out = [
        "## 3. Pore modules for the tier-2 units (D40)", "",
        "Five alignable superfamilies need a module rather than the full "
        "length (`Superfamily.module_rule`): P-loop, iGluR and Ca²⁺-release "
        "(`pore_loop`: the helix before each re-entrant pore loop through the "
        "helix after it), innexin clan and Hv (`tm_span`).", "",
        "**One extraction method for every sequence in a unit: "
        "`profile_projection`.** Each member is `hmmalign`ed to its own "
        "family's S3a profile and cut at that family's module span of match "
        "states. What differs between families is how the span was fixed, "
        "recorded in `module_spans.tsv`:", "",
        f"* **annotated** ({len(ann)} families) — the UniProt topology of the "
        f"family's annotated references ({sum(1 for r in refs if r['status'] == 'ok')} "
        "reference proteins), mapped to profile states, median over references;",
        f"* **unit HMM vote** ({len(vote)}: {', '.join(vote)}) — no reference in "
        "the family has an annotated pore loop, so a module HMM built from the "
        "annotated references' modules is searched over the members and those "
        "with the expected number of full hits vote. A vote boundary more than "
        "12 residues off every predicted helix moves outward to the bracketing "
        "helix: this fired once, for ITPR, whose RyR-seeded vote began in the "
        "luminal loop and would have dropped TM5.", "",
        _t(["unit", "seed modules", "seed families", "model states"],
           [[r["unit"], r["seeds"], r["families"], r["model_length"]] for r in models]),
        "",
        "### The rejected design, measured", "",
        "A single module HMM per unit was the first extraction instrument and "
        "was measured leave-one-family-out before use (`module_loo.tsv`): "
        "held-out families whose chains then carry exactly the expected number "
        "of modules —", "",
        _t(["family", "held out", "exact / members"],
           [[r["family"], r["held_out"], f"{r['exact']}/{r['members']}"]
            for r in loo if r["status"] == "ok"]), "",
        "It reaches the Kv, Nav/Cav and vertebrate iGluR families and fails on "
        "Kir, HCN, TRPM, TRPML and the whole innexin clan — so it is not the "
        "extractor. It survives only as the vote, where it is measured again "
        "(below).", "",
        "### Validation", "",
        f"**Extracted against annotated, on members never used as references: "
        f"{len(j)} modules, median Jaccard {statistics.median(j):.2f}, "
        f"{sum(x >= 0.8 for x in j)}/{len(j)} ≥ 0.8.** "
        f"**The vote, run with the held-out model on the annotated P-loop and "
        f"iGluR families, lands within {max(errs)} profile states of the "
        f"annotated span** (median {statistics.median(errs)}) — the region the "
        f"six vote families occupy. On the innexin clan it does not work "
        f"(80 states), and is not used there.", "",
        f"**Extraction**: {q['full']:,} full modules (≥ half the span's states "
        f"occupied), {q['partial']} partial, {q['absent']} absent. "
        + ", ".join(f"{u} {n:,}" for u, n in sorted(unit_full.items())) + ". "
        f"In the K⁺-filter families the TxGYG filter lies inside the extracted "
        f"module in **{kf['yes']:,}/{kf['yes'] + kf['no']:,}** "
        f"({kf['yes'] / (kf['yes'] + kf['no']):.1%}). Of the "
        f"{kf['no']} misses, {miss['no']} are chains with no canonical filter "
        "anywhere (degenerate or non-K⁺ filters, e.g. NaK's TVGDG) and "
        f"{miss['yes']} carry one elsewhere in the chain — "
        f"{miss_k2p} of them K2P modules, where one of the two pore domains "
        "carries a non-canonical filter and the other the canonical one.", ""]
    no = sorted({r["family"] for r in mods if r["quality"] == "no_span"})
    if no:
        out += [f"No span, no modules (D28): {', '.join(no)}.", ""]
    return out


def main() -> int:
    lines = ["# S6 — alignments and pore modules", "",
             "Rendered by `scripts/s6_report.py` from the tables in this "
             "directory (D13). Bulk FASTA, alignments, models and domtbls: "
             "`<data root>/alignments/s6/`.", "",
             "![S6](figures/alignments_modules.png)", ""]
    lines += sets_section() + align_section() + module_section()
    lines += ["## 4. What S7 and S8 take from here", "",
              "* **S7** roots each family alignment by adding outgroup sequences "
              "with `mafft --add --keeplength`, so the family columns do not move.",
              "* **S8** chooses representatives per clade × kingdom (D8) from the "
              "unit module FASTA and aligns the modules; every module carries "
              "`profile_projection` and its family's span basis in `modules.tsv`.",
              "* The profiles used are S3a's frozen library — the instrument that "
              "made the calls. The S2b R3 seed re-draw is a census revision, not "
              "an S6 step (roadmap emergent row).", ""]
    (D / "report.md").write_text("\n".join(lines))
    print(f"wrote {D / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
