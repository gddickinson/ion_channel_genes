"""s15_report.py — render results/method_contribution/report.md (D13).

    python3 scripts/s15_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv  # noqa: E402
from s15_contribution import HIGH, LOW, OUT, ROOT  # noqa: E402


def pct(x) -> str:
    return f"{float(x):.0%}" if x != "" else "—"


def main() -> int:
    sf = read_tsv(OUT / "curve_by_superfamily.tsv")
    fam = read_tsv(OUT / "curve_by_family.tsv")
    grp = read_tsv(OUT / "curve_by_group.tsv")
    hum = read_tsv(OUT / "human_genes.tsv")
    jh = read_tsv(ROOT / "results" / "panel_sweep" / "jackhmmer_completeness.tsv")
    n = sum(int(r["n"]) for r in sf)
    tot = {k: sum(int(r[f"n_{k}"]) for r in sf)
           for k in ("domain_call", "domain_enumeration", "profile", "genome")}
    prot = sum(int(r["proteome_rows"]) for r in sf)
    q_enum = Counter(r["q3_enumeration"] for r in sf)
    q_call = Counter(r["q3_call"] for r in sf)
    h = {k: sum(int(r[k]) for r in hum) for k in ("enumerated", "domain_call", "profile_call")}

    L = ["# S15 — what each method contributes (Q3)", "",
         "Rendered by `scripts/s15_report.py` from the tables in this directory (D13), "
         "which `scripts/s15_contribution.py` derives from census v4 "
         "(`<data root>/genomes/s5/census_v4.tsv.gz`) and the S2/S3b human-recall tables.", "",
         "![S15](figures/method_contribution.png)", "",
         "## The answer", "",
         f"**Domain search finds the channels; it cannot name most of them.** Of the "
         f"final census's {n:,} high-confidence census-family members (panel frame, "
         f"below), **{(tot['domain_call'] + tot['domain_enumeration']) / n:.1%} carry an "
         f"enumerated pore signature**, but S2's domain rules call only "
         f"**{tot['domain_call'] / n:.1%}** of them to the family the profiles do. "
         f"Profiles add {tot['profile']} proteome members domain search never enumerated "
         f"({tot['profile'] / n:.1%}); genomes add {tot['genome']} loci the proteomes "
         f"lack ({tot['genome'] / n:.1%}).", "",
         f"On the independent human frame (the {len(hum)} curated census genes): domain "
         f"search enumerates **{h['enumerated']}**, calls **{h['domain_call']}** to the "
         f"right family, and the profiles call **{h['profile_call']}**.", "",
         f"**Per superfamily (thresholds fixed in `s15_contribution.py` before any row "
         f"was read: ≥ {HIGH:.0%} / ≥ {LOW:.0%})** — enumeration: "
         + ", ".join(f"{k} {v}" for k, v in q_enum.most_common())
         + "; the family *call*: " + ", ".join(f"{k} {v}" for k, v in q_call.most_common())
         + ". Enumeration is not the bottleneck anywhere except "
         + ", ".join(f"{r['superfamily']} ({pct(r['domain_enumeration_frac'])})"
                     for r in sf if r["q3_enumeration"] != "domain_search")
         + ". The call is: in "
         + ", ".join(r["superfamily"] for r in sf if float(r["domain_call_frac"]) == 0)
         + " domain rules call **no** member to its family, because the families inside "
           "each share one architecture (D25) — S2 stops at the superfamily and only the "
           "profiles separate them.", "",
         "## Two frames", "",
         "* **Human frame — independent truth.** The catalogue's 320 human census genes, "
         "curated from the literature, never from a census method.",
         "* **Panel frame — the final census.** Every census v4 row the S3a profiles call "
         "to a census family at high confidence (D32), genome loci only with an intact "
         "frame (D37); contested families kept. The frame is defined by the best "
         "instrument, so the profile step reaches 100 % of proteome rows by construction: "
         "the curve measures what domain search misses *relative to it*. The profile "
         "step's own completeness is checked independently by jackhmmer (last column, "
         "S3b clean runs only, D10).", "",
         "## The curve, per superfamily", "",
         "| superfamily | members | domain call | + domain enumeration | + profile | + genome "
         "| Q3 (enumeration) | Q3 (call) | domain-only calls | human: enumerated / domain call "
         "/ profile | jackhmmer recall (clean runs) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in sf:
        clean = [j for j in jh if j["superfamily"] == r["superfamily"] and j["verdict"] == "clean"]
        jr = (f"{sum(int(j['recovered']) for j in clean)}/"
              f"{sum(int(j['profile_calls']) for j in clean)}" if clean else "no clean run")
        L.append(f"| {r['superfamily']} | {r['n']} | {pct(r['cum_domain_call'])} | "
                 f"{pct(r['cum_domain_enumeration'])} | {pct(r['cum_profile'])} | "
                 f"{pct(r['cum_genome'])} | {r['q3_enumeration']} | {r['q3_call']} | "
                 f"{r['domain_only_calls']} | {r['human_enumerated']} / "
                 f"{r['human_domain_call']} / {r['human_profile_call']} of "
                 f"{r['human_genes']} | {jr} |")
    L += ["", "*Domain-only calls*: S2 family calls on census v4 rows the profiles do not "
          "call (census v3 basis `s2_only`) — domain search's unique contribution, outside "
          "the frame because no profile confirms them.", "",
          "## Where domain search misses, by lineage", "",
          "Domain enumeration recall by panel group (proteome rows), superfamilies with a "
          "group below 95 %:", "",
          "| superfamily | group | members | enumerated | family call |", "|---|---|---|---|---|"]
    for r in grp:
        if float(r["domain_enumeration_frac"]) < HIGH and int(r["proteome_rows"]) >= 3:
            L.append(f"| {r['superfamily']} | {r['group']} | {r['proteome_rows']} | "
                     f"{pct(r['domain_enumeration_frac'])} | {pct(r['domain_call_frac'])} |")
    gen = sorted(fam, key=lambda r: -int(r["n_genome"]))[:8]
    L += ["", "## What the genome sweep adds", "",
          "Genome loci in the frame come from the two genome-only species (*Cornu*, "
          "*Torpedo*) and the proteome misses S5b confirmed. Largest by family: "
          + ", ".join(f"{r['family']} {r['n_genome']}" for r in gen if int(r["n_genome"])) + ".",
          "", "## Not done here", "",
          "* The two method designs the S3b and S5a emergent rows proposed for S15 — a "
          "superfamily-seeded jackhmmer and a six-frame profile scan of genomes — are new "
          "instruments, not measurements of the existing ones; they stay open rows.",
          "* This measures the census's own methods against each other and against the "
          "human gene list. What the *annotation databases* miss or mislabel per locus is "
          "S16.", ""]
    (OUT / "report.md").write_text("\n".join(L))
    (OUT / "summary.json").write_text(json.dumps(
        {"frame_members": n, **{f"n_{k}": v for k, v in tot.items()},
         "proteome_members": prot, "human": {"genes": len(hum), **h},
         "q3_enumeration": dict(q_enum), "q3_call": dict(q_call)}, indent=1))
    print(f"→ {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
