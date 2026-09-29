"""S2 — render `results/census_v2/report.md` purely from the committed tables (D13).

    python3 scripts/s2_report.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s0_lib import read_tsv          # noqa: E402
from scripts.s2_lib import OUT_DIR           # noqa: E402
from src.catalogue import CATALOGUE          # noqa: E402

CHANNELISH = ("channel", "channel_contested")

# Revisions of census v2, each rendered from its own transitions table.
REVISIONS = [
    ("s2b", "Revision r2 — hazard rules H2, H4, H13 as positive tests (S2b)",
     "S3a's calibration and its check against the IP3R project's census found "
     "three hazard rules that made a family call from what a protein *lacks* — "
     "which this project's conventions forbid. Each was rewritten so a family "
     "is called only on a domain it carries and its rivals do not; where no "
     "such domain exists (ITPR against RyR; AChBP against a receptor fragment) "
     "the architecture tier stops at the superfamily and a sequence-level tier "
     "makes the call (D33)."),
    ("s2c", "Revision r3 — hazard rules H11, H12 as positive tests (S2c)",
     "The two remaining absence rules. KCTD is now called on a KCTD C-terminal "
     "domain, and the T1 domain alone is superfamily-only; CFTR is still "
     "called on its R domain, and the SUR rule is gone — no domain identifies "
     "SUR, and none of the 980 calls the absence rule made was one — 918 "
     "bacterial ABC transporters and 62 eukaryotic fused gene models."),
]


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f} %" if b else "—"


def table(rows: list[dict], cols: list[str], limit: int | None = None) -> str:
    rows = rows[:limit] if limit else rows
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def main() -> int:
    S = json.loads((OUT_DIR / "summary.json").read_text())
    shards = read_tsv(OUT_DIR / "shard_counts.tsv")
    sigs = read_tsv(OUT_DIR / "signature_counts.tsv")
    fams = read_tsv(OUT_DIR / "census_families.tsv")
    sfo = read_tsv(OUT_DIR / "census_superfamily_only.tsv")
    tiers = read_tsv(OUT_DIR / "tier_attribution.tsv")
    unas = read_tsv(OUT_DIR / "unassigned_signatures.tsv")
    filt = read_tsv(OUT_DIR / "four_repeat_filters.tsv")
    hum = read_tsv(OUT_DIR / "human_recall.tsv")
    s1 = read_tsv(OUT_DIR / "s1_panel.tsv")
    man = read_tsv(OUT_DIR / "manifest.tsv")
    part = read_tsv(OUT_DIR / "partial_architectures.tsv")

    n = S["records"]
    s1_pos = [r for r in s1 if CATALOGUE[r["expected_family"]].status.value in CHANNELISH]
    s1_neg = [r for r in s1 if r not in s1_pos]
    s1_neg_in = [r for r in s1_neg if r["s2_status"] != "not_enumerated"]
    s1_neg_leak = [r for r in s1_neg_in if r["s2_status"] in CHANNELISH]
    union = next(r for r in shards if r["shard"] == "UNION")
    bad_shards = [r for r in shards if r["status"] != "ok"]
    bad_sigs = [r for r in sigs if r["status"] != "ok"]
    ip_diff = sorted(sigs, key=lambda r: -abs(int(r["uniprot_minus_interpro"] or 0)))
    chan = [r for r in fams if r["catalogue_status"] == "channel"]
    nonchan = [r for r in fams if r["catalogue_status"] != "channel"]
    hum_abs = [r for r in hum if r["s2_status"] == "absent"]
    hum_wrong = [r for r in hum if r["s2_status"] != "absent" and r["correct"] == "no"]

    L = [f"# S2 — census v2: every UniProtKB protein carrying a pore signature",
         "", f"*Rendered from the tables in this directory by `scripts/s2_report.py`. "
         f"UniProtKB release {S['uniprot_release']}. Revision "
         f"{S.get('revision', 'r1')} — earlier revisions are in git history.*", "",
         "## Headline", "",
         f"- **{n:,} records** carry at least one of the catalogue's "
         f"{len(sigs)} pore signatures (union count {int(union['uniprot_count']):,}; "
         f"fetched {int(union['fetched']):,}).",
         f"- **{S['family_calls']:,} ({pct(S['family_calls'], n)}) get a family call**, "
         f"{S['family_calls_channel']:,} of them to a family catalogued as a channel; "
         f"{S['superfamily_only']:,} ({pct(S['superfamily_only'], n)}) reach a superfamily only; "
         f"{S['unassigned']:,} ({pct(S['unassigned'], n)}) are `unassigned` and stay in the census "
         "with that label.",
         f"- **Every family call was made without the reference tier**: "
         f"{S['family_calls_by_architecture_or_hazard']:,} by architecture or hazard rules and "
         f"{S['family_calls_by_motif']:,} by the selectivity-filter motif "
         f"({S['four_repeat_projected']:,} four-repeat records projected onto Nav1.5).",
         f"- **Human census genes: {S['human_found']}/{S['human_genes']} enumerated, "
         f"{S['human_correct']}/{S['human_genes']} called to the right family** "
         "(reviewed human entries; the gene symbol is used here to *score*, never to classify — H15).",
         f"- **S1 panel: {S['s1_enumerated']}/{S['s1_panel']} enumerated; "
         f"{S['s1_correct_in_s2']} correct without the reference tier vs "
         f"{S['s1_correct_in_s1']} with it in S1.** Panel members from channel families "
         f"enumerated: **{sum(r['s2_status'] != 'not_enumerated' for r in s1_pos)}/{len(s1_pos)}**. "
         f"Non-channel members: {len(s1_neg) - len(s1_neg_in)}/{len(s1_neg)} carry no pore "
         f"signature and are excluded at enumeration; of the {len(s1_neg_in)} that are "
         f"enumerated, **{len(s1_neg_leak)} are called to a channel family** "
         f"({', '.join(r['gene'] + ' → ' + (r['s2_family'] or r['s2_status']) for r in s1_neg_in)}).",
         "",
         "## Completeness", "",
         "The union query was cut into taxonomic shards that partition it; the shard "
         "counts must sum to the union count, and each shard's fetched records must "
         "equal UniProt's own count for it.", "",
         table(shards, ["shard", "uniprot_count", "fetched", "status"]), "",
         f"Shards or union short: **{len(bad_shards)}**. Per-signature check "
         f"(records fetched carrying the signature == UniProt's count for it alone): "
         f"**{len(sigs) - len(bad_sigs)}/{len(sigs)} ok**.", ""]
    if bad_sigs:
        L += [table(bad_sigs, ["pfam", "uniprot_count", "fetched", "status"]), ""]
    L += ["UniProt and InterPro run on different release cycles, so their counts for "
          "the same Pfam entry differ. The census is UniProt's; the largest differences "
          "are listed so the gap is visible rather than assumed away.", "",
          table(ip_diff, ["pfam", "uniprot_count", "interpro_count",
                          "uniprot_minus_interpro"], 10), "",
          "## Calls by family — channels", "",
          table(chan, ["family", "superfamily", "records", "reviewed", "human",
                       "eukaryota", "bacteria", "archaea", "viruses", "gold",
                       "silver", "bronze"]), "",
          "## Calls by family — catalogued non-channels (controls, auxiliaries, "
          "transporters, out of scope)", "",
          "These records carry a pore signature and were positively called to a family "
          "the catalogue holds *in order to exclude it*. They are in the census so the "
          "exclusion is visible.", "",
          table(nonchan, ["family", "superfamily", "catalogue_status", "records",
                          "reviewed", "eukaryota", "bacteria", "archaea", "viruses"]), "",
          "## Superfamily only", "",
          "`ambiguous` counts the records whose architecture named two or more "
          "candidate families and could not choose (D7: a tie is not a call); the "
          "rest matched superfamily-level evidence only.", "",
          table(sfo, ["superfamily", "records", "ambiguous", "reviewed", "eukaryota", "bacteria",
                      "archaea", "viruses"]), "",
          "## Which tier decided", "",
          table(tiers, ["decisive_tier", "status", "records"]), "",
          "## Four-repeat selectivity filters", "",
          table(filt, ["filter_string", "family", "records"], 25), "",
          "## What brought the unassigned records in", "",
          "`records` counts every unassigned record carrying the signature; "
          "`solo_records` those for which it is the only pore signature.", "",
          table(unas, ["pfam", "records", "solo_records", "eukaryota", "bacteria",
                       "archaea", "viruses"], 30), "",
          "## Unassigned records carrying part of one family's architecture", "",
          "A derived family rule requires *every* `FAMILY`-level signature the "
          "catalogue declares for the family (`rules.derived_rules`). A record carrying "
          "some but not all of one family's signatures therefore goes unassigned. This "
          "table counts them: a measure of how far the declared architectures "
          "are from what the family's members actually carry, and the input to "
          "any rule revision — which must be re-benchmarked on S1 before it is "
          "adopted, not tuned here.", "",
          table(part, ["family", "declared_family_signatures", "unassigned_records"], 30), "",
          "## Human census genes", "",
          f"Not enumerated ({len(hum_abs)}): the reviewed human entry carries none of "
          "the pore signatures, or has no reviewed entry under that symbol.", "",
          table(hum_abs, ["gene", "expected_family"]), "",
          f"Enumerated but not called to the expected family ({len(hum_wrong)}):", "",
          table(hum_wrong, ["gene", "expected_family", "accession", "s2_family",
                            "s2_status", "confidence"]), "",
          "## S1 panel: S1 call vs S2 call", "",
          table([r for r in s1 if r["s1_correct"] != r["s2_correct"]],
                ["accession", "gene", "expected_family", "s1_family", "s2_family",
                 "s2_status", "s1_correct", "s2_correct"]), "",
          "(Only rows where S1 and S2 disagree on correctness are shown; all rows are "
          "in `s1_panel.tsv`.)", "",
          "## Bulk files (data root)", "",
          table(man, ["file", "bytes", "records", "sha256"]), ""]
    for tag, title, intro in REVISIONS:
        tt = OUT_DIR / f"{tag}_transitions.tsv"
        if not tt.exists():
            continue
        R = json.loads((OUT_DIR / f"{tag}_revision.json").read_text())
        L += ["", f"## {title}", "", intro + " "
              f"**{R['rechecked']:,} records** carry an accession the rewritten "
              f"rules consult ({', '.join('`' + x + '`' for x in R['accessions'])}) "
              f"and were re-classified; **{R['changed']:,} calls changed**, and "
              f"{R['changed_outside_recheck']} changed outside that set (checked, "
              f"not assumed). The previous call files are archived on the data "
              f"root under `calls_{R['base']}/`.", "",
              table(read_tsv(tt), ["before", "after", "hazards_after", "records"], 12), ""]
    (OUT_DIR / "report.md").write_text("\n".join(L))
    print(f"wrote {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
