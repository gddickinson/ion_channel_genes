"""S4 — render results/proteome_scope/report.md from the committed tables (D13).

Reads proteome_manifest.tsv, proteome_candidates.tsv, proteome_files.tsv and
summary.json; writes nothing else and computes no number the tables do not
already hold.

Run:  python3 scripts/s4_report.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.s3_hmm_lib import read_tsv  # noqa: E402
from scripts.s4_proteome_lib import OUT_DIR, RELEASE, RELEASE_DATE  # noqa: E402

BUSCO_FLOOR = 90.0      # reported, not applied: S5's D4 bar decides absences


def md_table(header: list[str], rows: list[list]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def fnum(x: str) -> str:
    return f"{int(x):,}" if x not in ("", None) else "—"


def render() -> str:
    man = read_tsv(OUT_DIR / "proteome_manifest.tsv")
    cands = read_tsv(OUT_DIR / "proteome_candidates.tsv")
    files = {f["upid"]: f for f in read_tsv(OUT_DIR / "proteome_files.tsv")}
    s = json.loads((OUT_DIR / "summary.json").read_text())
    prot = [m for m in man if m["status"] == "proteome"]
    other = [m for m in man if m["status"] != "proteome"]
    groups = sorted({m["group"] for m in man})
    cellular = [m for m in prot if m["group"] != "virus"]

    L = ["# S4 — Proteome scope (the denominator)", "",
         f"UniProt reference proteomes, release **{RELEASE}** (README dated "
         f"{RELEASE_DATE}) — the release census v2 was enumerated from (D31).", "",
         "## Headline", "",
         f"- **{s['species']} panel species** (`src/utils/species.py`) across "
         f"**{len(groups)} panel groups**: **{len(prot)} reference proteomes** "
         f"({len(cellular)} cellular, {len(prot) - len(cellular)} viral), "
         f"**{sum(m['status'] == 'genome_only' for m in other)} genome-only**, "
         f"**{sum(m['status'] == 'none' for m in other)} with neither**.",
         f"- **{s['entries']:,} canonical entries**, one per gene: the entry count "
         f"equals the proteome's declared gene count in "
         f"**{sum(int(f['entries']) == int(m['gene_count']) for m in prot for f in [files[m['upid']]])}"
         f"/{len(prot)}** proteomes.",
         f"- **Verified: {s['proteomes_exact']}/{s['proteomes']} exact** against the "
         f"release README (canonical entries = #(1), gene2acc rows = #(3)) **and** the "
         f"per-proteome `RELEASE.metalink` MD5; "
         f"{len(s['proteomes_reissued'])} reissued ({', '.join(s['proteomes_reissued']) or '—'}), "
         f"{len(s['proteomes_failed'])} failed.",
         f"- **Positive control:** {s['human_census_genes_in_proteome']}/"
         f"{s['human_census_genes_enumerated']} enumerated human census genes "
         "(census v2 `human_recall.tsv`) are entries of the human reference proteome.",
         f"- **Census v2 overlap:** {s['census_v2_in_panel']:,} of census v2's "
         f"{s['census_v2_records']:,} records are panel-proteome entries — the "
         "starting point S3b's profile sweep is measured against.",
         f"- Sweep DB for S3b: {s['sweep_db_seqs']:,} sequences, "
         f"{s['sweep_db_bytes'] / 1e6:.0f} MB, SHA-256 `{s['sweep_db_sha256'][:16]}…` "
         f"(`{s['sweep_db']}`); downloads {s['download_bytes'] / 1e6:.0f} MB.", ""]

    L += ["## The three rules", "",
          "1. **Candidates**: reference proteomes in the panel taxon's subtree "
          "(strains included) that the pinned release ships.",
          "2. **Selection**: most Swiss-Prot entries → BUSCO complete % → gene "
          "count → UPID (`s4_proteome_lib.select_proteome`).",
          "3. **No reference proteome** → the species' best current NCBI assembly "
          "(annotated > RefSeq > level > scaffold N50), status `genome_only`: its "
          "families can then be found only by S5's genomic sweep.", "",
          "Species with more than one candidate, and what the rule chose:", ""]
    by_sp = defaultdict(list)
    for c in cands:
        by_sp[c["species"]].append(c)
    rows = []
    for sp, cs in by_sp.items():
        if len(cs) > 1:
            for c in sorted(cs, key=lambda c: -int(c["selected"])):
                rows.append([sp if c is cs[0] or c["selected"] == "1" else "",
                             c["upid"], c["proteome_organism"], "**yes**" if c["selected"] == "1" else "",
                             fnum(c["reviewed"]), c["busco_complete_pct"] or "—", fnum(c["gene_count"])])
    L += md_table(["species", "UPID", "organism", "selected", "reviewed", "BUSCO %", "genes"], rows)
    L += [""]

    L += ["## Per group", ""]
    rows = []
    for g in groups:
        ms = [m for m in man if m["group"] == g]
        ps = [m for m in ms if m["status"] == "proteome"]
        rows.append([g, len(ms), len(ps), len(ms) - len(ps),
                     f"{sum(int(m['gene_count']) for m in ps):,}",
                     sum(int(files[m["upid"]]["in_census_v2"]) for m in ps)])
    L += md_table(["group", "species", "proteomes", "genome-only / none", "genes",
                   "census v2 records"], rows)
    L += [""]

    L += ["## Manifest", "",
          "Annotation source is carried per row (D9): a RefSeq gene set and a "
          "submitter gene set are not comparable evidence of absence.", ""]
    rows = []
    for m in man:
        f = files.get(m["upid"], {})
        asm = m["assembly_id"] + ("" if m["assembly_status"] in ("current", "") else
                                  f" ({m['assembly_status']}; now {m['current_assembly']})")
        rows.append([m["species"], m["group"], m["upid"] or f"*{m['status']}*",
                     fnum(m["gene_count"]) if m["upid"] else "—",
                     m["busco_complete_pct"] or "—", m["cpd_status"] or "—",
                     m["annotation_source"] or "—", asm, m["assembly_level"] or "—",
                     f.get("in_census_v2", "—"), f.get("verdict", "—")])
    L += md_table(["species", "group", "proteome", "genes", "BUSCO %", "CPD", "annotation",
                   "assembly", "level", "in census v2", "verified"], rows)
    L += [""]

    low = [m for m in prot if m["busco_complete_pct"] and float(m["busco_complete_pct"]) < BUSCO_FLOOR]
    prev = [m for m in man if m["assembly_status"] == "previous"]
    scaf = [m for m in man if m["assembly_level"] in ("Scaffold", "Contig")]
    L += ["## What S4 cannot decide — handed on", "",
          f"- **{len(low)} proteomes below {BUSCO_FLOOR:.0f} % BUSCO complete**: "
          + ", ".join(f"{m['species']} ({m['busco_complete_pct']})" for m in low)
          + ". An absence in one of these is weaker evidence; S5's D4 bars decide, "
          "this table only reports.",
          f"- **{len(prev)} proteomes are annotated on a superseded assembly version** "
          "(the manifest names the current one): "
          + ", ".join(m["species"] for m in prev)
          + ". S5 searches the assembly the gene set was built on, or says which it used.",
          f"- **{len(scaf)} assemblies are scaffold-level or below**: "
          + ", ".join(m["species"] for m in scaf) + ".",
          "- **Genome-only species**: "
          + "; ".join(f"{m['species']} — {m['assembly_id']} ({m['assembly_level']}, "
                      f"annotated={m['assembly_annotated']})" for m in other)
          + ". Neither has a reference proteome; both are catalogue exemplar sources.",
          "- **The viral proteomes hold no census v2 record**: the viroporin "
          "signatures are declared `SUBFAMILY` and so were never in the enumerated "
          "search space (D34). S3b's viroporin profile reaches these three proteomes; "
          "the census-space gap is an emergent row.",
          "- **The human proteome was reissued mid-release** (served file dated "
          + ", ".join(files[m["upid"]]["fasta_last_modified"] for m in prot
                      if m["species"] in s["proteomes_reissued"])
          + "): the README's counts describe the withdrawn file ("
          + ", ".join(f"{int(files[m['upid']]['entries_expected']):,} entries against "
                      f"{int(files[m['upid']]['entries']):,} served" for m in prot
                      if m["species"] in s["proteomes_reissued"])
          + " — the whole human UniProtKB set, not one per gene). The served file matches the "
          "metalink MD5 and the declared gene count, and is used (D35).", ""]
    return "\n".join(L) + "\n"


def main() -> int:
    (OUT_DIR / "report.md").write_text(render())
    print(f"[s4_report] wrote {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
