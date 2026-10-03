"""S10b — renders `results/repertoire/s10b_report.md` + `s10b_summary.json`
purely from the committed tables beside it (D13).

    python3 scripts/s10b_report.py
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv  # noqa: E402

OUT = ROOT / "results" / "repertoire"


def _t(rows: list[dict], cols: list[str], heads: list[str] | None = None) -> str:
    heads = heads or cols
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return "\n".join(out)


def summary() -> dict:
    q4 = json.loads((OUT / "q4.json").read_text())
    ct = read_tsv(OUT / "contig_test.tsv")
    dap = json.loads((OUT / "daphnia_summary.json").read_text())
    dcmp = read_tsv(OUT / "daphnia_compare.tsv")
    orth = read_tsv(OUT / "ortholog_check.tsv")
    disp = read_tsv(OUT / "ortholog_disputed.tsv") if (OUT / "ortholog_disputed.tsv").exists() else []
    k2p = read_tsv(OUT / "k2p_nonmetazoan.tsv")
    kk = read_tsv(OUT / "k2p_by_kingdom.tsv")
    return {
        "q4": q4["verdict"], "q4_why": q4["why"],
        "contig_status": {k: dict(Counter(r["status"] for r in ct if r["kind"] == k))
                          for k in ("q4", "prok")},
        "contig_status_by_family": {f: dict(Counter(r["status"] for r in ct if r["family"] == f))
                                    for f in sorted({r["family"] for r in ct})},
        "daphnia": dap,
        "daphnia_s5_absences": {r["family"]: r["new_verdict"] for r in dcmp
                                if r["old_verdict"] == "absent"},
        "orthologs": {f: dict(Counter(r["status"] for r in orth if r["family"] == f))
                      for f in sorted({r["family"] for r in orth})},
        "disputed_calls": dict(Counter(f"{r['p_family'] or r['p_call']}:{r['p_confidence']}"
                                       for r in disp)),
        "runner_screen": len(read_tsv(OUT / "absent_runner_screen.tsv")),
        "k2p_nonmetazoan_s4": len(k2p),
        "k2p_nonmetazoan_by_name": dict(Counter(r["protein_name"] for r in k2p)),
        "k2p_present_orders_by_kingdom": {r["kingdom"] or "protists": int(r["present"])
                                          for r in kk},
    }


def main() -> int:
    s = summary()
    (OUT / "s10b_summary.json").write_text(json.dumps(s, indent=1))
    ct = read_tsv(OUT / "contig_test.tsv")
    dcmp = read_tsv(OUT / "daphnia_compare.tsv")
    orth = read_tsv(OUT / "ortholog_check.tsv")
    disp = read_tsv(OUT / "ortholog_disputed.tsv") if (OUT / "ortholog_disputed.tsv").exists() else []
    k2p = read_tsv(OUT / "k2p_nonmetazoan.tsv")
    dap = s["daphnia"]
    q4rows = [r for r in ct if r["kind"] == "q4"]
    prow = [r for r in ct if r["kind"] == "prok"]
    L = ["# S10b — Q4 and the absence checks", "",
         "*Rendered by `scripts/s10b_report.py` from the tables beside it (D13). "
         "Rules: D51, fixed before any check was read; Q4 under D50 (6).*", "",
         "## Headline", "",
         f"**Q4 (animal MscS): {s['q4']}** — {s['q4_why']}.", "",
         "## 1. Q4 — the embedded-contig test (D51 (1)–(2))", "",
         "Each subject's contig was searched by `blastx` against the S4 panel "
         "(own species removed); a region's lineage is its best hit's only past a "
         "10 % bitscore margin. *embedded* = ≥ 1 region with a metazoan best hit.", "",
         _t(q4rows, ["subject", "species", "order", "confidence", "contig", "contig_len",
                     "regions", "metazoan_regions", "foreign_regions", "unclear_regions",
                     "locus_best_lineage", "locus_best_species", "status"]), "",
         "*Reading*: no high-confidence animal MscS sits on a contig with an animal "
         "gene — the *Nematostella* locus and the *Geodia* (sponge) protein are on "
         "contigs whose other genes are non-animal, the *Trichoplax* protein on a "
         "1.2 kb contig with nothing else. The contigs that are animal carry "
         "medium-only calls (bdelloid rotifers, a second *Geodia* protein), which "
         "D50 (2) never counts as a presence. So the presence arm fails; the "
         "absence arm fails too, because S5 never returned `absent` for an animal "
         "MscS cell (bait alignments the profiles do not call: `partial` / `gap`).",
         "",
         "S5's animal MscS cells (the absence arm): " +
         ", ".join(f"{k} {v}" for k, v in Counter(
             json.loads((OUT / "q4.json").read_text())["animal_cells"].values()).items()) + ".",
         "", "### Prokaryotic-family calls in animal proteomes (descriptive, S21)", "",
         _t([{"family": f, **{k: v for k, v in c.items()}}
             for f, c in s["contig_status_by_family"].items() if f != "mscs"],
            ["family", "embedded", "foreign", "no_evidence", "failed"]), "",
         _t(prow, ["subject", "family", "species", "order", "contig_len", "regions",
                   "metazoan_regions", "foreign_regions", "locus_best_lineage",
                   "locus_best_species", "status"]), "",
         "## 2. *Daphnia pulex* on its current assembly (D51 (3))", "",
         f"Assembly {dap['old_assembly']} (2011, scaffolds) → **{dap['new_assembly']}** "
         f"(chosen by S4's `assembly_rank`; N50 {int(dap['new_n50']):,} bp). Matched "
         f"detection {dap['old_matched_detection']} → **{dap['new_matched_detection']}** "
         f"({dap['new_matched_cells']} matched control cells), so the new verdicts "
         f"{'replace' if dap['replaces_old'] else 'do not replace'} the old for S10's reading. "
         f"{dap['changed_cells']} cells change.", "",
         _t([r for r in dcmp if r["changed"] == "1" or r["old_verdict"] == "absent"],
            ["family", "old_verdict", "new_verdict", "new_best_confidence", "new_found_loci"]), "",
         "*Reading*: of the 2011 assembly's absences, " +
         ", ".join(f for f, v in s["daphnia_s5_absences"].items() if v == "absent") +
         " stay `absent` on a chromosome-level genome with every matched control "
         "detected; " + ", ".join(f for f, v in s["daphnia_s5_absences"].items()
                                  if v != "absent") +
         " become `partial` (a single weak bait alignment no profile scores), so they "
         "are no longer controlled absences. The two S10a contradictions in *Daphnia* "
         "(HCN, GPHR) stay genome presences (HCN now medium).", "",
         "## 3. ZAC, PACC1, CLCC1 against two ortholog databases (D51 (4))", "",
         "Every S5 `absent` cell of the three families; NCBI Gene orthologs and "
         "Ensembl Compara queried with the human gene.", "",
         _t(orth, ["family", "species", "ncbi", "ncbi_ids", "ensembl", "ensembl_ids", "status"]),
         ""]
    if disp:
        L += ["Every protein a database lists where S5 found none, scored against the "
              "S3a profiles (D32) and placed on S5's genome by `tblastn`:", "",
              _t(disp, ["family", "species", "protein", "length", "p_call", "p_family",
                        "p_confidence", "win_score", "runner", "genome_hit", "genome_pident",
                        "s5_locus"]), ""]
    scr = read_tsv(OUT / "absent_runner_screen.tsv")
    L += ["*Reading*: the seven ZAC absences no database disputes stand as "
          "measured; the two fish ZAC disputes are proteins Ensembl's gene trees "
          "call ZACN orthologues and our profiles call 5-HT3 or Cys-loop "
          "superfamily-only — a disagreement between orthology inference and "
          "profile assignment that neither resolves. The *Takifugu* PACC1 dispute "
          "is ours: the NCBI protein is a high-confidence PACC1 by our own "
          "profile and lies inside an S5 locus whose majority bait was TMEM87, so "
          "S5's clustering absorbed it and rescue (which ignores hits inside a "
          "locus) never looked. CLCC1 in ecdysozoans: neither database covers it "
          "(Ensembl 116's main site carries no fly or worm; NCBI's CLCC1 set names "
          "no invertebrate of the panel).", "",
          "### Post-hoc screen: absences under a runner-up bait", "",
          "Found by the PACC1 check, so labelled post hoc: every S5 `absent` cell "
          "whose family is the runner-up bait of a locus in that genome. "
          f"{len(scr)} loci; {sum(r['same_superfamily'] == '0' for r in scr)} across "
          "superfamilies. Within a superfamily this is ordinary bait overlap (Nav "
          "under Cav, ENaC under ASIC); across superfamilies it is the PACC1 "
          "artefact and a no-call rat locus.", "",
          _t(scr, ["species", "family", "locus", "span", "bait_family", "p_family",
                   "p_confidence", "same_superfamily"]), ""]
    L += ["## 4. What the K2P profile calls outside animals (D51 (5))", "",
          f"{len(k2p)} non-metazoan high-confidence K2P calls on the S4 panel:", "",
          _t(k2p, ["species", "accession", "protein_name", "pfam", "tm_features",
                   "win_score", "runner"]), "",
          "S4b orders with a high-confidence K2P call, by kingdom: " +
          ", ".join(f"{k} {v}" for k, v in s["k2p_present_orders_by_kingdom"].items()) + ".",
          "", "*Reading*: on the S4 panel the K2P profile's non-animal calls are the "
          "plant TPK two-pore K⁺ channels (every one carries two `PF07885` pore "
          "domains and four TM helices) and one *Monosiga* protein; yeast TOK1 is not "
          "called. The K2P LECA presence of S10a therefore rests on TPK-type (and, in "
          "S4b, fungal and protist) two-pore channels — whether those are K2P by "
          "descent is S17's tree question.", "", "## Files", "",
          "`contig_subjects.tsv`, `contig_test.tsv`, `contig_regions.tsv`, `q4.json`; "
          "`daphnia_assembly.json`, `daphnia_cells.tsv`, `daphnia_compare.tsv`, "
          "`daphnia_summary.json`; `ortholog_check.tsv`, `ortholog_disputed.tsv`; "
          "`absent_runner_screen.tsv`; `k2p_nonmetazoan.tsv`, `k2p_by_kingdom.tsv`; `s10b_summary.json`; "
          "`figures/s10b_checks.png`. Raw fetches under `<data root>/raw_api/s10b/`."]
    (OUT / "s10b_report.md").write_text("\n".join(L) + "\n")
    print(json.dumps(s, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
