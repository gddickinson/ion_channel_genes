"""s5_census_v4.py — census v4: census v3 plus what only the genomes show.

    python3 scripts/s5_census_v4.py        # after s5_verdict.py

Census v4 is census v3 (the panel proteomes, S3b) **unchanged**, plus one
row per genome locus that the proteomes cannot supply. The rule, fixed
before any cell was read:

* a genome-only species (*Cornu*, *Torpedo*: no proteome) — every locus the
  S3a profiles call to a census family (high or medium, D32), basis
  `genome_only`. Nothing to merge with: there is no proteome call.
* a proteome species whose cell verdict is `genome_found` — only the loci
  called at **high confidence on an intact reading frame** (D37 (4)), basis
  `genome_found`: the proteome has no call for that family, the genome has a
  gene. A `genome_weak` cell adds nothing (retrocopies, chained models).

Genome rows are **never merged into a proteome call** and never change one.
A row is a *locus*, not a gene: copy number needs intact, exon-confirmed
loci (S11). Bulk: `<data root>/genomes/s5/census_v4.tsv.gz`. Committed:
`census_v4_status.tsv`, `family_by_species_v4.tsv`, `census_v4.json`.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, sha256, write_tsv  # noqa: E402
from s3b_lib import s3b_dir  # noqa: E402
from s5_ledger import CALLED, load_loci  # noqa: E402
from s5_lib import OUT_DIR, s5_dir  # noqa: E402

GENOME_FIELDS = ["assembly", "locus", "contig", "strand", "call_start", "call_end",
                 "call_identity", "call_intact", "call_bait"]
ADD_VERDICTS = ("genome_present", "genome_found")


def genome_rows(cells: list[dict], runs: dict) -> list[dict]:
    """Rows for the loci census v4 adds, by the rule in the docstring."""
    want = {(c["species"], c["family"]): c["verdict"] for c in cells
            if c["verdict"] in ADD_VERDICTS}
    out = []
    for sp in sorted({s for s, _ in want}):
        run = runs[sp]
        for L in load_loci(run["assembly"]):
            v = want.get((sp, L["p_family"]))
            if L["p_call"] != "family" or v is None or L["p_confidence"] not in CALLED:
                continue
            if v == "genome_found" and not (L["p_confidence"] == "high"
                                            and L.get("call_intact") == "1"):
                continue
            out.append({"target": f"{run['assembly']}:{L['locus']}", "accession": "",
                        "species": sp, "group": run["group"],
                        "length": len(L.get("call_translation", "")), "in_v2": 0,
                        "p_call": L["p_call"], "p_family": L["p_family"],
                        "p_superfamily": L["p_superfamily"],
                        "p_confidence": L["p_confidence"], "win_score": L["win_score"],
                        "win_coverage": L["win_coverage"], "runner": L["runner"],
                        "rel_margin": L["rel_margin"],
                        "v3_family": L["p_family"], "v3_superfamily": L["p_superfamily"],
                        "v3_status": "genome_locus",
                        "v3_basis": "genome_only" if v == "genome_present" else "genome_found",
                        "source": "genome", "assembly": run["assembly"],
                        **{k: L.get(k, "") for k in GENOME_FIELDS[1:]}})
    return out


def found_check(cells: list[dict], v3: list[dict], runs: dict) -> list[dict]:
    """Per `genome_found` cell: could the proteome hold the gene under a
    sister family's call? Lists the genome's loci of the family and the
    proteome entries whose runner-up profile is that family, with margins.
    A report on the verdict, never an input to it (added after S5b's cells
    were read, and labelled so)."""
    out = []
    for c in (c for c in cells if c["verdict"] == "genome_found"):
        near = [r for r in v3 if r["species"] == c["species"]
                and r["runner"] == c["family"] and r["p_call"] == "family"]
        loci = [L for L in load_loci(runs[c["species"]]["assembly"])
                if c["family"] in (L["p_family"], L["bait_family"])]
        out.append({"species": c["species"], "family": c["family"],
                    "genome_loci": len(loci),
                    "genome_loci_called": sum(L["p_family"] == c["family"] for L in loci),
                    "proteome_near": len(near),
                    "near_min_margin": min((float(r["rel_margin"]) for r in near), default=""),
                    "near_calls": ";".join(sorted(f"{r['accession']}:{r['p_family']}:{r['rel_margin']}"
                                                  for r in near))})
    return out


def main() -> int:
    runs = {r["species"]: r for r in read_tsv(OUT_DIR / "genome_runs.tsv")
            if not r.get("note", "").startswith("FAILED")}
    cells = read_tsv(OUT_DIR / "cells.tsv")
    v3_path = s3b_dir() / "census_v3.tsv.gz"
    v3 = read_tsv(v3_path)
    fields = list(v3[0].keys()) + ["source"] + [f for f in GENOME_FIELDS
                                                if f not in v3[0]]
    for r in v3:
        r["source"] = "proteome"
    add = genome_rows(cells, runs)
    out = s5_dir() / "census_v4.tsv.gz"
    write_tsv(out, fields, v3 + add)

    status = Counter((r["source"], r["v3_basis"], r["v3_status"]) for r in v3 + add)
    write_tsv(OUT_DIR / "census_v4_status.tsv", ["source", "v3_basis", "v3_status", "records"],
              [{"source": s, "v3_basis": b, "v3_status": st, "records": n}
               for (s, b, st), n in sorted(status.items(), key=lambda x: -x[1])])

    loci_n = Counter((r["species"], r["v3_family"]) for r in add)
    matrix = [{"family": c["family"], "superfamily": c["superfamily"],
               "species": c["species"], "group": c["group"],
               "proteome_records": c["proteome_records"],
               "proteome_high": c["proteome_high"],
               "genome_loci_added": loci_n[(c["species"], c["family"])],
               "verdict": c["verdict"]} for c in cells]
    write_tsv(OUT_DIR / "family_by_species_v4.tsv", list(matrix[0].keys()), matrix)

    chk = found_check(cells, v3, runs)
    write_tsv(OUT_DIR / "genome_found_check.tsv",
              list(chk[0].keys()) if chk else ["species", "family"], chk)

    summary = {"census_v3_rows": len(v3), "genome_rows": len(add),
               "genome_only_rows": sum(r["v3_basis"] == "genome_only" for r in add),
               "genome_found_rows": sum(r["v3_basis"] == "genome_found" for r in add),
               "genome_cells_added": len(loci_n),
               "census_v3_sha256": sha256(v3_path), "census_v4_sha256": sha256(out),
               "census_v4": str(out)}
    (OUT_DIR / "census_v4.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
