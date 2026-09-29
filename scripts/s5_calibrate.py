"""s5_calibrate.py — gene span and widest intron, measured from annotation.

    python3 scripts/s5_calibrate.py

For every locus the S3a profiles call to a census family at high confidence,
the annotated gene it overlaps (same strand, most shared bases) is read from
the assembly's own NCBI annotation (or its RefSeq twin's). Two thresholds
come out, neither from miniprot, whose `-G` shapes what it reports:

* `spans.tsv` — per family, the annotated gene spans; **D4's bar is the
  median** (the parent's rule: half the family's genes fit in a scaffold of
  that length). Pooled over the annotated genomes, with the per-group
  median and count kept beside it, because a nematode gene is not a mouse
  gene; `s5_verdict.d4_bar` chooses between them (D38).
* `introns.tsv` — per genome, the widest annotated intron at a called locus
  against the `-G` the sweep used. An annotated intron wider than `-G` is a
  gene miniprot can only report split.
"""

from __future__ import annotations

import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from src.catalogue import registry  # noqa: E402
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s5_annotation import GeneIndex, ensure_annotation, read_genes, seqid_map  # noqa: E402
from s5_ledger import load_loci  # noqa: E402
from s5_lib import OUT_DIR  # noqa: E402

LOCUS_FIELDS = ["species", "group", "family", "locus", "gene_id", "gene_name",
                "gene_span", "annotated_max_intron", "n_exons",
                "miniprot_max_intron", "annotation_accession"]
SPAN_FIELDS = ["family", "n_genes", "n_species", "median_span", "min_span",
               "max_span", "bar_bp", "by_group", "by_group_n"]
INTRON_FIELDS = ["species", "group", "annotation_accession", "called_loci",
                 "annotated_loci", "widest_annotated_intron", "widest_gene",
                 "max_intron_used", "over_G"]


def main() -> int:
    census = {f.key for f in registry.census_families()}
    rows, intron_rows = [], []
    for run in read_tsv(OUT_DIR / "genome_runs.tsv"):
        if run.get("note", "").startswith("FAILED"):
            continue
        gff, src = ensure_annotation(run["assembly"])
        called = [L for L in load_loci(run["assembly"])
                  if L["p_call"] == "family" and L["p_confidence"] == "high"
                  and L["p_family"] in census]
        rec = {"species": run["species"], "group": run["group"],
               "annotation_accession": src, "called_loci": len(called),
               "max_intron_used": run["max_intron"]}
        if gff is None:
            intron_rows.append({**rec, "annotated_loci": 0})
            print(f"{run['species']}: no annotation", flush=True)
            continue
        idx = GeneIndex(read_genes(gff, seqid_map(run["assembly"])))
        mine = []
        for L in called:
            cds = [tuple(map(int, x.split("-"))) for x in L["call_cds"].split(",") if x]
            g = idx.best_overlap(L["contig"], int(L["call_start"]),
                                 int(L["call_end"]), L["strand"], cds)
            if g is None:
                continue
            mine.append({"species": run["species"], "group": run["group"],
                         "family": L["p_family"], "locus": L["locus"],
                         "gene_id": g["gene_id"], "gene_name": g["name"],
                         "gene_span": g["end"] - g["start"] + 1,
                         "annotated_max_intron": g["max_intron"],
                         "n_exons": g["n_exons"],
                         "miniprot_max_intron": L["call_max_intron"],
                         "annotation_accession": src})
        rows += mine
        widest = max(mine, key=lambda r: r["annotated_max_intron"], default=None)
        rec.update(annotated_loci=len(mine),
                   widest_annotated_intron=widest["annotated_max_intron"] if widest else "",
                   widest_gene=f"{widest['gene_name']} ({widest['family']})" if widest else "",
                   over_G=sum(r["annotated_max_intron"] > int(run["max_intron"])
                              for r in mine))
        intron_rows.append(rec)
        print(f"{run['species']}: {len(mine)}/{len(called)} loci annotated; "
              f"widest intron {rec['widest_annotated_intron']}", flush=True)

    # one gene per (species, gene_id): several loci can hit one gene
    uniq = {(r["species"], r["gene_id"]): r for r in rows}.values()
    by_f = defaultdict(list)
    for r in uniq:
        by_f[r["family"]].append(r)
    span_rows = []
    for f in sorted(by_f):
        s = [r["gene_span"] for r in by_f[f]]
        grp = defaultdict(list)
        for r in by_f[f]:
            grp[r["group"]].append(r["gene_span"])
        span_rows.append({
            "family": f, "n_genes": len(s),
            "n_species": len({r["species"] for r in by_f[f]}),
            "median_span": int(statistics.median(s)), "min_span": min(s),
            "max_span": max(s), "bar_bp": int(statistics.median(s)),
            "by_group": ";".join(f"{g}:{int(statistics.median(v))}"
                                 for g, v in sorted(grp.items())),
            "by_group_n": ";".join(f"{g}:{len(v)}" for g, v in sorted(grp.items()))})
    write_tsv(OUT_DIR / "annotated_loci.tsv", LOCUS_FIELDS, rows)
    write_tsv(OUT_DIR / "spans.tsv", SPAN_FIELDS, span_rows)
    write_tsv(OUT_DIR / "introns.tsv", INTRON_FIELDS, intron_rows)
    print(f"{len(span_rows)} families with a measured span")
    return 0


if __name__ == "__main__":
    sys.exit(main())
