"""S4 — the proteome scope: one declared row per panel species.

Stages (run in order; each is resumable and reruns offline from the archive):

  manifest   candidates → selection → NCBI assembly for every row
             → results/proteome_scope/{proteome_candidates,proteome_manifest}.tsv
  download   canonical FASTA + gene2acc per selected proteome
             → <data root>/proteomes/s4/
  verify     every file counted against the release README (#1 canonical
             entries, #3 gene2acc rows) — the S2-style exactness check —
             SHA-256, gene groups, overlap with census v2; concatenates the
             sweep DB for S3b → proteome_files.tsv, summary.json

Run:  python3 scripts/s4_proteomes.py all
      python3 scripts/s4_proteomes.py manifest --refresh   # re-query the APIs
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import s4_proteome_lib as L  # noqa: E402
from scripts.s0_lib import live_progress  # noqa: E402
from scripts.s3_hmm_lib import iter_census_v2, read_tsv, sha256, write_tsv  # noqa: E402
from src.utils.species import all_species  # noqa: E402

CAND_FIELDS = ["species", "panel_taxid", "upid", "proteome_taxid", "proteome_organism",
               "in_release", "selected", "reviewed", "busco_complete_pct", "gene_count",
               "protein_count", "modified"]
MANIFEST_FIELDS = [
    "species", "common", "group", "panel_taxid", "status", "n_candidates",
    "n_candidates_not_in_release", "upid", "proteome_taxid", "proteome_organism",
    "superregnum", "modified", "gene_count", "protein_count", "reviewed",
    "readme_canonical", "readme_additional", "readme_gene2acc",
    "busco_complete_pct", "busco_fragmented", "busco_missing", "busco_lineage",
    "cpd_status", "annotation_source", "assembly_id", "assembly_source",
    "assembly_level", "refseq_category", "assembly_total_length",
    "scaffold_n50", "contig_n50", "assembly_annotated", "assembly_status",
    "current_assembly", "assembly_note"]
FILE_FIELDS = ["species", "upid", "fasta", "fasta_bytes", "fasta_sha256",
               "entries", "entries_expected", "gene2acc_rows", "gene2acc_expected",
               "gene_groups", "unknown_gene_rows", "in_census_v2",
               "in_census_v2_family", "md5_ok", "fasta_last_modified", "reissued", "verdict"]


def log(msg: str) -> None:
    print(f"[s4] {msg}", flush=True)


def progress(steps_done: dict) -> None:
    live_progress(L.LIVE, "S4", [(k, v) for k, v in steps_done.items()])


# ------------------------------------------------------------------ manifest
def assembly_columns(acc: str, refresh: bool) -> dict:
    rec = L.fetch_assembly(acc, refresh) if acc else None
    if not rec:
        return {"assembly_note": "no assembly accession" if not acc
                else "assembly not returned by NCBI Datasets"}
    a = L.flatten_assembly(rec)
    return {"assembly_level": a["level"], "refseq_category": a["refseq_category"],
            "assembly_total_length": a["total_length"], "scaffold_n50": a["scaffold_n50"],
            "contig_n50": a["contig_n50"], "assembly_annotated": int(a["annotated"]),
            "assembly_status": a["status"], "current_assembly": a["current_accession"]}


def genome_only_row(sp, refresh: bool) -> dict:
    """No reference proteome: the species' best current NCBI assembly, or none."""
    asms = [L.flatten_assembly(r) for r in L.fetch_taxon_assemblies(sp.taxon_id, refresh)]
    if not asms:
        return {"status": "none", "assembly_note": "no reference proteome, no assembly"}
    best = max(asms, key=L.assembly_rank)
    return {"status": "genome_only", "assembly_id": best["accession"],
            "assembly_source": "NCBI", "assembly_level": best["level"],
            "refseq_category": best["refseq_category"],
            "assembly_total_length": best["total_length"],
            "scaffold_n50": best["scaffold_n50"], "contig_n50": best["contig_n50"],
            "assembly_annotated": int(best["annotated"]),
            "assembly_status": best["status"], "current_assembly": best["current_accession"],
            "assembly_note": f"best of {len(asms)} current assemblies (assembly_rank)"}


def stage_manifest(refresh: bool) -> list[dict]:
    release, readme = L.parse_readme(L.fetch_release_readme(refresh))
    if L.RELEASE not in release:
        raise SystemExit(f"release README says {release!r}, expected {L.RELEASE} "
                         "(census v2's release); refusing to mix releases")
    log(f"{release}: {len(readme):,} reference proteomes in the release")
    cand_rows, rows = [], []
    for sp in all_species():
        cands = [L.flatten_proteome(r) for r in L.fetch_reference_proteomes(sp.taxon_id, refresh)]
        in_rel = [c for c in cands if c["upid"] in readme]
        chosen = L.select_proteome(in_rel)
        for c in cands:
            cand_rows.append({**c, "species": sp.scientific, "panel_taxid": sp.taxon_id,
                              "in_release": int(c["upid"] in readme),
                              "selected": int(chosen is not None and c["upid"] == chosen["upid"])})
        row = {"species": sp.scientific, "common": sp.common, "group": sp.group,
               "panel_taxid": sp.taxon_id, "n_candidates": len(cands),
               "n_candidates_not_in_release": len(cands) - len(in_rel)}
        if chosen:
            rd = readme[chosen["upid"]]
            if rd["taxid"] != chosen["proteome_taxid"]:
                raise SystemExit(f"{chosen['upid']}: README taxid {rd['taxid']} "
                                 f"!= API {chosen['proteome_taxid']}")
            row.update(chosen)
            row.update({"status": "proteome", "readme_canonical": rd["n_canonical"],
                        "readme_additional": rd["n_additional"],
                        "readme_gene2acc": rd["n_gene2acc"]})
            row.update(assembly_columns(chosen["assembly_id"], refresh))
        else:
            row.update(genome_only_row(sp, refresh))
        rows.append(row)
        log(f"  {sp.scientific[:34]:34s} {row['status']:11s} {row.get('upid', ''):12s} "
            f"{len(cands)} cand. {row.get('assembly_id', '')}")
    write_tsv(L.OUT_DIR / "proteome_candidates.tsv", CAND_FIELDS, cand_rows)
    write_tsv(L.OUT_DIR / "proteome_manifest.tsv", MANIFEST_FIELDS, rows)
    return rows


# ------------------------------------------------------------------ download
def stage_download() -> int:
    rows = [r for r in read_tsv(L.OUT_DIR / "proteome_manifest.tsv") if r["status"] == "proteome"]
    jobs = []
    for r in rows:
        for url, dest in L.proteome_files(r["upid"], int(r["proteome_taxid"]),
                                          r["superregnum"]).values():
            jobs.append((r["upid"], url, dest))
    with ThreadPoolExecutor(6) as pool:
        results = list(pool.map(lambda j: (j[0], j[2].name, L.download(j[1], j[2])), jobs))
    failed = [x for x in results if x[2].startswith("failed")]
    log(f"downloads: {sum(x[2] == 'ok' for x in results)} new, "
        f"{sum(x[2] == 'cached' for x in results)} cached, {len(failed)} failed")
    for f in failed:
        log(f"  FAILED {f}")
    return len(failed)


# ------------------------------------------------------------------ verify
def stage_verify() -> dict:
    rows = [r for r in read_tsv(L.OUT_DIR / "proteome_manifest.tsv") if r["status"] == "proteome"]
    files, acc_to_species = [], {}
    for r in rows:
        paths = L.proteome_files(r["upid"], int(r["proteome_taxid"]), r["superregnum"])
        fa, g2a = paths["fasta"][1], paths["gene2acc"][1]
        accs = L.fasta_accessions(fa)
        for a in accs:
            acc_to_species[a] = r["species"]
        gg = L.gene_groups(g2a)
        rec = {"species": r["species"], "upid": r["upid"], "fasta": fa.name,
               "fasta_bytes": fa.stat().st_size, "fasta_sha256": sha256(fa),
               "entries": len(accs), "entries_expected": int(r["readme_canonical"]),
               "gene2acc_rows": gg["rows"], "gene2acc_expected": int(r["readme_gene2acc"]),
               "gene_groups": gg["gene_groups"], "unknown_gene_rows": gg["unknown_gene_rows"],
               "in_census_v2": 0, "in_census_v2_family": 0}
        sums = L.fetch_metalink(r["upid"], r["superregnum"])
        md5_ok = all(sums.get(p.name) == L.md5(p) for p in (fa, g2a))
        modified = L.fetch_last_modified(paths["fasta"][0], fa.name)
        agrees = (rec["entries"] == rec["entries_expected"] == len(set(accs))
                  and rec["gene2acc_rows"] == rec["gene2acc_expected"])
        rec.update({"md5_ok": int(md5_ok), "fasta_last_modified": modified,
                    "reissued": int(modified > L.RELEASE_DATE)})
        rec["verdict"] = L.file_verdict(md5_ok, agrees, bool(rec["reissued"]),
                                        rec["entries"], int(r["gene_count"]))
        files.append(rec)
    log("census v2 overlap (one streaming pass)")
    by_sp = {f["species"]: f for f in files}
    n_census = 0
    for c in iter_census_v2(("accession", "family")):
        n_census += 1
        sp = acc_to_species.get(c["accession"])
        if sp:
            by_sp[sp]["in_census_v2"] += 1
            by_sp[sp]["in_census_v2_family"] += int(bool(c["family"]))
    write_tsv(L.OUT_DIR / "proteome_files.tsv", FILE_FIELDS, files)
    # positive control: every enumerated human census gene is in the proteome
    human = [r["accession"] for r in read_tsv(L.ROOT / "results" / "census_v2" /
                                              "human_recall.tsv") if r["accession"]]
    human_in = sum(acc_to_species.get(a) == "Homo sapiens" for a in human)

    db = L.proteome_dir() / "panel_refprot.fasta"
    n_seq = 0
    with db.open("w") as out:
        for f in files:
            with gzip.open(L.proteome_dir() / f["fasta"], "rt") as fh:
                for line in fh:
                    n_seq += line.startswith(">")
                    out.write(line)
    summary = {
        "release": L.RELEASE, "species": len(read_tsv(L.OUT_DIR / "proteome_manifest.tsv")),
        "proteomes": len(files), "proteomes_exact": sum(f["verdict"] == "exact" for f in files),
        "proteomes_reissued": [f["species"] for f in files if f["verdict"] == "reissued"],
        "proteomes_failed": [f["species"] for f in files if f["verdict"].startswith("FAIL")],
        "entries": sum(f["entries"] for f in files),
        "census_v2_records": n_census,
        "census_v2_in_panel": sum(f["in_census_v2"] for f in files),
        "human_census_genes_enumerated": len(human),
        "human_census_genes_in_proteome": human_in,
        "sweep_db": str(db), "sweep_db_seqs": n_seq, "sweep_db_bytes": db.stat().st_size,
        "sweep_db_sha256": sha256(db),
        "download_bytes": sum(f["fasta_bytes"] for f in files)}
    (L.OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1))
    log(f"{summary['proteomes_exact']}/{summary['proteomes']} proteomes exact, "
        f"reissued {summary['proteomes_reissued']}, failed {summary['proteomes_failed']}; "
        f"{summary['entries']:,} entries (one per gene); "
        f"sweep DB {n_seq:,} seqs")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["manifest", "download", "verify", "all"])
    ap.add_argument("--refresh", action="store_true", help="re-query the APIs")
    a = ap.parse_args()
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = {"manifest": False, "download": False, "verify": False}
    if a.stage in ("manifest", "all"):
        stage_manifest(a.refresh)
        steps["manifest"] = True
        progress(steps)
    if a.stage in ("download", "all"):
        if stage_download():
            return 1
        steps["download"] = True
        progress(steps)
    if a.stage in ("verify", "all"):
        s = stage_verify()
        steps["verify"] = True
        progress(steps)
        if s["proteomes_failed"]:
            log(f"FAILED {s['proteomes_failed']} — see proteome_files.tsv")
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
