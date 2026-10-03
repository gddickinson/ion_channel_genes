"""S10b — the embedded-contig test (D51 (1)) and the Q4 verdict (D51 (2)).

    bin/envpy scripts/s10b_contig.py subjects   # → results/repertoire/contig_subjects.tsv
    bin/envpy scripts/s10b_contig.py test       # fetch contigs, tblastn, blastx → regions
    bin/envpy scripts/s10b_contig.py q4         # D50 (6) / D51 (2) → q4.json

Subjects: every metazoan MscS locus S5 calls MscS (any confidence), every
metazoan S4b order's MscS call (any confidence), and — descriptively — the
metazoan S4b high-confidence calls of kcsa_prok, iglur_prok and mscl.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import s10b_lib as L  # noqa: E402
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s5_genome_io import build_fai, fetch_region, fna_path, read_fai  # noqa: E402
from s5_lib import manifest, s5_dir, sweep_assembly  # noqa: E402

OUT = ROOT / "results" / "repertoire"
SUBJECTS = OUT / "contig_subjects.tsv"
REGIONS = OUT / "contig_regions.tsv"
Q4 = OUT / "q4.json"
Q4_FAMILY = "mscs"
PROK_FAMILIES = ("kcsa_prok", "iglur_prok", "mscl")
SUBJECT_FIELDS = ["subject", "kind", "family", "confidence", "species", "taxid",
                  "order", "source", "accession", "contig"]


def _s4b(name: str):
    return L.data_root() / "proteomes" / "s4b" / name


def _s5_subjects() -> list[dict]:
    man, out = manifest(), []
    for sp, row in man.items():
        if L.lineage_of_group(row["group"]) != "Metazoa":
            continue
        p = s5_dir(sweep_assembly(row)) / "loci.tsv.gz"
        if not p.exists():
            continue
        for r in read_tsv(p):
            if r["p_call"] == "family" and r["p_family"] == Q4_FAMILY:
                out.append({"subject": f"{sp}:{r['locus']}", "kind": "q4",
                            "family": Q4_FAMILY, "confidence": r["p_confidence"],
                            "species": sp, "taxid": row["proteome_taxid"],
                            "order": "", "source": "s5_locus",
                            "accession": sweep_assembly(row),
                            "contig": f"{r['contig']}:{r['start']}-{r['end']}"})
    return out


def _s4b_subjects() -> list[dict]:
    uni = {r["target"]: r for r in read_tsv(_s4b("universe.tsv.gz"))
           if r["kingdom"] == "Metazoa"}
    out = []
    for r in read_tsv(_s4b("calls.tsv.gz")):
        u = uni.get(r["target"])
        if not u or r["p_call"] != "family":
            continue
        q4 = r["p_family"] == Q4_FAMILY and r["p_confidence"] in ("high", "medium")
        prok = r["p_family"] in PROK_FAMILIES and r["p_confidence"] == "high"
        if q4 or prok:
            out.append({"subject": r["target"].split("|")[1],
                        "kind": "q4" if q4 else "prok", "family": r["p_family"],
                        "confidence": r["p_confidence"], "species": "",
                        "taxid": "", "order": u["order"], "source": "s4b_proteome",
                        "accession": r["target"].split("|")[1], "contig": ""})
    return out


def cmd_subjects(_a) -> None:
    rows = _s5_subjects() + _s4b_subjects()
    write_tsv(SUBJECTS, SUBJECT_FIELDS, rows)
    print(len(rows), "subjects:",
          {k: sum(r["kind"] == k for r in rows) for k in ("q4", "prok")})


# ------------------------------------------------------------- contigs
def _local_contig(species: str, contig_id: str) -> tuple[str, str] | None:
    row = manifest().get(species)
    if not row:
        return None
    fna = fna_path(sweep_assembly(row))
    if not fna or not fna.exists():
        return None
    idx = read_fai(build_fai(fna))
    hit = next((k for k in idx if k == contig_id or k.split(".")[0] == contig_id), None)
    if not hit:
        return None
    return hit, fetch_region(fna, idx, hit, 1, idx[hit][0])


def _uniprot(acc: str) -> dict:
    return L.fetch_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                        L.raw_dir("uniprot") / f"{acc}.json")


def _protein_fasta(acc: str) -> str:
    return L.fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                   L.raw_dir("uniprot") / f"{acc}.fasta")


def _contig_for(s: dict) -> tuple[str, str, dict]:
    """(contig id, sequence, extra) for a subject — local genome first, INSDC else."""
    if s["source"] == "s5_locus":
        cid, span = s["contig"].split(":")
        local = _local_contig(s["species"], cid)
        lo, hi = map(int, span.split("-"))
        return local[0], local[1], {"locus": (lo, hi)}
    rec = _uniprot(s["accession"])
    s["taxid"] = str(rec["organism"]["taxonId"])
    s["species"] = rec["organism"]["scientificName"]
    xs = [x for x in rec.get("uniProtKBCrossReferences", []) if x["database"] == "EMBL"
          and any(p["value"] == "Genomic_DNA" for p in x.get("properties", []))]
    if not xs:
        raise RuntimeError(f"{s['accession']}: no genomic EMBL cross-reference")
    cid = xs[0]["id"]
    local = _local_contig(s["species"], cid)
    if local:
        return local[0], local[1], {}
    # INSDC record by NCBI efetch: ENA's browser API redirects a WGS contig
    # to its whole WGS set (tens of MB, gzipped)
    text = L.fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                   f"?db=nuccore&id={cid}&rettype=fasta&retmode=text",
                   L.raw_dir("insdc") / f"{cid}.fasta")
    seq = "".join(x.strip() for x in text.splitlines()[1:])
    if not seq:
        raise RuntimeError(f"{cid}: efetch returned no sequence")
    return cid, seq, {}


def _write_fa(path: Path, name: str, seq: str) -> Path:
    path.write_text(f">{name}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
    return path


def _locate(s: dict, work: Path, contig: Path) -> tuple[int, int]:
    """The subject protein on its contig: tblastn, best HSP's strand ± 20 kb."""
    prot = work / "protein.fasta"
    prot.write_text(_protein_fasta(s["accession"]))
    hs = [h for h in L.run_blast("tblastn", prot, work / "locate.tsv", subject=contig)
          if h["evalue"] <= 1e-5]
    if not hs:
        raise RuntimeError(f"{s['accession']}: protein not found on its contig")
    best = max(hs, key=lambda h: h["bits"])
    plus = best["ss"] < best["se"]
    near = [h for h in hs if (h["ss"] < h["se"]) == plus
            and abs(min(h["ss"], h["se"]) - min(best["ss"], best["se"])) < 20_000]
    return (min(min(h["ss"], h["se"]) for h in near),
            max(max(h["ss"], h["se"]) for h in near))


def test_one(s: dict, db: Path, ox: dict, taxa: dict) -> tuple[dict, list[dict]]:
    work = L.raw_dir("work", s["subject"].replace(":", "_").replace("|", "_"))
    cid, seq, extra = _contig_for(s)
    contig_fa = _write_fa(work / "contig.fasta", cid, seq)
    lo, hi = extra.get("locus") or _locate(s, work, contig_fa)
    off, query = 0, contig_fa
    if len(seq) > L.WINDOW_ABOVE:
        off = max(0, lo - 1 - L.WINDOW)
        query = _write_fa(work / "window.fasta", f"{cid}_w{off}",
                          seq[off:min(len(seq), hi + L.WINDOW)])
    hs = L.run_blast("blastx", query, work / "blastx.tsv", db=db,
                     extra=["-task", "blastx-fast", "-seg", "yes",
                            "-max_target_seqs", "250"])
    own = s["taxid"]
    hs = [{**h, "qs": h["qs"] + off, "qe": h["qe"] + off} for h in hs
          if h["evalue"] <= L.EVALUE and ox.get(h["s"], "") != own]
    inside = [h for h in hs if max(h["qs"], h["qe"]) >= lo - L.LOCUS_PAD
              and min(h["qs"], h["qe"]) <= hi + L.LOCUS_PAD]
    others = [h for h in hs if h not in inside]
    regions = [L.region_lineage(r, ox, taxa) for r in L.merge_regions(others)]
    own_locus = L.region_lineage(inside, ox, taxa) if inside else {}
    lins = {r["lineage"] for r in regions}
    status = ("embedded" if "Metazoa" in lins else
              "foreign" if lins - {"unclear"} else "no_evidence")
    row = {**s, "contig": cid, "contig_len": len(seq), "searched_bp": len(seq) - off
           if off == 0 else min(len(seq), hi + L.WINDOW) - off,
           "locus_start": lo, "locus_end": hi, "regions": len(regions),
           "metazoan_regions": sum(r["lineage"] == "Metazoa" for r in regions),
           "foreign_regions": sum(r["lineage"] in ("prokaryote", "other_eukaryote", "virus")
                                  for r in regions),
           "unclear_regions": sum(r["lineage"] == "unclear" for r in regions),
           "locus_best_lineage": own_locus.get("best_lineage", ""),
           "locus_best_species": own_locus.get("best_species", ""),
           "locus_best_pident": own_locus.get("best_pident", ""),
           "status": status}
    return row, [{"subject": s["subject"], "contig": cid, **r} for r in regions]


def cmd_test(_a) -> None:
    subs = read_tsv(SUBJECTS)
    db, ox, taxa = L.panel_db(), L.panel_ox(), L.panel_taxa()
    rows, regs = [], []
    for i, s in enumerate(subs, 1):
        try:
            row, rr = test_one(s, db, ox, taxa)
        except Exception as exc:           # recorded, never dropped
            row, rr = {**s, "status": "failed", "note": str(exc)[:200]}, []
        print(f"[{i}/{len(subs)}] {s['subject']} {s['family']} {row['status']}", flush=True)
        rows.append(row)
        regs += rr
    fields = SUBJECT_FIELDS + ["contig_len", "searched_bp", "locus_start", "locus_end",
                               "regions", "metazoan_regions", "foreign_regions",
                               "unclear_regions", "locus_best_lineage",
                               "locus_best_species", "locus_best_pident", "status", "note"]
    write_tsv(OUT / "contig_test.tsv", fields, rows)
    if regs:
        write_tsv(REGIONS, list(regs[0]), regs)


def cmd_q4(_a) -> None:
    """D50 (6) / D51 (2), verbatim."""
    rows = [r for r in read_tsv(OUT / "contig_test.tsv") if r["kind"] == "q4"]
    cells = [c for c in read_tsv(ROOT / "results" / "genome_sweep" / "cells.tsv")
             if c["family"] == Q4_FAMILY
             and L.lineage_of_group(c["group"]) == "Metazoa"]
    embedded_high = [r["subject"] for r in rows if r["status"] == "embedded"
                     and r["confidence"] == "high"]
    embedded_any = [r["subject"] for r in rows if r["status"] == "embedded"]
    not_absent = sorted({c["verdict"] for c in cells if c["verdict"] != "absent"})
    if embedded_high:
        verdict, why = "present", f"embedded high-confidence subject(s): {embedded_high}"
    elif not embedded_any and not not_absent:
        verdict, why = "absent", "no embedded subject and every animal S5 cell absent"
    else:
        reasons = []
        if embedded_any:
            reasons.append(f"embedded but medium-only: {embedded_any}")
        if not_absent:
            n = {v: sum(c["verdict"] == v for c in cells) for v in not_absent}
            reasons.append(f"animal S5 cells not `absent`: {n}")
        verdict, why = "unresolved", "; ".join(reasons)
    out = {"verdict": verdict, "why": why, "subjects": len(rows),
           "by_status": {s: [r["subject"] for r in rows if r["status"] == s]
                         for s in sorted({r["status"] for r in rows})},
           "animal_cells": {c["species"]: c["verdict"] for c in cells}}
    Q4.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ("verdict", "why", "by_status")}, indent=1))


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for n, fn in (("subjects", cmd_subjects), ("test", cmd_test), ("q4", cmd_q4)):
        sub.add_parser(n).set_defaults(fn=fn)
    a = ap.parse_args()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
