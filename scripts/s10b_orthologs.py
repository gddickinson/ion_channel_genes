"""S10b — S5's absences of ZAC, PACC1 and CLCC1 against two ortholog databases
(D51 (4)).

    bin/envpy scripts/s10b_orthologs.py fetch    # NCBI + Ensembl, archived
    bin/envpy scripts/s10b_orthologs.py check    # per cell: agrees / disputed / not covered
                                                 # + every disputed protein scored (D32)
                                                 #   and placed on the S5 genome
    bin/envpy scripts/s10b_orthologs.py screen   # absent cells under a runner-up bait

Targets are every S5 `absent` cell of the three families (cells.tsv), so the
list is S5's, not chosen here. Gene symbols are read from the catalogue and
used only to *query* the databases — never to call anything (H15).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import s10b_lib as L  # noqa: E402
from s3_assign import assign_one, collect_hits, fam_superfamily  # noqa: E402
from s3_hmm_lib import read_tsv, s3_dir, write_fasta, write_tsv  # noqa: E402
from s5_classify import PANEL_Z  # noqa: E402
from s5_genome_io import fna_path  # noqa: E402
from s5_ledger import load_loci  # noqa: E402
from s5_lib import manifest, sweep_assembly  # noqa: E402
from s5_rescue import read_hsps, run_tblastn  # noqa: E402
from src.catalogue import registry  # noqa: E402
from src.utils.species import ensembl_species_slug  # noqa: E402

FAMILIES = ("zac", "pacc", "clcc1")
OUT = ROOT / "results" / "repertoire"
NCBI = "https://api.ncbi.nlm.nih.gov/datasets/v2"
ENSEMBL = "https://rest.ensembl.org"


def targets() -> list[dict]:
    return [c for c in read_tsv(ROOT / "results" / "genome_sweep" / "cells.tsv")
            if c["family"] in FAMILIES and c["verdict"] == "absent"]


def human_gene(fam: str) -> str:
    genes = list(registry.CATALOGUE[fam].human_genes)
    if len(genes) != 1:
        raise RuntimeError(f"{fam}: expected one human gene, got {genes}")
    return genes[0]


def ncbi_orthologs(sym: str) -> dict:
    d = L.raw_dir("orthologs")
    g = L.fetch_json(f"{NCBI}/gene/symbol/{sym}/taxon/9606", d / f"ncbi_{sym}_human.json")
    gid = g["reports"][0]["gene"]["gene_id"]
    reps, token, page = [], "", 0
    while True:
        url = f"{NCBI}/gene/id/{gid}/orthologs?page_size=1000"
        url += f"&page_token={token}" if token else ""
        j = L.fetch_json(url, d / f"ncbi_{sym}_orthologs_{page}.json")
        reps += [r["gene"] for r in j.get("reports", [])]
        token, page = j.get("next_page_token", ""), page + 1
        if not token:
            break
    return {"gene_id": gid, "orthologs": reps}


def ensembl_orthologs(sym: str) -> list[dict]:
    j = L.fetch_json(f"{ENSEMBL}/homology/symbol/human/{sym}?type=orthologues;"
                     "format=condensed;content-type=application/json",
                     L.raw_dir("orthologs") / f"ensembl_{sym}.json")
    return [h for d in j.get("data", []) for h in d.get("homologies", [])]


def ensembl_species() -> set[str]:
    j = L.fetch_json(f"{ENSEMBL}/info/species?content-type=application/json",
                     L.raw_dir("orthologs") / "ensembl_species.json")
    return {s["name"] for s in j["species"]}


def cmd_fetch(_a) -> None:
    for fam in FAMILIES:
        sym = human_gene(fam)
        n = ncbi_orthologs(sym)
        e = ensembl_orthologs(sym)
        print(fam, sym, "NCBI gene", n["gene_id"], "orthologs", len(n["orthologs"]),
              "| Ensembl", len(e))
    print("Ensembl species:", len(ensembl_species()))


def _ncbi_protein(gene: dict) -> tuple[str, str]:
    d = L.raw_dir("orthologs")
    j = L.fetch_json(f"{NCBI}/gene/id/{gene['gene_id']}/product_report",
                     d / f"ncbi_product_{gene['gene_id']}.json")
    best = ("", 0)
    for rep in j.get("reports", []):
        for t in rep.get("product", {}).get("transcripts", []):
            p = t.get("protein", {})
            if p.get("accession_version") and p.get("length", 0) > best[1]:
                best = (p["accession_version"], p["length"])
    if not best[0]:
        return "", ""
    fa = L.fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
                 f"&id={best[0]}&rettype=fasta&retmode=text", d / f"{best[0]}.fasta")
    return best[0], "".join(fa.splitlines()[1:])


def _ensembl_protein(pid: str) -> str:
    fa = L.fetch(f"{ENSEMBL}/sequence/id/{pid}?type=protein;content-type=text/x-fasta",
                 L.raw_dir("orthologs") / f"{pid}.fasta")
    return "".join(fa.splitlines()[1:])


def score(items: list[tuple[str, str]], work: Path) -> dict[str, dict]:
    """D32 on each protein, with S5's hmmsearch settings."""
    faa, dom = work / "disputed.faa", work / "disputed.domtbl"
    write_fasta(faa, items)
    subprocess.run(["hmmsearch", "--cpu", "6", "--noali", "-E", "1e-3", "--domE", "1e-3",
                    "-Z", str(PANEL_Z), "--domZ", str(PANEL_Z), "-o", "/dev/null",
                    "--domtblout", str(dom), str(s3_dir("profiles") / "all.hmm"),
                    str(faa)], check=True)
    sf = fam_superfamily()
    return {t: assign_one(t, h, sf) for t, h in collect_hits([dom]).items()}


def place(species: str, items: list[tuple[str, str]], work: Path) -> dict[str, dict]:
    """Best tblastn HSP of each protein on S5's genome, and the S5 locus under it."""
    acc = sweep_assembly(manifest()[species])
    faa = work / f"{species.replace(' ', '_')}.faa"
    write_fasta(faa, items)
    hs = read_hsps(run_tblastn(acc, fna_path(acc), [i for i, _ in items],
                               "s10b_" + species.replace(" ", "_"), faa=faa))
    loci = load_loci(acc)
    out = {}
    for pid, _ in items:
        mine = [h for h in hs if h["qseqid"] == pid]
        if not mine:
            out[pid] = {"hit": "none"}
            continue
        b = max(mine, key=lambda h: float(h["bitscore"]))
        under = [x for x in loci if x["contig"] == b["contig"]
                 and int(x["start"]) <= b["end"] and b["start"] <= int(x["end"])]
        out[pid] = {"hit": f"{b['contig']}:{b['start']}-{b['end']}",
                    "bits": b["bitscore"], "pident": b["pident"],
                    "s5_locus": ";".join(f"{x['locus']}={x['p_call']}:{x['p_family']}"
                                         f":{x['p_confidence']}" for x in under) or "none"}
    return out


def cmd_check(_a) -> None:
    man, ens_sp = manifest(), ensembl_species()
    rows, disputed = [], []
    for c in targets():
        fam, sp = c["family"], c["species"]
        sym = human_gene(fam)
        n = ncbi_orthologs(sym)
        taxid = int(man[sp]["panel_taxid"])
        n_hits = [g for g in n["orthologs"] if int(g["tax_id"]) == taxid]
        slug = ensembl_species_slug(sp)
        e_hits = [h for h in ensembl_orthologs(sym) if h["species"] == slug]
        ncbi_cov = "listed" if n_hits else "none"
        ens_cov = ("listed" if e_hits else "none") if slug in ens_sp else "not_in_ensembl"
        # NCBI's coverage of a species is not knowable from one gene's list,
        # so "not covered" is reserved for a species Ensembl lacks and NCBI's
        # ortholog set for the gene names no species of its group at all
        group_taxa = {int(man[s]["panel_taxid"]) for s in man
                      if man[s]["group"] == c["group"]}
        ncbi_group = any(int(g["tax_id"]) in group_taxa for g in n["orthologs"])
        status = ("disputed" if n_hits or e_hits else
                  "agrees" if ens_cov != "not_in_ensembl" or ncbi_group else "not_covered")
        row = {"family": fam, "gene": sym, "species": sp, "group": c["group"],
               "s5_verdict": c["verdict"], "ncbi": ncbi_cov,
               "ncbi_ids": ",".join(str(g["gene_id"]) for g in n_hits),
               "ensembl": ens_cov, "ensembl_ids": ",".join(h["id"] for h in e_hits),
               "status": status}
        rows.append(row)
        for g in n_hits:
            pid, seq = _ncbi_protein(g)
            if seq:
                disputed.append((row, f"ncbi:{pid}", seq, g.get("symbol", "")))
        for h in e_hits:
            seq = _ensembl_protein(h["protein_id"])
            disputed.append((row, f"ensembl:{h['protein_id']}", seq, h["id"]))
    work = L.raw_dir("orthologs", "work")
    calls = score([(d[1], d[2]) for d in disputed], work) if disputed else {}
    placed = {}
    for sp in sorted({d[0]["species"] for d in disputed}):
        placed.update(place(sp, [(d[1], d[2]) for d in disputed
                                 if d[0]["species"] == sp], work))
    drows = []
    for row, pid, seq, label in disputed:
        cl, pl = calls.get(pid, {}), placed.get(pid, {})
        drows.append({"family": row["family"], "species": row["species"], "protein": pid,
                      "db_label": label, "length": len(seq),
                      "p_call": cl.get("p_call", "no_hit"), "p_family": cl.get("p_family", ""),
                      "p_confidence": cl.get("p_confidence", ""),
                      "win_score": cl.get("win_score", ""), "runner": cl.get("runner", ""),
                      "genome_hit": pl.get("hit", ""), "genome_bits": pl.get("bits", ""),
                      "genome_pident": pl.get("pident", ""),
                      "s5_locus": pl.get("s5_locus", "")})
    write_tsv(OUT / "ortholog_check.tsv", list(rows[0]), rows)
    if drows:
        write_tsv(OUT / "ortholog_disputed.tsv", list(drows[0]), drows)
    print(json.dumps({s: sum(r["status"] == s for r in rows)
                      for s in ("agrees", "disputed", "not_covered")}))


def cmd_screen(_a) -> None:
    """Post-hoc screen (found by the PACC1 check): every S5 `absent` cell whose
    family is the runner-up bait of a locus in that genome — a locus whose
    majority bait was another family may have absorbed it. A screen, not a
    verdict: S5's tables are not edited."""
    man = manifest()
    sf = {k: f.superfamily for k, f in registry.CATALOGUE.items()}
    rows = []
    for c in read_tsv(ROOT / "results" / "genome_sweep" / "cells.tsv"):
        if c["verdict"] != "absent":
            continue
        for x in load_loci(sweep_assembly(man[c["species"]])):
            if x["bait_runner"] == c["family"]:
                rows.append({"species": c["species"], "family": c["family"],
                             "locus": x["locus"], "span": x["span"],
                             "bait_family": x["bait_family"], "p_call": x["p_call"],
                             "p_family": x["p_family"], "p_confidence": x["p_confidence"],
                             "same_superfamily": int(sf.get(x["bait_family"]) ==
                                                     sf.get(c["family"]))})
    write_tsv(OUT / "absent_runner_screen.tsv", list(rows[0]), rows)
    print(len(rows), "loci;", sum(r["same_superfamily"] == 0 for r in rows),
          "across superfamilies")


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for n, fn in (("fetch", cmd_fetch), ("check", cmd_check), ("screen", cmd_screen)):
        sub.add_parser(n).set_defaults(fn=fn)
    a = ap.parse_args()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
