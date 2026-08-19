"""S0 — verify every claim the catalogue makes, against live databases.

The catalogue is a file full of assertions: this Pfam accession exists and
means that; this gene symbol resolves to that protein in that organism; this
family's members carry these domains; this taxon id is that species. Written
by hand, checked by nothing, such a file is a liability — the parent project
learned that when a `[lit]` claim about exon counts turned out to be false
and had propagated into three downstream documents.

So this script re-derives all of it and writes what it found:

    catalogue_validation.tsv   internal consistency (no network)
    pfam_verification.tsv      every accession: does it exist, is the short
                               name what we wrote, how many proteins carry it
    exemplars_resolved.tsv     gene+organism → accession, length, and whether
                               a declared accession actually is that gene
    exemplar_architecture.tsv  observed Pfam sets vs declared signatures
    signature_sharing.tsv      which families share which signature — the
                               hazard registry's empirical basis
    taxonomy_check.tsv         every taxon id in the species table
    reference_panel.fasta      the classifier's reference sequences
    verification_summary.json  headline counts for the report

Nothing here edits the catalogue. A mismatch is reported, and fixing it is a
decision a person makes — a script that silently rewrites its own inputs to
make its checks pass is not a check.

    python3 scripts/s0_catalogue_verify.py [--quick] [--no-sequences]
                                           [--limit N] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import (Fetcher, entry_by_accession, live_progress,
                            pfam_entry, pfam_protein_count, protein_pfams,
                            resolve_gene, sequence, taxon, write_tsv)
from src.catalogue import (CATALOGUE, HAZARDS, SUPERFAMILIES, census_families,
                           exemplars, shared_signatures, signature_index,
                           stats, validate)
from src.utils.species import all_species


def step(msg: str) -> None:
    print(f"[s0] {msg}", flush=True)


# ------------------------------------------------------------------ 1
def check_internal(out: Path) -> list[str]:
    problems = validate()
    rows = [(i + 1, p) for i, p in enumerate(problems)]
    write_tsv(out / "catalogue_validation.tsv", ("n", "problem"), rows)
    step(f"internal validation: {len(problems)} problem(s)")
    return problems


# ------------------------------------------------------------------ 2
def check_pfams(f: Fetcher, out: Path, with_counts: bool) -> dict:
    idx = signature_index()
    declared: dict[str, str] = {}
    for fam in CATALOGUE.values():
        for s in fam.signatures:
            declared.setdefault(s.accession, s.name)
    rows, bad, renamed = [], [], []
    for i, (acc, name) in enumerate(sorted(declared.items()), 1):
        e = pfam_entry(f, acc)
        n_prot = pfam_protein_count(f, acc) if (with_counts and e) else None
        if e is None:
            rows.append((acc, name, "", "", "", "MISSING", len(idx.get(acc, []))))
            bad.append(acc)
        else:
            ok = "ok" if e["short"] == name else "NAME_MISMATCH"
            if ok != "ok":
                renamed.append((acc, name, e["short"]))
            rows.append((acc, name, e["short"], e["name"][:70], n_prot, ok,
                         len(idx.get(acc, []))))
        if i % 20 == 0:
            step(f"  pfam {i}/{len(declared)}")
    write_tsv(out / "pfam_verification.tsv",
              ("accession", "declared_name", "interpro_short", "interpro_name",
               "n_uniprot_proteins", "status", "n_families_using"), rows)
    step(f"pfam: {len(declared)} checked, {len(bad)} missing, "
         f"{len(renamed)} name mismatch")
    return {"n": len(declared), "missing": bad, "renamed": renamed}


# ------------------------------------------------------------------ 3
def check_exemplars(f: Fetcher, out: Path, limit: int,
                    fetch_seqs: bool) -> dict:
    ex = exemplars(census_only=False)
    if limit:
        ex = ex[:limit]
    taxa = {s.scientific: s.taxon_id for s in all_species()}
    rows, arch_rows, fasta, resolved = [], [], [], {}
    n_ok = n_declared_ok = n_declared_bad = n_unresolved = 0

    for i, (fam_key, e) in enumerate(ex, 1):
        tax = taxa.get(e.species)
        rec = None
        route = ""
        if e.uniprot:
            rec = entry_by_accession(f, e.uniprot)
            route = "declared"
        if rec is None and tax:
            rec = resolve_gene(f, e.gene, tax)
            route = "resolved" if rec else ""
        if rec is None and tax:
            rec = resolve_gene(f, e.gene, tax, reviewed_only=False)
            route = "resolved_unreviewed" if rec else ""

        if rec is None:
            rows.append((fam_key, e.label, e.gene, e.species, tax or "",
                         e.uniprot, "", "", "", "", "UNRESOLVED"))
            n_unresolved += 1
            continue

        # UniProt carries no gene name at all for many prokaryotic, viral and
        # invertebrate entries — measured: 7 of this catalogue's exemplars,
        # including NavAb, ELIC, GluR0, FaNaC and the AChBPs. That is a
        # property of the record, not a disagreement with the catalogue, and
        # calling it a mismatch would bury the real mismatches in noise.
        has_gene = bool(rec["gene"])
        gene_ok = has_gene and rec["gene"].upper() == e.gene.upper()
        status = "ok"
        if e.uniprot and route == "declared":
            if gene_ok:
                n_declared_ok += 1
            elif not has_gene:
                status = "NO_GENE_NAME_IN_UNIPROT"
            else:
                status = "DECLARED_GENE_MISMATCH"
                n_declared_bad += 1
        elif not gene_ok:
            status = "NO_GENE_NAME_IN_UNIPROT" if not has_gene \
                else "GENE_NAME_DIFFERS"
        n_ok += 1
        resolved[e.gene] = rec["accession"]
        rows.append((fam_key, e.label, e.gene, e.species, tax or "",
                     e.uniprot, rec["accession"], rec["gene"], rec["length"],
                     rec["protein"][:60], status))

        observed = protein_pfams(f, rec["accession"])
        declared = {s.accession: s.copies for s in CATALOGUE[fam_key].signatures}
        missing = sorted(set(declared) - set(observed))
        extra = sorted(set(observed) - set(declared))
        arch_rows.append((
            fam_key, e.label, rec["accession"],
            ",".join(f"{a}x{n}" for a, n in sorted(observed.items())),
            ",".join(f"{a}x{n}" for a, n in sorted(declared.items())),
            ",".join(missing), ",".join(extra),
            "ok" if not missing else "DECLARED_NOT_OBSERVED"))

        if fetch_seqs:
            s = sequence(f, rec["accession"])
            if s:
                fasta.append(f">{e.label}|{fam_key}|{rec['accession']}\n{s}\n")
        if i % 10 == 0:
            step(f"  exemplar {i}/{len(ex)}  ({f.n_requests} requests)")

    write_tsv(out / "exemplars_resolved.tsv",
              ("family", "label", "gene", "species", "taxon_id",
               "declared_accession", "resolved_accession", "resolved_gene",
               "length_aa", "protein_name", "status"), rows)
    write_tsv(out / "exemplar_architecture.tsv",
              ("family", "label", "accession", "observed_pfam",
               "declared_pfam", "declared_not_observed", "observed_not_declared",
               "status"), arch_rows)
    if fasta:
        (out / "reference_panel.fasta").write_text("".join(fasta))
    step(f"exemplars: {n_ok}/{len(ex)} resolved, {n_unresolved} unresolved, "
         f"{n_declared_bad} declared accession(s) name the wrong gene, "
         f"{len(fasta)} sequences cached")
    return {"n": len(ex), "resolved": n_ok, "unresolved": n_unresolved,
            "declared_ok": n_declared_ok, "declared_bad": n_declared_bad,
            "sequences": len(fasta),
            "arch_mismatch": sum(1 for r in arch_rows if r[-1] != "ok")}


# ------------------------------------------------------------------ 4
def check_taxa(f: Fetcher, out: Path) -> dict:
    rows, bad = [], 0
    sp = all_species()
    for i, s in enumerate(sp, 1):
        t = taxon(f, s.taxon_id)
        if t is None:
            rows.append((s.taxon_id, s.scientific, "", "", "", "NOT_FOUND"))
            bad += 1
            continue
        got = (t["scientific"] or "").lower()
        want = s.scientific.lower()
        ok = "ok" if (want in got or got in want) else "NAME_MISMATCH"
        if ok != "ok":
            bad += 1
        rows.append((s.taxon_id, s.scientific, t["scientific"], t["rank"],
                     t["lineage"], ok))
        if i % 15 == 0:
            step(f"  taxon {i}/{len(sp)}")
    write_tsv(out / "taxonomy_check.tsv",
              ("taxon_id", "declared", "uniprot_scientific", "rank",
               "lineage", "status"), rows)
    step(f"taxonomy: {len(sp)} checked, {bad} problem(s)")
    return {"n": len(sp), "bad": bad}


# ------------------------------------------------------------------ 5
def write_sharing(out: Path) -> int:
    rows = []
    for acc, fams in sorted(shared_signatures().items(),
                            key=lambda kv: (-len(kv[1]), kv[0])):
        statuses = {CATALOGUE[k].status.value for k in fams}
        rows.append((acc, len(fams), ",".join(fams), ",".join(sorted(statuses)),
                     "yes" if len(statuses) > 1 else "no"))
    write_tsv(out / "signature_sharing.tsv",
              ("accession", "n_families", "families", "statuses",
               "crosses_channel_boundary"), rows)
    step(f"signature sharing: {len(rows)} shared signature(s)")
    return len(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "results" / "s0_baseline")
    ap.add_argument("--limit", type=int, default=0,
                    help="only the first N exemplars (smoke test)")
    ap.add_argument("--quick", action="store_true",
                    help="skip the per-Pfam protein counts (slowest part)")
    ap.add_argument("--no-sequences", action="store_true",
                    help="skip building reference_panel.fasta")
    args = ap.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    live = ROOT / "results" / "session_live.json"
    names = ["internal validation", "pfam accessions", "exemplars",
             "taxonomy", "signature sharing"]
    done = [False] * len(names)
    live_progress(live, "S0", list(zip(names, done)))

    t0 = time.time()
    f = Fetcher()
    summary: dict = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
                     "catalogue": stats()}

    problems = check_internal(out)
    summary["internal_problems"] = len(problems)
    done[0] = True
    live_progress(live, "S0", list(zip(names, done)))

    summary["pfam"] = check_pfams(f, out, with_counts=not args.quick)
    done[1] = True
    live_progress(live, "S0", list(zip(names, done)))

    summary["exemplars"] = check_exemplars(f, out, args.limit,
                                           not args.no_sequences)
    done[2] = True
    live_progress(live, "S0", list(zip(names, done)))

    summary["taxonomy"] = check_taxa(f, out)
    done[3] = True
    live_progress(live, "S0", list(zip(names, done)))

    summary["shared_signatures"] = write_sharing(out)
    done[4] = True
    live_progress(live, "S0", list(zip(names, done)))

    summary["requests"] = f.n_requests
    summary["request_failures"] = [{"url": u, "why": w} for u, w in f.failures]
    summary["elapsed_s"] = round(time.time() - t0, 1)
    summary["hazards"] = len(HAZARDS)
    summary["superfamilies_not_alignable"] = [
        k for k, sf in SUPERFAMILIES.items() if not sf.alignable]
    summary["census_families"] = len(census_families())
    (out / "verification_summary.json").write_text(json.dumps(summary, indent=2))

    step(f"done in {summary['elapsed_s']}s, {f.n_requests} requests, "
         f"{len(f.failures)} failure(s) → {out}")
    return 1 if (problems or f.failures) else 0


if __name__ == "__main__":
    sys.exit(main())
