"""s12_coverage.py — S12 (1): AlphaFold DB and PDB coverage of the final census (D54 (1)).

    bin/envpy scripts/s12_coverage.py [--jobs 8]
        → results/structures/coverage_members.tsv, coverage_families.tsv

Frame: S15's census v4 frame (high-confidence census-family profile calls);
one row per proteome member, genome loci counted per family as having no
accession. Per member: AFDB model of exactly that accession, whether its
sequence equals the census sequence (*exact*), global pLDDT and the share of
residues ≥ 70, *usable* (exact and global pLDDT ≥ 70); the PDB ids UniProt
cross-references.
"""

from __future__ import annotations

import argparse
import statistics
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

from s0_lib import live_progress
from s6_lib import LIVE, proteome_sequences
from s12_lib import (CATALOGUE, MIN_PLDDT, OUT, afdb_record, frame_rows,
                     uniprot_pdb, write_tsv)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = frame_rows()
    prot = [r for r in rows if r["source"] == "proteome"]
    genome = defaultdict(int)
    for r in rows:
        if r["source"] == "genome":
            genome[r["v3_family"]] += 1
    seqs = proteome_sequences({r["target"] for r in prot})
    accs = sorted({r["accession"] for r in prot})
    print(f"frame: {len(rows)} rows, {len(prot)} proteome, {len(accs)} accessions", flush=True)

    done = [0]

    def fetch(acc):
        rec = afdb_record(acc)
        done[0] += 1
        if done[0] % 250 == 0:
            live_progress(LIVE, "S12", [(f"AFDB records {done[0]}/{len(accs)}", False)])
            print(f"  {done[0]}/{len(accs)}", flush=True)
        return acc, rec
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        recs = dict(ex.map(fetch, accs))
    pdb = uniprot_pdb(accs)

    members = []
    for r in prot:
        rec, seq = recs[r["accession"]], seqs[r["target"]]
        model_seq = (rec or {}).get("uniprotSequence") or (rec or {}).get("sequence", "")
        exact = bool(rec) and model_seq == seq
        pl = float(rec["globalMetricValue"]) if rec else None
        f70 = (round(rec.get("fractionPlddtConfident", 0)
                     + rec.get("fractionPlddtVeryHigh", 0), 3) if rec else "")
        members.append({
            "family": r["v3_family"], "superfamily": CATALOGUE[r["v3_family"]].superfamily,
            "accession": r["accession"], "species": r["species"], "group": r["group"],
            "length": len(seq), "model": int(bool(rec)), "exact": int(exact),
            "model_version": (rec or {}).get("latestVersion", ""),
            "global_plddt": pl if pl is not None else "", "frac_ge70": f70,
            "usable": int(exact and pl is not None and pl >= MIN_PLDDT),
            "pdb_ids": ";".join(pdb.get(r["accession"], [])),
            "n_pdb": len(pdb.get(r["accession"], []))})
    write_tsv(OUT / "coverage_members.tsv", list(members[0]), members)

    by = defaultdict(list)
    for m in members:
        by[m["family"]].append(m)
    fams = []
    for fam in sorted(set(by) | set(genome),
                      key=lambda f: (CATALOGUE[f].superfamily, f)):
        ms = by.get(fam, [])
        n = len(ms)
        pls = [m["global_plddt"] for m in ms if m["global_plddt"] != ""]

        def frac(k):
            return round(sum(m[k] for m in ms) / n, 4) if n else ""
        fams.append({"family": fam, "superfamily": CATALOGUE[fam].superfamily,
                     "members": n, "genome_loci": genome.get(fam, 0),
                     "model": sum(m["model"] for m in ms), "exact": sum(m["exact"] for m in ms),
                     "usable": sum(m["usable"] for m in ms),
                     "with_pdb": sum(1 for m in ms if m["n_pdb"]),
                     "frac_model": frac("model"), "frac_usable": frac("usable"),
                     "median_plddt": round(statistics.median(pls), 1) if pls else "",
                     "human_with_pdb": sum(1 for m in ms if m["n_pdb"]
                                           and m["species"] == "Homo sapiens")})
    write_tsv(OUT / "coverage_families.tsv", list(fams[0]), fams)
    tot = {k: sum(f[k] for f in fams) for k in ("members", "model", "exact", "usable",
                                                "with_pdb", "genome_loci")}
    print(tot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
