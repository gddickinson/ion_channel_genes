"""s0_outgroups.py — S0 verification of the per-family tier-1 outgroups (S7c).

    python3 scripts/s0_outgroups.py

Every inline exemplar of a `RootSet` (`ChannelFamily.root_with`) is checked
live, the way `s0_catalogue_verify.py` checks family exemplars: the UniProt
accession resolves, the declared gene is one of the entry's gene names, the
declared species is in the entry's organism name, and InterPro places at
least one of the set's declared pore models (`RootSet.pfam`) on it. Family-
keyed root sets are only checked to name catalogue families — their
exemplars are the reference panel's, already verified by S0.

Writes `results/s0_baseline/outgroups.tsv` (one row per check) and
`outgroup_panel.fasta` (`>label|rootset|accession`), kept apart from
`reference_panel.fasta` so the classifier never scores against them.
Exits non-zero on any failed check or failed request; never edits the
catalogue.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from s0_lib import UNIPROT, Fetcher, protein_pfams, sequence  # noqa: E402
from s3_hmm_lib import write_tsv  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

OUT = ROOT / "results" / "s0_baseline"
FIELDS = ["family", "root_set", "label", "accession", "declared_gene",
          "entry_gene_names", "declared_species", "organism", "reviewed",
          "length", "pfams", "pore_model_found", "status"]


def entry(f: Fetcher, acc: str) -> dict | None:
    d = f.json(f"{UNIPROT}/uniprotkb/{acc}.json")
    if not d:
        return None
    names = []
    for g in d.get("genes") or []:
        for k in ("geneName",):
            if g.get(k):
                names.append(g[k]["value"])
        for k in ("synonyms", "orderedLocusNames", "orfNames"):
            names += [x["value"] for x in g.get(k, [])]
    org = d.get("organism") or {}
    return {"accession": d.get("primaryAccession", acc), "genes": names,
            "organism": " ".join([org.get("scientificName", ""),
                                  *(org.get("synonyms") or []),
                                  org.get("commonName", "")]),
            "reviewed": d.get("entryType", "").startswith("UniProtKB reviewed"),
            "length": (d.get("sequence") or {}).get("length")}


def main() -> int:
    f = Fetcher()
    rows, fasta, seen = [], [], {}
    for fam in sorted(CATALOGUE):
        for rs in CATALOGUE[fam].root_with:
            for o in rs.families:
                ok = o in CATALOGUE and o != fam
                rows.append({"family": fam, "root_set": rs.name, "label": o,
                             "status": "ok" if ok else "unknown_family"})
            for e in rs.exemplars:
                r = seen.get(e.uniprot) or entry(f, e.uniprot)
                seen[e.uniprot] = r
                row = {"family": fam, "root_set": rs.name, "label": e.label,
                       "accession": e.uniprot, "declared_gene": e.gene,
                       "declared_species": e.species}
                if r is None:
                    rows.append({**row, "status": "unresolved"})
                    continue
                pf = protein_pfams(f, e.uniprot)
                found = sorted(set(pf) & set(rs.pfam))
                problems = []
                if r["accession"] != e.uniprot:
                    problems.append("accession_changed")
                if e.gene not in r["genes"]:
                    problems.append("gene_mismatch")
                if e.species.lower() not in r["organism"].lower():
                    problems.append("species_mismatch")
                if rs.pfam and not found:
                    problems.append("no_pore_model")
                rows.append({**row, "entry_gene_names": ",".join(r["genes"]),
                             "organism": r["organism"].strip(),
                             "reviewed": r["reviewed"], "length": r["length"],
                             "pfams": ",".join(f"{k}x{v}" for k, v in sorted(pf.items())),
                             "pore_model_found": ",".join(found),
                             "status": ";".join(problems) or "ok"})
                if not problems and e.label not in {l for l, *_ in fasta}:
                    fasta.append((e.label, rs.name, e.uniprot, sequence(f, e.uniprot)))
    write_tsv(OUT / "outgroups.tsv", FIELDS, rows)
    with open(OUT / "outgroup_panel.fasta", "w") as fh:
        for label, rs, acc, seq in fasta:
            fh.write(f">{label}|{rs}|{acc}\n{seq}\n")
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"{len(rows)} checks, {len(bad)} failed, {f.n_requests} requests, "
          f"{len(f.failures)} request failures; {len(fasta)} sequences → "
          f"{OUT / 'outgroup_panel.fasta'}")
    for r in bad:
        print(f"  FAIL {r['family']} {r['root_set']} {r['label']}: {r['status']}")
    return 1 if bad or f.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
