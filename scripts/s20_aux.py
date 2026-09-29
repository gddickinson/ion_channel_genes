"""s20_aux.py — S20 step 3: are the auxiliary families homologous groups, and
how many auxiliaries does the panel carry?

    python3 scripts/s20_aux.py groups     # → aux_groups.tsv (+ aux_pairs.tsv)
    python3 scripts/s20_aux.py panel      # → aux_panel.tsv, aux_by_species.tsv

**groups.** The catalogue's `channel_associated` families are defined by
*which channel they serve*, not by descent: `assoc_k_beta` holds Kvβ
aldo-keto reductases, BK β subunits, LRRC γ subunits and DPP6/10. A profile
built across unrelated proteins is not a detector (S3a: their LOO decoys got
no hit). The test is direct: every human member against every other of its
family with `phmmer` (E ≤ 1e-3 either way = homologous); homology groups
are the connected components. A family with > 1 group pools unrelated
proteins.

**panel.** Every S3b panel entry the S3a profiles call to an auxiliary
family (any confidence) is searched with `phmmer` against the **whole human
reference proteome**. It is an auxiliary of group G only if its best human
hit (E ≤ 1e-5) is a human member of G — a reciprocal-best-hit test that an
LRR, Ig or aldo-keto-reductase protein fails, because its best human hit is
another LRR, Ig or AKR protein. Otherwise `not_auxiliary` (best human hit
named) or `unplaced` (no human hit), counted, never assigned.
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import iter_fasta, read_tsv, write_tsv  # noqa: E402
from s3b_lib import accession, panel_db, species_by_taxid  # noqa: E402
from s20_lib import OUT_DIR, ROOT, raw_dir  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

E_GROUP, E_PLACE = 1e-3, 1e-5
SEEDS = ROOT / "results" / "census_v3" / "seed_human_accessions.tsv"
AUX = sorted(f.key for f in CATALOGUE.values() if f.status.name == "CHANNEL_ASSOCIATED")


def human_aux() -> tuple[dict[str, dict], dict[str, str]]:
    """gene → {family, accession}; and accession → sequence from the panel DB."""
    genes = {r["gene"]: r for r in read_tsv(SEEDS) if r["family"] in AUX}
    accs = {r["accession"] for r in genes.values() if r["accession"]}
    seqs = {}
    for h, s in iter_fasta(panel_db()):
        a = accession(h.split()[0])
        if a in accs:
            seqs[a] = s
    return genes, seqs


def phmmer(query: Path, db: Path, cpu: int = 4) -> list[tuple[str, str, float, float]]:
    """(query, target, full-sequence E, bits) for every reported pair."""
    with tempfile.NamedTemporaryFile(suffix=".tbl", delete=False) as t:
        tbl = t.name
    p = subprocess.run(["phmmer", "--noali", "--cpu", str(cpu), "-E", "1",
                        "--tblout", tbl, str(query), str(db)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"phmmer failed: {p.stderr[:300]}")
    out = []
    for line in open(tbl):
        if line.startswith("#"):
            continue
        f = line.split()
        out.append((f[2], f[0], float(f[4]), float(f[5])))
    return out


def _components(nodes: list[str], edges: set[tuple[str, str]]) -> list[list[str]]:
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    groups = defaultdict(list)
    for n in nodes:
        groups[find(n)].append(n)
    return sorted((sorted(g) for g in groups.values()), key=lambda g: g[0])


def cmd_groups(_a) -> None:
    genes, seqs = human_aux()
    d = raw_dir("aux")
    fa = d / "human_aux.fasta"
    gene_of = {r["accession"]: g for g, r in genes.items() if r["accession"] in seqs}
    with open(fa, "w") as fh:
        for a, g in sorted(gene_of.items(), key=lambda x: x[1]):
            fh.write(f">{a}\n{seqs[a]}\n")
    hits = phmmer(fa, fa)
    fam = {a: genes[g]["family"] for a, g in gene_of.items()}
    pairs, edges = [], set()
    best: dict[tuple[str, str], float] = {}
    for q, t, e, bits in hits:
        if q == t:
            continue
        k = tuple(sorted((q, t)))
        best[k] = min(best.get(k, 1e9), e)
    for (a, b), e in sorted(best.items()):
        pairs.append({"gene_a": gene_of[a], "gene_b": gene_of[b], "family_a": fam[a],
                      "family_b": fam[b], "evalue": f"{e:.2g}"})
        if e <= E_GROUP and fam[a] == fam[b]:
            edges.add((a, b))
    write_tsv(OUT_DIR / "aux_pairs.tsv", ["gene_a", "gene_b", "family_a", "family_b",
                                          "evalue"], pairs)
    rows = []
    for f in AUX:
        members = [a for a in gene_of if fam[a] == f]
        comps = _components(members, {e for e in edges if e[0] in members})
        no_seq = [g for g, r in genes.items() if r["family"] == f and r["accession"] not in seqs]
        for i, comp in enumerate(comps, 1):
            rows.append({"family": f, "n_groups": len(comps), "group": f"{f}.{i}",
                         "genes": ",".join(sorted(gene_of[a] for a in comp)),
                         "accessions": ",".join(sorted(comp)),
                         "not_tested": ",".join(no_seq) if i == 1 else ""})
    write_tsv(OUT_DIR / "aux_groups.tsv", ["family", "n_groups", "group", "genes",
                                           "accessions", "not_tested"], rows)
    for r in rows:
        print(r["group"], r["n_groups"], r["genes"], r["not_tested"])


def cmd_panel(_a) -> None:
    groups = read_tsv(OUT_DIR / "aux_groups.tsv")
    group_of = {a: r["group"] for r in groups for a in r["accessions"].split(",") if a}
    calls = {}
    with gzip.open(raw_dir().parent.parent / "hmmer" / "s3b" / "panel_profile_calls.tsv.gz",
                   "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["p_call"] == "family" and r["p_family"] in AUX:
                calls[r["target"]] = r
    species = species_by_taxid()
    d = raw_dir("aux")
    fa = d / "panel_aux_calls.fasta"
    taxon = {}
    with open(fa, "w") as fh:
        for h, s in iter_fasta(panel_db()):
            t = h.split()[0]
            if t in calls:
                fh.write(f">{t}\n{s}\n")
                m = re.search(r"OX=(\d+)", h)
                taxon[t] = m.group(1) if m else ""
    human = d / "human_proteome.fasta"
    gene_name = {}
    with open(human, "w") as fh:
        for h, s in iter_fasta(panel_db()):
            if " OX=9606 " in h:
                a = accession(h.split()[0])
                m = re.search(r"GN=(\S+)", h)
                gene_name[a] = m.group(1) if m else ""
                fh.write(f">{a}\n{s}\n")
    hits = phmmer(fa, human, cpu=8)
    best: dict[str, tuple[float, str]] = {}
    for q, t, e, bits in hits:
        if e <= E_PLACE and (q not in best or bits > best[q][0]):
            best[q] = (bits, t)
    rows = []
    for t, r in sorted(calls.items()):
        sp = species.get(taxon.get(t, ""), {})
        hit = best.get(t)
        g = group_of.get(hit[1], "") if hit else ""
        placed = (g if g and g.split(".")[0] == r["p_family"] else
                  "not_auxiliary" if hit else "unplaced")
        rows.append({"target": t, "species": sp.get("species", ""),
                     "group": sp.get("group", ""), "family": r["p_family"],
                     "confidence": r["p_confidence"], "coverage": r["win_coverage"],
                     "homology_group": placed,
                     "best_human": hit[1] if hit else "",
                     "best_human_gene": gene_name.get(hit[1], "") if hit else "",
                     "bits": hit[0] if hit else ""})
    write_tsv(OUT_DIR / "aux_panel.tsv", list(rows[0]), rows)
    c = Counter((r["family"], r["homology_group"], r["confidence"]) for r in rows)
    for k, v in sorted(c.items()):
        print(k, v)
    by = Counter((r["homology_group"], r["species"], r["group"]) for r in rows
                 if r["confidence"] == "high")
    write_tsv(OUT_DIR / "aux_by_species.tsv", ["homology_group", "species", "panel_group",
                                               "high_confidence_calls"],
              [{"homology_group": g, "species": s, "panel_group": pg,
                "high_confidence_calls": n} for (g, s, pg), n in sorted(by.items())])


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("groups")
    sub.add_parser("panel")
    a = ap.parse_args()
    {"groups": cmd_groups, "panel": cmd_panel}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
