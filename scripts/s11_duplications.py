"""s11_duplications.py — S11a driver (D52): reconcile every tier-1 gene tree
with the NCBI species tree of the S4 panel.

    python3 scripts/s11_duplications.py reconcile   # → recon_trees, duplications, expansions
    python3 scripts/s11_duplications.py pairs       # → human_pairs, ohnolog_compare

All outputs under `results/duplication/`.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_lib as L                                   # noqa: E402
from s0_lib import live_progress                      # noqa: E402
from s7_newick import parse, root_partition           # noqa: E402
from s11_recon import Reconciler, edge_for_split, utree   # noqa: E402

TREE_FIELDS = ["family", "n_tips", "n_species", "root_mode", "s7b_root", "n_optimal_roots",
               "min_dups", "min_losses", "declared_dups", "declared_losses",
               "declared_optimal", "dups_stated", "dups_supported", "dups_unstated",
               "species_specific", "deep_supported", "root_split", "root_small_side", "note"]
DUP_FIELDS = ["family", "dup", "taxon", "taxon_depth", "species_specific", "n_tips",
              "n_species", "score", "supported", "child_ufboot", "n_genome_tips",
              "n_human", "human"]
PAIR_FIELDS = ["family", "acc1", "acc2", "gene1", "gene2", "hgnc1", "hgnc2", "node",
               "taxon", "supported", "window", "ohno_strict", "ohno_relaxed"]
TWO_R = {"Vertebrata", "Gnathostomata"}
BONY = "Euteleostomi"
THREE_R = "Clupeocephala"                 # LCA(Danio, Takifugu) on this species tree


def accepted(row: dict) -> bool:
    """D52 (3): an S7b root is used only where it is one clade and supported."""
    if row["family"] in L.S7D_UNDEFINED or row["root_rule"] != "rooted":
        return False
    if row["outgroup_monophyletic"] != "True":
        return False
    return row["n_outgroup"] == "1" or float(row["root_ufboot"] or 0) >= 95


def window(st, taxon: int) -> str:
    """D52 (6): the ohnologue window a duplication's mapped taxon names."""
    name = st.names[taxon]
    if name in TWO_R:
        return "2R"
    if name == BONY:
        return "bony_vertebrate"
    if name == THREE_R:
        return "3R"
    vert = st.names.index("Vertebrata")
    n = vert
    while n >= 0:
        if n == taxon:
            return "older"
        n = st.parent[n]
    return "younger"


def reconcile_family(row: dict, st, by_code: dict, source: dict, gn: dict) -> dict:
    fam = row["family"]
    tree = parse((L.TREES / f"{fam}.treefile").read_text())
    leaves = tree.leaves()
    og = {x for x in leaves if x.startswith(L.OG_PREFIX)}
    ingroup = [x for x in leaves if x not in og]
    species_of = {x: L.species_of_tip(x, by_code) for x in ingroup}
    t = utree(tree, drop=og)
    rc = Reconciler(t, st, species_of)
    opt, best = rc.optimal()
    out = {"family": fam, "n_tips": len(ingroup), "n_species": len(set(species_of.values())),
           "s7b_root": row["root_rule"], "n_optimal_roots": len(opt),
           "min_dups": best[0], "min_losses": best[1], "note": ""}
    if len(opt) == 1:                  # the reconciliation root's basal split, for S11b
        sides = sorted((rc.names(rc.sub(*k).leaves) for k in (opt[0], opt[0][::-1])), key=len)
        out["root_split"] = f"{len(sides[0])}|{len(sides[1])}"
        sp = Counter(species_of[x] for x in sides[0])
        out["root_small_side"] = "; ".join(f"{s} {n}" for s, n in sp.most_common(4))
    declared = None
    if og and row["root_rule"] == "rooted":
        rp = root_partition(tree, og)
        if rp:
            declared = edge_for_split(rc, set(rp[0][0]))
            r = rc.root(*declared)
            out.update(declared_dups=r.dups, declared_losses=r.losses,
                       declared_optimal=declared in opt)
    if accepted(row):
        roots, out["root_mode"] = [declared], "declared"
    else:
        roots, out["root_mode"] = opt, "reconciliation"
    per_root = [rc.dup_ids(e) for e in roots]
    stated_ids = set(per_root[0]).intersection(*per_root[1:])
    union_ids = set().union(*per_root)
    dups, pairs = [], []
    human = {x for x in ingroup if species_of[x] == "Homo sapiens"}
    for k, ident in enumerate(sorted(stated_ids, key=lambda i: -per_root[0][i].leaves.bit_length())):
        s = per_root[0][ident]
        names = rc.names(s.leaves)
        hum = sorted(x for x in names if x in human)
        dups.append({"family": fam, "dup": k, "taxon": st.names[s.m],
                     "taxon_depth": st.depth[s.m], "species_specific": bool(st.tip[s.m]),
                     "n_tips": len(names), "n_species": bin(s.species).count("1"),
                     "score": rc.score(s), "supported": rc.supported(ident, s),
                     "child_ufboot": ",".join("leaf" if x is None else f"{x:g}"
                                              for x in rc.child_support(ident, s)),
                     "n_genome_tips": sum(source.get(x) == "genome" for x in names),
                     "n_human": len(hum),
                     "human": ",".join(gn.get(_acc(x), _acc(x)) for x in hum)})
        kids = [set(rc.names(rc.sub(*c).leaves)) & human for c in s.kids]
        for a, b in itertools.combinations(range(len(kids)), 2):
            for x in kids[a]:
                for y in kids[b]:
                    pairs.append((x, y, k, s, rc.supported(ident, s)))
    out.update(dups_stated=len(dups), dups_supported=sum(d["supported"] for d in dups),
               dups_unstated=len(union_ids - stated_ids),
               species_specific=sum(d["supported"] and d["species_specific"] for d in dups),
               deep_supported=sum(d["supported"] and not d["species_specific"] for d in dups))
    if declared is None and row["root_rule"] == "rooted":
        out["note"] = "declared outgroup not one clade"
    out["_pairs"] = [(x, y, k, s.m, sup) for x, y, k, s, sup in pairs]
    out["_human"] = sorted(human)
    return out, dups


def _acc(tip: str) -> str:
    return tip.rsplit("__", 1)[0]


def cmd_reconcile(_a) -> None:
    _, st = L.species_tree()
    by_code = {L.code(r["species"]): r["species"] for r in L.panel()}
    source = {r["label"]: r["source"] for r in L.read_tsv(L.MEMBERS)}
    gn = L.human_gn()
    rows = L.read_tsv(L.TIER1)
    trees, dups, pairs = [], [], []
    for i, row in enumerate(rows):
        live_progress(L.LIVE, "S11a", [(f"reconcile {row['family']} ({i + 1}/{len(rows)})", False)])
        out, d = reconcile_family(row, st, by_code, source, gn)
        for x, y, k, m, sup in out.pop("_pairs"):
            pairs.append({"family": row["family"], "acc1": _acc(x), "acc2": _acc(y),
                          "node": k, "taxon": st.names[m], "supported": sup,
                          "window": window(st, m)})
        nh = out.pop("_human")
        # human pairs with no stated duplication between them
        seen = {(p["acc1"], p["acc2"]) for p in pairs if p["family"] == row["family"]}
        for x, y in itertools.combinations(nh, 2):
            a, b = _acc(x), _acc(y)
            if (a, b) not in seen and (b, a) not in seen:
                pairs.append({"family": row["family"], "acc1": a, "acc2": b, "node": "",
                              "taxon": "", "supported": "", "window": "none"})
        trees.append(out)
        dups.extend(d)
        print(f"{row['family']:22s} {out['root_mode']:15s} roots={out['n_optimal_roots']:4d} "
              f"dups={out['dups_stated']:4d} supported={out['dups_supported']:4d}", flush=True)
    L.write_tsv(L.OUT / "recon_trees.tsv", trees, TREE_FIELDS)
    L.write_tsv(L.OUT / "duplications.tsv", dups, DUP_FIELDS)
    exp = Counter((d["family"], d["taxon"], d["taxon_depth"]) for d in dups if d["supported"])
    L.write_tsv(L.OUT / "expansions.tsv",
                [{"family": f, "taxon": t, "taxon_depth": dp, "supported_dups": n}
                 for (f, t, dp), n in sorted(exp.items())],
                ["family", "taxon", "taxon_depth", "supported_dups"])
    L.write_tsv(L.OUT / "human_pairs_raw.tsv", pairs, PAIR_FIELDS)
    live_progress(L.LIVE, "S11a", [(f"reconcile {len(rows)} trees", True)])


def cmd_pairs(_a) -> None:
    """Human paralogue pairs → window, beside OHNOLOGS v2 (never an input)."""
    acc2h, ens2h, hsym = L.hgnc_maps()
    gn = L.human_gn()
    ohno = {}
    for crit in L.OHNO_CRITERIA:
        s = set()
        for p in L.ohnolog_pairs(crit):
            a, b = ens2h.get(p["ens1"]), ens2h.get(p["ens2"])
            if a and b:
                s.add(frozenset((a, b)))
        ohno[crit] = s
    rows = L.read_tsv(L.OUT / "human_pairs_raw.tsv")
    for r in rows:
        h1, h2 = acc2h.get(r["acc1"], ""), acc2h.get(r["acc2"], "")
        r.update(gene1=gn.get(r["acc1"], ""), gene2=gn.get(r["acc2"], ""), hgnc1=h1, hgnc2=h2)
        for crit in L.OHNO_CRITERIA:
            r[f"ohno_{crit}"] = (frozenset((h1, h2)) in ohno[crit]) if h1 and h2 else ""
    L.write_tsv(L.OUT / "human_pairs.tsv", rows, PAIR_FIELDS)
    comp = []
    for crit in L.OHNO_CRITERIA:
        c = Counter((r["window"], r[f"ohno_{crit}"]) for r in rows)
        for (w, o), n in sorted(c.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
            comp.append({"criterion": crit, "window": w, "in_ohnologs": o, "pairs": n})
    # OHNOLOGS pairs between two of our human genes in one family that we never pair
    ours = defaultdict(set)
    for r in rows:
        if r["hgnc1"]:
            ours[r["family"]].add(r["hgnc1"])
        if r["hgnc2"]:
            ours[r["family"]].add(r["hgnc2"])
    fam_of = {h: f for f, hs in ours.items() for h in hs}
    cross = Counter()
    for crit in L.OHNO_CRITERIA:
        for pr in ohno[crit]:
            a, b = sorted(pr)
            if a in fam_of and b in fam_of and fam_of[a] != fam_of[b]:
                cross[crit] += 1
    for crit in L.OHNO_CRITERIA:
        comp.append({"criterion": crit, "window": "across_families", "in_ohnologs": True,
                     "pairs": cross[crit]})
    L.write_tsv(L.OUT / "ohnolog_compare.tsv", comp,
                ["criterion", "window", "in_ohnologs", "pairs"])
    L.write_tsv(L.OUT / "human_genes.tsv", gene_rows(rows, ohno, hsym),
                ["family", "hgnc", "gene", "windows", "has_2R", "ohno_strict", "ohno_relaxed"])
    (L.OUT / "ohnolog_sources.json").write_text(json.dumps(
        {c: len(s) for c, s in ohno.items()}, indent=1) + "\n")
    print(json.dumps(comp, indent=0)[:2000])


def gene_rows(rows: list[dict], ohno: dict, hsym: dict) -> list[dict]:
    """Per human gene: the windows of its paralogue pairs, and whether OHNOLOGS
    lists it in any 2R pair (the per-gene frame, so large families do not
    dominate as they do in the pair counts)."""
    in_ohno = {c: {h for pr in s for h in pr} for c, s in ohno.items()}
    wins = defaultdict(set)
    fam = {}
    for r in rows:
        for h in (r["hgnc1"], r["hgnc2"]):
            if h:
                wins[h].add(r["window"])
                fam[h] = r["family"]
    return [{"family": fam[h], "hgnc": h, "gene": hsym.get(h, ""),
             "windows": ",".join(sorted(w)), "has_2R": "2R" in w,
             "ohno_strict": h in in_ohno["strict"], "ohno_relaxed": h in in_ohno["relaxed"]}
            for h, w in sorted(wins.items(), key=lambda kv: (fam[kv[0]], kv[0]))]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("reconcile").set_defaults(fn=cmd_reconcile)
    sub.add_parser("pairs").set_defaults(fn=cmd_pairs)
    a = ap.parse_args()
    a.fn(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
