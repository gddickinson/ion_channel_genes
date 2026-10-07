"""s12_network.py — S12 (4)–(5): S8a's fold network re-read with experimental
units and with TM-region units (D54). S8a's verdicts stay primary.

    bin/envpy scripts/s12_network.py [--jobs 8]
        → results/structures/network_pairs_<variant>.tsv, network_edges.tsv,
          network_units.tsv

Variants (D48's edge rule unchanged in each):

* ``experimental`` — every node with an experimental unit
  (`experimental.tsv`, D54 (2)) uses it; the rest keep S8a's AFDB unit.
  RyR enters only here.
* ``tm_region`` — every whole-model node cut to its exemplar's UniProt
  TRANSMEM/INTRAMEM span on S8a's AFDB model (pLDDT ≥ 70); module nodes keep
  S8a's unit.
"""

from __future__ import annotations

import argparse
import itertools
from concurrent.futures import ThreadPoolExecutor

from s8_fold_network import _tmalign, edge_rows, node_superfamily, read_ca, sf_medians
from s12_lib import (MIN_PLDDT, OUT, S8_NET, has_module, read_tsv, sdir, tm_region,
                     write_tsv)

MIN_RES = 30


def s8_units() -> dict:
    d = sdir().parent / "s8" / "units"
    return {r["node"]: (d / f"{r['node'].replace(':', '_')}.pdb", "afdb", int(r["residues"]))
            for r in read_tsv(S8_NET / "units.tsv")
            if r["residues"] and int(r["residues"]) >= MIN_RES}


def experimental_units() -> dict:
    units = s8_units()
    for r in read_tsv(OUT / "experimental.tsv"):
        if r["pdb"] and int(r["observed"]) >= MIN_RES:
            units[r["node"]] = (sdir("exp_units") / f"{r['node'].replace(':', '_')}.pdb",
                                f"pdb:{r['pdb']}_{r['chain']}", int(r["observed"]))
    return units


def tm_region_units() -> dict:
    units = s8_units()
    afdb = sdir().parent / "s8" / "afdb"
    for s in read_tsv(S8_NET / "structures.tsv"):
        node = s["node"]
        if not s["model"] or has_module(node):
            continue
        rng = tm_region(s["accession"])
        if rng is None:
            units.pop(node, None)
            continue
        keep = {n for n, _, pl in read_ca(afdb / s["model"])
                if rng[0] <= n <= rng[1] and pl >= MIN_PLDDT}
        lines = [l for l in (afdb / s["model"]).read_text().splitlines()
                 if l.startswith("ATOM") and int(l[22:26]) in keep]
        out = sdir("tm_units") / f"{node}.pdb"
        out.write_text("\n".join(lines) + "\nEND\n")
        n = len(keep)
        if n >= MIN_RES:
            units[node] = (out, f"tm_region:{rng[0]}-{rng[1]}", n)
        else:
            units.pop(node, None)
    return units


def measure(units: dict, jobs: int) -> list[dict]:
    pairs = list(itertools.combinations(sorted(units), 2))
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        res = list(ex.map(lambda p: _tmalign(units[p[0]][0], units[p[1]][0]), pairs))
    return [{"a": x, "b": y, "sf_a": node_superfamily(x), "sf_b": node_superfamily(y),
             "tm_avg": ta, "tm_short": max(t1, t2), "tm_long": min(t1, t2)}
            for (x, y), (t1, t2, ta) in zip(pairs, res)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    edges = [dict(r, variant="s8a_afdb") for r in read_tsv(S8_NET / "literature_edges.tsv")]
    unit_rows, medians = [], []
    for variant, build in (("experimental", experimental_units),
                           ("tm_region", tm_region_units)):
        units = build()
        unit_rows += [{"variant": variant, "node": n, "source": s, "residues": k}
                      for n, (_, s, k) in sorted(units.items())]
        pairs = measure(units, a.jobs)
        write_tsv(OUT / f"network_pairs_{variant}.tsv", list(pairs[0]), pairs)
        for r in edge_rows(pairs):
            edges.append(dict(r, variant=variant))
            print(f"{variant:13s} {r['a']:12s} {r['b']:12s} {r['median_tm']!s:7s} "
                  f"{r['rank_b_for_a']}/{r['rank_a_for_b']} {r['verdict']}", flush=True)
        medians += [{"variant": variant, "a": k[0], "b": k[1], "median_tm": round(v, 4)}
                    for k, v in sorted(sf_medians(pairs).items())]
    cols = ["variant"] + [c for c in edges[0] if c != "variant"]
    write_tsv(OUT / "network_edges.tsv", cols, edges)
    write_tsv(OUT / "network_units.tsv", list(unit_rows[0]), unit_rows)
    write_tsv(OUT / "network_medians.tsv", list(medians[0]), medians)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
