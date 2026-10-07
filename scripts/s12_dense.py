"""s12_dense.py — S12 (6): the dense fold network, one AFDB model per census
family × S4 group (D54 (6)).

    bin/envpy scripts/s12_dense.py reps       # → results/structures/dense_reps.tsv
    bin/envpy scripts/s12_dense.py measure    # TM-align + Foldseek all-vs-all
    bin/envpy scripts/s12_dense.py read       # edges, recovery, detectability

Representative: from S12's coverage frame, the *usable* member (exact
sequence, global pLDDT ≥ 70) with the highest global pLDDT, ties by
accession; cut by D48's unit rule (module 1 by S6's `hmmalign` + span, else
the whole model), residues with pLDDT ≥ 70 only, ≥ 30 residues kept.

Read: (a) every literature edge by D48's rule, a family pair's TM-score the
median over its representative pairs (the dense set has no VSD node, so the
Hv edge is unmeasured here); (b) superfamily recovery — a representative of
a superfamily with ≥ 2 measured families is recovered iff its best TM-score
partner outside its own family is in its own superfamily (Foldseek: best
E-value, beside); (c) Foldseek detectability per superfamily pair: share of
representative pairs with E ≤ 1e-3 (best of both directions).
"""

from __future__ import annotations

import argparse
import itertools
import statistics
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from s8_fold_network import _tmalign, edge_rows
from s12_lib import (CATALOGUE, MIN_PLDDT, OUT, afdb_model, afdb_record, has_module,
                     module_range, read_af_ca, read_tsv, sdir, write_tsv)

MIN_RES = 30
E_DETECT = 1e-3


def cmd_reps(_a) -> None:
    best: dict = {}
    for m in read_tsv(OUT / "coverage_members.tsv"):
        if m["usable"] != "1":
            continue
        key = (m["family"], m["group"])
        k = (-float(m["global_plddt"]), m["accession"])
        if key not in best or k < best[key][0]:
            best[key] = (k, m)
    rows = []
    for (fam, grp), (_, m) in sorted(best.items()):
        acc = m["accession"]
        pdb = afdb_model(afdb_record(acc))
        ca = read_af_ca(pdb)
        seq = "".join(a for _, a, _, _ in ca)
        row = {"family": fam, "superfamily": m["superfamily"], "group": grp,
               "accession": acc, "species": m["species"], "global_plddt": m["global_plddt"],
               "unit": "", "start": "", "end": "", "residues": "", "note": ""}
        if has_module(fam):
            rng = module_range(fam, seq, f"dense_{acc}")
            if rng is None:
                row["note"] = "module not found (< 50 % of span)"
                rows.append(row)
                continue
            rng, row["unit"] = (ca[rng[0] - 1][0], ca[rng[1] - 1][0]), "module 1"
        else:
            rng, row["unit"] = (ca[0][0], ca[-1][0]), "whole model"
        keep = {n for n, _, pl, _ in ca if rng[0] <= n <= rng[1] and pl >= MIN_PLDDT}
        lines = [l for l in pdb.read_text().splitlines()
                 if l.startswith("ATOM") and int(l[22:26]) in keep]
        (sdir("dense_units") / f"{fam}__{acc}.pdb").write_text("\n".join(lines) + "\nEND\n")
        row.update(start=rng[0], end=rng[1], residues=len(keep))
        if len(keep) < MIN_RES:
            row["note"] = f"only {len(keep)} residues at pLDDT ≥ {MIN_PLDDT:g}"
        rows.append(row)
        print(f"{fam:20s} {grp:16s} {acc:12s} {len(keep):5d}", flush=True)
    write_tsv(OUT / "dense_reps.tsv", list(rows[0]), rows)


def _reps() -> list[dict]:
    return [r for r in read_tsv(OUT / "dense_reps.tsv")
            if r["residues"] and int(r["residues"]) >= MIN_RES]


def _name(r: dict) -> str:
    return f"{r['family']}__{r['accession']}"


def cmd_measure(a) -> None:
    reps = _reps()
    d = sdir("dense_units")
    names = sorted(_name(r) for r in reps)
    for f in d.glob("*.pdb"):                      # only this run's units are searched
        if f.stem not in names:
            f.unlink()
    pairs = list(itertools.combinations(names, 2))
    print(f"{len(names)} units, {len(pairs)} pairs", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(lambda p: _tmalign(d / f"{p[0]}.pdb", d / f"{p[1]}.pdb"), pairs))
    fs = sdir("dense_foldseek")
    out = fs / "allvsall.m8"
    subprocess.run(["foldseek", "easy-search", str(d), str(d), str(out), str(fs / "tmp"),
                    "--exhaustive-search", "1", "-e", "10", "--threads", str(a.jobs),
                    "--format-output", "query,target,evalue,prob,bits", "-v", "1"],
                   check=True, capture_output=True, text=True)
    ev: dict = {}
    for line in out.read_text().splitlines():
        q, t, e, _, _ = line.split("\t")
        q, t = Path(q).stem, Path(t).stem
        if q != t:
            k = tuple(sorted((q, t)))
            ev[k] = min(float(e), ev.get(k, float("inf")))
    rows = [{"a": x, "b": y, "fam_a": x.split("__")[0], "fam_b": y.split("__")[0],
             "tm_avg": ta, "foldseek_evalue": ev.get((x, y), "")}
            for (x, y), (_, _, ta) in zip(pairs, res)]
    write_tsv(OUT / "dense_pairs.tsv", list(rows[0]), rows)


def cmd_read(_a) -> None:
    reps = {_name(r): r for r in _reps()}
    pairs = read_tsv(OUT / "dense_pairs.tsv")
    sf = {n: r["superfamily"] for n, r in reps.items()}

    # (a) family-pair medians → D48's edge rule
    fam_tm = defaultdict(list)
    for p in pairs:
        if p["fam_a"] != p["fam_b"]:
            fam_tm[tuple(sorted((p["fam_a"], p["fam_b"])))].append(float(p["tm_avg"]))
    fpairs = [{"a": a, "b": b, "sf_a": CATALOGUE[a].superfamily,
               "sf_b": CATALOGUE[b].superfamily, "tm_avg": statistics.median(v)}
              for (a, b), v in sorted(fam_tm.items())]
    edges = edge_rows(fpairs)
    write_tsv(OUT / "dense_edges.tsv", list(edges[0]), edges)
    for e in edges:
        print(f"{e['a']:12s} {e['b']:12s} {e['median_tm']!s:7s} "
              f"{e['rank_b_for_a']}/{e['rank_a_for_b']} {e['verdict']}")

    # (b) superfamily recovery
    best_tm, best_e = {}, {}
    for p in pairs:
        if p["fam_a"] == p["fam_b"]:
            continue
        t = float(p["tm_avg"])
        e = float(p["foldseek_evalue"]) if p["foldseek_evalue"] else None
        for q, s in ((p["a"], p["b"]), (p["b"], p["a"])):
            if q not in best_tm or t > best_tm[q][0]:
                best_tm[q] = (t, s)
            if e is not None and (q not in best_e or e < best_e[q][0]):
                best_e[q] = (e, s)
    fams_in = defaultdict(set)
    for n, r in reps.items():
        fams_in[r["superfamily"]].add(r["family"])
    rec = []
    for n, r in sorted(reps.items()):
        if len(fams_in[r["superfamily"]]) < 2:
            continue
        t, s = best_tm[n]
        e, se = best_e.get(n, (None, ""))
        rec.append({"rep": n, "family": r["family"], "superfamily": r["superfamily"],
                    "group": r["group"], "best_tm_partner": s, "best_tm": round(t, 4),
                    "tm_recovered": int(sf[s] == r["superfamily"]),
                    "best_foldseek_partner": se, "best_evalue": e if e is not None else "",
                    "foldseek_recovered": int(bool(se) and sf[se] == r["superfamily"])})
    write_tsv(OUT / "dense_recovery.tsv", list(rec[0]), rec)

    # (c) detectability per superfamily pair
    acc = defaultdict(lambda: [0, 0, []])
    for p in pairs:
        if p["fam_a"] == p["fam_b"]:
            continue
        k = tuple(sorted((sf[p["a"]], sf[p["b"]])))
        acc[k][0] += 1
        acc[k][1] += int(bool(p["foldseek_evalue"]) and float(p["foldseek_evalue"]) <= E_DETECT)
        acc[k][2].append(float(p["tm_avg"]))
    det = [{"a": a, "b": b, "rep_pairs": n, "foldseek_detected": d,
            "frac_detected": round(d / n, 4), "median_tm": round(statistics.median(t), 4)}
           for (a, b), (n, d, t) in sorted(acc.items())]
    write_tsv(OUT / "dense_detect.tsv", list(det[0]), det)
    by = defaultdict(lambda: [0, 0, 0])
    for r in rec:
        by[r["superfamily"]][0] += 1
        by[r["superfamily"]][1] += r["tm_recovered"]
        by[r["superfamily"]][2] += r["foldseek_recovered"]
    summ = [{"superfamily": k, "families": len(fams_in[k]), "reps": n, "tm_recovered": t,
             "foldseek_recovered": f} for k, (n, t, f) in sorted(by.items())]
    write_tsv(OUT / "dense_recovery_by_superfamily.tsv", list(summ[0]), summ)
    for s in summ:
        print(s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("reps")
    p = sub.add_parser("measure")
    p.add_argument("--jobs", type=int, default=8)
    sub.add_parser("read")
    a = ap.parse_args()
    {"reps": cmd_reps, "measure": cmd_measure, "read": cmd_read}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
