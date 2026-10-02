"""s8_fold_network.py — S8 tier 3: the fold network, measured (D27, D48).

    bin/envpy scripts/s8_fold_network.py structures   # one AFDB model per census family
    bin/envpy scripts/s8_fold_network.py units        # cut each to its comparison unit
    bin/envpy scripts/s8_fold_network.py measure      # TM-align + Foldseek all-vs-all
    bin/envpy scripts/s8_fold_network.py edges        # LITERATURE_EDGES → verdicts

Rules fixed before any structure was compared (**D48 (5)**):

* **Node = one census family**, represented by the first catalogue exemplar
  (S0-resolved order) with an AlphaFold DB model of exactly that accession;
  none → the node is *unmeasured* with the reason (never a PDB substitute
  picked after the fact).
* **Comparison unit**: for a superfamily with a `module_rule`, the family's
  first module cut by S6's method (`hmmalign` to the family profile, cut at
  its span); otherwise the whole model. Residues with pLDDT < 70 are dropped
  everywhere. One extra node, `kv_shaker:VSD` (S1–S4 by the exemplar's
  UniProt TRANSMEM features), exists because the Hv1 edge asserts a
  voltage-sensor homology, not a pore one.
* **Metric**: TM-align TM-score normalised by the average length
  (`-a T`); Foldseek (exhaustive, 3Di+AA) E-value reported beside it.
* **Edge verdict** for a literature edge A–B (superfamily level): the
  median TM-score over family pairs a ∈ A, b ∈ B (a ≠ b). `supported` if
  that median is ≥ 0.5 **and** B ranks first among all other superfamilies
  by A's median (and A first by B's); `not_distinguished` otherwise;
  `unmeasured` if a side has no structure. An edge is structural similarity,
  never phylogeny: no branch length, no support value.
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests  # noqa: E402

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s6_lib import s3_profiles  # noqa: E402
from s6_project import a2m_map, cut, load_spans  # noqa: E402
from s8_lib import OUT_DIR, ROOT, unit_families  # noqa: E402
from src.catalogue import CATALOGUE, SUPERFAMILIES  # noqa: E402
from src.phylo.network import LITERATURE_EDGES, FoldEdge, seed_network  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

MIN_PLDDT = 70.0
TM_SUPPORT = 0.5
EXEMPLARS = ROOT / "results" / "s0_baseline" / "exemplars_resolved.tsv"
NET_DIR = OUT_DIR / "fold_network"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
VSD_NODE = ("kv_shaker:VSD", "kv_shaker")


def sdir(*p: str) -> Path:
    """`<data root>/structures/s8/...` — raises without the drive (D1)."""
    d = require_data_root().joinpath("structures", "s8", *p)
    d.mkdir(parents=True, exist_ok=True)
    return d


def census_superfamilies() -> list[str]:
    return sorted(k for k in SUPERFAMILIES if unit_families(k))


def node_superfamily(node: str) -> str:
    return CATALOGUE[node.split(":")[0]].superfamily


# ------------------------------------------------------------ structures

def _afdb(acc: str) -> dict | None:
    f = sdir("raw") / f"afdb_{acc}.json"
    if not f.exists():
        r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60)
        if r.status_code in (400, 404, 422):
            f.write_text("[]")
        else:
            r.raise_for_status()
            f.write_text(r.text)
    hits = [e for e in json.loads(f.read_text()) if e.get("uniprotAccession") == acc]
    return hits[0] if hits else None


def cmd_structures(_a) -> None:
    rows = []
    ex = [r for r in read_tsv(EXEMPLARS) if r["resolved_accession"]]
    for sf in census_superfamilies():
        for fam in unit_families(sf):
            row = {"node": fam, "superfamily": sf, "exemplar": "", "accession": "",
                   "species": "", "model": "", "length": "", "note": ""}
            tried = []
            for r in (r for r in ex if r["family"] == fam):
                acc = r["resolved_accession"]
                tried.append(acc)
                e = _afdb(acc)
                if not e:
                    continue
                pdb = sdir("afdb") / Path(e["pdbUrl"]).name
                if not pdb.exists():
                    rr = requests.get(e["pdbUrl"], timeout=120)
                    rr.raise_for_status()
                    pdb.write_bytes(rr.content)
                row.update(exemplar=r["label"], accession=acc, species=r["species"],
                           model=pdb.name, length=e.get("uniprotEnd", ""))
                break
            else:
                row["note"] = ("no AlphaFold DB model for any exemplar: "
                               + ",".join(tried) if tried else "no resolved exemplar")
            rows.append(row)
            print(f"{fam:22s} {row['accession'] or '-':12s} {row['note']}", flush=True)
    write_tsv(NET_DIR / "structures.tsv", list(rows[0]), rows)


# ------------------------------------------------------------ units

def read_ca(pdb: Path) -> list[tuple[int, str, float]]:
    out = []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append((int(line[22:26]), AA3.get(line[17:20], "X"), float(line[60:66])))
    return out


def module_range(fam: str, seq: str) -> tuple[int, int] | None:
    d = sdir("hmm")
    fa, out = d / f"{fam}.fasta", d / f"{fam}.a2m"
    fa.write_text(f">q\n{seq}\n")
    p = subprocess.run(["hmmalign", "--amino", "--outformat", "A2M", "-o", str(out),
                        str(s3_profiles() / f"{fam}.hmm"), str(fa)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"hmmalign {fam}: {p.stderr.strip()[:300]}")
    a2m = "".join(l.strip() for l in out.read_text().splitlines()[1:])
    st, res = a2m_map(a2m)
    k0, k1, _ = load_spans()[fam][0]
    a, b, cover = cut(st, res, seq, k0, k1)
    return (a, b) if cover >= 0.5 else None


def vsd_range(acc: str) -> tuple[int, int]:
    f = sdir("raw") / f"uniprot_{acc}.json"
    if not f.exists():
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60)
        r.raise_for_status()
        f.write_text(r.text)
    tms = sorted((x["location"]["start"]["value"], x["location"]["end"]["value"])
                 for x in json.loads(f.read_text())["features"]
                 if x["type"] == "Transmembrane")
    return tms[0][0], tms[3][1]


def write_unit(pdb: Path, keep: set[int], out: Path) -> int:
    lines = [l for l in pdb.read_text().splitlines()
             if l.startswith("ATOM") and int(l[22:26]) in keep]
    out.write_text("\n".join(lines) + "\nEND\n")
    return len({int(l[22:26]) for l in lines})


def cmd_units(_a) -> None:
    rows = []
    structs = read_tsv(NET_DIR / "structures.tsv")
    todo = [(r["node"], r) for r in structs]
    todo += [(VSD_NODE[0], next(r for r in structs if r["node"] == VSD_NODE[1]))]
    for node, s in todo:
        fam = node.split(":")[0]
        row = {"node": node, "superfamily": s["superfamily"], "accession": s["accession"],
               "unit": "", "start": "", "end": "", "residues": "", "mean_plddt": "",
               "note": s["note"]}
        if not s["model"]:
            rows.append(row)
            continue
        ca = read_ca(sdir("afdb") / s["model"])
        seq = "".join(a for _, a, _ in ca)
        if node == VSD_NODE[0]:
            rng, unit = vsd_range(s["accession"]), "VSD S1-S4 (UniProt TRANSMEM)"
        elif SUPERFAMILIES[s["superfamily"]].module_rule:
            rng, unit = module_range(fam, seq), "module 1 (S6 span)"
            if rng is None:
                row["note"] = "module not found on the exemplar (< 50 % of span)"
                rows.append(row)
                continue
            rng = (ca[rng[0] - 1][0], ca[rng[1] - 1][0])
        else:
            rng, unit = (ca[0][0], ca[-1][0]), "whole model"
        keep = {n for n, _, pl in ca if rng[0] <= n <= rng[1] and pl >= MIN_PLDDT}
        pl = [p for n, _, p in ca if n in keep]
        out = sdir("units") / f"{node.replace(':', '_')}.pdb"
        n = write_unit(sdir("afdb") / s["model"], keep, out)
        row.update(unit=unit, start=rng[0], end=rng[1], residues=n,
                   mean_plddt=round(statistics.mean(pl), 1) if pl else "")
        if n < 30:
            row["note"] = f"only {n} residues at pLDDT ≥ {MIN_PLDDT:g}"
        rows.append(row)
        print(f"{node:22s} {unit:30s} {n:5d} res", flush=True)
    write_tsv(NET_DIR / "units.tsv", list(rows[0]), rows)


# ------------------------------------------------------------ measure

def _tmalign(a: Path, b: Path) -> tuple[float, float, float]:
    p = subprocess.run(["TMalign", str(a), str(b), "-a", "T"], capture_output=True,
                       text=True)
    tms = [float(x) for x in re.findall(r"TM-score=\s*([\d.]+)", p.stdout)]
    if len(tms) < 3:
        raise RuntimeError(f"TMalign {a.name} {b.name}: {p.stdout[-300:]}")
    return tms[0], tms[1], tms[2]   # by chain 1, by chain 2, by average


def cmd_measure(a) -> None:
    units = [r for r in read_tsv(NET_DIR / "units.tsv")
             if r["residues"] and int(r["residues"]) >= 30]
    path = {r["node"]: sdir("units") / f"{r['node'].replace(':', '_')}.pdb" for r in units}
    pairs = list(itertools.combinations(sorted(path), 2))
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(lambda p: _tmalign(path[p[0]], path[p[1]]), pairs))
    fs = _foldseek(sorted(path.values()))
    rows = []
    for (x, y), (t1, t2, ta) in zip(pairs, res):
        f = fs.get((x, y)) or fs.get((y, x)) or {}
        rows.append({"a": x, "b": y, "sf_a": node_superfamily(x),
                     "sf_b": node_superfamily(y), "tm_avg": ta,
                     "tm_short": max(t1, t2), "tm_long": min(t1, t2),
                     "foldseek_evalue": f.get("evalue", ""),
                     "foldseek_prob": f.get("prob", "")})
    write_tsv(NET_DIR / "pairs.tsv", list(rows[0]), rows)
    print(f"{len(rows)} pairs over {len(path)} nodes")


def _foldseek(pdbs: list[Path]) -> dict:
    d = sdir("foldseek")
    out, tmp = d / "allvsall.m8", d / "tmp"
    subprocess.run(["foldseek", "easy-search", str(sdir("units")), str(sdir("units")),
                    str(out), str(tmp), "--exhaustive-search", "1", "-e", "10",
                    "--format-output", "query,target,evalue,prob", "-v", "1"],
                   check=True, capture_output=True, text=True)
    best: dict = {}
    for line in out.read_text().splitlines():
        q, t, e, pr = line.split("\t")
        q, t = (Path(q).stem.replace("_VSD", ":VSD"), Path(t).stem.replace("_VSD", ":VSD"))
        if q == t:
            continue
        cur = best.get((q, t))
        if cur is None or float(e) < float(cur["evalue"]):
            best[(q, t)] = {"evalue": e, "prob": pr}
    return best


# ------------------------------------------------------------ edges

def sf_medians(pairs: list[dict]) -> dict[tuple[str, str], float]:
    """Median TM-score per unordered superfamily pair (VSD node kept apart)."""
    acc: dict = {}
    for p in pairs:
        a = "ploop:VSD" if p["a"] == VSD_NODE[0] else p["sf_a"]
        b = "ploop:VSD" if p["b"] == VSD_NODE[0] else p["sf_b"]
        acc.setdefault(tuple(sorted((a, b))), []).append(float(p["tm_avg"]))
    return {k: statistics.median(v) for k, v in acc.items()}


def cmd_edges(_a) -> None:
    pairs = read_tsv(NET_DIR / "pairs.tsv")
    med = sf_medians(pairs)
    nodes = {x for k in med for x in k}
    rows = []
    for a, b, claim in LITERATURE_EDGES:
        if a == "hv":
            b = "ploop:VSD"
        key = tuple(sorted((a, b)))
        if key not in med:
            rows.append({"a": a, "b": b, "claim": claim, "median_tm": "", "n_pairs": 0,
                         "rank_b_for_a": "", "rank_a_for_b": "", "best_other": "",
                         "best_other_tm": "", "verdict": "unmeasured"})
            continue

        def rank(x: str, y: str) -> int:
            """y's rank among x's other superfamilies (x itself too for a self-edge)."""
            others = sorted(((med[tuple(sorted((x, z)))], z) for z in nodes
                             if (z != x or y == x) and tuple(sorted((x, z))) in med),
                            reverse=True)
            return [z for _, z in others].index(y) + 1
        n = sum(1 for p in pairs if tuple(sorted((
            "ploop:VSD" if p["a"] == VSD_NODE[0] else p["sf_a"],
            "ploop:VSD" if p["b"] == VSD_NODE[0] else p["sf_b"]))) == key)
        ra, rb = rank(a, b), rank(b, a)
        bg = sorted(((v, k) for k, v in med.items()
                     if (a in k or b in k) and k != key and len(set(k)) == 2), reverse=True)
        ok = med[key] >= TM_SUPPORT and ra == 1 and rb == 1
        rows.append({"a": a, "b": b, "claim": claim, "median_tm": round(med[key], 4),
                     "n_pairs": n, "rank_b_for_a": ra, "rank_a_for_b": rb,
                     "best_other": "-".join(bg[0][1]) if bg else "",
                     "best_other_tm": round(bg[0][0], 4) if bg else "",
                     "verdict": "supported" if ok else "not_distinguished"})
        print(f"{a:14s} {b:12s} median TM {med[key]:.3f}  ranks {ra}/{rb}  "
              f"→ {rows[-1]['verdict']}", flush=True)
    write_tsv(NET_DIR / "literature_edges.tsv", list(rows[0]), rows)
    write_network(rows, pairs)
    write_tsv(NET_DIR / "superfamily_medians.tsv", ["a", "b", "median_tm"],
              [{"a": k[0], "b": k[1], "median_tm": round(v, 4)}
               for k, v in sorted(med.items())])


def write_network(edges: list[dict], pairs: list[dict]) -> None:
    """The tier-3 object: superfamily nodes, the literature edges measured."""
    net = seed_network("superfamily")
    for e in edges:
        if e["verdict"] == "unmeasured":
            net.note_unmeasured(e["a"], e["b"], "a side has no structure")
            continue
        net.add(FoldEdge(e["a"], e["b"], "tm_score", float(e["median_tm"]),
                         evidence=(f"median TM-align TM-score (avg length) over "
                                   f"{e['n_pairs']} family pairs; AFDB v6 models; "
                                   f"verdict {e['verdict']} (D48)"), measured=True))
    measured = {n for p in pairs for n in (p["sf_a"], p["sf_b"])}
    for n in net.nodes:
        if n not in measured:
            net.note_unmeasured(n, "*", "no family of this superfamily has a structure")
    net.write(NET_DIR / "network.json")
    print(net.summary())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("structures")
    sub.add_parser("units")
    p = sub.add_parser("measure")
    p.add_argument("--jobs", type=int, default=6)
    sub.add_parser("edges")
    a = ap.parse_args()
    NET_DIR.mkdir(parents=True, exist_ok=True)
    {"structures": cmd_structures, "units": cmd_units, "measure": cmd_measure,
     "edges": cmd_edges}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
