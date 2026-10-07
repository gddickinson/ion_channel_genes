"""s12_experimental.py — S12 (2)–(3): the experimental reference per family and
the AlphaFold-model check (D54).

    bin/envpy scripts/s12_experimental.py pick     # → results/structures/experimental.tsv
    bin/envpy scripts/s12_experimental.py check    # + model-vs-experiment TM-scores

Per census family (S8a's nodes, RyR included): the catalogue exemplars in S0
order; for each, its comparison unit (D48: module 1 by S6's `hmmalign` + span,
else the whole sequence) computed on the UniProt sequence; PDBe SIFTS
`best_structures` chains in PDBe's rank order; the first chain whose
*observed* residues (SIFTS per-residue UniProt numbering, updated mmCIF)
cover ≥ 50 % of the unit is the family's experimental reference. The unit is
cut from it, CA only. The VSD node takes the Shaker reference cut to S1–S4.

`check`: where the reference's accession is S8a's node accession, TM-align
S8a's AFDB unit against it; TM-score normalised by the experimental unit;
`consistent` iff ≥ 0.5.
"""

from __future__ import annotations

import argparse
import re
import subprocess

from s8_fold_network import VSD_NODE, vsd_range
from s12_lib import (EXEMPLARS, OUT, S8_NET, UNIT_COVER, experimental_ca, has_module,
                     module_range, read_tsv, sdir, sifts_best, uniprot_json, write_ca,
                     write_tsv)

CONSISTENT, HIGH = 0.5, 0.8


def unit_range(fam: str, acc: str, seq: str) -> tuple[int, int] | None:
    if has_module(fam):
        return module_range(fam, seq, f"{fam}_{acc}")
    return 1, len(seq)


def pick_one(fam: str, exemplars: list[dict]) -> dict:
    row = {"node": fam, "exemplar": "", "accession": "", "pdb": "", "chain": "",
           "method": "", "resolution": "", "unit_start": "", "unit_end": "",
           "unit_len": "", "observed": "", "unit_cover": "", "chains_tried": 0, "note": ""}
    tried = []
    for ex in exemplars:
        acc = ex["resolved_accession"]
        seq = uniprot_json(acc)["sequence"]["value"]
        rng = unit_range(fam, acc, seq)
        if rng is None:
            tried.append(f"{acc}:no_module")
            continue
        a, b = rng
        n_unit = b - a + 1
        chains = [c for c in sifts_best(acc)
                  if (min(b, c["unp_end"]) - max(a, c["unp_start"]) + 1) / n_unit >= UNIT_COVER]
        tried.append(f"{acc}:{len(chains)}")
        for c in chains:
            row["chains_tried"] += 1
            ca = experimental_ca(c["pdb_id"], c["chain_id"], acc)
            keep = [x for x in ca if a <= x[0] <= b]
            if len(keep) / n_unit < UNIT_COVER:
                continue
            out = sdir("exp_units") / f"{fam}.pdb"
            write_ca([x[2] for x in keep], out)
            row.update(exemplar=ex["label"], accession=acc, pdb=c["pdb_id"],
                       chain=c["chain_id"], method=c.get("experimental_method", ""),
                       resolution=c.get("resolution") or "", unit_start=a, unit_end=b,
                       unit_len=n_unit, observed=len(keep),
                       unit_cover=round(len(keep) / n_unit, 3))
            return row
    row["note"] = "no experimental chain covering ≥ 50 % of the unit: " + ",".join(tried)
    return row


def cmd_pick(_a) -> None:
    ex = [r for r in read_tsv(EXEMPLARS) if r["resolved_accession"]]
    rows = []
    for s in read_tsv(S8_NET / "structures.tsv"):
        fam = s["node"]
        r = pick_one(fam, [e for e in ex if e["family"] == fam])
        r["superfamily"] = s["superfamily"]
        r["node_accession"] = s["accession"]
        rows.append(r)
        print(f"{fam:22s} {r['accession'] or '-':11s} {r['pdb']:5s} {r['chain']:3s} "
              f"{r['observed']!s:>5s}/{r['unit_len']!s:5s} {r['note'][:60]}", flush=True)
    vs = next(r for r in rows if r["node"] == VSD_NODE[1])
    v = dict(vs, node=VSD_NODE[0], note="")
    if vs["pdb"]:
        a, b = vsd_range(vs["accession"])
        ca = [x for x in experimental_ca(vs["pdb"], vs["chain"], vs["accession"])
              if a <= x[0] <= b]
        write_ca([x[2] for x in ca], sdir("exp_units") / "kv_shaker_VSD.pdb")
        v.update(unit_start=a, unit_end=b, unit_len=b - a + 1, observed=len(ca),
                 unit_cover=round(len(ca) / (b - a + 1), 3))
    rows.append(v)
    cols = ["node", "superfamily", "node_accession", "exemplar", "accession", "pdb",
            "chain", "method", "resolution", "unit_start", "unit_end", "unit_len",
            "observed", "unit_cover", "chains_tried", "note"]
    write_tsv(OUT / "experimental.tsv", cols, rows)


def tmalign(a, b) -> dict:
    p = subprocess.run(["TMalign", str(a), str(b)], capture_output=True, text=True)
    tms = [float(x) for x in re.findall(r"TM-score=\s*([\d.]+)", p.stdout)]
    m = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([\d.]+)", p.stdout)
    if len(tms) < 2 or not m:
        raise RuntimeError(f"TMalign {a} {b}: {p.stdout[-300:]}")
    return {"tm_by_model": tms[0], "tm_by_exp": tms[1], "aligned": int(m.group(1)),
            "rmsd": float(m.group(2))}


def cmd_check(_a) -> None:
    rows = read_tsv(OUT / "experimental.tsv")
    units = {r["node"]: r for r in read_tsv(S8_NET / "units.tsv")}
    s8u = sdir().parent / "s8" / "units"
    for r in rows:
        r.update(model_residues="", tm_by_exp="", tm_by_model="", aligned="", rmsd="",
                 model_check="")
        u = units.get(r["node"], {})
        if not r["pdb"]:
            r["model_check"] = "no_experimental"
        elif r["accession"] != r["node_accession"] or not u.get("residues"):
            r["model_check"] = "different_accession" if u.get("residues") else "no_model"
        else:
            t = tmalign(s8u / f"{r['node'].replace(':', '_')}.pdb",
                        sdir("exp_units") / f"{r['node'].replace(':', '_')}.pdb")
            r.update(model_residues=u["residues"], **{k: round(v, 4) if isinstance(v, float)
                                                     else v for k, v in t.items()})
            r["model_check"] = ("consistent_high" if t["tm_by_exp"] >= HIGH else
                                "consistent" if t["tm_by_exp"] >= CONSISTENT else
                                "inconsistent")
        print(f"{r['node']:22s} {r['model_check']:20s} {r['tm_by_exp']}", flush=True)
    write_tsv(OUT / "experimental.tsv", list(rows[0]), rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("pick")
    sub.add_parser("check")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"pick": cmd_pick, "check": cmd_check}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
