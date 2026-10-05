"""s11b_roots.py — S11b (2): outgroup-free roots by a non-reversible model (D53).

    python3 scripts/s11b_roots.py prep                       # ingroup inputs + models
    python3 scripts/s11b_roots.py run [--jobs 5 --threads 2]  # NQ.pfam rootstrap + root test
    python3 scripts/s11b_roots.py parse                      # → roots_nonrev.tsv

For S7d's six families and the five positive controls: S7b's tier-1 input
with its outgroup rows removed (D41 mask unchanged — the tips S11a
reconciled). Model NQ.pfam with the rate-heterogeneity terms of S7b's BIC
model; `-B 1000` gives the rootstrap, `--root-test -zb 1000 -au` the AU
confidence set of root branches. The reversible check re-evaluates the NQ
ML tree under Q.pfam with the same terms (`-te`, same parameter count).

A root is resolved iff its rootstrap ≥ 95 and ΔlnL(NQ − Q.pfam) > 0; it
agrees with S11a iff its split is one of S11a's optimal reconciliation
edges on S7b's tree. Controls are recovered iff the NQ root split equals the
declared S7b root split. Rules fixed in D53 before any input was built.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import s11_lib as L  # noqa: E402
from s11_recon import Reconciler, utree  # noqa: E402
from s11b_lib import (CONTROL_FAMILIES, OUT, ROOT_FAMILIES, ROOTSTRAP_MIN,  # noqa: E402
                      iqtree, live, s11b_dir)
from s3_hmm_lib import read_tsv, sha256, write_tsv  # noqa: E402
from s7_lib import read_aln, s7_dir, write_aln  # noqa: E402
from s7_newick import parse, root_partition  # noqa: E402

FAMILIES = ROOT_FAMILIES + CONTROL_FAMILIES
AU_ALPHA = 0.05


def nonrev_model(bic: str) -> str:
    """S7b's BIC model → NQ.pfam with the same +I/+G/+R terms, no +F (D53 (5))."""
    terms = [t for t in bic.split("+")[1:] if not t.startswith("F")]
    return "+".join(["NQ.pfam", *terms])


def s7b_model(fam: str) -> str:
    text = (s7_dir("iqtree", fam) / f"{fam}.iqtree").read_text()
    return re.search(r"Best-fit model according to BIC:\s*(\S+)", text).group(1)


def cmd_prep(_a) -> None:
    rows = []
    for fam in FAMILIES:
        rows_in = read_aln(s7_dir("family") / f"{fam}.tier1.fasta")
        ingroup = [(l, s) for l, s in rows_in if not l.startswith(L.OG_PREFIX)]
        out = s11b_dir("roots", fam) / f"{fam}.ingroup.fasta"
        write_aln(out, ingroup)
        bic = s7b_model(fam)
        rows.append({"family": fam, "role": "root" if fam in ROOT_FAMILIES else "control",
                     "n_tips": len(ingroup), "cols": len(ingroup[0][1]),
                     "dropped_outgroup": len(rows_in) - len(ingroup), "s7b_model": bic,
                     "nq_model": nonrev_model(bic),
                     "rev_model": nonrev_model(bic).replace("NQ.pfam", "Q.pfam"),
                     "input_sha256": sha256(out)})
        print(f"{fam:22s} {len(ingroup):4d} tips  {bic:20s} → {rows[-1]['nq_model']}")
    write_tsv(OUT / "roots_inputs.tsv", list(rows[0]), rows)


def run_one(row: dict, threads: int) -> str:
    fam = row["family"]
    d = s11b_dir("roots", fam)
    src = d / f"{fam}.ingroup.fasta"
    iqtree(src, d / "nq", ["-m", row["nq_model"], "-B", "1000", "--root-test",
                           "-zb", "1000", "-au", "-seed", "1"], threads)
    iqtree(src, d / "rev", ["-m", row["rev_model"], "-te", str(d / "nq.treefile"),
                            "-seed", "1"], threads)
    return fam


def cmd_run(a) -> None:
    rows = read_tsv(OUT / "roots_inputs.tsv")
    if a.only:
        rows = [r for r in rows if r["family"] in a.only.split(",")]
    rows.sort(key=lambda r: -int(r["n_tips"]) ** 2 * int(r["cols"]))
    steps = [{"label": f"NQ.pfam rootstrap {r['family']}", "done": False} for r in rows]
    live(steps, a.jobs)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, r, a.threads) for r in rows]
        for fu in as_completed(futs):
            fam = fu.result()
            for st in steps:
                st["done"] |= st["label"].endswith(f" {fam}")
            live(steps, a.jobs)
            print(f"done {fam}", flush=True)


# ------------------------------------------------------------------ parse

_ANN = re.compile(r'([^(),:\[\];]*):([-0-9.eE]+)\[&id="(\d+)",rootstrap="([0-9.]+)"\]')


def rootstrap_tree(path: Path):
    """The rootstrap Nexus tree → (parsed tree with internal names '#id',
    {branch id: (leaf set below, rootstrap)})."""
    text = re.search(r"tree\s+\S+\s*=\s*(.*?);", path.read_text(), re.S).group(1)
    rs: dict[str, float] = {}

    def sub(m):
        name, length, bid, val = m.groups()
        rs[bid] = float(val)
        return f"{name or '#' + bid}:{length}" if name else f"#{bid}:{length}"
    tree = parse(_ANN.sub(sub, text) + ";")
    below: dict[str, tuple[frozenset, float]] = {}

    def walk(n) -> frozenset:
        s = frozenset([n.name]) if not n.children else frozenset().union(*map(walk, n.children))
        return s

    leaf_id = {m.group(1): m.group(3) for m in _ANN.finditer(text) if m.group(1)}
    stack = [tree]
    while stack:
        n = stack.pop()
        stack.extend(n.children)
        bid = n.name[1:] if n.name.startswith("#") else leaf_id.get(n.name)
        if bid is not None:
            below[bid] = (walk(n), rs[bid])
    return tree, below


def canon(side: frozenset, leaves: frozenset) -> frozenset:
    """One orientation for an unrooted split: the side without the first leaf."""
    first = min(leaves)
    return leaves - side if first in side else side


def roottest(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def lnl(iq: Path) -> float:
    return float(re.search(r"Log-likelihood of the tree:\s*(-?[\d.]+)", iq.read_text()).group(1))


def recon_splits(tree, st, by_code) -> tuple[set[frozenset], frozenset]:
    """Canonical splits of every optimal reconciliation root (D52 cost)."""
    leaves = tree.leaves()
    og = {x for x in leaves if x.startswith(L.OG_PREFIX)}
    ing = frozenset(x for x in leaves if x not in og)
    rc = Reconciler(utree(tree, drop=og), st, {x: L.species_of_tip(x, by_code) for x in ing})
    opt, _ = rc.optimal()
    return {canon(frozenset(rc.names(rc.sub(*e).leaves)), ing) for e in opt}, ing


def parse_family(row: dict, st, by_code) -> dict:
    fam = row["family"]
    d = s11b_dir("roots", fam)
    tree, below = rootstrap_tree(d / "nq.rootstrap.nex")
    leaves = frozenset(tree.leaves())
    rt = roottest(d / "nq.roottest.csv")
    best = min(rt, key=lambda r: float(r["deltaL"]))
    side, boot = below[best["ID"]]
    ml_split = canon(side, leaves)
    au_set = {canon(below[r["ID"]][0], leaves) for r in rt if float(r["p-AU"]) >= AU_ALPHA}
    d_lnl = lnl(d / "nq.iqtree") - lnl(d / "rev.iqtree")
    s7b = parse((L.TREES / f"{fam}.treefile").read_text())
    opt_s7b, _ = recon_splits(s7b, st, by_code)
    opt_nq, _ = recon_splits(parse((d / "nq.treefile").read_text()), st, by_code)
    small = min((ml_split, leaves - ml_split), key=len)
    species = sorted({L.species_of_tip(x, by_code) for x in small})
    out = {"family": fam, "role": row["role"], "n_tips": len(leaves), "model": row["nq_model"],
           "dlnl_nonrev": round(d_lnl, 2), "root_rootstrap": boot,
           "root_split": f"{len(small)}|{len(leaves) - len(small)}",
           "small_side_species": len(species),
           "small_side": "; ".join(species[:4]) + (" …" if len(species) > 4 else ""),
           "au_set": len(au_set), "branches_tested": len(rt),
           "resolved": boot >= ROOTSTRAP_MIN and d_lnl > 0,
           "recon_optima_s7b": len(opt_s7b), "agrees_s11a": ml_split in opt_s7b,
           "recon_optimum_nq_topology": ml_split in opt_nq,
           "s11a_root_in_au_set": bool(opt_s7b & au_set), "control_recovered": "",
           "ml_root_tip": next(iter(small)) if len(small) == 1 else "",
           "au_set_frac": round(len(au_set) / len(rt), 2),
           "max_rootstrap": max(v for s, v in below.values() if len(s) < len(leaves))}
    if row["role"] == "control":
        og = {x for x in s7b.leaves() if x.startswith(L.OG_PREFIX)}
        rp = root_partition(s7b, og)
        declared = canon(frozenset(rp[0][0]), leaves) if rp else None
        out["control_recovered"] = declared == ml_split
        out["declared_in_au_set"] = declared in au_set
    return out


FIELDS = ["family", "role", "n_tips", "model", "dlnl_nonrev", "root_rootstrap", "root_split",
          "small_side_species", "small_side", "au_set", "branches_tested", "resolved",
          "recon_optima_s7b", "agrees_s11a", "recon_optimum_nq_topology",
          "s11a_root_in_au_set", "control_recovered", "declared_in_au_set", "ml_root_tip",
          "au_set_frac", "max_rootstrap"]


def cmd_parse(_a) -> None:
    _, st = L.species_tree()
    by_code = {L.code(r["species"]): r["species"] for r in L.panel()}
    rows = []
    for row in read_tsv(OUT / "roots_inputs.tsv"):
        if not (s11b_dir("roots", row["family"]) / "rev.iqtree").exists():
            print(f"{row['family']}: not finished — skipped")
            continue
        r = parse_family(row, st, by_code)
        rows.append(r)
        print(" ".join(f"{k}={r.get(k, '')}" for k in FIELDS[:1] + FIELDS[4:]))
    write_tsv(OUT / "roots_nonrev.tsv", FIELDS, rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prep")
    sub.add_parser("parse")
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=5)
    r.add_argument("--threads", type=int, default=2)
    r.add_argument("--only", default="")
    a = ap.parse_args()
    {"prep": cmd_prep, "run": cmd_run, "parse": cmd_parse}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
