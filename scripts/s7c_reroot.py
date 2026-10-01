"""s7c_reroot.py — S7c: re-root the weakly rooted P-loop families.

    python3 scripts/s7c_reroot.py align [--jobs 3]   # D39 set + basal picks → L-INS-i
    python3 scripts/s7c_reroot.py prep               # one input per family × root set
    python3 scripts/s7c_reroot.py run [--jobs 5 --threads 2]
    python3 scripts/s7c_reroot.py parse              # → results/phylogeny/tier1_reroot.tsv

**What changes from S7b, and what does not.** The families are those with a
per-family `ChannelFamily.root_with` (two declared outgroups each, S0-
verified by `s0_outgroups.py`). Their sequence sets are S6's D39 sets plus
`s7c_basal.py`'s early-lineage picks (D47) — re-aligned from scratch with
MAFFT L-INS-i, the S6 method. Everything after is S7's method unchanged
(D41): each outgroup added with `mafft --add --keeplength`, trimAl `-gt 0.5`
computed on the ingroup alignment alone and applied to the outgroup rows,
IQ-TREE 2 with `s7_trees.IQ_ARGS`. Each family gets **one tree per root
set**; the ingroup alignment and mask are identical in both, so the two
trees differ only in which outgroup is attached.

**The root criterion (fixed before any re-rooted tree, D47).** A root is
accepted only if, in both trees, (a) the outgroup is one clade and (b) its
edge carries UFBoot ≥ 95, and (c) the ingroup's basal split is identical
under the two outgroups. Otherwise the family's root is `unresolved`, and
no root is picked by how it looks.

Bulk: `<data root>/trees/s7c/`. Committed: `results/phylogeny/tier1_reroot*.tsv`
and the treefiles under `results/phylogeny/tier1/reroot/`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s0_lib import live_progress  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_tsv  # noqa: E402
from s6_lib import s6_dir  # noqa: E402
from s7_lib import (LIVE, OUT_DIR, ROOT, apply_mask, metrics,  # noqa: E402
                    read_aln, trim_columns, write_aln)
from s7_trees import IQ_ARGS, OG_PREFIX, TRIM, _fasta_text  # noqa: E402
from s7c_basal import reroot_families, s7c_dir  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
OG_PANEL = ROOT / "results" / "s0_baseline" / "outgroup_panel.fasta"
MIN_ROOT_UFBOOT = 95
INPUT_FIELDS = ["family", "root_set", "n_ingroup", "n_d39", "n_basal",
                "outgroup", "og_dropped", "cols", "informative", "gap_frac",
                "og_occupancy", "input_sha256"]


# --- align -----------------------------------------------------------------
def align_family(fam: str, threads: int) -> dict:
    d = s7c_dir("family")
    d39 = [(h.split()[0], s) for h, s in iter_fasta(s6_dir("family") / f"{fam}.fasta")]
    basal = [(h.split()[0], s) for h, s in iter_fasta(s7c_dir() / f"{fam}.basal.fasta")]
    members = d / f"{fam}.members.fasta"
    write_aln(members, d39 + basal)
    out = d / f"{fam}.linsi.fasta"
    rec = d / f"{fam}.align.json"
    sha = sha256(members)
    if rec.exists() and json.loads(rec.read_text()).get("input_sha256") == sha:
        return {"family": fam, "skipped": True}
    cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--thread",
           str(threads), "--quiet", str(members)]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft failed on {fam}: {p.stderr.strip()[:300]}")
    rows = [(h.split()[0], s.upper()) for h, s in _fasta_text(p.stdout)]
    if [l for l, _ in rows] != [l for l, _ in d39 + basal]:
        raise RuntimeError(f"{fam}: L-INS-i returned different rows")
    if len({len(s) for _, s in rows}) != 1:
        raise RuntimeError(f"{fam}: ragged alignment")
    for (l, a), (_, s) in zip(rows, d39 + basal):
        if a.replace("-", "") != s.upper():
            raise RuntimeError(f"{fam}: residues changed in {l}")
    write_aln(out, rows)
    r = {"family": fam, "cmd": " ".join(cmd), "input_sha256": sha,
         "n_d39": len(d39), "n_basal": len(basal),
         "seconds": round(time.time() - t0, 1), "aln_sha256": sha256(out)}
    rec.write_text(json.dumps(r, indent=1))
    return r


def cmd_align(a) -> None:
    fams = a.only.split(",") if a.only else reroot_families()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(align_family, f, a.threads): f for f in fams}
        for fu in as_completed(futs):
            r = fu.result()
            print(f"{r['family']:10s} "
                  f"{'skipped' if r.get('skipped') else str(r['seconds']) + ' s'}",
                  flush=True)


# --- prep ------------------------------------------------------------------
def outgroup_rows(fam: str, rs, ingroup: list) -> tuple[list, list]:
    """The root set's sequences; one identical to an ingroup row is dropped."""
    inseqs = {s.replace("-", "").upper() for _, s in ingroup}
    want_fam = set(rs.families)
    want_ex = {e.label for e in rs.exemplars}
    got = []
    for h, s in iter_fasta(PANEL):
        lab, f, _ = h.split()[0].split("|")
        if f in want_fam:
            got.append((lab, s))
    for h, s in iter_fasta(OG_PANEL):
        lab, rsn, _ = h.split()[0].split("|")
        if lab in want_ex:
            got.append((lab, s))
            want_ex.discard(lab)
    if want_ex:
        raise RuntimeError(f"{fam}/{rs.name}: {sorted(want_ex)} not in "
                           f"{OG_PANEL.name} — run s0_outgroups.py")
    keep = [(OG_PREFIX + l, s) for l, s in got if s.upper() not in inseqs]
    dropped = [l for l, s in got if s.upper() in inseqs]
    return keep, dropped


def prep_one(fam: str, rs) -> dict:
    d = s7c_dir("family")
    base = d / f"{fam}.linsi.fasta"
    rec = json.loads((d / f"{fam}.align.json").read_text())
    ingroup = read_aln(base)
    og, dropped = outgroup_rows(fam, rs, ingroup)
    ogf = d / f"{fam}.{rs.name}.outgroup.fasta"
    write_aln(ogf, og)
    cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--add", str(ogf),
           "--keeplength", "--thread", "4", "--quiet", str(base)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft --add failed on {fam}/{rs.name}: {p.stderr[:300]}")
    rows = [(h.split()[0], s.upper()) for h, s in _fasta_text(p.stdout)]
    if rows[:len(ingroup)] != [(l, s.upper()) for l, s in ingroup]:
        raise RuntimeError(f"{fam}/{rs.name}: mafft --add moved an ingroup row")
    if len(rows) != len(ingroup) + len(og):
        raise RuntimeError(f"{fam}/{rs.name}: wrong row count after --add")
    cols = trim_columns(base, TRIM)
    trimmed = apply_mask(rows, cols)
    out = d / f"{fam}.{rs.name}.tier1.fasta"
    write_aln(out, trimmed)
    m = metrics(trimmed[:len(ingroup)], list(range(len(cols))))
    occ = [sum(ch != "-" for ch in s) / len(cols) for l, s in trimmed
           if l.startswith(OG_PREFIX)]
    return {"family": fam, "root_set": rs.name, "n_ingroup": len(ingroup),
            "n_d39": rec["n_d39"], "n_basal": rec["n_basal"],
            "outgroup": ",".join(l for l, _ in og), "og_dropped": ",".join(dropped),
            "cols": m["cols"], "informative": m["informative"],
            "gap_frac": m["gap_frac"],
            "og_occupancy": ",".join(f"{x:.2f}" for x in occ),
            "input_sha256": sha256(out)}


def cmd_prep(a) -> None:
    jobs = [(f, rs) for f in reroot_families() for rs in CATALOGUE[f].root_with]
    rows = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for r in ex.map(lambda j: prep_one(*j), jobs):
            rows.append(r)
            print(f"{r['family']:10s} {r['root_set']:24s} {r['cols']:>5} cols "
                  f"og={r['outgroup']} occ={r['og_occupancy']}", flush=True)
    write_tsv(OUT_DIR / "tier1_reroot_inputs.tsv", INPUT_FIELDS, rows)


# --- run -------------------------------------------------------------------
def tree_id(row: dict) -> str:
    return f"{row['family']}.{row['root_set']}"


def run_one(row: dict, threads: int) -> dict:
    tid = tree_id(row)
    d = s7c_dir("iqtree", tid)
    rec = d / "run.json"
    if rec.exists():
        r = json.loads(rec.read_text())
        if r.get("input_sha256") == row["input_sha256"] and (d / f"{tid}.treefile").exists():
            return {**r, "skipped": True}
    src = s7c_dir("family") / f"{tid}.tier1.fasta"
    cmd = [shutil.which("iqtree2"), "-s", str(src), *IQ_ARGS, "-T", str(threads),
           "--prefix", str(d / tid), "--quiet", "-o", row["outgroup"]]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"iqtree2 failed on {tid}: "
                           f"{(p.stderr or p.stdout).strip()[-400:]}")
    r = {"tree": tid, "cmd": " ".join(cmd), "seconds": round(time.time() - t0, 1),
         "input_sha256": row["input_sha256"]}
    rec.write_text(json.dumps(r, indent=1))
    return r


def cmd_run(a) -> None:
    rows = read_tsv(OUT_DIR / "tier1_reroot_inputs.tsv")
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["family"] in keep or tree_id(r) in keep]
    rows.sort(key=lambda r: -int(r["n_ingroup"]) ** 2 * int(r["cols"]))
    done = {tree_id(r): False for r in rows}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_one, r, a.threads): tree_id(r) for r in rows}
        for fu in as_completed(futs):
            r = fu.result()
            done[futs[fu]] = True
            print(f"{futs[fu]:34s} "
                  f"{'skipped' if r.get('skipped') else str(r['seconds']) + ' s'}",
                  flush=True)
            live_progress(LIVE, "S7c", list(done.items()), workers=a.jobs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("align")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--threads", type=int, default=3)
    p = sub.add_parser("prep")
    p.add_argument("--jobs", type=int, default=4)
    p = sub.add_parser("run")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=5)
    p.add_argument("--threads", type=int, default=2)
    sub.add_parser("parse")
    a = ap.parse_args()
    if a.cmd == "parse":
        from s7c_parse import cmd_parse
        cmd_parse(a)
    else:
        {"align": cmd_align, "prep": cmd_prep, "run": cmd_run}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
