"""s8_tier2.py — S8: one tier-2 tree per alignable multi-family superfamily.

    python3 scripts/s8_tier2.py reps                 # D8 representatives → tier2_reps.tsv
    python3 scripts/s8_tier2.py prep [--only ploop]  # root tips, L-INS-i, D41 mask → tier2_inputs.tsv
    python3 scripts/s8_tier2.py run [--jobs 4 --threads 2]
    python3 scripts/s8_tier2.py parse                # S8b → tier2_trees/families/placement.tsv

Every rule here was fixed before any tier-2 tree was built (**D48**):

* **Units** (`s8_lib.units()`): alignable superfamilies with ≥ 2 census
  families. Pore-module units (a declared `module_rule`) take S6's `full`
  modules, one tip per module; the others take full-length members.
  Single-family superfamilies are their tier-1 tree; refused ones go to the
  fold network (D27) — both listed in `tier2_units.tsv`.
* **Representatives (D8)**: per family × panel group × module, central-first
  greedy clustering at ≥ `THRESHOLD` identity (`s8_lib.greedy`). One
  threshold for every unit, chosen on tip counts alone (`tier2_threshold.tsv`).
* **Root**: the superfamily's `root_with` family (D41). Its tips are the
  outgroup **restricted to the kingdom of its catalogue exemplars** — a
  prokaryotic root family's animal members are ingroup tips, their placement
  a result (the S7a `plgic_prok` row). The root family's catalogue exemplars
  are added as tips unless a tip already carries the identical sequence
  (module units: cut by S6's method — `hmmalign` to the family profile, cut
  at its module span). No declared root → unrooted, never midpoint.
* **Method**: MAFFT L-INS-i over the unit; trimAl `-gt 0.5` (D41) on the
  unit alignment; IQ-TREE 2 with S7's `IQ_ARGS`. A run is skipped when its
  input SHA-256 matches the finished one.
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
import requests  # noqa: E402

from s0_lib import live_progress  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_tsv  # noqa: E402
from s6_lib import s3_profiles, s6_dir  # noqa: E402
from s6_project import a2m_map, cut, load_spans  # noqa: E402
from s7_lib import apply_mask, metrics, read_aln, trim_columns, write_aln  # noqa: E402
from s7_trees import IQ_ARGS  # noqa: E402
from s8_lib import (LIVE, OUT_DIR, ROOT, base_label, cells, greedy,  # noqa: E402
                    is_module_unit, matrix, members, s8_dir, single_family_units,
                    unit_families, units)
from src.catalogue import CATALOGUE, SUPERFAMILIES  # noqa: E402

THRESHOLD = 0.5          # D48: one identity threshold for every unit ...
#: ... unless the unit then has more than MAX_TIPS_PER_SITE tips per informative
#: site of its trimmed alignment at 0.5 (measured: P-loop 840 tips on 85 sites —
#: the ~105-residue pore module is the limit, not the trim). Such a unit takes
#: the highest threshold in `tier2_threshold.tsv` that meets the cap; fixed on
#: alignment properties before any tier-2 tree was built.
MAX_TIPS_PER_SITE = 4
UNIT_THRESHOLD = {"ploop": 0.3}   # 317 tips ≤ 4 × 85 (0.4: 558 > 340)
TRIM = "gt0.5"           # D41
PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
EXEMPLARS = ROOT / "results" / "s0_baseline" / "exemplars_resolved.tsv"
EX_PREFIX = "EX_"
REP_FIELDS = ["unit", "family", "group", "module", "tip", "label", "cluster_size",
              "cell_size"]
INPUT_FIELDS = ["unit", "kind", "n_tips", "n_families", "root_family", "n_outgroup",
                "outgroup_excluded", "exemplars_added", "cols", "informative",
                "gap_frac", "input_sha256"]


def tip_name(fam: str, label: str) -> str:
    return f"{fam}__{label}"


# ------------------------------------------------------------------ reps

def cmd_reps(_a) -> None:
    rows, thr = [], []
    for sf in units():
        c = cells(sf)
        counts = {t: 0 for t in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)}
        for (fam, grp, mod), seqs in sorted(c.items()):
            m = matrix(seqs)
            for t in counts:
                counts[t] += len(greedy(list(seqs), m, t))
            for rep, mem in greedy(list(seqs), m, UNIT_THRESHOLD.get(sf, THRESHOLD)):
                rows.append({"unit": sf, "family": fam, "group": grp,
                             "module": mod.split(":")[0], "tip": tip_name(fam, rep),
                             "label": base_label(rep), "cluster_size": len(mem),
                             "cell_size": len(seqs)})
        thr.append({"unit": sf, "kind": "module" if is_module_unit(sf) else "full_length",
                    "sequences": sum(len(v) for v in c.values()), "cells": len(c),
                    "threshold": UNIT_THRESHOLD.get(sf, THRESHOLD),
                    **{f"t{int(t * 100)}": n for t, n in counts.items()}})
        print(f"{sf:14s} {sum(1 for r in rows if r['unit'] == sf):5d} tips "
              f"from {thr[-1]['sequences']} sequences", flush=True)
    write_tsv(OUT_DIR / "tier2_reps.tsv", REP_FIELDS, rows)
    write_tsv(OUT_DIR / "tier2_threshold.tsv", list(thr[0]), thr)
    units_tab = [{"superfamily": k, "tier2": "tree", "families": len(unit_families(k)),
                  "kind": "module" if is_module_unit(k) else "full_length",
                  "note": ""} for k in units()]
    treed = {r["family"] for r in read_tsv(OUT_DIR / "tier1_trees.tsv")}
    units_tab += [{"superfamily": k, "tier2": "= tier 1", "families": 1,
                   "kind": "module" if is_module_unit(k) else "full_length",
                   "note": ("one census family: its tier-1 tree is the tier-2 tree"
                            if unit_families(k)[0] in treed else
                            "one census family, and no tier-1 tree "
                            "(< 4 sequences, or added after S6 — r4)")}
                  for k in single_family_units()]
    units_tab += [{"superfamily": k, "tier2": "refused (D27)",
                   "families": len(unit_families(k)), "kind": "",
                   "note": SUPERFAMILIES[k].notes or SUPERFAMILIES[k].anchor_module}
                  for k, sf in sorted(SUPERFAMILIES.items()) if not sf.alignable]
    write_tsv(OUT_DIR / "tier2_units.tsv",
              ["superfamily", "tier2", "families", "kind", "note"], units_tab)


# ------------------------------------------------------------------ root

def _kingdom(accession: str) -> str:
    """'prokaryote' | 'eukaryote' | 'virus' from the entry's UniProt lineage (archived)."""
    f = s8_dir("raw") / f"entry_{accession}.json"
    if not f.exists():
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json",
                         timeout=60)
        r.raise_for_status()
        f.write_text(r.text)
    lin = set(json.loads(f.read_text())["organism"].get("lineage", []))
    if "Viruses" in lin:
        return "virus"
    return "eukaryote" if "Eukaryota" in lin else "prokaryote"


def root_family(sf: str) -> str:
    rw = SUPERFAMILIES[sf].root_with
    return rw[0] if rw else ""


def root_exemplars(fam: str) -> tuple[list[tuple[str, str]], set[str]]:
    """([(label, sequence)], kingdoms) of the root family's catalogue exemplars."""
    seqs = {h.split()[0].split("|")[0]: s for h, s in iter_fasta(PANEL)
            if h.split()[0].split("|")[1] == fam}
    kingdoms, out = set(), []
    for r in read_tsv(EXEMPLARS):
        if r["family"] == fam and r["label"] in seqs:
            kingdoms.add(_kingdom(r["resolved_accession"]))
            out.append((r["label"], seqs[r["label"]]))
    return out, kingdoms


def _species_kingdom() -> dict[str, str]:
    mem = members()
    return {k: ("prokaryote" if r["group"] == "prokaryote" else
                "virus" if r["group"] == "virus" else "eukaryote")
            for k, r in mem.items()}


def cut_modules(fam: str, seqs: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """S6's method on extra sequences: hmmalign to `fam`'s profile, cut at its spans."""
    d = s8_dir("exemplars")
    fa, out = d / f"{fam}.fasta", d / f"{fam}.a2m"
    write_aln(fa, seqs)
    p = subprocess.run(["hmmalign", "--amino", "--outformat", "A2M", "-o", str(out),
                        str(s3_profiles() / f"{fam}.hmm"), str(fa)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"hmmalign {fam}: {p.stderr.strip()[:300]}")
    raw = dict(seqs)
    spans = load_spans()[fam]
    res = []
    for h, a2m in iter_fasta(out):
        lab = h.split()[0]
        st, rmap = a2m_map(a2m)
        for i, (k0, k1, _b) in enumerate(spans, 1):
            a, b, cover = cut(st, rmap, raw[lab], k0, k1)
            if cover >= 0.5:
                res.append((lab if len(spans) == 1 else f"{lab}_m{i}", raw[lab][a - 1:b]))
    return res


# ------------------------------------------------------------------ prep

def unit_sequences(sf: str, reps: list[dict]) -> dict[str, str]:
    if is_module_unit(sf):
        src = {h.split()[0]: s for h, s in
               iter_fasta(s6_dir("modules") / f"{sf}.modules.fasta")}
        return {r["tip"]: src[f"{r['family']}|{r['tip'].split('__', 1)[1]}"]
                for r in reps}
    out = {}
    for fam in sorted({r["family"] for r in reps}):
        src = {h.split()[0]: s for h, s in iter_fasta(s6_dir("family") / f"{fam}.fasta")}
        out.update({r["tip"]: src[r["label"]] for r in reps if r["family"] == fam})
    return out


def prep_unit(sf: str) -> dict:
    reps = [r for r in read_tsv(OUT_DIR / "tier2_reps.tsv") if r["unit"] == sf]
    seqs = unit_sequences(sf, reps)
    rfam = root_family(sf)
    og, excluded, added = [], [], []
    if rfam:
        ex, kingdoms = root_exemplars(rfam)
        kin = _species_kingdom()
        for r in reps:
            if r["family"] != rfam:
                continue
            (og if kin[(rfam, r["label"])] in kingdoms else excluded).append(r["tip"])
        if is_module_unit(sf):
            ex = cut_modules(rfam, ex)
        have = {s.upper() for s in seqs.values()}
        for lab, s in ex:
            if s.upper() in have:
                continue
            t = tip_name(rfam, EX_PREFIX + lab)
            seqs[t] = s
            og.append(t)
            added.append(t)
    d = s8_dir("unit")
    raw = d / f"{sf}.fasta"
    write_aln(raw, sorted(seqs.items()))
    aln = d / f"{sf}.linsi.fasta"
    p = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "4",
                        "--quiet", str(raw)], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft L-INS-i failed on {sf}: {p.stderr.strip()[:300]}")
    aln.write_text(p.stdout)
    rows = [(l, s.upper()) for l, s in read_aln(aln)]
    if len(rows) != len(seqs) or len({len(s) for _, s in rows}) != 1:
        raise RuntimeError(f"{sf}: alignment has {len(rows)} rows / ragged")
    for l, s in rows:
        if s.replace("-", "") != seqs[l].upper():
            raise RuntimeError(f"{sf}: mafft changed the residues of {l}")
    cols = trim_columns(aln, TRIM)
    out = d / f"{sf}.tier2.fasta"
    write_aln(out, apply_mask(rows, cols))
    m = metrics(apply_mask(rows, cols), list(range(len(cols))))
    (d / f"{sf}.outgroup.json").write_text(json.dumps(
        {"root_family": rfam, "outgroup": og, "excluded": excluded}, indent=1))
    return {"unit": sf, "kind": "module" if is_module_unit(sf) else "full_length",
            "n_tips": len(rows), "n_families": len({r["family"] for r in reps}),
            "root_family": rfam, "n_outgroup": len(og),
            "outgroup_excluded": ",".join(excluded), "exemplars_added": ",".join(added),
            "cols": m["cols"], "informative": m["informative"], "gap_frac": m["gap_frac"],
            "input_sha256": sha256(out)}


def cmd_prep(a) -> None:
    todo = a.only.split(",") if a.only else units()
    path = OUT_DIR / "tier2_inputs.tsv"
    prev = {r["unit"]: r for r in read_tsv(path)} if path.exists() and a.only else {}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(prep_unit, sf): sf for sf in todo}
        for fu in as_completed(futs):
            r = fu.result()
            prev[r["unit"]] = r
            print(f"{r['unit']:14s} {r['n_tips']:5} tips {r['cols']:5} cols "
                  f"{r['informative']:5} inf  og={r['n_outgroup']} "
                  f"(excluded {r['outgroup_excluded'] or '-'})", flush=True)
    write_tsv(path, INPUT_FIELDS, [prev[k] for k in sorted(prev)])


# ------------------------------------------------------------------ run

def run_unit(row: dict, threads: int) -> dict:
    sf = row["unit"]
    d = s8_dir("iqtree", sf)
    rec = d / "run.json"
    if rec.exists():
        r = json.loads(rec.read_text())
        if r.get("input_sha256") == row["input_sha256"] and (d / f"{sf}.treefile").exists():
            return {**r, "skipped": True}
    src = s8_dir("unit") / f"{sf}.tier2.fasta"
    og = json.loads((s8_dir("unit") / f"{sf}.outgroup.json").read_text())["outgroup"]
    cmd = [shutil.which("iqtree2"), "-s", str(src), *IQ_ARGS, "-T", str(threads),
           "--prefix", str(d / sf), "--quiet"]
    if og:
        cmd += ["-o", ",".join(og)]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"iqtree2 failed on {sf}: {(p.stderr or p.stdout).strip()[-400:]}")
    r = {"unit": sf, "cmd": " ".join(cmd), "seconds": round(time.time() - t0, 1),
         "input_sha256": row["input_sha256"]}
    rec.write_text(json.dumps(r, indent=1))
    return r


def cmd_run(a) -> None:
    rows = read_tsv(OUT_DIR / "tier2_inputs.tsv")
    if a.only:
        rows = [r for r in rows if r["unit"] in set(a.only.split(","))]
    rows.sort(key=lambda r: -int(r["n_tips"]) ** 2 * int(r["cols"]))
    done = {r["unit"]: False for r in rows}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_unit, r, a.threads): r["unit"] for r in rows}
        for fu in as_completed(futs):
            r = fu.result()
            done[futs[fu]] = True
            print(f"{futs[fu]:14s} {'skipped' if r.get('skipped') else str(r['seconds']) + ' s'}",
                  flush=True)
            live_progress(LIVE, "S8", list(done.items()), workers=a.jobs)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("reps")
    p = sub.add_parser("prep")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=3)
    p = sub.add_parser("run")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--threads", type=int, default=2)
    sub.add_parser("parse")
    a = ap.parse_args()
    if a.cmd == "parse":
        from s8_parse import cmd_parse
        cmd_parse(a)
        return 0
    {"reps": cmd_reps, "prep": cmd_prep, "run": cmd_run}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
