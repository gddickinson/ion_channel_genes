"""s7_trees.py — S7: one rooted ML tree per census family (tier 1).

    python3 scripts/s7_trees.py prep                 # outgroups + D41 trim → tier1_inputs.tsv
    python3 scripts/s7_trees.py run [--jobs 5 --threads 2] [--only nav,cav]
    python3 scripts/s7_trees.py parse                # → tier1_trees.tsv + treefiles

**Rooting (D41), from the catalogue only.** A family is rooted on the
catalogue exemplars of its superfamily's declared outgroup family
(`Superfamily.root_with`). The outgroup family itself, and every family whose
superfamily declares no outgroup, is built unrooted and reported so — never
midpoint-rooted, never rooted on an analyst's pick.

**Alignment.** The outgroup exemplars are added to S6's family L-INS-i
alignment with `mafft --localpair --maxiterate 1000 --add … --keeplength`,
so no family column moves (checked: every ingroup row comes back
unchanged). **Trimming (D41)**: trimAl `-gt 0.5` computed on the family
alignment alone, the same column mask then applied to the outgroup rows.

**Trees.** IQ-TREE 2, ModelFinder (`-m MFP -mset LG,WAG,JTT,Q.pfam`), 1000 UFBoot with `-bnni`,
`-seed 1`; one command for every family. A run is skipped when its input's
SHA-256 matches the one recorded with the finished tree.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s0_lib import live_progress  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_tsv  # noqa: E402
from s7_lib import (LIVE, OUT_DIR, ROOT, TREE_DIR, apply_mask,  # noqa: E402
                    family_alignment, metrics, read_aln, s7_dir,
                    trim_columns, write_aln)
from s7_newick import ingroup_support, parse, split_support  # noqa: E402
from src.catalogue import CATALOGUE, SUPERFAMILIES  # noqa: E402

TRIM = "gt0.5"                      # D41
OG_PREFIX = "OG_"
#: ModelFinder over four general empirical matrices, every family alike (D41):
#: the full set (~500 models) cost 20–40 s per model on the largest families.
IQ_ARGS = ["-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam", "-B", "1000", "-bnni",
           "-seed", "1"]
PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
INPUT_FIELDS = ["family", "superfamily", "n_ingroup", "root_rule",
                "outgroup_family", "outgroup", "og_dropped", "cols",
                "informative", "gap_frac", "og_occupancy", "input_sha256"]


def aligned_families() -> list[str]:
    return [r["family"] for r in read_tsv(ROOT / "results" / "alignments" /
                                          "alignments.tsv")
            if r["status"] == "aligned"]


def root_rule(fam: str) -> tuple[str, str]:
    """(rule, outgroup family) — the catalogue's declaration, nothing else."""
    sf = SUPERFAMILIES[CATALOGUE[fam].superfamily]
    if not sf.root_with:
        return "unrooted:no_declared_outgroup", ""
    if fam in sf.root_with:
        return "unrooted:is_superfamily_outgroup", ""
    return "rooted", sf.root_with[0]


def outgroup_seqs(og_family: str, ingroup: list[tuple[str, str]]) -> tuple[list, list]:
    """Catalogue exemplars of `og_family` from the S0 reference panel.

    An exemplar whose sequence is already an ingroup row is dropped and
    listed — it cannot be both the root and a member.
    """
    inseqs = {s.replace("-", "").upper() for _, s in ingroup}
    keep, dropped = [], []
    for h, s in iter_fasta(PANEL):
        label, fam, _acc = h.split()[0].split("|")
        if fam != og_family:
            continue
        (dropped if s.upper() in inseqs else keep).append((OG_PREFIX + label, s))
    return keep, [l for l, _ in dropped]


def _add_outgroup(fam: str, og: list, d: Path) -> list[tuple[str, str]]:
    base = family_alignment(fam)
    ogf = d / f"{fam}.outgroup.fasta"
    write_aln(ogf, og)
    cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--add", str(ogf),
           "--keeplength", "--thread", "4", "--quiet", str(base)]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft --add failed on {fam}: {p.stderr.strip()[:300]}")
    rows = [(h.split()[0], s.upper()) for h, s in _fasta_text(p.stdout)]
    before = [(l, s.upper()) for l, s in read_aln(base)]
    if rows[:len(before)] != before:
        raise RuntimeError(f"{fam}: mafft --add moved an ingroup row")
    if len(rows) != len(before) + len(og):
        raise RuntimeError(f"{fam}: {len(rows)} rows after --add, expected "
                           f"{len(before) + len(og)}")
    (d / f"{fam}.add.json").write_text(json.dumps({"cmd": " ".join(cmd)}))
    return rows


def _fasta_text(text: str):
    head, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if head is not None:
                yield head, "".join(buf)
            head, buf = line[1:], []
        else:
            buf.append(line.strip())
    if head is not None:
        yield head, "".join(buf)


def prep_family(fam: str) -> dict:
    d = s7_dir("family")
    base = family_alignment(fam)
    ingroup = read_aln(base)
    rule, ogfam = root_rule(fam)
    og, dropped = (outgroup_seqs(ogfam, ingroup) if ogfam else ([], []))
    if ogfam and not og:
        rule = "unrooted:outgroup_not_available"
    rows = _add_outgroup(fam, og, d) if og else [(l, s.upper()) for l, s in ingroup]
    cols = trim_columns(base, TRIM)
    trimmed = apply_mask(rows, cols)
    out = d / f"{fam}.tier1.fasta"
    write_aln(out, trimmed)
    m = metrics(trimmed[:len(ingroup)], list(range(len(cols))))
    occ = [sum(ch != "-" for ch in s) / len(cols) for l, s in trimmed
           if l.startswith(OG_PREFIX)]
    return {"family": fam, "superfamily": CATALOGUE[fam].superfamily,
            "n_ingroup": len(ingroup), "root_rule": rule,
            "outgroup_family": ogfam if og else "",
            "outgroup": ",".join(l for l, _ in og), "og_dropped": ",".join(dropped),
            "cols": m["cols"], "informative": m["informative"],
            "gap_frac": m["gap_frac"],
            "og_occupancy": ",".join(f"{x:.2f}" for x in occ),
            "input_sha256": sha256(out)}


def cmd_prep(a) -> None:
    fams = a.only.split(",") if a.only else aligned_families()
    rows, prev = [], {}
    path = OUT_DIR / "tier1_inputs.tsv"
    if a.only and path.exists():
        prev = {r["family"]: r for r in read_tsv(path)}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(prep_family, f): f for f in fams}
        for fu in as_completed(futs):
            r = fu.result()
            prev[r["family"]] = r
            print(f"{r['family']:20s} {r['root_rule']:34s} {r['cols']:>5} cols "
                  f"og={r['outgroup'] or '-'} occ={r['og_occupancy'] or '-'}",
                  flush=True)
    rows = [prev[f] for f in sorted(prev)]
    write_tsv(path, INPUT_FIELDS, rows)


def _iq_dir(fam: str) -> Path:
    return s7_dir("iqtree", fam)


def run_family(row: dict, threads: int) -> dict:
    fam = row["family"]
    d = _iq_dir(fam)
    rec = d / "run.json"
    if rec.exists():
        r = json.loads(rec.read_text())
        if r.get("input_sha256") == row["input_sha256"] and (d / f"{fam}.treefile").exists():
            return {**r, "skipped": True}
    src = s7_dir("family") / f"{fam}.tier1.fasta"
    cmd = [shutil.which("iqtree2"), "-s", str(src), *IQ_ARGS, "-T", str(threads),
           "--prefix", str(d / fam), "--quiet"]  # resumes from .ckp.gz
    if row["outgroup"]:
        cmd += ["-o", row["outgroup"]]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"iqtree2 failed on {fam}: "
                           f"{(p.stderr or p.stdout).strip()[-400:]}")
    r = {"family": fam, "cmd": " ".join(cmd), "seconds": round(time.time() - t0, 1),
         "input_sha256": row["input_sha256"]}
    rec.write_text(json.dumps(r, indent=1))
    return r


def cmd_run(a) -> None:
    rows = read_tsv(OUT_DIR / "tier1_inputs.tsv")
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["family"] in keep]
    rows.sort(key=lambda r: -int(r["n_ingroup"]) ** 2 * int(r["cols"]))
    done: dict[str, bool] = {r["family"]: False for r in rows}
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(run_family, r, a.threads): r["family"] for r in rows}
        for fu in as_completed(futs):
            fam = futs[fu]
            r = fu.result()
            done[fam] = True
            print(f"{fam:20s} {'skipped' if r.get('skipped') else str(r['seconds']) + ' s'}",
                  flush=True)
            live_progress(LIVE, "S7", [(f, v) for f, v in done.items()],
                          workers=a.jobs)


def _iqtree_facts(fam: str) -> dict:
    txt = (_iq_dir(fam) / f"{fam}.iqtree").read_text()

    def grab(pat: str) -> str:
        m = re.search(pat, txt)
        return m.group(1).strip() if m else ""
    return {"model": grab(r"Best-fit model according to BIC:\s*(\S+)"),
            "log_likelihood": grab(r"Log-likelihood of the tree:\s*(-?[\d.]+)"),
            "parsimony_informative": grab(r"Number of parsimony informative sites:\s*(\d+)")}


TREE_FIELDS = ["family", "superfamily", "n_ingroup", "cols", "informative",
               "model", "log_likelihood", "root_rule", "outgroup_family",
               "n_outgroup", "outgroup_monophyletic", "ingroup_monophyletic",
               "root_ufboot", "internal_edges", "ufboot_median", "frac_ge95",
               "frac_lt70", "seconds"]


def cmd_parse(_a) -> None:
    TREE_DIR.mkdir(parents=True, exist_ok=True)
    out = []
    for row in read_tsv(OUT_DIR / "tier1_inputs.tsv"):
        fam = row["family"]
        d = _iq_dir(fam)
        tf = d / f"{fam}.treefile"
        if not tf.exists():
            print(f"{fam}: no tree yet", flush=True)
            continue
        run = json.loads((d / "run.json").read_text())
        if run["input_sha256"] != row["input_sha256"]:
            raise RuntimeError(f"{fam}: tree built on a different input")
        tree = parse(tf.read_text())
        og = {l for l in row["outgroup"].split(",") if l}
        mono, sup = split_support(tree, og) if og else (None, None)
        s = ingroup_support(tree, og)
        shutil.copy(tf, TREE_DIR / f"{fam}.treefile")
        out.append({"family": fam, "superfamily": row["superfamily"],
                    "n_ingroup": row["n_ingroup"], "cols": row["cols"],
                    "informative": row["informative"], **_iqtree_facts(fam),
                    "root_rule": row["root_rule"],
                    "outgroup_family": row["outgroup_family"],
                    "n_outgroup": len(og),
                    # one split separates outgroup from family: both are clades
                    "outgroup_monophyletic": "" if not og else mono,
                    "ingroup_monophyletic": "" if not og else mono,
                    "root_ufboot": "" if sup is None else sup,
                    **{k: ("" if v is None else v) for k, v in s.items()},
                    "seconds": run["seconds"]})
    for r in out:
        r.pop("parsimony_informative", None)
    write_tsv(OUT_DIR / "tier1_trees.tsv", TREE_FIELDS, out)
    print(f"{len(out)} trees → {OUT_DIR / 'tier1_trees.tsv'}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prep")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=3)
    p = sub.add_parser("run")
    p.add_argument("--only", default="")
    p.add_argument("--jobs", type=int, default=5)
    p.add_argument("--threads", type=int, default=2)
    sub.add_parser("parse")
    a = ap.parse_args()
    {"prep": cmd_prep, "run": cmd_run, "parse": cmd_parse}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
