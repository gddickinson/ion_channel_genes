"""s13_jobs.py — D55 (7)–(8) as a job list: every codeml and HyPhy run S13 makes.

    python3 scripts/s13_jobs.py          # write inputs + jobs.tsv, run nothing

Site sets (`site_*`): M0, M1a, M2a, M7, M8, M8a, pairwise dS (`runmode -2`),
M0 without *Petromyzon* where the set holds it, HyPhy FEL.
Family contrast sets (`family_*`): M0 (branch dS for the stem gate), the
two-ratio null (one ω for every branch *within* an orthogroup, stems
background) and alternative (the anchor orthogroup's within-branches
separate), branch-site model A on the anchor orthogroup's stem (null ω₂ = 1;
alternative from each of `BS_OMEGA_STARTS`), HyPhy RELAX (test = anchor
within-branches, reference = the other orthogroups'; `--models Minimal`).
Orthogroups (`og_*`): M0 on the family alignment's rows for that orthogroup.

Labelled trees are drawn from `results/selection/trees/<set>.nwk` with
`s13_tree.paml_newick`, which raises if a labelled set is not a clade.
Inputs and job directories under `<data root>/selection/s13/jobs/`.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s3_hmm_lib import read_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402
from s13_codeml import Job  # noqa: E402
from s13_codon import write_phylip  # noqa: E402
from s13_lib import (ANCHORS, BS_OMEGA_STARTS, CYCLOSTOME, OUT,  # noqa: E402
                     anchor_label, s13_dir, species_code)
from s13_tree import UTree, drawn_for, paml_newick, prune  # noqa: E402

SITE_MODELS = {"m0": {}, "m1a": {"NSsites": 1}, "m2a": {"NSsites": 2},
               "m7": {"NSsites": 7}, "m8": {"NSsites": 8},
               "m8a": {"NSsites": 8, "fix_omega": 1, "omega": 1}}
COST = {"m0": 1, "pair": 1, "m1a": 2, "m2a": 3, "m7": 4, "m8": 5, "m8a": 5,
        "fel": 2, "relax": 6, "two": 2, "bs": 8}


def _codes() -> tuple[dict[str, str], dict[str, str]]:
    rows = read_tsv(OUT / "tip_codes.tsv")
    return ({r["label"]: r["code"] for r in rows}, {r["code"]: r["label"] for r in rows})


def _write_inputs(d: Path, rows: list[tuple[str, str]], nwk: str) -> tuple[Path, Path, Path]:
    d.mkdir(parents=True, exist_ok=True)
    write_phylip(d / "seq.phy", rows)
    write_fasta(d / "seq.fasta", rows)
    (d / "tree.nwk").write_text(nwk + "\n")
    return d / "seq.phy", d / "seq.fasta", d / "tree.nwk"


def _coded(set_name: str, lab2code: dict) -> list[tuple[str, str]]:
    return [(lab2code[k], s) for k, s in read_fasta(OUT / "codon" / f"{set_name}.fasta").items()]


def site_jobs(set_name: str, lab2code: dict, code2lab: dict) -> list[Job]:
    inp = s13_dir("jobs", "_inputs", set_name)
    rows = _coded(set_name, lab2code)
    nwk = (OUT / "trees" / f"{set_name}.nwk").read_text().strip()
    phy, fa, tree = _write_inputs(inp, rows, nwk)
    size = len(rows) * len(rows[0][1]) / 3
    jobs = [Job(f"{m}_{set_name}", "codeml", phy, tree, s13_dir("jobs", f"{m}_{set_name}"),
                dict(s), set_name, size * COST[m]) for m, s in SITE_MODELS.items()]
    jobs.append(Job(f"pair_{set_name}", "codeml", phy, tree, s13_dir("jobs", f"pair_{set_name}"),
                    {"runmode": -2}, set_name, size * COST["pair"]))
    jobs.append(Job(f"fel_{set_name}", "fel", fa, tree, s13_dir("jobs", f"fel_{set_name}"),
                    {"pvalue": 0.1, "ci": "No"}, set_name, size * COST["fel"]))
    cyc = [c for c, _ in rows if species_code(code2lab[c]) in CYCLOSTOME]
    if cyc:
        keep = [(c, s) for c, s in rows if c not in cyc]
        t = prune(UTree.from_newick(nwk), {c for c, _ in keep})
        sub = paml_newick(drawn_for(t, [frozenset(c for c, _ in keep)]), {})
        p2, _, t2 = _write_inputs(s13_dir("jobs", "_inputs", f"{set_name}_nocyc"), keep, sub)
        jobs.append(Job(f"m0_{set_name}_nocyc", "codeml", p2, t2,
                        s13_dir("jobs", f"m0_{set_name}_nocyc"), {}, set_name, size))
    return jobs


def family_jobs(fam: str, lab2code: dict) -> list[Job]:
    set_name = f"family_{fam}"
    sets = [r for r in read_tsv(OUT / "sets.tsv") if r["set"] == set_name]
    rows = _coded(set_name, lab2code)
    present = {c for c, _ in rows}
    og: dict[str, set[str]] = {}
    for r in sets:
        if lab2code[r["label"]] in present:
            og.setdefault(r["orthogroup"], set()).add(lab2code[r["label"]])
    anchor_og = next(k for k, v in og.items() if lab2code[anchor_label(fam)] in v)
    groups = {k: frozenset(v) for k, v in og.items()}
    t = UTree.from_newick((OUT / "trees" / f"{set_name}.nwk").read_text())
    drawn = drawn_for(t, list(groups.values()))
    multi = {k: g for k, g in groups.items() if len(g) > 1}   # a lone tip has no within-branch
    others = {g: "#1" for k, g in multi.items() if k != anchor_og}
    a = groups[anchor_og]
    trees = {
        "plain": paml_newick(drawn, {}),
        "two_null": paml_newick(drawn, {g: "#1" for g in multi.values()}),
        "two_alt": paml_newick(drawn, {**others, a: "#2"}),
        "bs": paml_newick(drawn, {}, stem={a: "#1"}),
        "relax": paml_newick(drawn, {**{g: "Reference" for g in others}, a: "Test"},
                             style="hyphy"),
    }
    base = s13_dir("jobs", "_inputs", set_name)
    phy, fa, plain = _write_inputs(base, rows, trees["plain"])
    paths = {}
    for k, nwk in trees.items():
        paths[k] = base / f"tree_{k}.nwk"
        paths[k].write_text(nwk + "\n")
    size = len(rows) * len(rows[0][1]) / 3

    def cj(name: str, tree: str, settings: dict, cost: str) -> Job:
        return Job(f"{name}_{set_name}", "codeml", phy, paths[tree],
                   s13_dir("jobs", f"{name}_{set_name}"), settings, set_name, size * COST[cost])

    jobs = [cj("m0", "plain", {}, "m0"),
            cj("two_null", "two_null", {"model": 2}, "two"),
            cj("two_alt", "two_alt", {"model": 2}, "two"),
            cj("bs_null", "bs", {"model": 2, "NSsites": 2, "fix_omega": 1, "omega": 1}, "bs")]
    jobs += [cj(f"bs_alt_w{w}", "bs", {"model": 2, "NSsites": 2, "fix_omega": 0, "omega": w}, "bs")
             for w in BS_OMEGA_STARTS]
    jobs.append(Job(f"relax_{set_name}", "relax", fa, paths["relax"],
                    s13_dir("jobs", f"relax_{set_name}"),
                    {"test": "Test", "reference": "Reference", "models": "Minimal"},
                    set_name, size * COST["relax"]))
    full = dict(rows)
    for k, g in groups.items():          # per-orthogroup M0 on the family columns
        sub = [(c, full[c]) for c in sorted(g)]
        if len(sub) < 2:
            continue
        nwk = (OUT / "trees" / f"og_{k}.nwk").read_text().strip()
        p, _, tr = _write_inputs(s13_dir("jobs", "_inputs", f"og_{k}"), sub, nwk)
        jobs.append(Job(f"m0_og_{k}", "codeml", p, tr, s13_dir("jobs", f"m0_og_{k}"),
                        {}, f"og_{k}", len(sub) * len(sub[0][1]) / 3))
    return jobs


def all_jobs() -> list[Job]:
    lab2code, code2lab = _codes()
    jobs: list[Job] = []
    for set_name in sorted({r["set"] for r in read_tsv(OUT / "sets.tsv") if r["kind"] == "site"}):
        jobs += site_jobs(set_name, lab2code, code2lab)
    for fam in ANCHORS:
        if (OUT / "codon" / f"family_{fam}.fasta").exists():
            jobs += family_jobs(fam, lab2code)
    return sorted(jobs, key=lambda j: -j.cost)


if __name__ == "__main__":
    js = all_jobs()
    write_tsv(OUT / "jobs.tsv", ["name", "tool", "group", "cost", "settings"],
              [{"name": j.name, "tool": j.tool, "group": j.group, "cost": int(j.cost),
                "settings": ";".join(f"{k}={v}" for k, v in j.settings.items())} for j in js])
    print(f"{len(js)} jobs; {sum(j.done() for j in js)} already done")
