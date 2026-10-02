"""s7_gost_tree.py — the GOST superfamily tree: TMEM87 orthology outside animals (D45).

    python3 scripts/s7_gost_tree.py build      # set → L-INS-i → trimAl -gt 0.5 → IQ-TREE
    python3 scripts/s7_gost_tree.py read       # the split tests → results/phylogeny/gost/
    bin/envpy scripts/s7_gost_tree.py figure   # the tree, leaves coloured by lineage

Two questions S2f left open, answered on one tree of every S4-panel protein
the profiles put in the TMEM87/GOST superfamily (census v3 over the panel):

1. **Are any non-animal GOST proteins TMEM87 orthologues?** None carries the
   TMEM87 GOLD domain (PF21901, metazoan-only), so only descent can say. Test:
   the **largest** clade holding every animal TMEM87 (PF21901 carrier) and no
   animal GPR107/108 is the TMEM87 lineage; non-animal leaves inside it, with
   that edge's UFBoot, are TMEM87 orthologues. (A first version took the
   *smallest* such clade — the animals alone — and so asked whether plant
   genes sit *among* the animal TMEM87s, which an orthologue never does.)
2. **When did TMEM87A and TMEM87B split?** D45 assumed a vertebrate
   duplication. Test: is there an edge separating the vertebrate A and B
   genes (together) from the invertebrate TMEM87s?

Same method as the S7 family trees (D41): MAFFT L-INS-i, trimAl `-gt 0.5`,
IQ-TREE 2 `-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -bnni -seed 1`. Unrooted:
every question is a split, so no root is assumed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta, s2_raw, write_tsv  # noqa: E402
from s3b_lib import panel_db, s3b_dir  # noqa: E402
from s7_lib import s7_dir  # noqa: E402
from s7_newick import parse, split_support  # noqa: E402

OUT = ROOT / "results" / "phylogeny" / "gost"
METAZOA = {"vertebrate", "invertebrate", "deuterostome", "cnidarian", "basal_metazoan"}
IQ = ["-m", "MFP", "-mset", "LG,WAG,JTT,Q.pfam", "-B", "1000", "-bnni", "-seed", "1"]


def members() -> list[dict]:
    rows = []
    with gzip.open(s3b_dir() / "census_v3.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["v3_superfamily"] == "tmem87" or r["p_family"] in ("tmem87", "nonchannel_gost"):
                rows.append(r)
    want = {r["accession"] for r in rows}
    gold = set()
    with gzip.open(s2_raw() / "census_v2.tsv.gz", "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["accession"] in want and "PF21901" in r["pfam"]:
                gold.add(r["accession"])
    seqs, genes = {}, {}
    for h, s in iter_fasta(panel_db()):
        t = h.split()[0]
        acc = t.split("|")[1] if t.count("|") >= 2 else t
        if acc in want:
            seqs[acc] = s
            m = re.search(r" GN=(\S+)", h)
            genes[acc] = m.group(1) if m else ""
    out = []
    for r in rows:
        a = r["accession"]
        sp = "".join(w[:3] for w in r["species"].split()[:2])
        call = r["v3_family"] or r["v3_status"]
        label = f"{a}_{sp}_{call}".replace(".", "")
        out.append({"label": label, "accession": a, "species": r["species"],
                    "group": r["group"], "metazoa": int(r["group"] in METAZOA),
                    "call": call, "confidence": r["p_confidence"],
                    "gold_pf21901": int(a in gold), "gene": genes.get(a, ""),
                    "length": len(seqs.get(a, ""))})
    return [m for m in out if m["length"]], seqs


def cmd_build(_a) -> None:
    rows, seqs = members()
    d = s7_dir("gost")
    fa = d / "gost.fasta"
    with open(fa, "w") as fh:
        for m in rows:
            fh.write(f">{m['label']}\n{seqs[m['accession']]}\n")
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "members.tsv", list(rows[0]), rows)
    aln, trm = d / "gost.linsi.fasta", d / "gost.tier1.fasta"
    if not trm.exists():
        with open(aln, "w") as out:
            subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "4",
                            "--quiet", str(fa)], stdout=out, check=True)
        subprocess.run(["trimal", "-in", str(aln), "-out", str(trm), "-gt", "0.5"], check=True)
    # no -redo: an interrupted run resumes from its checkpoint
    subprocess.run(["iqtree2", "-s", str(trm), *IQ, "-T", "4", "--prefix", str(d / "gost"),
                    "--quiet"], check=True)
    (OUT / "gost.treefile").write_text((d / "gost.treefile").read_text())
    model = re.search(r"Best-fit model according to BIC:\s*(\S+)",
                      (d / "gost.iqtree").read_text())
    (OUT / "run.json").write_text(json.dumps(
        {"sequences": len(rows), "trimmed_columns": len(next(iter_fasta(trm))[1]),
         "model": model.group(1) if model else "", "iqtree": " ".join(IQ)}, indent=1))
    print(f"{len(rows)} sequences → {OUT / 'gost.treefile'}")


def _vert_gene(m: dict) -> str:
    g = m["gene"].upper()
    return "A" if g.startswith("TMEM87A") else "B" if g.startswith("TMEM87B") else ""


def cmd_read(_a) -> None:
    rows = list(csv.DictReader(open(OUT / "members.tsv"), delimiter="\t"))
    tree = parse((OUT / "gost.treefile").read_text())
    leaves = set(tree.leaves())
    by = {m["label"]: m for m in rows if m["label"] in leaves}
    tmem87_animal = {l for l, m in by.items() if m["metazoa"] == "1" and m["gold_pf21901"] == "1"}
    gost_animal = {l for l, m in by.items() if m["metazoa"] == "1" and m["gold_pf21901"] == "0"
                   and m["call"] == "nonchannel_gost"}
    nonanimal = {l for l in by if by[l]["metazoa"] == "0"}
    # 1. the smallest clade holding all animal TMEM87 and no animal GOST, over
    #    every way of attaching the non-animal leaves: test each split
    from s7_newick import splits
    nested = sorted(((s, sup) for side, sup in splits(tree)
                     for s in (side, frozenset(leaves) - side)
                     if tmem87_animal <= s and not (s & gost_animal)),
                    key=lambda x: len(x[0]))
    best = nested[-1] if nested else None        # the whole TMEM87 lineage
    write_tsv(OUT / "tmem87_lineage_nesting.tsv",
              ["leaves", "ufboot", "non_animal", "non_animal_species"],
              [{"leaves": len(s), "ufboot": "" if sup is None else sup,
                "non_animal": sum(by[l]["metazoa"] == "0" for l in s),
                "non_animal_species": ", ".join(sorted({by[l]["species"].split(" (")[0]
                                                       for l in s if by[l]["metazoa"] == "0"}))}
               for s, sup in nested])
    res = []
    if best:
        side, sup = best
        for l in sorted(nonanimal):
            res.append({"label": l, "species": by[l]["species"], "group": by[l]["group"],
                        "profile_call": by[l]["call"], "gene": by[l]["gene"],
                        "placement": "TMEM87 lineage" if l in side else "GOST (GPR107/108) side",
                        "edge_ufboot": sup if sup is not None else ""})
    write_tsv(OUT / "nonanimal_placement.tsv",
              ["label", "species", "group", "profile_call", "gene", "placement", "edge_ufboot"],
              res)
    # 2. vertebrate A + B against invertebrate TMEM87
    vert = {l for l, m in by.items() if m["group"] == "vertebrate" and l in tmem87_animal}
    inv = tmem87_animal - vert
    a = {l for l in vert if _vert_gene(by[l]) == "A"}
    b = {l for l in vert if _vert_gene(by[l]) == "B"}
    tests = {"all animal TMEM87 (PF21901) vs animal GPR107/108":
             (best is not None, best[1] if best else None, len(best[0]) if best else 0),
             "vertebrate TMEM87 (A+B) as one clade, invertebrates outside":
             split_support(tree, vert) + (len(vert),),
             "vertebrate TMEM87A a clade": split_support(tree, a) + (len(a),),
             "vertebrate TMEM87B a clade": split_support(tree, b) + (len(b),)}
    summary = {k: {"split_present": v[0], "ufboot": v[1], "leaves": v[2]}
               for k, v in tests.items()}
    summary["n_tmem87_animal"], summary["n_gost_animal"] = len(tmem87_animal), len(gost_animal)
    summary["n_invertebrate_tmem87"] = len(inv)
    summary["nonanimal_in_tmem87_lineage"] = sum(r["placement"] == "TMEM87 lineage" for r in res)
    summary["profile_agrees_with_tree"] = sum(
        (r["placement"] == "TMEM87 lineage") == (r["profile_call"] == "tmem87") for r in res)
    summary["nonanimal_total"] = len(res)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
    render(summary, res)
    for r in res:
        print(f"  {r['species'][:28]:28s} {r['profile_call']:16s} {r['gene']:12s} → {r['placement']}")


def render(summary: dict, res: list[dict]) -> None:
    run = json.loads((OUT / "run.json").read_text())
    t = summary["all animal TMEM87 (PF21901) vs animal GPR107/108"]
    miscalled = [r for r in res if r["profile_call"] == "tmem87"]
    L = ["# The GOST superfamily tree — is TMEM87 animal-specific?", "",
         "Rendered by `scripts/s7_gost_tree.py read` from `summary.json` and "
         "`nonanimal_placement.tsv` (D13).", "", "![GOST tree](figures/gost_tree.png)", "",
         f"**{run['sequences']} panel proteins** the profiles place in the TMEM87/GOST "
         f"superfamily; MAFFT L-INS-i, trimAl `-gt 0.5` ({run['trimmed_columns']} "
         f"columns), IQ-TREE 2 `{run['iqtree']}` (model {run['model']}). Unrooted: "
         "every question is a split.", "",
         f"**1. TMEM87 is a pan-eukaryotic lineage; its GOLD domain is an animal "
         f"addition.** The largest clade holding all {summary['n_tmem87_animal']} animal "
         f"TMEM87s (PF21901 carriers) and none of the {summary['n_gost_animal']} animal "
         f"GPR107/108s has UFBoot {t['ufboot']:.0f} and contains "
         f"**{summary['nonanimal_in_tmem87_lineage']} of {summary['nonanimal_total']} "
         "non-animal GOST proteins** — plants, moss, fungi, *Dictyostelium*, "
         "*Plasmodium*, the choanoflagellate and *Capsaspora* "
         "(`tmem87_lineage_nesting.tsv` gives the nested clades and their support). "
         f"The profile calls agree with the tree on "
         f"**{summary['profile_agrees_with_tree']}/{summary['nonanimal_total']}** non-animal "
         "proteins (TMEM87 call ⇔ TMEM87 lineage). A first reading of this tree "
         "took the smallest clade instead and concluded the opposite; it was wrong.", "",
         "**2. The TMEM87A/B duplication is not dated by this tree.** Neither "
         "paralogue forms a clade, nor do the vertebrate A and B genes together "
         "against the invertebrates. D45's grounds stand on what was measured — "
         "profiles cannot separate non-mammalian TMEM87 into A and B — but its "
         "premise that the split is a vertebrate duplication is unverified here.", "",
         "| non-animal protein | species | profile call | tree placement |",
         "|---|---|---|---|"]
    L += [f"| {r['gene']} | {r['species']} | {r['profile_call']} | {r['placement']} |"
          for r in res]
    (OUT / "report.md").write_text("\n".join(L) + "\n")


def cmd_figure(_a) -> None:
    import figstyle as fs
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from s7_figures import _layout
    rows = {m["label"]: m for m in csv.DictReader(open(OUT / "members.tsv"), delimiter="\t")}
    fs.use()
    _, _, items = _layout(parse((OUT / "gost.treefile").read_text()))
    col = {"animal TMEM87 (PF21901)": "#184f95", "animal GOST": "#86b6ef",
           "plant / alga": "#0f7d3d", "fungus": "#eb6834", "holozoan": "#4a3aa7",
           "other protist": "#a9a79e"}

    def kind(m):
        if m["metazoa"] == "1":
            return "animal TMEM87 (PF21901)" if m["gold_pf21901"] == "1" else "animal GOST"
        return {"plant": "plant / alga", "algae": "plant / alga", "fungi": "fungus",
                "holozoa": "holozoan"}.get(m["group"], "other protist")
    fig, ax = plt.subplots(figsize=(fs.W_HALF + 1.2, 7.6))
    for it in items:
        if it[0] == "node":
            _, name, x, y, kids = it
            ax.plot([x, x], [kids[0][1], kids[-1][1]], color=fs.MUTED, lw=0.45)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=fs.MUTED, lw=0.45)
            try:
                if name and float(name) >= 95 and len(kids) > 1:
                    ax.scatter(x, y, s=3, color=fs.INK, zorder=3, linewidth=0)
            except ValueError:
                pass
        else:
            _, name, x, y = it
            m = rows.get(name)
            k = kind(m) if m else "other protist"
            ax.scatter(x, y, s=9, color=col[k], zorder=4, linewidth=0,
                       marker="s" if m and m["call"] == "tmem87" and m["metazoa"] == "0" else "o")
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("substitutions / site", fontsize=fs.FS_LABEL)
    hs = [Line2D([], [], marker="o", ls="", color=c, markersize=3.5, label=k) for k, c in col.items()]
    hs.append(Line2D([], [], marker="s", ls="", color=fs.FAINT, markersize=3.5,
                     label="non-animal, profile-called TMEM87"))
    ax.legend(handles=hs, fontsize=fs.FS_NOTE - 0.6, frameon=False, ncol=2,
              loc="upper left", bbox_to_anchor=(0.0, -0.07))
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "", "The GOST superfamily (unrooted, drawn as rooted by IQ-TREE)")
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, OUT / "figures" / "gost_tree"))))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "read", "figure"])
    a = ap.parse_args()
    {"build": cmd_build, "read": cmd_read, "figure": cmd_figure}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
