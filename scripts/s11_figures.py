"""s11_figures.py — S11a's headline figure, drawn from the committed tables (D13).

    bin/envpy scripts/s11_figures.py

`results/duplication/figures/duplications.png` (+ .pdf):

A. Supported duplications older than one species, by the taxon they map to.
B. Supported species-specific duplications, by species (genome-locus tips marked).
C. OHNOLOGS v2 strict 2R pairs inside each family, by the window our trees date them to.
D. How each gene tree was rooted, and whether the declared root is a reconciliation optimum.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import figstyle as fs
from s3_hmm_lib import read_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "duplication"
TOP = 14
WIN = [("2R", fs.SUPERFAMILY["ploop"], "2R window (Vertebrata, Gnathostomata)"),
       ("older", fs.SUPERFAMILY["cysloop"], "older than vertebrates"),
       ("other", fs.GRID, "bony-vertebrate, younger or none")]


def main() -> int:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
    fs.use()
    dups = [r for r in read_tsv(OUT / "duplications.tsv") if r["supported"] == "True"]
    pairs = read_tsv(OUT / "human_pairs.tsv")
    trees = read_tsv(OUT / "recon_trees.tsv")
    fig = plt.figure(figsize=(fs.W_FULL, 6.6))
    gs = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.62,
                  left=0.15, right=0.98, top=0.9, bottom=0.08)

    ax = fig.add_subplot(gs[0, 0])
    deep = Counter(r["taxon"] for r in dups if r["species_specific"] == "False").most_common(TOP)
    for i, (t, n) in enumerate(deep):
        col = fs.SUPERFAMILY["ploop"] if t in ("Vertebrata", "Gnathostomata") else \
            fs.SUPERFAMILY["iglur"] if t == "Clupeocephala" else fs.BLUES[2]
        ax.barh(i, n, color=col, height=0.66, lw=0)
        ax.text(n + 0.8, i, str(n), va="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.set_yticks(range(len(deep)), [t for t, _ in deep], fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("supported duplications (all families)")
    fs.despine(ax)
    ax.legend(handles=[Patch(facecolor=fs.SUPERFAMILY["ploop"], label="2R window"),
                       Patch(facecolor=fs.SUPERFAMILY["iglur"], label="teleost 3R window"),
                       Patch(facecolor=fs.BLUES[2], label="other lineage")],
              fontsize=fs.FS_NOTE, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.panel(ax, "A", "Duplications shared by more than one species")

    bx = fig.add_subplot(gs[0, 1])
    ss = [r for r in dups if r["species_specific"] == "True"]
    top = Counter(r["taxon"] for r in ss).most_common(TOP)
    gen = Counter(r["taxon"] for r in ss if int(r["n_genome_tips"]) > 0)
    for i, (sp, n) in enumerate(top):
        g = gen[sp]
        bx.barh(i, n - g, color=fs.BLUES[3], height=0.66, lw=0)
        if g:
            bx.barh(i, g, left=n - g, color=fs.FAINT, height=0.66, lw=0)
        bx.text(n + 1.5, i, str(n), va="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
    bx.set_yticks(range(len(top)), [s for s, _ in top], fontsize=fs.FS_TICK, style="italic")
    bx.invert_yaxis()
    bx.set_xlabel("supported species-specific duplications")
    fs.despine(bx)
    bx.legend(handles=[Patch(facecolor=fs.BLUES[3], label="proteome genes"),
                       Patch(facecolor=fs.FAINT, label="involves a genome locus")],
              fontsize=fs.FS_NOTE, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.panel(bx, "B", "Duplications within one species")

    cx = fig.add_subplot(gs[1, 0])
    strict = [r for r in pairs if r["ohno_strict"] == "True"]
    by = Counter((r["family"], r["window"] if r["window"] in ("2R", "older") else "other")
                 for r in strict)
    fams = sorted({f for f, _ in by}, key=lambda f: -sum(by[(f, w)] for w, _, _ in WIN))
    fams = [f for f in fams if sum(by[(f, w)] for w, _, _ in WIN) >= 4]
    for i, f in enumerate(fams):
        left = 0
        for w, col, _ in WIN:
            n = by[(f, w)]
            if n:
                cx.barh(i, n, left=left, color=col, height=0.66, lw=0)
                left += n
    cx.set_yticks(range(len(fams)), fams, fontsize=fs.FS_TICK)
    cx.invert_yaxis()
    cx.set_xlabel("human 2R ohnologue pairs (OHNOLOGS v2, strict)")
    fs.despine(cx)
    cx.legend(handles=[Patch(facecolor=c, label=l) for _, c, l in WIN],
              fontsize=fs.FS_NOTE, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.panel(cx, "C", "Synteny-supported 2R pairs: where the gene trees date them")

    dx = fig.add_subplot(gs[1, 1])
    cats = [("declared root, a reconciliation optimum",
             lambda r: r["root_mode"] == "declared" and r["declared_optimal"] == "True"),
            ("declared root, not an optimum",
             lambda r: r["root_mode"] == "declared" and r["declared_optimal"] != "True"),
            ("reconciliation root, unique",
             lambda r: r["root_mode"] == "reconciliation" and r["n_optimal_roots"] == "1"),
            ("reconciliation root, tied",
             lambda r: r["root_mode"] == "reconciliation" and r["n_optimal_roots"] != "1")]
    cols = [fs.SUPERFAMILY["ploop"], fs.SUPERFAMILY["cysloop"], fs.BLUES[2], fs.GRID]
    for i, ((lab, sel), col) in enumerate(zip(cats, cols)):
        n = sum(sel(r) for r in trees)
        dx.barh(i, n, color=col, height=0.66, lw=0)
        dx.text(n + 0.5, i, str(n), va="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
    dx.set_yticks(range(len(cats)), [c for c, _ in cats], fontsize=fs.FS_TICK)
    dx.invert_yaxis()
    dx.set_xlabel("tier-1 gene trees")
    fs.despine(dx)
    fs.panel(dx, "D", "How the 63 gene trees were rooted")
    fig.suptitle("S11a — duplication history of the ion-channel families by reconciliation",
                 fontsize=fs.FS_SUPTITLE, y=0.99)
    for p in fs.save(fig, OUT / "figures" / "duplications"):
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
