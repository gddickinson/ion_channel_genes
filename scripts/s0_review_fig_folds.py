"""Review figures 1, 6 and 7 — the folds, a real tree, and the forest.

**Figure 1 is a schematic and says so.** It draws membrane topology — helix
count, pore loops, subunit stoichiometry — from what the catalogue records
for each fold, so it is a diagram of the project's own data rather than a
rendering of any structure. Nothing about helix packing, tilt or the shape
of a pore should be read from it.

**Figure 6 is a measurement.** The Cys-loop tree in
`results/phylogeny/tier2_cysloop/` — MAFFT, trimAl, IQ-TREE 2 with
ModelFinder and 1000 ultrafast bootstraps, rooted on the two bacterial
channels — drawn from its Newick with its real support values.

**Figure 7 is a schematic of a rule.** It shows what the phylogeny protocol
allows: trees inside superfamilies, a network between them, and nothing
spanning the two. Decision **D27** made visible.

    python3 scripts/s0_review_fig_folds.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs
import review_figlib as R
from src.catalogue import CATALOGUE, SUPERFAMILIES
from src.phylo import refused_superfamilies

OUT = ROOT / "docs" / "figures"
TREE = ROOT / "results" / "phylogeny" / "tier2_cysloop" / "cysloop.treefile"

#: (title, subunits, TM per subunit, pore loops per subunit, kind, note).
#: Counts are the catalogue's; `kind` selects helices or strands.
FOLDS = [
    ("P-loop, 6TM/1P", 4, 6, 1, "helix", "Kv, TRP, CNG, HCN"),
    ("P-loop, 2TM/1P", 4, 2, 1, "helix", "Kir, KcsA"),
    ("P-loop, 4TM/2P", 2, 4, 2, "helix", "K2P"),
    ("P-loop, 24TM/4P", 1, 24, 4, "helix", "Nav, Cav, NALCN"),
    ("Cys-loop", 5, 4, 0, "helix", "nAChR, GABA-A, GlyR"),
    ("iGluR", 4, 3, 1, "helix", "AMPA, NMDA — inverted pore"),
    ("P2X", 3, 2, 0, "helix", "P2X1–7"),
    ("DEG/ENaC", 3, 2, 0, "helix", "ENaC, ASIC"),
    ("CLC", 2, 18, 0, "helix", "two independent pores"),
    ("TMEM16 clan", 2, 10, 0, "helix", "TMEM16, OSCA, TMC"),
    ("connexin", 6, 4, 0, "helix", "gap junctions"),
    ("β-barrel", 1, 19, 0, "strand", "VDAC"),
]


def fig_folds() -> Path:
    """One subunit drawn at its true transmembrane count, plus a top view.

    Drawing three subunits side by side and truncating the helix count — the
    first version of this figure — labels a panel "6 TM" while showing four,
    which is worse than no diagram. One subunit is drawn exactly, and
    stoichiometry is carried by the top view instead.
    """
    import matplotlib.pyplot as plt
    ncol, nrow = 4, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(fs.W_FULL, 5.0))
    for ax, (title, subunits, tm, loops, kind, note) in zip(axes.ravel(), FOLDS):
        R.membrane(ax, 0, 10, 0.0, 1.0)
        step = 9.4 / max(1, tm)
        width = min(0.36, step * 0.52)
        xs = [0.3 + step * (i + 0.5) for i in range(tm)]
        for i, x in enumerate(xs):
            if kind == "strand":
                R.strand(ax, x, 0.06, 0.94, width=width)
            else:
                R.helix(ax, x, 0.06, 0.94, width=width,
                        colour=fs.BLUES[3] if i % 2 == 0 else fs.BLUES[2])
        # pore loops between the last two helices of each repeat
        if loops:
            per = tm // loops
            for L in range(loops):
                a = xs[min(len(xs) - 1, per * L + per - 2)]
                b = xs[min(len(xs) - 1, per * L + per - 1)]
                R.pore_loop(ax, a, b, 1.0, 0.5)
        R.top_view(ax, 5.0, -0.55, subunits, r=0.33)
        ax.set_xlim(-0.4, 10.4)
        ax.set_ylim(-1.75, 2.05)
        ax.axis("off")
        ax.text(5, 1.80, title, ha="center", va="center",
                fontsize=fs.FS_LABEL - 0.6, color=fs.INK)
        ax.text(5, 1.40, f"{tm} TM per subunit"
                + (f" · {loops} pore loop{'s' if loops > 1 else ''}"
                   if loops else ""),
                ha="center", va="center", fontsize=fs.FS_TICK - 1.4,
                color=fs.MUTED)
        ax.text(6.6, -0.55, f"×{subunits}", ha="left", va="center",
                fontsize=fs.FS_TICK - 0.8, color=fs.MUTED)
        ax.text(5, -1.48, note, ha="center", va="center",
                fontsize=fs.FS_TICK - 1.4, color=fs.FAINT, style="italic")
    fig.suptitle("Schematic — one subunit at its recorded transmembrane count, "
                 "with the subunit stoichiometry below. Not structures: "
                 "nothing about helix packing or pore shape is implied.",
                 fontsize=fs.FS_TICK - 0.6, color=fs.MUTED, y=0.015)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return fs.save(fig, OUT / "fig1_folds")[0]


# ------------------------------------------------------------- newick
def parse_newick(text: str):
    """Minimal Newick → nested (children, label, support, length)."""
    text = text.strip().rstrip(";")
    pos = 0

    def node():
        nonlocal pos
        children = []
        if text[pos] == "(":
            pos += 1
            while True:
                children.append(node())
                if text[pos] == ",":
                    pos += 1
                    continue
                pos += 1          # ')'
                break
        start = pos
        while pos < len(text) and text[pos] not in ",):":
            pos += 1
        name = text[start:pos]
        length = 0.0
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",)":
                pos += 1
            length = float(text[start:pos] or 0)
        support = None
        if children and name:
            try:
                support = float(name)
                name = ""
            except ValueError:
                pass
        return {"children": children, "name": name, "support": support,
                "length": length}
    return node()


def fig_tree() -> Path | None:
    import matplotlib.pyplot as plt
    if not TREE.exists():
        print("[figs] no cysloop tree — run run.py --phylo cysloop --tier 2")
        return None
    root = parse_newick(TREE.read_text())

    leaves: list[str] = []

    def collect(n):
        if not n["children"]:
            leaves.append(n["name"])
        for c in n["children"]:
            collect(c)
    collect(root)

    ypos, xmax = {}, [0.0]

    def layout(n, x=0.0):
        x += n["length"]
        xmax[0] = max(xmax[0], x)
        if not n["children"]:
            ypos[id(n)] = leaves.index(n["name"])
            return x, ypos[id(n)]
        ys = [layout(c, x)[1] for c in n["children"]]
        ypos[id(n)] = sum(ys) / len(ys)
        return x, ypos[id(n)]
    layout(root)

    NICE = {"Hs_CHRNA1": "nAChR α1  (human)", "Tm_nAChR_a": "nAChR α  (Torpedo)",
            "Hs_HTR3A": "5-HT3A", "Ls_AChBP": "AChBP  (snail) — not a channel",
            "Hs_GABRA1": "GABA-A α1", "Hs_GLRA1": "glycine α1",
            "Gv_GLIC": "GLIC  (Gloeobacter)", "Ec_ELIC": "ELIC  (Dickeya)"}
    ANION = {"Hs_GABRA1", "Hs_GLRA1"}

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 2.9))

    def draw(n, x=0.0):
        x1 = x + n["length"]
        if n["children"]:
            ys = [ypos[id(c)] for c in n["children"]]
            ax.plot([x1, x1], [min(ys), max(ys)], color=fs.MUTED, lw=1.0)
            for c in n["children"]:
                ax.plot([x1, x1 + c["length"]], [ypos[id(c)], ypos[id(c)]],
                        color=fs.MUTED, lw=1.0)
                draw(c, x1)
            if n["support"] is not None and n is not root:
                ax.text(x1 - 0.012, ypos[id(n)] + 0.16,
                        f"{int(n['support'])}", ha="right", va="bottom",
                        fontsize=fs.FS_TICK - 1.8, color=fs.ACCENT)
        else:
            colour = (fs.SUPERFAMILY["cysloop"] if n["name"] in ANION
                      else (fs.FAINT if n["name"] == "Ls_AChBP" else fs.BLUES[4]))
            ax.plot([x1], [ypos[id(n)]], "o", ms=3.4, color=colour, zorder=4)
            ax.text(x1 + 0.02, ypos[id(n)], NICE.get(n["name"], n["name"]),
                    ha="left", va="center", fontsize=fs.FS_TICK - 0.6,
                    color=fs.INK)
    draw(root)
    ax.set_xlim(-0.05, xmax[0] * 1.62)
    ax.set_ylim(-0.7, len(leaves) - 0.3)
    ax.axis("off")
    ax.text(0, len(leaves) - 0.45,
            "tier-2 protein tree · MAFFT → trimAl → IQ-TREE 2 (LG+G4, 1000 "
            "UFBoot) · rooted on GLIC/ELIC",
            fontsize=fs.FS_TICK - 1.2, color=fs.MUTED, va="top")
    ax.text(0, -0.62, "numbers are ultrafast bootstrap support",
            fontsize=fs.FS_TICK - 1.6, color=fs.ACCENT, va="bottom")
    fig.tight_layout()
    return fs.save(fig, OUT / "fig6_cysloop_tree")[0]


def fig_forest() -> Path:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import FancyBboxPatch
    refused = {k for k, _ in refused_superfamilies()}
    keys = [k for k in SUPERFAMILIES
            if any(f.superfamily == k and f.census_member()
                   for f in CATALOGUE.values())]
    keys.sort(key=lambda k: (k in refused, k))

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 3.4))
    cols = 7
    rng = np.random.default_rng(11)
    centres = {}
    for i, k in enumerate(keys):
        cx, cy = (i % cols) * 1.5, -(i // cols) * 1.35
        centres[k] = (cx, cy)
        bad = k in refused
        ax.add_patch(FancyBboxPatch((cx - 0.62, cy - 0.42), 1.24, 0.92,
                                    boxstyle="round,pad=0.02,rounding_size=0.08",
                                    facecolor="#ffffff",
                                    edgecolor=(fs.STATUS["absent"] if bad
                                               else fs.GRID),
                                    linestyle=(":" if bad else "-"),
                                    linewidth=0.8, zorder=2))
        if bad:
            ax.text(cx, cy + 0.03, "no tree", ha="center", va="center",
                    fontsize=fs.FS_TICK - 1.6, color=fs.STATUS["absent"])
        else:
            # a small cartoon tree, deterministic per box
            n = 4
            xs = np.linspace(cx - 0.42, cx + 0.42, n)
            ys = cy - 0.30 + rng.random(n) * 0.06
            for x, y in zip(xs, ys):
                ax.plot([x, x], [y, cy + 0.16], color=fs.BLUES[3], lw=0.8)
                ax.plot([x], [y], "o", ms=1.9, color=fs.BLUES[4])
            ax.plot([xs[0], xs[-1]], [cy + 0.16, cy + 0.16],
                    color=fs.BLUES[3], lw=0.8)
            ax.plot([(xs[0] + xs[-1]) / 2] * 2, [cy + 0.16, cy + 0.30],
                    color=fs.BLUES[3], lw=0.8)
        ax.text(cx, cy - 0.52, R.superfamily_label(k), ha="center", va="top",
                fontsize=fs.FS_TICK - 2.0, color=fs.INK)
    # a few network edges between boxes — structural similarity, no ancestry
    for a, b in [("iglur", "ploop"), ("tmem16_like", "clc"),
                 ("innexin_like", "connexin"), ("hv", "ploop"),
                 ("ca_release", "ploop")]:
        if a in centres and b in centres:
            (x0, y0), (x1, y1) = centres[a], centres[b]
            # drawn behind the boxes, which are opaque, so an edge reads as
            # emerging from one box and entering another rather than crossing
            # the labels in between
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="-", color=fs.ACCENT,
                                        lw=0.9, alpha=0.7,
                                        connectionstyle="arc3,rad=0.34"),
                        zorder=1)
    ax.set_xlim(-1.0, cols * 1.5)
    ax.set_ylim(-(len(keys) // cols + 1) * 1.35 - 0.2, 1.55)
    ax.axis("off")
    ax.text(-1.0, 1.48,
            "Trees inside superfamilies. A network between them. Nothing "
            "spanning the two.", fontsize=fs.FS_LABEL - 0.4, color=fs.INK,
            va="top")
    ax.text(-1.0, 1.12,
            "dotted = the catalogue marks it non-alignable, so build_tier2() "
            "refuses;  violet = a structural relationship, drawn as an edge "
            "with no branch length and no ancestor implied",
            fontsize=fs.FS_TICK - 1.6, color=fs.MUTED, va="top")
    fig.tight_layout()
    return fs.save(fig, OUT / "fig7_forest")[0]


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    fs.use()
    for f in (fig_folds, fig_tree, fig_forest):
        p = f()
        if p:
            print(f"[figs] {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
