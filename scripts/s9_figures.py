"""s9_figures.py — S9's headline figure, drawn from the committed tables (D13).

    bin/envpy scripts/s9_figures.py

`results/filter_atlas/figures/filter_atlas.png` (+ .pdf):

A. The Cav tier-1 tree, each chain coloured by its four-repeat filter string —
   the Q5 test case (EEDD, the T-type locus).
B. The four-repeat filter strings of every Nav, Cav, NALCN, CatSper and TPC
   chain, as shares of each family.
C. Each filter position's retention index as a percentile of every
   parsimony-informative column of the same tree's alignment (the post-hoc
   column-background control).
"""

from __future__ import annotations

import sys
from collections import Counter

import figstyle as fs
from s3_hmm_lib import read_tsv
from s7_figures import _layout
from s7_newick import parse
from s9_lib import OUT_DIR, PHY_DIR, is_missing

CAV_COLOURS = [("EEEE", fs.SUPERFAMILY["ploop"], "EEEE (Cav1/Cav2-type)"),
               ("EEDD", fs.SUPERFAMILY["cysloop"], "EEDD (T-type)"),
               ("DDDD", fs.SUPERFAMILY["iglur"], "DDDD")]
OTHER, OUTGROUP = fs.ACCENT, fs.GRID
REPEAT_FAMS = [("nav", "Nav"), ("cav", "Cav"), ("nalcn", "NALCN"), ("catsper", "CatSper"),
               ("tpc", "TPC")]
FAM_LABEL = {"kv_shaker": "Kv (Shaker)", "kv_kcnq": "Kv (KCNQ)", "kv_eag": "EAG",
             "kv_modifier": "Kv modifier", "kca_slo": "Slo", "kca_sk": "SK", "kir": "Kir",
             "k2p": "K2P", "hcn": "HCN", "cng": "CNG", "trpa": "TRPA", "trpc": "TRPC",
             "trpm": "TRPM", "trpml": "TRPML", "trpn": "TRPN", "trpp": "TRPP",
             "trpv": "TRPV", "nav": "Nav", "cav": "Cav", "nalcn": "NALCN",
             "catsper": "CatSper", "tpc": "TPC", "kcsa_prok": "KcsA-like"}


def panel_tree(ax, chains) -> None:
    tree = parse((PHY_DIR / "tier1" / "cav.treefile").read_text())
    _, _, items = _layout(tree)
    xmax = max(it[2] for it in items)
    col = {s: c for s, c, _ in CAV_COLOURS}
    for it in items:
        if it[0] == "node":
            _, _, x, _, kids = it
            ax.plot([x, x], [kids[0][1], kids[-1][1]], color=fs.FAINT, lw=0.3)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=fs.FAINT, lw=0.3)
        else:
            _, name, x, y = it
            s = chains.get(name)
            c = OUTGROUP if name.startswith("OG_") else (
                fs.GRID if is_missing(s) else col.get(s, OTHER))
            ax.plot([x, xmax * 1.03], [y, y], color=c, lw=0.9, solid_capstyle="butt")
            ax.scatter(xmax * 1.05, y, s=3, color=c, lw=0)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlim(0, xmax * 1.1)
    ax.set_xlabel("substitutions / site", fontsize=fs.FS_LABEL)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, lw=2, label=l) for _, c, l in CAV_COLOURS]
    h += [Line2D([], [], color=OTHER, lw=2, label="other string"),
          Line2D([], [], color=fs.GRID, lw=2, label="a repeat unread, or outgroup")]
    ax.legend(handles=h, fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="lower left",
              bbox_to_anchor=(0.0, -0.2), ncol=2)
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "A", "Cav tree (201 chains), tips coloured by filter")


def panel_strings(ax, rows) -> None:
    shades = [fs.BLUES[4], fs.BLUES[2], fs.BLUES[1], fs.BLUES[0]]
    labels = []
    for i, (fam, lab) in enumerate(REPEAT_FAMS):
        st = [r["string"] for r in rows if r["family"] == fam]
        read = [s for s in st if not is_missing(s)]
        top = Counter(read).most_common(3)
        parts = top + [("other", len(read) - sum(n for _, n in top)),
                       ("unread", len(st) - len(read))]
        left = 0.0
        for j, (s, n) in enumerate(parts):
            if not n:
                continue
            w = n / len(st)
            c = shades[j] if j < 3 else (fs.HILITE if s == "other" else fs.GRID)
            ax.barh(i, w, left=left, color=c, edgecolor="white", lw=0.4, height=0.7)
            if w >= 0.1:
                ax.text(left + w / 2, i, s, ha="center", va="center",
                        fontsize=fs.FS_NOTE - 0.6,
                        color="white" if j == 0 else fs.INK)
            left += w
        labels.append(f"{lab} ({len(st)})")
    ax.set_yticks(range(len(REPEAT_FAMS)), labels, fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of chains (family, number of chains)", fontsize=fs.FS_LABEL)
    fs.despine(ax)
    fs.panel(ax, "B", "Filter strings, repeats I→IV")


def _kind(fam: str, pos: str) -> tuple[str, str]:
    if pos.startswith("repeat"):
        return "four-repeat locus", fs.SUPERFAMILY["cysloop"]
    if pos == "window 2":
        return "K window, x of TxGYG", fs.SUPERFAMILY["ploop"]
    if pos == "window 4":
        return "K window, Y/F of GYG", fs.SUPERFAMILY["iglur"]
    return "other window position", fs.FAINT


def panel_background(ax, rows) -> None:
    rows = [r for r in rows if r["tree"] == "tier1" and r["percentile"] != ""]
    fams = sorted({r["family"] for r in rows},
                  key=lambda f: (f not in dict(REPEAT_FAMS), FAM_LABEL[f]))
    seen = {}
    for i, fam in enumerate(fams):
        for r in rows:
            if r["family"] != fam:
                continue
            lab, c = _kind(fam, r["position"])
            h = ax.scatter(float(r["percentile"]), i, s=11, color=c, lw=0, zorder=3)
            seen.setdefault(lab, h)
    ax.axvline(50, color=fs.GRID, lw=0.8, zorder=1)
    ax.set_yticks(range(len(fams)), [FAM_LABEL[f] for f in fams], fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("percentile among the alignment's columns\n(retention index; "
                  "right = more tree-congruent)", fontsize=fs.FS_LABEL)
    order = ["four-repeat locus", "K window, x of TxGYG", "K window, Y/F of GYG",
             "other window position"]
    ax.legend([seen[k] for k in order if k in seen], [k for k in order if k in seen],
              fontsize=fs.FS_NOTE - 0.6, frameon=False, loc="upper left",
              bbox_to_anchor=(-0.05, -0.17), ncol=2)
    fs.hgrid(ax, "x")
    fs.despine(ax)
    fs.panel(ax, "C", "Filter positions vs every other column")


def main() -> None:
    fs.use()
    import matplotlib.pyplot as plt
    rows = read_tsv(OUT_DIR / "filter_chains.tsv")
    chains = {r["chain"]: r["string"] for r in rows if r["family"] == "cav"}
    fig = plt.figure(figsize=(fs.W_FULL, 6.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[0.8, 1.6],
                          hspace=0.55, wspace=0.55)
    panel_tree(fig.add_subplot(gs[:, 0]), chains)
    panel_strings(fig.add_subplot(gs[0, 1]), rows)
    panel_background(fig.add_subplot(gs[1, 1]), read_tsv(OUT_DIR / "column_background.tsv"))
    out = OUT_DIR / "figures"
    out.mkdir(parents=True, exist_ok=True)
    for p in fs.save(fig, out / "filter_atlas"):
        print(p)


if __name__ == "__main__":
    sys.exit(main())
