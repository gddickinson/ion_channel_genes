"""S8b figure — the tier-2 trees (`results/phylogeny/figures/tier2_trees.png`),
drawn from the committed treefiles and tables (D13). Called by `s8_figures.py`.

* **A** — the Cys-loop tree as IQ-TREE draws it (an arbitrary drawn root: the
  declared outgroup is split, so the tree is unrooted under D48), tips
  coloured by receptor group; the four animal `plgic_prok` tips marked.
* **B** — ingroup support per unit: share of internal edges at UFBoot ≥ 95,
  70–95 and < 70.
* **C** — families (per module for multi-module families) as clades, per unit:
  one clade with UFBoot ≥ 95, one clade below, not one clade with ≤ 5 other
  tips needing to join, not one clade with more.
"""

from __future__ import annotations

import figstyle as fs
from s3_hmm_lib import read_tsv
from s7_figures import _layout
from s7_newick import parse
from s8_report_trees import tested_groups

ANIMAL_PROK = {"plgic_prok__A0A9J7KP96__Braflo", "plgic_prok__A0ABM1AAJ3__Aplcal",
               "plgic_prok__V4A2P8__Lotgig", "plgic_prok__V4AGW4__Lotgig"}
CYS_GROUPS = [  # (label, families, colour) — five groups, within the palette's limit
    ("nicotinic ACh receptors", {"nachr"}, fs.SUPERFAMILY["ploop"]),
    ("5-HT3 and ZAC", {"ht3", "zac"}, fs.SUPERFAMILY["deg_enac"]),
    ("anion-selective (GABA-A, GlyR, invertebrate GluCl)",
     {"gabaa", "glyr", "plgic_invertebrate"}, fs.SUPERFAMILY["cysloop"]),
    ("prokaryotic GLIC, ELIC (declared outgroup)", {"plgic_prok"},
     fs.SUPERFAMILY["iglur"]),
]
UNIT_LABEL = {"ploop": "P-loop (module)", "cysloop": "Cys-loop", "deg_enac": "DEG/ENaC",
              "iglur": "iGluR (module)", "innexin_like": "innexin (module)",
              "p2x": "P2X", "ca_release": "ITPR/RyR (module)"}


def _colour(tip: str) -> str:
    if tip in ANIMAL_PROK:
        return fs.ACCENT
    fam = tip.split("__", 1)[0]
    return next((c for _, f, c in CYS_GROUPS if fam in f), fs.FAINT)


def panel_cysloop(ax, d) -> None:
    _, _, items = _layout(parse((d / "tier2" / "cysloop.treefile").read_text()))
    xmax = max(it[2] for it in items)
    for it in items:
        if it[0] == "node":
            _, name, x, y, kids = it
            ax.plot([x, x], [kids[0][1], kids[-1][1]], color=fs.FAINT, lw=0.35)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=fs.FAINT, lw=0.35)
        else:
            _, name, x, y = it
            ax.plot([x, xmax * 1.02], [y, y], color=_colour(name), lw=0.9,
                    solid_capstyle="butt", alpha=0.6 if name not in ANIMAL_PROK else 1)
            if name in ANIMAL_PROK:
                ax.scatter(xmax * 1.04, y, marker="<", s=10, color=fs.ACCENT, lw=0, zorder=4)
    ys = [it[3] for it in items if it[0] == "leaf" and it[1] in ANIMAL_PROK]
    ax.text(xmax * 1.07, sum(ys) / len(ys),
            "4 animal\n'prokaryotic'\npLGICs",
            fontsize=fs.FS_NOTE - 0.4, color=fs.ACCENT, va="center")
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlim(0, xmax * 1.4)
    ax.set_xlabel("substitutions / site (drawn root arbitrary)", fontsize=fs.FS_LABEL)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=c, lw=2, label=l) for l, _, c in CYS_GROUPS]
    h.append(Line2D([], [], color=fs.ACCENT, lw=2, label="animal 'prokaryotic' pLGIC"))
    ax.legend(handles=h, fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="lower left",
              bbox_to_anchor=(0.0, -0.17), ncol=1)
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "A", "Cys-loop receptors, 282 tips (unrooted)")


def _order(trees) -> list[dict]:
    return sorted(trees, key=lambda t: -int(t["n_tips"]))


def panel_support(ax, d) -> None:
    trees = _order(read_tsv(d / "tier2_trees.tsv"))
    cols = [(fs.BLUES[5], "UFBoot ≥ 95"), (fs.BLUES[2], "70–95"), (fs.GRID, "< 70")]
    for i, t in enumerate(trees):
        hi, lo = float(t["frac_ge95"]), float(t["frac_lt70"])
        left = 0.0
        for (c, _), w in zip(cols, (hi, 1 - hi - lo, lo)):
            ax.barh(i, w, left=left, color=c, height=0.7, lw=0)
            left += w
        ax.text(1.02, i, f"{t['n_tips']} tips / {t['informative']} sites",
                va="center", fontsize=fs.FS_NOTE - 0.4, color=fs.MUTED)
    ax.set_yticks(range(len(trees)), [UNIT_LABEL[t["unit"]] for t in trees],
                  fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.6)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("share of internal edges", fontsize=fs.FS_LABEL)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc=c, label=l) for c, l in cols], fontsize=fs.FS_NOTE - 0.4,
              frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.2), ncol=3)
    fs.despine(ax)
    fs.panel(ax, "B", "Branch support per tree")


def panel_clades(ax, d) -> None:
    trees = _order(read_tsv(d / "tier2_trees.tsv"))
    tg = tested_groups(read_tsv(d / "tier2_families.tsv"))
    cats = [("one clade, UFBoot ≥ 95", fs.BLUES[5]), ("one clade, < 95", fs.BLUES[2]),
            ("not one clade, ≤ 5 tips away", "#ef9a90"),
            ("not one clade, > 5 tips away", "#b3261e")]

    def cat(r):
        if r["one_clade"] == "True":
            return 0 if r["ufboot"] and float(r["ufboot"]) >= 95 else 1
        return 2 if int(r["n_intruders"]) <= 5 else 3
    for i, t in enumerate(trees):
        rows = [r for r in tg if r["unit"] == t["unit"]]
        n = [sum(cat(r) == k for r in rows) for k in range(4)]
        left = 0
        for (_, c), w in zip(cats, n):
            ax.barh(i, w, left=left, color=c, height=0.7, lw=0)
            left += w
        ax.text(left + 0.4, i, f"{n[0] + n[1]}/{len(rows)}", va="center",
                fontsize=fs.FS_NOTE - 0.4, color=fs.MUTED)
    ax.set_yticks(range(len(trees)), [UNIT_LABEL[t["unit"]] for t in trees],
                  fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("families (per module where a chain has several)", fontsize=fs.FS_LABEL)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(fc=c, label=l) for l, c in cats], fontsize=fs.FS_NOTE - 0.4,
              frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "C", "Is each family one clade?")


def figure(d, out) -> list:
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(fs.W_FULL, 6.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], hspace=0.45, wspace=0.75)
    panel_cysloop(fig.add_subplot(gs[:, 0]), d)
    panel_support(fig.add_subplot(gs[0, 1]), d)
    panel_clades(fig.add_subplot(gs[1, 1]), d)
    return fs.save(fig, out / "tier2_trees")
