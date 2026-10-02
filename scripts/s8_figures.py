"""S8 figure — the tier-3 fold network, measured (and, once S8b has parsed
them, the tier-2 pore-module trees).

Drawn from the committed tables (D13), through `figstyle.py`:

* **A** — superfamily × superfamily median TM-score (TM-align, average-length
  normalisation; one AFDB model per census family cut to its comparison unit),
  ordered by average-linkage clustering of the matrix, on the full 0–1 scale.
  Hatched cells have no pair to measure (the diagonal of a one-family
  superfamily). Literature edges boxed.
* **B** — each literature edge: every family-pair TM-score behind it, the
  0.5 bar fixed before measuring (D48), and the best median any *other*
  superfamily pair touching either end reaches — the background the rank
  test reads against.

    bin/envpy scripts/s8_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402

NET = ROOT / "results" / "phylogeny" / "fold_network"
OUT = ROOT / "results" / "phylogeny" / "figures"
VERDICT = {"supported": fs.BLUES[5], "not_distinguished": fs.FAINT,
           "unmeasured": fs.GRID}


def _sf(node: str, sfa: str) -> str:
    return "ploop:VSD" if node == "kv_shaker:VSD" else sfa


def matrix_order(med: dict, names: list[str]) -> list[str]:
    import numpy as np
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform
    n = len(names)
    d = np.zeros((n, n))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i != j:
                d[i, j] = 1 - med.get(tuple(sorted((a, b))), 0.0)
    return [names[i] for i in leaves_list(linkage(squareform(d), "average"))]


def panel_matrix(ax, fig) -> None:
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle
    rows = read_tsv(NET / "superfamily_medians.tsv")
    med = {tuple(sorted((r["a"], r["b"]))): float(r["median_tm"]) for r in rows}
    names = sorted({x for k in med for x in k})
    order = matrix_order(med, names)
    m = np.full((len(order), len(order)), np.nan)
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            m[i, j] = med.get(tuple(sorted((a, b))), np.nan)
    # Full 0–1 scale, no clipping; the ramp starts at a pale tint, never white,
    # so an unmeasured cell (no within-superfamily pair: one census family) is
    # drawn apart from every measured value.
    cmap = LinearSegmentedColormap.from_list("tm", ["#f4f8fd", *fs.BLUES, "#0b2a52"])
    im = ax.imshow(np.ma.masked_invalid(m), cmap=cmap, vmin=0.0, vmax=1.0,
                   interpolation="nearest")
    for i, j in zip(*np.where(np.isnan(m))):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc=fs.GRID, ec=fs.FAINT,
                               lw=0, hatch="////"))

    ax.set_xticks(range(len(order)), order, rotation=90, fontsize=5)
    ax.set_yticks(range(len(order)), order, fontsize=5)
    ax.tick_params(length=0)
    for e in read_tsv(NET / "literature_edges.tsv"):
        if e["verdict"] == "unmeasured":
            continue
        i, j = order.index(e["a"]), order.index(e["b"])
        for (x, y) in {(i, j), (j, i)}:
            ax.add_patch(Rectangle(
                (y - 0.5, x - 0.5), 1, 1, fill=False, lw=1.0,
                ec=fs.ACCENT if e["verdict"] == "supported" else fs.INK))
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("median TM-score", fontsize=6)
    cb.ax.tick_params(labelsize=5)
    from matplotlib.patches import Patch
    cb.ax.legend(handles=[Patch(fc=fs.GRID, ec=fs.FAINT, lw=0, hatch="////",
                                label="no pair:\none family")],
                 fontsize=5, frameon=False, loc="upper left",
                 bbox_to_anchor=(-0.3, -0.04), handlelength=1.2)
    fs.panel(ax, "A", "Fold similarity between superfamilies (tier 3)")


def panel_edges(ax) -> None:
    pairs = read_tsv(NET / "pairs.tsv")
    edges = read_tsv(NET / "literature_edges.tsv")
    for i, e in enumerate(edges):
        key = tuple(sorted((e["a"], e["b"])))
        vals = [float(p["tm_avg"]) for p in pairs
                if tuple(sorted((_sf(p["a"], p["sf_a"]), _sf(p["b"], p["sf_b"])))) == key]
        xs = [i + 0.25 * ((k % 7) / 6 - 0.5) for k in range(len(vals))]
        ax.scatter(xs, vals, s=5, color=VERDICT[e["verdict"]], lw=0, zorder=3)
        if e["median_tm"]:
            ax.hlines(float(e["median_tm"]), i - 0.25, i + 0.25, color=fs.INK, lw=1.2)
        if e["best_other_tm"]:
            ax.hlines(float(e["best_other_tm"]), i - 0.3, i + 0.3, color=fs.ACCENT,
                      lw=0.8, ls=(0, (2, 1.5)))
    ax.axhline(0.5, color=fs.MUTED, lw=0.6, ls=(0, (4, 2)))
    short = {"tmem16_like": "TMEM16", "iglur": "iGluR", "ploop": "P-loop",
             "innexin_like": "innexin", "connexin": "connexin", "ca_release": "ITPR",
             "hv": "Hv1", "ploop:VSD": "Kv VSD"}
    verdict = {"supported": "supported", "not_distinguished": "not dist.",
               "unmeasured": "unmeasured"}
    labels = [f"{short.get(e['a'], e['a'])}–\n{short.get(e['b'], e['b'])}\n"
              f"{verdict[e['verdict']]}" for e in edges]
    ax.set_xticks(range(len(edges)), labels, fontsize=5)
    ax.set_ylabel("TM-score, family pairs", fontsize=6)
    ax.set_ylim(0.2, 1.0)
    ax.set_xlim(-0.6, len(edges) - 0.4)
    ax.tick_params(labelsize=5)
    ax.text(len(edges) - 0.45, 0.51, "0.5 bar (D48)", fontsize=5, color=fs.MUTED,
            ha="right", va="bottom")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=fs.INK, lw=1.2, label="edge median"),
                       Line2D([], [], color=fs.ACCENT, lw=0.8, ls=(0, (2, 1.5)),
                              label="best other pair at either end")],
              fontsize=5, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "B", "The literature edges, measured")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(fs.W_FULL, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.35)
    panel_matrix(fig.add_subplot(gs[0]), fig)
    panel_edges(fig.add_subplot(gs[1]))
    print("\n".join(map(str, fs.save(fig, OUT / "fold_network"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
