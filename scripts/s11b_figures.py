"""s11b_figures.py — S11b's headline figure, drawn from the committed tables (D13).

    bin/envpy scripts/s11b_figures.py

`results/duplication/figures/s11b_repeats_roots.png` (+ .pdf):

A. The whole-repeat ML tree (unrooted, drawn root arbitrary), tips coloured by
   repeat class I–IV; TPC's two repeats in grey with markers.
B. AU test p-values of the three repeat pairings, whole repeat and pore module,
   with the 0.05 rejection line.
C. Rootstrap of the non-reversible ML root per family — S7d's six and the five
   controls — marked by agreement with S11a's reconciliation root (six) or with
   the declared outgroup root (controls).
"""
from __future__ import annotations

import sys
from pathlib import Path

import figstyle as fs
from s3_hmm_lib import read_tsv
from s7_figures import _layout
from s7_newick import parse

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "duplication"
CLASS_COL = {"I": fs.SUPERFAMILY["ploop"], "II": fs.SUPERFAMILY["cysloop"],
             "III": fs.SUPERFAMILY["iglur"], "IV": fs.ACCENT,
             "tI": fs.FAINT, "tII": fs.MUTED}
HYP_LABEL = {"H13": "{I,III} | {II,IV}", "H12": "{I,II} | {III,IV}",
             "H14": "{I,IV} | {II,III}"}


def panel_tree(ax) -> None:
    tf = OUT / "repeats" / "repeat.treefile"
    _, _, items = _layout(parse(tf.read_text()))
    xmax = max(it[2] for it in items)
    for it in items:
        if it[0] == "node":
            _, _, x, y, kids = it
            ax.plot([x, x], [kids[0][1], kids[-1][1]], color=fs.GRID, lw=0.5)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=fs.GRID, lw=0.5)
        else:
            _, name, x, y = it
            cls, fam = name.split("__")[:2]
            ax.plot([x, xmax * 1.02], [y, y], color=CLASS_COL[cls], lw=1.1,
                    solid_capstyle="butt")
            if cls.startswith("t"):
                ax.scatter(xmax * 1.05, y, s=9, lw=0, color=CLASS_COL[cls],
                           marker="o" if cls == "tI" else "s", zorder=4)
            else:
                ax.text(xmax * 1.05, y, fam, fontsize=4.2, va="center", color=fs.MUTED)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlim(0, xmax * 1.25)
    ax.set_xlabel("substitutions / site (drawn root arbitrary)", fontsize=fs.FS_LABEL)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=CLASS_COL[c], lw=2, label=f"repeat {c}") for c in
         ("I", "II", "III", "IV")]
    h += [Line2D([], [], color=CLASS_COL["tI"], marker="o", lw=2, ms=3, label="TPC repeat I"),
          Line2D([], [], color=CLASS_COL["tII"], marker="s", lw=2, ms=3, label="TPC repeat II")]
    ax.legend(handles=h, fontsize=fs.FS_NOTE - 0.4, frameon=True, framealpha=1,
              edgecolor=fs.GRID, loc="lower left", handlelength=1.4)
    fs.despine(ax, keep=("bottom",))
    n = sum(1 for it in items if it[0] == "leaf")
    fs.panel(ax, "A", f"Whole S1–S6 repeats, {n} tips (unrooted)")


def panel_au(ax) -> None:
    au = [r for r in read_tsv(OUT / "repeat_au.tsv") if r["tree"] != "ML"]
    units = [("repeat", "whole repeat", fs.BLUES[4], "o"),
             ("pore", "pore module", fs.BLUES[2], "D")]
    hyps = list(HYP_LABEL)
    for k, (u, lab, col, mk) in enumerate(units):
        for i, h in enumerate(hyps):
            r = next(x for x in au if x["unit"] == u and x["tree"] == h)
            p = max(float(r["p_AU"]), 1e-5)
            ax.scatter(p, i + (k - 0.5) * 0.25, color=col, marker=mk, s=22, lw=0,
                       label=lab if i == 0 else None, zorder=3)
            ax.text(p * 1.35, i + (k - 0.5) * 0.25, f"{float(r['p_AU']):.3g}",
                    va="center", fontsize=fs.FS_NOTE - 0.4, color=fs.MUTED)
    ax.axvline(0.05, color=fs.ACCENT, lw=0.8, ls="--")
    ax.text(0.05, -0.75, " p = 0.05", color=fs.ACCENT, fontsize=fs.FS_NOTE, va="center")
    ax.set_xscale("log")
    ax.set_xlim(1e-5, 3)
    ax.set_ylim(len(hyps) - 0.4, -1.0)
    ax.set_yticks(range(len(hyps)), [f"{h}  {HYP_LABEL[h]}" for h in hyps],
                  fontsize=fs.FS_TICK)
    ax.set_xlabel("AU test p-value (log scale; left of the line = rejected)")
    ax.legend(fontsize=fs.FS_NOTE, frameon=True, framealpha=1, edgecolor=fs.GRID,
              loc="lower left")
    fs.despine(ax)
    fs.panel(ax, "B", "Which repeats pair: AU test")


def panel_roots(ax) -> None:
    p = OUT / "roots_nonrev.tsv"
    inputs = read_tsv(OUT / "roots_inputs.tsv")
    res = {r["family"]: r for r in read_tsv(p)} if p.exists() else {}
    for i, row in enumerate(inputs):
        fam, r = row["family"], res.get(row["family"])
        if r is None:
            ax.text(1, i, "running", va="center", fontsize=fs.FS_NOTE, color=fs.FAINT,
                    style="italic")
            continue
        ok = (r["control_recovered"] if row["role"] == "control" else r["agrees_s11a"]) == "True"
        v = float(r["root_rootstrap"])
        ax.barh(i, v, color=fs.BLUES[4] if ok else fs.FAINT, height=0.62, lw=0)
        worse = float(r["dlnl_nonrev"]) <= 0
        ax.text(v + 1.5, i, f"{v:.0f}" + (" (Q.pfam fits better)" if worse else ""),
                va="center", fontsize=fs.FS_NOTE - 0.4, color=fs.MUTED)
    ax.axvline(95, color=fs.ACCENT, lw=0.8, ls="--")
    n6 = sum(r["role"] == "root" for r in inputs)
    ax.axhline(n6 - 0.5, color=fs.GRID, lw=0.8)
    ax.set_yticks(range(len(inputs)),
                  [f"{r['family']} ({r['n_tips']})" for r in inputs], fontsize=fs.FS_TICK)
    ax.text(97, -0.9, "95", color=fs.ACCENT, fontsize=fs.FS_NOTE)
    ax.set_xlim(0, 118)
    ax.set_ylim(len(inputs) + 1.4, -1.2)
    ax.set_xlabel("rootstrap of the non-reversible ML root (%)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=fs.BLUES[4], label="same root as reconciliation "
                             "(six) / declared outgroup (controls)"),
                       Patch(facecolor=fs.FAINT, label="different root")],
              fontsize=fs.FS_NOTE - 0.4, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.despine(ax)
    fs.panel(ax, "C", "Outgroup-free roots: S7d's six above, controls below")


def main() -> int:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    fs.use()
    fig = plt.figure(figsize=(fs.W_FULL, 7.4))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1], hspace=0.42, wspace=0.55,
                  left=0.04, right=0.97, top=0.93, bottom=0.07)
    panel_tree(fig.add_subplot(gs[:, 0]))
    panel_au(fig.add_subplot(gs[0, 1]))
    panel_roots(fig.add_subplot(gs[1, 1]))
    (OUT / "figures").mkdir(exist_ok=True)
    fs.save(fig, OUT / "figures" / "s11b_repeats_roots")
    return 0


if __name__ == "__main__":
    sys.exit(main())
