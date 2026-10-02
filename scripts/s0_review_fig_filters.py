"""Review figures 2 and 3 — the two selectivity-filter alignments.

Both are rendered from `results/s0_baseline/filter_*.tsv`, which
`scripts/s0_filter_atlas.py` computed from real sequences. Neither is a
schematic.

**Figure 2 — the potassium signature.** Nineteen proteins spanning the
potassium branch, from the *Streptomyces* prototype to a cyanobacterial
glutamate receptor, each showing the TxGYG motif found in its own sequence
with seven residues of context on each side. No alignment was performed and
none is needed: the motif is located directly. The one panel member that
fails is *Bacillus* NaK, whose filter reads TVGDG — one substitution from the
signature, and not potassium-selective.

**Figure 3 — the four-repeat locus.** The four residues that separate the
sodium, calcium and leak channels, projected from human Nav1.5 by MAFFT
(**D26**), each shown with its context in the query. Two-repeat channels are
included so the figure shows the method declining to call rather than only
its successes.

    python3 scripts/s0_review_fig_filters.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs
import review_figlib as R
from scripts.s0_lib import read_tsv

DATA = ROOT / "results" / "s0_baseline"
OUT = ROOT / "docs" / "figures"

NICE = {
    "Sl_KcsA": "KcsA  (S. lividans)", "Mt_MthK": "MthK  (archaeal)",
    "Bc_NaK": "NaK  (B. cereus)", "Dm_Shaker": "Shaker  (Drosophila)",
    "Hs_KCNA1": "Kv1.1  KCNA1", "Hs_KCNB1": "Kv2.1  KCNB1",
    "Hs_KCND2": "Kv4.2  KCND2", "Hs_KCNQ1": "Kv7.1  KCNQ1",
    "Hs_KCNQ2": "Kv7.2  KCNQ2", "Hs_KCNH2": "hERG  KCNH2",
    "Hs_KCNH1": "Kv10.1  KCNH1", "Hs_KCNMA1": "BK  KCNMA1",
    "Hs_KCNJ2": "Kir2.1  KCNJ2", "Hs_KCNJ11": "Kir6.2  KCNJ11",
    "Hs_KCNK2": "TREK-1  KCNK2", "Hs_KCNK3": "TASK-1  KCNK3",
    "Hs_HCN1": "HCN1", "Hs_HCN4": "HCN4", "Ss_GluR0": "GluR0  (Synechocystis)",
    "Hs_SCN5A": "Nav1.5  SCN5A", "Hs_SCN1A": "Nav1.1  SCN1A",
    "Hs_CACNA1C": "Cav1.2  CACNA1C", "Hs_CACNA1G": "Cav3.1  CACNA1G",
    "Hs_NALCN": "NALCN", "Hs_TPCN1": "TPC1  TPCN1",
    "Hs_CATSPER1": "CatSper1",
}
FAMILY_NOTE = {
    "Hs_SCN5A": "Na+", "Hs_SCN1A": "Na+", "Hs_CACNA1C": "Ca2+",
    "Hs_CACNA1G": "Ca2+ (T-type)", "Hs_NALCN": "leak",
    "Hs_TPCN1": "two repeats", "Hs_CATSPER1": "one repeat",
}


def fig_k_filter(rows: list[dict]) -> Path:
    import matplotlib.pyplot as plt
    rows = [r for r in rows if r["status"] in ("ok", "NEAR_MISS")]
    ncols = max(len(r["left"]) + len(r["motif"]) + len(r["right"]) for r in rows)
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 0.30 * len(rows) + 1.15))
    for i, r in enumerate(rows):
        y = len(rows) - i - 1
        seq = r["left"] + r["motif"] + r["right"]
        hl = (len(r["left"]), len(r["left"]) + len(r["motif"]))
        R.draw_sequence_row(ax, y, seq, highlight=hl)
        ax.text(-0.7, y + 0.5, NICE.get(r["label"], r["label"]),
                ha="right", va="center", fontsize=fs.FS_TICK - 0.4,
                color=fs.INK)
        if r["status"] == "NEAR_MISS":
            ax.text(ncols + 0.7, y + 0.5,
                    "not K+-selective — Y→D", ha="left", va="center",
                    fontsize=fs.FS_TICK - 1.0, color=fs.STATUS["absent"])
    ax.set_xlim(-9.5, ncols + 8.5)
    ax.set_ylim(-1.5, len(rows) + 0.2)
    ax.axis("off")
    left = max(len(r["left"]) for r in rows)
    width = max(len(r["motif"]) for r in rows)
    ax.text(left + width / 2, len(rows) + 0.05, "selectivity filter",
            ha="center", va="bottom", fontsize=fs.FS_LABEL - 0.6,
            color=fs.MUTED)
    R.residue_legend(ax, loc="lower center", ncol=4)
    fig.tight_layout()
    return fs.save(fig, OUT / "fig2_potassium_filter", legend_in="docs/review (numbered review figure)")[0]


def fig_four_repeat(rows: list[dict]) -> Path:
    import matplotlib.pyplot as plt
    rows = [r for r in rows if r["status"] in ("ok", "PARTIAL")]
    fig, axes = plt.subplots(1, 4, figsize=(fs.W_FULL, 0.40 * len(rows) + 1.35),
                             sharey=True)
    for rep, ax in enumerate(axes):
        for i, r in enumerate(rows):
            y = len(rows) - i - 1
            ctx = (r["contexts"].split("|") + ["", "", "", ""])[rep]
            if not ctx:
                ax.text(5.5, y + 0.5, "no equivalent position",
                        ha="center", va="center", fontsize=fs.FS_TICK - 1.4,
                        color=fs.FAINT, style="italic")
                continue
            mid = len(ctx) // 2
            R.draw_sequence_row(ax, y, ctx, highlight=(mid, mid + 1),
                                fontsize=fs.FS_TICK - 2.0)
        ax.set_xlim(-0.5, 11.5)
        ax.set_ylim(-0.4, len(rows) + 0.9)
        ax.axis("off")
        ax.text(5.5, len(rows) + 0.15, f"repeat {'I II III IV'.split()[rep]}",
                ha="center", va="bottom", fontsize=fs.FS_LABEL - 0.6,
                color=fs.MUTED)
    for i, r in enumerate(rows):
        y = len(rows) - i - 1
        axes[0].text(-1.2, y + 0.5, NICE.get(r["label"], r["label"]),
                     ha="right", va="center", fontsize=fs.FS_TICK - 0.4)
        sig = r["signature"].replace("-", "·")
        axes[-1].text(12.3, y + 0.5, f"{sig}   {FAMILY_NOTE.get(r['label'],'')}",
                      ha="left", va="center", fontsize=fs.FS_TICK - 0.2,
                      color=(fs.INK if r["status"] == "ok" else fs.FAINT),
                      fontweight="bold" if r["status"] == "ok" else "normal")
    fig.tight_layout(rect=(0.10, 0, 0.88, 1))
    return fs.save(fig, OUT / "fig3_four_repeat_filter", legend_in="docs/review (numbered review figure)")[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    fs.use()
    k = read_tsv(DATA / "filter_k.tsv")
    f4 = read_tsv(DATA / "filter_four_repeat.tsv")
    if not k or not f4:
        print("[figs] run scripts/s0_filter_atlas.py first")
        return 1
    print(f"[figs] {fig_k_filter(k)}")
    print(f"[figs] {fig_four_repeat(f4)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
