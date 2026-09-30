"""S5b figure — the presence matrix: every census family in every panel genome.

Species (rows, panel order) × census families (columns, grouped by
superfamily), each cell coloured by its S5 verdict on the same ordered
evidence scale as the r4 figure (`s5r4_figures.VERDICT`). The matrix S10's
repertoire reconstruction starts from. Drawn from
`results/genome_sweep/cells.tsv`.

    python3 scripts/s5_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402
from s5_lib import manifest  # noqa: E402
from s5r4_figures import VERDICT  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

D = ROOT / "results" / "genome_sweep"


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    cells = {(c["species"], c["family"]): c["verdict"] for c in read_tsv(D / "cells.tsv")}
    man = manifest()
    species = [s for s in man if any(k[0] == s for k in cells)]
    order = {sf: i for i, sf in enumerate(fs.SUPERFAMILY_ORDER)}
    fams = sorted({k[1] for k in cells},
                  key=lambda f: (order.get(CATALOGUE[f].superfamily, 99),
                                 CATALOGUE[f].superfamily, f))
    col = {v: c for v, _, c in VERDICT}
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 7.2))
    for i, s in enumerate(species):
        for j, f in enumerate(fams):
            ax.add_patch(plt.Rectangle((j, i), 0.9, 0.86,
                                       color=col.get(cells.get((s, f), ""), "white"),
                                       linewidth=0))
    prev = None
    for j, f in enumerate(fams):
        sf = CATALOGUE[f].superfamily
        if sf != prev and j:
            ax.axvline(j - 0.05, color=fs.INK, lw=0.4)
        prev = sf
    prev = None
    for i, s in enumerate(species):
        if man[s]["group"] != prev and i:
            ax.axhline(i - 0.07, color=fs.GRID, lw=0.5)
        prev = man[s]["group"]
    ax.set_xlim(0, len(fams))
    ax.set_ylim(len(species), 0)
    ax.set_xticks([j + 0.45 for j in range(len(fams))], fams, rotation=90,
                  fontsize=fs.FS_TICK - 2.6)
    ax.xaxis.tick_top()
    ax.set_yticks([i + 0.43 for i in range(len(species))],
                  [s if len(s) < 26 else s[:24] + "…" for s in species],
                  fontsize=fs.FS_TICK - 2.4, style="italic")
    for side in ("left", "right", "bottom", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(handles=[Patch(color=c, label=l) for _, l, c in VERDICT if l],
              fontsize=fs.FS_NOTE - 0.6, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.005), ncol=4)
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "presence_matrix"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
