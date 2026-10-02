"""S5 r4 figure — the seven families added after S20, across the 52 genomes.

One panel: species (rows, in panel order, grouped) × new census family
(columns), each cell coloured by its S5 verdict on figstyle's ordered
evidence scale — in the proteome, found in the genome only, partial / gap,
trace, controlled absence; near-white where no in-group bait makes the cell
informative. Drawn from `results/genome_sweep/r4_cells.tsv`.

    python3 scripts/s5r4_figures.py
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

D = ROOT / "results" / "genome_sweep"
FAMS = ["pacc", "tmco1", "tmem87", "tmem109", "clcc1", "mitok", "gphr"]
VERDICT = [  # (verdict, legend label, colour) — best evidence first
    ("present", "in the proteome", fs.STATUS["found_annotated"]),
    ("genome_found", "genome only: proteome miss", fs.STATUS["found_unannotated"]),
    ("genome_present", "genome-only species", fs.STATUS["found_no_annotation"]),
    ("genome_weak", "genome, weak", fs.STATUS["tblastn_trace_ambiguous"]),
    ("partial", "partial / gap", fs.STATUS["assembly_gap"]),
    ("gap", None, fs.STATUS["assembly_gap"]),
    ("trace", "trace", fs.STATUS["tblastn_trace"]),
    ("absent", "controlled absence", fs.STATUS["absent"]),
    ("no_locus_unrescued", "not informative", fs.HILITE),
]


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    cells = {(c["species"], c["family"]): c["verdict"]
             for c in read_tsv(D / "r4_cells.tsv")}
    man = manifest()
    species = [s for s in man if any((s, f) in cells for f in FAMS)]
    col = {v: c for v, _, c in VERDICT}
    fig, ax = plt.subplots(figsize=(fs.W_HALF + 0.9, 7.4))
    for i, s in enumerate(species):
        for j, f in enumerate(FAMS):
            v = cells.get((s, f), "")
            ax.add_patch(plt.Rectangle((j, i), 0.92, 0.86, color=col.get(v, "white"),
                                       linewidth=0))
    prev = None
    for i, s in enumerate(species):
        g = man[s]["group"]
        if g != prev and i:
            ax.axhline(i - 0.07, color=fs.GRID, lw=0.6)
        prev = g
    ax.set_xlim(0, len(FAMS))
    ax.set_ylim(len(species), 0)
    ax.set_xticks([j + 0.46 for j in range(len(FAMS))], FAMS, rotation=45,
                  ha="left", fontsize=fs.FS_TICK - 0.6)
    ax.xaxis.tick_top()
    ax.set_yticks([i + 0.43 for i in range(len(species))],
                  [s if len(s) < 30 else s[:28] + "…" for s in species],
                  fontsize=fs.FS_TICK - 2.0, style="italic")
    for side in ("left", "right", "bottom", "top"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(handles=[Patch(color=c, label=l) for _, l, c in VERDICT if l],
              fontsize=fs.FS_NOTE - 0.6, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.005), ncol=2)
    fig.suptitle("Families added after S20, in the 52 genomes (S5 r4)",
                 fontsize=fs.FS_TITLE, x=0.02, ha="left", fontweight="normal")
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "genome_r4"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
