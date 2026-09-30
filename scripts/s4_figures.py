"""S4 figure — the denominator: 52 panel species and the quality of what was searched.

Each species' reference proteome by BUSCO completeness (x) against its
assembly's scaffold N50 (y, log), coloured by panel group; genome-only
species (no proteome) marked on the axis. Drawn from
`results/proteome_scope/proteome_manifest.tsv`.

    python3 scripts/s4_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402

D = ROOT / "results" / "proteome_scope"
GROUPS = {"vertebrate": "#184f95", "deuterostome": "#6da7ec", "invertebrate": "#eb6834",
          "cnidarian": "#f3a683", "basal_metazoan": "#1baf7a", "holozoa": "#4a3aa7",
          "plant": "#0f7d3d", "prokaryote": "#8a897f", "virus": "#d6d5cf"}


def f(x: str) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    rows = read_tsv(D / "proteome_manifest.tsv")
    fig, ax = plt.subplots(figsize=(fs.W_FULL * 0.62, 3.8))
    seen = set()
    for r in rows:
        b, n = f(r["busco_complete_pct"]), f(r["scaffold_n50"])
        g = r["group"] if r["group"] in GROUPS else "fungi & protists"
        col = GROUPS.get(g, "#c9c3b0")
        if r["status"] != "proteome" or b is None or not n:
            continue
        ax.scatter(b, n, s=16, color=col, linewidth=0, alpha=0.9,
                   label=g if g not in seen else None)
        seen.add(g)
        if b < 85 or n < 2e5:
            ax.annotate(r["common"] or r["species"], (b, n), fontsize=fs.FS_NOTE - 1.2,
                        color=fs.MUTED, xytext=(3, 2), textcoords="offset points")
    genome_only = [r["species"] for r in rows if r["status"] != "proteome"]
    ax.set_yscale("log")
    ax.set_xlabel("BUSCO complete (%) of the reference proteome", fontsize=fs.FS_LABEL)
    ax.set_ylabel("assembly scaffold N50 (bp)", fontsize=fs.FS_LABEL)
    ax.axvline(90, color=fs.GRID, lw=0.6, ls="--")
    n_prot = sum(r["status"] == "proteome" for r in rows)
    ax.legend(fontsize=fs.FS_NOTE - 0.6, frameon=False, loc="center left",
              bbox_to_anchor=(1.0, 0.5))
    fs.despine(ax)
    fs.panel(ax, "", f"The denominator: {n_prot} proteomes + genome-only "
                     + ", ".join(g.split()[0] for g in genome_only))
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "panel"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
