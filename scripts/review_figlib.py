"""Drawing primitives shared by the review's figures.

Everything here is a *drawing* helper — sequence rows, membrane cartoons,
domain bars. The data always comes from a committed table; nothing in this
module fetches, computes or decides anything.

Two conventions the whole figure set follows:

**A schematic is labelled a schematic.** The topology cartoons are drawings
of what the catalogue records (helix count, pore loops, stoichiometry), not
renderings of structures, and every panel that is one says so in its
caption. Mixing a schematic and a measurement in one figure without marking
which is which is the fastest way to make a review misleading.

**Residues are coloured by chemistry, one scheme, everywhere.** The
alignment figures colour by residue class rather than by conservation, so
the same amino acid is the same colour in every panel and a reader can
compare across figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import figstyle as fs

#: Residue classes. Muted enough that a whole alignment does not vibrate,
#: separated enough to read at 6 pt.
RESIDUE_CLASS = {
    **{a: "hydrophobic" for a in "AVLIMFWCP"},
    **{a: "polar" for a in "STNQGY"},
    **{a: "acidic" for a in "DE"},
    **{a: "basic" for a in "KRH"},
}
CLASS_COLOUR = {
    "hydrophobic": "#e8e6df",
    "polar": "#cde2fb",
    "acidic": "#f6c7c2",
    "basic": "#cfe6d5",
    "gap": "#ffffff",
}
CLASS_LABEL = {
    "hydrophobic": "hydrophobic (AVLIMFWCP)",
    "polar": "polar (STNQGY)",
    "acidic": "acidic (DE)",
    "basic": "basic (KRH)",
}


#: Short display names for the superfamilies. The catalogue's names are
#: precise and too long for an axis; these are used in every figure so a
#: reader tracks the same label across panels.
SUPERFAMILY_SHORT = {
    "ploop": "P-loop (VGIC)", "cysloop": "Cys-loop (pLGIC)",
    "iglur": "iGluR", "tmem16_like": "TMEM16 / OSCA / TMC",
    "connexin": "connexins", "deg_enac": "DEG / ENaC",
    "innexin_like": "innexin / pannexin / LRRC8", "p2x": "P2X",
    "ca_release": "ITPR / RYR", "calhm": "CALHM", "clic": "CLIC",
    "bestrophin": "bestrophins", "clc": "CLC", "orai": "ORAI",
    "otopetrin": "otopetrins", "porin": "VDAC / porins",
    "tweety": "tweety", "mcu": "MCU", "piezo": "Piezo", "tric": "TRIC",
    "abc_channel": "CFTR (ABC fold)", "hv": "Hv1", "tmem175": "TMEM175",
    "msc": "MscL / MscS", "viroporin": "viroporins", "": "(no superfamily)",
}


def superfamily_label(key: str) -> str:
    from src.catalogue import SUPERFAMILIES
    if key in SUPERFAMILY_SHORT:
        return SUPERFAMILY_SHORT[key]
    sf = SUPERFAMILIES.get(key)
    if sf is None:
        return key or "(no superfamily)"
    name = sf.name.split("(")[0].strip()
    return (name[:26] + "…") if len(name) > 27 else name


def residue_colour(ch: str) -> str:
    return CLASS_COLOUR.get(RESIDUE_CLASS.get(ch.upper(), ""), "#ffffff")


def draw_sequence_row(ax, y: float, seq: str, x0: float = 0.0,
                      cell: float = 1.0, highlight: tuple[int, int] | None = None,
                      fontsize: float | None = None) -> None:
    """One row of residues as coloured cells with the letter on top."""
    from matplotlib.patches import Rectangle
    fontsize = fontsize or fs.FS_TICK - 0.8
    for i, ch in enumerate(seq):
        x = x0 + i * cell
        inside = highlight and highlight[0] <= i < highlight[1]
        ax.add_patch(Rectangle((x, y), cell, 1.0,
                               facecolor=residue_colour(ch),
                               edgecolor="#ffffff", linewidth=0.4, zorder=2))
        if inside:
            ax.add_patch(Rectangle((x, y), cell, 1.0, facecolor="none",
                                   edgecolor=fs.INK, linewidth=0.9, zorder=4))
        ax.text(x + cell / 2, y + 0.5, ch, ha="center", va="center",
                fontsize=fontsize, color=fs.INK, zorder=5,
                fontweight="bold" if inside else "normal")


def residue_legend(ax, loc: str = "lower right", ncol: int = 4) -> None:
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=CLASS_COLOUR[k], edgecolor="#ffffff",
                             label=CLASS_LABEL[k])
                       for k in ("hydrophobic", "polar", "acidic", "basic")],
              fontsize=fs.FS_TICK - 1.2, frameon=False, loc=loc, ncol=ncol,
              handlelength=1.0, handleheight=0.9, columnspacing=1.0)


def membrane(ax, x0: float, x1: float, y0: float = 0.0, y1: float = 1.0) -> None:
    """The lipid bilayer as a pale band — the reference frame for topology."""
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#f2f1ec",
                           edgecolor="#e5e4df", linewidth=0.6, zorder=1))


def helix(ax, x: float, y0: float, y1: float, width: float = 0.34,
          colour: str | None = None, label: str = "") -> None:
    """One transmembrane helix."""
    from matplotlib.patches import FancyBboxPatch
    ax.add_patch(FancyBboxPatch((x - width / 2, y0), width, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=0.12",
                                facecolor=colour or fs.BLUES[2],
                                edgecolor=fs.MUTED, linewidth=0.5, zorder=3))
    if label:
        ax.text(x, y1 + 0.08, label, ha="center", va="bottom",
                fontsize=fs.FS_TICK - 1.4, color=fs.MUTED, zorder=4)


def pore_loop(ax, x0: float, x1: float, y_top: float, depth: float = 0.45,
              colour: str | None = None) -> None:
    """A re-entrant loop: in from the top, down, and back out."""
    import numpy as np
    xs = np.linspace(x0, x1, 40)
    mid = (x0 + x1) / 2
    ys = y_top - depth * np.exp(-((xs - mid) / ((x1 - x0) / 3.2)) ** 2)
    ax.plot(xs, ys, color=colour or fs.ACCENT, linewidth=1.5,
            solid_capstyle="round", zorder=4)


def strand(ax, x: float, y0: float, y1: float, width: float = 0.26,
           colour: str | None = None) -> None:
    """One β-strand — flat, so a barrel is never drawn as a helix bundle."""
    from matplotlib.patches import Polygon
    c = colour or fs.SUPERFAMILY["clc"]
    ax.add_patch(Polygon([[x - width / 2, y0], [x + width / 2, y0],
                          [x + width / 2, y1 - 0.10], [x, y1],
                          [x - width / 2, y1 - 0.10]],
                         closed=True, facecolor=c, edgecolor=fs.MUTED,
                         linewidth=0.4, zorder=3))


def top_view(ax, cx: float, cy: float, n: int, r: float = 0.30,
             colour: str | None = None, pore: bool = True) -> None:
    """Stoichiometry as a top view: n subunits around a central pore."""
    import numpy as np
    from matplotlib.patches import Circle
    for k in range(n):
        a = 2 * np.pi * k / n - np.pi / 2
        ax.add_patch(Circle((cx + r * np.cos(a), cy + r * np.sin(a)), r * 0.42,
                            facecolor=colour or fs.BLUES[2],
                            edgecolor=fs.MUTED, linewidth=0.4, zorder=3))
    if pore:
        ax.add_patch(Circle((cx, cy), r * 0.20, facecolor="#ffffff",
                            edgecolor=fs.MUTED, linewidth=0.4, zorder=4))


def domain_bar(ax, y: float, length: int, domains: list[tuple],
               height: float = 0.62, colours: dict | None = None,
               label_min_frac: float = 0.055, scale: int | None = None) -> None:
    """A protein as a scale bar with its domains drawn where they occur.

    `domains` is `[(start, end, short_name, accession)]`, 1-based inclusive,
    exactly as `domain_positions.tsv` stores them.

    `scale` is the width of the *panel*, not of this protein, and it is what
    the label threshold is measured against. Using the protein's own length
    labels a 150-residue domain on a 580-residue protein drawn beside a
    4,303-residue one, where the box is a few pixels wide and the text lands
    on its neighbours.
    """
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((0, y + height / 2 - 0.045), length, 0.09,
                           facecolor="#dddcd8", edgecolor="none", zorder=2))
    for start, end, short, acc in domains:
        c = (colours or {}).get(acc, fs.BLUES[2])
        ax.add_patch(Rectangle((start, y), end - start, height,
                               facecolor=c, edgecolor="#ffffff",
                               linewidth=0.5, zorder=3))
        if (end - start) / max(1, scale or length) >= label_min_frac:
            ax.text((start + end) / 2, y + height / 2, short,
                    ha="center", va="center", fontsize=fs.FS_TICK - 1.6,
                    color=fs.INK, zorder=4)
