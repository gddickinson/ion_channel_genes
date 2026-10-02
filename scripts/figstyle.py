"""One figure style for every publication figure in this project.

Import this and call `figstyle.use()` before building a figure. It fixes the
three things that made the figure set inconsistent:

**Size.** Figures used to be drawn 11-17 inches wide and then scaled to the
6.7-inch text block of the page, so 7 pt labels arrived at 3-4 pt. Everything is
now drawn at its final printed width (`W_FULL`), so a point is a point and the
PDF does no scaling.

**Colour.** One validated palette. Superfamily identity is categorical and
always the same hues; evidence, annotation quality and classification
confidence are *ordered*, so they use a diverging blue -> grey -> red scale
instead of more categorical hues. The values are the reference palette's slots
1/2/3/7 and its blue/red ramps; the set clears all-pairs CVD separation and
each arm of the diverging scale is a monotone single-hue ramp (checked with
the palette validator, not by eye). With twenty-five superfamilies in the
catalogue there are far more groups than safe hues, so most figures group or
facet rather than colour — see `SUPERFAMILY`.

**Chrome.** Hairline solid axes, no top/right spines, recessive grid, sans face
matching the manuscript's headings, and text saved as text (`pdf.fonttype 42`).
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import matplotlib
from matplotlib import font_manager

# ----------------------------------------------------------------- geometry

# Page geometry has one definition, in the manuscript module that also places
# the figures; importing it here keeps drawing and placement in step.
try:
    from s14_lib import H_MAX, W_FULL, W_HALF
except ImportError:                                # standalone use
    W_FULL, W_HALF, H_MAX = 6.7, 3.25, 8.6

# ------------------------------------------------------------------- colour

#: Categorical: superfamily identity. Reference palette slots 1, 2, 3 and the
#: accent — the four that clear the all-pairs CVD floor at 7 pt. **There are
#: twenty-five superfamilies and four safe categorical hues**, which is not a
#: palette problem, it is a figure-design constraint: a figure that needs to
#: distinguish more than about six groups is the wrong figure. Everything
#: outside the named set is `OTHER` grey, and a figure that greys out its
#: subject has to be redrawn, not recoloured.
SUPERFAMILY = {
    "ploop": "#2a78d6",         # blue   — the voltage-gated-like superfamily
    "cysloop": "#eb6834",       # orange — pentameric ligand-gated
    "iglur": "#1baf7a",         # aqua   — tetrameric ligand-gated
    "clc": "#4a3aa7",           # violet — the chloride channel/transporter fold
    "deg_enac": "#8a897f",      # warm grey
    "p2x": "#a9a79e",
    "tmem16_like": "#c9c3b0",
    "other": "#d6d5cf",
}
SUPERFAMILY_ORDER = ["ploop", "cysloop", "iglur", "clc", "deg_enac", "p2x",
                     "tmem16_like", "other"]
SUPERFAMILY_LABEL = {
    "ploop": "P-loop (VGIC)", "cysloop": "Cys-loop", "iglur": "iGluR",
    "clc": "CLC", "deg_enac": "DEG/ENaC", "p2x": "P2X",
    "tmem16_like": "TMEM16/OSCA/TMC", "other": "other superfamilies",
}

#: Categorical: what the channel conducts. Selectivity is a *literature*
#: attribute of a family, never a per-sequence prediction (hazard H14), so
#: this palette is only ever used to colour catalogue-level statements.
SELECTIVITY = {
    "K+": "#2a78d6", "Na+": "#eb6834", "Ca2+": "#1baf7a",
    "anion": "#4a3aa7", "cation_nonselective": "#8a897f",
    "H+": "#b3261e", "large_pore": "#a9a79e", "water": "#c9c3b0",
    "unknown": "#d6d5cf",
}
SELECTIVITY_ORDER = ["K+", "Na+", "Ca2+", "anion", "cation_nonselective",
                     "H+", "large_pore", "water", "unknown"]

#: Ordered: how confident a classification is. Same semantics as STATUS —
#: dark blue is best-evidenced, red is nothing — so the two scales can sit in
#: one figure without a reader having to relearn the direction.
CONFIDENCE = {
    "gold": "#184f95",       # two or more independent tiers agree
    "silver": "#3987e5",     # one tier calls it, nothing contradicts
    "bronze": "#d6d5cf",     # superfamily only, or a derived rule alone
    "unassigned": "#b3261e",
}
CONFIDENCE_ORDER = ["gold", "silver", "bronze", "unassigned"]

#: Diverging: how good the evidence for the gene is, present -> absent.
#: Blue arm (found), neutral greys (cannot tell), red arm (dead or gone).
STATUS = {
    "found_annotated": "#184f95",
    "found_unannotated": "#2a78d6",
    "found_no_annotation": "#86b6ef",
    "assembly_gap": "#d6d5cf",
    "tblastn_trace_ambiguous": "#a9a79e",
    "tblastn_trace": "#ef9a90",
    "absent": "#b3261e",
}
#: Reading order for stacked bars and legends: best evidence first.
STATUS_ORDER = [
    "found_annotated", "found_unannotated", "found_no_annotation",
    "assembly_gap", "tblastn_trace_ambiguous", "tblastn_trace", "absent",
]
#: Short labels — the raw status strings are too long for a legend at 7 pt.
STATUS_LABEL = {
    "found_annotated": "found, annotated",
    "found_unannotated": "found, not annotated",
    "found_no_annotation": "found, no gene set",
    "assembly_gap": "assembly gap",
    "tblastn_trace_ambiguous": "trace, ambiguous",
    "tblastn_trace": "remnant",
    "absent": "absent",
}

#: Diverging: annotation quality, same semantic direction as STATUS.
QUALITY = {
    "complete": "#184f95",
    "split": "#86b6ef",
    "fragmentary": "#d6d5cf",
    "noncoding": "#ef9a90",
    "unannotated": "#b3261e",
}

#: Clinical significance — a status encoding, never reused for a series.
#: The channelopathies are the largest Mendelian disease class attached to
#: any protein family — long-QT, Dravet, cystic fibrosis, myotonia,
#: malignant hyperthermia, deafness, polycystic kidney disease — so this
#: scale carries weight rather than being decorative.
CLINICAL = {"pathogenic": "#b3261e", "uncertain": "#a9a79e",
            "benign": "#0f7d3d"}

#: Sequential blue ramp (light -> dark) for magnitude.
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]

INK = "#14140f"        # primary text
MUTED = "#52514e"      # secondary text
FAINT = "#8a897f"      # annotation, de-emphasised marks
GRID = "#e5e4df"       # gridlines and hairline rules
SURFACE = "#ffffff"
ACCENT = "#4a3aa7"     # callouts, highlight boxes (violet, slot 7)
HILITE = "#f2f1ec"     # region shading behind a highlighted span

# ------------------------------------------------------------------ typography

FS_TICK = 6.6
FS_LABEL = 7.4
FS_TITLE = 8.2         # panel title
FS_LETTER = 9.0        # panel letter
FS_SUPTITLE = 9.4
FS_NOTE = 6.4          # in-plot annotation

_FONT_STACK = ["TeX Gyre Heros", "Helvetica Neue", "Helvetica", "Arial",
               "DejaVu Sans"]


def _register_document_sans() -> None:
    """Make the manuscript's sans face available to matplotlib.

    The PDF's headings are TeX Gyre Heros, which lives in the TeX tree and is
    not on matplotlib's font path. Registering it means figure text and heading
    text are the same typeface; if TeX is absent the stack falls back.
    """
    try:
        out = subprocess.run(
            ["kpsewhich", "texgyreheros-regular.otf", "texgyreheros-bold.otf",
             "texgyreheros-italic.otf", "texgyreheros-bolditalic.otf"],
            capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return
    for line in out.stdout.split("\n"):
        path = line.strip()
        if path and Path(path).exists():
            try:
                font_manager.fontManager.addfont(path)
            except Exception:                       # noqa: BLE001 - optional
                pass


_REGISTERED = False


def use() -> None:
    """Apply the project figure style. Safe to call more than once."""
    global _REGISTERED
    if not _REGISTERED:
        _register_document_sans()
        _REGISTERED = True
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": _FONT_STACK,
        "font.size": FS_LABEL,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.titlepad": 4.0,
        "axes.labelsize": FS_LABEL,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.labelpad": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "grid.linestyle": "-",
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.major.pad": 1.8,
        "ytick.major.pad": 1.8,
        "legend.fontsize": FS_TICK,
        "legend.frameon": False,
        "legend.handlelength": 1.1,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.0,
        "legend.labelspacing": 0.35,
        "legend.borderaxespad": 0.2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.0,
        "patch.linewidth": 0.0,
        "figure.titlesize": FS_SUPTITLE,
        "figure.titleweight": "bold",
        "figure.dpi": 120,
        "figure.facecolor": SURFACE,
        "savefig.dpi": 400,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "text.color": INK,
    })


# ------------------------------------------------------------------- helpers

def panel(ax, letter: str, title: str = "", pad: float = 4.0) -> None:
    """Bold panel letter at the top-left, description in roman beside it."""
    ax.set_title("  " + title if title else "", loc="left", pad=pad,
                 fontsize=FS_TITLE, color=INK)
    ax.annotate(letter, xy=(0.0, 1.0), xycoords="axes fraction",
                xytext=(-1.0, pad + 1.0), textcoords="offset points",
                fontsize=FS_LETTER, fontweight="bold", color=INK,
                ha="right", va="baseline", annotation_clip=False)


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def hgrid(ax, axis: str = "y") -> None:
    """A recessive solid grid, drawn behind the marks."""
    ax.set_axisbelow(True)
    ax.grid(True, axis=axis, color=GRID, linewidth=0.5, linestyle="-")


def superfamily_handles(labels=None):
    """Legend handles for the named superfamilies, in fixed order.

    Anything not in `SUPERFAMILY` gets the `other` grey — see the palette
    note: there are twenty-five superfamilies and four safe categorical
    hues, so a figure that needs more groups than this needs faceting, not
    a bigger palette.
    """
    from matplotlib.patches import Patch
    labels = labels or SUPERFAMILY_ORDER
    return [Patch(facecolor=SUPERFAMILY.get(k, SUPERFAMILY["other"]),
                  label=SUPERFAMILY_LABEL.get(k, k)) for k in labels]


def confidence_handles(labels=None):
    """Legend handles for the classification confidence scale."""
    from matplotlib.patches import Patch
    labels = labels or CONFIDENCE_ORDER
    return [Patch(facecolor=CONFIDENCE[k], label=k) for k in labels]


def check_titles(fig, name: str, legend_in: str = "") -> None:
    """Every figure is titled: a figure title, or a title on every panel.

    Project rule (CLAUDE.md, *Figures*): a reader must know what a figure is
    without hunting for its caption. A colour bar is not a panel. The one
    exception is a figure whose title and legend live in a document beside
    it (the review's numbered figures) — `legend_in` names that document,
    so the exception is stated at the call, never silent.
    """
    if legend_in:
        return
    sup = fig._suptitle.get_text().strip() if fig._suptitle else ""
    panels = [ax for ax in fig.axes if ax.get_visible() and ax.get_label() != "<colorbar>"]
    untitled = [ax for ax in panels
                if not any(ax.get_title(loc=l).strip() for l in ("left", "center", "right"))]
    if not sup and untitled:
        raise RuntimeError(
            f"{name}: {len(untitled)} of {len(panels)} panel(s) have no title and "
            f"the figure has none — give each panel a title (figstyle.panel) or the "
            f"figure a suptitle; pass legend_in= only if the legend lives elsewhere")


def save(fig, stem, formats=("png", "pdf"), dpi: int = 400,
         legend_in: str = "") -> list[Path]:
    """Save one figure to every format, and report if it exceeds the page.

    Refuses an untitled figure (`check_titles`); `legend_in` names the
    document holding the legend of a figure deliberately drawn without one.
    """
    stem = Path(stem)
    check_titles(fig, stem.name, legend_in)
    stem.parent.mkdir(parents=True, exist_ok=True)
    w, h = fig.get_size_inches()
    if w > W_FULL + 0.01:
        print(f"  [figstyle] {stem.name}: {w:.2f} in wide — wider than the "
              f"{W_FULL} in text block, it will be scaled down in the PDF")
    if h > H_MAX:
        print(f"  [figstyle] {stem.name}: {h:.2f} in tall — taller than "
              f"{H_MAX} in, it will be scaled down in the PDF")
    # A single artist placed in the wrong coordinate system can push the
    # tight bounding box hundreds of axes-heights off-canvas, and savefig
    # will happily write a several-hundred-megapixel file. Catch it here.
    try:
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        if tb.width > 3 * w or tb.height > 3 * h:
            raise RuntimeError(
                f"{stem.name}: tight bbox is {tb.width:.1f}x{tb.height:.1f} in "
                f"for a {w:.1f}x{h:.1f} in figure — an artist is drawn far "
                f"off-canvas (check transform= on any annotation)")
    except RuntimeError:
        raise
    except Exception:                              # noqa: BLE001 - best effort
        pass
    out = []
    for ext in formats:
        path = stem.with_suffix(f".{ext}")
        # A character the figure font lacks is dropped from the output with
        # only a warning, so a label can silently lose glyphs. Promote it.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fig.savefig(path, dpi=dpi)
        missing = {str(w.message) for w in caught
                   if "missing from font" in str(w.message)}
        if missing:
            raise RuntimeError(f"{stem.name}: {len(missing)} glyph(s) have no "
                               f"outline in the figure font and would be "
                               f"dropped — {sorted(missing)[0]}")
        out.append(path)
    return out
