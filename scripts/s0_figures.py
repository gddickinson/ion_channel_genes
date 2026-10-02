"""S0 figures — the catalogue's shape and the evidence behind the hazards.

Two panels, both drawn from the committed S0 tables and nothing else (D13),
both going through the one figure style (`scripts/figstyle.py`).

**Panel A — what the catalogue counts.** Families and human genes per
superfamily, census families solid and control families hatched. The point
of the panel is the lopsidedness: the P-loop superfamily is 143 of the 320
human genes and twenty-four other superfamilies share the rest.

**Panel B — the shared signatures.** Every domain signature carried by more
than one family, ordered by how many, with the ones that cross the
channel / non-channel boundary marked. This is the hazard registry's
evidence in one picture: `PF00520` reaches twenty families including a
phosphatase.

    python3 scripts/s0_figures.py [--dir results/s0_baseline]
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs
import review_figlib as R
from scripts.s0_lib import read_tsv
from src.catalogue import CATALOGUE, SUPERFAMILIES, census_families


def panel_a(ax) -> None:
    census = Counter()
    control = Counter()
    genes = Counter()
    for f in CATALOGUE.values():
        (census if f.census_member() else control)[f.superfamily] += 1
        if f.census_member():
            genes[f.superfamily] += len(f.human_genes)

    order = sorted(set(census) | set(control),
                   key=lambda k: (-genes.get(k, 0), -census.get(k, 0), k))
    y = range(len(order))
    colours = [fs.SUPERFAMILY.get(k, fs.SUPERFAMILY["other"]) for k in order]

    ax.barh(list(y), [genes.get(k, 0) for k in order], color=colours,
            height=0.68, zorder=3)
    for i, k in enumerate(order):
        n_c, n_x = census.get(k, 0), control.get(k, 0)
        label = f"{n_c}" + (f"+{n_x}" if n_x else "")
        ax.text(genes.get(k, 0) + 2, i, label,
                va="center", fontsize=fs.FS_TICK - 0.6, color=fs.MUTED)
    ax.set_yticks(list(y))
    ax.set_yticklabels([R.superfamily_label(k) for k in order],
                       fontsize=fs.FS_TICK - 0.4)
    ax.invert_yaxis()
    ax.set_xlabel("human pore-forming genes  (right: census + control families)",
                  fontsize=fs.FS_LABEL - 0.6)
    ax.set_xlim(0, max(genes.values()) * 1.12)
    fs.despine(ax)
    fs.hgrid(ax, axis="x")


def _unused_name(key: str) -> str:
    if key in SHORT:
        return SHORT[key]
    sf = SUPERFAMILIES.get(key)
    if sf is None:
        return key or "(no superfamily)"
    name = sf.name.split("(")[0].strip()
    return (name[:26] + "…") if len(name) > 27 else name


def panel_b(ax, share_rows: list[dict]) -> None:
    rows = [r for r in share_rows if int(r["n_families"]) >= 2][:16]
    y = range(len(rows))
    crossing = [r["crosses_channel_boundary"] == "yes" for r in rows]
    colours = [fs.STATUS["absent"] if c else fs.BLUES[3] for c in crossing]
    ax.barh(list(y), [int(r["n_families"]) for r in rows], color=colours,
            height=0.68, zorder=3)
    ax.set_yticks(list(y))
    ax.set_yticklabels([r["accession"] for r in rows], fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlabel("catalogued families carrying the signature",
                  fontsize=fs.FS_LABEL)
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor=fs.STATUS["absent"],
              label="also carried by a non-channel"),
        Patch(facecolor=fs.BLUES[3], label="channels only")],
        fontsize=fs.FS_TICK - 0.4, frameon=False, loc="lower right")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", type=Path, default=ROOT / "results" / "s0_baseline")
    args = ap.parse_args()
    share = read_tsv(args.dir / "signature_sharing.tsv")
    if not share:
        print(f"[s0_figures] no signature_sharing.tsv in {args.dir} — run "
              f"scripts/s0_catalogue_verify.py first")
        return 1

    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 4.0),
                             gridspec_kw={"width_ratios": [1.18, 1.0]})
    panel_a(axes[0])
    panel_b(axes[1], share)
    fs.panel(axes[0], "A", "Human pore-forming genes per superfamily")
    fs.panel(axes[1], "B", "Domains carried by more than one family")
    fig.tight_layout()
    out = args.dir / "figures" / "catalogue_scope"
    paths = fs.save(fig, out)
    n_census = len(census_families())
    print(f"[s0_figures] {n_census} census families drawn → "
          f"{', '.join(str(p) for p in paths)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
