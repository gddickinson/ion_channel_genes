"""S6 figure — what went into the alignments, and whether the pore modules are right.

Four panels, all drawn from the committed S6 tables (D13) through
`scripts/figstyle.py`:

* **A** — the D39 sets per superfamily: included (high-confidence profile
  calls, intact genome loci) against every exclusion reason.
* **B** — the family alignments: L-INS-i columns against columns kept by
  trimAl, per family (log scale), coloured by superfamily.
* **C** — module extraction checked against members' own UniProt topology
  (never a reference's): Jaccard overlap of extracted and annotated module,
  one dot per module, per family.
* **D** — the vote that fixes the module span of the six unannotated
  families, measured on the annotated ones with a held-out model: distance
  from the annotated span in profile states.

    python3 scripts/s6_figures.py
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

D = ROOT / "results" / "alignments"
REASONS = [("include", "included", "#184f95"),
           ("profile_medium", "medium-confidence call", "#86b6ef"),
           ("no_profile_call", "no profile call (S2 only)", "#a9a79e"),
           ("genome_not_intact", "genome locus not intact", "#ef9a90"),
           ("duplicate_of", "identical sequence", "#d6d5cf")]


def sf_colour(fam: str) -> str:
    sf = CATALOGUE[fam].superfamily
    return fs.SUPERFAMILY.get(sf, fs.SUPERFAMILY["other"])


def sf_label(sf: str) -> str:
    return fs.SUPERFAMILY_LABEL.get(sf, sf)


def panel_a(ax) -> None:
    c: dict[str, Counter] = defaultdict(Counter)
    for m in read_tsv(D / "members.tsv"):
        c[CATALOGUE[m["family"]].superfamily][m["verdict"].split(":")[0]] += 1
    order = sorted(c, key=lambda s: sum(c[s].values()))
    y = range(len(order))
    left = [0] * len(order)
    for key, lab, col in REASONS:
        v = [c[s][key] for s in order]
        ax.barh(list(y), v, left=left, color=col, height=0.72, label=lab,
                linewidth=0)
        left = [a + b for a, b in zip(left, v)]
    ax.set_yticks(list(y), [sf_label(s) for s in order], fontsize=fs.FS_TICK)
    ax.set_xlabel("census v4 rows called to a family", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.hgrid(ax, "x")
    fs.panel(ax, "A", "Alignment sets (D39)")


def panel_b(ax) -> None:
    rows = [r for r in read_tsv(D / "alignments.tsv") if r["status"] == "aligned"]
    if not rows:
        ax.text(0.5, 0.5, "alignments not yet run", ha="center",
                transform=ax.transAxes, color=fs.FAINT)
        fs.panel(ax, "B", "Family alignments")
        return
    for r in rows:
        ax.scatter(int(r["aln_cols"]), int(r["trim_cols"]), s=8 + int(r["included"]) / 8,
                   color=sf_colour(r["family"]), alpha=0.85, linewidth=0)
    lo, hi = 50, max(int(r["aln_cols"]) for r in rows) * 1.2
    ax.plot([lo, hi], [lo, hi], color=fs.FAINT, lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("L-INS-i columns", fontsize=fs.FS_LABEL)
    ax.set_ylabel("columns kept by trimAl", fontsize=fs.FS_LABEL)
    ax.text(0.03, 0.95, f"{len(rows)} families · dot area ~ sequences",
            transform=ax.transAxes, fontsize=fs.FS_NOTE, color=fs.MUTED, va="top")
    ax.legend(handles=fs.superfamily_handles(["ploop", "cysloop", "iglur", "other"]),
              fontsize=fs.FS_NOTE, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, 0.9))
    fs.despine(ax)
    fs.panel(ax, "B", "Family alignments")


def _families_in_unit_order() -> list[str]:
    return [r["family"] for r in read_tsv(D / "module_spans.tsv") if r["module"] in ("1", "")]


def panel_c(ax) -> None:
    by: dict[str, list[float]] = defaultdict(list)
    for r in read_tsv(D / "module_validation.tsv"):
        if r["status"] == "ok":
            by[r["family"]].append(float(r["jaccard"]))
    fams = [f for f in _families_in_unit_order() if f in by][::-1]
    for i, f in enumerate(fams):
        v = by[f]
        ax.scatter(v, [i] * len(v), s=6, color=sf_colour(f), alpha=0.6, linewidth=0)
    ax.axvline(0.8, color=fs.FAINT, lw=0.6, ls="--")
    ax.set_yticks(range(len(fams)), fams, fontsize=fs.FS_TICK - 0.8)
    ax.set_xlim(0, 1.03)
    ax.set_xlabel("Jaccard, extracted vs UniProt-annotated module", fontsize=fs.FS_LABEL)
    n = sum(map(len, by.values()))
    good = sum(x >= 0.8 for v in by.values() for x in v)
    ax.text(0.03, 0.02, f"{good}/{n} modules ≥ 0.8", transform=ax.transAxes,
            fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.despine(ax)
    fs.panel(ax, "C", "Module extraction, validated")


def panel_d(ax) -> None:
    err: dict[str, int] = {}
    vote_fams = []
    for r in read_tsv(D / "module_spans.tsv"):
        if r["basis"].startswith("unit_hmm_vote") and r["module"] == "1":
            vote_fams.append(r["family"])
        if r["vote_error"]:
            err[r["family"]] = max(err.get(r["family"], 0), int(r["vote_error"]))
    fams = [f for f in _families_in_unit_order() if f in err][::-1]
    cap = 30
    ax.barh(range(len(fams)), [min(err[f], cap) for f in fams],
            color=[sf_colour(f) for f in fams], height=0.7, linewidth=0)
    for i, f in enumerate(fams):
        if err[f] > cap:
            ax.text(cap - 0.5, i, f"{err[f]} →", ha="right", va="center",
                    fontsize=fs.FS_NOTE, color=fs.INK)
    ax.set_xlim(0, cap)
    ax.axvline(12, color=fs.ACCENT, lw=0.8, ls="--")
    ax.set_yticks(range(len(fams)), fams, fontsize=fs.FS_TICK - 0.8)
    ax.set_xlabel("held-out vote error (profile states)", fontsize=fs.FS_LABEL)
    ax.text(0.97, 0.98, "bound 12 states\n\nvote used for:\n" + "\n".join(vote_fams),
            transform=ax.transAxes, ha="right", va="top",
            fontsize=fs.FS_NOTE, color=fs.ACCENT)
    fs.despine(ax)
    fs.hgrid(ax, "x")
    fs.panel(ax, "D", "The span vote, measured")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(fs.W_FULL, 7.6),
                             gridspec_kw={"height_ratios": [1, 1.6]})
    panel_a(axes[0][0])
    panel_b(axes[0][1])
    panel_c(axes[1][0])
    panel_d(axes[1][1])
    fig.tight_layout()
    paths = fs.save(fig, D / "figures" / "alignments_modules")
    print("\n".join(map(str, paths)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
