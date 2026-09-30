"""S15 figure — what each method contributes, per superfamily (Q3).

Two panels from the committed S15 tables (D13), through `figstyle.py`:

* **A** — the panel-frame curve: for each superfamily, the share of its
  final-census members first found by the domain *call*, by domain
  enumeration without a call, by the profiles alone, and by the genome sweep.
* **B** — the human frame: per superfamily, the curated human genes domain
  search enumerates, calls to the right family, and the profiles call.

    python3 scripts/s15_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402

D = ROOT / "results" / "method_contribution"
STEPS = [("domain_call", "domain search: right family", "#184f95"),
         ("domain_enumeration", "domain search: found, family not called", "#86b6ef"),
         ("profile", "profile HMM only", "#eb6834"),
         ("genome", "genome sweep only", "#1baf7a")]


def panel_a(ax, sf) -> None:
    for i, r in enumerate(sf):
        n, left = int(r["n"]), 0.0
        for key, lab, col in STEPS:
            v = int(r[f"n_{key}"]) / n
            ax.barh(i, v, left=left, color=col, height=0.72, linewidth=0,
                    label=lab if i == 0 else None)
            left += v
        ax.text(1.01, i, f"{n:,}", va="center", fontsize=fs.FS_NOTE - 0.6, color=fs.MUTED)
    ax.set_yticks(range(len(sf)), [fs.SUPERFAMILY_LABEL.get(r["superfamily"], r["superfamily"])
                                   for r in sf], fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.1)
    ax.set_xticks([0, 0.5, 1.0], ["0", "50 %", "100 %"])
    ax.set_xlabel("share of the superfamily's final-census members", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.07), ncol=2)
    fs.despine(ax)
    fs.panel(ax, "A", "Who finds each channel (panel frame)")


def panel_b(ax, sf) -> None:
    rows = [r for r in sf if int(r["human_genes"])]
    keys = [("human_enumerated", "enumerated by domain search", "#86b6ef"),
            ("human_domain_call", "right family by domain rules", "#184f95"),
            ("human_profile_call", "right family by profiles", "#eb6834")]
    h = 0.26
    for j, (k, lab, col) in enumerate(keys):
        ax.barh([i + (j - 1) * h for i in range(len(rows))],
                [int(r[k]) / int(r["human_genes"]) for r in rows], height=h,
                color=col, linewidth=0, label=lab)
    for i, r in enumerate(rows):
        ax.text(1.01, i, r["human_genes"], va="center", fontsize=fs.FS_NOTE - 0.6,
                color=fs.MUTED)
    ax.set_yticks(range(len(rows)), [fs.SUPERFAMILY_LABEL.get(r["superfamily"],
                                                              r["superfamily"])
                                     for r in rows], fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.1)
    ax.set_xticks([0, 0.5, 1.0], ["0", "50 %", "100 %"])
    ax.set_xlabel("share of the curated human genes", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.07), ncol=1)
    fs.despine(ax)
    fs.panel(ax, "B", "The human genes (independent frame)")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    sf = read_tsv(D / "curve_by_superfamily.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 5.4))
    panel_a(axes[0], sf)
    panel_b(axes[1], sf)
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "method_contribution"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
