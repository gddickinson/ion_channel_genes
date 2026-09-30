"""S1 figure — the classifier benchmark: who made each call, and the hazards.

* **A** — the 72 positive-panel proteins by the tier that decided a correct
  call (architecture/hazard rule, filter motif, reference identity) or by
  failure, per superfamily. The question S1 had to answer: is the classifier
  a nearest-neighbour lookup? (No — most correct calls are not reference.)
* **B** — every hazard: panel proteins it touches, called right and wrong.

Drawn from `results/benchmark_controls/calls.tsv` and `hazards.tsv`.

    python3 scripts/s1_figures.py
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

D = ROOT / "results" / "benchmark_controls"
TIERS = [("architecture", "architecture / hazard rule", "#184f95"),
         ("motif", "selectivity-filter motif", "#1baf7a"),
         ("reference", "reference identity", "#86b6ef"),
         ("wrong", "wrong family", "#b3261e"),
         ("none", "no family call", "#d6d5cf")]


def outcome(r: dict) -> str:
    if r["called_family"] != r["expected_family"]:
        return "wrong" if r["called_family"] else "none"
    t = r["decisive_tier"]
    return "architecture" if t in ("architecture", "hazard") else t or "none"


def panel_a(ax) -> None:
    rows = [r for r in read_tsv(D / "calls.tsv")
            if r["expected_family"] in CATALOGUE and CATALOGUE[r["expected_family"]].census_member()]
    by = defaultdict(Counter)
    for r in rows:
        by[CATALOGUE[r["expected_family"]].superfamily][outcome(r)] += 1
    sfs = sorted(by, key=lambda s: -sum(by[s].values()))
    for i, sf in enumerate(sfs):
        left = 0
        for k, lab, col in TIERS:
            ax.barh(i, by[sf][k], left=left, color=col, height=0.7, linewidth=0,
                    label=lab if i == 0 else None)
            left += by[sf][k]
    ax.set_yticks(range(len(sfs)), [fs.SUPERFAMILY_LABEL.get(s, s) for s in sfs],
                  fontsize=fs.FS_TICK - 0.8)
    ax.invert_yaxis()
    tot = Counter(outcome(r) for r in rows)
    ax.set_xlabel("positive-panel proteins", fontsize=fs.FS_LABEL)
    ax.text(0.97, 0.03, f"correct {tot['architecture'] + tot['motif'] + tot['reference']}"
            f"/{len(rows)}: rules {tot['architecture']}, motif {tot['motif']}, "
            f"reference {tot['reference']}", transform=ax.transAxes, ha="right",
            fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.legend(fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.1), ncol=2)
    fs.despine(ax)
    fs.panel(ax, "A", "Which tier made each call")


def panel_b(ax) -> None:
    rows = read_tsv(D / "hazards.tsv")
    y = range(len(rows))
    ok = [int(r["n_correct"]) for r in rows]
    bad = [int(r["n_wrong"]) for r in rows]
    ax.barh(list(y), ok, color="#184f95", height=0.7, linewidth=0, label="called right")
    ax.barh(list(y), bad, left=ok, color="#b3261e", height=0.7, linewidth=0,
            label="called wrong")
    ax.set_yticks(list(y), [r["hazard"] for r in rows], fontsize=fs.FS_TICK - 0.8)
    ax.invert_yaxis()
    ax.set_xlabel("panel proteins the hazard touches", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "B", f"All {len(rows)} hazards exercised")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 4.4),
                             gridspec_kw={"width_ratios": [1.3, 1]})
    panel_a(axes[0])
    panel_b(axes[1])
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "benchmark"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
