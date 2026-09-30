"""Census revision r4 figure — what the new families found, and whether their
look-alikes separate.

* **A** — the 26,783 r4 records by their census v3a call and confidence.
* **B** — every reviewed r4 record's D32 margin over the runner-up profile,
  per called family; the look-alike pairs (H17–H20) sit side by side, and a
  margin under 0.30 is a medium-confidence call.

Drawn from `results/census_v3/r4_new_records.tsv` and `r4_reviewed_calls.tsv`.

    python3 scripts/s3r4_figures.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402

D = ROOT / "results" / "census_v3"
CONF = [("high", "#184f95"), ("medium", "#86b6ef"), ("low", "#d6d5cf"), ("none", "#e5e4df")]
PAIRS = ["pacc", "clcc1", "gphr", "tmco1", "nonchannel_emc3", "tmem109",
         "nonchannel_bri3bp", "tmem87", "nonchannel_tmem87b", "nonchannel_gost"]


def panel_a(ax) -> None:
    by = defaultdict(lambda: defaultdict(int))
    for r in read_tsv(D / "r4_new_records.tsv"):
        by[r["v3_call"]][r["p_confidence"] or "none"] += int(r["records"])
    fams = sorted(by, key=lambda f: -sum(by[f].values()))
    for i, f in enumerate(fams):
        left = 0
        for c, col in CONF:
            v = by[f][c]
            ax.barh(i, v, left=left, color=col, height=0.7, linewidth=0,
                    label=c if i == 0 else None)
            left += v
        ax.text(left + 60, i, f"{left:,}", va="center", fontsize=fs.FS_NOTE - 0.6,
                color=fs.MUTED)
    ax.set_yticks(range(len(fams)), [f.replace("nonchannel_", "decoy: ") for f in fams],
                  fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlim(0, max(sum(v.values()) for v in by.values()) * 1.22)
    ax.set_xlabel("r4 records (census v3a call)", fontsize=fs.FS_LABEL)
    ax.legend(title="profile confidence", title_fontsize=fs.FS_NOTE, fontsize=fs.FS_NOTE,
              frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "A", "What r4 added")


def panel_b(ax) -> None:
    rows = [r for r in read_tsv(D / "r4_reviewed_calls.tsv") if r["p_family"] in PAIRS]
    y = {f: i for i, f in enumerate(PAIRS)}
    for r in rows:
        m = float(r["rel_margin"]) if r["rel_margin"] else 1.0
        ax.scatter(m, y[r["p_family"]], s=12, linewidth=0.6,
                   facecolor="white" if r["seed"] == "1" else fs.BLUES[5],
                   edgecolor=fs.BLUES[5], zorder=3)
    ax.axvline(0.30, color=fs.ACCENT, lw=0.7, ls="--")
    ax.text(0.31, len(PAIRS) - 0.4, "high ≥ 0.30", fontsize=fs.FS_NOTE, color=fs.ACCENT)
    for k, (a, b) in enumerate([(3, 4), (5, 6), (7, 9)]):
        ax.axhspan(a - 0.45, b + 0.45, color=fs.HILITE, zorder=0)
        ax.text(1.02, (a + b) / 2, ["H17", "H19", "H18"][k], va="center",
                fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.set_yticks(range(len(PAIRS)), [p.replace("nonchannel_", "decoy: ") for p in PAIRS],
                  fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.1)
    ax.set_xlabel("margin over the runner-up profile (reviewed records)", fontsize=fs.FS_LABEL)
    ax.text(0.02, 0.98, "open: seed\nfilled: held out", transform=ax.transAxes, va="top",
            fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.despine(ax)
    fs.panel(ax, "B", "The look-alikes (H17–H19)")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.6))
    panel_a(axes[0])
    panel_b(axes[1])
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "census_r4"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
