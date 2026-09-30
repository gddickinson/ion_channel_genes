"""S3a figure — the profile library, calibrated against S2 and put to work.

* **A** — per family, how often the S3a profile call agrees with S2's domain
  call where both call (seeds excluded), against the number of S2 calls.
* **B** — S2's superfamily-only records per superfamily, by what the
  profiles made of them: a family call, still superfamily-only, or other.

Drawn from `results/census_v3/calibration.tsv` and
`superfamily_only_resolved.tsv`.

    python3 scripts/s3_figures.py
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

D = ROOT / "results" / "census_v3"


def panel_a(ax) -> None:
    rows = [r for r in read_tsv(D / "calibration.tsv") if int(r["v2_calls"]) >= 20]
    for r in rows:
        sf = CATALOGUE[r["family"]].superfamily if r["family"] in CATALOGUE else "other"
        ax.scatter(int(r["v2_calls"]), float(r["agreement"]), s=12, linewidth=0, alpha=0.85,
                   color=fs.SUPERFAMILY.get(sf, fs.SUPERFAMILY["other"]))
        if float(r["agreement"]) < 0.995:
            ax.annotate(r["family"], (int(r["v2_calls"]), float(r["agreement"])),
                        fontsize=fs.FS_NOTE - 0.8, color=fs.MUTED, xytext=(3, -3),
                        textcoords="offset points")
    ax.set_xscale("log")
    ax.set_ylim(min(float(r["agreement"]) for r in rows) - 0.003, 1.001)
    ax.axhline(1, color=fs.GRID, lw=0.6)
    agree = sum(int(r["agree"]) for r in rows)
    both = agree + sum(int(r["disagree"]) for r in rows)
    ax.text(0.03, 0.05, f"agreement where both call: {agree / both:.1%}",
            transform=ax.transAxes, fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.set_xlabel("S2 domain-rule calls (seeds excluded)", fontsize=fs.FS_LABEL)
    ax.set_ylabel("profile agrees", fontsize=fs.FS_LABEL)
    ax.legend(handles=fs.superfamily_handles(["ploop", "cysloop", "iglur", "other"]),
              fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "A", "Profiles calibrated against S2")


def panel_b(ax) -> None:
    by = defaultdict(Counter)
    for r in read_tsv(D / "superfamily_only_resolved.tsv"):
        o = r["v3_outcome"]
        k = ("still superfamily-only" if o == "superfamily_only" else
             "other / unassigned" if o in ("unassigned", "conflict", "ambiguous", "")
             else "family call")
        by[r["v2_superfamily"]][k] += int(r["records"])
    sfs = sorted(by, key=lambda s: -sum(by[s].values()))[:10]
    cols = [("family call", "#184f95"), ("still superfamily-only", "#86b6ef"),
            ("other / unassigned", "#d6d5cf")]
    for i, s in enumerate(sfs):
        left, tot = 0, sum(by[s].values())
        for k, c in cols:
            v = by[s][k] / tot
            ax.barh(i, v, left=left, color=c, height=0.7, linewidth=0,
                    label=k if i == 0 else None)
            left += v
        ax.text(1.01, i, f"{tot:,}", va="center", fontsize=fs.FS_NOTE - 0.6, color=fs.MUTED)
    ax.set_yticks(range(len(sfs)), [fs.SUPERFAMILY_LABEL.get(s, s) for s in sfs],
                  fontsize=fs.FS_TICK - 0.8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.5, 1], ["0", "50 %", "100 %"])
    ax.set_xlabel("S2 superfamily-only records", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.14), ncol=2)
    fs.despine(ax)
    fs.panel(ax, "B", "What the profiles resolved")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.7))
    panel_a(axes[0])
    panel_b(axes[1])
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "profiles"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
