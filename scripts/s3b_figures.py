"""S3b figure — what domain search missed across the 50 panel proteomes.

* **A** — census families with high-confidence panel members that no
  enumerated pore signature reaches (outside census v2), by family.
* **B** — the same misses by panel group, as a share of each group's
  channel-family members.

Drawn from `results/panel_sweep/domain_search_missed.tsv` and
`missed_by_group.tsv`.

    python3 scripts/s3b_figures.py
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

D = ROOT / "results" / "panel_sweep"
CHANNEL = ("channel", "channel_contested")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    rows = [r for r in read_tsv(D / "domain_search_missed.tsv")
            if r["status"] in CHANNEL and int(r["missed_high"])]
    rows.sort(key=lambda r: -int(r["missed_high"]))
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.9),
                             gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    ax.barh(range(len(rows)), [int(r["missed_high"]) for r in rows],
            color=[fs.SUPERFAMILY.get(r["superfamily"], fs.SUPERFAMILY["other"]) for r in rows],
            height=0.7, linewidth=0)
    for i, r in enumerate(rows):
        ax.text(int(r["missed_high"]) + 0.4, i, f"{float(r['missed_frac']):.0%} of "
                f"{r['called']}", va="center", fontsize=fs.FS_NOTE - 1.0, color=fs.MUTED)
    ax.set_yticks(range(len(rows)), [r["family"] for r in rows], fontsize=fs.FS_TICK - 1.0)
    ax.invert_yaxis()
    ax.set_xlim(0, max(int(r["missed_high"]) for r in rows) * 1.45)
    ax.set_xlabel("high-confidence members outside census v2", fontsize=fs.FS_LABEL)
    ax.legend(handles=fs.superfamily_handles(["ploop", "cysloop", "iglur", "other"]),
              fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "A", "Channel members no domain search reaches")
    ax = axes[1]
    tot, miss = Counter(), Counter()
    for r in read_tsv(D / "missed_by_group.tsv"):
        # census families only: controls (aquaporins, SLC26, gasdermins …)
        # carry no pore signature by design, so "missed" means nothing there
        if not (r["family"] in CATALOGUE and CATALOGUE[r["family"]].census_member()):
            continue
        tot[r["group"]] += int(r["in_census_v2"]) + int(r["missed_by_domain_search"])
        miss[r["group"]] += int(r["missed_by_domain_search"])
    groups = sorted((g for g in tot if tot[g] >= 20), key=lambda g: -miss[g] / tot[g])
    ax.barh(range(len(groups)), [miss[g] / tot[g] for g in groups], color="#184f95",
            height=0.7, linewidth=0)
    for i, g in enumerate(groups):
        ax.text(miss[g] / tot[g] + 0.003, i, f"{miss[g]}/{tot[g]}", va="center",
                fontsize=fs.FS_NOTE - 1.0, color=fs.MUTED)
    ax.set_yticks(range(len(groups)), groups, fontsize=fs.FS_TICK - 1.0)
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(lambda x, _: f"{x:.0%}")
    ax.set_xlabel("members missed by domain search", fontsize=fs.FS_LABEL)
    fs.despine(ax)
    fs.panel(ax, "B", "By lineage (census families)")
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "domain_search_missed"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
