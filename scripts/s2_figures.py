"""S2 figure — census v2: every record enumerated, and how far domain rules get.

* **A** — records per superfamily, split by S2's call: a family, the
  superfamily only, or unassigned. Families with identical architectures
  (Cys-loop, iGluR, DEG/ENaC, P2X, CLC) stop at the superfamily by design.
* **B** — the largest families by record count, called by domain rules.

Drawn from `results/census_v2/census_families.tsv`, `census_superfamily_only.tsv`
and `census_status.tsv`.

    python3 scripts/s2_figures.py
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

D = ROOT / "results" / "census_v2"


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fam = read_tsv(D / "census_families.tsv")
    sfo = read_tsv(D / "census_superfamily_only.tsv")
    status = {r["status"]: 0 for r in read_tsv(D / "census_status.tsv")}
    for r in read_tsv(D / "census_status.tsv"):
        status[r["status"]] += int(r["records"])
    sf_col = [c for c in sfo[0] if c != "records"][0]
    by_fam, by_sfo = Counter(), Counter()
    for r in fam:
        by_fam[r["superfamily"] or "other"] += int(r["records"])
    for r in sfo:
        by_sfo[r[sf_col]] += int(r["records"])
    sfs = sorted(set(by_fam) | set(by_sfo), key=lambda s: -(by_fam[s] + by_sfo[s]))[:14]
    fig, axes = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.8))
    ax = axes[0]
    for i, s in enumerate(sfs):
        ax.barh(i, by_fam[s], color="#184f95", height=0.7, linewidth=0,
                label="family call" if i == 0 else None)
        ax.barh(i, by_sfo[s], left=by_fam[s], color="#86b6ef", height=0.7, linewidth=0,
                label="superfamily only" if i == 0 else None)
    ax.set_yticks(range(len(sfs)), [fs.SUPERFAMILY_LABEL.get(s, s) for s in sfs],
                  fontsize=fs.FS_TICK - 0.8)
    ax.invert_yaxis()
    ax.set_xlabel("census v2 records (thousands)", fontsize=fs.FS_LABEL)
    ax.xaxis.set_major_formatter(lambda x, _: f"{x / 1000:.0f}")
    n = sum(status.values())
    ax.text(0.97, 0.03, f"{n:,} records; unassigned {status.get('unassigned', 0):,}",
            transform=ax.transAxes, ha="right", fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="center right")
    fs.despine(ax)
    fs.panel(ax, "A", "Called by domain rules, per superfamily")
    ax = axes[1]
    top = sorted((r for r in fam if r["catalogue_status"] in ("channel", "channel_contested")),
                 key=lambda r: -int(r["records"]))[:16]
    ax.barh(range(len(top)), [int(r["records"]) for r in top],
            color=[fs.SUPERFAMILY.get(r["superfamily"], fs.SUPERFAMILY["other"]) for r in top],
            height=0.7, linewidth=0)
    ax.set_yticks(range(len(top)), [r["family"] for r in top], fontsize=fs.FS_TICK - 0.8)
    ax.invert_yaxis()
    ax.set_xlabel("records called to the census family (thousands)", fontsize=fs.FS_LABEL)
    ax.xaxis.set_major_formatter(lambda x, _: f"{x / 1000:.0f}")
    ax.legend(handles=fs.superfamily_handles(["ploop", "cysloop", "iglur", "other"]),
              fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "B", "The largest families")
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "census_v2"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
