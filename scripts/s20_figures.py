"""S20 figure — what the published channelome is made of, and the auxiliaries.

Four panels from the committed S20 tables (D13), through `figstyle.py`:

* **A** — each database channelome (GtoPdb, HGNC, UniProt KW-0407, their
  union) split into pore-forming census genes and everything else, by reason.
* **B** — census pore genes each list leaves out, by family.
* **C** — the auxiliary families against their homology groups: grouped by
  the channel they serve, not by descent.
* **D** — the panel census of auxiliaries per homology group: S3b profile
  calls that pass the reciprocal-best-human-hit test, and those that fail.

    python3 scripts/s20_figures.py
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

D = ROOT / "results" / "auxiliary"
LISTS = [("gtopdb", "GtoPdb"), ("hgnc", "HGNC"), ("uniprot", "UniProt\nKW-0407"),
         ("any_list", "union")]
PARTS = [("pore_census", "pore-forming (census)", "#184f95"),
         ("auxiliary", "auxiliary subunit", "#eb6834"),
         ("auxiliary_uncatalogued", "auxiliary, uncatalogued", "#f3a683"),
         ("out_of_scope", "aquaporin", "#4a3aa7"),
         ("out_of_scope_uncatalogued", None, "#4a3aa7"),
         ("transporter", "CLC / SLC26 transporter", "#8a897f"),
         ("transporter_or_enzyme", "transporter or enzyme", "#a9a79e"),
         ("pore_candidate", "proposed pore, uncatalogued", "#1baf7a"),
         ("paracellular", "claudin", "#d6d5cf"),
         ("pseudogene", "pseudogene", "#e5e4df")]


def panel_a(ax) -> None:
    dec = {r["list"]: r for r in read_tsv(D / "list_decomposition.tsv")}
    for i, (k, lab) in enumerate(LISTS):
        left = 0
        for key, name, col in PARTS:
            v = int(dec[k][key])
            ax.barh(i, v, left=left, color=col, height=0.66, linewidth=0,
                    label=name if (i == 0 and name) else None)
            left += v
        ax.text(left + 4, i, str(left), va="center", fontsize=fs.FS_NOTE, color=fs.INK)
    ax.set_yticks(range(len(LISTS)), [l for _, l in LISTS], fontsize=fs.FS_TICK)
    ax.invert_yaxis()
    ax.set_xlim(0, 460)
    ax.set_xlabel("human genes", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE - 0.6, frameon=False, loc="upper center", ncol=2,
              bbox_to_anchor=(0.45, -0.16), handlelength=1.0, columnspacing=0.8)
    fs.despine(ax)
    fs.panel(ax, "A", "Three database channelomes")


def panel_b(ax) -> None:
    rows = read_tsv(D / "census_missed_by_list.tsv")[:14]
    keys = [("gtopdb", "GtoPdb", "#184f95"), ("hgnc", "HGNC", "#6da7ec"),
            ("uniprot", "UniProt", "#eb6834")]
    h = 0.26
    for j, (k, lab, col) in enumerate(keys):
        ax.barh([i + (j - 1) * h for i in range(len(rows))],
                [int(r[f"{k}_missed"]) for r in rows], height=h, color=col,
                linewidth=0, label=lab)
    ax.set_yticks(range(len(rows)), [r["family"] for r in rows], fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlabel("census pore genes not on the list", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.hgrid(ax, "x")
    fs.panel(ax, "B", "What each list leaves out")


def panel_c(ax) -> None:
    groups = read_tsv(D / "aux_groups.tsv")
    fams = sorted({r["family"] for r in groups},
                  key=lambda f: (-sum(r["family"] == f for r in groups), f))
    for i, f in enumerate(fams):
        gs = [r for r in groups if r["family"] == f]
        left = 0
        for j, g in enumerate(gs):
            n = len(g["genes"].split(","))
            ax.barh(i, n, left=left, height=0.66, linewidth=0.6, edgecolor="white",
                    color=[fs.BLUES[5], fs.BLUES[3], fs.ACCENT, "#eb6834"][j % 4] if len(gs) > 1 else fs.BLUES[1])
            if n >= 2:
                ax.text(left + n / 2, i, g["genes"].split(",")[0].rstrip("0123456789")
                        [:7], ha="center", va="center", fontsize=fs.FS_NOTE - 1.4,
                        color="white" if len(gs) > 1 else fs.INK)
            left += n
    ax.set_yticks(range(len(fams)), [f.replace("assoc_", "") for f in fams],
                  fontsize=fs.FS_TICK - 0.6)
    ax.invert_yaxis()
    ax.set_xlabel("human genes (segments = homology groups)", fontsize=fs.FS_LABEL)
    pooled = sum(1 for f in fams if sum(r["family"] == f for r in groups) > 1)
    ax.text(0.97, 0.03, f"{pooled}/{len(fams)} families pool\nunrelated proteins",
            transform=ax.transAxes, ha="right", fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.despine(ax)
    fs.panel(ax, "C", "Auxiliary families are grouped by partner")


def panel_d(ax) -> None:
    p = D / "aux_panel.tsv"
    if not p.exists():
        ax.text(0.5, 0.5, "panel search not yet run", ha="center", transform=ax.transAxes)
        fs.panel(ax, "D", "Auxiliaries across the panel")
        return
    rows = read_tsv(p)
    groups = read_tsv(D / "aux_groups.tsv")
    fail = Counter(r["family"] for r in rows if not r["homology_group"].startswith("assoc_"))
    labels, ok, bad = [], [], []
    for g in groups:
        labels.append(f"{g['group'].replace('assoc_', '')} {g['genes'].split(',')[0]}")
        ok.append(sum(r["homology_group"] == g["group"] for r in rows))
        bad.append(fail[g["family"]] if g["group"].endswith(".1") else 0)
    y = range(len(labels))
    ax.barh(list(y), ok, color="#184f95", height=0.66, linewidth=0,
            label="best human hit is the group")
    ax.barh(list(y), bad, left=ok, color="#ef9a90", height=0.66, linewidth=0,
            label="best human hit elsewhere (family total)")
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yticks(list(y), labels, fontsize=fs.FS_TICK - 1.2)
    ax.invert_yaxis()
    ax.set_xlabel("S3b panel profile calls", fontsize=fs.FS_LABEL)
    ax.legend(fontsize=fs.FS_NOTE - 0.4, frameon=False, loc="upper center",
              bbox_to_anchor=(0.4, -0.1), ncol=1)
    fs.despine(ax)
    fs.panel(ax, "D", "Auxiliaries across the panel")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(fs.W_FULL, 7.8),
                             gridspec_kw={"height_ratios": [1, 1.25]})
    panel_a(axes[0][0])
    panel_b(axes[0][1])
    panel_c(axes[1][0])
    panel_d(axes[1][1])
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "auxiliary"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
