"""S12 figure — structures (`results/structures/figures/structures.png`).

Drawn from the committed S12 tables (D13), through `figstyle.py`:

* **A** — AlphaFold DB coverage of the final census per superfamily: share of
  members with a usable model (exact sequence, global pLDDT ≥ 70), a model
  that is not usable, and no model; the share with any PDB entry marked.
* **B** — the AlphaFold models checked against experiment: TM-score of each
  family's S8a model unit against the experimental unit of the same protein
  (normalised by the experimental unit), by experimental method.
* **C** — superfamily recovery in the dense set: share of representatives
  whose closest structure outside their own family lies in their own
  superfamily, by TM-align and by Foldseek.
* **D** — the five literature edges read four ways: S8a's primary reading,
  experimental units, TM-region units, and the dense set.

    bin/envpy scripts/s12_figures.py
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

RES = ROOT / "results" / "structures"
OUT = RES / "figures"
SHORT = {"tmem16_like": "TMEM16", "iglur": "iGluR", "ploop": "P-loop",
         "innexin_like": "innexin", "connexin": "connexin", "ca_release": "Ca-release",
         "hv": "Hv1", "ploop:VSD": "Kv VSD"}
METHOD = {"Electron Microscopy": (fs.BLUES[5], "o", "cryo-EM"),
          "X-ray diffraction": (fs.SUPERFAMILY["cysloop"], "s", "X-ray"),
          "Solution NMR": (fs.ACCENT, "^", "solution NMR")}
READINGS = [("s8a_afdb", "S8a (AFDB, primary)"), ("experimental", "experimental units"),
            ("tm_region", "TM-region units"), ("dense", "dense set")]


def panel_coverage(ax) -> None:
    rows = read_tsv(RES / "coverage_members.tsv")
    by = defaultdict(lambda: [0, 0, 0, 0])
    for r in rows:
        b = by[r["superfamily"]]
        b[0] += 1
        b[1] += int(r["usable"])
        b[2] += int(r["model"]) - int(r["usable"])
        b[3] += int(r["n_pdb"] != "0")
    order = sorted(by, key=lambda k: by[k][1] / by[k][0])
    for i, k in enumerate(order):
        n, u, m, p = by[k]
        ax.barh(i, u / n, color=fs.BLUES[4], height=0.75, lw=0)
        ax.barh(i, m / n, left=u / n, color=fs.BLUES[1], height=0.75, lw=0)
        ax.barh(i, 1 - (u + m) / n, left=(u + m) / n, color=fs.GRID, height=0.75, lw=0)
        ax.plot([p / n], [i], marker="|", color=fs.INK, ms=5, mew=1.0)
        ax.text(1.01, i, f"{n:,}", fontsize=4.6, va="center", color=fs.MUTED)
    ax.set_yticks(range(len(order)), [k.replace("_", " ") for k in order], fontsize=4.6)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("share of census members", fontsize=6)
    ax.tick_params(axis="x", labelsize=5)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=fs.BLUES[4], label="usable model"),
                       Patch(color=fs.BLUES[1], label="model, not usable"),
                       Patch(color=fs.GRID, label="no model"),
                       Line2D([], [], marker="|", color=fs.INK, ls="", ms=5,
                              label="share with a PDB entry")],
              fontsize=4.8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.12), ncol=2)
    fs.despine(ax)
    fs.panel(ax, "A", "AlphaFold DB coverage of the census")


def panel_check(ax) -> None:
    rows = [r for r in read_tsv(RES / "experimental.tsv") if r["tm_by_exp"]]
    for meth, (col, mk, lab) in METHOD.items():
        xs = [int(r["observed"]) for r in rows if r["method"] == meth]
        ys = [float(r["tm_by_exp"]) for r in rows if r["method"] == meth]
        ax.scatter(xs, ys, s=10, color=col, marker=mk, lw=0, label=f"{lab} ({len(xs)})",
                   zorder=3)
    for r in rows:
        if float(r["tm_by_exp"]) < 0.5:
            ax.annotate(r["node"].replace("_", " "), (int(r["observed"]),
                        float(r["tm_by_exp"])), xytext=(4, -9), textcoords="offset points",
                        fontsize=4.8, color=fs.MUTED)
    for y, lab in ((0.5, "0.5 consistent"), (0.8, "0.8")):
        ax.axhline(y, color=fs.MUTED, lw=0.5, ls=(0, (4, 2)))
        ax.text(1550, y + 0.01, lab, fontsize=4.8, color=fs.MUTED, ha="right")
    ax.set_xscale("log")
    ax.set_xlim(40, 1600)
    ax.set_ylim(0.25, 1.02)
    ax.set_xlabel("experimental unit, observed residues (log scale)", fontsize=6)
    ax.set_ylabel("TM-score, model vs experiment", fontsize=6)
    ax.tick_params(labelsize=5)
    ax.legend(fontsize=4.8, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "B", "AlphaFold models against experiment")


def panel_edges(ax) -> None:
    edges = read_tsv(RES / "network_edges.tsv")
    dense = read_tsv(RES / "dense_edges.tsv")
    keys = [(e["a"], e["b"]) for e in edges if e["variant"] == "s8a_afdb"]
    look = {(e["variant"], e["a"], e["b"]): e for e in edges}
    look.update({("dense", e["a"], e["b"]): e for e in dense})
    cols = [fs.INK, fs.BLUES[5], fs.SUPERFAMILY["iglur"], fs.SUPERFAMILY["cysloop"]]
    for j, ((v, lab), col) in enumerate(zip(READINGS, cols)):
        for i, (a, b) in enumerate(keys):
            e = look.get((v, a, b))
            if not e or not e["median_tm"]:
                ax.text(i + (j - 1.5) * 0.17, 0.33, "×", fontsize=6, color=col,
                        ha="center", va="center")
                continue
            sup = e["verdict"] == "supported"
            ax.scatter(i + (j - 1.5) * 0.17, float(e["median_tm"]), s=16, marker="o",
                       facecolor=col if sup else "white", edgecolor=col, lw=0.9,
                       zorder=3)
    ax.axhline(0.5, color=fs.MUTED, lw=0.5, ls=(0, (4, 2)))
    ax.set_xticks(range(len(keys)),
                  [f"{SHORT.get(a, a)}–\n{SHORT.get(b, b)}" for a, b in keys], fontsize=5)
    ax.set_ylim(0.3, 0.7)
    ax.set_xlim(-0.6, len(keys) - 0.4)
    ax.set_ylabel("median TM-score of the edge", fontsize=6)
    ax.tick_params(labelsize=5)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, ms=4, label=lab)
                       for (_, lab), c in zip(READINGS, cols)],
              fontsize=4.6, frameon=False, loc="upper right", ncol=2)
    ax.text(len(keys) - 0.45, 0.305, "filled = supported (D48); open = not distinguished;"
            " × = unmeasured", fontsize=4.4, color=fs.MUTED, ha="right")
    fs.despine(ax)
    fs.panel(ax, "D", "The literature edges, read four ways")


def panel_recovery(ax) -> None:
    rows = read_tsv(RES / "dense_recovery_by_superfamily.tsv")
    rows.sort(key=lambda r: int(r["tm_recovered"]) / int(r["reps"]))
    for i, r in enumerate(rows):
        n = int(r["reps"])
        ax.barh(i + 0.18, int(r["tm_recovered"]) / n, height=0.34, color=fs.BLUES[4], lw=0)
        ax.barh(i - 0.18, int(r["foldseek_recovered"]) / n, height=0.34,
                color=fs.SUPERFAMILY["cysloop"], lw=0)
        ax.text(1.01, i, f"{n} reps, {r['families']} fam.", fontsize=4.6, va="center",
                color=fs.MUTED)
    ax.set_yticks(range(len(rows)), [r["superfamily"].replace("_", " ") for r in rows],
                  fontsize=5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of representatives recovered", fontsize=6)
    ax.tick_params(axis="x", labelsize=5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=fs.BLUES[4], label="best TM-align partner"),
                       Patch(color=fs.SUPERFAMILY["cysloop"],
                             label="best Foldseek E-value")],
              fontsize=4.8, frameon=False, loc="upper center",
              bbox_to_anchor=(0.45, -0.12), ncol=2)
    fs.despine(ax)
    fs.panel(ax, "C", "Does structure find the superfamily?")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(fs.W_FULL, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1.25, 1],
                          wspace=0.55, hspace=0.5)
    panel_coverage(fig.add_subplot(gs[0, 0]))
    panel_check(fig.add_subplot(gs[0, 1]))
    panel_recovery(fig.add_subplot(gs[1, 0]))
    panel_edges(fig.add_subplot(gs[1, 1]))
    print("\n".join(map(str, fs.save(fig, OUT / "structures"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
