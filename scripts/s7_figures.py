"""S7 figure — the tier-1 forest: how it was trimmed, how it is rooted, how
well it is supported, and one tree drawn.

Four panels from the committed S7 tables (D13), through `figstyle.py`:

* **A** — trimming (D41): informative sites kept by trimAl `-automated1`
  against `-gt 0.5`, one dot per family.
* **B** — support: share of a family's internal edges at UFBoot ≥ 95
  against its size, coloured by superfamily.
* **C** — rooting: the root edge's UFBoot per rooted family, with the
  outgroup's column occupancy; families whose outgroup is not one clade are
  marked, unrooted families counted by reason. After S7c, each re-rooted
  family shows its S7b root and its two declared outgroups' roots, with the
  D47 verdict.
* **D** — one family tree drawn: the ryanodine receptors rooted on the
  ITPR exemplars, leaves coloured by panel group.

    python3 scripts/s7_figures.py
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
from s7_newick import parse  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

D = ROOT / "results" / "phylogeny"
EXAMPLE = "ryr"
GROUP_COLOURS = {"vertebrate": "#2a78d6", "invertebrate": "#eb6834",
                 "outgroup": fs.INK}


def sf_colour(fam: str) -> str:
    return fs.SUPERFAMILY.get(CATALOGUE[fam].superfamily, fs.SUPERFAMILY["other"])


def panel_a(ax) -> None:
    rows = read_tsv(D / "tier1_trim_compare.tsv")
    a = {r["family"]: int(r["informative"]) for r in rows if r["method"] == "automated1"}
    g = {r["family"]: int(r["informative"]) for r in rows if r["method"] == "gt0.5"}
    for f in a:
        ax.scatter(a[f], g[f], s=10, color=sf_colour(f), alpha=0.85, linewidth=0)
    for f in ("k2p", "nachr", "cng"):
        ax.annotate(f, (a[f], g[f]), fontsize=fs.FS_NOTE, color=fs.MUTED,
                    xytext=(3, -8), textcoords="offset points")
    lo, hi = 5, max(g.values()) * 1.3
    ax.plot([lo, hi], [lo, hi], color=fs.FAINT, lw=0.6, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("informative sites, trimAl -automated1 (S6)", fontsize=fs.FS_LABEL)
    ax.set_ylabel("informative sites, trimAl -gt 0.5 (D41)", fontsize=fs.FS_LABEL)
    ge = sum(g[f] >= a[f] for f in a)
    ax.text(0.03, 0.95, f"≥ automated1 in {ge}/{len(a)} families\nchosen before any tree",
            transform=ax.transAxes, fontsize=fs.FS_NOTE, color=fs.MUTED, va="top")
    ax.legend(handles=fs.superfamily_handles(["ploop", "cysloop", "iglur", "other"]),
              fontsize=fs.FS_NOTE, frameon=False, loc="lower right")
    fs.despine(ax)
    fs.panel(ax, "A", "Trimming, measured (D41)")


def panel_b(ax, trees) -> None:
    for t in trees:
        if t["frac_ge95"] == "":
            continue
        ax.scatter(int(t["n_ingroup"]), float(t["frac_ge95"]), s=6 + int(t["cols"]) / 60,
                   color=sf_colour(t["family"]), alpha=0.85, linewidth=0)
        if int(t["n_ingroup"]) >= 250 or (float(t["frac_ge95"]) < 0.3 and int(t["n_ingroup"]) >= 10):
            ax.annotate(t["family"], (int(t["n_ingroup"]), float(t["frac_ge95"])),
                        fontsize=fs.FS_NOTE, color=fs.MUTED, xytext=(3, 2),
                        textcoords="offset points")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("sequences in the family tree", fontsize=fs.FS_LABEL)
    ax.set_ylabel("internal edges with UFBoot ≥ 95", fontsize=fs.FS_LABEL)
    ax.text(0.97, 0.03, "dot area ~ alignment columns", transform=ax.transAxes,
            ha="right", fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.despine(ax)
    fs.hgrid(ax, "y")
    fs.panel(ax, "B", f"Support across {len(trees)} family trees")


def _reroot() -> tuple[dict, dict]:
    """S7c: {family: [tree rows]} and {family: verdict row}, empty before S7c parse."""
    tp, fp = D / "tier1_reroot_trees.tsv", D / "tier1_reroot.tsv"
    if not fp.exists():
        return {}, {}
    by = {}
    for t in read_tsv(tp):
        by.setdefault(t["family"], []).append(t)
    return by, {r["family"]: r for r in read_tsv(fp)}


def _bar(ax, y, t, colour, h) -> None:
    clade = str(t["outgroup_monophyletic"]) == "True"
    v = float(t["root_ufboot"]) if clade and t["root_ufboot"] else 0
    ax.barh(y, v if clade else 100, color=colour if clade else "none",
            edgecolor=fs.FAINT if not clade else "none", hatch=None if clade else "////",
            height=h, linewidth=0.5)


def panel_c(ax, trees, inputs) -> None:
    occ = {r["family"]: max(float(x) for x in r["og_occupancy"].split(","))
           for r in inputs if r["og_occupancy"]}
    rr_trees, rr = _reroot()
    rooted = sorted((t for t in trees if t["outgroup_family"]),
                    key=lambda t: (CATALOGUE[t["family"]].superfamily, t["family"]))
    for i, t in enumerate(rooted):
        clade = t["outgroup_monophyletic"] == "True"
        single = t["n_outgroup"] == "1"
        if single:      # one outgroup sequence: a terminal edge, no UFBoot exists
            ax.text(2, i, "single-sequence outgroup: root edge has no support value",
                    va="center", fontsize=fs.FS_NOTE - 1.2, color=fs.MUTED)
            continue
        if t["family"] in rr:          # S7c: old root (grey) + the two declared outgroups
            _bar(ax, i - 0.27, t, fs.FAINT, 0.24)
            for k, nt in enumerate(rr_trees[t["family"]]):
                _bar(ax, i + 0.27 * k, nt, sf_colour(t["family"]), 0.24)
            ok = rr[t["family"]]["verdict"] == "resolved"
            ax.text(102, i, "resolved" if ok else "unresolved", va="center",
                    fontsize=fs.FS_NOTE - 0.8, color=fs.INK if ok else fs.ACCENT)
            continue
        _bar(ax, i, t, sf_colour(t["family"]), 0.72)
        ax.text(102, i, f"{occ.get(t['family'], 0):.2f}", va="center",
                fontsize=fs.FS_NOTE - 0.6, color=fs.MUTED)
    ax.axvline(95, color=fs.ACCENT, lw=0.7, ls="--")
    ax.set_yticks(range(len(rooted)), [t["family"] for t in rooted],
                  fontsize=fs.FS_TICK - 1.2)
    ax.set_xlim(0, 125 if rr else 115)
    ax.set_xticks([0, 50, 95])
    ax.set_xlabel("root-edge UFBoot (hatched: outgroup not one clade)"
                  + ("\nS7c rows: grey = KcsA root (S7b), then outgroups 1, 2" if rr else ""),
                  fontsize=fs.FS_LABEL)
    ax.text(102, len(rooted) - 0.2, "occ.", fontsize=fs.FS_NOTE - 0.6, color=fs.MUTED)
    c = Counter(t["root_rule"].split(":")[-1] for t in trees if not t["outgroup_family"])
    ax.text(0.02, -0.075 if len(rooted) > 30 else -0.12,
            "unrooted: " + ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in c.items()),
            transform=ax.transAxes, fontsize=fs.FS_NOTE, color=fs.MUTED, va="top")
    fs.despine(ax)
    fs.panel(ax, "C", "Roots from the catalogue's outgroups"
             + (" (S7c re-rooted)" if rr else ""))


def _layout(node, depth=0.0, ys=None, out=None):
    """Rectangular phylogram coordinates: x = root distance, y = leaf order."""
    if out is None:
        out, ys = [], [0]
    x = depth + (node.length or 0.0)
    if not node.children:
        y = ys[0]
        ys[0] += 1
        out.append(("leaf", node.name, x, y))
        return x, y, out
    kids = [_layout(c, x, ys, out)[:2] for c in node.children]
    y = (kids[0][1] + kids[-1][1]) / 2
    out.append(("node", node.name, x, y, kids))
    return x, y, out


def panel_d(ax) -> None:
    tf = D / "tier1" / f"{EXAMPLE}.treefile"
    if not tf.exists():
        ax.text(0.5, 0.5, "tree not yet built", ha="center", transform=ax.transAxes)
        fs.panel(ax, "D", EXAMPLE)
        return
    group = {m["label"]: m["group"] for m in read_tsv(ROOT / "results" / "alignments" /
                                                      "members.tsv")}
    _, _, items = _layout(parse(tf.read_text()))
    cnt = Counter()
    for it in items:
        if it[0] == "node":
            _, name, x, y, kids = it
            ax.plot([x, x], [kids[0][1], kids[-1][1]], color=fs.MUTED, lw=0.5)
            for kx, ky in kids:
                ax.plot([x, kx], [ky, ky], color=fs.MUTED, lw=0.5)
            try:
                if name and float(name) >= 95 and len(kids) > 1:
                    ax.scatter(x, y, s=4, color=fs.INK, zorder=3, linewidth=0)
            except ValueError:
                pass
        else:
            _, name, x, y = it
            g = "outgroup" if name.startswith("OG_") else group.get(name, "other")
            cnt[g] += 1
            ax.scatter(x, y, s=5, color=GROUP_COLOURS.get(g, fs.FAINT), linewidth=0, zorder=3)
            if g == "outgroup":
                ax.text(x, y, "  " + name[3:], va="center", fontsize=fs.FS_NOTE - 0.6)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("substitutions / site", fontsize=fs.FS_LABEL)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", ls="", color=GROUP_COLOURS.get(g, fs.FAINT),
                      markersize=3, label=f"{g} ({n})") for g, n in cnt.most_common()]
    handles.append(Line2D([], [], marker="o", ls="", color=fs.INK, markersize=2,
                          label="node UFBoot ≥ 95"))
    ax.legend(handles=handles, fontsize=fs.FS_NOTE, frameon=False, loc="upper left")
    fs.despine(ax, keep=("bottom",))
    fs.panel(ax, "D", "Ryanodine receptors, rooted on ITPR")


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    tp = D / "tier1_trees.tsv"
    trees = read_tsv(tp) if tp.exists() else []
    inputs = read_tsv(D / "tier1_inputs.tsv")
    if not trees:                       # S7a: the trimming decision alone
        fig, ax = plt.subplots(figsize=(fs.W_HALF, 3.0))
        panel_a(ax)
        fig.tight_layout()
        print("\n".join(map(str, fs.save(fig, D / "figures" / "tier1_trim"))))
        return 0
    fig, axes = plt.subplots(2, 2, figsize=(fs.W_FULL, 8.2),
                             gridspec_kw={"height_ratios": [1, 1.7]})
    panel_a(axes[0][0])
    panel_b(axes[0][1], trees)
    panel_c(axes[1][0], trees, inputs)
    panel_d(axes[1][1])
    fig.tight_layout()
    paths = fs.save(fig, D / "figures" / "tier1_trees")
    print("\n".join(map(str, paths)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
