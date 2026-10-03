"""s10_figures.py — S10a's headline figure, drawn from the committed tables (D13).

    bin/envpy scripts/s10_figures.py

`results/repertoire/figures/repertoire.png` (+ .pdf):

A. The reconstructed state of every census family at fifteen ancestral nodes
   of the NCBI species tree (primary cost, gain = 2 losses); a dot marks a
   node whose state changes under the other two costs.
B. Stated gains and losses per family on the 437-order tree (primary cost),
   the robust part (stated under all three costs) solid.
C. Every stated loss by absence strength (D46/D50): proteome-only,
   controlled by S5's genome sweep, or contradicted by it.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict

import figstyle as fs
from s3_hmm_lib import read_tsv
from s10_lib import GAINS, OUT, PRIMARY

sys.path.insert(0, str(OUT.parents[1]))
from src.catalogue import registry  # noqa: E402

NODES = [("Eukaryota", "LECA"), ("Opisthokonta", "Opisthokonta"), ("Fungi", "Fungi"),
         ("Metazoa", "Metazoa"), ("Bilateria", "Bilateria"), ("Ecdysozoa", "Ecdysozoa"),
         ("Lophotrochozoa", "Lophotrochozoa"), ("Deuterostomia", "Deuterostomia"),
         ("Vertebrata", "Vertebrata"), ("Mammalia", "Mammalia"),
         ("Viridiplantae", "Viridiplantae"), ("Embryophyta", "Land plants"),
         ("Sar", "SAR"), ("Amoebozoa", "Amoebozoa"), ("Discoba", "Discoba")]
SHORT = {
    "kv_shaker": "Kv (Shaker)", "kv_kcnq": "Kv (KCNQ)", "kv_eag": "EAG", "kv_modifier": "Kv modifier",
    "kca_slo": "Slo", "kca_sk": "SK/IK", "kir": "Kir", "k2p": "K2P", "hcn": "HCN", "cng": "CNG",
    "nav": "Nav", "cav": "Cav", "nalcn": "NALCN", "catsper": "CatSper", "tpc": "TPC",
    "trpa": "TRPA", "trpc": "TRPC", "trpm": "TRPM", "trpml": "TRPML", "trpn": "TRPN",
    "trpp": "TRPP", "trpv": "TRPV", "kcsa_prok": "KcsA-like (prok.)", "nachr": "nAChR",
    "gabaa": "GABA-A", "glyr": "GlyR", "ht3": "5-HT3", "zac": "ZAC",
    "plgic_invertebrate": "pLGIC (invert.)", "plgic_prok": "pLGIC (prok.)", "ampa": "AMPA",
    "kainate": "Kainate", "nmda": "NMDA", "delta_glur": "Delta GluR",
    "iglur_nonvertebrate": "iGluR (non-vert.)", "iglur_prok": "GluR0 (prok.)",
    "p2x": "P2X (metazoan)", "p2x_nonmetazoan": "P2X (non-metazoan)", "asic": "ASIC",
    "enac": "ENaC", "deg_invertebrate": "Degenerin (invert.)", "clc_channel": "CLC",
    "bestrophin": "Bestrophin", "cftr": "CFTR", "tweety": "Tweety",
    "ano_channel": "ANO channel", "ano_scramblase": "ANO scramblase", "osca_tmem63": "OSCA/TMEM63",
    "tmc": "TMC", "piezo": "Piezo", "mscl": "MscL", "mscs": "MscS", "connexin": "Connexin",
    "innexin": "Innexin", "pannexin": "Pannexin", "lrrc8": "LRRC8", "calhm": "CALHM",
    "itpr": "IP3R", "ryr": "RyR", "tric": "TRIC", "mcu": "MCU", "vdac": "VDAC",
    "tmem175": "TMEM175", "orai": "Orai", "hv1": "Hv1", "otop": "Otopetrin", "clic": "CLIC",
    "viroporin": "Viroporin", "pacc": "PACC1", "tmco1": "TMCO1", "tmem87": "TMEM87",
    "tmem109": "TMEM109", "clcc1": "CLCC1", "mitok": "MITOK", "gphr": "GPHR"}
PRESENT, ABSENT, AMBIG = fs.SUPERFAMILY["ploop"], "#eef2f8", fs.GRID
STRENGTH = [("proteome-only", fs.FAINT), ("controlled", fs.SUPERFAMILY["iglur"]),
            ("contradicted", fs.SUPERFAMILY["cysloop"])]


def family_order(fams: set[str]) -> list[str]:
    """Superfamily blocks in catalogue order, families in catalogue order."""
    return [f.key for f in registry.census_families() if f.key in fams]


def main() -> int:
    fs.use()
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    nodes = defaultdict(dict)
    for r in read_tsv(OUT / "nodes.tsv"):
        nodes[(r["family"], r["node"])][r["cost"]] = r["state"]
    fam = {(r["family"], r["cost"]): r for r in read_tsv(OUT / "families.tsv")}
    ev = read_tsv(OUT / "events.tsv")
    losses = read_tsv(OUT / "losses.tsv")
    fams = family_order({f for f, _ in fam})
    sf_of = {f.key: f.superfamily for f in registry.census_families()}

    fig = plt.figure(figsize=(fs.W_FULL, 8.4))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.55, 1.0], height_ratios=[1.0, 0.06],
                  wspace=0.06, hspace=0.20, left=0.13, right=0.98, top=0.935, bottom=0.03)
    ax = fig.add_subplot(gs[0, 0])
    n = len(fams)
    for i, f in enumerate(fams):
        for j, (node, _) in enumerate(NODES):
            st = nodes[(f, node)]
            s = st.get(PRIMARY, "?")
            c = PRESENT if s == "1" else ABSENT if s == "0" else AMBIG
            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=c, edgecolor="white", lw=0.4))
            if len(set(st.values())) > 1:
                ax.scatter(j + 0.5, i + 0.5, s=2.2, color=fs.INK if s != "1" else "white", lw=0)
    prev = None
    for i, f in enumerate(fams):
        if prev is not None and sf_of[f] != prev:
            ax.axhline(i, color=fs.INK, lw=0.5)
        prev = sf_of[f]
    ax.set_xlim(0, len(NODES))
    ax.set_ylim(n, 0)
    ax.set_yticks([i + 0.5 for i in range(n)])
    ax.set_yticklabels([SHORT.get(f, f) for f in fams], fontsize=4.6)
    ax.set_xticks([j + 0.5 for j in range(len(NODES))])
    ax.set_xticklabels([lab for _, lab in NODES], rotation=55, ha="right",
                       rotation_mode="anchor", fontsize=5.4)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    fs.panel(ax, "A", "Reconstructed state at ancestral nodes")

    bx = fig.add_subplot(gs[0, 1], sharey=ax)
    rob = Counter((e["family"], e["event"]) for e in ev
                  if e["cost"] == PRIMARY and e["stated"] == "1" and e["robust"] == "1")
    for i, f in enumerate(fams):
        r = fam[(f, PRIMARY)]
        g, lo = int(r["gains"]), int(r["losses"])
        bx.barh(i + 0.5, -g, height=0.7, color=PRESENT, alpha=0.35, lw=0)
        bx.barh(i + 0.5, -rob[(f, "gain")], height=0.7, color=PRESENT, lw=0)
        bx.barh(i + 0.5, lo, height=0.7, color=fs.SUPERFAMILY["cysloop"], alpha=0.35, lw=0)
        bx.barh(i + 0.5, rob[(f, "loss")], height=0.7, color=fs.SUPERFAMILY["cysloop"], lw=0)
    bx.axvline(0, color=fs.INK, lw=0.5)
    bx.set_xlim(-14, 60)
    bx.set_xticks([-10, 0, 20, 40, 60])
    bx.set_xticklabels(["10", "0", "20", "40", "60"])
    bx.tick_params(axis="y", left=False, labelleft=False)
    bx.text(-7, 1.0, "gains", ha="center", fontsize=fs.FS_NOTE, color=PRESENT)
    bx.text(25, 1.0, "losses", ha="center", fontsize=fs.FS_NOTE, color=fs.SUPERFAMILY["cysloop"])
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    keys = [Patch(facecolor=PRESENT, label="A: present"), Patch(facecolor=ABSENT, edgecolor=fs.GRID,
            lw=0.4, label="A: absent"), Patch(facecolor=AMBIG, label="A: ambiguous (tied optima)"),
            Line2D([], [], marker="o", ls="", ms=2, color=fs.INK, label="A: state differs under g = 1 or Dollo"),
            Patch(facecolor=PRESENT, alpha=0.35, label="B: gains, cost-dependent"),
            Patch(facecolor=fs.SUPERFAMILY["cysloop"], alpha=0.35, label="B: losses, cost-dependent")]
    bx.legend(handles=keys, loc="upper right", bbox_to_anchor=(1.0, 0.985), fontsize=4.8,
              frameon=True, framealpha=1, edgecolor=fs.GRID, handlelength=1.2)
    fs.despine(bx, keep=("bottom",))
    bx.grid(axis="x", color=fs.GRID, lw=0.4)
    bx.set_axisbelow(True)
    bx.set_xlabel("stated events on the order tree (solid: under all three costs)",
                  fontsize=fs.FS_TICK)
    fs.panel(bx, "B", "Gains and losses per family")

    cx = fig.add_subplot(gs[1, :])
    c = Counter(r["strength"] for r in losses)
    left, total = 0, sum(c.values())
    for k, col in STRENGTH:
        cx.barh(0, c[k], left=left, color=col, height=0.8, lw=0)
        left += c[k]
    x = 0
    for k, col in STRENGTH:
        cx.text(x, -0.75, f"{k} {c[k]}", va="top", ha="left", fontsize=fs.FS_NOTE, color=col)
        x += total * 0.22
    cx.set_xlim(0, total)
    cx.set_ylim(-1.6, 0.5)
    cx.set_yticks([])
    cx.set_xticks([])
    for s in cx.spines.values():
        s.set_visible(False)
    fs.panel(cx, "C", f"All {len(losses)} stated losses by absence strength (S5 genome overlay)")

    fig.suptitle("Ion-channel repertoire across 437 eukaryotic orders (parsimony on NCBI taxonomy)",
                 fontsize=fs.FS_SUPTITLE, y=0.99)
    out = fs.save(fig, OUT / "figures" / "repertoire")
    print("\n".join(str(p) for p in out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
