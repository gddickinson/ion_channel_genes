"""Review figures 4, 5 and 8 — length, domain coverage, and the architecture traps.

All three are rendered from committed tables and none is a schematic.

**Figure 4 — the length range.** Every reference protein, by superfamily, on
a log axis. MscL is 136 residues and RYR1 is 5,038; the 37-fold range inside
a single superfamily is the reason a tier-2 tree has to be built on the pore
module rather than on full-length sequence (**D27**).

**Figure 5 — which families carry which domain.** A presence matrix from
`exemplar_architecture.tsv`, ordered so the shared signatures come first.
This is section 12 in one panel: `PF00520` reaching twenty families
including a phosphatase, and the TRP families that carry no copy of it.

**Figure 8 — the architecture traps.** Scale diagrams from measured InterPro
coordinates for the four pairs domain composition cannot separate. The
proteins in each pair are drawn to the same scale, so "identical
architecture" is something the reader can check rather than take on trust.

    python3 scripts/s0_review_fig_domains.py
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs
import review_figlib as R
from scripts.s0_lib import read_tsv
from src.catalogue import CATALOGUE, SUPERFAMILIES

DATA = ROOT / "results" / "s0_baseline"
OUT = ROOT / "docs" / "figures"

#: The pairs and trios section 12.4 tabulates, with the point of each.
TRAPS: list[tuple[str, list[str], str]] = [
    ("Four repeats, three different ions",
     ["SCN5A", "CACNA1C", "CACNA1G", "NALCN"],
     "identical composition; separated only by the filter (fig. 3)"),
    ("An ABC transporter that is a channel, and one that is not",
     ["CFTR", "ABCC8"],
     "one accession apart: PF14396, CFTR's R domain"),
    ("A chloride channel and a lipid scramblase",
     ["ANO1", "ANO6"],
     "no sequence-level difference is known"),
    ("One domain, three families",
     ["MCOLN1", "PKD2", "PKD1"],
     "PF08016 spans TRPML, TRPP and a protein that is not a pore"),
    ("Domains that exist outside channels",
     ["KCNA1", "KCTD1", "GRIA1", "GRM1", "CHRNA1", "AChBP", "HVCN1", "TPTE"],
     "each pair shares its first-listed domain and only one is a channel"),
]


def short_names() -> dict[str, str]:
    """`{accession: Pfam short name}` from the S0 verification table.

    The protein-level InterPro endpoint returns the *long* name
    ("Ion transport protein"), which does not fit in a domain box and turns
    a scale diagram into overlapping text. The short names ("Ion_trans")
    were already resolved and committed by `s0_catalogue_verify.py`.
    """
    return {r["accession"]: (r["interpro_short"] or r["accession"])
            for r in read_tsv(DATA / "pfam_verification.tsv")}


def domain_colours(rows: list[dict]) -> dict[str, str]:
    """Stable colours: the shared/decoy domains get the accent, the rest a ramp."""
    accs, seen = [], set()
    for r in rows:
        if r["pfam"] not in seen:
            seen.add(r["pfam"])
            accs.append(r["pfam"])
    ramp = [fs.BLUES[1], fs.BLUES[3], fs.BLUES[5], fs.SUPERFAMILY["cysloop"],
            fs.SUPERFAMILY["iglur"], fs.ACCENT, fs.MUTED, fs.FAINT]
    return {a: ramp[i % len(ramp)] for i, a in enumerate(accs)}


def fig_lengths() -> Path:
    import matplotlib.pyplot as plt
    seqs: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for r in read_tsv(DATA / "exemplars_resolved.tsv"):
        if not r.get("length_aa", "").isdigit():
            continue
        fam = CATALOGUE.get(r["family"])
        if fam is None or not fam.census_member():
            continue
        seqs[fam.superfamily].append((r["label"], int(r["length_aa"])))

    order = sorted(seqs, key=lambda k: -max(v for _, v in seqs[k]))
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 3.5))
    for i, sf in enumerate(order):
        vals = [v for _, v in seqs[sf]]
        colour = fs.SUPERFAMILY.get(sf, fs.SUPERFAMILY["other"])
        ax.scatter(vals, [i] * len(vals), s=13, color=colour,
                   edgecolor="white", linewidth=0.35, zorder=3)
        ax.plot([min(vals), max(vals)], [i, i], color=colour, linewidth=1.0,
                alpha=0.45, zorder=2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([R.superfamily_label(k) for k in order],
                       fontsize=fs.FS_TICK - 1.0)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("length of the reference protein (residues, log scale)",
                  fontsize=fs.FS_LABEL)
    ax.set_xlim(90, 7000)
    for x, label in ((136, "MscL 136"), (5038, "RYR1 5,038")):
        ax.axvline(x, color=fs.FAINT, linewidth=0.6, linestyle=":", zorder=1)
        ax.text(x, len(order) - 0.3, label, rotation=90, ha="right",
                va="bottom", fontsize=fs.FS_TICK - 1.4, color=fs.MUTED)
    fs.despine(ax)
    fs.hgrid(ax, axis="x")
    fig.tight_layout()
    return fs.save(fig, OUT / "fig4_length_range", legend_in="docs/review (numbered review figure)")[0]


def fig_domain_matrix() -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    rows = read_tsv(DATA / "exemplar_architecture.tsv")
    fam_domains: dict[str, set] = defaultdict(set)
    for r in rows:
        for token in (r["observed_pfam"] or "").split(","):
            if token:
                fam_domains[r["family"]].add(token.split("x")[0])
    counts: dict[str, int] = defaultdict(int)
    for fam, ds in fam_domains.items():
        for d in ds:
            counts[d] += 1
    shared = [d for d, n in sorted(counts.items(), key=lambda kv: -kv[1])
              if n >= 2][:14]
    fams = [f for f in CATALOGUE
            if f in fam_domains and (fam_domains[f] & set(shared))]
    fams.sort(key=lambda f: (CATALOGUE[f].superfamily,
                             not CATALOGUE[f].census_member(), f))

    fig, ax = plt.subplots(figsize=(fs.W_FULL, 0.155 * len(fams) + 1.5))
    for j, d in enumerate(shared):
        for i, f in enumerate(fams):
            if d in fam_domains[f]:
                census = CATALOGUE[f].census_member()
                ax.add_patch(Rectangle((j + 0.1, i + 0.12), 0.8, 0.76,
                                       facecolor=(fs.BLUES[3] if census
                                                  else fs.STATUS["absent"]),
                                       edgecolor="none", zorder=3))
    ax.set_xlim(0, len(shared)); ax.set_ylim(0, len(fams))
    ax.set_xticks([j + 0.5 for j in range(len(shared))])
    ax.set_xticklabels([f"{d}\n({counts[d]})" for d in shared],
                       fontsize=fs.FS_TICK - 1.6)
    ax.set_yticks([i + 0.5 for i in range(len(fams))])
    ax.set_yticklabels(fams, fontsize=fs.FS_TICK - 1.8)
    ax.invert_yaxis()
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("domain accession (number of families carrying it)",
                  fontsize=fs.FS_LABEL - 0.6)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=fs.BLUES[3], label="ion-channel family"),
                       Patch(facecolor=fs.STATUS["absent"],
                             label="catalogued, not an ion channel")],
              fontsize=fs.FS_TICK - 1.2, frameon=False,
              loc="upper left", bbox_to_anchor=(0, -0.055), ncol=2)
    fig.tight_layout()
    return fs.save(fig, OUT / "fig5_domain_matrix", legend_in="docs/review (numbered review figure)")[0]


def fig_traps() -> Path:
    import matplotlib.pyplot as plt
    rows = read_tsv(DATA / "domain_positions.tsv")
    by_gene: dict[str, list[dict]] = defaultdict(list)
    lengths: dict[str, int] = {}
    for r in rows:
        by_gene[r["gene"]].append(r)
        lengths[r["gene"]] = int(r["length_aa"])
    colours = domain_colours(rows)
    shorts = short_names()
    seen_domains: dict[str, str] = {}

    n = sum(len(g) for _, g, _ in TRAPS)
    fig, axes = plt.subplots(len(TRAPS), 1, figsize=(fs.W_FULL, 0.30 * n + 2.4),
                             gridspec_kw={"height_ratios": [len(g) for _, g, _ in TRAPS]})
    for ax, (title, genes, note) in zip(axes, TRAPS):
        widest = max(lengths.get(g, 1) for g in genes)
        for i, g in enumerate(genes):
            y = len(genes) - i - 1
            doms = [(int(r["start"]), int(r["end"]),
                     shorts.get(r["pfam"], r["pfam"]), r["pfam"])
                    for r in by_gene.get(g, [])]
            for _, _, nm, acc in doms:
                seen_domains[acc] = nm
            R.domain_bar(ax, y, lengths.get(g, 0), doms, colours=colours,
                         label_min_frac=0.085, scale=widest)
            ax.text(-widest * 0.035, y + 0.31, g, ha="right", va="center",
                    fontsize=fs.FS_TICK - 0.6)
            ax.text(lengths.get(g, 0) + widest * 0.012, y + 0.31,
                    f"{lengths.get(g, 0):,} aa", ha="left", va="center",
                    fontsize=fs.FS_TICK - 1.6, color=fs.MUTED)
        ax.set_xlim(-widest * 0.34, widest * 1.16)
        ax.set_ylim(-0.25, len(genes) + 0.55)
        ax.axis("off")
        ax.text(-widest * 0.20, len(genes) + 0.10, title, ha="left",
                va="bottom", fontsize=fs.FS_LABEL - 0.4, color=fs.INK)
        ax.text(widest * 1.16, len(genes) + 0.10, note, ha="right",
                va="bottom", fontsize=fs.FS_TICK - 1.4, color=fs.MUTED,
                style="italic")
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=colours[a],
                     label=(a if nm == a else f"{a}  {nm}"))
               for a, nm in sorted(seen_domains.items(), key=lambda kv: kv[1])]
    axes[-1].legend(handles=handles, fontsize=fs.FS_TICK - 2.0, frameon=False,
                    loc="upper left", bbox_to_anchor=(0.0, -0.06), ncol=4,
                    handlelength=1.1, handleheight=0.9, columnspacing=1.0)
    fig.tight_layout()
    return fs.save(fig, OUT / "fig8_architecture_traps", legend_in="docs/review (numbered review figure)")[0]


def main() -> int:
    argparse.ArgumentParser(description=__doc__).parse_args()
    fs.use()
    for f in (fig_lengths, fig_domain_matrix, fig_traps):
        print(f"[figs] {f()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
