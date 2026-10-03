"""s10b_figures.py — S10b's headline figure, drawn from the committed tables (D13).

    bin/envpy scripts/s10b_figures.py

`results/repertoire/figures/s10b_checks.png` (+ .pdf):

A. The embedded-contig test (D51 (1)): every animal locus or protein called
   MscS, MscL, KcsA-like or GluR0, by what the other genes on its contig
   look like.
B. *Daphnia pulex*: the 2011 and the current assembly's S5 verdict for every
   cell that was `absent` or changed.
C. S5's ZAC, PACC1 and CLCC1 absences against NCBI and Ensembl orthologs.
D. K2P presence per kingdom in S4b's orders, and what the plant calls are.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import figstyle as fs
from s3_hmm_lib import read_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "repertoire"
SHORT = {"mscs": "MscS", "mscl": "MscL", "kcsa_prok": "KcsA-like", "iglur_prok": "GluR0",
         "zac": "ZAC", "pacc": "PACC1", "clcc1": "CLCC1"}
STATUS = [("embedded", fs.SUPERFAMILY["ploop"], "embedded (a neighbour is animal)"),
          ("foreign", fs.SUPERFAMILY["cysloop"], "foreign (neighbours non-animal)"),
          ("no_evidence", fs.GRID, "no other gene on the contig"),
          ("failed", fs.INK, "contig not retrieved")]
VERDICT = {"absent": fs.SUPERFAMILY["cysloop"], "partial": fs.SUPERFAMILY["deg_enac"],
           "genome_found": fs.SUPERFAMILY["ploop"], "genome_weak": fs.BLUES[1],
           "no_locus_unrescued": fs.GRID, "gap": fs.SUPERFAMILY["p2x"]}
ORTH = {"agrees": fs.SUPERFAMILY["iglur"], "disputed": fs.SUPERFAMILY["cysloop"],
        "not_covered": fs.GRID}


def main() -> int:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Patch
    fs.use()
    ct = read_tsv(OUT / "contig_test.tsv")
    dc = read_tsv(OUT / "daphnia_compare.tsv")
    orth = read_tsv(OUT / "ortholog_check.tsv")
    kk = read_tsv(OUT / "k2p_by_kingdom.tsv")
    fig = plt.figure(figsize=(fs.W_FULL, 6.4))
    gs = GridSpec(2, 2, figure=fig, hspace=0.62, wspace=0.55,
                  left=0.13, right=0.98, top=0.9, bottom=0.08)

    ax = fig.add_subplot(gs[0, 0])
    groups = [("MscS, high", lambda r: r["family"] == "mscs" and r["confidence"] == "high"),
              ("MscS, medium", lambda r: r["family"] == "mscs" and r["confidence"] == "medium")]
    groups += [(SHORT[f], (lambda f: lambda r: r["family"] == f)(f))
               for f in ("mscl", "kcsa_prok", "iglur_prok") if any(r["family"] == f for r in ct)]
    for i, (_, sel) in enumerate(groups):
        c, left = Counter(r["status"] for r in ct if sel(r)), 0
        for key, col, _ in STATUS:
            if c[key]:
                ax.barh(i, c[key], left=left, color=col, height=0.62, lw=0)
                ax.text(left + c[key] / 2, i, str(c[key]), ha="center", va="center",
                        fontsize=fs.FS_NOTE, color="white" if key in ("embedded", "foreign",
                                                                     "failed") else fs.INK)
                left += c[key]
    ax.set_yticks(range(len(groups)), [g for g, _ in groups])
    ax.invert_yaxis()
    ax.set_xlabel("animal loci / proteins tested")
    fs.despine(ax)
    ax.legend(handles=[Patch(facecolor=c, label=l) for _, c, l in STATUS],
              fontsize=fs.FS_NOTE, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.panel(ax, "A", "Prokaryote-type channels in animals: what is on their contig")

    bx = fig.add_subplot(gs[0, 1])
    rows = [r for r in dc if r["old_verdict"] == "absent" or r["changed"] == "1"]
    for i, r in enumerate(rows):
        for j, k in enumerate(("old_verdict", "new_verdict")):
            bx.scatter(j, i, s=46, marker="s", color=VERDICT.get(r[k], fs.FAINT), lw=0)
    bx.set_yticks(range(len(rows)), [r["family"] for r in rows], fontsize=fs.FS_TICK)
    bx.set_xticks([0, 1], ["2011 assembly", "2021 assembly"])
    bx.set_xlim(-0.6, 1.6)
    bx.invert_yaxis()
    fs.despine(bx, keep=())
    bx.tick_params(length=0)
    used = sorted({r[k] for r in rows for k in ("old_verdict", "new_verdict")})
    bx.legend(handles=[Patch(facecolor=VERDICT.get(v, fs.FAINT), label=v) for v in used],
              fontsize=fs.FS_NOTE, loc="center left", bbox_to_anchor=(1.0, 0.5),
              frameon=False, handlelength=1.0)
    fs.panel(bx, "B", "Daphnia pulex: S5 verdicts on two assemblies")

    cx = fig.add_subplot(gs[1, 0])
    for i, r in enumerate(orth):
        cx.barh(i, 1, color=ORTH[r["status"]], height=0.8, lw=0)
        cx.text(1.05, i, f"NCBI {r['ncbi']}, Ensembl {r['ensembl'].replace('_', ' ')}",
                va="center", fontsize=fs.FS_NOTE, color=fs.MUTED)
    cx.set_yticks(range(len(orth)), [f"{SHORT[r['family']]} · {r['species']}" for r in orth],
                  fontsize=fs.FS_TICK)
    cx.set_xlim(0, 3.2)
    cx.set_xticks([])
    cx.invert_yaxis()
    fs.despine(cx, keep=())
    cx.legend(handles=[Patch(facecolor=c, label=k.replace("_", " ")) for k, c in ORTH.items()],
              fontsize=fs.FS_NOTE, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
              frameon=False, handlelength=1.0)
    fs.panel(cx, "C", "S5 absences against two ortholog databases")

    dx = fig.add_subplot(gs[1, 1])
    for i, r in enumerate(kk):
        left = 0
        for key, col in (("present", fs.SUPERFAMILY["ploop"]), ("absent", fs.SUPERFAMILY["cysloop"]),
                         ("missing", fs.GRID)):
            n = int(r[key])
            if n:
                dx.barh(i, n, left=left, color=col, height=0.62, lw=0)
            left += n
    dx.set_yticks(range(len(kk)), [r["kingdom"] or "protists" for r in kk])
    dx.invert_yaxis()
    dx.set_xlabel("S4b orders (one proteome each)")
    fs.despine(dx)
    dx.legend(handles=[Patch(facecolor=fs.SUPERFAMILY["ploop"], label="K2P called (high)"),
                       Patch(facecolor=fs.SUPERFAMILY["cysloop"], label="no call"),
                       Patch(facecolor=fs.GRID, label="missing (medium or low BUSCO)")],
              fontsize=fs.FS_NOTE, loc="lower right", frameon=True, framealpha=1,
              edgecolor=fs.GRID, handlelength=1.0)
    fs.panel(dx, "D", "K2P calls by kingdom")
    fig.suptitle("S10b — the animal MscS question and five absence checks",
                 fontsize=fs.FS_SUPTITLE, y=0.99)
    for p in fs.save(fig, OUT / "figures" / "s10b_checks"):
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
