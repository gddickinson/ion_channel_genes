"""s4b_report.py — S4b's report, instrument check and figure (D13).

    bin/envpy scripts/s4b_report.py

* **Instrument check**: the 35 orders whose S4b proteome is also an S4 panel
  proteome — high-confidence family presence on the dense panel against the
  S3b panel profile calls (`<data root>/hmmer/s3b/panel_profile_calls.tsv.gz`, same
  proteome, same profiles; only `-Z` differs) — S3b's high-confidence
  profile calls per species, S4b's definition exactly.
* **The matrix**: census family × order presence (high-confidence proteome
  calls), drawn with orders grouped by kingdom and phylum.

→ `results/panel_density/report.md`, `instrument_check.tsv`,
`figures/order_matrix.png`.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

D = ROOT / "results" / "panel_density"
KINGDOM_ORDER = ["Metazoa", "Fungi", "Viridiplantae", ""]


def instrument_check(mx: list[dict]) -> dict:
    s4 = {r["upid"]: r["species"] for r in read_tsv(ROOT / "results" / "proteome_scope" /
                                                    "proteome_manifest.tsv")
          if r["status"] == "proteome"}
    panel = {r["upid"]: r for r in read_tsv(D / "order_panel.tsv") if r["in_s4_panel"] == "1"}
    # S3b's own high-confidence profile calls per species — the same
    # definition S4b's matrix uses (family call, high confidence)
    import gzip
    from s3b_lib import s3b_dir
    sp_of = {}
    with gzip.open(s3b_dir() / "panel_universe.tsv.gz", "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp_of[f[0]] = f[3]
    s3b = defaultdict(bool)
    with gzip.open(s3b_dir() / "panel_profile_calls.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["p_call"] == "family" and r["p_confidence"] == "high":
                s3b[(sp_of.get(r["target"], ""), r["p_family"])] = True
    rows, agree, total = [], 0, 0
    for r in mx:
        if r["upid"] not in panel:
            continue
        sp = s4[r["upid"]]
        a, b = r["present"] == "1", s3b[(sp, r["family"])]
        total += 1
        agree += a == b
        if a != b:
            rows.append({"species": sp, "family": r["family"], "s4b_present": int(a),
                         "s3b_present": int(b)})
    write_tsv(D / "instrument_check.tsv", ["species", "family", "s4b_present", "s3b_present"],
              rows)
    return {"cells": total, "agree": agree, "disagree": total - agree,
            "proteomes": len({r["upid"] for r in mx if r["upid"] in panel})}


def figure(mx: list[dict]) -> None:
    import figstyle as fs
    import matplotlib.pyplot as plt
    fs.use()
    orders = {}
    for r in mx:
        orders[r["order"]] = (r["kingdom"], r["phylum"], r["class"])
    okeys = sorted(orders, key=lambda o: (KINGDOM_ORDER.index(orders[o][0])
                                          if orders[o][0] in KINGDOM_ORDER else 9,
                                          orders[o][1], orders[o][2], o))
    sfo = {sf: i for i, sf in enumerate(fs.SUPERFAMILY_ORDER)}
    fams = sorted({r["family"] for r in mx},
                  key=lambda f: (sfo.get(CATALOGUE[f].superfamily, 99),
                                 CATALOGUE[f].superfamily, f))
    pres = {(r["order"], r["family"]) for r in mx if r["present"] == "1"}
    fig, ax = plt.subplots(figsize=(fs.W_FULL, 6.2))
    xs, ys = [], []
    for j, o in enumerate(okeys):
        for i, f in enumerate(fams):
            if (o, f) in pres:
                xs.append(j)
                ys.append(i)
    ax.scatter(xs, ys, s=1.1, marker="s", color=fs.BLUES[5], linewidth=0)
    prev = None
    for j, o in enumerate(okeys):
        k = orders[o][0] or "protists"
        if k != prev:
            ax.axvline(j - 0.5, color=fs.INK, lw=0.5)
            ax.text(j + 2, -1.5, k, fontsize=fs.FS_NOTE, va="bottom")
            prev = k
    ax.set_xlim(-0.5, len(okeys) - 0.5)
    ax.set_ylim(len(fams) - 0.5, -0.5)
    ax.set_yticks(range(len(fams)), fams, fontsize=fs.FS_TICK - 2.6)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(okeys)} eukaryotic orders (one reference proteome each), "
                  "grouped by kingdom and phylum", fontsize=fs.FS_LABEL)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    print("\n".join(map(str, fs.save(fig, D / "figures" / "order_matrix"))))


def main() -> int:
    mx = read_tsv(D / "order_matrix.tsv")
    sel = json.loads((D / "order_panel.json").read_text())
    db = json.loads((D / "dense_db.json").read_text())
    summ = json.loads((D / "order_matrix.json").read_text())
    chk = instrument_check(mx)
    by_k = Counter(r["kingdom"] or "protists" for r in read_tsv(D / "order_panel.tsv"))
    fam_orders = summ["families_by_orders_present"]
    absent_everywhere = sorted(f for f in {r["family"] for r in mx} if f not in fam_orders)
    L = ["# S4b — the dense panel: one reference proteome per eukaryotic order", "",
         "Rendered by `scripts/s4b_report.py` from `order_panel.tsv`, "
         "`order_panel_files.tsv`, `dense_db.json`, `order_matrix.tsv` and "
         "`instrument_check.tsv` (D13). Bulk: `<data root>/proteomes/s4b/`.", "",
         "![S4b](figures/order_matrix.png)", "",
         f"**{sel['orders']} orders** from {sel['release']}'s "
         f"{sel['eukaryotic_reference_proteomes']:,} eukaryotic reference proteomes "
         f"({', '.join(f'{k} {v}' for k, v in by_k.most_common())}); "
         f"{sel['from_s4_panel']} keep their S4-panel proteome. "
         f"{sel['without_order_rank']} proteomes have no NCBI order rank (mostly protists) "
         "and are outside the rule. All files MD5-verified against the release metalink; "
         f"**{db['sequences']:,} sequences**, entry counts equal to the README's for all "
         "but human (the D35 reissue).", "",
         f"**Instrument check** — the {chk['proteomes']} proteomes shared with the S4 "
         f"panel: high-confidence family presence agrees with S3b on "
         f"**{chk['agree']:,} / {chk['cells']:,}** cells (`instrument_check.tsv` lists "
         "the rest; same proteome and profiles, only `-Z` and r4–r6's later profile "
         "versions differ).", "",
         f"**The matrix**: {summ['present_cells']:,} present cells (high-confidence "
         f"proteome calls) over {summ['orders']} orders × {summ['families']} census "
         f"families. Absent from every eukaryotic order: {', '.join(absent_everywhere) or 'none'}.",
         "", "| family | orders present | | family | orders present |",
         "|---|---|---|---|---|"]
    items = sorted(fam_orders.items(), key=lambda x: -x[1])
    half = (len(items) + 1) // 2
    for a, b in zip(items[:half], items[half:] + [("", "")]):
        L.append(f"| {a[0]} | {a[1]} | | {b[0]} | {b[1]} |")
    L += ["", "**What these cells are (D46).** Each is a proteome call: present = a "
          "high-confidence S3a profile call in that order's reference proteome. An "
          "empty cell is an **annotation-level absence** — no S5 genome check, no "
          "matched-bait control. S10 states absences as *controlled* (S5, the 52 "
          "species) or *proteome-only* (here), never merged.", ""]
    (D / "report.md").write_text("\n".join(L))
    figure(mx)
    print(json.dumps(chk))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
