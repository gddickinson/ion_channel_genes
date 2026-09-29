"""s20_report.py — render results/auxiliary/report.md from S20's tables (D13).

    python3 scripts/s20_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv  # noqa: E402
from s20_curate import CLASSES, PARTS  # noqa: E402
from s20_lib import LISTS, OUT_DIR  # noqa: E402

NAMES = {"gtopdb": "GtoPdb", "hgnc": "HGNC", "uniprot": "UniProt KW-0407",
         "any_list": "union"}


def sec_lists(L: list[str]) -> None:
    meta = json.loads((OUT_DIR / "lists.json").read_text())
    dec = {r["list"]: r for r in read_tsv(OUT_DIR / "list_decomposition.tsv")}
    genes = read_tsv(OUT_DIR / "channelome_genes.tsv")
    n_pore = sum(g["category"] == "pore_census" for g in genes)
    u = dec["any_list"]
    L += ["## 1. Three database channelomes, one denominator", "",
          "Three human ion-channel lists maintained by databases, each archived with "
          f"its release under `<data root>/raw_api/s20/`: **GtoPdb** ({meta['gtopdb']}; "
          "target types vgic, lgic, other_ic), **HGNC** (gene group 177 'Ion channels' "
          f"and every group below it; {meta['hgnc']}) and **UniProt** (reviewed human "
          f"entries with keyword KW-0407 'Ion channel', {meta['uniprot']} — the release "
          "census v2 was enumerated from). Every gene is joined on its HGNC id; the "
          "catalogue's own symbols resolved 468/468 (renamed: "
          + ", ".join(f"{k} → {v}" for k, v in meta["catalogue_symbols_renamed"].items())
          + ").", "",
          f"**The three lists hold {dec['gtopdb']['total']}, {dec['hgnc']['total']} and "
          f"{dec['uniprot']['total']} genes, and {u['total']} together** — the published "
          "range of the human channelome, reproduced by three curated databases. "
          f"Of the union, **{u['pore_census']} are this catalogue's pore-forming census "
          f"genes** ({n_pore} in the catalogue; {u['census_missed']} are on no list) and "
          f"**{u['not_pore_census']} are not**:", "",
          "| list | total | pore (census) | " + " | ".join(p.replace("_", " ") for p in PARTS[1:])
          + " | census pore genes missed |", "|---|" + "---|" * (len(PARTS) + 2)]
    for k in (*LISTS, "any_list"):
        r = dec[k]
        L.append(f"| {NAMES[k]} | {r['total']} | " + " | ".join(r[p] for p in PARTS)
                 + f" | {r['census_missed']} |")
    L += ["", "Column meanings — *auxiliary*: a catalogued `channel_associated` subunit; "
          "*auxiliary uncatalogued*, *pore candidate*, *transporter or enzyme*, "
          "*paracellular*, *pseudogene*, *out of scope uncatalogued*: the hand-curated "
          "classes of §4; *out of scope*: catalogued aquaporins; *transporter*: "
          "catalogued CLC and SLC26 transporters.", ""]


def sec_spread(L: list[str]) -> None:
    dec = {r["list"]: r for r in read_tsv(OUT_DIR / "list_decomposition.tsv")}
    miss = read_tsv(OUT_DIR / "census_missed_by_list.tsv")
    genes = read_tsv(OUT_DIR / "channelome_genes.tsv")
    core = sum(g["n_lists"] == "3" for g in genes)
    core_pore = sum(g["n_lists"] == "3" and g["category"] == "pore_census" for g in genes)
    L += ["## 2. Where the spread comes from", "",
          f"**{core} genes are on all three lists, and all {core_pore} of them are "
          "pore-forming census genes** — the agreed core is pure. The lists differ on "
          "two independent axes, and both are scope, not biology:", "",
          "* **What they add beyond the pore-forming genes**: "
          + "; ".join(f"{NAMES[k]} {dec[k]['not_pore_census']}" for k in LISTS)
          + " — auxiliary subunits, aquaporins, CLC transporters, and (UniProt only) "
          "enzymes, transporters and claudins carrying the 'Ion channel' keyword.",
          "* **Which contested or large-pore families they leave out**: "
          + "; ".join(f"{NAMES[k]} misses {dec[k]['census_missed']} census genes"
                      for k in LISTS) + ".", "",
          "| family | " + " | ".join(f"{NAMES[k]} missed" for k in LISTS) + " | genes |",
          "|---|---|---|---|---|"]
    for r in miss:
        genes = sorted({g for k in LISTS for g in r[f"{k}_genes"].split(",") if g})
        L.append(f"| {r['family']} | " + " | ".join(r[f"{k}_missed"] for k in LISTS)
                 + f" | {', '.join(genes)} |")
    L += ["", "Two of these are single annotation choices with large effects: "
          "**UniProt's KW-0407 is on none of the 21 connexins** (their keyword is 'Gap "
          "junction'), and **GtoPdb lists no TMEM16 scramblase, TMC, CLIC, tweety, "
          "OSCA or otopetrin** — the contested families of the scope document.", ""]


def sec_aux(L: list[str]) -> None:
    dec = {r["list"]: r for r in read_tsv(OUT_DIR / "list_decomposition.tsv")}
    genes = read_tsv(OUT_DIR / "channelome_genes.tsv")
    aux = [g for g in genes if g["category"] == "auxiliary"]
    unc = [r for r in read_tsv(OUT_DIR / "uncatalogued.tsv")
           if r["class"] == "auxiliary_uncatalogued"]
    n_pore = sum(g["category"] == "pore_census" for g in genes)
    fams = Counter(g["family"] for g in aux)
    n_all = len(aux) + len(unc)
    L += ["## 3. The auxiliary subunits — how many, and who counts them", "",
          f"**The catalogue names {len(aux)} human auxiliary genes in {len(fams)} "
          f"families; the lists add {len(unc)} it does not** ("
          + ", ".join(r["symbol"] for r in unc) + "). "
          f"Counted as channels, all {n_all} would inflate the {n_pore} pore-forming "
          f"genes by **{n_all / n_pore:.0%}** — not the 'roughly 15 %' the scope "
          "document stated before this measurement.", "",
          "**How often the lists count them**: " + "; ".join(
              f"{NAMES[k]} {int(dec[k]['auxiliary']) + int(dec[k]['auxiliary_uncatalogued'])}"
              f" ({float(dec[k]['auxiliary_share']):.1%} of its total)"
              for k in (*LISTS, "any_list")) + ". Each list counts a different set:", "",
          "| auxiliary family | genes | GtoPdb | HGNC | UniProt |", "|---|---|---|---|---|"]
    by = defaultdict(list)
    for g in aux:
        by[g["family"]].append(g)
    for f, gs in sorted(by.items()):
        L.append(f"| {f} | {len(gs)} | " + " | ".join(
            str(sum(g[k] == "1" for g in gs)) for k in LISTS) + " |")
    L.append("")


def sec_uncat(L: list[str]) -> None:
    rows = read_tsv(OUT_DIR / "uncatalogued.tsv")
    c = Counter(r["class"] for r in rows)
    L += ["## 4. Genes a list carries and the catalogue does not name", "",
          f"**{len(rows)} genes**, each classified by hand (`s20_curate.CURATED`, "
          "provenance CURATED — read from the UniProt record name or HGNC locus type, "
          "not yet checked against primary literature):", "",
          "| class | genes | meaning |", "|---|---|---|"]
    for k, meaning in CLASSES.items():
        L.append(f"| {k} | {c[k]} | {meaning} |")
    L += ["", "| gene | lists | class | reason |", "|---|---|---|---|"]
    for r in rows:
        lists = ", ".join(NAMES[k] for k in LISTS if r[k] == "1")
        L.append(f"| {r['symbol']} | {lists} | {r['class']} | {r['reason']} |")
    L += ["", "The two catalogue gaps are emergent rows, not edits made here: the "
          "`auxiliary_uncatalogued` genes (KChIP1–4, TMEM37) belong in "
          "`channel_associated` families, and the `pore_candidate` genes — PACC1 above "
          "all, a proton-activated Cl⁻ channel with solved structures — are candidate "
          "census families, which is a change to the census search space (D34).", ""]


def sec_groups(L: list[str]) -> None:
    groups = read_tsv(OUT_DIR / "aux_groups.tsv")
    fams = Counter(r["family"] for r in groups)
    pooled = [f for f, n in fams.items() if n > 1]
    L += ["## 5. The auxiliary families are grouped by partner, not by descent", "",
          f"All-against-all `phmmer` among each family's human members (E ≤ 1e-3 = "
          f"homologous; groups are connected components): **{len(pooled)} of "
          f"{len(fams)} auxiliary families pool unrelated proteins**, "
          f"{len(groups)} homology groups in all. A profile built across a pooled "
          "family is not a detector of any of its parts (S3a: the LOO decoys got no "
          "hit), so the panel census below counts homology groups, not families.", "",
          "| family | groups | homology groups (human genes) |", "|---|---|---|"]
    for f in sorted(fams):
        gs = [r for r in groups if r["family"] == f]
        txt = " · ".join(r["genes"].replace(",", ", ") for r in gs)
        nt = gs[0]["not_tested"]
        L.append(f"| {f} | {len(gs)} | {txt}{f' (not tested: {nt})' if nt else ''} |")
    L.append("")


def sec_panel(L: list[str]) -> None:
    p = OUT_DIR / "aux_panel.tsv"
    if not p.exists():
        L += ["## 6. Auxiliaries across the panel", "", "*Not yet run.*", ""]
        return
    rows = read_tsv(p)
    c = Counter(r["homology_group"] for r in rows)
    placed = [r for r in rows if r["homology_group"].startswith("assoc_")]
    L += ["## 6. Auxiliaries across the panel", "",
          f"S3b's panel profile calls to an auxiliary family: **{len(rows)}**. Each was "
          "searched against the whole human proteome; it counts as an auxiliary of a "
          "group only if its **best human hit is a human member of that group** "
          "(reciprocal best hit, E ≤ 1e-5). "
          f"**{len(placed)} pass**; {c['not_auxiliary']} have a better human hit outside "
          f"the group, {c['unplaced']} no human hit.", "",
          "| homology group | human genes | calls | pass (high) | pass (medium) | "
          "not auxiliary | species (high) | panel groups (high) |",
          "|---|---|---|---|---|---|---|---|"]
    groups = {r["group"]: r for r in read_tsv(OUT_DIR / "aux_groups.tsv")}
    fam_rows = defaultdict(list)
    for r in rows:
        fam_rows[r["family"]].append(r)
    for g, gr in groups.items():
        fam = gr["family"]
        mine = [r for r in rows if r["homology_group"] == g]
        hi = [r for r in mine if r["confidence"] == "high"]
        L.append(f"| {g} | {gr['genes'].replace(',', ', ')} | "
                 f"{len(fam_rows[fam]) if g.endswith('.1') else ''} | {len(hi)} | "
                 f"{len(mine) - len(hi)} | "
                 f"{sum(r['homology_group'] != fam and not r['homology_group'].startswith(fam) for r in fam_rows[fam]) if g.endswith('.1') else ''} | "
                 f"{len({r['species'] for r in hi})} | "
                 f"{', '.join(sorted({r['group'] for r in hi}))} |")
    zero = [g for g, gr in groups.items()
            if not any(r["homology_group"] == g for r in rows)]
    vert = [g for g in groups if (lambda hi: hi and {r["group"] for r in hi} == {"vertebrate"})(
        [r for r in rows if r["homology_group"] == g and r["confidence"] == "high"])]
    L += ["", f"**Vertebrate-only in the panel ({len(vert)} groups)**: "
          + ", ".join(groups[g]["genes"].split(",")[0] for g in vert)
          + ". **Wider**: " + ", ".join(
              f"{groups[g]['genes'].split(',')[0]} ("
              + ", ".join(sorted({r['group'] for r in rows if r['homology_group'] == g
                                  and r['confidence'] == 'high'})) + ")"
              for g in groups if g not in vert and g not in zero) + ".",
          "", "Two readings these rows do **not** support. *Vertebrate-only* is where a "
          "human-seeded instrument found the group, not an absence claim: no zero cell "
          "here has passed D37's genome test (CatSper auxiliaries, for one, are reported "
          "outside vertebrates *(pending: S10)*). And a reciprocal best human hit shows "
          "**homology, not function**: the Kvβ group's plant and prokaryotic members are "
          "aldo-keto reductases whose closest human relative happens to be Kvβ, and "
          "nothing here says they serve a channel.",
          "", f"**No call at all for {len(zero)} groups** ("
          + ", ".join(groups[g]["genes"] for g in zero) + "): the pooled family profile "
          "does not reach them. That is a detection limit of a profile built across "
          "unrelated proteins (§5), not an absence — including in human, where these "
          "genes exist.", ""]
    worst = Counter(r["best_human_gene"] for r in rows if r["homology_group"] == "not_auxiliary")
    L += ["", "The calls that fail are the shared-domain drift S3b found for TRPN and "
          "LRRC8: their best human hits are "
          + ", ".join(f"{g} {n}" for g, n in worst.most_common(8) if g) + " ….", ""]


def main() -> int:
    L = ["# S20 — the auxiliary subunits and the published channelome", "",
         "Rendered by `scripts/s20_report.py` from the tables in this directory (D13). "
         "Raw lists, HGNC lookups and phmmer inputs: `<data root>/raw_api/s20/`.", "",
         "![S20](figures/auxiliary.png)", ""]
    sec_lists(L)
    sec_spread(L)
    sec_aux(L)
    sec_uncat(L)
    sec_groups(L)
    sec_panel(L)
    (OUT_DIR / "report.md").write_text("\n".join(L))
    print(f"→ {OUT_DIR / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
