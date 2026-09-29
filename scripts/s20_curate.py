"""s20_curate.py — S20 step 2: what the uncatalogued list genes are, and what
each database channelome is made of.

    python3 scripts/s20_curate.py

**The curation is by hand and says so.** Every gene a database list carries
and the catalogue does not name gets one class from `CURATED` below, with a
one-line reason drawn from its UniProt record name or HGNC locus type
(archived by `s20_lists.py`). Provenance `CURATED` (catalogue grammar): read,
not yet checked against primary literature. A gene with no entry here is a
hard error — no uncatalogued gene is classified by default.

→ `uncatalogued.tsv`, `list_decomposition.tsv` (each list's total split into
the reasons a gene is on it), `census_missed_by_list.tsv` (census pore genes
a list leaves out, by family).
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s20_lib import LISTS, OUT_DIR  # noqa: E402

#: class → meaning. `auxiliary_uncatalogued` and `pore_candidate` are
#: catalogue gaps; the rest are outside the D23 definition.
CLASSES = {
    "pseudogene": "not a protein-coding gene (HGNC locus type) or a UniProt 'putative' pseudogene product",
    "auxiliary_uncatalogued": "an auxiliary subunit of a catalogued channel that the catalogue does not list",
    "pore_candidate": "a proposed pore-forming channel not in the catalogue",
    "transporter_or_enzyme": "a transporter or enzyme with a reported channel-like activity",
    "paracellular": "forms paracellular (tight-junction) pores, not a transmembrane channel",
    "out_of_scope_uncatalogued": "out of scope by D23 (water channel) and not listed in the catalogue",
}

CURATED: dict[str, tuple[str, str]] = {
    "TRPC2": ("pseudogene", "HGNC locus type pseudogene; functional in rodents"),
    "GJA6P": ("pseudogene", "HGNC locus type pseudogene"),
    "FXYD6P3": ("pseudogene", "UniProt 'putative FXYD domain-containing ion transport regulator 8'"),
    "KCNIP1": ("auxiliary_uncatalogued", "Kv channel-interacting protein (KChIP), Kv4 auxiliary"),
    "KCNIP2": ("auxiliary_uncatalogued", "KChIP2, Kv4 auxiliary"),
    "KCNIP3": ("auxiliary_uncatalogued", "KChIP3, Kv4 auxiliary"),
    "KCNIP4": ("auxiliary_uncatalogued", "KChIP4, Kv4 auxiliary"),
    "TMEM37": ("auxiliary_uncatalogued", "voltage-dependent calcium channel gamma-like subunit"),
    "PACC1": ("pore_candidate", "proton-activated chloride channel (PAC / ASOR, TMEM206)"),
    "TMCO1": ("pore_candidate", "calcium load-activated calcium channel (CLAC)"),
    "TMEM87A": ("pore_candidate", "Golgi-pH regulating cation channel (Elkin1)"),
    "TMEM109": ("pore_candidate", "voltage-gated cation channel TMEM109 (mitsugumin-23)"),
    "CLCC1": ("pore_candidate", "ER anion channel 1 (chloride channel CLIC-like 1)"),
    "CCDC51": ("pore_candidate", "mitochondrial potassium channel (MITOK)"),
    "GPHRA": ("pore_candidate", "Golgi pH regulator A (GPR89A), reported anion channel"),
    "GPHRB": ("pore_candidate", "Golgi pH regulator B (GPR89B), reported anion channel"),
    "MFSD8": ("transporter_or_enzyme", "MFS transporter CLN7, reported lysosomal Cl- channel"),
    "UCP1": ("transporter_or_enzyme", "mitochondrial carrier SLC25A7, H+ leak"),
    "SLC17A6": ("transporter_or_enzyme", "vesicular glutamate transporter 2, Cl- conductance"),
    "SLC17A7": ("transporter_or_enzyme", "vesicular glutamate transporter 1, Cl- conductance"),
    "SLC17A8": ("transporter_or_enzyme", "vesicular glutamate transporter 3, Cl- conductance"),
    "NOX5": ("transporter_or_enzyme", "NADPH oxidase 5, H+ conductance"),
    "CYBB": ("transporter_or_enzyme", "NADPH oxidase 2 (gp91phox), H+ conductance"),
    "STING1": ("transporter_or_enzyme", "innate-immunity adaptor STING, reported H+ channel"),
    "CLDN4": ("paracellular", "claudin-4, tight-junction pore"),
    "CLDN17": ("paracellular", "claudin-17, tight-junction anion pore"),
    "AQP12B": ("out_of_scope_uncatalogued", "aquaporin-12B; the catalogue lists AQP12A only"),
}

#: The categories each list's total decomposes into, in reporting order.
PARTS = ["pore_census", "auxiliary", "auxiliary_uncatalogued", "out_of_scope",
         "out_of_scope_uncatalogued", "transporter", "transporter_or_enzyme",
         "pore_candidate", "paracellular", "pseudogene"]


def main() -> int:
    genes = read_tsv(OUT_DIR / "channelome_genes.tsv")
    unc = [g for g in genes if g["category"] == "uncatalogued"]
    missing = sorted(g["symbol"] for g in unc if g["symbol"] not in CURATED)
    if missing:
        raise SystemExit(f"uncatalogued genes with no curated class: {missing}")
    stale = sorted(set(CURATED) - {g["symbol"] for g in unc})
    if stale:
        raise SystemExit(f"curated entries no list carries any more: {stale}")
    rows = [{"symbol": g["symbol"], "hgnc_id": g["hgnc_id"],
             **{k: g[k] for k in LISTS}, "class": CURATED[g["symbol"]][0],
             "reason": CURATED[g["symbol"]][1], "provenance": "CURATED"}
            for g in sorted(unc, key=lambda g: (CURATED[g["symbol"]][0], g["symbol"]))]
    write_tsv(OUT_DIR / "uncatalogued.tsv", list(rows[0]), rows)

    def part(g: dict) -> str:
        return CURATED[g["symbol"]][0] if g["category"] == "uncatalogued" else g["category"]

    dec = []
    for name in (*LISTS, "any_list"):
        sel = [g for g in genes if (g["n_lists"] != "0" if name == "any_list" else g[name] == "1")]
        c = Counter(part(g) for g in sel)
        unknown = set(c) - set(PARTS)
        if unknown:
            raise SystemExit(f"{name}: unplaced categories {unknown}")
        n_pore = sum(1 for g in genes if g["category"] == "pore_census")
        dec.append({"list": name, "total": len(sel), **{k: c[k] for k in PARTS},
                    "census_missed": n_pore - c["pore_census"],
                    "not_pore_census": len(sel) - c["pore_census"],
                    "auxiliary_share": round((c["auxiliary"] + c["auxiliary_uncatalogued"])
                                             / len(sel), 4)})
    write_tsv(OUT_DIR / "list_decomposition.tsv", list(dec[0]), dec)

    miss = defaultdict(lambda: {k: [] for k in LISTS})
    for g in genes:
        if g["category"] == "pore_census":
            for k in LISTS:
                if g[k] == "0":
                    miss[g["family"]][k].append(g["symbol"])
    mrows = [{"family": f, **{f"{k}_missed": len(v[k]) for k in LISTS},
              **{f"{k}_genes": ",".join(sorted(v[k])) for k in LISTS}}
             for f, v in sorted(miss.items(), key=lambda x: -sum(map(len, x[1].values())))]
    write_tsv(OUT_DIR / "census_missed_by_list.tsv", list(mrows[0]), mrows)
    for d in dec:
        print(d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
