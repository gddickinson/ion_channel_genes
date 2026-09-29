"""s20_lists.py — S20 step 1: the three database channelomes against the catalogue.

    python3 scripts/s20_lists.py [--refresh]

→ `results/auxiliary/channelome_genes.tsv` (one row per HGNC id in any list
or in the catalogue: which lists carry it, its catalogue family and
category) and `list_composition.tsv` (per list × category counts).

The catalogue's category is the only reading applied (`s20_lib.category`).
Genes a list carries and the catalogue does not name are `uncatalogued`;
what they are is S20 step 2 (`uncatalogued.tsv`, curated by hand, each with
a reason), never inferred from the list's own label.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import write_tsv  # noqa: E402
from s20_lib import (LISTS, OUT_DIR, catalogue_genes, category,  # noqa: E402
                     gtopdb_list, hgnc_list, hgnc_resolve, uniprot_list)

GENE_FIELDS = ["hgnc_id", "symbol", "catalogue_symbol", "resolved_by", "family",
               "superfamily", "status", "category", *LISTS, "n_lists",
               "gtopdb_group", "hgnc_group"]
CATS = ["pore_census", "auxiliary", "transporter", "non_channel_homolog",
        "out_of_scope", "uncatalogued"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    g, gv = gtopdb_list(a.refresh)
    h, hv = hgnc_list(a.refresh)
    u, uv, u_nohgnc = uniprot_list(a.refresh)
    lists = {"gtopdb": g, "hgnc": h, "uniprot": u}

    cat = catalogue_genes()
    res = hgnc_resolve(sorted({c["symbol"] for c in cat}), a.refresh)
    unresolved = sorted(s for s, r in res.items() if not r["hgnc_id"])
    by_id: dict[str, dict] = {}
    for c in cat:
        r = res[c["symbol"]]
        if not r["hgnc_id"]:
            continue
        if r["hgnc_id"] in by_id:
            raise RuntimeError(f"{c['symbol']} and {by_id[r['hgnc_id']]['symbol']} "
                               f"are one HGNC gene ({r['hgnc_id']})")
        by_id[r["hgnc_id"]] = {**c, "approved": r["approved"], "how": r["how"]}

    rows = []
    for hid in sorted(set(by_id) | set(g) | set(h) | set(u)):
        c = by_id.get(hid)
        sym = (c["approved"] if c else "") or next(
            (L[hid]["symbol"] for L in (h, g, u) if hid in L), "")
        member = {k: int(hid in L) for k, L in lists.items()}
        rows.append({"hgnc_id": hid, "symbol": sym,
                     "catalogue_symbol": c["symbol"] if c else "",
                     "resolved_by": c["how"] if c else "",
                     "family": c["family"] if c else "",
                     "superfamily": c["superfamily"] if c else "",
                     "status": c["status"] if c else "",
                     "category": category(c["status"], c["census"]) if c else "uncatalogued",
                     **member, "n_lists": sum(member.values()),
                     "gtopdb_group": g.get(hid, {}).get("group", ""),
                     "hgnc_group": h.get(hid, {}).get("group", "")})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT_DIR / "channelome_genes.tsv", GENE_FIELDS, rows)

    comp = []
    for name in (*LISTS, "any_list", "catalogue"):
        if name == "any_list":
            sel = [r for r in rows if r["n_lists"]]
        elif name == "catalogue":
            sel = [r for r in rows if r["category"] != "uncatalogued"]
        else:
            sel = [r for r in rows if r[name]]
        c = Counter(r["category"] for r in sel)
        comp.append({"list": name, "total": len(sel), **{k: c[k] for k in CATS}})
    write_tsv(OUT_DIR / "list_composition.tsv", ["list", "total", *CATS], comp)
    (OUT_DIR / "lists.json").write_text(json.dumps({
        "gtopdb": gv, "hgnc": hv, "uniprot": uv,
        "uniprot_without_hgnc": u_nohgnc,
        "catalogue_symbols_unresolved": unresolved,
        "catalogue_symbols_renamed": {s: r["approved"] for s, r in res.items()
                                      if r["hgnc_id"] and r["approved"] != s}},
        indent=1))
    for r in comp:
        print(r)
    print("unresolved catalogue symbols:", unresolved)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
