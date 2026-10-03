"""s10_repertoire.py — S10a driver: species tree, characters, reconstruction, overlay (D50).

    python3 scripts/s10_repertoire.py tree          # NCBI lineages → species tree
    python3 scripts/s10_repertoire.py characters    # S4b matrix → 0 / 1 / missing
    python3 scripts/s10_repertoire.py reconstruct   # Sankoff at g = 1, 2, Dollo-like
    python3 scripts/s10_repertoire.py overlay       # S5's controlled cells on every loss
    python3 scripts/s10_repertoire.py all

Every rule is D50's, fixed before any family was reconstructed. Output →
`results/repertoire/`; the archived taxonomy → `<data root>/raw_api/s10/taxonomy/`.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import date

from s10_lib import (BUSCO_FLOOR, GAINS, LIVE, OUT, PRIMARY, ROOT, ROOT_TAXON,
                     build_tree, events, fetch_lineages, newick, sankoff)

sys.path.insert(0, str(ROOT / "scripts"))
from s0_lib import live_progress, write_tsv  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402

DENSE = ROOT / "results" / "panel_density"
S4 = ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv"
CELLS = ROOT / "results" / "genome_sweep" / "cells.tsv"
GENOME_EVIDENCE = {"present", "genome_found", "genome_present"}


def _order_lineage(lin, order):
    for i, (_, nm, rk) in enumerate(lin):
        if rk == "order" and nm == order:
            return lin[:i + 1]
    return None


def load_tree():
    panel = read_tsv(DENSE / "order_panel.tsv")
    s4 = [r for r in read_tsv(S4) if r["group"] not in ("prokaryote", "virus")]
    s4_ids = [r["proteome_taxid"] or r["panel_taxid"] for r in s4]
    lins = fetch_lineages([r["taxid"] for r in panel] + s4_ids)
    tips, problems = {}, []
    for r in panel:
        lin = _order_lineage(lins[r["taxid"]], r["order"])
        if lin is None:
            problems.append((r["order"], r["taxid"], "order not on NCBI lineage"))
            continue
        tips[r["order"]] = lin
    tree = build_tree(tips)
    return tree, tips, lins, panel, s4, problems


def cmd_tree(_a) -> None:
    tree, tips, lins, panel, s4, problems = load_tree()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "species_tree.nwk").write_text(newick(tree) + "\n")
    rows = []
    for n in tree.walk():
        rows.append([n.name, n.taxid, n.rank or "", "tip" if n.tip else "internal",
                     len(n.children), len(n.tips())])
    write_tsv(OUT / "species_tree.tsv", ["node", "taxid", "rank", "kind", "children", "tips"], rows)
    # the S4 species placed on the order tips (for the overlay)
    place = []
    order_of = {o: o for o in tips}
    for r in s4:
        lin = lins[r["proteome_taxid"] or r["panel_taxid"]]
        orders = [nm for _, nm, rk in lin if rk == "order"]
        o = orders[0] if orders else ""
        place.append([r["species"], r["group"], o, "1" if o in order_of else "0"])
    write_tsv(OUT / "s4_species_orders.tsv", ["species", "group", "order", "on_tree"], place)
    internal = [n for n in tree.walk() if not n.tip]
    summ = {"retrieved": date.today().isoformat(), "source": "NCBI Taxonomy efetch LineageEx",
            "tips": len(tree.tips()), "internal_nodes": len(internal),
            "polytomies": sum(1 for n in internal if len(n.children) > 2),
            "root": tree.name, "root_children": sorted(c.name for c in tree.children),
            "max_children": max(len(n.children) for n in internal),
            "orders_unplaced": problems,
            "s4_species_on_tree": sum(1 for p in place if p[3] == "1"), "s4_species": len(place)}
    (OUT / "species_tree.json").write_text(json.dumps(summ, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summ.items() if k != "root_children"}, indent=1))
    print("root children:", len(summ["root_children"]))


def cmd_characters(_a) -> None:
    panel = {r["order"]: r for r in read_tsv(DENSE / "order_panel.tsv")}
    rows, tally = [], Counter()
    for r in read_tsv(DENSE / "order_matrix.tsv"):
        busco = float(panel[r["order"]]["busco_c"] or 0)
        hi, anyc = int(r["high"]), int(r["any"])
        if hi > 0:
            st, why = "1", "high-confidence call"
        elif anyc > 0:
            st, why = "?", "medium-only call"
        elif busco < BUSCO_FLOOR:
            st, why = "?", f"no call, BUSCO {busco:.1f} < {BUSCO_FLOOR:.0f}"
        else:
            st, why = "0", "no call"
        tally[(r["family"], st)] += 1
        rows.append([r["family"], r["order"], r["kingdom"], hi, anyc, f"{busco:.1f}", st, why])
    write_tsv(OUT / "characters.tsv",
              ["family", "order", "kingdom", "high", "any", "busco_c", "state", "why"], rows)
    fams = sorted({f for f, _ in tally})
    print(len(fams), "families;", Counter(r[6] for r in rows))


def _characters():
    ch = defaultdict(dict)
    for r in read_tsv(OUT / "characters.tsv"):
        ch[r["family"]][r["order"]] = None if r["state"] == "?" else int(r["state"])
    return ch


def cmd_reconstruct(_a) -> None:
    tree = load_tree()[0]
    ch = _characters()
    fam_rows, ev_rows, node_rows = [], [], []
    stated = defaultdict(set)                   # (family, child, event) → costs stating it
    steps = [(f"reconstruct {k}", False) for k in GAINS]
    for i, (key, g) in enumerate(GAINS.items()):
        for fam in sorted(ch):
            res = sankoff(tree, ch[fam], g)
            ev = events(res)
            for e in ev:
                if e["stated"]:
                    stated[(fam, e["child"], e["event"])].add(key)
                n = e["node"]
                ev_rows.append([fam, key, e["parent"], e["child"], n.rank or "", len(n.tips()),
                                e["event"], e["stated"]])
            c = Counter((e["event"], e["stated"]) for e in ev)
            st = ch[fam]
            fam_rows.append([fam, key, res["score"] if g < 10 ** 6 else res["score"] % 10 ** 6,
                             {1: "present", 0: "absent", None: "ambiguous"}[res["root"]],
                             sum(1 for v in st.values() if v == 1), sum(1 for v in st.values() if v == 0),
                             sum(1 for v in st.values() if v is None),
                             c[("gain", 1)], c[("loss", 1)],
                             sum(v for (k, s), v in c.items() if not s)])
            for n in tree.walk():
                if not n.tip:
                    node_rows.append([fam, key, n.name, n.rank or "", len(n.tips()),
                                      {1: "1", 0: "0", None: "?"}[res["node"][id(n)]]])
        steps[i] = (steps[i][0], True)
        live_progress(LIVE, "S10a", steps)
    write_tsv(OUT / "families.tsv", ["family", "cost", "score", "leca", "present", "absent", "missing",
                                     "gains", "losses", "ambiguous_edges"], fam_rows)
    for r in ev_rows:
        r.append("1" if r[7] and len(stated[(r[0], r[3], r[6])]) == len(GAINS) else "0")
    write_tsv(OUT / "events.tsv", ["family", "cost", "parent", "child", "child_rank", "child_tips",
                                   "event", "stated", "robust"], ev_rows)
    write_tsv(OUT / "nodes.tsv", ["family", "cost", "node", "rank", "tips", "state"], node_rows)
    print(len(fam_rows), "family×cost rows;", len(ev_rows), "events")


def _s5_verdict(r: dict) -> str:
    """S5's `present` is any-confidence; D50's presence is high-confidence
    (D50 addendum (a)), so a medium-only proteome cell reads `present_medium`."""
    if r["verdict"] == "present" and int(r["proteome_high"] or 0) == 0:
        return "present_medium"
    return r["verdict"]


def absent_region(node, family: str, state: dict) -> set[str]:
    """Tips below a loss that the reconstruction leaves absent: descend until a
    node whose primary state is not 0 (a regain, or ambiguous) — D50 addendum (b)."""
    if node.tip and not node.children:
        return {node.tip}
    if node is not None and state.get((family, node.name), "0") != "0":
        return set()
    return set().union(*(absent_region(c, family, state) for c in node.children))


def cmd_overlay(_a) -> None:
    tree = load_tree()[0]
    by_name = {n.name: n for n in tree.walk()}
    sp_order = {r["species"]: r["order"] for r in read_tsv(OUT / "s4_species_orders.tsv")
                if r["on_tree"] == "1"}
    cells = {(r["species"], r["family"]): _s5_verdict(r) for r in read_tsv(CELLS)}
    state = {(r["family"], r["node"]): r["state"] for r in read_tsv(OUT / "nodes.tsv")
             if r["cost"] == PRIMARY}
    rows = []
    for e in read_tsv(OUT / "events.tsv"):
        if e["cost"] != PRIMARY or e["event"] != "loss" or e["stated"] != "1":
            continue
        tips = absent_region(by_name[e["child"]], e["family"], state)
        sp = sorted(s for s, o in sp_order.items() if o in tips)
        v = {s: cells.get((s, e["family"]), "") for s in sp}
        if any(x in GENOME_EVIDENCE for x in v.values()):
            strength = "contradicted"
        elif any(x == "absent" for x in v.values()):
            strength = "controlled"
        else:
            strength = "proteome-only"
        rows.append([e["family"], e["parent"], e["child"], e["child_rank"], e["child_tips"],
                     e["robust"], strength, len(sp),
                     "; ".join(f"{s}={x or 'none'}" for s, x in v.items())])
    write_tsv(OUT / "losses.tsv", ["family", "parent", "child", "child_rank", "child_tips", "robust",
                                   "strength", "s4_species", "s5_verdicts"], rows)
    print(len(rows), "stated losses (primary);", Counter(r[6] for r in rows))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["tree", "characters", "reconstruct", "overlay", "all"])
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cmds = [cmd_tree, cmd_characters, cmd_reconstruct, cmd_overlay]
    for c in (cmds if a.cmd == "all" else [globals()[f"cmd_{a.cmd}"]]):
        c(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
