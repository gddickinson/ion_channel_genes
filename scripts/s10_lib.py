"""s10_lib.py — S10a shared pieces: the species tree and the reconstruction (D50).

**The species tree is an input (D15).** NCBI Taxonomy lineages (`efetch`
`LineageEx`) for the S4b order proteomes, archived under
`<data root>/raw_api/s10/`; tips are orders, internal nodes every named taxon
on the lineages with unary chains collapsed, polytomies kept.

**Sankoff parsimony, exact on multifurcations.** States 0/1; loss costs 1,
gain costs `g`, presence at the root costs one gain. A down pass gives each
node's subtree cost per state, an up pass the cost of everything outside it;
their sum is the best total with the node fixed in that state, so a node's
state (or an edge's state pair) is *stated* only when it is the unique
optimum (D50 (3)). Missing tips cost nothing in either state.
"""

from __future__ import annotations

import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "repertoire"
LIVE = ROOT / "results" / "session_live.json"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ROOT_TAXON = "Eukaryota"
BUSCO_FLOOR = 70.0
GAINS = {"g1": 1, "g2": 2, "dollo": 10 ** 6}
PRIMARY = "g2"
INF = 10 ** 15


def raw(*parts: str) -> Path:
    """`<data root>/raw_api/s10/...` — raises without the drive (D1)."""
    p = require_data_root() / "raw_api" / "s10"
    for part in parts:
        p = p / part
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------- taxonomy
def fetch_lineages(taxids: list[str], name: str = "ncbi_taxonomy") -> dict[str, list[tuple[str, str, str]]]:
    """taxid → [(taxid, name, rank), …] root-first, ending with the taxon itself.

    Each efetch batch is archived as `<name>_<i>.xml`; a rerun reads the archive.
    """
    out: dict[str, list[tuple[str, str, str]]] = {}
    ids = sorted(set(taxids), key=int)
    for i in range(0, len(ids), 150):
        path = raw("taxonomy", f"{name}_{i // 150:03d}.xml")
        chunk = ids[i:i + 150]
        if not path.exists() or not _covers(path, chunk):
            for attempt in range(5):
                r = requests.post(EFETCH, data={"db": "taxonomy", "id": ",".join(chunk),
                                                "retmode": "xml", "tool": "ion_channel_census"},
                                  timeout=120)
                if r.status_code == 200 and "<TaxaSet" in r.text:
                    break
                time.sleep(2 ** attempt)
            r.raise_for_status()
            path.write_text(r.text)
            time.sleep(0.4)
        out.update(parse_taxa(path.read_text()))
    missing = [t for t in ids if t not in out]
    if missing:
        raise RuntimeError(f"NCBI Taxonomy returned no record for {missing}")
    return out


def _covers(path: Path, chunk: list[str]) -> bool:
    return set(chunk) <= set(parse_taxa(path.read_text()))


def parse_taxa(xml: str) -> dict[str, list[tuple[str, str, str]]]:
    out = {}
    for tx in ET.fromstring(xml).findall("Taxon"):
        tid = tx.findtext("TaxId")
        lin = [(a.findtext("TaxId"), a.findtext("ScientificName"), a.findtext("Rank"))
               for a in tx.findall("LineageEx/Taxon")]
        lin.append((tid, tx.findtext("ScientificName"), tx.findtext("Rank")))
        out[tid] = lin
        for a in tx.findall("AkaTaxIds/TaxId"):   # merged ids answer for the new one
            out[a.text] = lin
    return out


# ---------------------------------------------------------------- tree
@dataclass
class TNode:
    name: str
    rank: str = ""
    taxid: str = ""
    children: list["TNode"] = field(default_factory=list)
    tip: str | None = None                     # the order a tip stands for

    def tips(self) -> list[str]:
        return [self.tip] if self.tip else [t for c in self.children for t in c.tips()]

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()


def build_tree(tips: dict[str, list[tuple[str, str, str]]], root: str = ROOT_TAXON) -> TNode:
    """tips: order name → its lineage (root-first, ending at the order).

    Nodes are taxa; the tree starts at `root`; a node with one child is
    collapsed into the child (the deeper name is kept), polytomies stay.
    """
    nodes: dict[str, TNode] = {}
    top = None
    for order, lin in sorted(tips.items()):
        names = [x[1] for x in lin]
        if root not in names:
            raise ValueError(f"{order}: {root} not on its lineage")
        lin = lin[names.index(root):]
        parent = None
        for tid, nm, rk in lin:
            n = nodes.get(tid)
            if n is None:
                n = nodes[tid] = TNode(nm, rk, tid)
                if parent is not None:
                    parent.children.append(n)
            parent = n
        if parent.tip and parent.tip != order:
            raise ValueError(f"two tips on one taxon: {parent.tip}, {order}")
        if parent.children:
            # an order whose lineage is an ancestor of another tip: give it its own leaf
            leaf = TNode(order, "order", parent.taxid, tip=order)
            parent.children.append(leaf)
        else:
            parent.tip = order
        top = top or nodes[lin[0][0]]
    return collapse(top)


def collapse(n: TNode) -> TNode:
    n.children = [collapse(c) for c in n.children]
    while len(n.children) == 1 and not n.tip:
        n = n.children[0]
    return n


def newick(n: TNode) -> str:
    def q(s: str) -> str:
        return "'" + s.replace("'", "") + "'"
    def w(x: TNode) -> str:
        if x.tip:
            return q(x.tip)
        return "(" + ",".join(w(c) for c in x.children) + ")" + q(x.name)
    return w(n) + ";"


# ---------------------------------------------------------------- Sankoff
def _cost(s: int, t: int, g: int) -> int:
    return 0 if s == t else (g if t == 1 else 1)


def sankoff(tree: TNode, states: dict[str, int | None], g: int) -> dict:
    """Marginal Sankoff. Returns {'score', 'node': {id: state|None},
    'edges': [(parent, child, pair|None, changes:set)], 'root': state|None}."""
    down: dict[int, list[int]] = {}
    post = list(tree.walk())[::-1]
    for n in post:
        if n.tip is not None and not n.children:
            st = states.get(n.tip)
            down[id(n)] = [0, 0] if st is None else ([0, INF] if st == 0 else [INF, 0])
        else:
            down[id(n)] = [sum(min(down[id(c)][t] + _cost(s, t, g) for t in (0, 1))
                               for c in n.children) for s in (0, 1)]
    up: dict[int, list[int]] = {id(tree): [0, g]}
    edges = []
    for p in tree.walk():
        for c in p.children:
            best_c = [min(down[id(c)][t] + _cost(s, t, g) for t in (0, 1)) for s in (0, 1)]
            joint = {(s, t): up[id(p)][s] + down[id(p)][s] - best_c[s] + _cost(s, t, g) + down[id(c)][t]
                     for s in (0, 1) for t in (0, 1)}
            up[id(c)] = [min(joint[(s, t)] - down[id(c)][t] for s in (0, 1)) for t in (0, 1)]
            edges.append((p, c, joint))
    total = [down[id(tree)][s] + up[id(tree)][s] for s in (0, 1)]
    score = min(total)
    node = {}
    for n in tree.walk():
        m = [down[id(n)][s] + up[id(n)][s] for s in (0, 1)]
        opt = [s for s in (0, 1) if m[s] == score]
        node[id(n)] = opt[0] if len(opt) == 1 else None
    out_edges = []
    for p, c, joint in edges:
        opt = [k for k, v in joint.items() if v == score]
        out_edges.append((p, c, opt[0] if len(opt) == 1 else None,
                          {("gain" if t else "loss") for s, t in opt if s != t}))
    rs = [s for s in (0, 1) if total[s] == score]
    return {"score": score, "node": node, "edges": out_edges,
            "root": rs[0] if len(rs) == 1 else None}


def events(res: dict) -> list[dict]:
    """Every edge where a change is possible: stated gain/loss or ambiguous."""
    out = []
    for p, c, pair, changes in res["edges"]:
        if pair is not None and pair[0] != pair[1]:
            out.append({"parent": p.name, "child": c.name, "event": "gain" if pair[1] else "loss",
                        "stated": 1, "node": c})
        elif pair is None and changes:
            out.append({"parent": p.name, "child": c.name, "event": "/".join(sorted(changes)),
                        "stated": 0, "node": c})
    return out
