"""s13_tree.py — the unrooted-tree operations S13's sets and labels need (D55).

S7b's trees are unrooted (IQ-TREE draws a trifurcating root), and every set
S13 builds is a **side of an edge**, which needs no root:

* `prune()` — restrict to a leaf set, suppress the degree-2 nodes left
  behind (their two edge lengths summed, the **larger** UFBoot kept: every
  edge on the merged path induces the same split of the kept leaves, so a
  bootstrap tree holding any one of them holds that split — the larger
  value is a lower bound on its support, the smaller is not);
* `sides()` — every edge's two leaf sides with the edge's UFBoot;
* `orthogroup()` — D55 (2): of the sides holding `h` and no other anchor
  tip, the largest whose bounding edge has UFBoot ≥ 95. The sides that hold
  `h` and exclude another fixed tip are nested, so "the largest" is unique;
* `drawn_for()` — the tree redrawn from a node outside every given set, so
  that each set is a clade and PAML's node labels mean what they say;
* `paml_newick()` — node labels `#k` on the branches *within* each set
  (its stem excluded) and `#s` on a named set's stem, HyPhy `{Name}` the
  same way. A set that is not a clade in the drawn tree raises — codeml
  marks a node, so a non-clade foreground would silently mark a larger one.
"""

from __future__ import annotations

from collections import defaultdict

from s7_newick import Node, parse, support

SUPPORT_MIN = 95.0


class UTree:
    """An undirected tree: adjacency with edge lengths and UFBoot."""

    def __init__(self) -> None:
        self.adj: dict[int, dict[int, tuple[float, float | None]]] = defaultdict(dict)
        self.name: dict[int, str] = {}

    @classmethod
    def from_newick(cls, text: str) -> "UTree":
        t, counter = cls(), [0]

        def add(n: Node, parent: int | None) -> None:
            nid = counter[0]
            counter[0] += 1
            if not n.children:
                t.name[nid] = n.name
            t.adj[nid]
            if parent is not None:
                sup = support(n.name) if n.children else None
                t.link(parent, nid, n.length or 0.0, sup)
            for c in n.children:
                add(c, nid)

        add(parse(text), None)
        return t

    def link(self, a: int, b: int, length: float, sup: float | None) -> None:
        self.adj[a][b] = (length, sup)
        self.adj[b][a] = (length, sup)

    def leaves(self) -> set[str]:
        return set(self.name.values())

    def leaf_id(self, label: str) -> int:
        return next(i for i, n in self.name.items() if n == label)

    def side(self, a: int, b: int) -> frozenset[str]:
        """Leaves reached from b without crossing back to a."""
        out, stack, seen = set(), [b], {a, b}
        while stack:
            x = stack.pop()
            if x in self.name:
                out.add(self.name[x])
            for y in self.adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        return frozenset(out)


def prune(t: UTree, keep: set[str]) -> UTree:
    """A copy restricted to `keep`, degree-2 nodes suppressed."""
    missing = keep - t.leaves()
    if missing:
        raise ValueError(f"prune: {len(missing)} tips not in tree, e.g. {sorted(missing)[:3]}")
    u = UTree()
    for a, nb in t.adj.items():
        u.adj[a] = dict(nb)
    u.name = {i: n for i, n in t.name.items() if n in keep}
    changed = True
    while changed:
        changed = False
        for x in list(u.adj):
            deg = len(u.adj[x])
            if x in t.name and x not in u.name and deg <= 1:      # dropped leaf
                for y in u.adj.pop(x):
                    u.adj[y].pop(x, None)
                changed = True
            elif x not in t.name and deg == 1:                    # bare internal
                for y in u.adj.pop(x):
                    u.adj[y].pop(x, None)
                changed = True
            elif x not in u.name and deg == 2 and x not in t.name:
                (y1, (l1, s1)), (y2, (l2, s2)) = u.adj[x].items()
                sup = max((s for s in (s1, s2) if s is not None), default=None)
                del u.adj[x]
                u.adj[y1].pop(x)
                u.adj[y2].pop(x)
                u.link(y1, y2, l1 + l2, sup)
                changed = True
    return u


def sides(t: UTree) -> list[tuple[frozenset[str], float | None, tuple[int, int]]]:
    """Both sides of every edge: (leaf side beyond b, UFBoot, (a, b))."""
    out, done = [], set()
    for a, nb in t.adj.items():
        for b, (_, sup) in nb.items():
            if (b, a) in done:
                continue
            done.add((a, b))
            out.append((t.side(a, b), sup, (a, b)))
            out.append((t.side(b, a), sup, (b, a)))
    return out


def orthogroup(t: UTree, h: str, others: set[str]) -> tuple[frozenset[str], float | None]:
    """D55 (2): the largest side holding `h`, none of `others`, UFBoot ≥ 95."""
    if not others:
        return frozenset(t.leaves()), None
    if len(others) == 1:
        # With two anchor tips the complement of the other one's terminal
        # edge (trivially supported) holds h and excludes it, so "the largest
        # side" would be everything but that one tip. Not a defined case.
        raise ValueError("orthogroup rule needs 0 or >= 2 other anchor tips")
    ok = [(s, sup) for s, sup, _ in sides(t)
          if h in s and not (s & others)
          and (sup is None or sup >= SUPPORT_MIN or len(s) == 1)]
    best = max(ok, key=lambda x: len(x[0]))
    return best[0], (best[1] if len(best[0]) > 1 else None)


def bounding_edge(t: UTree, group: frozenset[str]) -> tuple[int, int]:
    """(outer, inner) node pair of the edge whose inner side is `group`."""
    for s, _, (a, b) in sides(t):
        if s == group:
            return a, b
    raise ValueError(f"set of {len(group)} tips is not one side of any edge")


def drawn_for(t: UTree, groups: list[frozenset[str]]) -> Node:
    """Redraw rooted at a node outside every group, so each is a clade."""
    if len(groups) == 1 and groups[0] == frozenset(t.leaves()):
        start = next(x for x in t.adj if x not in t.name)
    else:
        outer, _ = bounding_edge(t, groups[0])
        start = outer
        for g in groups[1:]:
            if start in _inner_nodes(t, g):
                raise ValueError("groups overlap: no node lies outside all of them")

    def build(x: int, frm: int | None) -> Node:
        kids = [build(y, x) for y in t.adj[x] if y != frm]
        length = t.adj[x][frm][0] if frm is not None else None
        return Node(name=t.name.get(x, ""), length=length, children=kids)

    return build(start, None)


def _inner_nodes(t: UTree, group: frozenset[str]) -> set[int]:
    outer, inner = bounding_edge(t, group)
    out, stack, seen = set(), [inner], {outer, inner}
    while stack:
        x = stack.pop()
        out.add(x)
        for y in t.adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return out


def _clade_node(root: Node, group: frozenset[str]) -> Node:
    best = None
    stack = [root]
    while stack:
        n = stack.pop()
        if frozenset(n.leaves()) == group:
            best = n
        stack.extend(n.children)
    if best is None:
        raise ValueError(f"set of {len(group)} tips is not a clade in the drawn tree")
    return best


def paml_newick(root: Node, within: dict[frozenset[str], str],
                stem: dict[frozenset[str], str] | None = None,
                lengths: bool = False, style: str = "paml") -> str:
    """Newick with labels: every branch inside each `within` set (stem
    excluded) gets its mark; each `stem` set's stem branch gets its own.
    style 'paml' writes `#1`; 'hyphy' writes `{Name}`."""
    stem = stem or {}
    mark: dict[int, str] = {}
    for g, tag in within.items():
        top = _clade_node(root, g)
        stack = list(top.children)
        while stack:
            n = stack.pop()
            mark[id(n)] = tag
            stack.extend(n.children)
    for g, tag in stem.items():
        mark[id(_clade_node(root, g))] = tag

    def fmt(tag: str) -> str:
        return (tag if style == "paml" else f"{{{tag}}}") if tag else ""

    def rec(n: Node, top: bool = False) -> str:
        body = n.name if not n.children else "(" + ",".join(rec(c) for c in n.children) + ")"
        out = body + ("" if top else fmt(mark.get(id(n), "")))
        if lengths and n.length is not None and not top:
            out += f":{n.length:.6f}"
        return out

    return rec(root, top=True) + ";"
