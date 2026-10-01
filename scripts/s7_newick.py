"""s7_newick.py — a stdlib Newick reader and the split tests S7 reads trees with.

IQ-TREE writes UFBoot support as the internal-node label. Every question S7
asks of a tree is a question about a **split** (a bipartition of the leaves
made by one edge), which does not depend on where the tree is drawn rooted:

* is the outgroup set separated from the family by one edge? (rootable)
* what UFBoot support does that edge carry? (root support)
* what share of the ingroup's internal edges carry ≥ 95 support?
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Node:
    name: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)

    def leaves(self) -> list[str]:
        if not self.children:
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def parse(text: str) -> Node:
    """Parse one Newick tree (labels without quotes or whitespace)."""
    s = text.strip().rstrip(";")
    pos = 0

    def label() -> tuple[str, float | None]:
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in ",():":
            pos += 1
        name = s[start:pos]
        length = None
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            length = float(s[start:pos])
        return name, length

    def node() -> Node:
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            n.children.append(node())
            while s[pos] == ",":
                pos += 1
                n.children.append(node())
            if s[pos] != ")":
                raise ValueError(f"expected ')' at {pos}")
            pos += 1
        n.name, n.length = label()
        return n

    root = node()
    if pos != len(s):
        raise ValueError(f"trailing text at {pos}")
    return root


def support(name: str) -> float | None:
    """UFBoot from an IQ-TREE internal label ('98', or 'SH/UFBoot' '85.1/98')."""
    if not name:
        return None
    try:
        return float(name.split("/")[-1])
    except ValueError:
        return None


def splits(root: Node) -> list[tuple[frozenset[str], float | None]]:
    """Every non-trivial split below the drawn root, with its support."""
    out = []

    def walk(n: Node) -> frozenset[str]:
        if not n.children:
            return frozenset([n.name])
        below = frozenset().union(*(walk(c) for c in n.children))
        if n is not root:
            out.append((below, support(n.name)))
        return below

    walk(root)
    return out


def split_support(root: Node, group: set[str]) -> tuple[bool, float | None]:
    """Is `group` one side of some edge? Returns (present, UFBoot)."""
    all_leaves = frozenset(root.leaves())
    g = frozenset(group)
    comp = all_leaves - g
    if len(g) == 1 or len(comp) == 1:
        return True, None          # a terminal edge: always present, no support
    for side, sup in splits(root):
        if side == g or side == comp:
            return True, sup
    return False, None


def ingroup_support(root: Node, outgroup: set[str]) -> dict:
    """Support summary over internal edges that split the ingroup only."""
    all_leaves = frozenset(root.leaves())
    ingroup = all_leaves - frozenset(outgroup)
    vals = []
    for side, sup in splits(root):
        other = all_leaves - side
        inner = side if side <= ingroup else other if other <= ingroup else None
        if inner is None or len(inner) < 2 or len(ingroup - inner) < 1:
            continue
        if sup is not None:
            vals.append(sup)
    vals.sort()
    n = len(vals)
    return {"internal_edges": n,
            "ufboot_median": vals[n // 2] if n else None,
            "frac_ge95": round(sum(v >= 95 for v in vals) / n, 4) if n else None,
            "frac_lt70": round(sum(v < 70 for v in vals) / n, 4) if n else None}


def root_partition(root: Node, outgroup: set[str]) -> list[tuple[frozenset[str], float | None]] | None:
    """The ingroup's basal split when the tree is rooted on `outgroup` (S7c).

    The outgroup must be one side of an edge; the ingroup-side node of that
    edge is the ingroup's root, and its other neighbours' leaf sets are the
    root's clades (two for a bifurcating root), each with the UFBoot of its
    own edge. None when the outgroup is not a clade — the root is undefined.
    """
    adj: dict[int, list[tuple[Node, float | None]]] = {}
    nodes: dict[int, Node] = {}

    def link(a: Node, b: Node, sup: float | None) -> None:
        adj.setdefault(id(a), []).append((b, sup))
        adj.setdefault(id(b), []).append((a, sup))
        nodes[id(a)], nodes[id(b)] = a, b

    def walk(n: Node) -> None:
        for c in n.children:
            link(n, c, support(c.name) if c.children else None)
            walk(c)
    walk(root)

    def leaves_away(start: Node, frm: Node) -> frozenset[str]:
        out, stack = set(), [(start, frm)]
        while stack:
            n, p = stack.pop()
            nbrs = [m for m, _ in adj.get(id(n), []) if m is not p]
            if not nbrs:
                out.add(n.name)
            stack.extend((m, n) for m in nbrs)
        return frozenset(out)

    og = frozenset(outgroup)
    for a_id, nbrs in adj.items():
        a = nodes[a_id]
        for b, _ in nbrs:
            if leaves_away(b, a) == og:     # edge a–b, outgroup beyond b
                return sorted(((leaves_away(m, a), s) for m, s in adj[a_id]
                               if m is not b), key=lambda x: (len(x[0]), sorted(x[0])))
    return None
