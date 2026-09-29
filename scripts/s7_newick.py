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
