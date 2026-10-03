"""s11_recon.py — gene-tree / species-tree reconciliation under D52 (S11a).

Pure logic, no IO. A gene tree is held **unrooted** (`UTree`): every edge
carries its length and UFBoot, and every possible root is an edge. For a
directed edge p → w ("node w hanging from p") the subtree below w does not
depend on where the root is, so its leaf set, species set, species-tree
mapping, duplication flag and (duplications, losses) cost are computed once
and memoised; the cost of rooting on edge (u, v) is then
cost(u → v) + cost(v → u) + the root node's own term — all roots in O(n).

* **Duplication** (D52 (4)): a node whose children's species sets
  intersect — species overlap, independent of the species tree.
* **Losses** (D52 (5)): LCA mapping onto the species tree; one loss per
  species-tree node skipped between a node's mapping and a child's
  (a speciation node is expected to skip none), floored at 0 where a
  polytomy maps parent and child to the same node.
* **Identity**: a non-root node is its directed edge (p, w) — the same leaf
  set and child split under every root on the far side of p. The root node
  is ("root", edge). A duplication is *stated* iff it is one under every
  optimal root.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

from s7_newick import Node, support

sys.setrecursionlimit(20000)
SUPPORTED = 95.0                       # D52 (4): each child a leaf or UFBoot ≥ 95


# ---------------------------------------------------------------- species tree
@dataclass
class STree:
    """Species tree with species bitmasks; `lca(mask)` memoised."""
    names: list[str]                    # node index → taxon name
    tip: list[str | None]               # node index → species (tips only)
    parent: list[int]
    children: list[list[int]]
    depth: list[int]
    mask: list[int]                     # node → species bitmask below it
    species_index: dict[str, int]       # species → bit
    _lca: dict[int, int] = field(default_factory=dict)

    def lca(self, m: int) -> int:
        hit = self._lca.get(m)
        if hit is not None:
            return hit
        n = 0
        while True:
            nxt = [c for c in self.children[n] if self.mask[c] & m == m]
            if not nxt:
                break
            n = nxt[0]
        self._lca[m] = n
        return n


def stree_from(root, tip_attr: str = "tip") -> STree:
    """From an S10 `TNode` tree (name, children, tip)."""
    names, tip, parent, children, depth = [], [], [], [], []

    def walk(n, p: int, d: int) -> int:
        i = len(names)
        names.append(n.name)
        tip.append(getattr(n, tip_attr))
        parent.append(p)
        children.append([])
        depth.append(d)
        for c in n.children:
            children[i].append(walk(c, i, d + 1))
        return i

    walk(root, -1, 0)
    species = sorted(t for t in tip if t)
    idx = {s: k for k, s in enumerate(species)}
    mask = [0] * len(names)
    for i in sorted(range(len(names)), key=lambda k: -depth[k]):
        if tip[i]:
            mask[i] |= 1 << idx[tip[i]]
        if parent[i] >= 0:
            mask[parent[i]] |= mask[i]
    return STree(names, tip, parent, children, depth, mask, idx)


# ---------------------------------------------------------------- gene tree
@dataclass
class UTree:
    leaf: dict[int, str]                          # node id → leaf name
    adj: dict[int, list[int]]
    sup: dict[frozenset, float | None]            # edge → UFBoot (None: terminal)

    def edges(self) -> list[tuple[int, int]]:
        return sorted({tuple(sorted(e)) for e in self.sup})


def utree(root: Node, drop: set[str] = frozenset()) -> UTree:
    """Unrooted tree from a parsed Newick, `drop` leaves pruned.

    A node left with two neighbours is suppressed and its two edges merged;
    the merged edge keeps the **smaller** support (conservative: either
    original edge now induces the same ingroup split).
    """
    leaf, adj, sup = {}, {}, {}
    nid = 0

    def add(n: Node) -> int:
        nonlocal nid
        i = nid
        nid += 1
        adj[i] = []
        if not n.children:
            leaf[i] = n.name
        for c in n.children:
            j = add(c)
            adj[i].append(j)
            adj[j].append(i)
            sup[frozenset((i, j))] = support(c.name) if c.children else None
        return i

    add(root)
    for i in [k for k, v in leaf.items() if v in drop]:
        _remove(i, leaf, adj, sup)
    for i in list(adj):                                # dangling internal nodes
        _prune_dead(i, leaf, adj, sup)
    for i in list(adj):
        if i in adj and i not in leaf and len(adj[i]) == 2:
            _suppress(i, adj, sup)
    return UTree(leaf, adj, sup)


def _remove(i, leaf, adj, sup):
    for j in adj.pop(i):
        adj[j].remove(i)
        sup.pop(frozenset((i, j)))
    leaf.pop(i, None)


def _prune_dead(i, leaf, adj, sup):
    while i in adj and i not in leaf and len(adj[i]) <= 1:
        nb = list(adj[i])
        _remove(i, leaf, adj, sup)
        if not nb:
            return
        i = nb[0]


def _suppress(i, adj, sup):
    a, b = adj.pop(i)
    sa, sb = sup.pop(frozenset((i, a))), sup.pop(frozenset((i, b)))
    vals = [s for s in (sa, sb) if s is not None]
    adj[a].remove(i)
    adj[b].remove(i)
    adj[a].append(b)
    adj[b].append(a)
    sup[frozenset((a, b))] = min(vals) if vals else None


# ---------------------------------------------------------------- reconciliation
@dataclass
class Sub:
    """The subtree of node w hanging from p (directed edge p → w)."""
    leaves: int            # bitmask over gene-tree leaf ids
    species: int           # bitmask over species
    m: int                 # species-tree node it maps to
    dup: bool
    kids: tuple            # child directed edges
    dups: int              # duplications in the subtree (incl. w)
    losses: int


class Reconciler:
    def __init__(self, t: UTree, st: STree, species_of: dict[str, str]):
        self.t, self.st = t, st
        self.lid = {n: k for k, n in enumerate(sorted(t.leaf))}
        self.sp = {n: 1 << st.species_index[species_of[t.leaf[n]]] for n in t.leaf}
        self.memo: dict[tuple[int, int], Sub] = {}

    def sub(self, p: int, w: int) -> Sub:
        key = (p, w)
        hit = self.memo.get(key)
        if hit is not None:
            return hit
        if w in self.t.leaf:
            s = Sub(1 << self.lid[w], self.sp[w], self.st.lca(self.sp[w]), False, (), 0, 0)
        else:
            kids = tuple((w, c) for c in self.t.adj[w] if c != p)
            s = self._join(kids)
        self.memo[key] = s
        return s

    def _join(self, kids: tuple) -> Sub:
        subs = [self.sub(*k) for k in kids]
        leaves = species = 0
        dup = False
        for s in subs:
            dup = dup or bool(species & s.species)
            leaves |= s.leaves
            species |= s.species
        m = self.st.lca(species)
        d = self.st.depth
        losses = sum(max(0, d[s.m] - d[m] - (0 if dup else 1)) for s in subs)
        return Sub(leaves, species, m, dup, kids, sum(s.dups for s in subs) + dup,
                   sum(s.losses for s in subs) + losses)

    def root(self, u: int, v: int) -> Sub:
        return self._join(((u, v), (v, u)))

    def optimal(self) -> tuple[list[tuple[int, int]], tuple[int, int]]:
        """Every edge whose rooting minimises (duplications, losses)."""
        costs = {e: (r.dups, r.losses) for e in self.t.edges() for r in [self.root(*e)]}
        best = min(costs.values())
        return [e for e, c in costs.items() if c == best], best

    def dup_ids(self, e: tuple[int, int]) -> dict:
        """Identity → (kind, directed edge or root edge) of every duplication under root e."""
        out = {}
        r = self.root(*e)
        if r.dup:
            out[("root", e)] = r
        stack = list(r.kids)
        while stack:
            k = stack.pop()
            s = self.sub(*k)
            if s.dup:
                out[k] = s
            stack.extend(s.kids)
        return out

    def child_support(self, ident, s: Sub) -> list[float | None]:
        """UFBoot of each child edge; None for a leaf child."""
        out = []
        for k in s.kids:
            sup = self.t.sup[frozenset(k)]
            out.append(None if k[1] in self.t.leaf else sup)
        return out

    def supported(self, ident, s: Sub) -> bool:
        return all(x is None or x >= SUPPORTED for x in self.child_support(ident, s))

    def score(self, s: Sub) -> float:
        a, b = (self.sub(*k).species for k in s.kids[:2])
        return round(bin(a & b).count("1") / bin(a | b).count("1"), 4)

    def names(self, mask: int) -> list[str]:
        inv = {k: self.t.leaf[n] for n, k in self.lid.items()}
        return [inv[k] for k in range(len(inv)) if mask >> k & 1]


def edge_for_split(rc: Reconciler, side: set[str]) -> tuple[int, int] | None:
    """The edge whose one side is exactly `side` (leaf names), if any."""
    target = 0
    name_to = {v: k for k, v in rc.t.leaf.items()}
    for n in side:
        target |= 1 << rc.lid[name_to[n]]
    for u, v in rc.t.edges():
        if rc.sub(u, v).leaves == target or rc.sub(v, u).leaves == target:
            return (u, v)
    return None
