"""Tier 3 — the fold-similarity network that replaces the tree nobody can build.

Between superfamilies there is no alignment and therefore no phylogeny. What
there *is* is measurable structural similarity: TMEM16, OSCA and TMC share a
ten-transmembrane fold with a conduction groove in the same place; the iGluR
pore is an inverted Kir pore; pannexin, innexin and LRRC8 build the same
large pore. Those relationships are real, they are probably ancestral, and
the evidence for them is structural rather than sequential.

So tier 3 is a **network**, and it is drawn and named as one: nodes are
families, edges are structural similarity (Foldseek TM-score or E-value
between representative predicted structures), and no edge is ever given a
branch length or a bootstrap value. The output is explicitly not a
phylogram, and the module has no code that could produce one.

Where structures are missing the edge is absent and recorded as absent — a
missing edge in this network means "not measured", never "not related", and
the manifest keeps the two apart.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from ..catalogue import CATALOGUE, SUPERFAMILIES


@dataclass
class FoldEdge:
    source: str            # family or superfamily key
    target: str
    metric: str            # "tm_score" | "foldseek_evalue" | "literature"
    value: float
    evidence: str = ""     # structure ids, tool version, or citation
    measured: bool = True

    def key(self) -> tuple[str, str]:
        return tuple(sorted((self.source, self.target)))


@dataclass
class FoldNetwork:
    level: str = "superfamily"     # or "family"
    nodes: list[str] = field(default_factory=list)
    edges: list[FoldEdge] = field(default_factory=list)
    unmeasured: list[tuple[str, str, str]] = field(default_factory=list)

    def add(self, edge: FoldEdge) -> None:
        self.edges.append(edge)

    def note_unmeasured(self, a: str, b: str, why: str) -> None:
        self.unmeasured.append((a, b, why))

    def neighbours(self, node: str) -> list[FoldEdge]:
        return [e for e in self.edges if node in (e.source, e.target)]

    def summary(self) -> str:
        return (f"fold network at {self.level} level: {len(self.nodes)} nodes, "
                f"{len(self.edges)} measured edge(s), "
                f"{len(self.unmeasured)} pair(s) not measured")

    def to_dict(self) -> dict:
        return {
            "level": self.level,
            "nodes": self.nodes,
            "edges": [e.__dict__ for e in self.edges],
            "unmeasured": self.unmeasured,
            "caveat": ("Edges are structural similarity, not phylogeny. No "
                       "branch lengths, no support values, no common "
                       "ancestor is implied by an edge (D27)."),
        }

    def write(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2))
        return path

    def edge_table(self) -> list[tuple]:
        return [(e.source, e.target, e.metric, f"{e.value:.4g}",
                 "measured" if e.measured else "asserted", e.evidence)
                for e in sorted(self.edges, key=lambda e: -e.value)]


def seed_network(level: str = "superfamily") -> FoldNetwork:
    """A network with every node and no edges — the honest starting state.

    The literature-asserted relationships are added as `measured=False`
    edges by `scripts/s8_fold_network.py` only after it has said so; nothing
    is asserted here, because a network that ships with its conclusions
    already drawn is a diagram, not a measurement.
    """
    if level == "superfamily":
        nodes = sorted(SUPERFAMILIES)
    else:
        nodes = sorted(f.key for f in CATALOGUE.values() if f.census_member())
    return FoldNetwork(level=level, nodes=nodes)


#: Relationships the literature asserts and this project must test rather
#: than inherit. Each becomes a `measured=False` edge until S8 measures it.
LITERATURE_EDGES: list[tuple[str, str, str]] = [
    ("tmem16_like", "tmem16_like",
     "TMEM16, OSCA/TMEM63 and TMC share a 10-TM fold with the conduction "
     "groove between TM4 and TM6 (structural studies)"),
    ("iglur", "ploop",
     "the iGluR pore module is an inverted P-loop, structurally equivalent "
     "to Kir"),
    ("innexin_like", "connexin",
     "innexin/pannexin/LRRC8 and connexin build convergent large pores with "
     "no detectable sequence relationship"),
    ("ca_release", "ploop",
     "the ITPR/RYR pore module is a P-loop channel carrying PF00520"),
    ("hv", "ploop",
     "Hv1's voltage-sensor domain is homologous to the VSD of the P-loop "
     "channels, without the pore domain"),
]
