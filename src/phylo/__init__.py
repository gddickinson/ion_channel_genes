"""Hierarchical phylogeny: a forest of trees plus a fold network.

Contents:
    modules.py  — pore-module extraction (envelope-based or reference-projected)
    forest.py   — tier-1 (within family) and tier-2 (within superfamily,
                  pore-module) trees; refuses non-alignable superfamilies
    network.py  — tier-3 fold-similarity network, which is explicitly not a tree

    from src.phylo import build_tier1, build_tier2, NotAlignable, seed_network
"""

from .forest import (NotAlignable, TreeRun, alignable_superfamilies,
                     build_tier1, build_tier2, build_tree,
                     needs_modules, refused_superfamilies, toolchain,
                     write_run)
from .modules import (ModuleSet, PoreModule, module_by_projection,
                      module_from_envelope, modules_from_envelopes,
                      whole_sequence_module)
from .network import (LITERATURE_EDGES, FoldEdge, FoldNetwork, seed_network)

__all__ = [
    "build_tree", "build_tier1", "build_tier2", "NotAlignable", "TreeRun",
    "toolchain", "write_run", "alignable_superfamilies",
    "refused_superfamilies", "needs_modules",
    "PoreModule", "ModuleSet", "module_from_envelope",
    "modules_from_envelopes", "module_by_projection", "whole_sequence_module",
    "FoldNetwork", "FoldEdge", "seed_network", "LITERATURE_EDGES",
]
