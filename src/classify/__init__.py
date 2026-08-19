"""Classification of a protein into the ion-channel catalogue.

Contents:
    rules.py       — the architecture tier: domain-composition rules derived
                     from the catalogue, plus the hand-written hazard rules
    motifs.py      — the motif tier: the K+ filter regex and the four-repeat
                     selectivity-filter projection (Nav/Cav/NALCN)
    reference.py   — the reference tier: best-hit identity to the exemplar
                     panel with a margin (D7)
    classifier.py  — combines the tiers into one `ChannelCall` with an audit
                     trail and a confidence tier
    report.py      — benchmark tables: confusion matrix, per-family recall,
                     per-hazard precision

    from src.classify import ChannelQuery, classify, ReferenceSet
"""

from .classifier import (ChannelCall, ChannelQuery, Evidence, classify,
                         classify_all)
from .motifs import ANCHORS, FOUR_REPEAT_ANCHOR, MotifResult, k_filter_hits
from .reference import ReferenceResult, ReferenceSet, fetch_uniprot_sequence
from .rules import (ArchitectureRule, HAZARD_RULES, RuleMatch, all_rules,
                    derived_rules, match_architecture, topology_check)

__all__ = [
    "ChannelQuery", "ChannelCall", "Evidence", "classify", "classify_all",
    "ArchitectureRule", "HAZARD_RULES", "RuleMatch", "all_rules",
    "derived_rules", "match_architecture", "topology_check",
    "MotifResult", "ANCHORS", "FOUR_REPEAT_ANCHOR", "k_filter_hits",
    "ReferenceSet", "ReferenceResult", "fetch_uniprot_sequence",
]
