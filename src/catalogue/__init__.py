"""The ion-channel catalogue — this project's definition of its own subject.

Contents:
    schema.py         — ChannelFamily / Superfamily / Signature / Hazard types
    vgic_k.py         — P-loop superfamily, potassium branch (Kv, Kir, K2P, Slo, SK)
    vgic_cation.py    — P-loop superfamily, Nav / Cav / TPC / CNG / HCN / TRP
    lgic.py           — Cys-loop, iGluR, P2X and DEG/ENaC superfamilies
    anion.py          — CLC, bestrophin, CFTR, tweety
    tmem16_like.py    — anoctamins, OSCA/TMEM63, TMC (one fold, no sequence homology)
    mechano.py        — Piezo, MscL, MscS
    largepore.py      — connexins, pannexins/innexins, LRRC8, CALHM
    intracellular.py  — ITPR, RYR, TRIC, MCU, VDAC, TMEM175
    other.py          — ORAI, Hv1, otopetrins, CLIC, viroporins
    proposed.py       — proposed channels added after S20 (PACC1 …) + their look-alikes
    controls.py       — auxiliary subunits, non-channel homologues, out-of-scope
    hazards.py        — the twenty recorded ways to get a classification wrong
    registry.py       — assembly, lookups, `validate()`, `stats()`

Import from `registry` (or from this package), never from a division file.

    from src.catalogue import CATALOGUE, census_families, validate
"""

from .hazards import HAZARD_BY_ID, HAZARDS
from .registry import (CATALOGUE, SUPERFAMILIES, census_families,
                       control_families, exemplars, families, family,
                       family_for_gene, hazards_for, human_genes,
                       pore_signatures, reference_panel, shared_signatures,
                       signature_index, stats, superfamily, validate,
                       with_resolved_exemplars)
from .schema import (CENSUS_STATUSES, ChannelFamily, Exemplar, Fold, Gating,
                     Hazard, Level, Provenance, Selectivity, Signature, Status,
                     Superfamily)

__all__ = [
    "CATALOGUE", "SUPERFAMILIES", "HAZARDS", "HAZARD_BY_ID", "CENSUS_STATUSES",
    "ChannelFamily", "Superfamily", "Signature", "Exemplar", "Hazard",
    "Level", "Provenance", "Gating", "Selectivity", "Fold", "Status",
    "family", "superfamily", "families", "census_families", "control_families",
    "human_genes", "family_for_gene", "signature_index", "shared_signatures",
    "pore_signatures", "exemplars", "reference_panel", "hazards_for",
    "with_resolved_exemplars", "validate", "stats",
]
