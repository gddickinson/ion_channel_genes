"""Mechanosensitive channels that are not TMEM16-like: Piezo and the
prokaryotic MscL/MscS families.

Three unrelated solutions to the same problem — turning membrane tension
into an open pore — and none of them shares a residue with the others.
Piezo is a 2,500-residue three-bladed propeller of 38 transmembrane helices
per subunit that curves the membrane into a dome; MscL is a 136-residue
pentamer that acts as a tension-gated safety valve; MscS is a heptamer with
a lipid-facing gate. Force-from-lipid gating is a mechanism, not a clade,
and the catalogue keeps them apart so that a "mechanosensitive channel"
count is never mistaken for a family.

Piezo is also the reason the length band matters: at 2,521 aa (PIEZO1,
measured) it sits in the same band as the ryanodine receptors and above
every other plasma-membrane channel, so a size-only filter designed for one
family excludes the other.
"""

from __future__ import annotations

from .schema import ChannelFamily as CF
from .schema import Exemplar as Ex
from .schema import Fold
from .schema import Gating as G
from .schema import Level as L
from .schema import Provenance as P
from .schema import Selectivity as Sel
from .schema import Signature as Sig
from .schema import Status as St
from .schema import Superfamily as SF

SUPERFAMILIES: list[SF] = [
    SF("piezo", "Piezo mechanosensitive channels", Fold.PIEZO, True,
       "full-length (only two paralogues)",
       (Sig("pfam", "PF12166", "Piezo_cap", L.SUPERFAMILY, 1, P.DB),)),
    SF("msc", "Prokaryotic mechanosensitive channels (MscL, MscS)", Fold.MSC,
       False, "no cross-family alignment; MscL and MscS are unrelated",
       (), (),
       notes=("Two families under one heading for convenience, explicitly not "
              "one superfamily: `alignable=False` stops the driver from "
              "building a tree across them.")),
]

FAMILIES: list[CF] = [
    CF(
        key="piezo",
        name="Piezo channels (PIEZO1, PIEZO2)",
        superfamily="piezo", status=St.CHANNEL, fold=Fold.PIEZO,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL,),
        stoichiometry="trimer", tm_per_subunit=38,
        signatures=(Sig("pfam", "PF12166", "Piezo_cap", L.FAMILY, 1, P.DB,
                        "the extracellular cap over the pore"),
                    Sig("pfam", "PF15917", "Piezo_TM25-28", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF24871", "Piezo_TM1-24", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF23188", "THU_Piezo1", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF24874", "Piezo_THU9_anchor", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_PIEZO1", "PIEZO1", "Homo sapiens", "Q92508"),
                   Ex("Hs_PIEZO2", "PIEZO2", "Homo sapiens")),
        human_genes=("PIEZO1", "PIEZO2"),
        length_band_aa=(2500, 2850),
        notes=("38 TM helices per subunit — the largest TM count of any known "
               "channel. PIEZO1 loss-of-function causes lymphatic dysplasia and "
               "gain-of-function causes hereditary xerocytosis; PIEZO2 carries "
               "touch and proprioception. The parent project of this codebase "
               "(`../piezo_genes`) is a full census of exactly this family, so "
               "its results are the natural external check on ours."),
    ),
    CF(
        key="mscl",
        name="Large-conductance mechanosensitive channel (MscL)",
        superfamily="msc", status=St.CHANNEL, fold=Fold.MSC,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL,),
        stoichiometry="pentamer", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF01741", "MscL", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Ec_MscL", "mscL", "Escherichia coli", "P0A742"),),
        other_genes=("mscL",),
        length_band_aa=(120, 160),
        notes=("136 aa (measured) — the smallest channel in the catalogue, and "
               "the emergency valve that saves a bacterium from osmotic "
               "downshock. Its 3 nS conductance is the largest."),
    ),
    CF(
        key="mscs",
        name="Small-conductance mechanosensitive channels (MscS, plant MSL)",
        superfamily="msc", status=St.CHANNEL, fold=Fold.MSC,
        selectivity=Sel.ANION, gating=(G.MECHANICAL, G.VOLTAGE),
        stoichiometry="heptamer", tm_per_subunit=3,
        signatures=(Sig("pfam", "PF00924", "MS_channel_2nd", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF05552", "MS_channel_1st_1", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF21082", "MS_channel_3rd", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Ec_MscS", "mscS", "Escherichia coli", "P0C0S1"),
                   Ex("At_MSL10", "MSL10", "Arabidopsis thaliana")),
        other_genes=("mscS", "mscK", "MSL1", "MSL10"),
        length_band_aa=(250, 750),
        notes=("A superfamily in bacteria, archaea, plants, fungi and protists "
               "and absent from animals — the sharpest presence/absence claim "
               "the census makes, and one the roadmap tests rather than "
               "assumes (Q4)."),
    ),
]
