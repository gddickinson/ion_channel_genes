"""Large-pore channels: connexins, pannexins/innexins, CALHM and LRRC8.

These conduct ions, but a pore wide enough for ATP or glutathione is a
different object from a selectivity filter that discriminates K+ from Na+ by
a fraction of an ångström. They are in the census because they are channels
by every operational definition — a gated aqueous pathway across the
membrane — and flagged `LARGE_PORE` so that any statement about ion
selectivity can exclude them cleanly.

The interesting structure here is a convergence with a twist. Connexins
(vertebrate gap junctions), innexins (invertebrate gap junctions),
pannexins (vertebrate, non-junctional) and LRRC8 (the volume-regulated
anion channel) build the same kind of hexa/octameric large pore, and
*pannexin, innexin and LRRC8 are genuinely homologous* — measured: PANX1 and
PANX2 carry `PF00876` (Innexin), and LRRC8A/B carry `PF12534`, whose Pfam
name is literally "Pannexin-like TM region of LRRC8". Connexins (`PF00029`)
are not detectably related to any of them: vertebrates run two independent
gap-junction systems, one of which they share with invertebrates and one of
which they invented. CALHM is a fourth, unrelated solution.

That makes this the catalogue's best-supported case of convergent evolution
at the fold level, and it is stated as a hypothesis with the evidence
attached rather than as a fact, because "not detectably related" is a
statement about detection methods (**Q6**).
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

INNEXIN = Sig("pfam", "PF00876", "Innexin", L.SUPERFAMILY, 1, P.DB,
              "measured on both invertebrate innexins and vertebrate pannexins")

SUPERFAMILIES: list[SF] = [
    SF("connexin", "Connexins (vertebrate gap junctions)", Fold.CONNEXIN, True,
       "full-length", (Sig("pfam", "PF00029", "Connexin", L.SUPERFAMILY, 1, P.DB),),
       notes="No detectable relationship to innexin/pannexin/LRRC8."),
    SF("innexin_like", "Innexin / pannexin / LRRC8 clan", Fold.INNEXIN, True,
       "TM region", (INNEXIN,),
       notes=("Alignable across pannexin and innexin; LRRC8's membrane region "
              "is homologous (`PF12534`) but its LRR domain is not, so "
              "cross-family alignment is restricted to the TM region."),
       module_rule="tm_span"),
    SF("calhm", "Calcium homeostasis modulators", Fold.CALHM, True, "full-length",
       (Sig("pfam", "PF14798", "Ca_hom_mod", L.SUPERFAMILY, 1, P.DB),)),
]

FAMILIES: list[CF] = [
    CF(
        key="connexin",
        name="Connexins (gap-junction hemichannels)",
        superfamily="connexin", status=St.CHANNEL, fold=Fold.CONNEXIN,
        selectivity=Sel.LARGE_PORE, gating=(G.VOLTAGE, G.CALCIUM, G.PROTON),
        stoichiometry="hexamer (hemichannel); dodecamer across a junction",
        tm_per_subunit=4,
        signatures=(Sig("pfam", "PF00029", "Connexin", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF03508", "Connexin43", L.SUBFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_GJA1", "GJA1", "Homo sapiens", "P17302", "Cx43"),
                   Ex("Hs_GJB2", "GJB2", "Homo sapiens", "P29033", "Cx26")),
        human_genes=("GJA1", "GJA3", "GJA4", "GJA5", "GJA8", "GJA9", "GJA10",
                     "GJB1", "GJB2", "GJB3", "GJB4", "GJB5", "GJB6", "GJB7",
                     "GJC1", "GJC2", "GJC3", "GJD2", "GJD3", "GJD4", "GJE1"),
        length_band_aa=(200, 550),
        notes=("21 human genes named twice over — by Greek-letter subfamily "
               "(GJA1) and by predicted molecular weight (Cx43) — with the two "
               "systems not in register. Pure nomenclature hazard; GJB2 "
               "mutations are the commonest cause of inherited deafness."),
    ),
    CF(
        key="pannexin",
        name="Pannexins (PANX1–PANX3)",
        superfamily="innexin_like", status=St.CHANNEL, fold=Fold.INNEXIN,
        selectivity=Sel.LARGE_PORE, gating=(G.VOLTAGE, G.MECHANICAL, G.CALCIUM),
        stoichiometry="heptamer", tm_per_subunit=4,
        signatures=(INNEXIN,),
        exemplars=(Ex("Hs_PANX1", "PANX1", "Homo sapiens", "Q96RD7"),
                   Ex("Hs_PANX2", "PANX2", "Homo sapiens", "Q96RD6")),
        human_genes=("PANX1", "PANX2", "PANX3"),
        length_band_aa=(390, 680),
        notes=("Glycosylated, so they do not form junctions: a released-ATP "
               "conduit rather than a cell-to-cell channel. Their `PF00876` "
               "annotation is what proves the innexin relationship."),
    ),
    CF(
        key="innexin",
        name="Innexins (invertebrate gap junctions)",
        superfamily="innexin_like", status=St.CHANNEL, fold=Fold.INNEXIN,
        selectivity=Sel.LARGE_PORE, gating=(G.VOLTAGE, G.CALCIUM),
        stoichiometry="octamer", tm_per_subunit=4,
        signatures=(INNEXIN,),
        exemplars=(Ex("Ce_unc-7", "unc-7", "Caenorhabditis elegans"),
                   Ex("Dm_ogre", "ogre", "Drosophila melanogaster")),
        other_genes=("unc-7", "unc-9", "inx-6", "ogre", "shakB"),
        length_band_aa=(350, 550),
        notes=("~25 in *C. elegans*, ~8 in *Drosophila*. Invertebrates have no "
               "connexins at all, so the two gap-junction systems are a clean "
               "lineage split."),
    ),
    CF(
        key="lrrc8",
        name="LRRC8 volume-regulated anion channels (VRAC)",
        superfamily="innexin_like", status=St.CHANNEL, fold=Fold.LRRC8,
        selectivity=Sel.ANION, gating=(G.VOLUME,),
        stoichiometry="hexamer, obligate heteromer with LRRC8A", tm_per_subunit=4,
        signatures=(Sig("pfam", "PF12534", "Pannexin_like", L.FAMILY, 1, P.DB,
                        "Pfam's own name for it is 'Pannexin-like TM region of "
                        "LRRC8' — the homology is in the annotation"),
                    Sig("pfam", "PF13855", "LRR_8", L.ACCESSORY, 1, P.DB)),
        exemplars=(Ex("Hs_LRRC8A", "LRRC8A", "Homo sapiens", "Q8IWT6"),
                   Ex("Hs_LRRC8B", "LRRC8B", "Homo sapiens", "Q6P9F7")),
        human_genes=("LRRC8A", "LRRC8B", "LRRC8C", "LRRC8D", "LRRC8E"),
        length_band_aa=(700, 860),
        notes=("The molecular identity of VRAC, unresolved for thirty years and "
               "settled in 2014 by two parallel genome-wide screens. A pannexin "
               "fold with a leucine-rich-repeat domain bolted on; substrate "
               "selectivity (including cisplatin uptake) is set by which "
               "LRRC8B–E subunits join LRRC8A."),
    ),
    CF(
        key="calhm",
        name="Calcium homeostasis modulators (CALHM1–CALHM6)",
        superfamily="calhm", status=St.CHANNEL, fold=Fold.CALHM,
        selectivity=Sel.LARGE_PORE, gating=(G.VOLTAGE, G.CALCIUM),
        stoichiometry="octamer to undecamer", tm_per_subunit=4,
        signatures=(Sig("pfam", "PF14798", "Ca_hom_mod", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_CALHM1", "CALHM1", "Homo sapiens", "Q8IU99"),
                   Ex("Hs_CALHM2", "CALHM2", "Homo sapiens", "Q9HA72")),
        human_genes=("CALHM1", "CALHM2", "CALHM3", "CALHM4", "CALHM5", "CALHM6"),
        length_band_aa=(300, 400),
        notes=("The ATP-release channel of sweet, bitter and umami taste cells "
               "(CALHM1/CALHM3). Oligomeric state varies by paralogue — CALHM2 "
               "is an undecamer — which makes stoichiometry a per-protein "
               "property rather than a family constant."),
    ),
]
