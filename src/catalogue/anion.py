"""Anion channels that are not TMEM16 and not large-pore: CLC, bestrophin,
CFTR and tweety.

Two of the catalogue's four "is it even a channel?" decisions live here.

**CLC.** One fold, one gene family, two transport mechanisms. CLC-1, CLC-2,
CLC-Ka and CLC-Kb are genuine chloride channels; CLC-3 through CLC-7 are
2Cl⁻/H⁺ antiporters, and so is the bacterial structural prototype ClC-ec1.
The difference is a single conserved glutamate on the intracellular side —
present in the transporters, absent (or non-functional) in the channels.
No sequence-family method separates them, because they *are* one family.
The catalogue therefore splits CLC at the sub-family level with an explicit
mechanism field, and the classification pipeline reports the family call and
the mechanism call separately (**D24**: membership and mechanism are
different questions, and conflating them is how the channel literature
ended up with two incompatible member counts).

**CFTR.** An ABC transporter — nucleotide-binding domains, transmembrane
domains, the whole architecture — that conducts anions down their gradient.
It is the counter-example to every rule that says "ABC fold ⇒ not a
channel", and it is why the classifier's exclusion tests are positive tests
on the pore, never on the presence of a transporter domain.

**Tweety** is the honest open case: TTYH1–3 were proposed as the
volume-regulated anion channel, LRRC8 turned out to be that channel, and
recent structures suggest the tweety proteins shape membranes rather than
conduct. Catalogued as contested, counted with a flag, and named in the
roadmap as a question rather than quietly dropped.
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

CLC_SIG = Sig("pfam", "PF00654", "Voltage_CLC", L.SUPERFAMILY, 1, P.DB,
              "carried by CLC channels and CLC antiporters alike (H5)")
CBS = Sig("pfam", "PF00571", "CBS", L.ACCESSORY, 2, P.DB,
          "eukaryotic CLCs only; a nucleotide-sensing pair, not a pore feature")

SUPERFAMILIES: list[SF] = [
    SF("clc", "CLC chloride channel/transporter superfamily", Fold.CLC, True,
       "full-length (one family, so tier 1 and tier 2 coincide)",
       (CLC_SIG,), ("clc_prokaryotic",),
       notes=("The double-barrelled architecture — two pores, one per subunit "
              "of an antiparallel dimer — was inferred from noise analysis "
              "thirty years before the structure confirmed it.")),
    SF("bestrophin", "Bestrophins", Fold.BESTROPHIN, True, "full-length",
       (Sig("pfam", "PF01062", "Bestrophin", L.SUPERFAMILY, 1, P.CURATED),)),
    SF("abc_channel", "ABC-fold anion channel (CFTR)", Fold.ABC, False,
       "not alignable to any other channel; aligns to ABC transporters",
       (), (),
       notes=("A superfamily of one, whose relatives are all transporters. "
              "Its tier-2 'phylogeny' is a phylogeny of the ABCC subfamily, "
              "which is the correct answer and not a channel tree.")),
    SF("tweety", "Tweety family", Fold.UNKNOWN, True, "full-length",
       (Sig("pfam", "PF04906", "Tweety", L.SUPERFAMILY, 1, P.CURATED),)),
]

FAMILIES: list[CF] = [
    CF(
        key="clc_channel",
        name="CLC chloride channels (CLC-1, CLC-2, CLC-Ka, CLC-Kb)",
        superfamily="clc", status=St.CHANNEL, fold=Fold.CLC,
        selectivity=Sel.ANION, gating=(G.VOLTAGE, G.VOLUME, G.PROTON),
        stoichiometry="homodimer, two independent pores", tm_per_subunit=18,
        signatures=(CLC_SIG, CBS),
        exemplars=(Ex("Hs_CLCN1", "CLCN1", "Homo sapiens", "P35523"),
                   Ex("Hs_CLCN2", "CLCN2", "Homo sapiens")),
        human_genes=("CLCN1", "CLCN2", "CLCNKA", "CLCNKB"),
        length_band_aa=(680, 1000), confusable_with=("clc_transporter", "H5"),
        notes=("CLC-K needs barttin (BSND) to traffic — an accessory subunit "
               "whose loss causes Bartter syndrome type IV, catalogued in "
               "`controls.py` and not counted as a channel."),
    ),
    CF(
        key="clc_transporter",
        name="CLC 2Cl-/H+ antiporters (CLC-3 to CLC-7)",
        superfamily="clc", status=St.TRANSPORTER, fold=Fold.CLC,
        selectivity=Sel.ANION, gating=(G.VOLTAGE,),
        stoichiometry="homodimer", tm_per_subunit=18,
        signatures=(CLC_SIG, CBS),
        exemplars=(Ex("Hs_CLCN7", "CLCN7", "Homo sapiens", "P51798"),
                   Ex("Hs_CLCN5", "CLCN5", "Homo sapiens")),
        human_genes=("CLCN3", "CLCN4", "CLCN5", "CLCN6", "CLCN7"),
        length_band_aa=(760, 900), confusable_with=("clc_channel", "H5"),
        notes=("Excluded from the channel census by mechanism, included in the "
               "CLC family by descent. The proton-glutamate is the "
               "discriminating residue; the classifier reports it as evidence "
               "and the roadmap treats 'how well does the glutamate predict "
               "mechanism across the whole family?' as an open question."),
    ),
    CF(
        key="clc_prokaryotic",
        name="Prokaryotic CLC homologues (ClC-ec1)",
        superfamily="clc", status=St.TRANSPORTER, fold=Fold.CLC,
        selectivity=Sel.ANION, gating=(G.PROTON,),
        stoichiometry="homodimer", tm_per_subunit=18,
        signatures=(CLC_SIG,),
        exemplars=(Ex("Ec_ClC-ec1", "clcA", "Escherichia coli", "P37019",
                      "older literature calls the gene eriC"),),
        other_genes=("clcA", "eriC"),
        length_band_aa=(420, 500),
        notes=("The structural prototype of the whole family is an antiporter, "
               "not a channel — the fold's ancestral function is transport, and "
               "channel behaviour is the derived state. The rooting outgroup "
               "for the CLC tree."),
    ),
    CF(
        key="bestrophin",
        name="Bestrophins (BEST1–BEST4)",
        superfamily="bestrophin", status=St.CHANNEL, fold=Fold.BESTROPHIN,
        selectivity=Sel.ANION, gating=(G.CALCIUM, G.VOLUME),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(Sig("pfam", "PF01062", "Bestrophin", L.FAMILY, 1, P.CURATED),),
        exemplars=(Ex("Hs_BEST1", "BEST1", "Homo sapiens", "O76090"),
                   Ex("Kp_BEST", "best", "Klebsiella pneumoniae", "W9BH30",
                      "the bacterial homologue that gave the first structure")),
        human_genes=("BEST1", "BEST2", "BEST3", "BEST4"),
        length_band_aa=(400, 700),
        notes=("Pentameric with a single central pore — unusual for an anion "
               "channel, and unrelated to every other family here. BEST1 "
               "mutations cause Best vitelliform macular dystrophy."),
    ),
    CF(
        key="cftr",
        name="CFTR (ABCC7) — the ABC-fold anion channel",
        superfamily="abc_channel", status=St.CHANNEL, fold=Fold.ABC,
        selectivity=Sel.ANION, gating=(G.PHOSPHORYLATION, G.LIGAND_INTRACELLULAR),
        stoichiometry="monomer", tm_per_subunit=12,
        signatures=(Sig("pfam", "PF00664", "ABC_membrane", L.SHARED_WITH_DECOY, 2, P.DB,
                        "shared with every ABC transporter — the exclusion trap"),
                    Sig("pfam", "PF00005", "ABC_tran", L.SHARED_WITH_DECOY, 2, P.DB),
                    Sig("pfam", "PF14396", "CFTR_R", L.FAMILY, 1, P.CURATED,
                        "the regulatory R domain — the one CFTR-specific feature")),
        exemplars=(Ex("Hs_CFTR", "CFTR", "Homo sapiens", "P13569"),),
        human_genes=("CFTR",),
        length_band_aa=(1450, 1500),
        confusable_with=("H11",),
        notes=("Gated by ATP binding at the NBD dimer interface and licensed by "
               "PKA phosphorylation of R. The catalogue's proof that fold does "
               "not determine function: the same architecture as its ABCC "
               "siblings, which are transporters, plus SUR1/SUR2, which "
               "regulate a channel without conducting."),
    ),
    CF(
        key="tweety",
        name="Tweety homologues (TTYH1–TTYH3)",
        superfamily="tweety", status=St.CHANNEL_CONTESTED, fold=Fold.UNKNOWN,
        selectivity=Sel.ANION, gating=(G.VOLUME, G.CALCIUM),
        stoichiometry="dimer", tm_per_subunit=5,
        signatures=(Sig("pfam", "PF04906", "Tweety", L.FAMILY, 1, P.CURATED),),
        exemplars=(Ex("Hs_TTYH1", "TTYH1", "Homo sapiens", "Q9H313"),),
        human_genes=("TTYH1", "TTYH2", "TTYH3"),
        length_band_aa=(450, 550),
        notes=("Proposed as the volume-regulated and Ca2+-activated anion "
               "channels before LRRC8 and TMEM16A took those roles; recent "
               "structures show no obvious ion-conduction pathway. Kept, "
               "flagged contested, and named as an open question — dropping a "
               "family because the literature moved would make the census a "
               "record of opinion."),
    ),
]
