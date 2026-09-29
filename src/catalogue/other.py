"""Single-family superfamilies: ORAI, Hv1, otopetrins, CLIC and the viroporins.

Five families with nothing in common except that each is the only one of its
kind, and each breaks a rule the rest of the catalogue relies on.

* **ORAI** is the store-operated calcium channel: a hexamer with a single
  central pore, gated not by voltage or ligand but by physical contact with
  STIM1 in the ER membrane across a junction. Its gate is in another
  organelle.
* **Hv1 (HVCN1)** is a voltage sensor with no pore domain at all — protons
  move through the S4 helix bundle itself. Measured: InterPro annotates
  HVCN1 with `PF00520`, the pore-module model, which it does not have. Any
  rule of the form "`PF00520` ⇒ P-loop pore" is therefore false, and its
  nearest relative by sequence is TPTE/TPTE2, a *phosphatase* (hazard H9).
* **Otopetrins** are proton channels built from two structurally similar
  12-TM halves with no relationship to anything else here; OTOP1 is the
  sour-taste receptor.
* **CLIC** proteins are soluble glutathione-S-transferase-fold monomers
  that insert into membranes and conduct chloride. A protein that is both
  cytosolic enzyme fold and channel breaks the assumption that membership is
  a property of the sequence rather than of the state, so CLIC is catalogued
  as contested with the reason attached.
* **Viroporins** are 60–100-residue viral proteins that oligomerise into
  ion-conducting pores. They are in the catalogue as a division of their
  own because they are unambiguously channels, they are the target of
  amantadine, and any pipeline that claims to find "all ion channels" in a
  sequence database will meet them.
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
    SF("orai", "ORAI store-operated Ca2+ channels", Fold.ORAI, True, "full-length",
       (Sig("pfam", "PF07856", "Orai-1", L.SUPERFAMILY, 1, P.DB),)),
    SF("hv", "Voltage-gated proton channels", Fold.VSD_ONLY, True,
       "the voltage-sensor domain",
       (Sig("pfam", "PF00520", "Ion_trans", L.SHARED_WITH_DECOY, 1, P.DB),),
       notes=("Sequence-alignable to the VSD of P-loop channels and to the VSD "
              "of the voltage-sensing phosphatases, which are not channels. "
              "The one place in the catalogue where the tier-2 tree has to "
              "include a non-channel outgroup to be rootable at all."),
       module_rule="tm_span"),
    SF("otopetrin", "Otopetrin proton channels", Fold.OTOPETRIN, True, "full-length",
       (Sig("pfam", "PF03189", "Otopetrin", L.SUPERFAMILY, 1, P.DB),)),
    SF("clic", "Chloride intracellular channels", Fold.GST, True, "full-length",
       (Sig("pfam", "PF22441", "CLIC-like_N", L.SUPERFAMILY, 1, P.DB),)),
    SF("viroporin", "Viroporins", Fold.VIROPORIN, False,
       "none — the member families are unrelated to each other", (), (),
       notes="A functional grouping, marked non-alignable so it is never treed."),
]

FAMILIES: list[CF] = [
    CF(
        key="orai",
        name="ORAI / CRAC channels (ORAI1–3)",
        superfamily="orai", status=St.CHANNEL, fold=Fold.ORAI,
        selectivity=Sel.CA, gating=(G.LIGAND_INTRACELLULAR,),
        stoichiometry="hexamer", tm_per_subunit=4,
        signatures=(Sig("pfam", "PF07856", "Orai-1", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_ORAI1", "ORAI1", "Homo sapiens", "Q96D31"),
                   Ex("Hs_ORAI2", "ORAI2", "Homo sapiens", "Q96SN7"),
                   Ex("Dm_Orai", "Orai", "Drosophila melanogaster")),
        human_genes=("ORAI1", "ORAI2", "ORAI3"),
        length_band_aa=(250, 310),
        notes=("Gated by STIM1 across an ER–plasma-membrane junction; ORAI1 "
               "loss causes a severe combined immunodeficiency. STIM1/STIM2 "
               "are catalogued as channel-associated — they are the gate, not "
               "the pore, and counting them would double-count the channel."),
    ),
    CF(
        key="hv1",
        name="Voltage-gated proton channel (HVCN1 / Hv1)",
        superfamily="hv", status=St.CHANNEL, fold=Fold.VSD_ONLY,
        selectivity=Sel.PROTON, gating=(G.VOLTAGE, G.PROTON),
        stoichiometry="dimer, one conduction path per subunit", tm_per_subunit=4,
        pore_loops_per_subunit=0,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on HVCN1 — which has no pore domain (H9)"),
                    Sig("pfam", "PF16799", "VGPC1_C", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_HVCN1", "HVCN1", "Homo sapiens", "Q96D96"),
                   Ex("Ci_Hv1", "HVCN1", "Ciona intestinalis", "Q1JV40")),
        human_genes=("HVCN1",),
        length_band_aa=(260, 280), filter_motif="D112 (the selectivity aspartate)",
        confusable_with=("H9",),
        notes=("Zero pore loops: the proton pathway is the voltage sensor. "
               "Sustains the respiratory burst in neutrophils by compensating "
               "the charge NOX2 moves. Its closest sequence relatives are the "
               "voltage-sensing phosphatases (TPTE, TPTE2, Ci-VSP), which "
               "carry the same VSD and a PTEN phosphatase domain instead of a "
               "conduction path — so 'has a voltage-sensor domain' is not "
               "evidence of a channel."),
    ),
    CF(
        key="otop",
        name="Otopetrin proton channels (OTOP1–3)",
        superfamily="otopetrin", status=St.CHANNEL, fold=Fold.OTOPETRIN,
        selectivity=Sel.PROTON, gating=(G.PROTON,),
        stoichiometry="dimer", tm_per_subunit=12,
        signatures=(Sig("pfam", "PF03189", "Otopetrin", L.FAMILY, 3, P.DB,
                        "three copies of the model per subunit (measured)"),),
        exemplars=(Ex("Hs_OTOP1", "OTOP1", "Homo sapiens", "Q7RTM1"),
                   Ex("Hs_OTOP2", "OTOP2", "Homo sapiens", "Q7RTS6")),
        human_genes=("OTOP1", "OTOP2", "OTOP3"),
        length_band_aa=(550, 620),
        notes=("Named for otoconia — the vestibular stones whose formation "
               "OTOP1 supports — a decade before anyone showed it was a "
               "channel. OTOP1 is the sour taste receptor. A reminder that "
               "gene names encode the first phenotype found, not the "
               "molecular function."),
    ),
    CF(
        key="clic",
        name="Chloride intracellular channels (CLIC1–6)",
        superfamily="clic", status=St.CHANNEL_CONTESTED, fold=Fold.GST,
        selectivity=Sel.ANION, gating=(G.PROTON, G.LIPID),
        stoichiometry="monomer (soluble) → oligomer (membrane-inserted)",
        tm_per_subunit=1,
        signatures=(Sig("pfam", "PF22441", "CLIC-like_N", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF13410", "GST_C_2", L.SHARED_WITH_DECOY, 1, P.DB,
                        "the glutathione-S-transferase fold; measured on CLIC1, "
                        "absent from CLIC4 — Pfam's coverage is uneven here")),
        exemplars=(Ex("Hs_CLIC1", "CLIC1", "Homo sapiens", "O00299"),
                   Ex("Hs_CLIC4", "CLIC4", "Homo sapiens", "Q9Y696")),
        human_genes=("CLIC1", "CLIC2", "CLIC3", "CLIC4", "CLIC5", "CLIC6"),
        length_band_aa=(230, 640),
        notes=("Metamorphic: a soluble GST-fold enzyme that refolds and inserts "
               "into membranes. Recent work argues the physiological activity "
               "is glutathione-dependent oxidoreductase, not conduction. Kept "
               "as contested because both claims have direct evidence, and "
               "because a catalogue that resolves live disputes by omission is "
               "not a catalogue."),
    ),
    CF(
        key="viroporin",
        name="Viroporins (influenza M2, HIV-1 Vpu, coronavirus E and 3a)",
        superfamily="viroporin", status=St.CHANNEL, fold=Fold.VIROPORIN,
        selectivity=Sel.PROTON, gating=(G.PROTON,),
        stoichiometry="tetramer (M2) / pentamer (Vpu, E)", tm_per_subunit=1,
        signatures=(Sig("pfam", "PF00599", "Flu_M2", L.SUBFAMILY, 1, P.DB),
                    Sig("pfam", "PF00558", "Vpu", L.SUBFAMILY, 1, P.DB),
                    Sig("pfam", "PF02723", "CoV_E", L.SUBFAMILY, 1, P.DB),
                    Sig("pfam", "PF11289", "bCoV_viroporin", L.SUBFAMILY, 1, P.DB)),
        exemplars=(Ex("IAV_M2", "M", "Influenza A virus"),
                   Ex("HIV1_Vpu", "vpu", "Human immunodeficiency virus type 1", "P05919"),
                   Ex("SARS2_E", "E", "Severe acute respiratory syndrome coronavirus 2")),
        other_genes=("M2", "vpu", "E", "ORF3a", "p7", "2B"),
        length_band_aa=(60, 275),
        notes=("The amantadine target. In the catalogue with four subfamily "
               "signatures rather than one family signature, because the "
               "viroporins are a convergent functional class: M2, Vpu and E "
               "are no more related to each other than any of them is to a "
               "human channel."),
    ),
]
