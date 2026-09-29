"""The ligand-gated superfamilies: Cys-loop, iGluR, P2X and DEG/ENaC.

Four superfamilies, four independent origins, four different oligomeric
states — pentamer, tetramer, trimer, trimer. None of them is homologous to
the P-loop channels or to each other, which is the single most important
structural fact in this project: **there is no alignment that contains a
nicotinic receptor and a Kv channel**, so there is no tree either. The
phylogeny protocol handles each as its own tier-2 unit and compares them
only through the fold network (**D27**).

Each superfamily has its own trap, and each trap is a *ligand-binding domain
shared with something that is not a channel*:

* **Cys-loop.** `PF02931` (Neur_chan_LBD) is a complete soluble protein in
  the molluscan acetylcholine-binding proteins. A hit on the LBD alone
  means "cholinergic ligand-binding fold", not "channel"; the pentameric
  channel call requires the transmembrane domain `PF02932` as well (H2).
* **iGluR.** The clamshell `PF01094` (ANF_receptor) is the ancestral
  bacterial periplasmic binding protein, and it is also the ligand-binding
  domain of the class C GPCRs — every metabotropic glutamate receptor, both
  GABA-B subunits, the calcium-sensing receptor and the sweet/umami taste
  receptors. Searching for iGluRs by clamshell alone returns a receptor
  family that is not ionotropic at all (H3).
* **P2X and DEG/ENaC** are cleaner: `PF00864` and `PF00858` are effectively
  family-specific. Their difficulty is at the other end — assigning the
  invertebrate and non-metazoan members, where the human-anchored reference
  panel has nothing close.
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

NEUR_LBD = Sig("pfam", "PF02931", "Neur_chan_LBD", L.SUPERFAMILY, 1, P.DB,
               "also a complete soluble protein in AChBP — never a channel "
               "call on its own (H2)")
NEUR_MEMB = Sig("pfam", "PF02932", "Neur_chan_memb", L.SUPERFAMILY, 1, P.DB,
                "the four-helix TM bundle; this is what makes it a channel")
LIG_CHAN = Sig("pfam", "PF00060", "Lig_chan", L.SUPERFAMILY, 1, P.DB,
               "iGluR transmembrane + pore region (an inverted P-loop)")
LIG_CHAN_GLU = Sig("pfam", "PF10613", "Lig_chan-Glu_bd", L.SUPERFAMILY, 1, P.CURATED)
ANF_RECEPTOR = Sig("pfam", "PF01094", "ANF_receptor", L.SHARED_WITH_DECOY, 1, P.DB,
                   "the Venus flytrap clamshell — shared with every class C "
                   "GPCR and with bacterial periplasmic binding proteins (H3)")

SUPERFAMILIES: list[SF] = [
    SF("cysloop", "Pentameric ligand-gated ion channels (Cys-loop receptors)",
       Fold.CYS_LOOP, True, "LBD + M1–M4 (full-length alignable within superfamily)",
       (NEUR_LBD, NEUR_MEMB), ("plgic_prok",),
       notes=("Alignable end to end across the whole superfamily — the "
              "friendliest tier-2 unit in the catalogue, and the one where a "
              "conventional phylogeny is straightforwardly correct.")),
    SF("iglur", "Ionotropic glutamate receptors", Fold.IGLUR, True,
       "M1–P–M3 pore module for tier 2 (ATD, S1S2 and M4 flank it; the "
       "GluR0 root has no ATD)",
       (ANF_RECEPTOR, LIG_CHAN, LIG_CHAN_GLU), ("iglur_prok",),
       notes=("The pore module is an inverted P-loop: structurally related to "
              "Kir, sequence-undetectably so. A fold-network edge, never a "
              "tree edge (D27)."),
       module_rule="pore_loop"),
    SF("p2x", "P2X purinergic receptors", Fold.P2X, True,
       "full-length", (Sig("pfam", "PF00864", "P2X_receptor", L.SUPERFAMILY, 1, P.DB),),
       ("p2x_nonmetazoan",)),
    SF("deg_enac", "Degenerin / epithelial sodium channel superfamily",
       Fold.ENAC_DEG, True, "full-length",
       (Sig("pfam", "PF00858", "ASC", L.SUPERFAMILY, 1, P.DB),),
       ("deg_invertebrate",)),
]

FAMILIES: list[CF] = [
    # ---------------------------------------------------------- Cys-loop
    CF(
        key="nachr",
        name="Nicotinic acetylcholine receptors",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Hs_CHRNA1", "CHRNA1", "Homo sapiens", "P02708"),
                   Ex("Hs_CHRNA7", "CHRNA7", "Homo sapiens", "", "homopentameric"),
                   Ex("Tm_nAChR_a", "CHRNA1", "Torpedo marmorata", "P02711",
                      "the electric-ray receptor, first channel ever purified")),
        human_genes=("CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5",
                     "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10",
                     "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4",
                     "CHRND", "CHRNE", "CHRNG"),
        length_band_aa=(430, 630), iuphar_class="Nicotinic acetylcholine receptors",
        confusable_with=("H2",),
        notes="Cation-selective; the muscle receptor is the α2βδε/γ pentamer.",
    ),
    CF(
        key="gabaa",
        name="GABA-A receptors",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.ANION, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Hs_GABRA1", "GABRA1", "Homo sapiens", "P14867"),
                   Ex("Hs_GABRB3", "GABRB3", "Homo sapiens")),
        human_genes=("GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6",
                     "GABRB1", "GABRB2", "GABRB3", "GABRG1", "GABRG2", "GABRG3",
                     "GABRD", "GABRE", "GABRP", "GABRQ",
                     "GABRR1", "GABRR2", "GABRR3"),
        length_band_aa=(430, 640),
        confusable_with=("glyr", "H2", "H14"),
        notes=("Anion-selective, unlike the nicotinic receptors: the charge "
               "selectivity of a Cys-loop channel is set by a short "
               "intracellular ring at the M1–M2 boundary and flips with a few "
               "substitutions — so selectivity here is a *predicted* label, "
               "and the catalogue marks it as such. GABRR1-3 are the rho "
               "subunits once called GABA-C."),
    ),
    CF(
        key="glyr",
        name="Glycine receptors",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.ANION, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Hs_GLRA1", "GLRA1", "Homo sapiens", "P23415"),),
        human_genes=("GLRA1", "GLRA2", "GLRA3", "GLRA4", "GLRB"),
        length_band_aa=(440, 530),
        notes="GLRA4 is a pseudogene in human and functional in mouse.",
    ),
    CF(
        key="ht3",
        name="5-HT3 serotonin receptors",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Hs_HTR3A", "HTR3A", "Homo sapiens", "P46098"),),
        human_genes=("HTR3A", "HTR3B", "HTR3C", "HTR3D", "HTR3E"),
        length_band_aa=(440, 490),
        notes=("The only ionotropic serotonin receptor; the other thirteen "
               "5-HT receptors are GPCRs sharing the name and nothing else — a "
               "pure nomenclature hazard, and the reason gene-symbol matching "
               "is never a classification route in this project."),
    ),
    CF(
        key="zac",
        name="Zinc-activated channel (ZACN)",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Hs_ZACN", "ZACN", "Homo sapiens"),),
        human_genes=("ZACN",),
        length_band_aa=(400, 450), confusable_with=("H2",),
        notes=("Functional in human, pseudogenised in rodents. Measured "
               "2026-08-19: InterPro annotates ZACN with `PF02931` and **no** "
               "`PF02932` — so the H2 rule 'a pentameric channel call requires "
               "the TM domain as well' rejects a real channel. The rule keeps "
               "its priority and the exception is recorded, because a rule that "
               "is quietly relaxed to fit one protein stops being a test."),
    ),
    CF(
        key="plgic_invertebrate",
        name="Invertebrate anion-selective Cys-loop channels (GluCl, HisCl, pHCl)",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.ANION, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Ce_glc-1", "glc-1", "Caenorhabditis elegans"),
                   Ex("Dm_GluClalpha", "GluClalpha", "Drosophila melanogaster")),
        other_genes=("glc-1", "avr-14", "GluClalpha", "HisCl1"),
        length_band_aa=(400, 500),
        notes=("Glutamate-gated chloride channels: the ivermectin target, "
               "absent from vertebrates. A superfamily-level presence/absence "
               "result that only a non-human census can see."),
    ),
    CF(
        key="plgic_prok",
        name="Prokaryotic pentameric channels (GLIC, ELIC)",
        superfamily="cysloop", status=St.CHANNEL, fold=Fold.CYS_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.PROTON, G.LIGAND_EXTRACELLULAR),
        stoichiometry="pentamer", tm_per_subunit=4,
        signatures=(NEUR_LBD, NEUR_MEMB),
        exemplars=(Ex("Gv_GLIC", "glvI", "Gloeobacter violaceus", "Q7NDN8"),
                   Ex("Ec_ELIC", "elic", "Dickeya chrysanthemi", "P0C7B7")),
        other_genes=("GLIC", "ELIC"),
        length_band_aa=(300, 360),
        notes=("They lack the eponymous Cys-loop disulfide entirely — the "
               "superfamily's defining feature is not universal within it. "
               "The rooting outgroup for every Cys-loop tree."),
    ),
    # ------------------------------------------------------------- iGluR
    CF(
        key="ampa",
        name="AMPA receptors (GluA1–4)",
        superfamily="iglur", status=St.CHANNEL, fold=Fold.IGLUR,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="tetramer (dimer of dimers)", tm_per_subunit=3,
        pore_loops_per_subunit=1,
        signatures=(ANF_RECEPTOR, LIG_CHAN, LIG_CHAN_GLU),
        exemplars=(Ex("Hs_GRIA1", "GRIA1", "Homo sapiens", "P42261"),
                   Ex("Hs_GRIA2", "GRIA2", "Homo sapiens", "", "Q/R-edited")),
        human_genes=("GRIA1", "GRIA2", "GRIA3", "GRIA4"),
        length_band_aa=(880, 910),
        confusable_with=("H3",),
        notes=("GRIA2 Ca2+ permeability is set by RNA editing, not by sequence: "
               "a genomic census cannot see it. Recorded as a limit, not "
               "silently ignored."),
    ),
    CF(
        key="kainate",
        name="Kainate receptors (GluK1–5)",
        superfamily="iglur", status=St.CHANNEL, fold=Fold.IGLUR,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="tetramer", tm_per_subunit=3, pore_loops_per_subunit=1,
        signatures=(ANF_RECEPTOR, LIG_CHAN, LIG_CHAN_GLU),
        exemplars=(Ex("Hs_GRIK1", "GRIK1", "Homo sapiens", "P39086"),),
        human_genes=("GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5"),
        length_band_aa=(900, 980), confusable_with=("H3",),
    ),
    CF(
        key="nmda",
        name="NMDA receptors (GluN1, GluN2A–D, GluN3A–B)",
        superfamily="iglur", status=St.CHANNEL, fold=Fold.IGLUR,
        selectivity=Sel.CA, gating=(G.LIGAND_EXTRACELLULAR, G.VOLTAGE),
        stoichiometry="obligate heterotetramer", tm_per_subunit=3,
        pore_loops_per_subunit=1,
        signatures=(ANF_RECEPTOR, LIG_CHAN, LIG_CHAN_GLU),
        exemplars=(Ex("Hs_GRIN1", "GRIN1", "Homo sapiens", "Q05586"),
                   Ex("Hs_GRIN2B", "GRIN2B", "Homo sapiens")),
        human_genes=("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D",
                     "GRIN3A", "GRIN3B"),
        length_band_aa=(900, 1500), confusable_with=("H3",),
        notes=("Coincidence detector: glutamate plus glycine plus depolarisation "
               "to relieve the Mg2+ block. The voltage dependence comes from a "
               "blocking ion, not from a voltage sensor — the catalogue records "
               "gating mechanism, not just gating stimulus."),
    ),
    CF(
        key="delta_glur",
        name="Delta receptors (GluD1, GluD2)",
        superfamily="iglur", status=St.CHANNEL_CONTESTED, fold=Fold.IGLUR,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.UNKNOWN,),
        stoichiometry="tetramer", tm_per_subunit=3, pore_loops_per_subunit=1,
        signatures=(ANF_RECEPTOR, LIG_CHAN, LIG_CHAN_GLU),
        exemplars=(Ex("Hs_GRID2", "GRID2", "Homo sapiens"),),
        human_genes=("GRID1", "GRID2"),
        length_band_aa=(950, 1030),
        notes=("Orphan receptors that do not gate to glutamate and may work as "
               "synaptic adhesion molecules; the Lurcher mutation makes GRID2 "
               "constitutively open, which is what shows the pore is real. "
               "Counted, flagged contested."),
    ),
    CF(
        key="iglur_nonvertebrate",
        name="Non-vertebrate iGluRs (plant GLR, insect ionotropic receptors)",
        superfamily="iglur", status=St.CHANNEL, fold=Fold.IGLUR,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="tetramer", tm_per_subunit=3, pore_loops_per_subunit=1,
        signatures=(LIG_CHAN,),
        exemplars=(Ex("At_GLR3.4", "GLR3.4", "Arabidopsis thaliana"),
                   Ex("Dm_IR25a", "Ir25a", "Drosophila melanogaster")),
        other_genes=("GLR3.3", "GLR3.4", "Ir8a", "Ir25a", "Ir76b"),
        length_band_aa=(600, 1000),
        notes=("The insect IRs are chemosensory receptors derived from iGluRs "
               "that mostly lost the clamshell — searching for them by "
               "`PF01094` finds nothing. Recall for this family is the honest "
               "test of whether the pipeline is domain-driven or "
               "reference-driven."),
    ),
    CF(
        key="iglur_prok",
        name="Prokaryotic glutamate receptor (GluR0)",
        superfamily="iglur", status=St.CHANNEL, fold=Fold.IGLUR,
        selectivity=Sel.K, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="tetramer", tm_per_subunit=2, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF07885", "Ion_trans_2", L.FAMILY, 1, P.DB,
                        "measured: GluR0 carries the *potassium-channel* pore "
                        "model and no `PF00060` at all"),
                    Sig("pfam", "PF00497", "SBP_bac_3", L.FAMILY, 1, P.DB,
                        "the bacterial solute-binding domain the iGluR "
                        "clamshell descends from")),
        exemplars=(Ex("Ss_GluR0", "glr0", "Synechocystis sp. PCC 6803", "P73797"),),
        other_genes=("GluR0",),
        length_band_aa=(300, 400), filter_motif="TVGYG",
        notes=("A K+-selective glutamate receptor. **Measured 2026-08-19: its "
               "Pfam architecture is `PF07885` + `PF00497` — a potassium "
               "channel pore fused to a bacterial binding protein, with none "
               "of the eukaryotic iGluR models.** That is the inverted-P-loop "
               "story visible directly in the annotation, and it means the "
               "superfamily's own outgroup is unreachable by any search built "
               "on iGluR signatures — S3's profile methods have to find it."),
    ),
    # --------------------------------------------------------------- P2X
    CF(
        key="p2x",
        name="P2X receptors",
        superfamily="p2x", status=St.CHANNEL, fold=Fold.P2X,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="trimer", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF00864", "P2X_receptor", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_P2RX7", "P2RX7", "Homo sapiens", "Q99572"),
                   Ex("Hs_P2RX2", "P2RX2", "Homo sapiens", "Q9UBL9")),
        human_genes=("P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX5", "P2RX6", "P2RX7"),
        length_band_aa=(380, 600),
        notes=("ATP-gated, trimeric, with both termini intracellular and a "
               "large disulfide-rich ectodomain. P2RX7's long C-terminus "
               "drives the macropore/inflammasome phenotype."),
    ),
    CF(
        key="p2x_nonmetazoan",
        name="Non-metazoan P2X receptors",
        superfamily="p2x", status=St.CHANNEL, fold=Fold.P2X,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIGAND_EXTRACELLULAR,),
        stoichiometry="trimer", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF00864", "P2X_receptor", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Dd_P2XA", "p2xA", "Dictyostelium discoideum"),),
        other_genes=("p2xA", "p2xB", "p2xE"),
        length_band_aa=(350, 450),
        notes=("Present in amoebozoa and algae, absent from fungi, plants, "
               "*C. elegans* and *Drosophila* — a patchy distribution that "
               "makes P2X the best test case for distinguishing gene loss from "
               "database absence."),
    ),
    # ---------------------------------------------------------- DEG/ENaC
    CF(
        key="enac",
        name="Epithelial sodium channels (ENaC)",
        superfamily="deg_enac", status=St.CHANNEL, fold=Fold.ENAC_DEG,
        selectivity=Sel.NA, gating=(G.LEAK, G.MECHANICAL),
        stoichiometry="heterotrimer (α/β/γ)", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF00858", "ASC", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_SCNN1A", "SCNN1A", "Homo sapiens", "P37088"),),
        human_genes=("SCNN1A", "SCNN1B", "SCNN1G", "SCNN1D"),
        length_band_aa=(630, 900),
        notes=("Amiloride-sensitive; constitutively open and regulated by "
               "proteolysis rather than by a gate — the reason `Gating.LEAK` "
               "and `Gating.MECHANICAL` are both recorded."),
    ),
    CF(
        key="asic",
        name="Acid-sensing ion channels (ASIC)",
        superfamily="deg_enac", status=St.CHANNEL, fold=Fold.ENAC_DEG,
        selectivity=Sel.NA, gating=(G.PROTON,),
        stoichiometry="trimer", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF00858", "ASC", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_ASIC1", "ASIC1", "Homo sapiens", "P78348"),),
        human_genes=("ASIC1", "ASIC2", "ASIC3", "ASIC4", "ASIC5"),
        length_band_aa=(500, 600),
        notes="Formerly ACCN1–5; the symbol change is itself a census hazard.",
    ),
    CF(
        key="deg_invertebrate",
        name="Invertebrate degenerins and peptide-gated channels (MEC, FaNaC, ppk)",
        superfamily="deg_enac", status=St.CHANNEL, fold=Fold.ENAC_DEG,
        selectivity=Sel.NA, gating=(G.MECHANICAL, G.LIGAND_EXTRACELLULAR),
        stoichiometry="trimer", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF00858", "ASC", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Ce_mec-4", "mec-4", "Caenorhabditis elegans", "P24612"),
                   Ex("Ha_FaNaC", "fanac", "Cornu aspersum", "Q25011")),
        other_genes=("mec-4", "mec-10", "deg-1", "unc-8", "ppk", "FaNaC"),
        length_band_aa=(600, 950),
        notes=("The gain-of-function alleles that kill the neuron they are in "
               "gave the superfamily its name. *C. elegans* has ~30 members and "
               "*Drosophila* ~30 pickpocket genes, so the superfamily's size is "
               "strongly lineage-dependent — a real result the census must not "
               "flatten by counting human genes only."),
    ),
]
