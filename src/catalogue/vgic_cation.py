"""P-loop superfamily, non-potassium branch — Nav, Cav, NALCN, CatSper, TPC,
CNG, HCN and the TRP families.

This is where domain architecture runs out. `PF00520×4` describes the sodium
channels, the calcium channels, NALCN and the CatSper subunits equally well:
four homologous 6TM repeats in one polypeptide, a pseudo-tetramer built by
two rounds of internal duplication from a Kv-like ancestor. Pfam cannot tell
them apart, InterPro's family entries largely restate the gene names they
were built from, and the length bands overlap. **The discrimination is the
selectivity filter** — one residue contributed by each of the four repeats:

    Nav      D / E / K / A      (the DEKA locus)
    Cav      E / E / E / E
    NALCN    E / E / K / E      (a Nav-like lysine in repeat III)
    CatSper  degenerate, per-subunit; not resolved here            [open]

The four-residue motif is read off the pore-module alignment, not off the
raw sequence, so it is a positive test that works on an unnamed gene model
in a newly assembled genome (`src/classify/motifs.py`; decision **D26**).

The TRP families are a second problem of a different kind, and a worse one
than it looks. Their N-terminal architectures (ankyrin repeats, TRPM
homology regions, the polycystin domain) differ so much that a full-length
alignment is dominated by insertions, so tier-1 TRP trees are built on the
pore module as well. But the deeper problem is that **Pfam does not cover
the TRP pore uniformly**: measured against InterPro on 2026-08-19, TRPV1,
TRPV5, TRPA1 and TRPM2 carry `PF00520`, while TRPC3, TRPM8, MCOLN1/2 and
PKD2 carry none of it — their pores are modelled by `PF08344`, `PF18139`
and `PF08016` instead. A census that enumerates P-loop channels by
`PF00520` therefore loses most of the TRP division silently, which is
exactly the kind of failure the S1 control benchmark exists to catch.
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
from .vgic_k import ION_TRANS

FAMILIES: list[CF] = [
    CF(
        key="nav",
        name="Voltage-gated sodium channels (Nav1.1–1.9, Nax)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.NA, gating=(G.VOLTAGE,),
        stoichiometry="pseudo-tetramer (four repeats in one chain)",
        tm_per_subunit=24, pore_loops_per_subunit=4,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 4, P.DB),
                    Sig("pfam", "PF06512", "Na_trans_assoc", L.FAMILY, 1, P.DB,
                        "the one architectural feature that is Nav-specific"),
                    Sig("pfam", "PF11933", "Na_trans_cytopl", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_SCN5A", "SCN5A", "Homo sapiens", "Q14524"),
                   Ex("Hs_SCN1A", "SCN1A", "Homo sapiens"),
                   Ex("Ab_NavAb", "navAb", "Aliarcobacter butzleri", "A8EVM5",
                      "prokaryotic homotetrameric Nav — the rooting outgroup")),
        human_genes=("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN7A",
                     "SCN8A", "SCN9A", "SCN10A", "SCN11A"),
        length_band_aa=(1600, 2100), filter_motif="DEKA",
        iuphar_class="Voltage-gated sodium channels",
        confusable_with=("cav", "nalcn", "catsper", "H1"),
        notes=("SCN7A (Nax) is not voltage-gated — a concentration sensor that "
               "keeps the architecture. Included, flagged: the census counts "
               "genes by descent, not by measured gating. Measured "
               "2026-08-19: the prokaryotic NavAb carries neither `PF06512` "
               "nor `PF11933` — it is a 6TM homotetramer with none of the "
               "cytoplasmic apparatus — so the Nav family rule, which "
               "requires both, cannot reach the family's own outgroup."),
    ),
    CF(
        key="cav",
        name="Voltage-gated calcium channels (Cav1–Cav3)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CA, gating=(G.VOLTAGE,),
        stoichiometry="pseudo-tetramer (four repeats in one chain)",
        tm_per_subunit=24, pore_loops_per_subunit=4,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 4, P.DB),
                    Sig("pfam", "PF08763", "Ca_chan_IQ", L.SUBFAMILY, 1, P.DB,
                        "Cav1/Cav2 only — measured absent from CACNA1G (Cav3), "
                        "which carries nothing but Ion_trans×4"),
                    Sig("pfam", "PF16905", "GPHH", L.SUBFAMILY, 1, P.DB),
                    Sig("pfam", "PF16885", "CAC1F_C", L.SUBFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_CACNA1C", "CACNA1C", "Homo sapiens", "Q13936"),
                   Ex("Hs_CACNA1G", "CACNA1G", "Homo sapiens", "", "Cav3.1, T-type")),
        human_genes=("CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E",
                     "CACNA1F", "CACNA1G", "CACNA1H", "CACNA1I", "CACNA1S"),
        length_band_aa=(1750, 2500), filter_motif="EEEE",
        iuphar_class="Voltage-gated calcium channels",
        confusable_with=("nav", "nalcn", "H1"),
        notes=("Cav3 (T-type) is closer to Nav than to Cav1/Cav2 in the "
               "N-terminal repeats but keeps the EEEE filter — a case where "
               "motif and tree disagree, and both get reported."),
    ),
    CF(
        key="nalcn",
        name="NALCN sodium leak channel",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LEAK,),
        stoichiometry="pseudo-tetramer + FAM155/UNC79/UNC80",
        tm_per_subunit=24, pore_loops_per_subunit=4,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 4, P.DB),),
        exemplars=(Ex("Hs_NALCN", "NALCN", "Homo sapiens", "Q8IZF0"),),
        human_genes=("NALCN",),
        length_band_aa=(1600, 1800), filter_motif="EEKE",
        confusable_with=("nav", "cav", "H1"),
        notes=("A single-gene family with a filter that is neither Nav's nor "
               "Cav's; the strongest single case for filter-based classification."),
    ),
    CF(
        key="catsper",
        name="CatSper sperm-specific Ca2+ channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CA, gating=(G.VOLTAGE, G.PROTON),
        stoichiometry="heterotetramer of four distinct 6TM subunits",
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,),
        exemplars=(Ex("Hs_CATSPER1", "CATSPER1", "Homo sapiens"),
                   Ex("Hs_CATSPER2", "CATSPER2", "Homo sapiens")),
        human_genes=("CATSPER1", "CATSPER2", "CATSPER3", "CATSPER4"),
        length_band_aa=(400, 800), filter_motif="[open]",
        confusable_with=("cav", "H1"),
        notes=("The only heterotetramer in the superfamily assembled from four "
               "*different* genes; each subunit is 6TM, so a naive architecture "
               "rule files them with the Kv branch. Auxiliary subunits "
               "(CATSPERB/G/D/E/Z) are in `controls.py`, not here."),
    ),
    CF(
        key="tpc",
        name="Two-pore channels (TPC1, TPC2)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.NA, gating=(G.VOLTAGE, G.LIGAND_INTRACELLULAR),
        stoichiometry="dimer of two-repeat subunits",
        tm_per_subunit=12, pore_loops_per_subunit=2,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 2, P.DB),),
        exemplars=(Ex("Hs_TPCN1", "TPCN1", "Homo sapiens"),
                   Ex("Hs_TPCN2", "TPCN2", "Homo sapiens"),
                   Ex("At_TPC1", "TPC1", "Arabidopsis thaliana")),
        human_genes=("TPCN1", "TPCN2"),
        length_band_aa=(750, 900), filter_motif="[open]",
        notes=("The evolutionary intermediate the 6TM→24TM story predicts: two "
               "repeats, not one and not four. Endolysosomal, PI(3,5)P2- and "
               "NAADP-regulated. Their placement is a tier-2 result, not an "
               "assumption."),
    ),
    CF(
        key="cng",
        name="Cyclic-nucleotide-gated channels (CNG)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.CYCLIC_NUCLEOTIDE,),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF00027", "cNMP_binding", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF16526", "CLZ", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_CNGA1", "CNGA1", "Homo sapiens", "P29973"),
                   Ex("Hs_CNGB1", "CNGB1", "Homo sapiens")),
        human_genes=("CNGA1", "CNGA2", "CNGA3", "CNGA4", "CNGB1", "CNGB3"),
        length_band_aa=(550, 1400), filter_motif="[open]",
        confusable_with=("hcn", "kv_eag", "H4"),
        notes="Voltage-insensitive despite a full S4. Ligand binds the CNBD.",
    ),
    CF(
        key="hcn",
        name="Hyperpolarisation-activated cyclic-nucleotide-gated channels (HCN)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.VOLTAGE, G.CYCLIC_NUCLEOTIDE),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF00027", "cNMP_binding", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF08412", "Ion_trans_N", L.FAMILY, 1, P.DB,
                        "the HCN-specific N-terminal extension")),
        exemplars=(Ex("Hs_HCN1", "HCN1", "Homo sapiens", "O60741"),
                   Ex("Hs_HCN4", "HCN4", "Homo sapiens")),
        human_genes=("HCN1", "HCN2", "HCN3", "HCN4"),
        length_band_aa=(750, 1250), filter_motif="CIGYG",
        confusable_with=("cng", "kv_eag", "H4"),
        notes=("Opens on hyperpolarisation — the inverted voltage dependence is "
               "a gating property, invisible to every sequence-level test, and "
               "one of the clearest limits on what this project can classify."),
    ),
    CF(
        key="trpc",
        name="TRPC — canonical transient receptor potential channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE,
        gating=(G.LIGAND_INTRACELLULAR, G.LIPID, G.MECHANICAL),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF08344", "TRP_2", L.FAMILY, 1, P.DB,
                        "'Transient receptor ion channel II' — TRPC's pore model; "
                        "measured: TRPC3 carries NO PF00520 at all"),
                    Sig("pfam", "PF23317", "YVC1_C", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF12796", "Ank_2", L.ACCESSORY, 1, P.DB)),
        exemplars=(Ex("Hs_TRPC3", "TRPC3", "Homo sapiens"),
                   Ex("Hs_TRPC6", "TRPC6", "Homo sapiens"),
                   Ex("Dm_trp", "trp", "Drosophila melanogaster", "P19334",
                      "the founding member of the whole TRP division")),
        human_genes=("TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6", "TRPC7"),
        other_genes=("TRPC2",),
        length_band_aa=(750, 1100), filter_motif="[open]",
        notes="TRPC2 is a unitary pseudogene in human and a functional "
              "vomeronasal channel in mouse — a per-lineage presence call the "
              "census must make correctly.",
    ),
    CF(
        key="trpv",
        name="TRPV — vanilloid receptor channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE,
        gating=(G.THERMAL, G.LIGAND_EXTRACELLULAR, G.PROTON),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF00023", "Ank", L.ACCESSORY, 1, P.DB),
                    Sig("pfam", "PF12796", "Ank_2", L.ACCESSORY, 1, P.DB)),
        exemplars=(Ex("Hs_TRPV1", "TRPV1", "Homo sapiens", "Q8NER1"),
                   Ex("Hs_TRPV5", "TRPV5", "Homo sapiens", "", "Ca2+-selective")),
        human_genes=("TRPV1", "TRPV2", "TRPV3", "TRPV4", "TRPV5", "TRPV6"),
        length_band_aa=(700, 900), filter_motif="[open]",
        notes="TRPV5/6 are the only strongly Ca2+-selective TRPs.",
    ),
    CF(
        key="trpm",
        name="TRPM — melastatin-related channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE,
        gating=(G.THERMAL, G.LIGAND_INTRACELLULAR, G.CALCIUM),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF18139", "LSDAT_euk", L.FAMILY, 1, P.DB,
                        "the TRPM homology region; present across the family"),
                    Sig("pfam", "PF25508", "TRPM2", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF16519", "TRPM_tetra", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF23317", "YVC1_C", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 1, P.DB,
                        "measured present on TRPM2, absent from TRPM8 — the "
                        "pore model does not cover the family uniformly")),
        exemplars=(Ex("Hs_TRPM8", "TRPM8", "Homo sapiens", "Q7Z2W7"),
                   Ex("Hs_TRPM2", "TRPM2", "Homo sapiens", "", "NUDT9 kinase-fused"),
                   Ex("Hs_TRPM7", "TRPM7", "Homo sapiens", "", "alpha-kinase-fused")),
        human_genes=("TRPM1", "TRPM2", "TRPM3", "TRPM4", "TRPM5",
                     "TRPM6", "TRPM7", "TRPM8"),
        length_band_aa=(1100, 2100), filter_motif="[open]",
        notes=("TRPM2, TRPM6 and TRPM7 are chanzymes — a channel fused to an "
               "enzyme. Any classifier that assumes one domain-of-interest per "
               "protein misfiles them."),
    ),
    CF(
        key="trpa",
        name="TRPA1 — ankyrin transient receptor potential channel",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE,
        gating=(G.LIGAND_EXTRACELLULAR, G.THERMAL, G.MECHANICAL),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF12796", "Ank_2", L.ACCESSORY, 4, P.DB),
                    Sig("pfam", "PF13637", "Ank_4", L.ACCESSORY, 1, P.DB)),
        exemplars=(Ex("Hs_TRPA1", "TRPA1", "Homo sapiens", "O75762"),),
        human_genes=("TRPA1",),
        length_band_aa=(1100, 1200), filter_motif="[open]",
        notes="~14–16 N-terminal ankyrin repeats — the longest in the catalogue.",
    ),
    CF(
        key="trpml",
        name="TRPML — mucolipins",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.LIPID, G.PROTON),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF08016", "PKD_channel", L.SHARED_WITH_DECOY, 1, P.DB,
                        "shared with TRPP *and* with PKD1, which is not a pore (H13)"),
                    Sig("pfam", "PF21381", "MCLN_ECD", L.FAMILY, 1, P.DB,
                        "the extracytosolic domain — what separates TRPML from TRPP")),
        exemplars=(Ex("Hs_MCOLN1", "MCOLN1", "Homo sapiens", "Q9GZU1"),
                   Ex("Hs_MCOLN2", "MCOLN2", "Homo sapiens", "Q8IZK6")),
        human_genes=("MCOLN1", "MCOLN2", "MCOLN3"),
        length_band_aa=(550, 600), filter_motif="[open]",
        notes="Endolysosomal, PI(3,5)P2-gated; MCOLN1 loss causes mucolipidosis IV.",
    ),
    CF(
        key="trpp",
        name="TRPP / polycystin-2 channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL, G.CALCIUM),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF08016", "PKD_channel", L.SHARED_WITH_DECOY, 1, P.DB,
                        "shared with TRPML and with PKD1 (H13)"),
                    Sig("pfam", "PF20519", "Polycystin_dom", L.FAMILY, 1, P.DB,
                        "also on PKD1 — the discriminator is PF21381's absence"),
                    Sig("pfam", "PF18109", "Fer4_24", L.SUBFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_PKD2", "PKD2", "Homo sapiens", "Q13563"),
                   Ex("Hs_PKD2L1", "PKD2L1", "Homo sapiens")),
        human_genes=("PKD2", "PKD2L1", "PKD2L2"),
        length_band_aa=(600, 900), filter_motif="[open]",
        confusable_with=("H13",),
        notes=("Assembles 3:1 with PKD1, an 11TM protein that contributes one "
               "pore-like segment but conducts nothing alone; PKD1 and its "
               "relatives are catalogued as channel-associated (H13)."),
    ),
    CF(
        key="trpn",
        name="TRPN / NOMPC — mechanotransduction channels (no human member)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL,),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS, Sig("pfam", "PF12796", "Ank_2", L.ACCESSORY, 9, P.CURATED)),
        exemplars=(Ex("Dm_nompC", "nompC", "Drosophila melanogaster", "Q7KIQ2"),
                   Ex("Ce_trp-4", "trp-4", "Caenorhabditis elegans")),
        other_genes=("nompC", "trp-4"),
        length_band_aa=(1500, 1800), filter_motif="[open]",
        notes=("Lost in mammals. In the catalogue because a census that only "
               "counts human genes cannot make an evolutionary claim, and "
               "because its 29-ankyrin gating spring is the clearest "
               "mechanotransduction model in the superfamily."),
    ),
]
