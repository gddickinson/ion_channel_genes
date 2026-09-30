"""The things that are not ion channels — catalogued anyway, and marked.

Every one of these will turn up in a search for ion channels: they carry a
channel domain, or a channel gene name, or they are physically part of a
channel complex. Leaving them out of the catalogue would not keep them out
of the results; it would only mean the pipeline meets them with no record of
what they are. So they are here, with a `Status` that keeps them out of the
census and a note saying who mistakes them for what.

Four kinds:

`CHANNEL_ASSOCIATED` — real subunits of real channel complexes that do not
line a pore. KCNE1 is 129 residues that change KCNQ1's kinetics beyond
recognition; SUR1 is an ABC transporter that makes Kir6.2 ATP-sensitive;
barttin is required for CLC-K to reach the membrane at all. They are
essential and they are not channels, and a census that counts them inflates
by 24 % (S20: 76 human auxiliary genes against 320 pore genes).

`NON_CHANNEL_HOMOLOG` — proteins that carry a channel's diagnostic domain
and have nothing to do with conduction. Measured here, not assumed: KCTD1
carries `PF02214` (the Kv T1 domain), GRM1 carries `PF01094` (the iGluR
clamshell), POMT1 carries `PF02815` (the IP3R/RyR MIR domain), TPTE carries
`PF00520` (the pore module) while being a phosphatase.

`TRANSPORTER` — coupled transport through a channel-like fold.

`OUT_OF_SCOPE` — genuine channels that are not *ion* channels. Aquaporins
are the boundary case that defines the project's scope: a gated aqueous
pore, selective, ancient, well characterised, and it conducts water. The
scope decision is recorded (**D23**) rather than assumed, and aquaporins
stay in the catalogue as the control that proves the boundary is enforced.
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

_NONE: tuple = ()

FAMILIES: list[CF] = [
    # ------------------------------------------------ auxiliary subunits
    CF(
        key="assoc_kcne",
        name="KCNE (MinK-related) accessory subunits",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="1–4 per channel",
        tm_per_subunit=1,
        signatures=(Sig("pfam", "PF02060", "ISK_Channel", L.FAMILY, 1, P.DB,
                        "Pfam calls it a 'slow voltage-gated potassium channel' "
                        "— the name is wrong and the annotation is the trap"),),
        exemplars=(Ex("Hs_KCNE1", "KCNE1", "Homo sapiens", "P15382"),),
        human_genes=("KCNE1", "KCNE2", "KCNE3", "KCNE4", "KCNE5"),
        length_band_aa=(100, 180), confusable_with=("kv_kcnq", "H11"),
        notes=("129 aa, one TM helix, no pore — and it converts KCNQ1 into the "
               "cardiac IKs channel. Both the gene symbol and the Pfam entry "
               "name it a channel."),
    ),
    CF(
        key="assoc_k_beta",
        name="K+ channel beta and gamma subunits (KCNAB, KCNMB, LRRC26, DPP6/10)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="4 per channel",
        signatures=(Sig("pfam", "PF00248", "Aldo_ket_red", L.FAMILY, 1, P.DB,
                        "KCNAB is an aldo-keto reductase bound to the T1 domain"),
                    Sig("pfam", "PF03185", "CaKB", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF00930", "DPPIV_N", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNAB1", "KCNAB1", "Homo sapiens", "Q14722"),
                   Ex("Hs_KCNMB1", "KCNMB1", "Homo sapiens", "Q16558"),
                   Ex("Hs_DPP6", "DPP6", "Homo sapiens", "P42658")),
        human_genes=("KCNAB1", "KCNAB2", "KCNAB3", "KCNMB1", "KCNMB2",
                     "KCNMB3", "KCNMB4", "LRRC26", "LRRC38", "LRRC52",
                     "LRRC55", "DPP6", "DPP10"),
        notes=("DPP6 and DPP10 are catalytically dead peptidase homologues "
               "that set Kv4 kinetics — an enzyme fold recruited as a channel "
               "subunit, the mirror image of CLIC."),
    ),
    CF(
        key="assoc_cav_aux",
        name="Cav auxiliary subunits (beta, alpha2delta, gamma/TARP)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="1 each per channel",
        signatures=(Sig("pfam", "PF00625", "Guanylate_kin", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF12052", "VGCC_beta4Aa_N", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF08473", "VGCC_alpha2", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF00822", "PMP22_Claudin", L.FAMILY, 1, P.DB,
                        "the gamma/TARP subunits are claudins"),
                    Sig("pfam", "PF15108", "TMEM37", L.FAMILY, 1, P.DB,
                        "TMEM37 (PR1, the gamma-like subunit): its own model, "
                        "not the CACNG one (measured 2026-09-29)")),
        exemplars=(Ex("Hs_CACNB1", "CACNB1", "Homo sapiens", "Q02641"),
                   Ex("Hs_CACNA2D1", "CACNA2D1", "Homo sapiens", "P54289"),
                   Ex("Hs_CACNG2", "CACNG2", "Homo sapiens", "Q9Y698", "stargazin")),
        human_genes=("CACNB1", "CACNB2", "CACNB3", "CACNB4",
                     "CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNA2D4",
                     "CACNG1", "CACNG2", "CACNG3", "CACNG4", "CACNG5",
                     "CACNG6", "CACNG7", "CACNG8", "TMEM37"),
        notes=("The `CACNG` symbols are the worst false-friend set in the "
               "catalogue: CACNG2–8 are TARPs, AMPA-receptor subunits with a "
               "claudin fold and no role in calcium channels at all. "
               "α2δ (CACNA2D1) is the gabapentin/pregabalin target. TMEM37 "
               "(voltage-dependent calcium channel γ-like subunit) added after "
               "S20 found it on UniProt KW-0407."),
    ),
    CF(
        key="assoc_kchip",
        name="Kv channel-interacting proteins (KChIP1–4)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="4 per Kv4 channel",
        signatures=(Sig("pfam", "PF13499", "EF-hand_7", L.SHARED_WITH_DECOY, 1, P.DB,
                        "EF-hands of the neuronal calcium sensor family "
                        "(recoverin, hippocalcin, NCS1) — never evidence"),
                    Sig("pfam", "PF13833", "EF-hand_8", L.SHARED_WITH_DECOY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNIP1", "KCNIP1", "Homo sapiens", "Q9NZI2"),),
        human_genes=("KCNIP1", "KCNIP2", "KCNIP3", "KCNIP4"),
        length_band_aa=(200, 290), provenance=P.CURATED,
        notes=("Cytosolic Ca2+ sensors that bind the Kv4 T1 domain and set "
               "the A-type current; KCNIP3 is also DREAM/calsenilin, a "
               "transcriptional repressor. Added after S20 found all four on "
               "UniProt KW-0407."),
    ),
    CF(
        key="assoc_nav_beta",
        name="Nav beta subunits (SCN1B–SCN4B)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, tm_per_subunit=1,
        signatures=(Sig("pfam", "PF07686", "V-set", L.FAMILY, 1, P.DB,
                        "an immunoglobulin V-set domain — a cell-adhesion "
                        "molecule doing a channel job"),),
        exemplars=(Ex("Hs_SCN1B", "SCN1B", "Homo sapiens", "Q07699"),),
        human_genes=("SCN1B", "SCN2B", "SCN3B", "SCN4B"),
        notes="`SCNxB` symbols sort next to `SCNxA` in every gene list.",
    ),
    CF(
        key="assoc_sur",
        name="Sulfonylurea receptors (ABCC8, ABCC9)",
        superfamily="abc_channel", status=St.CHANNEL_ASSOCIATED, fold=Fold.ABC,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="4 per Kir6 channel",
        tm_per_subunit=17,
        # Generic ABC domains (~1.7 M UniProt proteins): no level at which
        # they identify SUR. There is no positive SUR test (H11, S2c).
        signatures=(Sig("pfam", "PF00664", "ABC_membrane", L.SHARED_WITH_DECOY, 2, P.DB),
                    Sig("pfam", "PF00005", "ABC_tran", L.SHARED_WITH_DECOY, 2, P.DB)),
        exemplars=(Ex("Hs_ABCC8", "ABCC8", "Homo sapiens", "Q09428"),),
        human_genes=("ABCC8", "ABCC9"),
        confusable_with=("cftr", "kir", "H11"),
        notes=("Measured: ABCC8 and CFTR carry the same two ABC Pfam domains "
               "in the same copy numbers. The only architectural difference is "
               "CFTR's R domain (`PF14396`) — so 'ABC fold' cannot decide "
               "channel from non-channel in either direction, and the "
               "discriminator is a single accession."),
    ),
    CF(
        # A decoy for profile methods, not for domain rules (S3a2). The S3a
        # best-profile assignment scores a record against every catalogued
        # family and calls the best one: with no generic ABC transporter on
        # the menu, the two-seed SUR profile won 525 census records, 508 of
        # them bacterial peptide-exporting ABC transporters. Agreement
        # between instruments is evidence only when every family that could
        # score has a profile.
        key="nonchannel_abc_transporter",
        name="ABC transporters that are neither CFTR nor a sulfonylurea receptor",
        superfamily="abc_channel", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.ABC,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="monomer or dimer",
        signatures=(Sig("pfam", "PF00664", "ABC_membrane", L.SHARED_WITH_DECOY, 1, P.DB),
                    Sig("pfam", "PF00005", "ABC_tran", L.SHARED_WITH_DECOY, 1, P.DB),
                    Sig("pfam", "PF03412", "Peptidase_C39", L.ACCESSORY, 1, P.DB,
                        "the bacterial peptide-processing exporters (HlyB, SunT) "
                        "the SUR profile was winning")),
        exemplars=(Ex("Hs_ABCC1", "ABCC1", "Homo sapiens", "P33527",
                      "MRP1 — the closest human relative of SUR, TMD0 included"),
                   Ex("Ec_HlyB", "hlyB", "Escherichia coli", "P08716",
                      "C39-peptidase ABC exporter"),
                   Ex("Bs_SunT", "sunT", "Bacillus subtilis", "P68579",
                      "C39-peptidase ABC exporter"),
                   Ex("Ec_MsbA", "msbA", "Escherichia coli", "P60752",
                      "bacterial ABC half-transporter")),
        human_genes=("ABCC1", "ABCC2", "ABCC3", "ABCC4", "ABCC5", "ABCC6",
                     "ABCC10", "ABCC11", "ABCC12"),
        confusable_with=("cftr", "assoc_sur", "H11"),
        provenance=P.DB,
        notes=("ABCC family members other than CFTR (ABCC7) and SUR1/2 "
               "(ABCC8/9), plus bacterial exporters. Exemplar Pfam "
               "architectures confirmed live 2026-09-28."),
    ),
    CF(
        key="assoc_stim",
        name="STIM proteins (STIM1, STIM2)",
        superfamily="orai", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, tm_per_subunit=1,
        signatures=(Sig("pfam", "PF07647", "SAM_2", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF16533", "SOAR", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_STIM1", "STIM1", "Homo sapiens", "Q13586"),),
        human_genes=("STIM1", "STIM2"),
        notes="The ER calcium sensor that gates ORAI. A gate is not a pore.",
    ),
    CF(
        key="assoc_mcu_reg",
        name="MCU regulators (MICU1–3, EMRE/SMDT1, MCUR1)",
        superfamily="mcu", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF13833", "EF-hand_8", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF10161", "DDDD", L.FAMILY, 1, P.DB,
                        "EMRE; Pfam still calls it a putative mitochondrial "
                        "precursor protein")),
        exemplars=(Ex("Hs_MICU1", "MICU1", "Homo sapiens", "Q9BPX6"),
                   Ex("Hs_SMDT1", "SMDT1", "Homo sapiens", "Q9H4I9")),
        human_genes=("MICU1", "MICU2", "MICU3", "SMDT1", "MCUR1"),
        notes="EMRE is required for conduction and contributes no pore lining.",
    ),
    CF(
        key="assoc_iglur_aux",
        name="iGluR auxiliary subunits (NETO, CNIH, SHISA, GSG1L)",
        superfamily="iglur", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF00431", "CUB", L.FAMILY, 2, P.DB),),
        exemplars=(Ex("Hs_NETO1", "NETO1", "Homo sapiens", "Q8TDF5"),),
        human_genes=("NETO1", "NETO2", "CNIH2", "CNIH3",
                     "SHISA6", "SHISA7", "SHISA8", "SHISA9", "GSG1L"),
        notes="With the TARPs, these set nearly every measurable AMPA property.",
    ),
    CF(
        key="assoc_catsper_aux",
        name="CatSper auxiliary subunits (CATSPERB/G/D/E/Z)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF15149", "CATSPERB_C", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_CATSPERB", "CATSPERB", "Homo sapiens", "Q9H7T0"),),
        human_genes=("CATSPERB", "CATSPERG", "CATSPERD", "CATSPERE",
                     "CATSPERZ", "TMEM249"),
        notes=("Share the `CATSPER` symbol root with the four pore subunits. "
               "Symbol-prefix matching would count nine CatSper channels where "
               "there are four pore-forming genes and one channel."),
    ),
    CF(
        key="assoc_clc_aux",
        name="CLC accessory proteins (barttin, CLCA)",
        superfamily="clc", status=St.CHANNEL_ASSOCIATED, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF15462", "Barttin", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_BSND", "BSND", "Homo sapiens", "Q8WZ55"),),
        human_genes=("BSND", "CLCA1", "CLCA2", "CLCA3P", "CLCA4"),
        notes=("The CLCA proteins are secreted metalloprotease-fold activators "
               "of TMEM16A — named 'chloride channel accessory' and related to "
               "neither CLC nor any channel."),
    ),
    CF(
        key="assoc_polycystin1",
        name="Polycystin-1 family (PKD1, PKD1L1–3)",
        superfamily="ploop", status=St.CHANNEL_ASSOCIATED, fold=Fold.P_LOOP,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        stoichiometry="1 per PKD2 trimer", tm_per_subunit=11,
        signatures=(Sig("pfam", "PF08016", "PKD_channel", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on PKD1 — the same model as TRPP and TRPML"),
                    Sig("pfam", "PF20519", "Polycystin_dom", L.SHARED_WITH_DECOY, 1, P.DB),
                    Sig("pfam", "PF00801", "PKD", L.ACCESSORY, 15, P.DB)),
        exemplars=(Ex("Hs_PKD1", "PKD1", "Homo sapiens", "P98161"),),
        human_genes=("PKD1", "PKD1L1", "PKD1L2", "PKD1L3"),
        confusable_with=("trpp", "H13"),
        notes=("4,303 aa with a 3,000-residue ectodomain and a channel-like "
               "C-terminus that contributes one subunit to a PKD2 trimer "
               "without conducting on its own. The commonest cause of "
               "autosomal dominant polycystic kidney disease. Whether the "
               "PKD1:PKD2 1:3 complex is one channel or a channel plus a "
               "receptor is the reason this sits in `controls.py` with the "
               "argument written down."),
    ),
    # ------------------------------------------- non-channel homologues
    CF(
        key="nonchannel_kctd",
        name="KCTD proteins — the T1 domain without a channel",
        superfamily="ploop", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF02214", "BTB_2", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on KCTD1 — the same model as the Kv T1 domain"),
                    # The positive H12 tests (S2c). Measured in census v2 on
                    # PF02214 carriers without a pore module: 14,354 carry one
                    # of these, and the S3a profiles call 14,317 of them KCTD
                    # and 3 anything else. Each covers one KCTD clade.
                    Sig("pfam", "PF31093", "KCTD5_C-like", L.FAMILY, 1, P.DB,
                        "KCTD2/5/17 C-terminal domain"),
                    Sig("pfam", "PF23110", "H1_KCTD8_12_16", L.FAMILY, 1, P.DB,
                        "KCTD8/12/16 H1 domain"),
                    Sig("pfam", "PF31104", "KCTD10_C-like", L.FAMILY, 1, P.DB,
                        "KCTD10/13/TNFAIP1 C-terminal domain"),
                    Sig("pfam", "PF20871", "KCTD1-15_CTD", L.FAMILY, 1, P.DB,
                        "KCTD1/15 C-terminal domain"),
                    Sig("pfam", "PF31099", "SHKBP1_KCTD3_C", L.FAMILY, 1, P.DB,
                        "SHKBP1/KCTD3 C-terminal beta-propeller")),
        exemplars=(Ex("Hs_KCTD1", "KCTD1", "Homo sapiens", "Q719H9"),),
        human_genes=("KCTD1", "KCTD2", "KCTD3", "KCTD5", "KCTD6", "KCTD7",
                     "KCTD8", "KCTD10", "KCTD11", "KCTD12", "KCTD13",
                     "KCTD15", "KCTD16", "KCTD17", "KCTD20", "KCTD21"),
        confusable_with=("kv_shaker", "H12"),
        notes=("'Potassium channel tetramerisation domain containing' — 25 "
               "genes named after a channel domain, mostly Cullin3 ubiquitin "
               "ligase adaptors, several of them GABA-B receptor subunits. "
               "None is a channel. Symbol and domain both mislead."),
    ),
    CF(
        key="nonchannel_class_c_gpcr",
        name="Class C GPCRs — the iGluR clamshell without a pore",
        superfamily="iglur", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE, tm_per_subunit=7,
        signatures=(Sig("pfam", "PF01094", "ANF_receptor", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on GRM1 — the same model as the iGluR ATD"),
                    Sig("pfam", "PF00003", "7tm_3", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_GRM1", "GRM1", "Homo sapiens", "Q13255"),),
        human_genes=("GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7",
                     "GRM8", "GABBR1", "GABBR2", "CASR",
                     "TAS1R1", "TAS1R2", "TAS1R3"),
        confusable_with=("ampa", "nmda", "H3"),
        notes=("Metabotropic, seven-TM, G-protein-coupled. `GRM` and `GRIN` "
               "differ by two letters and by a whole signalling mechanism; the "
               "shared clamshell is real homology, which is what makes a "
               "domain-only search return them."),
    ),
    CF(
        key="nonchannel_achbp",
        name="Acetylcholine-binding proteins — a soluble Cys-loop LBD",
        superfamily="cysloop", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.CYS_LOOP,
        selectivity=Sel.UNKNOWN, gating=_NONE, stoichiometry="soluble pentamer",
        tm_per_subunit=0,
        signatures=(Sig("pfam", "PF02931", "Neur_chan_LBD", L.SHARED_WITH_DECOY, 1, P.DB),),
        exemplars=(Ex("Ls_AChBP", "achbp", "Lymnaea stagnalis", "P58154"),),
        other_genes=("AChBP",),
        length_band_aa=(210, 240), confusable_with=("nachr", "H2"),
        notes=("Secreted by snail glia into the synapse; the structure that "
               "made nicotinic pharmacology interpretable. Zero TM helices, so "
               "the H2 test (LBD *and* TM domain) rejects it — while the same "
               "test wrongly rejects human ZACN, which InterPro annotates with "
               "no TM domain either. One test, one true positive and one false "
               "negative, both recorded."),
    ),
    CF(
        key="nonchannel_pomt",
        name="Protein O-mannosyltransferases (POMT1, POMT2)",
        superfamily="ca_release", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.UNKNOWN,
        selectivity=Sel.UNKNOWN, gating=_NONE,
        signatures=(Sig("pfam", "PF02815", "MIR", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on POMT1 — the same model as ITPR/RYR"),
                    Sig("pfam", "PF02366", "PMT", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_POMT1", "POMT1", "Homo sapiens", "Q9Y6A1"),),
        human_genes=("POMT1", "POMT2"),
        confusable_with=("itpr", "ryr", "H10"),
        notes=("Glycosyltransferases of the ER membrane carrying the MIR "
               "domain. Inherited as a control from the parent project, where "
               "they were the MIR-sharing decoy in the S1 benchmark."),
    ),
    CF(
        key="nonchannel_vsp",
        name="Voltage-sensing phosphatases (TPTE, TPTE2, Ci-VSP)",
        superfamily="hv", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.VSD_ONLY,
        selectivity=Sel.UNKNOWN, gating=_NONE, tm_per_subunit=4,
        signatures=(Sig("pfam", "PF00520", "Ion_trans", L.SHARED_WITH_DECOY, 1, P.DB,
                        "measured on TPTE — a pore-module annotation on a "
                        "protein with no pore"),
                    Sig("pfam", "PF10409", "PTEN_C2", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_TPTE", "TPTE", "Homo sapiens", "P56180"),
                   Ex("Ci_VSP", "Ci-VSP", "Ciona intestinalis", "Q4W8A1")),
        human_genes=("TPTE", "TPTE2"),
        confusable_with=("hv1", "H9"),
        notes=("A voltage sensor coupled to a PTEN-like phosphatase: "
               "depolarisation dephosphorylates PIP2. The single best "
               "counter-example to 'voltage sensor ⇒ voltage-gated channel', "
               "and the nearest relative of Hv1, which is a channel."),
    ),
    CF(
        key="nonchannel_gasdermin",
        name="Gasdermins — pore-forming, not ion-conducting",
        superfamily="", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.UNKNOWN,
        selectivity=Sel.LARGE_PORE, gating=_NONE,
        signatures=(Sig("pfam", "PF04598", "Gasdermin", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_GSDMD", "GSDMD", "Homo sapiens", "P57764"),),
        human_genes=("GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME", "PJVK"),
        notes=("A 27-mer β-barrel of 21 nm that kills the cell it opens in. "
               "Excluded because the census is of gated, selective conduction "
               "pathways, not of every hole in a membrane — the same rule that "
               "excludes perforin, MLKL and the bacterial toxins. Recorded so "
               "that the exclusion is a decision (**D23**) with a boundary."),
    ),
    # ----------------------------------------- transporters & out of scope
    CF(
        key="transporter_slc26",
        name="SLC26 anion transporters and prestin",
        superfamily="", status=St.TRANSPORTER, fold=Fold.UNKNOWN,
        selectivity=Sel.ANION, gating=_NONE,
        signatures=(Sig("pfam", "PF00916", "Sulfate_transp", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF01740", "STAS", L.ACCESSORY, 1, P.DB)),
        exemplars=(Ex("Hs_SLC26A5", "SLC26A5", "Homo sapiens", "P58743", "prestin"),),
        human_genes=("SLC26A1", "SLC26A2", "SLC26A3", "SLC26A4", "SLC26A5",
                     "SLC26A6", "SLC26A7", "SLC26A8", "SLC26A9", "SLC26A11"),
        notes=("Prestin is the outer-hair-cell motor: a transporter that moves "
               "charge without completing a transport cycle. SLC26A9 has "
               "measurable channel-like conductance. The family sits exactly "
               "on the channel/transporter boundary and is catalogued as "
               "transporter with the dissent noted."),
    ),
    CF(
        key="outofscope_aquaporin",
        name="Aquaporins — channels that are not ion channels",
        superfamily="", status=St.OUT_OF_SCOPE, fold=Fold.AQUAPORIN,
        selectivity=Sel.WATER, gating=(G.PROTON, G.PHOSPHORYLATION),
        stoichiometry="tetramer, one pore per subunit", tm_per_subunit=6,
        signatures=(Sig("pfam", "PF00230", "MIP", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_AQP1", "AQP1", "Homo sapiens", "P29972"),),
        human_genes=("AQP0", "AQP1", "AQP2", "AQP3", "AQP4", "AQP5", "AQP6",
                     "AQP7", "AQP8", "AQP9", "AQP10", "AQP11", "AQP12A"),
        notes=("The scope boundary made concrete. Aquaporins exclude protons "
               "by design — the exact opposite of an ion channel's job — yet "
               "AQP1 conducts cations when cGMP-gated and AQP6 conducts anions. "
               "Out of scope by decision **D23**, kept as the control that "
               "shows the boundary is enforced rather than assumed."),
    ),
    CF(
        key="outofscope_bacterial_porin",
        name="Bacterial outer-membrane porins (OmpF, OmpC)",
        superfamily="porin", status=St.OUT_OF_SCOPE, fold=Fold.BETA_BARREL,
        selectivity=Sel.CATION_NONSELECTIVE, gating=_NONE, tm_per_subunit=16,
        signatures=(Sig("pfam", "PF00267", "Porin_1", L.FAMILY, 1, P.CURATED),),
        exemplars=(Ex("Ec_OmpF", "ompF", "Escherichia coli"),),
        other_genes=("ompF", "ompC", "phoE"),
        notes=("Ungated diffusion pores. Excluded on the gating criterion, "
               "which is what keeps VDAC — voltage-gated, eukaryotic, "
               "regulated — inside the census while its architectural "
               "relatives stay out."),
    ),
]
