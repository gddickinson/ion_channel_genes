"""P-loop superfamily, potassium branch — the largest division of the catalogue.

Roughly 78 of the ~330 human pore-forming ion-channel genes are potassium
channels, and they are the reason the catalogue records signatures with a
*level*: every family here carries `PF00520` (Ion_trans) or its 2TM cousin
`PF07885` (Ion_trans_2), and so does every sodium, calcium, TRP, CNG and
HCN channel. A Pfam hit places a protein in the P-loop superfamily and
tells you nothing about which family it belongs to (**D25**).

What *does* separate them, in the order the classifier tries it:

1. **Topology.** 6TM/1P (Kv, KCNQ, EAG, Slo, SK), 2TM/1P (Kir), 4TM/2P
   (K2P). The pore-loop count is the cleanest split in the whole
   superfamily and survives at any evolutionary distance.
2. **Accessory architecture.** The T1 tetramerisation domain (`PF02214`)
   for Shaker-related Kv; the cyclic-nucleotide-binding homology domain
   (`PF00027`) plus PAS for EAG/ERG/ELK; the RCK domains for Slo.
3. **The selectivity filter itself.** T-x-G-Y-G is the K+ signature and is
   the one motif in this superfamily that is both diagnostic and ancient —
   it is as recognisable in *KcsA* as in Kv1.1. Where the filter degenerates
   (the silent Kv5/6/8/9 modifiers, `KCNK7`, `KCNJ18`) the family call falls
   back to reference identity, and the record is marked.

Signatures are declared only where confidence is high. Anything the S0
verifier finds on an exemplar that is not declared here is reported as an
*observed* architecture rather than silently adopted, so this file stays a
set of falsifiable claims instead of a copy of InterPro.
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

# --- shared superfamily signatures ---------------------------------------
ION_TRANS = Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 1, P.DB,
                "6TM/1P pore module; shared by Kv, Nav, Cav, TRP, CNG, HCN")
ION_TRANS_2 = Sig("pfam", "PF07885", "Ion_trans_2", L.SUPERFAMILY, 1, P.DB,
                  "2TM/1P pore module; Kir, K2P, and the prokaryotic KcsA fold")

PLOOP_SUPERFAMILY = SF(
    key="ploop",
    name="P-loop (voltage-gated-like) channel superfamily",
    fold=Fold.P_LOOP,
    alignable=True,
    anchor_module="S5–P-loop–S6 (the pore module, ~120 aa)",
    shared_signatures=(ION_TRANS, ION_TRANS_2),
    root_with=("kcsa_prok",),
    module_rule="pore_loop",
    notes=("The one superfamily large enough to need its own within-superfamily "
           "tree (tier 2 of the phylogeny protocol). Full-length alignment is "
           "meaningless across it — a Nav subunit is 2,000 aa and a Kir is 400 — "
           "so tier 2 aligns the pore module only, and the tree is explicitly a "
           "pore-module tree, not a protein tree."),
)

# --- families -------------------------------------------------------------
FAMILIES: list[CF] = [
    CF(
        key="kv_shaker",
        name="Shaker-related voltage-gated K+ channels (Kv1–Kv4)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.VOLTAGE,),
        stoichiometry="tetramer", tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF02214", "BTB_2", L.FAMILY, 1, P.DB,
                        "the T1 tetramerisation domain, which Pfam files as a "
                        "BTB/POZ domain — also carried by the 25 KCTD proteins, "
                        "which are not channels (hazard H12, measured)"),
                    Sig("pfam", "PF03521", "Kv2channel", L.SUBFAMILY, 2, P.DB,
                        "KCNB (Kv2) only")),
        exemplars=(Ex("Hs_KCNA1", "KCNA1", "Homo sapiens", "Q09470"),
                   Ex("Hs_KCNB1", "KCNB1", "Homo sapiens", "Q14721"),
                   Ex("Hs_KCND2", "KCND2", "Homo sapiens", "Q9NZV8"),
                   Ex("Dm_Shaker", "Sh", "Drosophila melanogaster", "P08510",
                      "the founding member")),
        human_genes=("KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6",
                     "KCNA7", "KCNA10", "KCNB1", "KCNB2", "KCNC1", "KCNC2",
                     "KCNC3", "KCNC4", "KCND1", "KCND2", "KCND3"),
        length_band_aa=(400, 1100), filter_motif="TVGYG",
        iuphar_class="Voltage-gated potassium channels",
        confusable_with=("kv_modifier", "H12"),
        notes="Kv1 = KCNA, Kv2 = KCNB, Kv3 = KCNC, Kv4 = KCND.",
    ),
    CF(
        key="kv_modifier",
        name="Electrically silent Kv modifier subunits (Kv5, Kv6, Kv8, Kv9)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.VOLTAGE,),
        stoichiometry="obligate heterotetramer with Kv2",
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS, Sig("pfam", "PF02214", "BTB_2", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNS3", "KCNS3", "Homo sapiens", "Q9BQ31"),
                   Ex("Hs_KCNV2", "KCNV2", "Homo sapiens", "Q8TDN2")),
        human_genes=("KCNF1", "KCNG1", "KCNG2", "KCNG3", "KCNG4",
                     "KCNS1", "KCNS2", "KCNS3", "KCNV1", "KCNV2"),
        length_band_aa=(400, 600), filter_motif="degenerate",
        confusable_with=("kv_shaker",),
        notes=("Pore-forming by fold but not by function: they cannot conduct "
               "alone and only reach the membrane with Kv2. Counted in the "
               "census as channels with a flag, because excluding them would "
               "make the census a claim about function, not about genes."),
    ),
    CF(
        key="kv_kcnq",
        name="KCNQ (Kv7) channels",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.VOLTAGE, G.LIPID),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF03520", "KCNQ_channel", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNQ1", "KCNQ1", "Homo sapiens", "P51787"),
                   Ex("Hs_KCNQ2", "KCNQ2", "Homo sapiens", "O43526")),
        human_genes=("KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5"),
        length_band_aa=(600, 1000), filter_motif="TIGYG",
        notes="PIP2-dependent; KCNQ1 needs KCNE1 (hazard H11).",
    ),
    CF(
        key="kv_eag",
        name="EAG/ERG/ELK channels (Kv10–Kv12)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.VOLTAGE,),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF00027", "cNMP_binding", L.FAMILY, 1, P.CURATED,
                        "cyclic-nucleotide-binding *homology* domain — does not "
                        "bind cNMP here; shared with CNG/HCN, which do (H4)"),
                    Sig("pfam", "PF13426", "PAS_9", L.FAMILY, 1, P.DB,
                        "the N-terminal PAS/eag domain")),
        exemplars=(Ex("Hs_KCNH2", "KCNH2", "Homo sapiens", "Q12809", "hERG"),
                   Ex("Hs_KCNH1", "KCNH1", "Homo sapiens", "O95259")),
        human_genes=("KCNH1", "KCNH2", "KCNH3", "KCNH4", "KCNH5",
                     "KCNH6", "KCNH7", "KCNH8"),
        length_band_aa=(900, 1200), filter_motif="SVGFG",
        confusable_with=("cng", "hcn"),
        notes=("The cNBHD is the sharpest architecture trap in the K+ branch: "
               "PF00027 + PF00520 fits EAG, CNG and HCN equally."),
    ),
    CF(
        key="kca_slo",
        name="Slo channels (BK, Slo3, Slack/Slick)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.VOLTAGE, G.CALCIUM),
        tm_per_subunit=7, pore_loops_per_subunit=1,
        signatures=(ION_TRANS,
                    Sig("pfam", "PF03493", "BK_channel_a", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF22614", "Slo-like_RCK", L.FAMILY, 2, P.DB,
                        "the RCK gating ring; two per subunit"),
                    Sig("pfam", "PF21014", "Slowpoke_C", L.SUBFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNMA1", "KCNMA1", "Homo sapiens", "Q12791"),
                   Ex("Hs_KCNT1", "KCNT1", "Homo sapiens", "Q5JUK3")),
        human_genes=("KCNMA1", "KCNU1", "KCNT1", "KCNT2"),
        length_band_aa=(1100, 1400), filter_motif="TVGYG",
        notes=("BK is the only human channel with a seventh TM (S0) and an "
               "extracellular N-terminus. KCNT1/2 are Na+-activated, not "
               "Ca2+-activated, despite the family name. Measured 2026-08-19: "
               "KCNMA1 carries `PF00520` and KCNT1 carries `PF07885` instead "
               "— one family, two different pore models, which is hazard H8 "
               "inside a single family."),
    ),
    CF(
        key="kca_sk",
        name="Small/intermediate-conductance Ca2+-activated K+ channels (SK, IK)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.CALCIUM,),
        tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(ION_TRANS_2,
                    Sig("pfam", "PF03530", "SK_channel", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF02888", "CaMBD", L.FAMILY, 1, P.DB,
                        "the constitutively bound calmodulin is the Ca2+ sensor")),
        exemplars=(Ex("Hs_KCNN2", "KCNN2", "Homo sapiens", "Q9H2S1"),
                   Ex("Hs_KCNN4", "KCNN4", "Homo sapiens", "O15554")),
        human_genes=("KCNN1", "KCNN2", "KCNN3", "KCNN4"),
        length_band_aa=(400, 750), filter_motif="TVGYG",
        notes=("Voltage-insensitive despite an intact S4: gated by calmodulin "
               "bound constitutively to the C-terminus. A worked example of why "
               "gating cannot be read off the architecture. Measured surprise: "
               "Pfam annotates the SK pore with `PF07885` (the 2TM model) even "
               "though the subunit is 6TM, so a rule keyed on `PF00520` for "
               "'6TM K+ channel' misses the whole family."),
    ),
    CF(
        key="kir",
        name="Inward-rectifier K+ channels (Kir1–Kir7)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.LIGAND_INTRACELLULAR, G.LIPID),
        tm_per_subunit=2, pore_loops_per_subunit=1,
        signatures=(Sig("pfam", "PF01007", "IRK", L.FAMILY, 1, P.DB,
                        "inward-rectifier transmembrane domain. Measured "
                        "2026-08-19: KCNJ2 and KCNJ11 carry NO `PF07885` at "
                        "all, so Kir and K2P are not distinguished by the "
                        "copy number of one model — the K2P test is `PF07885`"
                        "×2 and the Kir test is this accession (H16)"),
                    Sig("pfam", "PF08466", "IRK_N", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF17655", "IRK_C", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_KCNJ2", "KCNJ2", "Homo sapiens", "P63252"),
                   Ex("Hs_KCNJ11", "KCNJ11", "Homo sapiens", "Q14654", "Kir6.2, the SUR partner")),
        human_genes=("KCNJ1", "KCNJ2", "KCNJ3", "KCNJ4", "KCNJ5", "KCNJ6",
                     "KCNJ8", "KCNJ9", "KCNJ10", "KCNJ11", "KCNJ12", "KCNJ13",
                     "KCNJ14", "KCNJ15", "KCNJ16", "KCNJ18"),
        length_band_aa=(350, 500), filter_motif="TIGYG",
        confusable_with=("k2p", "H11"),
        notes=("Kir6.1/6.2 are pores whose ATP sensitivity lives in SUR1/SUR2 "
               "(ABCC8/ABCC9) — an ABC transporter that is a channel subunit "
               "and not a channel (hazard H11)."),
    ),
    CF(
        key="k2p",
        name="Two-pore-domain K+ channels (K2P)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.LEAK, G.MECHANICAL, G.PROTON, G.LIPID),
        stoichiometry="dimer", tm_per_subunit=4, pore_loops_per_subunit=2,
        signatures=(Sig("pfam", "PF07885", "Ion_trans_2", L.SUPERFAMILY, 2, P.DB,
                        "two copies per subunit — the topological signature"),),
        exemplars=(Ex("Hs_KCNK2", "KCNK2", "Homo sapiens", "O95069", "TREK-1"),
                   Ex("Hs_KCNK3", "KCNK3", "Homo sapiens", "O14649", "TASK-1")),
        human_genes=("KCNK1", "KCNK2", "KCNK3", "KCNK4", "KCNK5", "KCNK6",
                     "KCNK7", "KCNK9", "KCNK10", "KCNK12", "KCNK13",
                     "KCNK15", "KCNK16", "KCNK17", "KCNK18"),
        length_band_aa=(300, 550), filter_motif="TxGYG + TxGFG",
        confusable_with=("kir",),
        notes=("Dimeric, not tetrameric: two subunits × two pore loops still "
               "makes the fourfold-symmetric filter. The copy count of "
               "PF07885 separates K2P from Kir without any reference "
               "sequence — the cheapest reliable test in the catalogue."),
    ),
    CF(
        key="kcsa_prok",
        name="Prokaryotic K+ channels (KcsA, MthK, NaK, KirBac)",
        superfamily="ploop", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.K, gating=(G.PROTON, G.CALCIUM, G.LIGAND_INTRACELLULAR),
        tm_per_subunit=2, pore_loops_per_subunit=1,
        signatures=(ION_TRANS_2,),
        exemplars=(Ex("Sl_KcsA", "kcsA", "Streptomyces lividans", "P0A334"),
                   Ex("Mt_MthK", "mthK", "Methanothermobacter thermautotrophicus", "O27564"),
                   Ex("Bc_NaK", "nak", "Bacillus cereus", "Q81HW2")),
        other_genes=("kcsA", "mthK", "nak", "kirBac1.1"),
        length_band_aa=(150, 350), filter_motif="TVGYG",
        notes=("The rooting outgroup for the whole potassium branch and the "
               "structural reference for the selectivity filter. Keeping the "
               "prokaryotic channels in the catalogue is what makes the "
               "superfamily tree rootable at all."),
    ),
]
