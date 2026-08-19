"""Channels of intracellular membranes: IP3 receptors, ryanodine receptors,
TRIC, the mitochondrial calcium uniporter, VDAC and TMEM175.

The organelle channels are systematically under-represented in every
channel census, for a reason worth naming: they were found by physiology on
isolated membranes rather than by patch clamp on cells, so they entered the
literature later, with less standard nomenclature, and several of them still
have contested molecular identities. A census that works from
plasma-membrane-shaped assumptions — a signal peptide, an extracellular
ligand-binding domain, a voltage sensor — finds fewer of them than exist.

**ITPR and RYR are inherited work.** The parent project of this codebase
(`../ip3r_genes`) is a full census of the IP3 receptors, and it established
by measurement that the two families share all four of their diagnostic Pfam
domains (`PF08709`, `PF01365`, `PF08454`, `PF02815`), that 49 % of zebrafish
`PF08709` records are ryanodine receptors, and that a discovery scorer
without a positive ITPR/RYR test promotes every RyR decoy. That result is
carried here as hazard **H4** and as decision **D14**, restated for a
catalogue where the same problem recurs a dozen times over.
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

INS145 = Sig("pfam", "PF08709", "Ins145_P3_rec", L.SUPERFAMILY, 1, P.DB,
             "carried by ITPR *and* RYR — Pfam's name says 'IP3/ryanodine' (H4)")
RIH = Sig("pfam", "PF01365", "RYDR_ITPR", L.SUPERFAMILY, 2, P.DB)
RIH_ASSOC = Sig("pfam", "PF08454", "RIH_assoc", L.SUPERFAMILY, 1, P.DB)
MIR = Sig("pfam", "PF02815", "MIR", L.SHARED_WITH_DECOY, 1, P.DB,
          "also carried by the protein O-mannosyltransferases POMT1/2 (H10)")

SUPERFAMILIES: list[SF] = [
    SF("ca_release", "Intracellular Ca2+-release channels (ITPR, RYR)",
       Fold.P_LOOP, True, "the C-terminal pore module + RIH domains",
       (INS145, RIH, RIH_ASSOC, MIR), ("itpr",),
       notes=("Their pore module is a P-loop channel — so `ca_release` is a "
              "division *inside* the P-loop clan structurally, but it is kept "
              "as its own superfamily because the 2,700–5,000 aa cytosolic "
              "solenoid makes full-length alignment to a 400 aa Kir "
              "meaningless. Tier 2 aligns the pore module and says so.")),
    SF("tric", "Trimeric intracellular cation channels", Fold.TRIC, True,
       "full-length", (Sig("pfam", "PF05197", "TRIC", L.SUPERFAMILY, 1, P.DB),)),
    SF("mcu", "Mitochondrial calcium uniporter", Fold.MCU, True, "full-length",
       (Sig("pfam", "PF04678", "MCU", L.SUPERFAMILY, 1, P.DB),)),
    SF("porin", "Eukaryotic porins (VDAC)", Fold.BETA_BARREL, True, "full-length",
       (Sig("pfam", "PF01459", "Porin_3", L.SUPERFAMILY, 1, P.DB),),
       notes=("The only β-barrel in the census. Bacterial porins (OmpF and "
              "relatives) share the architecture and are catalogued as "
              "out-of-scope: a diffusion pore with no gating and no "
              "selectivity is not what this project counts.")),
    SF("tmem175", "TMEM175 lysosomal channels", Fold.TMEM175, True, "full-length",
       (Sig("pfam", "PF06736", "TMEM175", L.SUPERFAMILY, 2, P.DB),)),
]

FAMILIES: list[CF] = [
    CF(
        key="itpr",
        name="Inositol 1,4,5-trisphosphate receptors (ITPR1–3)",
        superfamily="ca_release", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CA, gating=(G.LIGAND_INTRACELLULAR, G.CALCIUM),
        stoichiometry="tetramer", tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(INS145, RIH, RIH_ASSOC, MIR,
                    Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_ITPR1", "ITPR1", "Homo sapiens", "Q14643"),
                   Ex("Hs_ITPR2", "ITPR2", "Homo sapiens", "Q14571"),
                   Ex("Hs_ITPR3", "ITPR3", "Homo sapiens", "Q14573")),
        human_genes=("ITPR1", "ITPR2", "ITPR3"),
        length_band_aa=(2000, 3600), confusable_with=("ryr", "H4", "H10"),
        provenance=P.DB,
        notes=("Measured lengths 2,758 / 2,701 / 2,671 aa. The ER's "
               "ligand-gated calcium-release channel; ITPR1 variants cause "
               "SCA15/SCA29 and Gillespie syndrome. Fully censused by the "
               "parent project, whose results are this project's external "
               "check on the `ca_release` branch — never copied into it."),
    ),
    CF(
        key="ryr",
        name="Ryanodine receptors (RYR1–3)",
        superfamily="ca_release", status=St.CHANNEL, fold=Fold.P_LOOP,
        selectivity=Sel.CA, gating=(G.CALCIUM, G.VOLTAGE, G.LIGAND_INTRACELLULAR),
        stoichiometry="tetramer", tm_per_subunit=6, pore_loops_per_subunit=1,
        signatures=(INS145, RIH, RIH_ASSOC, MIR,
                    Sig("pfam", "PF02026", "RyR", L.FAMILY, 4, P.DB,
                        "the RyR repeat — the one domain ITPR does not carry"),
                    Sig("pfam", "PF06459", "RR_TM4-6", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF00622", "SPRY", L.ACCESSORY, 3, P.DB),
                    Sig("pfam", "PF00520", "Ion_trans", L.SUPERFAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_RYR1", "RYR1", "Homo sapiens", "P21817"),
                   Ex("Hs_RYR2", "RYR2", "Homo sapiens", "Q92736"),
                   Ex("Hs_RYR3", "RYR3", "Homo sapiens", "Q15413")),
        human_genes=("RYR1", "RYR2", "RYR3"),
        length_band_aa=(4800, 5100), confusable_with=("itpr", "H4"),
        provenance=P.DB,
        notes=("5,038 aa (measured) — the largest ion channel known, and the "
               "largest single-chain protein in the catalogue. RYR1 variants "
               "cause malignant hyperthermia and central core disease; RYR2 "
               "causes CPVT. `PF02026` is the positive discriminator against "
               "ITPR that the parent project's benchmark was missing."),
    ),
    CF(
        key="tric",
        name="Trimeric intracellular cation channels (TMEM38A/B)",
        superfamily="tric", status=St.CHANNEL, fold=Fold.TRIC,
        selectivity=Sel.K, gating=(G.VOLTAGE, G.CALCIUM),
        stoichiometry="trimer", tm_per_subunit=7,
        signatures=(Sig("pfam", "PF05197", "TRIC", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_TMEM38A", "TMEM38A", "Homo sapiens", "Q9H6F2"),
                   Ex("Hs_TMEM38B", "TMEM38B", "Homo sapiens")),
        human_genes=("TMEM38A", "TMEM38B"),
        length_band_aa=(280, 310),
        notes=("The counter-ion channel that lets the SR release calcium "
               "without building an opposing voltage — a channel whose whole "
               "function is to make another channel work. TMEM38B mutations "
               "cause a recessive osteogenesis imperfecta."),
    ),
    CF(
        key="mcu",
        name="Mitochondrial calcium uniporter (MCU, MCUB)",
        superfamily="mcu", status=St.CHANNEL, fold=Fold.MCU,
        selectivity=Sel.CA, gating=(G.CALCIUM,),
        stoichiometry="tetramer + MICU regulators", tm_per_subunit=2,
        signatures=(Sig("pfam", "PF04678", "MCU", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_MCU", "MCU", "Homo sapiens", "Q8NE86"),
                   Ex("Hs_MCUB", "MCUB", "Homo sapiens", "Q9NWR8")),
        human_genes=("MCU", "MCUB"),
        length_band_aa=(330, 360),
        notes=("The most Ca2+-selective channel measured, and molecularly "
               "unidentified until 2011 — fifty years after the current it "
               "carries was described. Its regulators (MICU1–3, EMRE/SMDT1, "
               "MCUR1) are catalogued as channel-associated: MCUB is a pore "
               "subunit, EMRE is not."),
    ),
    CF(
        key="vdac",
        name="Voltage-dependent anion channels (VDAC1–3)",
        superfamily="porin", status=St.CHANNEL, fold=Fold.BETA_BARREL,
        selectivity=Sel.ANION, gating=(G.VOLTAGE,),
        stoichiometry="monomer (19-stranded β-barrel)", tm_per_subunit=19,
        signatures=(Sig("pfam", "PF01459", "Porin_3", L.FAMILY, 1, P.DB),),
        exemplars=(Ex("Hs_VDAC1", "VDAC1", "Homo sapiens", "P21796"),
                   Ex("Hs_VDAC2", "VDAC2", "Homo sapiens", "P45880")),
        human_genes=("VDAC1", "VDAC2", "VDAC3"),
        length_band_aa=(275, 300),
        notes=("The main conduit across the outer mitochondrial membrane and "
               "the only all-β channel in the census. An odd strand number "
               "(19) means one parallel β-pair, unlike every bacterial porin. "
               "Its role in the permeability transition pore is contested and "
               "is catalogued as an open question, not as a function."),
    ),
    CF(
        key="tmem175",
        name="TMEM175 endolysosomal K+/H+ channel",
        superfamily="tmem175", status=St.CHANNEL, fold=Fold.TMEM175,
        selectivity=Sel.K, gating=(G.PROTON, G.LIGAND_INTRACELLULAR),
        stoichiometry="dimer of two-repeat subunits", tm_per_subunit=12,
        signatures=(Sig("pfam", "PF06736", "TMEM175", L.FAMILY, 2, P.DB),),
        exemplars=(Ex("Hs_TMEM175", "TMEM175", "Homo sapiens", "Q9BSA9"),),
        human_genes=("TMEM175",),
        length_band_aa=(490, 520), filter_motif="none — no GYG",
        notes=("K+-selective without a selectivity filter: the conduction "
               "pathway is lined by isoleucines, not by the backbone carbonyl "
               "cage every other K+ channel uses. Proof that potassium "
               "selectivity evolved at least twice, and the reason the "
               "classifier's motif tier can never be its only K+ evidence. "
               "A Parkinson's-disease risk locus."),
    ),
]
