"""Proposed channels added after S20, and the look-alikes they must be told from.

S20 compared the catalogue against three database channelomes (GtoPdb,
HGNC, UniProt KW-0407) and found eight human proteins that UniProt calls
ion channels and the catalogue did not name. They are added here, each as
its own superfamily — no relationship between them or to anything else is
claimed — with the evidence stated and every assertion `CURATED`:

* **PACC1** (TMEM206, PAC / ASOR) — the proton-activated chloride channel;
  solved structures and knockout phenotypes. `CHANNEL`.
* **TMCO1** (CLAC), **TMEM87A** (Elkin1 / GolpHCat), **TMEM109**
  (mitsugumin-23), **CLCC1** (MCLC / ERAC1), **CCDC51** (MITOK) and
  **GPHRA/GPHRB** (GPR89A/B) — each reported to conduct ions, none with the
  weight of evidence the census families carry. `CHANNEL_CONTESTED`, the
  reason attached, as the scope document requires of a live dispute.

**Search space (D34, D43).** Declared first with `enumerate=False` (D42),
then brought into the census by revision r4 (D43): census v2's 67
signatures plus these eight, walked as a delta, with the old records and
their calls kept unchanged.

**Look-alikes.** Three of the eight share their only domain with a protein
nobody calls a channel, so the domain is `SHARED_WITH_DECOY` and the decoy
is catalogued here beside it (hazards H17–H19): TMCO1 with EMC3 (an ER
membrane-complex insertase subunit), TMEM87A with TMEM87B and the GOST
proteins GPR107/GPR108, TMEM109 with BRI3BP. CCDC51 carries no Pfam domain
at all (measured 2026-09-29), so like Hv1 outside mammals it is reachable
only by profile.
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

_NOTE = "added after S20; enumerated from census v2 r4 (D43)"


def _sig(acc: str, name: str, level: L = L.FAMILY, note: str = _NOTE,
         enumerate: bool = True) -> Sig:
    """Census-family signatures are enumerated from r4 on (D43) — including
    the SHARED_WITH_DECOY ones, because for TMCO1, TMEM87A and TMEM109 the
    shared domain is the only one there is; the decoy then arrives in the
    census too and is told apart by its own profile."""
    return Sig("pfam", acc, name, level, 1, P.DB, note, enumerate=enumerate)


SUPERFAMILIES: list[SF] = [
    SF("pac", "Proton-activated chloride channels", Fold.PAC, True, "full-length",
       (_sig("PF15122", "TMEM206", L.SUPERFAMILY),)),
    SF("tmco1", "TMCO1 / EMC3 membrane-insertase fold", Fold.OXA1_LIKE, True,
       "full-length", (_sig("PF01956", "EMC3_TMCO1", L.SHARED_WITH_DECOY),),
       notes="Oxa1/YidC insertase relatives; only TMCO1 is reported to conduct."),
    SF("tmem87", "TMEM87 / GOST seven-TM proteins", Fold.GOST, True, "full-length",
       (_sig("PF06814", "GOST_TM", L.SHARED_WITH_DECOY),)),
    SF("tmem109", "TMEM109 / BRI3BP", Fold.TMEM109, True, "full-length",
       (_sig("PF14965", "BRI3BP", L.SHARED_WITH_DECOY),)),
    SF("mclc", "Mid-1-related chloride channels", Fold.MCLC, True, "full-length",
       (_sig("PF05934", "MCLC", L.SUPERFAMILY),)),
    SF("mitok", "Mitochondrial K+ channel MITOK", Fold.UNKNOWN, True, "full-length", ()),
    SF("gphr", "Golgi pH regulators", Fold.GPHR, True, "full-length",
       (_sig("PF12537", "GPHR_N", L.SUPERFAMILY),)),
]

FAMILIES: list[CF] = [
    CF(key="pacc", name="Proton-activated chloride channel (PACC1 / TMEM206)",
       superfamily="pac", status=St.CHANNEL, fold=Fold.PAC,
       selectivity=Sel.ANION, gating=(G.PROTON,), stoichiometry="trimer",
       tm_per_subunit=2, signatures=(_sig("PF15122", "TMEM206"),),
       exemplars=(Ex("Hs_PACC1", "PACC1", "Homo sapiens", "Q9H813"),),
       human_genes=("PACC1",), length_band_aa=(300, 400),
       notes=("Opens at acidic extracellular pH; mediates acid-induced cell "
              "death. On UniProt KW-0407, not on GtoPdb or HGNC 177 (S20)."),
       provenance=P.CURATED),
    CF(key="tmco1", name="Calcium load-activated calcium channel (TMCO1)",
       superfamily="tmco1", status=St.CHANNEL_CONTESTED, fold=Fold.OXA1_LIKE,
       selectivity=Sel.CA, gating=(G.CALCIUM,), stoichiometry="unknown",
       signatures=(_sig("PF01956", "EMC3_TMCO1", L.SHARED_WITH_DECOY),),
       exemplars=(Ex("Hs_TMCO1", "TMCO1", "Homo sapiens", "Q9UM00"),
                  Ex("Mm_Tmco1", "Tmco1", "Mus musculus", "Q921L3"),
                  Ex("Dr_tmco1", "tmco1", "Danio rerio", "Q6DGW9"),
                  Ex("Dd_tmco1", "tmco1", "Dictyostelium discoideum", "Q54TU8")),
       human_genes=("TMCO1",), length_band_aa=(180, 280), confusable_with=("H17",),
       notes=("Reported as an ER Ca2+ leak channel that tetramerises on "
              "store overload; also a subunit of the GEL insertion complex. "
              "Contested: an insertase-fold protein with a second, "
              "better-supported job."),
       provenance=P.CURATED),
    CF(key="tmem87", name="TMEM87 Golgi cation channels (TMEM87A, TMEM87B), by descent",
       superfamily="tmem87", status=St.CHANNEL_CONTESTED, fold=Fold.GOST,
       selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL, G.VOLTAGE),
       stoichiometry="unknown", tm_per_subunit=7,
       signatures=(_sig("PF06814", "GOST_TM", L.SHARED_WITH_DECOY),
                   _sig("PF21901", "TMEM87A-B_GOLD", L.FAMILY,
                        "TMEM87 GOLD domain: all 1,660 UniProt carriers are "
                        "metazoan TMEM87s, none a GOST protein (measured "
                        "2026-09-30) — the family's positive test (H18)")),
       exemplars=(Ex("Hs_TMEM87A", "TMEM87A", "Homo sapiens", "Q8NBN3"),
                  Ex("Mm_Tmem87a", "Tmem87a", "Mus musculus", "Q8BXN9"),
                  Ex("Hs_TMEM87B", "TMEM87B", "Homo sapiens", "Q96K49"),
                  Ex("Mm_Tmem87b", "Tmem87b", "Mus musculus", "Q8BKU8")),
       human_genes=("TMEM87A", "TMEM87B"), length_band_aa=(480, 620),
       confusable_with=("H18",),
       notes=("Family by descent (D24, D45): TMEM87A and TMEM87B are "
              "vertebrate paralogues and an invertebrate's single TMEM87 is "
              "co-orthologous to both, so the family holds all three. "
              "Mechanism is a separate, literature field: channel activity "
              "is reported for TMEM87A only — a mechanosensitive cation "
              "channel (Elkin1) and a voltage-gated Golgi pH-regulating "
              "channel (GolpHCat); TMEM87B has none reported, and is counted "
              "by descent, not because it is shown to conduct."),
       provenance=P.CURATED),
    CF(key="tmem109", name="Mitsugumin-23 cation channel (TMEM109)",
       superfamily="tmem109", status=St.CHANNEL_CONTESTED, fold=Fold.TMEM109,
       selectivity=Sel.CATION_NONSELECTIVE, gating=(G.VOLTAGE,),
       stoichiometry="unknown",
       signatures=(_sig("PF14965", "BRI3BP", L.SHARED_WITH_DECOY),),
       exemplars=(Ex("Hs_TMEM109", "TMEM109", "Homo sapiens", "Q9BVC6"),
                  Ex("Mm_Tmem109", "Tmem109", "Mus musculus", "Q3UBX0")),
       human_genes=("TMEM109",), length_band_aa=(200, 290), confusable_with=("H19",),
       notes="ER/SR membrane protein reported to form a voltage-gated cation pore.",
       provenance=P.CURATED),
    CF(key="clcc1", name="ER anion channel CLCC1 (MCLC / ERAC1)",
       superfamily="mclc", status=St.CHANNEL_CONTESTED, fold=Fold.MCLC,
       selectivity=Sel.ANION, gating=(G.UNKNOWN,), stoichiometry="unknown",
       signatures=(_sig("PF05934", "MCLC"),),
       exemplars=(Ex("Hs_CLCC1", "CLCC1", "Homo sapiens", "Q96S66"),),
       human_genes=("CLCC1",), length_band_aa=(480, 620),
       notes=("Named by analogy with yeast Mid-1; reported as an ER Cl- "
              "channel. Not a CLC and not a CLIC, despite both names."),
       provenance=P.CURATED),
    CF(key="mitok", name="Mitochondrial ATP-sensitive K+ channel MITOK (CCDC51)",
       superfamily="mitok", status=St.CHANNEL_CONTESTED, fold=Fold.UNKNOWN,
       selectivity=Sel.K, gating=(G.LIGAND_INTRACELLULAR,), stoichiometry="unknown",
       signatures=(),
       exemplars=(Ex("Hs_CCDC51", "CCDC51", "Homo sapiens", "Q96ER9"),),
       human_genes=("CCDC51",), length_band_aa=(350, 470),
       notes=("Proposed as the pore of mitoKATP with MITOSUR (ABCB8). Carries "
              "no Pfam domain (measured 2026-09-29): no domain search can "
              "find it, only a profile."),
       provenance=P.CURATED),
    CF(key="gphr", name="Golgi pH regulator anion channels (GPHRA / GPHRB, GPR89)",
       superfamily="gphr", status=St.CHANNEL_CONTESTED, fold=Fold.GPHR,
       selectivity=Sel.ANION, gating=(G.VOLTAGE,), stoichiometry="unknown",
       tm_per_subunit=9,
       signatures=(_sig("PF12537", "GPHR_N"),
                   _sig("PF12430", "ABA_GPCR", note="Pfam names the plant "
                        "homologues (GTG1/2) abscisic-acid GPCRs — the name is "
                        "a claim about the plant proteins, not a GPCR fold")),
       exemplars=(Ex("Hs_GPHRA", "GPHRA", "Homo sapiens", "B7ZAQ6"),
                  Ex("Hs_GPHRB", "GPHRB", "Homo sapiens", "P0CG08")),
       human_genes=("GPHRA", "GPHRB"), length_band_aa=(400, 520),
       notes=("Reported as a voltage-dependent Golgi anion channel that "
              "counters luminal acidification. GPHRA and GPHRB are "
              "near-identical duplicates."),
       provenance=P.CURATED),
    # ------------------------------------------------ the look-alikes
    CF(key="nonchannel_emc3", name="ER membrane complex subunit EMC3",
       superfamily="tmco1", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.OXA1_LIKE,
       selectivity=Sel.UNKNOWN, gating=(), stoichiometry="1 per EMC",
       signatures=(_sig("PF01956", "EMC3_TMCO1", L.SHARED_WITH_DECOY),),
       exemplars=(Ex("Hs_EMC3", "EMC3", "Homo sapiens", "Q9P0I2"),
                  Ex("Mm_Emc3", "Emc3", "Mus musculus", "Q99KI3"),
                  Ex("Dr_emc3", "emc3", "Danio rerio", "Q7SXW4"),
                  Ex("Sc_EMC3", "EMC3", "Saccharomyces cerevisiae", "P36039"),
                  Ex("Dd_emc3", "emc3", "Dictyostelium discoideum", "Q54YN3")),
       human_genes=("EMC3",), length_band_aa=(220, 300), confusable_with=("H17",),
       notes="The insertase core of the ER membrane complex; carries TMCO1's only domain.",
       provenance=P.CURATED),
    CF(key="nonchannel_gost", name="GOST seven-TM proteins GPR107 and GPR108",
       superfamily="tmem87", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.GOST,
       selectivity=Sel.UNKNOWN, gating=(),
       signatures=(_sig("PF06814", "GOST_TM", L.SHARED_WITH_DECOY),),
       exemplars=(Ex("Hs_GPR107", "GPR107", "Homo sapiens", "Q5VW38"),
                  Ex("Mm_Gpr107", "Gpr107", "Mus musculus", "Q8BUV8"),
                  Ex("Mm_Gpr108", "Gpr108", "Mus musculus", "Q91WD0"),
                  Ex("At_CAND6", "CAND6", "Arabidopsis thaliana", "Q9LZ39",
                     "plant GOST protein; no TMEM87 GOLD domain"),
                  Ex("Sc_PTM1", "PTM1", "Saccharomyces cerevisiae", "P32857",
                     "yeast GOST protein; no TMEM87 GOLD domain")),
       human_genes=("GPR107", "GPR108"), length_band_aa=(480, 650),
       confusable_with=("H18",),
       notes=("Golgi seven-TM proteins; GPR107/108 are not GPCRs despite the "
              "symbol. They — and the fungal and plant GOST proteins, seeded "
              "here after TMEM87-by-descent drew them in (D45) — lack TMEM87's "
              "GOLD domain (PF21901), whose 1,660 carriers are all metazoan."),
       provenance=P.CURATED),
    CF(key="nonchannel_bri3bp", name="BRI3-binding protein (BRI3BP)",
       superfamily="tmem109", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.TMEM109,
       selectivity=Sel.UNKNOWN, gating=(),
       signatures=(_sig("PF14965", "BRI3BP", L.SHARED_WITH_DECOY),),
       exemplars=(Ex("Hs_BRI3BP", "BRI3BP", "Homo sapiens", "Q8WY22"),),
       human_genes=("BRI3BP",), length_band_aa=(210, 290), confusable_with=("H19",),
       notes="Carries TMEM109's only domain; no channel activity reported.",
       provenance=P.CURATED),
    CF(key="nonchannel_ncs",
       name="Neuronal calcium sensors (NCS1, hippocalcin, recoverin, VILIPs, GCAPs)",
       superfamily="ploop", status=St.NON_CHANNEL_HOMOLOG, fold=Fold.UNKNOWN,
       selectivity=Sel.UNKNOWN, gating=(),
       signatures=(Sig("pfam", "PF13499", "EF-hand_7", L.SHARED_WITH_DECOY, 1, P.DB),
                   Sig("pfam", "PF00036", "EF-hand_1", L.SHARED_WITH_DECOY, 1, P.DB,
                       "measured on NCS1, HPCA, GUCA1A (S0, 2026-09-29); the "
                       "KChIPs carry PF13833 instead")),
       exemplars=(Ex("Hs_NCS1", "NCS1", "Homo sapiens", "P62166"),
                  Ex("Hs_HPCA", "HPCA", "Homo sapiens", "P84074"),
                  Ex("Hs_RCVRN", "RCVRN", "Homo sapiens", "P35243"),
                  Ex("Hs_GUCA1A", "GUCA1A", "Homo sapiens", "P43080")),
       human_genes=("NCS1", "HPCA", "HPCAL1", "HPCAL4", "NCALD", "VSNL1",
                    "RCVRN", "GUCA1A", "GUCA1B", "GUCA1C"),
       length_band_aa=(180, 220), confusable_with=("H20",),
       notes=("The KChIPs' own family: EF-hand calcium sensors with no channel "
              "role. Added in r4 after the one-family KChIP profile called "
              "recoverin, hippocalcin, NCS-1 and the GCAPs at high confidence — "
              "a profile is evidence only when every family that could score "
              "has one (S2c)."),
       provenance=P.CURATED),
]
