"""The hazard registry — every recorded way to classify an ion channel wrongly.

The parent project (`../ip3r_genes`) had exactly one of these. ITPR and RYR
share all four of their diagnostic Pfam domains, and when that project's
discovery scorer was first benchmarked it promoted **six out of six**
ryanodine-receptor decoys as candidate novel IP3 receptors. The fix — a
positive test on distances to a labelled sister-family bait — became its
decision D14 and is carried here as D14 too.

Across all ion channels the same shape of problem recurs at least twenty
times, so it is data rather than prose. Every hazard names the families it
confuses, the evidence that makes them look alike, the positive test that
separates them, and the module that owns that test. `scripts/s1_benchmark.py`
turns each row into a measurement: how many decoys does the discriminator
actually reject, and does it reject anything real?

Two rules govern this file.

1. **A hazard is only closed by a test, never by a caveat.** "Analysts should
   be careful of X" is not an entry; `discriminator` has to be something a
   script can run.
2. **The discriminator must be a positive test.** Never "it is not a Nav
   because its symbol says CACNA1C" — always "the four filter residues read
   E/E/E/E, which is Cav's signature and not Nav's". Names are evidence
   about annotators, not about proteins, and an unnamed gene model in a new
   genome has none.

Severity is the cost of getting it wrong at census scale: `high` means the
hazard silently changes the headline count.
"""

from __future__ import annotations

from .schema import Hazard
from .schema import Provenance as P

HAZARDS: list[Hazard] = [
    Hazard(
        "H1", "Nav, Cav, NALCN and CatSper share one domain architecture",
        ("nav", "cav", "nalcn", "catsper"),
        shared_evidence=("`PF00520`×4 in a single chain, 1,600–2,500 aa, same "
                         "24-TM pseudo-tetramer. Measured: CACNA1G (Cav3) "
                         "carries Ion_trans×4 and nothing else — no Cav-specific "
                         "domain at all."),
        discriminator=("Read the four selectivity-filter residues off the "
                       "pore-module alignment: DEKA = Nav, EEEE = Cav, "
                       "EEKE = NALCN. Confirm with reference identity margin "
                       "≥ 0.10 (D7)."),
        test_owner="src/classify/motifs.py:filter_signature()",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H2", "The Cys-loop ligand-binding domain exists without a channel",
        ("nachr", "nonchannel_achbp", "zac"),
        shared_evidence=("`PF02931` (Neur_chan_LBD) is a complete soluble "
                         "protein in the molluscan acetylcholine-binding "
                         "proteins."),
        discriminator=("Require `PF02932` (the TM region) as well as the LBD "
                       "for a channel. **There is no positive AChBP test at "
                       "the architecture tier (S2b, measured):** the absence "
                       "rule (LBD without `PF02932` ⇒ AChBP) called 6,633 "
                       "receptor-length proteins AChBP, and its best positive "
                       "replacement — a complete, soluble, AChBP-sized LBD — "
                       "agreed with the S3a profiles on 53 of 1,263 records, "
                       "~75 % of its calls in lineages with no known AChBP. "
                       "An LBD without the TM region is superfamily-only; "
                       "AChBP is called by a sequence-level tier (reference "
                       "margin, or the S3a profile margin — D32). ZACN, GLIC "
                       "and ELIC, which carry the LBD with no `PF02932`, fall "
                       "to the superfamily the same way."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H3", "The iGluR clamshell is the class C GPCR ligand-binding domain",
        ("ampa", "kainate", "nmda", "nonchannel_class_c_gpcr"),
        shared_evidence=("`PF01094` (ANF_receptor) — measured on GRM1 — is the "
                         "ancestral bacterial periplasmic binding protein, "
                         "shared by every mGluR, both GABA-B subunits, CaSR and "
                         "the sweet/umami taste receptors."),
        discriminator=("Require `PF00060` (Lig_chan, the pore region). "
                       "`PF01094` alone classifies as class C GPCR when "
                       "`PF00003` (7tm_3) is also present."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H4", "IP3 receptors and ryanodine receptors share every diagnostic domain",
        ("itpr", "ryr"),
        shared_evidence=("`PF08709`, `PF01365`, `PF08454`, `PF02815` are all "
                         "carried by both. The parent project measured 49 % of "
                         "zebrafish `PF08709` records to be RyRs."),
        discriminator=("`PF02026` (the RyR repeat, ×4) or `PF06459` with the "
                       "shared core ⇒ RYR (positive). **There is no positive "
                       "ITPR test at the architecture tier** — ITPR carries "
                       "no domain RyR lacks — so the core alone is "
                       "superfamily-only and the ITPR call comes from the "
                       "labelled-bait identity margin (D14a) or the S3a "
                       "profile margin (D32). The absence rule it replaced "
                       "(S2b) called N-terminal RyR fragments ITPR: 502 "
                       "records the IP3R project calls RYR. Length (2.7 kaa "
                       "vs 5.0 kaa) supports the call and never makes it."),
        test_owner="src/classify/rules.py + src/discovery/candidates.py",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H5", "Half the CLC family are antiporters, not channels",
        ("clc_channel", "clc_transporter", "clc_prokaryotic"),
        shared_evidence="One family, one fold, one Pfam model (`PF00654`).",
        discriminator=("Presence and context of the intracellular proton "
                       "glutamate (Glu_in); plus the tier-1 CLC tree, in which "
                       "the channel branch is nested inside a transporter "
                       "family. Reported as a mechanism call separate from the "
                       "family call (D24)."),
        test_owner="src/classify/motifs.py + src/phylo/forest.py",
        severity="medium", provenance=P.LIT,
    ),
    Hazard(
        "H6", "Anoctamin channels and scramblases are architecturally identical",
        ("ano_channel", "ano_scramblase"),
        shared_evidence=("Measured: ANO1 and ANO6 both carry exactly `PF04547` "
                         "+ `PF16178`, at similar lengths."),
        discriminator=("None at the domain level. Assignment is by reference "
                       "identity plus tier-1 clade membership, and the split is "
                       "flagged as literature-derived in the catalogue rather "
                       "than presented as measured."),
        test_owner="src/classify/reference.py",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H7", "Most TRP families carry no `PF00520` at all",
        ("trpc", "trpm", "trpml", "trpp", "trpv", "trpa"),
        shared_evidence=("Measured 2026-08-19: TRPV1, TRPV5, TRPA1 and TRPM2 "
                         "carry `PF00520`; TRPC3, TRPM8, MCOLN1/2 and PKD2 "
                         "carry none. Their pores are modelled by `PF08344`, "
                         "`PF18139`/`PF25508` and `PF08016` instead."),
        discriminator=("The census enumerates the P-loop superfamily from the "
                       "union of all seven pore models, never from `PF00520` "
                       "alone; recall per TRP family is a reported S1 metric."),
        test_owner="src/catalogue/registry.py:pore_signatures()",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H8", "Pfam's pore model does not track the number of transmembrane helices",
        ("kca_sk", "kir", "k2p", "kca_slo"),
        shared_evidence=("Measured: the 6TM SK channels (KCNN2) are annotated "
                         "with `PF07885` — the 2TM model — and KCNT1 likewise, "
                         "while KCNMA1 gets `PF00520`."),
        discriminator=("Topology is taken from UniProt TM features or a "
                       "hydropathy prediction, never inferred from which Pfam "
                       "pore model matched."),
        test_owner="src/classify/rules.py:topology_check()",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H9", "A voltage-sensor domain is not evidence of a channel",
        ("hv1", "nonchannel_vsp"),
        shared_evidence=("Measured: both HVCN1 (a channel with no pore domain) "
                         "and TPTE (a phosphatase) are annotated `PF00520`."),
        discriminator=("Presence of `PF10409`/PTEN phosphatase domain "
                       "classifies as VSP; `PF16799` classifies as Hv1. "
                       "Neither call may be made from `PF00520`. Measured "
                       "2026-08-19: `PF16799` is on human HVCN1 and **not** on "
                       "*Ciona* Hv1, so the positive Hv1 test is currently a "
                       "mammal-only instrument."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H10", "The MIR domain is shared with the O-mannosyltransferases",
        ("itpr", "ryr", "nonchannel_pomt"),
        shared_evidence="Measured: POMT1 carries `PF02815`, as ITPR and RYR do.",
        discriminator=("`PF02366` (PMT) present ⇒ mannosyltransferase. "
                       "Inherited from the parent project's S1 decoy panel."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="low", provenance=P.DB,
    ),
    Hazard(
        "H11", "CFTR and the sulfonylurea receptors are the same ABC architecture",
        ("cftr", "assoc_sur", "nonchannel_abc_transporter"),
        shared_evidence=("Measured: CFTR and ABCC8 both carry `PF00664`×2 + "
                         "`PF00005`×2. One is a chloride channel, the other "
                         "regulates a potassium channel it does not conduct "
                         "through."),
        discriminator=("`PF14396` (CFTR_R) present ⇒ CFTR (positive). **There "
                       "is no positive SUR test** (S2c): ABCC8/9 carry only "
                       "the generic ABC domains, and TMD0 (`PF24357`) is on "
                       "the long MRPs too. Without the R domain there is no "
                       "architecture-tier call. The absence rule it replaced "
                       "(no R domain ⇒ SUR) made 980 census calls, none a SUR: "
                       "918 bacterial ABC transporters with a cNMP and a C39 "
                       "peptidase domain, 62 eukaryotic fused gene models. "
                       "At the profile tier the generic ABC transporters "
                       "have their own decoy profile (S3a2), without which "
                       "the SUR profile won them."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H12", "The Kv tetramerisation domain is a generic BTB/POZ domain",
        ("kv_shaker", "kv_modifier", "nonchannel_kctd"),
        shared_evidence=("Measured: KCTD1 carries `PF02214`, the same model as "
                         "the Kv1 T1 domain. 25 human KCTD genes are named "
                         "after a channel domain and none is a channel."),
        discriminator=("KCTD is called on a KCTD C-terminal domain (S2c: "
                       "`PF31093`, `PF23110`, `PF31104`, `PF20871`, `PF31099` "
                       "— on 14,354 census records, 14,317 profile-KCTD). A "
                       "Kv needs a pore model (`PF00520`/`PF07885`) with T1. "
                       "T1 with neither is superfamily-only: the absence rule "
                       "it replaced (T1 without a pore ⇒ KCTD) called "
                       "N-terminal Kv gene models KCTD, and the shape test "
                       "(T1, 0 TM, complete) was contradicted by the S3a "
                       "profiles on 160 records."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H13", "`PF08016` covers TRPML, TRPP *and* polycystin-1",
        ("trpml", "trpp", "assoc_polycystin1"),
        shared_evidence=("Measured: MCOLN1, PKD2 and PKD1 all carry "
                         "`PF08016` (Polycystin cation channel)."),
        discriminator=("Each family by a domain it carries (S2b): `PF21381` "
                       "(MCLN_ECD) ⇒ TRPML; `PF18109` with `PF20519` ⇒ TRPP; "
                       "PLAT (`PF01477`), REJ (`PF02010`), GPS (`PF01825`) or "
                       "`PF00801`×≥5 with the channel domain ⇒ polycystin-1 "
                       "(channel-associated). The channel domain alone is "
                       "superfamily-only. The absence rule it replaced "
                       "(`PF20519` without the mucolipin domain ⇒ TRPP) was "
                       "measured in S3a calling 2,989 polycystin-1-like "
                       "proteins (median 2,263 aa) TRPP."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H14", "Ion selectivity is not recoverable from sequence family",
        ("gabaa", "nachr", "ampa", "trpv"),
        shared_evidence=("Cys-loop charge selectivity is set by a short ring at "
                         "the M1–M2 boundary and flips with three "
                         "substitutions; AMPA-receptor Ca2+ permeability is set "
                         "by RNA editing, which the genome does not show."),
        discriminator=("None. Selectivity is recorded as a *family-level "
                       "literature attribute*, never predicted per sequence, "
                       "and the classifier emits no selectivity call."),
        test_owner="(no test — a declared limit of the method)",
        severity="high", provenance=P.LIT,
    ),
    Hazard(
        "H15", "Gene-symbol prefixes group pores with their accessory subunits",
        ("assoc_kcne", "assoc_cav_aux", "assoc_nav_beta",
         "assoc_catsper_aux", "nonchannel_kctd"),
        shared_evidence=("`KCNE*`, `KCNMB*`, `CACNB*`, `CACNG*`, `CACNA2D*`, "
                         "`SCN*B`, `CATSPERB/G/D/E/Z` and `KCTD*` all sort "
                         "inside the channel symbol space."),
        discriminator=("Membership is decided by pore evidence only. Symbols "
                       "are never used to include or exclude — measured cost of "
                       "getting this wrong: +24 % on the human count (S20: 76 "
                       "auxiliary genes against 320 pore genes)."),
        test_owner="src/classify/classifier.py (symbol never consulted)",
        severity="high", provenance=P.DB,
    ),
    Hazard(
        "H16", "Kir and K2P look like one 2TM/1P group and are modelled differently",
        ("kir", "k2p"),
        shared_evidence=("Both are 2TM/1P potassium channels with the same "
                         "pore architecture, and the literature groups them "
                         "together."),
        discriminator=("Two `PF07885` per chain ⇒ K2P; `PF01007` ⇒ Kir. "
                       "Measured 2026-08-19: KCNK2 carries `PF07885`×2 and "
                       "**KCNJ2/KCNJ11 carry no `PF07885` at all** — so this "
                       "is not one model at two copy numbers, it is two "
                       "different models, and a rule written the first way "
                       "would have failed on every Kir."),
        test_owner="src/classify/rules.py:ARCHITECTURE_RULES",
        severity="low", provenance=P.DB,
    ),
    Hazard(
        "H17", "TMCO1 shares its only domain with the EMC3 insertase",
        ("tmco1", "nonchannel_emc3"),
        shared_evidence=("`PF01956` (EMC3/TMCO1-like) is the whole Pfam "
                         "architecture of both; human reviewed carriers are "
                         "exactly TMCO1 and EMC3 (measured 2026-09-29)."),
        discriminator=("None at the architecture tier: `PF01956` is "
                       "SHARED_WITH_DECOY and stops at the superfamily (D33). "
                       "The family call needs a sequence-level margin (S3a "
                       "profile or reference, D7/D32) against an EMC3 "
                       "profile — neither built yet."),
        test_owner="(open — needs profiles for tmco1 and nonchannel_emc3)",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H18", "TMEM87A shares its domains with TMEM87B and the GOST proteins",
        ("tmem87", "nonchannel_gost"),
        shared_evidence=("`PF06814` (GOST seven-TM) on TMEM87A/B, GPR107, "
                         "GPR108; `PF21901` (TMEM87 GOLD) on TMEM87A and "
                         "TMEM87B (measured 2026-09-29)."),
        discriminator=("None at the architecture tier; both signatures are "
                       "SHARED_WITH_DECOY. Family call by profile margin only "
                       "— and TMEM87A vs TMEM87B is a paralogue split no "
                       "domain can make."),
        test_owner="(open — needs profiles for tmem87 and nonchannel_gost)",
        severity="medium", provenance=P.DB,
    ),
    Hazard(
        "H19", "TMEM109 shares its only domain with BRI3BP",
        ("tmem109", "nonchannel_bri3bp"),
        shared_evidence=("`PF14965` is the whole Pfam architecture of both "
                         "human reviewed carriers, TMEM109 and BRI3BP "
                         "(measured 2026-09-29)."),
        discriminator=("None at the architecture tier; family call by "
                       "profile margin only."),
        test_owner="(open — needs profiles for tmem109 and nonchannel_bri3bp)",
        severity="low", provenance=P.DB,
    ),
    Hazard(
        "H20", "KChIPs are neuronal calcium sensors",
        ("assoc_kchip", "nonchannel_ncs"),
        shared_evidence=("KCNIP1–4 belong to the neuronal calcium sensor "
                         "family and carry its EF-hand models (`PF13499`, "
                         "`PF13833`), as do NCS1, hippocalcin, recoverin, the "
                         "VILIPs and the GCAPs. Measured in r4: the KChIP "
                         "profile alone called 22 human NCS/EF-hand proteins."),
        discriminator=("None at the architecture tier. Profile margin between "
                       "`assoc_kchip` and `nonchannel_ncs` (D32)."),
        test_owner="scripts/s3r4_compare.py (r4_reviewed_calls.tsv)",
        severity="low", provenance=P.DB,
    ),
]

HAZARD_BY_ID: dict[str, Hazard] = {h.hid: h for h in HAZARDS}
