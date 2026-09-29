# S20 — the auxiliary subunits and the published channelome

Rendered by `scripts/s20_report.py` from the tables in this directory (D13). Raw lists, HGNC lookups and phmmer inputs: `<data root>/raw_api/s20/`.

![S20](figures/auxiliary.png)

## 1. Three database channelomes, one denominator

Three human ion-channel lists maintained by databases, each archived with its release under `<data root>/raw_api/s20/`: **GtoPdb** (GtoPdb Version: 2026.3 - published: 2026-09-16; target types vgic, lgic, other_ic), **HGNC** (gene group 177 'Ion channels' and every group below it; HGNC group 177 branch, fetched Tue, 29 Sep 2026 23:24:34 GMT) and **UniProt** (reviewed human entries with keyword KW-0407 'Ion channel', UniProtKB 2026_03 — the release census v2 was enumerated from). Every gene is joined on its HGNC id; the catalogue's own symbols resolved 468/468 (renamed: AQP0 → MIP, TMEM249 → CATSPERQ).

**The three lists hold 285, 331 and 338 genes, and 400 together** — the published range of the human channelome, reproduced by three curated databases. Of the union, **314 are this catalogue's pore-forming census genes** (320 in the catalogue; 6 are on no list) and **86 are not**:

| list | total | pore (census) | auxiliary | auxiliary uncatalogued | out of scope | out of scope uncatalogued | transporter | transporter or enzyme | pore candidate | paracellular | pseudogene | census pore genes missed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GtoPdb | 285 | 263 | 4 | 0 | 12 | 0 | 5 | 0 | 0 | 0 | 1 | 57 |
| HGNC | 331 | 290 | 20 | 0 | 13 | 1 | 5 | 0 | 0 | 0 | 2 | 30 |
| UniProt KW-0407 | 338 | 281 | 31 | 5 | 0 | 0 | 2 | 8 | 8 | 2 | 1 | 39 |
| union | 400 | 314 | 39 | 5 | 13 | 1 | 7 | 8 | 8 | 2 | 3 | 6 |

Column meanings — *auxiliary*: a catalogued `channel_associated` subunit; *auxiliary uncatalogued*, *pore candidate*, *transporter or enzyme*, *paracellular*, *pseudogene*, *out of scope uncatalogued*: the hand-curated classes of §4; *out of scope*: catalogued aquaporins; *transporter*: catalogued CLC and SLC26 transporters.

## 2. Where the spread comes from

**238 genes are on all three lists, and all 238 of them are pore-forming census genes** — the agreed core is pure. The lists differ on two independent axes, and both are scope, not biology:

* **What they add beyond the pore-forming genes**: GtoPdb 22; HGNC 41; UniProt KW-0407 57 — auxiliary subunits, aquaporins, CLC transporters, and (UniProt only) enzymes, transporters and claudins carrying the 'Ion channel' keyword.
* **Which contested or large-pore families they leave out**: GtoPdb misses 57 census genes; HGNC misses 30 census genes; UniProt KW-0407 misses 39 census genes.

| family | GtoPdb missed | HGNC missed | UniProt KW-0407 missed | genes |
|---|---|---|---|---|
| tmc | 8 | 8 | 5 | TMC1, TMC2, TMC3, TMC4, TMC5, TMC6, TMC7, TMC8 |
| connexin | 0 | 0 | 21 | GJA1, GJA10, GJA3, GJA4, GJA5, GJA8, GJA9, GJB1, GJB2, GJB3, GJB4, GJB5, GJB6, GJB7, GJC1, GJC2, GJC3, GJD2, GJD3, GJD4, GJE1 |
| ano_scramblase | 8 | 0 | 7 | ANO10, ANO3, ANO4, ANO5, ANO6, ANO7, ANO8, ANO9 |
| clic | 6 | 6 | 0 | CLIC1, CLIC2, CLIC3, CLIC4, CLIC5, CLIC6 |
| vdac | 3 | 0 | 3 | VDAC1, VDAC2, VDAC3 |
| tweety | 3 | 3 | 0 | TTYH1, TTYH2, TTYH3 |
| osca_tmem63 | 3 | 3 | 0 | TMEM63A, TMEM63B, TMEM63C |
| otop | 3 | 3 | 0 | OTOP1, OTOP2, OTOP3 |
| calhm | 6 | 0 | 0 | CALHM1, CALHM2, CALHM3, CALHM4, CALHM5, CALHM6 |
| lrrc8 | 5 | 0 | 0 | LRRC8A, LRRC8B, LRRC8C, LRRC8D, LRRC8E |
| bestrophin | 4 | 0 | 0 | BEST1, BEST2, BEST3, BEST4 |
| tric | 2 | 2 | 0 | TMEM38A, TMEM38B |
| nav | 1 | 1 | 1 | SCN7A |
| mcu | 0 | 2 | 1 | MCU, MCUB |
| asic | 2 | 0 | 0 | ASIC4, ASIC5 |
| tmem175 | 1 | 1 | 0 | TMEM175 |
| ano_channel | 1 | 0 | 0 | ANO2 |
| trpm | 0 | 1 | 0 | TRPM8 |
| glyr | 0 | 0 | 1 | GLRA4 |
| kir | 1 | 0 | 0 | KCNJ18 |

Two of these are single annotation choices with large effects: **UniProt's KW-0407 is on none of the 21 connexins** (their keyword is 'Gap junction'), and **GtoPdb lists no TMEM16 scramblase, TMC, CLIC, tweety, OSCA or otopetrin** — the contested families of the scope document.

## 3. The auxiliary subunits — how many, and who counts them

**The catalogue names 71 human auxiliary genes in 11 families; the lists add 5 it does not** (KCNIP1, KCNIP2, KCNIP3, KCNIP4, TMEM37). Counted as channels, all 76 would inflate the 320 pore-forming genes by **24%** — not the 'roughly 15 %' the scope document stated before this measurement.

**How often the lists count them**: GtoPdb 4 (1.4% of its total); HGNC 20 (6.0% of its total); UniProt KW-0407 36 (10.7% of its total); union 44 (11.0% of its total). Each list counts a different set:

| auxiliary family | genes | GtoPdb | HGNC | UniProt |
|---|---|---|---|---|
| assoc_catsper_aux | 6 | 0 | 0 | 0 |
| assoc_cav_aux | 16 | 0 | 16 | 15 |
| assoc_clc_aux | 5 | 0 | 0 | 0 |
| assoc_iglur_aux | 9 | 0 | 0 | 0 |
| assoc_k_beta | 13 | 0 | 0 | 9 |
| assoc_kcne | 5 | 0 | 0 | 3 |
| assoc_mcu_reg | 5 | 4 | 0 | 0 |
| assoc_nav_beta | 4 | 0 | 4 | 1 |
| assoc_polycystin1 | 4 | 0 | 0 | 3 |
| assoc_stim | 2 | 0 | 0 | 0 |
| assoc_sur | 2 | 0 | 0 | 0 |

## 4. Genes a list carries and the catalogue does not name

**27 genes**, each classified by hand (`s20_curate.CURATED`, provenance CURATED — read from the UniProt record name or HGNC locus type, not yet checked against primary literature):

| class | genes | meaning |
|---|---|---|
| pseudogene | 3 | not a protein-coding gene (HGNC locus type) or a UniProt 'putative' pseudogene product |
| auxiliary_uncatalogued | 5 | an auxiliary subunit of a catalogued channel that the catalogue does not list |
| pore_candidate | 8 | a proposed pore-forming channel not in the catalogue |
| transporter_or_enzyme | 8 | a transporter or enzyme with a reported channel-like activity |
| paracellular | 2 | forms paracellular (tight-junction) pores, not a transmembrane channel |
| out_of_scope_uncatalogued | 1 | out of scope by D23 (water channel) and not listed in the catalogue |

| gene | lists | class | reason |
|---|---|---|---|
| KCNIP1 | UniProt KW-0407 | auxiliary_uncatalogued | Kv channel-interacting protein (KChIP), Kv4 auxiliary |
| KCNIP2 | UniProt KW-0407 | auxiliary_uncatalogued | KChIP2, Kv4 auxiliary |
| KCNIP3 | UniProt KW-0407 | auxiliary_uncatalogued | KChIP3, Kv4 auxiliary |
| KCNIP4 | UniProt KW-0407 | auxiliary_uncatalogued | KChIP4, Kv4 auxiliary |
| TMEM37 | UniProt KW-0407 | auxiliary_uncatalogued | voltage-dependent calcium channel gamma-like subunit |
| AQP12B | HGNC | out_of_scope_uncatalogued | aquaporin-12B; the catalogue lists AQP12A only |
| CLDN17 | UniProt KW-0407 | paracellular | claudin-17, tight-junction anion pore |
| CLDN4 | UniProt KW-0407 | paracellular | claudin-4, tight-junction pore |
| CCDC51 | UniProt KW-0407 | pore_candidate | mitochondrial potassium channel (MITOK) |
| CLCC1 | UniProt KW-0407 | pore_candidate | ER anion channel 1 (chloride channel CLIC-like 1) |
| GPHRA | UniProt KW-0407 | pore_candidate | Golgi pH regulator A (GPR89A), reported anion channel |
| GPHRB | UniProt KW-0407 | pore_candidate | Golgi pH regulator B (GPR89B), reported anion channel |
| PACC1 | UniProt KW-0407 | pore_candidate | proton-activated chloride channel (PAC / ASOR, TMEM206) |
| TMCO1 | UniProt KW-0407 | pore_candidate | calcium load-activated calcium channel (CLAC) |
| TMEM109 | UniProt KW-0407 | pore_candidate | voltage-gated cation channel TMEM109 (mitsugumin-23) |
| TMEM87A | UniProt KW-0407 | pore_candidate | Golgi-pH regulating cation channel (Elkin1) |
| FXYD6P3 | UniProt KW-0407 | pseudogene | UniProt 'putative FXYD domain-containing ion transport regulator 8' |
| GJA6P | HGNC | pseudogene | HGNC locus type pseudogene |
| TRPC2 | GtoPdb, HGNC | pseudogene | HGNC locus type pseudogene; functional in rodents |
| CYBB | UniProt KW-0407 | transporter_or_enzyme | NADPH oxidase 2 (gp91phox), H+ conductance |
| MFSD8 | UniProt KW-0407 | transporter_or_enzyme | MFS transporter CLN7, reported lysosomal Cl- channel |
| NOX5 | UniProt KW-0407 | transporter_or_enzyme | NADPH oxidase 5, H+ conductance |
| SLC17A6 | UniProt KW-0407 | transporter_or_enzyme | vesicular glutamate transporter 2, Cl- conductance |
| SLC17A7 | UniProt KW-0407 | transporter_or_enzyme | vesicular glutamate transporter 1, Cl- conductance |
| SLC17A8 | UniProt KW-0407 | transporter_or_enzyme | vesicular glutamate transporter 3, Cl- conductance |
| STING1 | UniProt KW-0407 | transporter_or_enzyme | innate-immunity adaptor STING, reported H+ channel |
| UCP1 | UniProt KW-0407 | transporter_or_enzyme | mitochondrial carrier SLC25A7, H+ leak |

The two catalogue gaps are emergent rows, not edits made here: the `auxiliary_uncatalogued` genes (KChIP1–4, TMEM37) belong in `channel_associated` families, and the `pore_candidate` genes — PACC1 above all, a proton-activated Cl⁻ channel with solved structures — are candidate census families, which is a change to the census search space (D34).

## 5. The auxiliary families are grouped by partner, not by descent

All-against-all `phmmer` among each family's human members (E ≤ 1e-3 = homologous; groups are connected components): **6 of 11 auxiliary families pool unrelated proteins**, 25 homology groups in all. A profile built across a pooled family is not a detector of any of its parts (S3a: the LOO decoys got no hit), so the panel census below counts homology groups, not families.

| family | groups | homology groups (human genes) |
|---|---|---|
| assoc_catsper_aux | 4 | CATSPERD, CATSPERE · CATSPERG · CATSPERB · CATSPERZ (not tested: TMEM249) |
| assoc_cav_aux | 3 | CACNB1, CACNB2, CACNB3, CACNB4 · CACNG1, CACNG2, CACNG3, CACNG4, CACNG5, CACNG6, CACNG7, CACNG8 · CACNA2D1, CACNA2D2, CACNA2D3, CACNA2D4 |
| assoc_clc_aux | 2 | CLCA1, CLCA2, CLCA4 · BSND (not tested: CLCA3P) |
| assoc_iglur_aux | 4 | SHISA6, SHISA7, SHISA8, SHISA9 · CNIH2, CNIH3 · GSG1L · NETO1, NETO2 |
| assoc_k_beta | 4 | KCNAB1, KCNAB2, KCNAB3 · DPP10, DPP6 · KCNMB1, KCNMB2, KCNMB3, KCNMB4 · LRRC26, LRRC38, LRRC52, LRRC55 |
| assoc_kcne | 1 | KCNE1, KCNE2, KCNE3, KCNE4, KCNE5 |
| assoc_mcu_reg | 3 | MICU1, MICU2, MICU3 · MCUR1 · SMDT1 |
| assoc_nav_beta | 1 | SCN1B, SCN2B, SCN3B, SCN4B |
| assoc_polycystin1 | 1 | PKD1, PKD1L1, PKD1L2, PKD1L3 |
| assoc_stim | 1 | STIM1, STIM2 |
| assoc_sur | 1 | ABCC8, ABCC9 |

## 6. Auxiliaries across the panel

S3b's panel profile calls to an auxiliary family: **3053**. Each was searched against the whole human proteome; it counts as an auxiliary of a group only if its **best human hit is a human member of that group** (reciprocal best hit, E ≤ 1e-5). **1103 pass**; 1949 have a better human hit outside the group, 1 no human hit.

| homology group | human genes | calls | pass (high) | pass (medium) | not auxiliary | species (high) | panel groups (high) |
|---|---|---|---|---|---|---|---|
| assoc_catsper_aux.1 | CATSPERD, CATSPERE | 38 | 12 | 6 | 0 | 7 | vertebrate |
| assoc_catsper_aux.2 | CATSPERG |  | 11 | 0 |  | 8 | vertebrate |
| assoc_catsper_aux.3 | CATSPERB |  | 8 | 1 |  | 7 | vertebrate |
| assoc_catsper_aux.4 | CATSPERZ |  | 0 | 0 |  | 0 |  |
| assoc_cav_aux.1 | CACNB1, CACNB2, CACNB3, CACNB4 | 265 | 59 | 8 | 0 | 19 | deuterostome, invertebrate, vertebrate |
| assoc_cav_aux.2 | CACNG1, CACNG2, CACNG3, CACNG4, CACNG5, CACNG6, CACNG7, CACNG8 |  | 63 | 48 |  | 13 | vertebrate |
| assoc_cav_aux.3 | CACNA2D1, CACNA2D2, CACNA2D3, CACNA2D4 |  | 77 | 10 |  | 23 | cnidarian, deuterostome, invertebrate, vertebrate |
| assoc_clc_aux.1 | CLCA1, CLCA2, CLCA4 | 118 | 108 | 9 | 1 | 20 | basal_metazoan, deuterostome, invertebrate, vertebrate |
| assoc_clc_aux.2 | BSND |  | 0 | 0 |  | 0 |  |
| assoc_iglur_aux.1 | SHISA6, SHISA7, SHISA8, SHISA9 | 138 | 48 | 11 | 24 | 13 | vertebrate |
| assoc_iglur_aux.2 | CNIH2, CNIH3 |  | 0 | 0 |  | 0 |  |
| assoc_iglur_aux.3 | GSG1L |  | 23 | 12 |  | 13 | vertebrate |
| assoc_iglur_aux.4 | NETO1, NETO2 |  | 20 | 0 |  | 10 | vertebrate |
| assoc_k_beta.1 | KCNAB1, KCNAB2, KCNAB3 | 1574 | 60 | 3 | 1397 | 27 | basal_metazoan, ciliate, deuterostome, holozoa, invertebrate, plant, prokaryote, vertebrate |
| assoc_k_beta.2 | DPP10, DPP6 |  | 39 | 2 |  | 18 | deuterostome, invertebrate, vertebrate |
| assoc_k_beta.3 | KCNMB1, KCNMB2, KCNMB3, KCNMB4 |  | 36 | 5 |  | 13 | vertebrate |
| assoc_k_beta.4 | LRRC26, LRRC38, LRRC52, LRRC55 |  | 32 | 0 |  | 12 | vertebrate |
| assoc_kcne.1 | KCNE1, KCNE2, KCNE3, KCNE4, KCNE5 | 49 | 43 | 6 | 0 | 11 | vertebrate |
| assoc_mcu_reg.1 | MICU1, MICU2, MICU3 | 148 | 77 | 1 | 53 | 34 | algae, amoebozoa, basal_metazoan, cnidarian, deuterostome, excavate, holozoa, invertebrate, plant, vertebrate |
| assoc_mcu_reg.2 | MCUR1 |  | 4 | 13 |  | 4 | vertebrate |
| assoc_mcu_reg.3 | SMDT1 |  | 0 | 0 |  | 0 |  |
| assoc_nav_beta.1 | SCN1B, SCN2B, SCN3B, SCN4B | 448 | 59 | 1 | 388 | 13 | vertebrate |
| assoc_polycystin1.1 | PKD1, PKD1L1, PKD1L2, PKD1L3 | 199 | 89 | 23 | 87 | 21 | basal_metazoan, cnidarian, deuterostome, invertebrate, vertebrate |
| assoc_stim.1 | STIM1, STIM2 | 45 | 32 | 13 | 0 | 15 | deuterostome, invertebrate, vertebrate |
| assoc_sur.1 | ABCC8, ABCC9 | 31 | 26 | 5 | 0 | 13 | vertebrate |

**Vertebrate-only in the panel (13 groups)**: CATSPERD, CATSPERG, CATSPERB, CACNG1, SHISA6, GSG1L, NETO1, KCNMB1, LRRC26, KCNE1, MCUR1, SCN1B, ABCC8. **Wider**: CACNB1 (deuterostome, invertebrate, vertebrate), CACNA2D1 (cnidarian, deuterostome, invertebrate, vertebrate), CLCA1 (basal_metazoan, deuterostome, invertebrate, vertebrate), KCNAB1 (basal_metazoan, ciliate, deuterostome, holozoa, invertebrate, plant, prokaryote, vertebrate), DPP10 (deuterostome, invertebrate, vertebrate), MICU1 (algae, amoebozoa, basal_metazoan, cnidarian, deuterostome, excavate, holozoa, invertebrate, plant, vertebrate), PKD1 (basal_metazoan, cnidarian, deuterostome, invertebrate, vertebrate), STIM1 (deuterostome, invertebrate, vertebrate).

Two readings these rows do **not** support. *Vertebrate-only* is where a human-seeded instrument found the group, not an absence claim: no zero cell here has passed D37's genome test (CatSper auxiliaries, for one, are reported outside vertebrates *(pending: S10)*). And a reciprocal best human hit shows **homology, not function**: the Kvβ group's plant and prokaryotic members are aldo-keto reductases whose closest human relative happens to be Kvβ, and nothing here says they serve a channel.

**No call at all for 4 groups** (CATSPERZ, BSND, CNIH2,CNIH3, SMDT1): the pooled family profile does not reach them. That is a detection limit of a profile built across unrelated proteins (§5), not an absence — including in human, where these genes exist.


The calls that fail are the shared-domain drift S3b found for TRPN and LRRC8: their best human hits are SLIT1 223, SLIT3 163, SLIT2 99, HMCN1 96, TTN 76, LOXHD1 63, HMCN2 43, BTNL2 34 ….
