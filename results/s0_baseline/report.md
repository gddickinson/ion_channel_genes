# S0 — catalogue verification baseline

*Generated 2026-09-30 09:37:58 by `scripts/s0_report.py` from the tables in `results/s0_baseline`. Nothing here is hand-written (D13). 901 live requests, 0 failure(s), 569.0 s.*

## What the catalogue claims

| quantity | n |
|---|---|
| superfamilies | 32 |
| families total | 103 |
| families census | 75 |
| families control | 28 |
| human genes census | 328 |
| human genes catalogued | 496 |
| exemplars | 195 |
| exemplars with accession | 195 |
| signatures | 132 |
| shared signatures | 31 |
| hazards | 20 |

The census denominator is **328 human pore-forming genes** across **75 families** in **32 superfamilies**, with **28 further families catalogued specifically so they can be excluded** — auxiliary subunits, domain-sharing non-channels and channels that are not ion channels.

## 1. Internal consistency

`src.catalogue.validate()` is clean: every family points at a superfamily that exists, no human gene is claimed by two census families, every exemplar is resolvable, every `confusable_with` and every hazard family reference matches a catalogue key.

## 2. Pfam accessions

132 distinct accessions declared across the catalogue; **132 verified (100 %)**, 0 missing from InterPro, 0 carrying a different short name.

**The ten largest search spaces** — how many UniProt proteins carry each signature. This is the census's raw denominator, and the reason S2 is a task of its own:

| accession | name | UniProt proteins | families |
|---|---|---|---|
| PF00005 | ABC_tran | 1,659,432 | 3 |
| PF12796 | Ank_2 | 542,452 | 4 |
| PF13855 | LRR_8 | 364,110 | 1 |
| PF00664 | ABC_membrane | 291,689 | 3 |
| PF00571 | CBS | 250,156 | 2 |
| PF13499 | EF-hand_7 | 239,985 | 2 |
| PF00248 | Aldo_ket_red | 239,151 | 1 |
| PF07686 | V-set | 217,115 | 1 |
| PF00023 | Ank | 210,260 | 1 |
| PF00520 | Ion_trans | 206,358 | 20 |

## 3. Exemplars

195 exemplars declared; 188 resolved cleanly, 0 unresolvable, 0 with a declared accession that names a different gene, 0 where UniProt's primary symbol differs from the catalogue's.

| status | n |
|---|---|
| ok | 188 |
| NO_GENE_NAME_IN_UNIPROT | 7 |

Every exemplar that did not resolve cleanly:

| family | label | gene | species | declared | resolved | resolved gene | status |
|---|---|---|---|---|---|---|---|
| kcsa_prok | Bc_NaK | nak | Bacillus cereus | Q81HW2 | Q81HW2 | — | NO_GENE_NAME_IN_UNIPROT |
| nav | Ab_NavAb | navAb | Aliarcobacter butzleri | A8EVM5 | A8EVM5 | — | NO_GENE_NAME_IN_UNIPROT |
| plgic_prok | Ec_ELIC | elic | Dickeya chrysanthemi | P0C7B7 | P0C7B7 | — | NO_GENE_NAME_IN_UNIPROT |
| iglur_prok | Ss_GluR0 | glr0 | Synechocystis sp. PCC 6803 | P73797 | P73797 | — | NO_GENE_NAME_IN_UNIPROT |
| deg_invertebrate | Ha_FaNaC | fanac | Cornu aspersum | Q25011 | Q25011 | — | NO_GENE_NAME_IN_UNIPROT |
| bestrophin | Kp_BEST | best | Klebsiella pneumoniae | W9BH30 | W9BH30 | — | NO_GENE_NAME_IN_UNIPROT |
| nonchannel_achbp | Ls_AChBP | achbp | Lymnaea stagnalis | P58154 | P58154 | — | NO_GENE_NAME_IN_UNIPROT |

`reference_panel.fasta` holds 195 sequences — the panel the classifier's reference tier scores against.

## 4. Declared architecture vs observed

195 exemplars had their Pfam architecture re-derived. **157 match the catalogue's declaration; 38 do not.**

A mismatch is not automatically an error in the catalogue — a signature declared for a family need not be present on every member, and that is exactly what hazard **H7** records for the TRP families. It is, however, always something a person should look at.

| family | exemplar | declared but not observed | observed but not declared |
|---|---|---|---|
| kv_shaker | Hs_KCNA1 | PF03521 | — |
| kv_shaker | Hs_KCND2 | PF03521 | PF11601,PF11879 |
| kv_shaker | Dm_Shaker | PF03521 | — |
| kca_slo | Hs_KCNT1 | PF00520,PF21014 | PF07885 |
| kir | Hs_KCNJ11 | PF08466 | — |
| nav | Ab_NavAb | PF06512,PF11933 | — |
| cav | Hs_CACNA1G | PF08763,PF16885,PF16905 | — |
| cng | Hs_CNGB1 | PF16526 | — |
| trpc | Dm_trp | PF23317 | PF00023,PF00520 |
| trpm | Hs_TRPM8 | PF00520,PF16519 | — |
| trpm | Hs_TRPM2 | PF16519,PF23317 | PF25969 |
| trpm | Hs_TRPM7 | PF00520 | PF02816 |
| trpn | Ce_trp-4 | PF00520 | PF00023,PF13637,PF23317 |
| zac | Hs_ZACN | PF02932 | — |
| plgic_prok | Gv_GLIC | PF02932 | — |
| plgic_prok | Ec_ELIC | PF02932 | — |
| clc_channel | Hs_CLCN1 | PF00571 | — |
| clc_channel | Hs_CLCN2 | PF00571 | — |
| bestrophin | Kp_BEST | PF01062 | PF25539 |
| mscs | At_MSL10 | PF05552,PF21082 | PF25886 |
| connexin | Hs_GJB2 | PF03508 | — |
| hv1 | Ci_Hv1 | PF16799 | — |
| clic | Hs_CLIC4 | PF13410 | — |
| viroporin | IAV_M2 | PF00558,PF02723,PF11289 | — |
| viroporin | HIV1_Vpu | PF00599,PF02723,PF11289 | — |
| viroporin | SARS2_E | PF00558,PF00599,PF11289 | — |
| nonchannel_ncs | Hs_RCVRN | PF00036 | PF13833 |
| assoc_k_beta | Hs_KCNAB1 | PF00930,PF03185 | — |
| assoc_k_beta | Hs_KCNMB1 | PF00248,PF00930 | — |
| assoc_k_beta | Hs_DPP6 | PF00248,PF03185 | PF00326 |
| assoc_cav_aux | Hs_CACNB1 | PF00822,PF08473,PF15108 | — |
| assoc_cav_aux | Hs_CACNA2D1 | PF00625,PF00822,PF12052,PF15108 | PF00092,PF08399,PF30670 |
| assoc_cav_aux | Hs_CACNG2 | PF00625,PF08473,PF12052,PF15108 | — |
| nonchannel_abc_transporter | Hs_ABCC1 | PF03412 | PF24357 |
| nonchannel_abc_transporter | Ec_MsbA | PF03412 | — |
| assoc_mcu_reg | Hs_MICU1 | PF10161 | PF13202 |
| assoc_mcu_reg | Hs_SMDT1 | PF13833 | — |
| nonchannel_kctd | Hs_KCTD1 | PF23110,PF31093,PF31099,PF31104 | — |

## 5. Shared signatures — the hazard registry's evidence

31 signatures are carried by more than one family. **19 of them cross the channel / non-channel boundary** — a domain that is evidence for a channel family and is also carried by something the catalogue does not count as a channel.

| accession | families | statuses |
|---|---|---|
| PF00520 | kv_shaker,kv_modifier,kv_kcnq,kv_eag,kca_slo,nav,cav,nalcn,catsper,tpc,cng,hcn,t | channel,non_channel_homolog |
| PF02931 | nachr,gabaa,glyr,ht3,zac,plgic_invertebrate,plgic_prok,nonchannel_achbp | channel,non_channel_homolog |
| PF02932 | nachr,gabaa,glyr,ht3,zac,plgic_invertebrate,plgic_prok | channel |
| PF00060 | ampa,kainate,nmda,delta_glur,iglur_nonvertebrate | channel,channel_contested |
| PF01094 | ampa,kainate,nmda,delta_glur,nonchannel_class_c_gpcr | channel,channel_contested,non_channel_homolog |
| PF07885 | kca_sk,k2p,kcsa_prok,iglur_prok | channel |
| PF10613 | ampa,kainate,nmda,delta_glur | channel,channel_contested |
| PF12796 | trpc,trpv,trpa,trpn | channel |
| PF00005 | cftr,assoc_sur,nonchannel_abc_transporter | channel,channel_associated,non_channel_homolog |
| PF00027 | kv_eag,cng,hcn | channel |
| PF00654 | clc_channel,clc_transporter,clc_prokaryotic | channel,transporter |
| PF00664 | cftr,assoc_sur,nonchannel_abc_transporter | channel,channel_associated,non_channel_homolog |
| PF00858 | enac,asic,deg_invertebrate | channel |
| PF02214 | kv_shaker,kv_modifier,nonchannel_kctd | channel,non_channel_homolog |
| PF02815 | itpr,ryr,nonchannel_pomt | channel,non_channel_homolog |

## 6. Taxonomy

52 taxon ids checked against the UniProt taxonomy service; 52 confirmed, 0 mismatched.

## 7. What this does not verify

- **Family membership.** That `KCNA1` belongs to `kv_shaker` is a literature assignment. S1 tests whether the classifier reproduces it; nothing here tests whether it is right.
- **Selectivity and gating.** Both are literature attributes of a family (hazard **H14**) and are not recoverable from sequence, so no database check can confirm them.
- **Completeness.** A verified catalogue is not a complete one. Whether 320 is the right number of human pore-forming genes is the question S2 answers, and it is answered against a declared search space, never against "all ion channels".

