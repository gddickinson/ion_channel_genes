# S3a — profile library and best-profile assignment → census v3a

*Rendered from the tables in this directory by `scripts/s3_report.py`. Census v2 = UniProtKB 2026_03.*

Sweep database: census v2 collapsed to **1,187,702 unique sequences** (1,245,200 records, 768,363,148 residues); every profile searched with `-Z` fixed to that size. Wall-clock per profile in `sweep_runs.tsv` (total 6.9 profile-hours).

## 1. The profile library

**90 profiles**, one per catalogue family (census and control families alike, so a decoy has a profile to win), from **834 seeds**: R1 curated human genes 455, R2 catalogue exemplars 36, R3 reviewed S2 family calls (one per species, round-robin over groups) 343. Seeds are disjoint across families by construction (the build aborts otherwise). MAFFT L-INS-i single-threaded, hmmbuild; SHA-256 of every seed set, alignment and profile in `profile_build.tsv`. Match states total 71,421.

**18 profiles rest on one or two sequences** — the families with no human genes and no S2 family calls. They are reported as thin, not padded:

| family | superfamily | n_seeds | match_states |
|---|---|---|---|
| ano_channel | tmem16_like | 2 | 1035 |
| assoc_stim | orai | 2 | 827 |
| assoc_sur | abc_channel | 2 | 1582 |
| clc_prokaryotic | clc | 1 | 473 |
| deg_invertebrate | deg_enac | 2 | 898 |
| delta_glur | iglur | 2 | 1021 |
| iglur_nonvertebrate | iglur | 2 | 1100 |
| iglur_prok | iglur | 1 | 397 |
| innexin | innexin_like | 2 | 546 |
| nonchannel_achbp | cysloop | 1 | 229 |
| nonchannel_pomt | ca_release | 2 | 812 |
| outofscope_bacterial_porin | porin | 1 | 362 |
| p2x_nonmetazoan | p2x | 1 | 378 |
| plgic_invertebrate | cysloop | 2 | 498 |
| plgic_prok | cysloop | 2 | 371 |
| trpa | ploop | 1 | 1119 |
| trpn | ploop | 2 | 1972 |
| zac | cysloop | 1 | 412 |

## 2. The instrument, measured before use

Gates fixed before any census result: best profile ≥ 30 bits, spanning ≥ 30% of its match states (D30's number, applied to profiles), and beating the runner-up by ≥ 10% of its own score (D7).

**A. S1 control panel, leave-one-out** (95 rebuilt profiles; 2 queries are their family's only seed and have no LOO profile): **69/71 positives correct**, 26/26 decoys not called to a channel family. For comparison, S1's three-tier classifier scored 50/72 and S2's reference-free calls 45 on the same proteins. (S1 counted 72 positives and 25 decoys; here the CLC transporter, catalogue status `transporter`, is scored as a decoy.) The stronger decoy test — called to its **own** control family, a positive call rather than a mere non-channel one — holds for **22/26**; the other 4 (4 with no profile hit at all) come from auxiliary-subunit families that pool *unrelated* proteins under one catalogue key — `assoc_k_beta` holds an aldo-keto reductase (Kvβ), the BK β subunits, the LRRC γ subunits and DPP6/10 — so once the query is left out, no remaining seed is its homologue. One profile cannot represent a family that is not a family; the census is unaffected (none of these proteins carries a pore signature).

**B. Held-out orthologues** — reviewed, non-human, non-seed census records whose gene symbol is a catalogue human gene (the symbol scores the call, it never makes it — H15): **534/534 correct**, 49 control-family orthologues, 0 called to a channel. **This set is almost entirely vertebrate** (see `benchmark_calls.tsv`), so it tests the profiles on the lineage they were seeded from; section 3 is the non-vertebrate test.

Every benchmark miss:

| set | accession | expected_family | loo | outcome | p_call | p_family | win_score | win_coverage | runner | rel_margin |
|---|---|---|---|---|---|---|---|---|---|---|
| A_s1_panel | O75762 | trpa | sole_seed | no_call | module |  | 183.7 | 0.2449 | trpv | 0.5972 |
| A_s1_panel | Q401N2 | zac | sole_seed | wrong_family_same_sf | family | ht3 | 122.1 | 0.7285 | nachr | 0.4341 |

## 3. Calibration against S2 (seeds excluded)

On the 346,082 records S2 called to a family (seed accessions removed), the profile makes the **same call on 309,846** and a **different family call on 10,121** — agreement 96.8 % where both instruments speak. The two read different evidence (Pfam architecture and filter motif vs a full-length profile margin), so agreement is a measurement, not a tautology.

| domain | v2_calls | agree | disagree | profile_superfamily_only | profile_module | profile_low_score | profile_no_hit | agreement |
|---|---|---|---|---|---|---|---|---|
| Archaea | 271 | 130 | 8 | 0 | 75 | 16 | 42 | 0.942 |
| Bacteria | 28380 | 22728 | 177 | 26 | 3289 | 770 | 1390 | 0.9923 |
| Eukaryota | 317366 | 286973 | 9936 | 960 | 17198 | 617 | 1670 | 0.9665 |
| Viruses | 65 | 15 | 0 | 0 | 19 | 9 | 22 | 1.0 |

Per family (largest 30; full table `calibration.tsv`):

| family | v2_calls | agree | disagree | profile_superfamily_only | profile_ambiguous | profile_module | agreement |
|---|---|---|---|---|---|---|---|
| nonchannel_kctd | 31747 | 26680 | 362 | 242 | 0 | 4299 | 0.9866 |
| k2p | 25978 | 25466 | 100 | 133 | 1 | 275 | 0.9961 |
| cav | 20057 | 19959 | 44 | 1 | 0 | 53 | 0.9978 |
| connexin | 19810 | 19646 | 0 | 0 | 0 | 107 | 1.0 |
| osca_tmem63 | 18504 | 18499 | 0 | 0 | 3 | 2 | 1.0 |
| kir | 17481 | 17298 | 0 | 0 | 0 | 178 | 1.0 |
| vdac | 17029 | 16373 | 0 | 0 | 0 | 332 | 1.0 |
| mscl | 16560 | 16524 | 2 | 0 | 0 | 15 | 0.9999 |
| kca_slo | 12988 | 12699 | 0 | 0 | 0 | 286 | 1.0 |
| nonchannel_achbp | 12566 | 59 | 6588 | 519 | 1 | 4880 | 0.0089 |
| nav | 12186 | 12165 | 0 | 0 | 0 | 21 | 1.0 |
| kv_kcnq | 8829 | 8196 | 0 | 0 | 0 | 592 | 1.0 |
| tmc | 8480 | 8221 | 0 | 0 | 0 | 256 | 1.0 |
| bestrophin | 8224 | 6018 | 1 | 0 | 0 | 2067 | 0.9998 |
| trpp | 7566 | 3601 | 2987 | 68 | 5 | 820 | 0.5466 |
| kv_eag | 7477 | 7468 | 0 | 0 | 0 | 2 | 1.0 |
| itpr | 7199 | 5768 | 0 | 0 | 0 | 1375 | 1.0 |
| tmem175 | 6981 | 1803 | 0 | 0 | 0 | 2994 | 1.0 |
| mcu | 6468 | 6188 | 3 | 0 | 0 | 218 | 0.9995 |
| clic | 6341 | 6278 | 0 | 0 | 0 | 56 | 1.0 |
| trpc | 6210 | 6210 | 0 | 0 | 0 | 0 | 1.0 |
| otop | 6100 | 5359 | 0 | 0 | 1 | 335 | 1.0 |
| mscs | 5982 | 5982 | 0 | 0 | 0 | 0 | 1.0 |
| ryr | 5610 | 5473 | 0 | 0 | 0 | 137 | 1.0 |
| calhm | 4934 | 4767 | 0 | 0 | 0 | 51 | 1.0 |
| cng | 4700 | 4670 | 0 | 0 | 0 | 30 | 1.0 |
| trpml | 4685 | 4495 | 0 | 0 | 0 | 186 | 1.0 |
| tweety | 4666 | 4441 | 1 | 0 | 0 | 216 | 0.9998 |
| kca_sk | 4544 | 4543 | 0 | 0 | 0 | 1 | 1.0 |
| lrrc8 | 4347 | 4099 | 0 | 0 | 0 | 239 | 1.0 |

## 4. Census v3a

**1,245,200 records**; **730,790 (58.7 %) now carry a family call** — both instruments 310,386, S2 only 26,115, profile only 394,289 — against 346,627 (27.8 %) in S2. Superfamily only 93,766; unassigned 407,787; **conflict 12,857**, kept and counted.

**What the profile-only calls have not been tested on.** The superfamily splits below (Cys-loop, iGluR, DEG/ENaC, P2X, CLC, TMEM16) are families S2 never called, so section 3 cannot check them, and section 2's orthologue test is vertebrate. Where those calls fall in invertebrates, plants and protists they rest on profiles seeded mostly from human genes and have no independent check yet.

| status | records |
|---|---|
| channel | 613270 |
| unassigned | 407787 |
| superfamily_only | 93766 |
| channel_contested | 42154 |
| non_channel_homolog | 38165 |
| transporter | 35214 |
| conflict | 12857 |
| channel_associated | 1987 |

### What happened to S2's superfamily-only calls

| superfamily | v2_superfamily_only | resolved_to_family | outcomes |
|---|---|---|---|
| ploop | 164,252 | 101,224 (61.6 %) | superfamily_only 60,311, kv_shaker 21,577, kv_modifier 10,432, cng 7,861 |
| cysloop | 70,846 | 66,512 (93.9 %) | nachr 30,613, gabaa 24,782, glyr 5,820, superfamily_only 4,330 |
| iglur | 56,369 | 46,401 (82.3 %) | kainate 12,315, iglur_nonvertebrate 10,429, nmda 10,096, superfamily_only 9,966 |
| clc | 47,655 | 41,140 (86.3 %) | clc_transporter 18,351, clc_prokaryotic 16,475, superfamily_only 6,515, clc_channel 6,314 |
| deg_enac | 25,925 | 16,542 (63.8 %) | asic 11,984, superfamily_only 9,377, enac 3,742, deg_invertebrate 816 |
| tmem16_like | 13,287 | 13,059 (98.3 %) | ano_scramblase 9,720, ano_channel 3,339, superfamily_only 228 |
| innexin_like | 11,653 | 9,741 (83.6 %) | innexin 7,290, pannexin 2,451, superfamily_only 1,910, conflict 2 |
| p2x | 8,295 | 7,902 (95.3 %) | p2x 7,736, superfamily_only 393, p2x_nonmetazoan 166 |

### What happened to S2's unassigned records

500,291 records S2 left unassigned: `no_hit` 320,251 (64.0 %); `family` 91,768 (18.3 %); `module` 70,499 (14.1 %); `low_score` 17,033 (3.4 %); `superfamily_only` 736 (0.1 %); `ambiguous` 4 (0.0 %). A `module` verdict is a profile match over less than 30% of the family profile and is deliberately not a call; `no_hit` means no profile reported the record at all. S2 traced most of its unassigned records to three co-domain signatures (cNMP, SBP_bac_3, PAS); which of these two verdicts those records received is not broken down here.

| p_call | v3_family | records |
|---|---|---|
| no_hit | - | 320251 |
| module | - | 70499 |
| family | mscs | 61863 |
| low_score | - | 17033 |
| family | ano_scramblase | 8732 |
| family | trpm | 2831 |
| family | osca_tmem63 | 2818 |
| family | trpv | 2499 |
| family | cng | 1839 |
| family | piezo | 1552 |
| family | gabaa | 1420 |
| family | nachr | 1163 |
| family | trpc | 1028 |
| superfamily_only | - | 736 |
| family | hcn | 709 |

### Calls that rest on S2 alone

`s2_only` is a family call S2 made that the profile did not confirm — usually because the profile abstained (`module`: the record matches under 30 % of the family profile). These calls are only as good as S2's rule, and section 6 shows one rule that is not good enough. Largest groups (`call_support.tsv`):

| v3_family | s2_tier | p_call | records |
|---|---|---|---|
| nonchannel_achbp | hazard | module | 4880 |
| nonchannel_kctd | hazard | module | 4299 |
| tmem175 | architecture | module | 2994 |
| bestrophin | architecture | module | 2067 |
| tmem175 | architecture | no_hit | 1407 |
| itpr | hazard | module | 1375 |
| tmem175 | architecture | low_score | 777 |
| kv_kcnq | architecture | module | 592 |
| nonchannel_achbp | hazard | superfamily_only | 519 |
| trpp | hazard | module | 439 |
| trpp | architecture | module | 381 |
| assoc_sur | hazard | module | 377 |

### Conflicts

Where S2 and the profile name different families. S2's deciding tier is shown because a filter-motif call and a profile call disagreeing is a different finding from an architecture call and a profile call disagreeing.

| s2_call | profile_family | s2_tier | records |
|---|---|---|---|
| nonchannel_achbp | nachr | hazard | 3726 |
| trpp | assoc_polycystin1 | hazard | 2750 |
| [ploop] | hv1 | - | 1973 |
| nonchannel_achbp | gabaa | hazard | 1566 |
| nonchannel_achbp | ht3 | hazard | 437 |
| nonchannel_achbp | plgic_prok | hazard | 433 |
| [ploop] | itpr | - | 375 |
| trpp | assoc_polycystin1 | architecture | 239 |
| [ploop] | iglur_prok | - | 196 |
| nonchannel_kctd | kv_shaker | hazard | 190 |
| nonchannel_achbp | plgic_invertebrate | hazard | 174 |
| nonchannel_achbp | glyr | hazard | 171 |
| nonchannel_kctd | kv_modifier | hazard | 168 |
| [ploop] | ryr | - | 83 |
| [ploop] | nonchannel_vsp | - | 81 |
| nonchannel_achbp | zac | hazard | 76 |
| cav | nalcn | motif | 44 |
| assoc_sur | nalcn | hazard | 27 |
| k2p | kca_sk | hazard | 23 |
| k2p | cng | hazard | 22 |
| k2p | kcsa_prok | hazard | 19 |
| k2p | kv_eag | hazard | 14 |
| k2p | kca_slo | hazard | 6 |
| k2p | kir | hazard | 5 |
| [ploop] | mscs | - | 4 |

## 5. Human census genes

**318/320 called to the right family in v3a**, against 173/320 in S2. Most human genes are R1 seeds of their own profile, so this row measures the merge, not the profiles' generalisation — section 2 is the generalisation test. Every human gene not called right:

| gene | expected_family | s2_family | p_call | p_family | v3_family | v3_basis |
|---|---|---|---|---|---|---|
| GLRA4 | glyr |  |  |  |  | not_enumerated |
| ZACN | zac | nonchannel_achbp | family | zac |  | conflict |

## 6. External check — the parent projects' censuses

The IP3R and PIEZO projects censused three of this catalogue's families with their own seeds, profiles and search spaces (vertebrate reference proteomes and genome sweeps). Their results are **compared, never imported** (CLAUDE.md: port the method, never a result). Only UniProt accessions can be matched; genome-derived models and Ensembl proteins are counted as unmatched.

**IP3R census v6, call against call:** of 15,601 parent ITPR/RYR calls on accessions census v2 also holds, census v3a makes the **same call on 12,204 (78.2 %)**; ITPR and RYR swapped on 502; another family on 39; the rest carry no family call here (verdicts beginning `v3_`).

| project | parent_call | verdict | records |
|---|---|---|---|
| ip3r_genes/census_v6 | ITPR | agree | 6426 |
| ip3r_genes/census_v6 | ITPR | not_a_uniprot_accession | 671 |
| ip3r_genes/census_v6 | ITPR | not_in_census_v2 | 222 |
| ip3r_genes/census_v6 | ITPR | v3_conflict | 375 |
| ip3r_genes/census_v6 | ITPR | v3_superfamily_only | 448 |
| ip3r_genes/census_v6 | ITPR | v3_unassigned | 848 |
| ip3r_genes/census_v6 | RYR | agree | 5778 |
| ip3r_genes/census_v6 | RYR | called_other_family | 39 |
| ip3r_genes/census_v6 | RYR | not_a_uniprot_accession | 570 |
| ip3r_genes/census_v6 | RYR | not_in_census_v2 | 394 |
| ip3r_genes/census_v6 | RYR | swapped_itpr_ryr | 502 |
| ip3r_genes/census_v6 | RYR | v3_conflict | 83 |
| ip3r_genes/census_v6 | RYR | v3_superfamily_only | 298 |
| ip3r_genes/census_v6 | RYR | v3_unassigned | 804 |
| ip3r_genes/census_v6 | conflict | parent_conflict__ours_unassigned | 2 |
| ip3r_genes/census_v6 | unassigned | not_in_census_v2 | 4 |
| ip3r_genes/census_v6 | unassigned | parent_unassigned__ours_itpr | 409 |
| ip3r_genes/census_v6 | unassigned | parent_unassigned__ours_unassigned | 192 |
| piezo_genes/census_v5 | piezo_like | called_other_family | 1 |
| piezo_genes/census_v5 | piezo_like | called_piezo | 5055 |
| piezo_genes/census_v5 | piezo_like | not_a_uniprot_accession | 168 |
| piezo_genes/census_v5 | piezo_like | not_in_census_v2 | 847 |
| piezo_genes/census_v5 | piezo_like | v3_unassigned | 2258 |

**The 502 swapped records are all one error, and it is S2's.** 502 of them were called ITPR by S2's hazard rule `H4-itpr` — `PF08709` *without* `PF02026`/`PF06459` — while the profile abstained (median length 76 aa). `PF08709` is a domain ITPR and RyR *share* (H4: they share every diagnostic domain), so the rule is an absence test, which this project's conventions forbid, and it calls N-terminal RyR fragments IP3 receptors. On evidence both families carry, the right answer is no family call — not ITPR, and not necessarily RYR either. The parent project's basis for each RYR call is in `parent_reason`; 384 of the 502 are its architecture-only calls. The profile does not repeat the error; the merge inherits it because an abstaining profile leaves S2's call standing.

PIEZO census v5 is a membership list that includes fragments and short-motif hits; 2,258 of its records are unassigned here. Whether those are fragments below this census's gates or real PIEZOs missed has not been checked.

Reverse direction — of this census's own itpr / ryr / piezo calls, how many the parent census holds at all. A `no` is expected where this census reaches beyond the parent's search space (non-vertebrates for IP3R before its S20 sweep, UniProt entries outside reference proteomes, a newer UniProt release):

| our_family | v3_basis | in_parent_census | records |
|---|---|---|---|
| itpr | both | yes | 5777 |
| itpr | profile_only | yes | 129 |
| itpr | s2_only | yes | 1431 |
| piezo | both | yes | 3512 |
| piezo | profile_only | no | 10 |
| piezo | profile_only | yes | 1542 |
| piezo | s2_only | yes | 1 |
| ryr | both | yes | 5481 |
| ryr | profile_only | yes | 160 |
| ryr | s2_only | yes | 137 |

5,656 individual disagreements are listed in `external_disagreements.tsv` with both sides' evidence. Parent files and their SHA-256:

| file | rows | sha256 |
|---|---|---|
| ip3r_genes/results/census_v6/census_v6.tsv | 18065 | ab0b4f80da7633201f74b546b1750f570745ae4bef77c7438f474ff667adfff1 |
| piezo_genes/results/census_v5/piezo_like_census.csv | 8329 | a59112fa96e8f977eb01c825caec9ce5ac826ebf91e490651a5927274be937db |

## Files

| file | bytes | sha256 | records |
|---|---|---|---|
| hmmer/s3/census_v3.tsv.gz | 21250688 | 71c5f87807707ee9399cc2917af31851d26dface6b9f074e06750c4336583d4c | 1245200 |
| hmmer/s3/profile_calls.tsv.gz | 32820224 | fb8ca242e9eddb9a7e1dcb6980a60b74aecef309184b1e875cfc994a81b6575f | 857837 |
