# S3a — profile library and best-profile assignment → census v3a

*Rendered from the tables in this directory by `scripts/s3_report.py`. Census v2 = UniProtKB 2026_03.*

Sweep database: census v2 collapsed to **1,187,702 unique sequences** (1,245,200 records, 768,363,148 residues); every profile searched with `-Z` fixed to that size. Wall-clock per profile in `sweep_runs.tsv` (total 6.9 profile-hours).

## 1. The profile library

**104 profiles**, one per catalogue family (census and control families alike, so a decoy has a profile to win), from **908 seeds**: R1 curated human genes 491, R2 catalogue exemplars 51, R3 reviewed S2 family calls (one per species, round-robin over groups) 366. Seeds are disjoint across families by construction (the build aborts otherwise). MAFFT L-INS-i single-threaded, hmmbuild; SHA-256 of every seed set, alignment and profile in `profile_build.tsv`. Match states total 77,544.

**23 profiles rest on one or two sequences** — the families with no human genes and no S2 family calls. They are reported as thin, not padded:

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
| mitok | mitok | 1 | 411 |
| nonchannel_achbp | cysloop | 1 | 229 |
| nonchannel_bri3bp | tmem109 | 1 | 251 |
| nonchannel_pomt | ca_release | 2 | 812 |
| nonchannel_tmem87b | tmem87 | 2 | 556 |
| outofscope_bacterial_porin | porin | 1 | 362 |
| p2x_nonmetazoan | p2x | 1 | 378 |
| plgic_invertebrate | cysloop | 2 | 498 |
| plgic_prok | cysloop | 2 | 371 |
| tmem109 | tmem109 | 2 | 243 |
| tmem87 | tmem87 | 2 | 555 |
| trpa | ploop | 1 | 1119 |
| trpn | ploop | 2 | 1972 |
| zac | cysloop | 1 | 412 |

## 2. The instrument, measured before use

Gates fixed before any census result: best profile ≥ 30 bits, spanning ≥ 30% of its match states (D30's number, applied to profiles), and beating the runner-up by ≥ 10% of its own score (D7).

**A. S1 control panel, leave-one-out** (95 rebuilt profiles; 2 queries are their family's only seed and have no LOO profile): **69/71 positives correct**, 26/26 decoys not called to a channel family. For comparison, S1's three-tier classifier scored 50/72 and S2's reference-free calls 45 on the same proteins. (S1 counted 72 positives and 25 decoys; here the CLC transporter, catalogue status `transporter`, is scored as a decoy.) The stronger decoy test — called to its **own** control family, a positive call rather than a mere non-channel one — holds for **22/26**; the other 4 (4 with no profile hit at all) come from auxiliary-subunit families that pool *unrelated* proteins under one catalogue key — `assoc_k_beta` holds an aldo-keto reductase (Kvβ), the BK β subunits, the LRRC γ subunits and DPP6/10 — so once the query is left out, no remaining seed is its homologue. One profile cannot represent a family that is not a family; the census is unaffected (none of these proteins carries a pore signature).

**B. Held-out orthologues** — reviewed, non-human, non-seed census records whose gene symbol is a catalogue human gene (the symbol scores the call, it never makes it — H15): **744/745 correct**, 59 control-family orthologues, 0 called to a channel. **This set is almost entirely vertebrate** (see `benchmark_calls.tsv`), so it tests the profiles on the lineage they were seeded from; section 3 is the non-vertebrate test.

Every benchmark miss:

| set | accession | expected_family | loo | outcome | p_call | p_family | win_score | win_coverage | runner | rel_margin |
|---|---|---|---|---|---|---|---|---|---|---|
| A_s1_panel | O75762 | trpa | sole_seed | no_call | module |  | 183.7 | 0.2449 | trpv | 0.5972 |
| A_s1_panel | Q401N2 | zac | sole_seed | wrong_family_same_sf | family | ht3 | 122.1 | 0.7285 | nachr | 0.4341 |
| B_orthologues | Q28EW0 | tmem87 |  | wrong_family_same_sf | family | nonchannel_tmem87b | 764.5 | 0.9299 | tmem87 | 0.2736 |

## 3. Calibration against S2 (seeds excluded)

On the 311,997 records S2 called to a family (seed accessions removed), the profile makes the **same call on 298,772** and a **different family call on 174** — agreement 99.9 % where both instruments speak. The two read different evidence (Pfam architecture and filter motif vs a full-length profile margin), so agreement is a measurement, not a tautology.

| domain | v2_calls | agree | disagree | profile_superfamily_only | profile_module | profile_low_score | profile_no_hit | agreement |
|---|---|---|---|---|---|---|---|---|
| Archaea | 263 | 130 | 0 | 0 | 75 | 16 | 42 | 1.0 |
| Bacteria | 27279 | 22220 | 0 | 1 | 2922 | 770 | 1366 | 1.0 |
| Eukaryota | 284449 | 276419 | 174 | 139 | 6195 | 318 | 1196 | 0.9994 |
| Viruses | 6 | 3 | 0 | 0 | 1 | 2 | 0 | 1.0 |

Per family (largest 30; full table `calibration.tsv`):

| family | v2_calls | agree | disagree | profile_superfamily_only | profile_ambiguous | profile_module | agreement |
|---|---|---|---|---|---|---|---|
| k2p | 25980 | 25466 | 101 | 133 | 1 | 276 | 0.996 |
| cav | 20058 | 19959 | 44 | 1 | 1 | 53 | 0.9978 |
| connexin | 19809 | 19645 | 0 | 0 | 0 | 107 | 1.0 |
| osca_tmem63 | 18504 | 18498 | 0 | 0 | 4 | 2 | 1.0 |
| kir | 17487 | 17300 | 4 | 0 | 0 | 178 | 0.9998 |
| vdac | 17030 | 16373 | 1 | 0 | 0 | 332 | 0.9999 |
| mscl | 16560 | 16524 | 2 | 0 | 0 | 15 | 0.9999 |
| nonchannel_kctd | 14335 | 14297 | 4 | 0 | 0 | 34 | 0.9997 |
| kca_slo | 12989 | 12699 | 1 | 0 | 0 | 286 | 0.9999 |
| nav | 12183 | 12162 | 0 | 0 | 0 | 21 | 1.0 |
| kv_kcnq | 8829 | 8196 | 0 | 0 | 0 | 592 | 1.0 |
| tmc | 8480 | 8221 | 0 | 0 | 0 | 256 | 1.0 |
| bestrophin | 8224 | 6018 | 1 | 0 | 0 | 2067 | 0.9998 |
| kv_eag | 7477 | 7468 | 0 | 0 | 0 | 2 | 1.0 |
| tmem175 | 6981 | 1803 | 0 | 0 | 0 | 2994 | 1.0 |
| mcu | 6468 | 6188 | 3 | 0 | 0 | 218 | 0.9995 |
| clic | 6341 | 6277 | 1 | 0 | 0 | 56 | 0.9998 |
| trpc | 6211 | 6210 | 1 | 0 | 0 | 0 | 0.9998 |
| otop | 6100 | 5359 | 0 | 0 | 1 | 335 | 1.0 |
| mscs | 5982 | 5982 | 0 | 0 | 0 | 0 | 1.0 |
| ryr | 5610 | 5473 | 0 | 0 | 0 | 137 | 1.0 |
| calhm | 4934 | 4767 | 0 | 0 | 0 | 51 | 1.0 |
| cng | 4700 | 4670 | 0 | 0 | 0 | 30 | 1.0 |
| trpml | 4686 | 4494 | 2 | 0 | 0 | 186 | 0.9996 |
| tweety | 4666 | 4441 | 1 | 0 | 0 | 216 | 0.9998 |
| kca_sk | 4544 | 4543 | 0 | 0 | 0 | 1 | 1.0 |
| lrrc8 | 4347 | 4099 | 0 | 0 | 0 | 239 | 1.0 |
| assoc_polycystin1 | 4134 | 3969 | 5 | 4 | 1 | 153 | 0.9987 |
| gphr | 3601 | 3596 | 0 | 0 | 0 | 3 | 1.0 |
| orai | 3575 | 3247 | 0 | 0 | 0 | 125 | 1.0 |

## 4. Census v3a

**1,271,983 records**; **749,734 (58.9 %) now carry a family call** — both instruments 299,320, S2 only 13,051, profile only 437,363 — against 312,545 (24.6 %) in S2 (r3 (S2c, 2026-09-28: H2/H4/H11/H12/H13 positive tests)). Superfamily only 109,749; unassigned 409,578; **conflict 2,922**, kept and counted.

**What the profile-only calls have not been tested on.** The superfamily splits below (Cys-loop, iGluR, DEG/ENaC, P2X, CLC, TMEM16) are families S2 never called, so section 3 cannot check them, and section 2's orthologue test is vertebrate. Where those calls fall in invertebrates, plants and protists they rest on profiles seeded mostly from human genes and have no independent check yet.

| status | records |
|---|---|
| channel | 618815 |
| unassigned | 409578 |
| superfamily_only | 109749 |
| channel_contested | 51634 |
| non_channel_homolog | 39891 |
| transporter | 35213 |
| channel_associated | 4181 |
| conflict | 2922 |

### What happened to S2's superfamily-only calls

| superfamily | v2_superfamily_only | resolved_to_family | outcomes |
|---|---|---|---|
| ploop | 183,801 | 115,273 (62.7 %) | superfamily_only 65,803, kv_shaker 21,767, nonchannel_kctd 12,394, kv_modifier 10,600 |
| cysloop | 83,413 | 73,154 (87.7 %) | nachr 34,338, gabaa 26,348, superfamily_only 10,250, glyr 5,991 |
| iglur | 56,369 | 46,400 (82.3 %) | kainate 12,315, iglur_nonvertebrate 10,428, nmda 10,096, superfamily_only 9,966 |
| clc | 47,655 | 41,139 (86.3 %) | clc_transporter 18,351, clc_prokaryotic 16,474, superfamily_only 6,516, clc_channel 6,314 |
| deg_enac | 25,927 | 16,541 (63.8 %) | asic 11,983, superfamily_only 9,377, enac 3,742, deg_invertebrate 816 |
| tmem16_like | 13,287 | 13,059 (98.3 %) | ano_scramblase 9,720, ano_channel 3,339, superfamily_only 228 |
| innexin_like | 11,653 | 9,741 (83.6 %) | innexin 7,290, pannexin 2,451, superfamily_only 1,910, conflict 2 |
| p2x | 8,295 | 7,902 (95.3 %) | p2x 7,736, superfamily_only 393, p2x_nonmetazoan 166 |
| ca_release | 7,208 | 5,777 (80.1 %) | itpr 5,777, superfamily_only 1,431 |

### What happened to S2's unassigned records

521,830 records S2 left unassigned: `no_hit` 320,879 (61.5 %); `family` 108,377 (20.8 %); `module` 71,556 (13.7 %); `low_score` 17,140 (3.3 %); `superfamily_only` 3,875 (0.7 %); `ambiguous` 3 (0.0 %). A `module` verdict is a profile match over less than 30% of the family profile and is deliberately not a call; `no_hit` means no profile reported the record at all. S2 traced most of its unassigned records to three co-domain signatures (cNMP, SBP_bac_3, PAS); which of these two verdicts those records received is not broken down here.

| p_call | v3_family | records |
|---|---|---|
| no_hit | - | 320879 |
| module | - | 71556 |
| family | mscs | 61863 |
| low_score | - | 17140 |
| family | ano_scramblase | 8732 |
| family | nonchannel_emc3 | 4294 |
| family | nonchannel_gost | 4242 |
| superfamily_only | - | 3875 |
| family | trpm | 2831 |
| family | osca_tmem63 | 2817 |
| family | trpv | 2499 |
| family | tmco1 | 2259 |
| family | nonchannel_tmem87b | 1949 |
| family | cng | 1839 |
| family | piezo | 1552 |

### Calls that rest on S2 alone

`s2_only` is a family call S2 made that the profile did not confirm — usually because the profile abstained (`module`: the record matches under 30 % of the family profile). These calls are only as good as S2's rule: the first S3a merge found three S2 rules that were not (H2, H4, H13 — rewritten in S2b, see the census v2 report). Largest groups (`call_support.tsv`):

| v3_family | s2_tier | p_call | records |
|---|---|---|---|
| tmem175 | architecture | module | 2994 |
| bestrophin | architecture | module | 2067 |
| tmem175 | architecture | no_hit | 1407 |
| tmem175 | architecture | low_score | 777 |
| kv_kcnq | architecture | module | 592 |
| otop | architecture | no_hit | 360 |
| otop | architecture | module | 335 |
| vdac | architecture | module | 332 |
| kca_slo | architecture | module | 286 |
| k2p | hazard | module | 276 |
| tmc | architecture | module | 256 |
| vdac | architecture | no_hit | 254 |

### Conflicts

Where S2 and the profile name different families. S2's deciding tier is shown because a filter-motif call and a profile call disagreeing is a different finding from an architecture call and a profile call disagreeing.

| s2_call | profile_family | s2_tier | records |
|---|---|---|---|
| [ploop] | hv1 | - | 1973 |
| [ploop] | itpr | - | 375 |
| [ploop] | iglur_prok | - | 196 |
| [ploop] | ryr | - | 83 |
| [ploop] | nonchannel_vsp | - | 81 |
| cav | nalcn | motif | 44 |
| k2p | kca_sk | hazard | 23 |
| k2p | cng | hazard | 22 |
| k2p | kcsa_prok | hazard | 19 |
| k2p | kv_eag | hazard | 14 |
| [ploop] | nonchannel_abc_transporter | - | 8 |
| k2p | kca_slo | hazard | 6 |
| [cysloop] | kv_eag | - | 5 |
| assoc_polycystin1 | trpp | hazard | 5 |
| k2p | kir | hazard | 5 |
| [ploop] | mscs | - | 4 |
| kir | assoc_sur | hazard | 4 |
| k2p | iglur_prok | hazard | 4 |
| nalcn | cav | motif | 3 |
| [deg_enac] | ano_scramblase | - | 3 |
| mcu | assoc_mcu_reg | architecture | 3 |
| [ploop] | ano_scramblase | - | 3 |
| k2p | assoc_catsper_aux | hazard | 2 |
| nonchannel_kctd | outofscope_aquaporin | hazard | 2 |
| [deg_enac] | nonchannel_abc_transporter | - | 2 |

## 5. Human census genes

**326/328 called to the right family in v3a**, against 173/328 in S2. Most human genes are R1 seeds of their own profile, so this row measures the merge, not the profiles' generalisation — section 2 is the generalisation test. Every human gene not called right:

| gene | expected_family | s2_family | p_call | p_family | v3_family | v3_basis |
|---|---|---|---|---|---|---|
| GLRA4 | glyr |  |  |  |  | not_enumerated |
| CCDC51 | mitok |  |  |  |  | not_enumerated |

## 6. External check — the parent projects' censuses

The IP3R and PIEZO projects censused three of this catalogue's families with their own seeds, profiles and search spaces (vertebrate reference proteomes and genome sweeps). Their results are **compared, never imported** (CLAUDE.md: port the method, never a result). Only UniProt accessions can be matched; genome-derived models and Ensembl proteins are counted as unmatched.

**IP3R census v6, call against call:** of 15,601 parent ITPR/RYR calls on accessions census v2 also holds, census v3a makes the **same call on 11,684 (74.9 %)**; ITPR and RYR swapped on 0; another family on 32; the rest carry no family call here (verdicts beginning `v3_`).

| project | parent_call | verdict | records |
|---|---|---|---|
| ip3r_genes/census_v6 | ITPR | agree | 5906 |
| ip3r_genes/census_v6 | ITPR | not_a_uniprot_accession | 671 |
| ip3r_genes/census_v6 | ITPR | not_in_census_v2 | 222 |
| ip3r_genes/census_v6 | ITPR | v3_conflict | 375 |
| ip3r_genes/census_v6 | ITPR | v3_superfamily_only | 968 |
| ip3r_genes/census_v6 | ITPR | v3_unassigned | 848 |
| ip3r_genes/census_v6 | RYR | agree | 5778 |
| ip3r_genes/census_v6 | RYR | called_other_family | 32 |
| ip3r_genes/census_v6 | RYR | not_a_uniprot_accession | 570 |
| ip3r_genes/census_v6 | RYR | not_in_census_v2 | 394 |
| ip3r_genes/census_v6 | RYR | v3_conflict | 83 |
| ip3r_genes/census_v6 | RYR | v3_superfamily_only | 810 |
| ip3r_genes/census_v6 | RYR | v3_unassigned | 801 |
| ip3r_genes/census_v6 | conflict | parent_conflict__ours_unassigned | 2 |
| ip3r_genes/census_v6 | unassigned | not_in_census_v2 | 4 |
| ip3r_genes/census_v6 | unassigned | parent_unassigned__ours_superfamily_only | 409 |
| ip3r_genes/census_v6 | unassigned | parent_unassigned__ours_unassigned | 192 |
| piezo_genes/census_v5 | piezo_like | called_other_family | 1 |
| piezo_genes/census_v5 | piezo_like | called_piezo | 5055 |
| piezo_genes/census_v5 | piezo_like | not_a_uniprot_accession | 168 |
| piezo_genes/census_v5 | piezo_like | not_in_census_v2 | 847 |
| piezo_genes/census_v5 | piezo_like | v3_unassigned | 2258 |

**No ITPR/RYR swap remains.** The first S3a merge found every swap to be one S2 error — the absence rule `H4-itpr` (`PF08709`, which both families carry, without the RyR domains ⇒ ITPR) calling N-terminal RyR fragments IP3 receptors. S2b removed it: ITPR has no positive architectural test, so those records are superfamily-only unless the profile calls them. That is why fewer parent ITPR calls are matched than before — the ones withdrawn had agreed by chance, on evidence that could not decide.

PIEZO census v5 is a membership list that includes fragments and short-motif hits; 2,258 of its records are unassigned here. Whether those are fragments below this census's gates or real PIEZOs missed has not been checked.

Reverse direction — of this census's own itpr / ryr / piezo calls, how many the parent census holds at all. A `no` is expected where this census reaches beyond the parent's search space (non-vertebrates for IP3R before its S20 sweep, UniProt entries outside reference proteomes, a newer UniProt release):

| our_family | v3_basis | in_parent_census | records |
|---|---|---|---|
| itpr | profile_only | yes | 5906 |
| piezo | both | yes | 3512 |
| piezo | profile_only | no | 10 |
| piezo | profile_only | yes | 1542 |
| piezo | s2_only | yes | 1 |
| ryr | both | yes | 5481 |
| ryr | profile_only | yes | 160 |
| ryr | s2_only | yes | 137 |

6,176 individual disagreements are listed in `external_disagreements.tsv` with both sides' evidence. Parent files and their SHA-256:

| file | rows | sha256 |
|---|---|---|
| ip3r_genes/results/census_v6/census_v6.tsv | 18065 | ab0b4f80da7633201f74b546b1750f570745ae4bef77c7438f474ff667adfff1 |
| piezo_genes/results/census_v5/piezo_like_census.csv | 8329 | a59112fa96e8f977eb01c825caec9ce5ac826ebf91e490651a5927274be937db |

## Files

| file | bytes | sha256 | records |
|---|---|---|---|
| hmmer/s3/census_v3.tsv.gz | 21744212 | 98d819275a11b98813eee5f2c2cb979462ae33858979876d79438b69f978dd87 | 1271983 |
| hmmer/s3/profile_calls.tsv.gz | 33847189 | b0295a0a021f40692704809365133c22784e5d0047011ba9d8b127d97c2171c4 | 881229 |
