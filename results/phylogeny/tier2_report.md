# S8 — tier-2 pore-module trees and the tier-3 fold network

Rendered by `scripts/s8_report.py` from the committed tables (D13). Rules: **D48** (fixed before any tier-2 tree or structure comparison).

## 1. What gets a tier-2 tree

**7 units get a tier-2 tree** (alignable, ≥ 2 census families); **21 single-family superfamilies** are their tier-1 tree (8 of them have none: < 4 sequences or added after S6); **4 are refused by D27** (`build_tier2` raises `NotAlignable`) and go to the fold network — a result, not an error.

| superfamily | tier2 | families | kind | note |
|---|---|---|---|---|
| ca_release | tree | 2 | module |  |
| cysloop | tree | 7 | full_length |  |
| deg_enac | tree | 3 | full_length |  |
| iglur | tree | 6 | module |  |
| innexin_like | tree | 3 | module |  |
| p2x | tree | 2 | full_length |  |
| ploop | tree | 23 | module |  |
| bestrophin | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| calhm | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| clc | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| clic | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| connexin | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| gphr | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| hv | = tier 1 | 1 | module | one census family: its tier-1 tree is the tier-2 tree |
| mclc | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| mcu | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| mitok | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| orai | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| otopetrin | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| pac | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| piezo | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| porin | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| tmco1 | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| tmem109 | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| tmem175 | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| tmem87 | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| tric | = tier 1 | 1 | full_length | one census family: its tier-1 tree is the tier-2 tree |
| tweety | = tier 1 | 1 | full_length | one census family, and no tier-1 tree (< 4 sequences, or added after S6 — r4) |
| abc_channel | refused (D27) | 1 |  | A superfamily of one, whose relatives are all transporters. Its tier-2 'phylogeny' is a phylogeny of the ABCC subfamily, which is the correct answer and not a channel tree. |
| msc | refused (D27) | 2 |  | Two families under one heading for convenience, explicitly not one superfamily: `alignable=False` stops the driver from building a tree across them. |
| tmem16_like | refused (D27) | 4 |  | Homology asserted from structure, undetectable in sequence. `alignable=False` is enforced by the phylogeny driver: a tier-2 tree across this superfamily is refused, not attempted and caveated. |
| viroporin | refused (D27) | 1 |  | A functional grouping, marked non-alignable so it is never treed. |

## 2. Representatives (D8, D48)

Cell = family × panel group × module (repeat). Within a cell, members are ordered by centrality (mean within-family identity) and clustered greedily; a member joins the first representative at ≥ the threshold. One threshold (0.5) for every unit, except where it leaves more than 4 tips per informative site: the P-loop pore module is ~105 residues and trims to 85 columns, so 840 tips at 0.5 became 317 at 0.3 — fixed on alignment properties before any tree. Tip counts at every candidate threshold:

| unit | kind | sequences | cells | threshold | t30 | t40 | t50 | t60 | t70 | t80 | t90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ca_release | module | 142 | 12 | 0.5 | 13 | 20 | 27 | 36 | 46 | 70 | 89 |
| cysloop | full_length | 1136 | 17 | 0.5 | 61 | 173 | 281 | 392 | 519 | 711 | 924 |
| deg_enac | full_length | 240 | 7 | 0.5 | 42 | 71 | 100 | 131 | 169 | 198 | 216 |
| iglur | module | 312 | 13 | 0.5 | 18 | 25 | 30 | 40 | 49 | 61 | 98 |
| innexin_like | module | 225 | 5 | 0.5 | 13 | 33 | 59 | 85 | 114 | 140 | 184 |
| p2x | full_length | 124 | 10 | 0.5 | 10 | 15 | 32 | 53 | 86 | 98 | 109 |
| ploop | module | 4523 | 192 | 0.3 | 317 | 558 | 839 | 1186 | 1562 | 2060 | 2878 |

## 3. Inputs

MAFFT L-INS-i over each unit (rows, raggedness and residues checked), trimAl `-gt 0.5` (D41), IQ-TREE 2 with S7's `IQ_ARGS`. **Root**: the superfamily's `root_with` family; its tips are the outgroup only within the kingdom of its catalogue exemplars, and the exemplars are added as tips (module units: cut by S6's method) unless a tip already carries the sequence. The four animal `plgic_prok` members are therefore *ingroup* tips in the Cys-loop tree — where they fall is the answer to the S7a row. innexin clan: no declared outgroup, unrooted.

| unit | kind | n_tips | n_families | root_family | n_outgroup | cols | informative | gap_frac |
|---|---|---|---|---|---|---|---|---|
| ca_release | module | 30 | 2 | itpr | 25 | 104 | 96 | 0.0391 |
| cysloop | full_length | 282 | 7 | plgic_prok | 2 | 415 | 412 | 0.1355 |
| deg_enac | full_length | 101 | 3 | deg_invertebrate | 7 | 467 | 451 | 0.0837 |
| iglur | module | 30 | 6 | iglur_prok | 1 | 88 | 85 | 0.0716 |
| innexin_like | module | 59 | 3 |  | 0 | 265 | 258 | 0.09 |
| p2x | full_length | 33 | 2 | p2x_nonmetazoan | 8 | 387 | 358 | 0.086 |
| ploop | module | 319 | 23 | kcsa_prok | 4 | 90 | 90 | 0.0998 |

Outgroup tips excluded by the kingdom rule: cysloop: plgic_prok__A0A9J7KP96__Braflo,plgic_prok__V4A2P8__Lotgig,plgic_prok__V4AGW4__Lotgig,plgic_prok__A0ABM1AAJ3__Aplcal.

## 4. ITPR's module span, checked against structure (S6 row)

S6 fixed ITPR's span by a RyR-seeded vote plus the helix snap. On the three human ITPRs the extracted module starts at UniProt's TM5 and ends 2–3 residues past TM6 — the convention of every annotated family — and contains the GVGD filter; on cryo-EM 6DQJ (human ITPR3) it starts inside the TM5 helix, spans the pore helix and filter, and ends inside TM6.

| accession | module_start | module_end | uniprot_tm5 | uniprot_tm6 | start_minus_tm5 | end_minus_tm6 | filter_GVGD | structure | tm5_helix | pore_helix | tm6_helix | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Q14571 | 2395 | 2544 | 2395-2415 | 2522-2542 | 0 | 2 | 2499 |  |  |  |  | consistent |
| Q14573 | 2369 | 2520 | 2369-2389 | 2497-2517 | 0 | 3 | 2475 | 6DQJ | 2366-2390 | 2458-2468 | 2490-2536 | consistent |
| Q14643 | 2449 | 2601 | 2449-2469 | 2578-2598 | 0 | 3 | 2556 |  |  |  |  | consistent |

## 5. The tier-3 fold network (D27, D48)

One AlphaFold DB model per census family (first S0 exemplar with a model of exactly that accession), cut to its comparison unit (module 1 for module superfamilies, whole model otherwise; pLDDT ≥ 70 only), plus a Kv voltage-sensor node. **75 nodes, 2775 pairs**, TM-align (average-length TM-score) and Foldseek (exhaustive). Unmeasured: ryr (no AlphaFold DB model for any exemplar: P21817,Q92736,Q15413).

**Verdict rule, fixed first**: an edge is `supported` if the median family-pair TM-score is ≥ 0.5 *and* each side's best other superfamily is the other side; `not_distinguished` otherwise. An edge is structural similarity, never phylogeny — `network.json` carries no branch lengths or support values.

| a | b | median_tm | n_pairs | rank_b_for_a | rank_a_for_b | best_other | best_other_tm | verdict |
|---|---|---|---|---|---|---|---|---|
| tmem16_like | tmem16_like | 0.4839 | 6 | 1 | 1 | tmem16_like-tmem175 | 0.242 | not_distinguished |
| iglur | ploop | 0.6364 | 138 | 1 | 1 | ca_release-iglur | 0.5573 | supported |
| innexin_like | connexin | 0.5648 | 3 | 1 | 1 | innexin_like-tweety | 0.3938 | supported |
| ca_release | ploop | 0.556 | 23 | 2 | 2 | iglur-ploop | 0.6364 | not_distinguished |
| hv | ploop:VSD | 0.5922 | 1 | 1 | 1 | hv-tmco1 | 0.479 | supported |

Reading: small helical pore modules score 0.4–0.6 against almost any helical unit of their size, which is why the rank test, not the bar alone, carries the verdict. TMEM16/OSCA/TMC misses the bar (0.484) while sitting far above its best outside partner — `not_distinguished` as written. ITPR's module (65 confident residues) is as close to the iGluR pore as to the P-loop pore: it is a P-loop-like pore, but this measurement cannot place it nearer the P-loop than the iGluR.

**Strongest pairs no literature edge asserts** (for S18, not claims):

| a | b | median_tm |
|---|---|---|
| ca_release | iglur | 0.5573 |
| deg_enac | pac | 0.5077 |
| hv | tmco1 | 0.479 |
| ca_release | hv | 0.4707 |
| tmco1 | viroporin | 0.4409 |
| ca_release | viroporin | 0.4247 |
| ploop:VSD | tmco1 | 0.4211 |
| hv | orai | 0.4204 |

Figure: `results/phylogeny/figures/fold_network.png`.

## 6. Tier-2 trees (S8b)

**7 trees, 0 failures** (IQ-TREE 2, S7's `IQ_ARGS`; 9.1 run-hours summed, Cys-loop the longest at 5.5 h). Read with `scripts/s8_tier2.py parse` (`s8_parse.py`) under the D48 rules, none revisited after a tree was seen. Treefiles: `results/phylogeny/tier2/`.

### 6.1 Roots

**3 of 7 trees are rooted** on their declared outgroup (ca_release, iglur, p2x); **3 have an outgroup the ingroup splits** (cysloop, deg_enac, ploop) and are read unrooted, never re-rooted on a substitute; the innexin clan declares none. A single-tip outgroup (iGluR on GluR0) is a terminal edge with no support value.

| unit | n_tips | informative | model | root_family | n_outgroup | outgroup_monophyletic | root_ufboot | root_clade_sizes | root_clade_ufboot | root_small_clade |
|---|---|---|---|---|---|---|---|---|---|---|
| ca_release | 30 | 96 | LG+G4 | itpr | 25 | True | 100 | 1,4 | ,61 | ryr:1 |
| cysloop | 282 | 412 | LG+F+I+R8 | plgic_prok | 2 | False |  |  |  |  |
| deg_enac | 101 | 451 | LG+R6 | deg_invertebrate | 7 | False |  |  |  |  |
| iglur | 30 | 85 | Q.pfam+R3 | iglur_prok | 1 | True |  | 7,22 | 96,82 | iglur_nonvertebrate:7 |
| innexin_like | 59 | 258 | Q.pfam+I+R5 |  | 0 |  |  |  |  |  |
| p2x | 33 | 358 | LG+G4 | p2x_nonmetazoan | 8 | True | 100 | 1,24 | ,59 | p2x:1 |
| ploop | 319 | 90 | LG+F+R6 | kcsa_prok | 4 | False |  |  |  |  |

Where a root exists, the P2X and Ca²⁺-release trees' basal split isolates a single tip (one metazoan P2X against 24; one RyR module against four), the pattern S7d recorded at tier 1; the iGluR tree's separates seven of the eight non-vertebrate iGluR modules from the rest (UFBoot 96 / 82).

### 6.2 Support

Share of ingroup internal edges at UFBoot ≥ 95 and < 70:

| unit | internal_edges | ufboot_median | frac_ge95 | frac_lt70 |
|---|---|---|---|---|
| ca_release | 3 | 60.0 | 0.0 | 1.0 |
| cysloop | 278 | 86.0 | 0.3741 | 0.3058 |
| deg_enac | 92 | 95.0 | 0.5109 | 0.2935 |
| iglur | 27 | 73.0 | 0.1852 | 0.4815 |
| innexin_like | 56 | 87.0 | 0.375 | 0.3214 |
| p2x | 23 | 74.0 | 0.1304 | 0.4348 |
| ploop | 309 | 64.0 | 0.1521 | 0.5728 |

The P-loop pore-module tree is the weakest by construction — 319 tips on 90 informative columns (the D48 cap of 4 tips per site, met at the 0.3 threshold): its deep branching order is not resolved, and nothing in § 6.3 should be read as a statement about the order of the P-loop families' divergence.

### 6.3 Families as clades

Each family's tips tested as one clade — per module for families with several modules (four-repeat chains, K2P, TPC), whose repeats are not expected to group. Rooted sense where the tree is rooted, unrooted (one side of any edge) otherwise. **34 of 54 groups are one clade.** `n_intruders` is the fewest other tips that would have to join a group to make it one — a descriptive near-miss measure, never a pass.

| unit | family | module | n_tips | sense | one_clade | ufboot | n_intruders | intruder_families |
|---|---|---|---|---|---|---|---|---|
| ca_release | ryr |  | 5 | rooted | True | 100 | 0 |  |
| cysloop | gabaa |  | 22 | unrooted | False |  | 4 | glyr:2,plgic_invertebrate:2 |
| cysloop | glyr |  | 2 | unrooted | True | 100 | 0 |  |
| cysloop | ht3 |  | 11 | unrooted | True | 100 | 0 |  |
| cysloop | nachr |  | 237 | unrooted | False |  | 13 | ht3:11,zac:2 |
| cysloop | plgic_invertebrate |  | 2 | unrooted | True | 100 | 0 |  |
| cysloop | plgic_prok |  | 6 | unrooted | True | 87 | 0 |  |
| cysloop | zac |  | 2 | unrooted | True | 100 | 0 |  |
| deg_enac | asic |  | 84 | unrooted | True | 100 | 0 |  |
| deg_enac | deg_invertebrate |  | 7 | unrooted | False |  | 10 | enac:10 |
| deg_enac | enac |  | 10 | unrooted | True | 100 | 0 |  |
| iglur | ampa |  | 2 | rooted | True | 84 | 0 |  |
| iglur | delta_glur |  | 1 | rooted | True |  | 0 |  |
| iglur | iglur_nonvertebrate |  | 8 | rooted | False |  | 21 | nmda:15,kainate:3,ampa:2,delta_glur:1 |
| iglur | kainate |  | 3 | rooted | True | 89 | 0 |  |
| iglur | nmda |  | 15 | rooted | False |  | 7 | kainate:3,ampa:2,delta_glur:1,iglur_nonvertebrate:1 |
| innexin_like | innexin |  | 39 | unrooted | True | 99 | 0 |  |
| innexin_like | lrrc8 |  | 12 | unrooted | True | 100 | 0 |  |
| innexin_like | pannexin |  | 8 | unrooted | True | 70 | 0 |  |
| p2x | p2x |  | 25 | rooted | True | 100 | 0 |  |
| ploop | catsper |  | 7 | unrooted | False |  | 12 | tpc:12 |
| ploop | cav | m1 | 12 | unrooted | False |  | 55 | nav:15,tpc:12,cav:11,nalcn:10,catsper:7 |
| ploop | cav | m2 | 14 | unrooted | False |  | 172 | cav:35,tpc:31,nav:20,nalcn:19,trpc:16,trpm:9,trpp:9,trpv:8,catsper:7,trpa:7,trpml:6,trpn:5 |
| ploop | cav | m3 | 12 | unrooted | False |  | 113 | cav:36,tpc:31,nav:20,nalcn:19,catsper:7 |
| ploop | cav | m4 | 11 | unrooted | False |  | 95 | tpc:31,cav:23,nalcn:19,nav:15,catsper:7 |
| ploop | cng |  | 26 | unrooted | False |  | 23 | kir:11,kv_eag:8,hcn:4 |
| ploop | hcn |  | 4 | unrooted | False |  | 1 | cng:1 |
| ploop | k2p | m1 | 22 | unrooted | False |  | 111 | k2p:27,cng:26,kca_slo:14,kir:11,kv_eag:8,kv_shaker:6,kv_kcnq:5,hcn:4,kca_sk:4,kcsa_prok:4,kv_modifier:1,trpc:1 |
| ploop | k2p | m2 | 27 | unrooted | False |  | 35 | kca_slo:14,kv_shaker:6,kv_kcnq:5,kca_sk:4,kcsa_prok:4,k2p:1,kv_modifier:1 |
| ploop | kca_sk |  | 4 | unrooted | True | 89 | 0 |  |
| ploop | kca_slo |  | 14 | unrooted | True | 97 | 0 |  |
| ploop | kcsa_prok |  | 4 | unrooted | False |  | 32 | kca_slo:14,kv_shaker:6,kv_kcnq:5,kca_sk:4,k2p:2,kv_modifier:1 |
| ploop | kir |  | 11 | unrooted | True | 88 | 0 |  |
| ploop | kv_eag |  | 8 | unrooted | False |  | 18 | cng:14,hcn:4 |
| ploop | kv_kcnq |  | 5 | unrooted | True | 63 | 0 |  |
| ploop | kv_modifier |  | 1 | unrooted | True |  | 0 |  |
| ploop | kv_shaker |  | 6 | unrooted | False |  | 1 | kv_modifier:1 |
| ploop | nalcn | m1 | 5 | unrooted | True | 100 | 0 |  |
| ploop | nalcn | m2 | 5 | unrooted | True | 100 | 0 |  |
| ploop | nalcn | m3 | 5 | unrooted | True | 100 | 0 |  |
| ploop | nalcn | m4 | 4 | unrooted | True | 100 | 0 |  |
| ploop | nav | m1 | 5 | unrooted | True | 100 | 0 |  |
| ploop | nav | m2 | 5 | unrooted | True | 95 | 0 |  |
| ploop | nav | m3 | 5 | unrooted | True | 98 | 0 |  |
| ploop | nav | m4 | 5 | unrooted | True | 99 | 0 |  |
| ploop | tpc | m1 | 12 | unrooted | True | 95 | 0 |  |
| ploop | tpc | m2 | 19 | unrooted | False |  | 6 | nalcn:5,cav:1 |
| ploop | trpa |  | 7 | unrooted | True | 98 | 0 |  |
| ploop | trpc |  | 17 | unrooted | False |  | 171 | cav:49,tpc:31,nav:20,nalcn:19,trpm:9,trpp:9,trpv:8,catsper:7,trpa:7,trpml:6,trpn:5,k2p:1 |
| ploop | trpm |  | 9 | unrooted | True | 100 | 0 |  |
| ploop | trpml |  | 6 | unrooted | True | 100 | 0 |  |
| ploop | trpn |  | 5 | unrooted | True | 83 | 0 |  |
| ploop | trpp |  | 9 | unrooted | False |  | 51 | trpc:16,trpm:9,trpv:8,trpa:7,trpml:6,trpn:5 |
| ploop | trpv |  | 8 | unrooted | True | 88 | 0 |  |

Reading, from the table: every four-repeat chain's repeats group by repeat within Nav and NALCN (each repeat one clade at UFBoot ≥ 95), while Cav's repeats are scattered across the four-repeat and TPC tips — at 90 columns that is a resolution limit, not a finding (S11 owns the repeat order). Near misses are mostly a related family nested inside: the silent Kv modifier inside Shaker, ENaC inside the invertebrate degenerins, HCN with one CNG tip, 5-HT3 and ZAC inside nAChR, and in the iGluR module tree the vertebrate non-NMDA families beside NMDA tips.

### 6.4 The four animal `plgic_prok` tips (S7a row)

D48's kingdom rule kept the four animal members the prokaryotic-pLGIC profile calls (*Branchiostoma*, *Aplysia*, *Lottia* ×2) out of the Cys-loop outgroup as ingroup tips. Each placed tip's nested clades, smallest first:

| tip | level | clade_size | ufboot | composition |
|---|---|---|---|---|
| plgic_prok__A0A9J7KP96__Braflo | 1 | 4 | 100 | plgic_prok:3 |
| plgic_prok__A0A9J7KP96__Braflo | 2 | 5 | 69 | plgic_prok:4 |
| plgic_prok__A0A9J7KP96__Braflo | 3 | 6 | 87 | plgic_prok:5 |
| plgic_prok__A0A9J7KP96__Braflo | 4 | 32 | 100 | gabaa:22,plgic_prok:5,glyr:2,plgic_invertebrate:2 |
| plgic_prok__A0ABM1AAJ3__Aplcal | 1 | 3 | 99 | plgic_prok:2 |
| plgic_prok__A0ABM1AAJ3__Aplcal | 2 | 4 | 100 | plgic_prok:3 |
| plgic_prok__A0ABM1AAJ3__Aplcal | 3 | 5 | 69 | plgic_prok:4 |
| plgic_prok__A0ABM1AAJ3__Aplcal | 4 | 6 | 87 | plgic_prok:5 |
| plgic_prok__A0ABM1AAJ3__Aplcal | 5 | 32 | 100 | gabaa:22,plgic_prok:5,glyr:2,plgic_invertebrate:2 |
| plgic_prok__V4A2P8__Lotgig | 1 | 2 | 58 | plgic_prok:1 |
| plgic_prok__V4A2P8__Lotgig | 2 | 3 | 99 | plgic_prok:2 |
| plgic_prok__V4A2P8__Lotgig | 3 | 4 | 100 | plgic_prok:3 |
| plgic_prok__V4A2P8__Lotgig | 4 | 5 | 69 | plgic_prok:4 |
| plgic_prok__V4A2P8__Lotgig | 5 | 6 | 87 | plgic_prok:5 |
| plgic_prok__V4A2P8__Lotgig | 6 | 32 | 100 | gabaa:22,plgic_prok:5,glyr:2,plgic_invertebrate:2 |
| plgic_prok__V4AGW4__Lotgig | 1 | 2 | 58 | plgic_prok:1 |
| plgic_prok__V4AGW4__Lotgig | 2 | 3 | 99 | plgic_prok:2 |
| plgic_prok__V4AGW4__Lotgig | 3 | 4 | 100 | plgic_prok:3 |
| plgic_prok__V4AGW4__Lotgig | 4 | 5 | 69 | plgic_prok:4 |
| plgic_prok__V4AGW4__Lotgig | 5 | 6 | 87 | plgic_prok:5 |
| plgic_prok__V4AGW4__Lotgig | 6 | 32 | 100 | gabaa:22,plgic_prok:5,glyr:2,plgic_invertebrate:2 |

**They are nested in no Cys-loop family.** The four form one clade, which joins the two prokaryotic tips (GLIC, ELIC) — the whole `plgic_prok` set (6 tips) is one side of an edge at UFBoot 87. Because they sit between GLIC and ELIC, the declared two-tip prokaryotic outgroup is split and the Cys-loop tree stays unrooted under D48 (not re-rooted on the six-tip clade, which would be choosing a root after reading the tree). Either a bacterial-type pLGIC lineage in these animals (horizontal transfer, or retained from an ancient lineage) or divergent sequences drawn to the long prokaryotic branches; the tree cannot tell the two apart — S21 (genomic context, contamination check) carries the question.

Figure: `results/phylogeny/figures/tier2_trees.png`.
