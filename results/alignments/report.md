# S6 — alignments and pore modules

Rendered by `scripts/s6_report.py` from the tables in this directory (D13). Bulk FASTA, alignments, models and domtbls: `<data root>/alignments/s6/`.

![S6](figures/alignments_modules.png)

## 1. The alignment sets (D39)

**8,526 census v4 rows are called to a census family; 6,658 enter an alignment** (323 of them genome loci) across 68 families. Excluded, each counted in `members.tsv`:

| reason | rows |
|---|---|
| profile_medium | 1,697 |
| no_profile_call | 91 |
| genome_not_intact | 71 |
| duplicate_of | 9 |

The rule was fixed before any alignment ran: a sequence enters its family's alignment only on a **high-confidence** S3a profile call (D32 margin ≥ 0.30, ≥ half the profile covered), and a genome locus only on an intact reading frame. Medium calls are where S3b found the ankyrin/LRR upper bound and where sister families sit inside the margin; a tier-1 tree is the wrong place to carry either.

## 2. Family alignments — MAFFT L-INS-i + trimAl

**63 families aligned**, 5 with fewer than 4 sequences (no alignment), 0 not run. One method for every family — `mafft --localpair --maxiterate 1000`, then `trimal -automated1`; each run checked for row count, raggedness and unchanged residues, and recorded by SHA-256 (input, alignment, trimmed) so a rerun resumes by family.

L-INS-i wall time summed over families: 4.4 h. trimAl keeps a median 22% of columns (range 3%–92%).

| family | sequences | L-INS-i cols | gap frac | trimmed cols | trim gap | L-INS-i s |
|---|---|---|---|---|---|---|
| nachr | 726 | 8159 | 0.9386 | 208 | 0.1304 | 6514.8 |
| connexin | 402 | 1997 | 0.832 | 154 | 0.0524 | 307.7 |
| k2p | 380 | 5149 | 0.9184 | 145 | 0.0698 | 791.5 |
| cng | 313 | 8643 | 0.8899 | 365 | 0.0705 | 2869.8 |
| kv_shaker | 312 | 3152 | 0.8291 | 278 | 0.0272 | 515.0 |
| kir | 298 | 3241 | 0.8672 | 257 | 0.0583 | 306.3 |
| gabaa | 279 | 3152 | 0.8484 | 303 | 0.0636 | 196.7 |
| cav | 201 | 9686 | 0.7812 | 833 | 0.0305 | 1700.9 |
| asic | 187 | 3187 | 0.8358 | 244 | 0.0668 | 75.3 |
| trpc | 174 | 3631 | 0.7561 | 569 | 0.0701 | 230.0 |
| kv_eag | 158 | 3525 | 0.7109 | 612 | 0.0695 | 221.3 |
| vdac | 146 | 1289 | 0.7609 | 151 | 0.0535 | 16.4 |
| trpm | 145 | 5131 | 0.7043 | 580 | 0.0274 | 309.7 |
| osca_tmem63 | 143 | 6383 | 0.8603 | 472 | 0.0511 | 230.5 |
| nav | 138 | 4526 | 0.5872 | 1381 | 0.0735 | 289.9 |
| nmda | 131 | 4080 | 0.7045 | 530 | 0.0344 | 158.0 |
| trpv | 115 | 2907 | 0.7364 | 388 | 0.0331 | 37.6 |
| p2x | 114 | 2292 | 0.806 | 263 | 0.0634 | 22.5 |
| kca_slo | 106 | 4017 | 0.7111 | 411 | 0.0337 | 79.4 |
| bestrophin | 103 | 1952 | 0.7332 | 279 | 0.0359 | 16.4 |
| lrrc8 | 102 | 1326 | 0.39 | 650 | 0.0368 | 17.3 |
| tpc | 100 | 4689 | 0.8143 | 453 | 0.0315 | 109.3 |
| kv_kcnq | 97 | 2399 | 0.6807 | 291 | 0.0128 | 39.8 |
| kv_modifier | 90 | 942 | 0.4735 | 400 | 0.0622 | 5.4 |
| itpr | 89 | 7365 | 0.6272 | 1501 | 0.0225 | 426.1 |
| calhm | 82 | 690 | 0.532 | 207 | 0.0303 | 4.4 |
| innexin | 82 | 1022 | 0.5906 | 238 | 0.0116 | 9.6 |
| otop | 81 | 2798 | 0.7766 | 338 | 0.0262 | 43.9 |
| kca_sk | 80 | 1973 | 0.6895 | 280 | 0.0146 | 13.9 |
| kainate | 78 | 1676 | 0.4472 | 773 | 0.0626 | 12.7 |
| clc_channel | 72 | 2323 | 0.6255 | 510 | 0.0298 | 15.2 |
| mscs | 67 | 3158 | 0.7871 | 156 | 0.0349 | 49.5 |
| glyr | 65 | 854 | 0.447 | 532 | 0.1321 | 3.5 |
| trpp | 65 | 1809 | 0.5791 | 549 | 0.0875 | 9.3 |
| hcn | 64 | 4806 | 0.7882 | 431 | 0.0297 | 36.7 |
| ampa | 62 | 1963 | 0.5369 | 853 | 0.0085 | 7.1 |
| trpml | 62 | 1231 | 0.5293 | 349 | 0.0229 | 7.1 |
| ryr | 56 | 8742 | 0.4377 | 3664 | 0.0362 | 176.5 |
| orai | 54 | 1851 | 0.8444 | 154 | 0.018 | 3.1 |
| ht3 | 52 | 614 | 0.2843 | 353 | 0.0622 | 2.3 |
| trpa | 49 | 3566 | 0.6561 | 550 | 0.0249 | 21.8 |
| tric | 48 | 441 | 0.3425 | 198 | 0.0462 | 1.4 |
| pannexin | 46 | 1528 | 0.6733 | 297 | 0.0231 | 3.3 |
| piezo | 45 | 5009 | 0.493 | 1220 | 0.0279 | 47.3 |
| tmem175 | 45 | 871 | 0.3971 | 348 | 0.0155 | 2.9 |
| mcu | 43 | 994 | 0.6513 | 306 | 0.0147 | 1.5 |
| catsper | 41 | 1362 | 0.6429 | 232 | 0.0258 | 3.9 |
| enac | 41 | 1053 | 0.3955 | 415 | 0.0309 | 3.0 |
| iglur_nonvertebrate | 38 | 1290 | 0.2923 | 551 | 0.0331 | 4.5 |
| ano_channel | 31 | 1800 | 0.4541 | 789 | 0.048 | 2.8 |
| hv1 | 26 | 1375 | 0.6945 | 144 | 0.0059 | 2.3 |
| nalcn | 24 | 2501 | 0.3193 | 1211 | 0.0343 | 4.1 |
| trpn | 19 | 2182 | 0.259 | 1600 | 0.0327 | 3.1 |
| cftr | 14 | 1605 | 0.0776 | 1471 | 0.0088 | 1.6 |
| deg_invertebrate | 12 | 1290 | 0.458 | 641 | 0.178 | 1.6 |
| p2x_nonmetazoan | 10 | 587 | 0.3448 | 264 | 0.0394 | 0.8 |
| ano_scramblase | 9 | 1639 | 0.4848 | 645 | 0.0722 | 1.1 |
| clic | 8 | 727 | 0.5275 | 191 | 0.0183 | 0.8 |
| plgic_invertebrate | 7 | 754 | 0.3579 | 402 | 0.0057 | 0.9 |
| mscl | 6 | 174 | 0.1964 | 97 | 0.0189 | 0.8 |
| plgic_prok | 5 | 1281 | 0.5716 | 278 | 0.0165 | 1.0 |
| kcsa_prok | 4 | 401 | 0.2431 | 307 | 0.158 | 0.8 |
| tmc | 4 | 1113 | 0.2727 | 594 | 0.0446 | 1.0 |

Too few sequences for an alignment: delta_glur, iglur_prok, tweety, viroporin, zac.

## 3. Pore modules for the tier-2 units (D40)

Five alignable superfamilies need a module rather than the full length (`Superfamily.module_rule`): P-loop, iGluR and Ca²⁺-release (`pore_loop`: the helix before each re-entrant pore loop through the helix after it), innexin clan and Hv (`tm_span`).

**One extraction method for every sequence in a unit: `profile_projection`.** Each member is `hmmalign`ed to its own family's S3a profile and cut at that family's module span of match states. What differs between families is how the span was fixed, recorded in `module_spans.tsv`:

* **annotated** (29 families) — the UniProt topology of the family's annotated references (58 reference proteins), mapped to profile states, median over references;
* **unit HMM vote** (6: cng, iglur_nonvertebrate, iglur_prok, itpr, trpc, trpn) — no reference in the family has an annotated pore loop, so a module HMM built from the annotated references' modules is searched over the members and those with the expected number of full hits vote. A vote boundary more than 12 residues off every predicted helix moves outward to the bracketing helix: this fired once, for ITPR, whose RyR-seeded vote began in the luminal loop and would have dropped TM5.

| unit | seed modules | seed families | model states |
|---|---|---|---|
| ca_release | 3 | 1 | 102 |
| hv | 2 | 1 | 120 |
| iglur | 6 | 4 | 88 |
| innexin_like | 6 | 3 | 265 |
| ploop | 58 | 20 | 99 |

### The rejected design, measured

A single module HMM per unit was the first extraction instrument and was measured leave-one-family-out before use (`module_loo.tsv`): held-out families whose chains then carry exactly the expected number of modules —

| family | held out | exact / members |
|---|---|---|
| itpr | False | 80/89 |
| ampa | True | 62/62 |
| kainate | True | 78/78 |
| nmda | True | 127/131 |
| delta_glur | True | 3/3 |
| iglur_nonvertebrate | False | 27/38 |
| iglur_prok | False | 1/1 |
| pannexin | True | 0/46 |
| innexin | True | 5/82 |
| lrrc8 | True | 0/102 |
| kv_shaker | True | 305/312 |
| kv_modifier | True | 90/90 |
| kv_kcnq | True | 95/97 |
| kv_eag | True | 145/158 |
| kca_slo | True | 98/106 |
| kca_sk | True | 80/80 |
| kir | True | 22/298 |
| k2p | True | 207/380 |
| kcsa_prok | True | 3/4 |
| nav | True | 125/138 |
| cav | True | 190/201 |
| nalcn | True | 22/24 |
| catsper | True | 41/41 |
| tpc | True | 82/100 |
| cng | False | 163/313 |
| hcn | True | 0/64 |
| trpc | False | 102/174 |
| trpv | True | 74/115 |
| trpm | True | 0/145 |
| trpa | True | 18/49 |
| trpml | True | 2/62 |
| trpp | True | 60/65 |
| trpn | False | 10/19 |

It reaches the Kv, Nav/Cav and vertebrate iGluR families and fails on Kir, HCN, TRPM, TRPML and the whole innexin clan — so it is not the extractor. It survives only as the vote, where it is measured again (below).

### Validation

**Extracted against annotated, on members never used as references: 441 modules, median Jaccard 0.97, 436/441 ≥ 0.8.** **The vote, run with the held-out model on the annotated P-loop and iGluR families, lands within 12 profile states of the annotated span** (median 4) — the region the six vote families occupy. On the innexin clan it does not work (80 states), and is not used there.

**Extraction**: 5,228 full modules (≥ half the span's states occupied), 60 partial, 30 absent. ca_release 142, hv 26, iglur 312, innexin_like 225, ploop 4,523. In the K⁺-filter families the TxGYG filter lies inside the extracted module in **1,772/1,834** (96.6%). Of the 62 misses, 45 are chains with no canonical filter anywhere (degenerate or non-K⁺ filters, e.g. NaK's TVGDG) and 17 carry one elsewhere in the chain — 16 of them K2P modules, where one of the two pore domains carries a non-canonical filter and the other the canonical one.

## 4. What S7 and S8 take from here

* **S7** roots each family alignment by adding outgroup sequences with `mafft --add --keeplength`, so the family columns do not move.
* **S8** chooses representatives per clade × kingdom (D8) from the unit module FASTA and aligns the modules; every module carries `profile_projection` and its family's span basis in `modules.tsv`.
* The profiles used are S3a's frozen library — the instrument that made the calls. The S2b R3 seed re-draw is a census revision, not an S6 step (roadmap emergent row).
