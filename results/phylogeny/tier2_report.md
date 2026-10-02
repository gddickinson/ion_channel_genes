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

## 6. Tier-2 trees

Running detached (`scripts/s8_tier2.py run`, log `<data root>/trees/s8/run.log`); parsed in S8b.
