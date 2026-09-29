# S7 — tier-1 phylogenies

Rendered by `scripts/s7_report.py` from `tier1_trim_compare.tsv`, `tier1_inputs.tsv` and `tier1_trees.tsv` (D13). Treefiles: `tier1/`. Bulk (alignments with outgroups, IQ-TREE output): `<data root>/trees/s7/`.

![S7](figures/tier1_trim.png)

## 1. Trimming, decided before any tree (D41)

S6 trimmed with trimAl `-automated1`, which on low-identity families switches to its strict mode and removes the variable-but-aligned columns that carry the within-family signal. Four candidates were measured on all 63 family alignments on alignment properties alone (`tier1_trim_compare.tsv`); **no tree was built or read to choose.**

| method | median columns | median informative sites | median gap fraction | median share of a member's residues kept | median of the per-family minimum |
|---|---|---|---|---|---|
| none | 1963 | 869 | 0.627 | 1.00 | 1.00 |
| automated1 | 388 | 324 | 0.034 | 0.61 | 0.32 |
| gappyout | 474 | 408 | 0.021 | 0.80 | 0.43 |
| gt0.5 | 558 | 504 | 0.073 | 0.91 | 0.53 |

**Chosen: `-gt 0.5`** — keep every column in which at least half the family has a residue, for every family. It keeps a median 504 informative sites against 324 for automated1, at least as many in **62/63** families (the exception, glyr: 281 vs 325), and keeps a median 91% of each member's own residues against 61%, at a median gap fraction of 7.3%. The three families S6 flagged: k2p 145 → 302 columns, nachr 208 → 440 columns, cng 365 → 768 columns.

## 2. Rooting — from the catalogue only (D41)

A family is rooted on the catalogue exemplars of its superfamily's declared outgroup family (`Superfamily.root_with`), added to its S6 alignment with `mafft --localpair --maxiterate 1000 --add --keeplength` (every ingroup row checked unchanged) and trimmed by the same D41 column mask. Nothing is midpoint-rooted.

| rooting | families |
|---|---|
| rooted (36) | ampa, asic, catsper, cav, clc_channel, cng, enac, gabaa, glyr, hcn, ht3, iglur_nonvertebrate, k2p, kainate, kca_sk, kca_slo, kir, kv_eag, kv_kcnq, kv_modifier, kv_shaker, nachr, nalcn, nav, nmda, p2x, plgic_invertebrate, ryr, tpc, trpa, trpc, trpm, trpml, trpn, trpp, trpv |
| unrooted:is_superfamily_outgroup (5) | deg_invertebrate, itpr, kcsa_prok, p2x_nonmetazoan, plgic_prok |
| unrooted:no_declared_outgroup (22) | ano_channel, ano_scramblase, bestrophin, calhm, cftr, clic, connexin, hv1, innexin, lrrc8, mcu, mscl, mscs, orai, osca_tmem63, otop, pannexin, piezo, tmc, tmem175, tric, vdac |

No tree at all (< 4 sequences in the D39 set): delta_glur, iglur_prok, tweety, viroporin, zac.

**Outgroup occupancy** is the share of the family's trimmed columns an outgroup residue fills. The P-loop families are rooted on KcsA, MthK and NaK — two-helix pores — so in a four-repeat or TRP alignment the outgroup is mostly gap:

| family | outgroup | occupancy |
|---|---|---|
| ampa | Ss_GluR0 | 0.43 |
| asic | Ce_mec-4,Ha_FaNaC | 0.95,0.97 |
| catsper | Sl_KcsA,Mt_MthK,Bc_NaK | 0.33,0.48,0.25 |
| cav | Sl_KcsA,Mt_MthK,Bc_NaK | 0.09,0.16,0.06 |
| clc_channel | Ec_ClC-ec1 | 0.55 |
| cng | Sl_KcsA,Mt_MthK,Bc_NaK | 0.19,0.42,0.15 |
| enac | Ce_mec-4,Ha_FaNaC | 0.86,0.76 |
| gabaa | Gv_GLIC,Ec_ELIC | 0.75,0.69 |
| glyr | Gv_GLIC,Ec_ELIC | 0.75,0.68 |
| hcn | Sl_KcsA,Mt_MthK,Bc_NaK | 0.16,0.30,0.11 |
| ht3 | Gv_GLIC,Ec_ELIC | 0.71,0.67 |
| iglur_nonvertebrate | Ss_GluR0 | 0.40 |
| k2p | Sl_KcsA,Mt_MthK,Bc_NaK | 0.43,0.68,0.32 |
| kainate | Ss_GluR0 | 0.42 |
| kca_sk | Sl_KcsA,Mt_MthK,Bc_NaK | 0.28,0.46,0.19 |
| kca_slo | Sl_KcsA,Mt_MthK,Bc_NaK | 0.15,0.30,0.10 |
| kir | Sl_KcsA,Mt_MthK,Bc_NaK | 0.37,0.60,0.26 |
| kv_eag | Sl_KcsA,Mt_MthK,Bc_NaK | 0.17,0.30,0.12 |
| kv_kcnq | Sl_KcsA,Mt_MthK,Bc_NaK | 0.23,0.42,0.17 |
| kv_modifier | Sl_KcsA,Mt_MthK,Bc_NaK | 0.30,0.30,0.21 |
| kv_shaker | Sl_KcsA,Mt_MthK,Bc_NaK | 0.32,0.34,0.22 |
| nachr | Gv_GLIC,Ec_ELIC | 0.71,0.67 |
| nalcn | Sl_KcsA,Mt_MthK,Bc_NaK | 0.09,0.18,0.06 |
| nav | Sl_KcsA,Mt_MthK,Bc_NaK | 0.07,0.16,0.05 |
| nmda | Ss_GluR0 | 0.38 |
| p2x | Dd_P2XA | 0.83 |
| plgic_invertebrate | Gv_GLIC,Ec_ELIC | 0.71,0.68 |
| ryr | Hs_ITPR1,Hs_ITPR2,Hs_ITPR3 | 0.47,0.46,0.46 |
| tpc | Sl_KcsA,Mt_MthK,Bc_NaK | 0.17,0.23,0.15 |
| trpa | Sl_KcsA,Mt_MthK,Bc_NaK | 0.12,0.19,0.10 |
| trpc | Sl_KcsA,Mt_MthK,Bc_NaK | 0.19,0.31,0.14 |
| trpm | Sl_KcsA,Mt_MthK,Bc_NaK | 0.10,0.22,0.07 |
| trpml | Sl_KcsA,Mt_MthK,Bc_NaK | 0.27,0.49,0.19 |
| trpn | Sl_KcsA,Mt_MthK,Bc_NaK | 0.09,0.19,0.07 |
| trpp | Sl_KcsA,Mt_MthK,Bc_NaK | 0.18,0.37,0.12 |
| trpv | Sl_KcsA,Mt_MthK,Bc_NaK | 0.21,0.36,0.15 |

## 3. Trees

IQ-TREE 2 (`-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -bnni -seed 1`) on every family; the model is ModelFinder's BIC choice among four general empirical matrices with their rate and frequency variants. Support is summarised over the family's own internal edges (the outgroup excluded).

*No trees yet.*

## 4. What these trees are and are not

* A **protein tree** of one family (tier 1, full length after D41 trimming). Cross-family structure is S8's pore-module trees; cross-superfamily relationships are the fold network, never a tree (D27).
* Sets are S6's D39 sets: high-confidence profile calls and intact genome loci only. A family's tree therefore holds what the panel's 52 species carry, not every sequence in UniProt.
* Where the outgroup is a two-helix prokaryotic pore and the family is a 24-helix chain, the root rests on the pore columns alone; the occupancy table above says how few. Root support is reported per tree, never assumed.
