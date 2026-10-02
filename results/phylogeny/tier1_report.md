# S7 — tier-1 phylogenies

Rendered by `scripts/s7_report.py` from `tier1_trim_compare.tsv`, `tier1_inputs.tsv` and `tier1_trees.tsv` (D13). Treefiles: `tier1/`. Bulk (alignments with outgroups, IQ-TREE output): `<data root>/trees/s7/`.

![S7](figures/tier1_trees.png)

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

**63 trees.** 36 have an outgroup. In 6 it is a single sequence (ampa, clc_channel, iglur_nonvertebrate, kainate, nmda, p2x) — the root edge is then a terminal edge with no support value, and the root is as good as the outgroup's placement. Of the 30 with several outgroup sequences, the outgroup forms one clade — the family is monophyletic with respect to it and the root is defined — in **26**, with UFBoot ≥ 95 on the root edge in **22**. Median share of a family's internal edges at UFBoot ≥ 95: **0.54**. Substitution matrices chosen: Q.pfam 23, JTT 22, LG 18. Wall time summed over families: 147.8 h.

**Outgroup not one clade** (root undefined — reported, not repaired): cng (outgroup kcsa_prok), k2p (outgroup kcsa_prok), kv_kcnq (outgroup kcsa_prok), kv_shaker (outgroup kcsa_prok).

| family | seqs | cols | model | rooting | root UFBoot | internal edges | median UFBoot | ≥ 95 | < 70 |
|---|---|---|---|---|---|---|---|---|---|
| nachr | 726 | 440 | Q.pfam+R10 | on plgic_prok | 100.0 | 724 | 96.0 | 0.5276 | 0.2003 |
| connexin | 402 | 278 | JTT+F+R9 | no declared outgroup | — | 399 | 87.0 | 0.391 | 0.3008 |
| k2p | 380 | 302 | Q.pfam+R9 | on kcsa_prok | not a clade | 365 | 88.0 | 0.3973 | 0.2986 |
| cng | 313 | 768 | Q.pfam+F+R10 | on kcsa_prok | not a clade | 304 | 97.0 | 0.5362 | 0.2368 |
| kv_shaker | 312 | 458 | Q.pfam+I+R9 | on kcsa_prok | not a clade | 303 | 92.0 | 0.4587 | 0.2277 |
| kir | 298 | 377 | Q.pfam+R8 | on kcsa_prok | 100.0 | 296 | 88.0 | 0.4122 | 0.2331 |
| gabaa | 279 | 441 | JTT+R7 | on plgic_prok | 100.0 | 277 | 93.0 | 0.4729 | 0.1949 |
| cav | 201 | 1786 | Q.pfam+I+R8 | on kcsa_prok | 100.0 | 199 | 100.0 | 0.6884 | 0.1005 |
| asic | 187 | 460 | Q.pfam+R8 | on deg_invertebrate | 89.0 | 185 | 98.0 | 0.5838 | 0.1838 |
| trpc | 174 | 772 | Q.pfam+R8 | on kcsa_prok | 100.0 | 172 | 97.0 | 0.5349 | 0.2093 |
| kv_eag | 158 | 909 | JTT+F+I+R7 | on kcsa_prok | 69.0 | 156 | 98.0 | 0.6026 | 0.2179 |
| vdac | 146 | 282 | Q.pfam+F+I+R4 | no declared outgroup | — | 143 | 92.0 | 0.4336 | 0.2238 |
| trpm | 145 | 1355 | JTT+F+R7 | on kcsa_prok | 100.0 | 143 | 100.0 | 0.7413 | 0.0839 |
| osca_tmem63 | 143 | 718 | LG+F+R6 | no declared outgroup | — | 140 | 100.0 | 0.6429 | 0.1571 |
| nav | 138 | 1947 | JTT+F+R8 | on kcsa_prok | 100.0 | 136 | 100.0 | 0.6691 | 0.1618 |
| nmda | 131 | 1009 | JTT+F+R6 | on iglur_prok | — | 129 | 100.0 | 0.7054 | 0.1085 |
| trpv | 115 | 683 | Q.pfam+R8 | on kcsa_prok | 100.0 | 113 | 100.0 | 0.6637 | 0.0973 |
| p2x | 114 | 398 | Q.pfam+R6 | on p2x_nonmetazoan | — | 112 | 93.0 | 0.4732 | 0.1964 |
| kca_slo | 106 | 1033 | LG+F+I+R7 | on kcsa_prok | 31.0 | 104 | 100.0 | 0.6731 | 0.1731 |
| bestrophin | 103 | 447 | Q.pfam+R7 | no declared outgroup | — | 100 | 96.0 | 0.52 | 0.23 |
| lrrc8 | 102 | 806 | JTT+F+R6 | no declared outgroup | — | 99 | 99.0 | 0.5859 | 0.1919 |
| tpc | 100 | 732 | Q.pfam+R7 | on kcsa_prok | 99.0 | 98 | 100.0 | 0.6429 | 0.1735 |
| kv_kcnq | 97 | 667 | JTT+I+R5 | on kcsa_prok | not a clade | 86 | 96.0 | 0.5349 | 0.1977 |
| kv_modifier | 90 | 499 | JTT+F+I+R5 | on kcsa_prok | 92.0 | 88 | 99.0 | 0.6023 | 0.2045 |
| itpr | 89 | 2643 | Q.pfam+F+I+R8 | is superfamily outgroup | — | 86 | 100.0 | 0.7907 | 0.0581 |
| calhm | 82 | 311 | Q.pfam+I+R5 | no declared outgroup | — | 79 | 95.0 | 0.5063 | 0.2278 |
| innexin | 82 | 401 | LG+I+R6 | no declared outgroup | — | 79 | 94.0 | 0.481 | 0.2025 |
| otop | 81 | 522 | JTT+F+I+R5 | no declared outgroup | — | 78 | 100.0 | 0.6795 | 0.1282 |
| kca_sk | 80 | 544 | JTT+I+R5 | on kcsa_prok | 100.0 | 78 | 92.0 | 0.4103 | 0.2179 |
| kainate | 78 | 902 | JTT+I+R4 | on iglur_prok | — | 76 | 96.0 | 0.5132 | 0.2895 |
| clc_channel | 72 | 835 | Q.pfam+R6 | on clc_prokaryotic | — | 70 | 98.0 | 0.5286 | 0.1714 |
| mscs | 67 | 430 | Q.pfam+R6 | no declared outgroup | — | 64 | 97.0 | 0.5625 | 0.3125 |
| glyr | 65 | 454 | JTT+R4 | on plgic_prok | 100.0 | 63 | 94.0 | 0.4762 | 0.3175 |
| trpp | 65 | 749 | Q.pfam+I+R5 | on kcsa_prok | 100.0 | 63 | 96.0 | 0.5397 | 0.3175 |
| hcn | 64 | 937 | JTT+F+R5 | on kcsa_prok | 100.0 | 62 | 97.0 | 0.5645 | 0.2097 |
| ampa | 62 | 902 | JTT+G4 | on iglur_prok | — | 60 | 92.0 | 0.45 | 0.1833 |
| trpml | 62 | 558 | LG+R6 | on kcsa_prok | 100.0 | 60 | 96.0 | 0.5 | 0.2833 |
| ryr | 56 | 4964 | JTT+F+I+R6 | on itpr | 100.0 | 54 | 100.0 | 0.7037 | 0.0926 |
| orai | 54 | 256 | Q.pfam+R5 | no declared outgroup | — | 51 | 84.0 | 0.3725 | 0.3529 |
| ht3 | 52 | 457 | JTT+R5 | on plgic_prok | 100.0 | 50 | 97.0 | 0.54 | 0.2 |
| trpa | 49 | 1125 | LG+R6 | on kcsa_prok | 100.0 | 47 | 97.0 | 0.5532 | 0.1915 |
| tric | 48 | 292 | LG+F+I+G4 | no declared outgroup | — | 45 | 93.0 | 0.4889 | 0.3111 |
| pannexin | 46 | 437 | JTT+I+G4 | no declared outgroup | — | 43 | 90.0 | 0.4651 | 0.2326 |
| piezo | 45 | 2566 | Q.pfam+R6 | no declared outgroup | — | 42 | 100.0 | 0.6905 | 0.1429 |
| tmem175 | 45 | 524 | JTT+F+R5 | no declared outgroup | — | 42 | 99.0 | 0.5714 | 0.2381 |
| mcu | 43 | 326 | LG+R5 | no declared outgroup | — | 40 | 82.0 | 0.35 | 0.4 |
| catsper | 41 | 421 | LG+I+G4 | on kcsa_prok | 100.0 | 39 | 99.0 | 0.5897 | 0.0513 |
| enac | 41 | 654 | JTT+R4 | on deg_invertebrate | 99.0 | 39 | 100.0 | 0.8205 | 0.0 |
| iglur_nonvertebrate | 38 | 916 | Q.pfam+R5 | on iglur_prok | — | 36 | 96.0 | 0.5278 | 0.1389 |
| ano_channel | 31 | 981 | JTT+R5 | no declared outgroup | — | 28 | 100.0 | 0.6429 | 0.1429 |
| hv1 | 26 | 308 | LG+G4 | no declared outgroup | — | 23 | 95.0 | 0.5217 | 0.087 |
| nalcn | 24 | 1738 | LG+F+R4 | on kcsa_prok | 100.0 | 22 | 98.0 | 0.5909 | 0.2273 |
| trpn | 19 | 1634 | LG+I+R4 | on kcsa_prok | 100.0 | 17 | 82.0 | 0.4118 | 0.2941 |
| cftr | 14 | 1486 | JTT+F+R3 | no declared outgroup | — | 11 | 100.0 | 0.7273 | 0.0909 |
| deg_invertebrate | 12 | 762 | LG+I+G4 | is superfamily outgroup | — | 9 | 81.0 | 0.2222 | 0.3333 |
| p2x_nonmetazoan | 10 | 376 | LG+G4 | is superfamily outgroup | — | 7 | 95.0 | 0.5714 | 0.2857 |
| ano_scramblase | 9 | 796 | LG+G4 | no declared outgroup | — | 6 | 99.0 | 0.6667 | 0.1667 |
| clic | 8 | 242 | LG+G4 | no declared outgroup | — | 5 | 97.0 | 0.6 | 0.2 |
| plgic_invertebrate | 7 | 446 | Q.pfam+I+R2 | on plgic_prok | 100.0 | 5 | 95.0 | 0.6 | 0.2 |
| mscl | 6 | 143 | Q.pfam+R2 | no declared outgroup | — | 3 | 89.0 | 0.3333 | 0.0 |
| plgic_prok | 5 | 375 | LG+G4 | is superfamily outgroup | — | 2 | 49.0 | 0.0 | 1.0 |
| kcsa_prok | 4 | 350 | LG+G4 | is superfamily outgroup | — | 1 | 78.0 | 0.0 | 0.0 |
| tmc | 4 | 857 | LG+G4 | no declared outgroup | — | 1 | 53.0 | 0.0 | 1.0 |

## 4. Re-rooting the weakly rooted P-loop families (S7c, D47)

S7b's KcsA/MthK/NaK outgroup was not one clade in four P-loop families and weakly supported in two more. Each of the six now declares **two outgroups in the catalogue** (`ChannelFamily.root_with`), fixed before any re-rooted tree: a prokaryotic relative with the family's architecture (inline exemplars verified live by `s0_outgroups.py` → `s0_baseline/outgroups.tsv`) and two sister eukaryotic families. K2P has no prokaryotic relative of its own, so its prokaryotic choice remains the superfamily outgroup.

| family | outgroup 1 | outgroup 2 |
|---|---|---|
| cng | prok_cng: Ml_MloK1, St_SthK | sister_hcn_kv_eag: hcn, kv_eag |
| k2p | prok_2tm_k: kcsa_prok | sister_kir_kv_shaker: kir, kv_shaker |
| kca_slo | prok_6tm_k: Ap_KvAP, Mj_MVP, Ec_Kch | sister_kv_shaker_kv_kcnq: kv_shaker, kv_kcnq |
| kv_eag | prok_6tm_k: Ap_KvAP, Mj_MVP, Ec_Kch | sister_hcn_cng: hcn, cng |
| kv_kcnq | prok_6tm_k: Ap_KvAP, Mj_MVP, Ec_Kch | sister_kv_shaker_kv_eag: kv_shaker, kv_eag |
| kv_shaker | prok_6tm_k: Ap_KvAP, Mj_MVP, Ec_Kch | sister_kv_kcnq_kv_eag: kv_kcnq, kv_eag |

**Denser basal sampling (D47).** Each family's D39 set gains one member per early-diverging clade of S4b's dense panel (kingdom × phylum × class; outside Bilateria), the clade's highest-scoring high-confidence profile call, complete and inside the length band ±25 % (`reroot_basal.tsv`). The six sets therefore differ from S6's D39 sets — a stated exception. Re-aligned with L-INS-i; everything after is S7's method unchanged, and the ingroup alignment and trim mask are identical in a family's two trees.

| family | D39 set | basal added | clades already in D39 | trimmed cols | outgroup occupancy (1 / 2) |
|---|---|---|---|---|---|
| cng | 313 | 23 | 3 | 765 | 0.13,0.53 / 0.78,0.76,0.71,0.75 |
| k2p | 380 | 28 | 3 | 303 | 0.43,0.67,0.29 / 0.78,0.83,0.85,0.82,0.79,0.76 |
| kca_slo | 106 | 39 | 4 | 946 | 0.28,0.21,0.40 / 0.33,0.59,0.45,0.42,0.48,0.58 |
| kv_eag | 158 | 10 | 2 | 911 | 0.26,0.20,0.37 / 0.62,0.75,0.65,0.79 |
| kv_kcnq | 97 | 2 | 0 | 676 | 0.42,0.31,0.51 / 0.47,0.80,0.60,0.58,0.81,0.81 |
| kv_shaker | 312 | 5 | 1 | 445 | 0.63,0.46,0.59 / 0.84,0.78,0.87,0.82 |

**The root criterion, fixed before any re-rooted tree:** in both trees (a) the outgroup is one clade and (b) its edge has UFBoot ≥ 95, and (c) the ingroup's basal split is *identical* under the two outgroups. Anything else is `unresolved`; no root is picked by how it looks.

| family | tree | seqs | model | outgroup one clade | root UFBoot | basal split (sizes) | split UFBoot | basal picks in smaller clade | ingroup inside the outgroup's clade (basal picks) | ≥ 95 |
|---|---|---|---|---|---|---|---|---|---|---|
| cng | prok_cng | 336 | Q.pfam+F+I+R10 | False | — | — | — | — | 157 (18) | 0.5566 |
| cng | sister_hcn_kv_eag | 336 | Q.pfam+F+I+R10 | False | — | — | — | — | 152 (13) | 0.5424 |
| k2p | prok_2tm_k | 408 | Q.pfam+R9 | False | — | — | — | — | 24 (13) | 0.3756 |
| k2p | sister_kir_kv_shaker | 408 | Q.pfam+R9 | False | — | — | — | — | 1 (0) | 0.3571 |
| kca_slo | prok_6tm_k | 145 | LG+F+R8 | False | — | — | — | — | 80 (27) | 0.6176 |
| kca_slo | sister_kv_shaker_kv_kcnq | 145 | LG+F+R7 | True | 100.0 | 1,144 | —,69 | 1/1 | 0 (0) | 0.5944 |
| kv_eag | prok_6tm_k | 168 | JTT+F+R8 | True | 100.0 | 1,167 | —,98 | 0/1 | 0 (0) | 0.5904 |
| kv_eag | sister_hcn_cng | 168 | JTT+F+I+R8 | True | 100.0 | 2,166 | 53,44 | 1/2 | 0 (0) | 0.5964 |
| kv_kcnq | prok_6tm_k | 99 | JTT+I+R6 | False | — | — | — | — | 41 (2) | 0.5543 |
| kv_kcnq | sister_kv_shaker_kv_eag | 99 | JTT+R6 | True | 68.0 | 4,95 | 57,29 | 0/4 | 0 (0) | 0.5258 |
| kv_shaker | prok_6tm_k | 317 | Q.pfam+I+R8 | False | — | — | — | — | 4 (1) | 0.4345 |
| kv_shaker | sister_kv_kcnq_kv_eag | 317 | Q.pfam+I+R9 | False | — | — | — | — | 151 (5) | 0.4704 |

*Ingroup inside the outgroup's clade* is a near-miss measure, not a test: the fewest ingroup sequences sharing one side of an edge with the whole outgroup (0 when it is one clade; `s7_newick.outgroup_intruders`), with how many are S7c basal picks. It never makes an outgroup count as one clade.

**Near misses:** k2p / sister_kir_kv_shaker — 1 ingroup sequence(s) inside the outgroup clade (A0ABM4D6C1__Hydvul); kv_shaker / prok_6tm_k — 4 ingroup sequence(s) inside the outgroup clade (A0A4D9CT78__Nansal,A0ABM4CGI6__Hydvul,A0ABM4CGM1__Hydvul,A0ABM4CLI9__Hydvul).

**Where a root is defined (4 trees), it splits off a handful of sequences, never two substantial clades:** kca_slo / sister_kv_shaker_kv_kcnq 1,144 (smaller side L8HDA6__Acacas; UFBoot —,69); kv_eag / prok_6tm_k 1,167 (smaller side I7MAP2__Tetthe; UFBoot —,98); kv_eag / sister_hcn_cng 2,166 (smaller side A0AAW1PKB9___Mybis,L1K4G9__Guithe; UFBoot 53,44); kv_kcnq / sister_kv_shaker_kv_eag 4,95 (smaller side A1Z856__Dromel,GCA_964187895.1_OZ076453.1_21552645-22401613-__Corasp; UFBoot 57,29). A root that isolates one to four sequences is the pattern long-branch attraction to a distant outgroup produces; the criterion does not depend on that reading.

| family | (a) one clade | (b) UFBoot ≥ 95 | (c) same root | smaller-clade Jaccard | verdict | S7b (KcsA): one clade / root UFBoot | S7b split vs new |
|---|---|---|---|---|---|---|---|
| cng | False | False | False |  | **unresolved** — outgroup not one clade | False / — | — |
| k2p | False | False | False |  | **unresolved** — outgroup not one clade | False / — | — |
| kca_slo | False | False | False |  | **unresolved** — outgroup not one clade | True / 31.0 | — |
| kv_eag | True | True | False | 0.0 | **unresolved** — root moves with the outgroup | True / 69.0 | — |
| kv_kcnq | False | False | False |  | **unresolved** — outgroup not one clade | False / — | — |
| kv_shaker | False | False | False |  | **unresolved** — outgroup not one clade | False / — | — |

**Resolved: 0/6** (none). **Unresolved: cng, k2p, kca_slo, kv_eag, kv_kcnq, kv_shaker** — S11 must read these families' duplication order without a root, or from reconciliation with the species tree (its own job), and say so.

## 5. What these trees are and are not

* A **protein tree** of one family (tier 1, full length after D41 trimming). Cross-family structure is S8's pore-module trees; cross-superfamily relationships are the fold network, never a tree (D27).
* Sets are S6's D39 sets: high-confidence profile calls and intact genome loci only. A family's tree therefore holds what the panel's 52 species carry, not every sequence in UniProt.
* Where the outgroup is a two-helix prokaryotic pore and the family is a 24-helix chain, the root rests on the pore columns alone; the occupancy table above says how few. Root support is reported per tree, never assumed; § 4 re-roots the six P-loop families where it failed.
