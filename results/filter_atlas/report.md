# S9 — Selectivity-filter atlas (Q5)

Rendered by `scripts/s9_report.py` from the tables beside it (D13). Rules: D49.

## 1. The instrument

4,604 P-loop pore modules (S6, D39 members) in 3,035 chains were added to S8's untrimmed P-loop tier-2 L-INS-i alignment (`mafft --localpair --add --keeplength --mapout`, one run per family chunk), so every module is read in one coordinate system. Anchors, located through MAFFT's residue→column map (counting residues along an aligned row is wrong wherever `--keeplength` deleted an insertion upstream — the first read made that error on Nav1.5 repeat III and was corrected before any tree was read):

| anchor | reference | position | residue | repeat | column |
|---|---|---|---|---|---|
| k_window | KcsA P0A334 | 75 | T |  | 433 |
| k_window | KcsA P0A334 | 76 | V |  | 434 |
| k_window | KcsA P0A334 | 77 | G |  | 449 |
| k_window | KcsA P0A334 | 78 | Y |  | 450 |
| k_window | KcsA P0A334 | 79 | G |  | 451 |
| repeat_locus | Nav1.5 Q14524 | 372 | D | 1 | 434 |
| repeat_locus | Nav1.5 Q14524 | 898 | E | 2 | 434 |
| repeat_locus | Nav1.5 Q14524 | 1419 | K | 3 | 434 |
| repeat_locus | Nav1.5 Q14524 | 1711 | A | 4 | 434 |

**All four Nav1.5 locus residues (D372, E898, K1419, A1711) fall in one column (434), the column KcsA's V76 — the x of T-x-G-Y-G — occupies.** The four-repeat locus and the K⁺ signature are read at the same alignment position (an alignment statement, not a structural one).

**Validation (D49 (4)).**

* Classifier's pairwise projection (MAFFT to Nav1.5) vs the shared-alignment read, four-repeat chains: **356 / 363 agree** (98.1 %; identical, or identical wherever both read a residue). Per family: nav {'yes': 128, 'yes_read': 10}; cav {'yes': 187, 'yes_read': 8, 'no': 6}; nalcn {'yes': 20, 'yes_read': 3, 'no': 1}; catsper {'not_comparable': 41}; tpc {'not_comparable': 100}. Every disagreement differs at one position, all in divergent chains (`projection_check.tsv`). CatSper (one repeat per chain) and TPC (two) are outside the projection's design and are not compared.
* K window vs `K_FILTER_RE` in the module: of 2,058 modules carrying a TxGYG-like motif, **1,962 read the same motif in the window**, 9 a different K motif (two in the module), 87 a non-K window (mostly CNG); 0 windows read a K motif the regex did not find. 81 modules have no sequence (S6 `absent`).
* Every tier-2 tip re-added reads what its own row reads: **317 / 317**.

## 2. The atlas: four-repeat strings

| family | chains | unread | strings |
|---|---|---|---|
| nav | 138 | 11 | DEKA 104, DEEA 12, DEKG 3, DKEA 2, DENS 2, DDQA 1 |
| cav | 201 | 6 | EEEE 112, EEDD 49, DDDD 11, EEDE 5, NDDN 3, QEEE 2 |
| nalcn | 24 | 2 | EEKE 18, EKEE 2, EDEE 1, EEEE 1 |
| catsper | 41 | 0 | D 39, E 2 |
| tpc | 100 | 5 | AN 69, SN 6, DA 6, SG 4, ES 3, DS 3 |

Strings with no `FILTER_CALLS` entry (DEEA, DKEA, EKEE, DDDD, …) are recorded as read, with their placement below; S9 adds no call (D49 (6)). K-window families: `filter_modules.tsv` / `filter_chains.tsv`.

## 3. Congruence (D49 (5))

Fitch changes on each tree against 1,000 tip-label permutations; RI = retention index (1 = every change once, 0 = no more structure than a star tree). Rooted (S7b outgroup one clade): catsper, cav, hcn, kca_sk, kir, kv_modifier, nalcn, nav, tpc, trpa, trpc, trpm, trpml, trpn, trpp, trpv; the rest, including the six S7d left undefined, are read unrooted.

| tree | family | character | rooted | n_tips | n_states | changes | min_changes | star_changes | ri | p_perm |
|---|---|---|---|---|---|---|---|---|---|---|
| tier1 | catsper | repeat_locus | True | 41 | 2 | 2 | 1 | 2 | 0.0 | 1.0 |
| tier1 | cav | repeat_locus | True | 195 | 19 | 20 | 18 | 83 | 0.9692 | 0.001 |
| tier1 | cng | k_window | False | 165 | 41 | 54 | 40 | 88 | 0.7083 | 0.001 |
| tier1 | hcn | k_window | True | 62 | 2 | 1 | 1 | 1 |  |  |
| tier1 | k2p | k_window | False | 344 | 40 | 82 | 39 | 209 | 0.7471 | 0.001 |
| tier1 | kca_sk | k_window | True | 80 | 7 | 15 | 6 | 33 | 0.6667 | 0.001 |
| tier1 | kca_slo | k_window | False | 104 | 6 | 15 | 5 | 40 | 0.7143 | 0.001 |
| tier1 | kcsa_prok | k_window | False | 4 | 3 | 2 | 2 | 2 |  |  |
| tier1 | kir | k_window | True | 294 | 14 | 18 | 13 | 60 | 0.8936 | 0.001 |
| tier1 | kv_eag | k_window | False | 157 | 9 | 19 | 8 | 45 | 0.7027 | 0.001 |
| tier1 | kv_kcnq | k_window | False | 96 | 3 | 5 | 2 | 9 | 0.5714 | 0.001 |
| tier1 | kv_modifier | k_window | True | 90 | 1 | 0 | 0 | 0 |  |  |
| tier1 | kv_shaker | k_window | False | 304 | 11 | 20 | 10 | 125 | 0.913 | 0.001 |
| tier1 | nalcn | repeat_locus | True | 22 | 4 | 4 | 3 | 4 | 0.0 | 1.0 |
| tier1 | nav | repeat_locus | True | 127 | 9 | 10 | 8 | 23 | 0.8667 | 0.001 |
| tier1 | tpc | repeat_locus | True | 95 | 8 | 12 | 7 | 26 | 0.7368 | 0.001 |
| tier1 | trpa | k_window | True | 48 | 28 | 30 | 27 | 42 | 0.8 | 0.001 |
| tier1 | trpc | k_window | True | 164 | 76 | 76 | 75 | 139 | 0.9844 | 0.001 |
| tier1 | trpm | k_window | True | 129 | 24 | 29 | 23 | 97 | 0.9189 | 0.001 |
| tier1 | trpml | k_window | True | 62 | 7 | 12 | 6 | 13 | 0.1429 | 0.1958 |
| tier1 | trpn | k_window | True | 15 | 10 | 10 | 9 | 11 | 0.5 | 0.2927 |
| tier1 | trpp | k_window | True | 65 | 34 | 38 | 33 | 57 | 0.7917 | 0.001 |
| tier1 | trpv | k_window | True | 94 | 48 | 50 | 47 | 69 | 0.8636 | 0.001 |
| tier2 | ploop | locus | False | 312 | 16 | 89 | 15 | 240 | 0.6711 | 0.001 |
| tier2 | ploop | window | False | 293 | 167 | 192 | 166 | 265 | 0.7374 | 0.001 |

**16 of 20 testable tier-1 trees: p ≤ 0.001.** Not significant: catsper, nalcn, trpml, trpn — each a family whose minority string is carried by one or two scattered chains.

### Per string (≥ 2 carriers), four-repeat families

| family | string | n_carriers | sense | one_clade | ufboot | n_intruders | intruder_strings | origins |
|---|---|---|---|---|---|---|---|---|
| catsper | D | 39 | rooted | False |  | 2 | E:2 | 6 |
| catsper | E | 2 | rooted | False |  | 39 | D:39 | 2 |
| cav | EEEE | 112 | rooted | False |  | 35 | DDDD:11,EEDE:5,NDDN:3,EEDD:2,QEEE:2,EFEE:1 | 31 |
| cav | EEDD | 49 | rooted | False |  | 146 | EEEE:112,DDDD:11,EEDE:5,NDDN:3,QEEE:2,IEEE:1 | 12 |
| cav | DDDD | 11 | rooted | False |  | 3 | NDDN:3 | 2 |
| cav | EEDE | 5 | rooted | False |  | 9 | EEEE:2,EEDD:2,DEDD:1,EVDE:1,EDEE:1,EESD:1 | 3 |
| cav | NDDN | 3 | rooted | True | 91 | 0 |  | 1 |
| cav | QEEE | 2 | rooted | True | 100 | 0 |  | 1 |
| nalcn | EEKE | 18 | rooted | False |  | 3 | EKEE:2,EEEE:1 | 2 |
| nalcn | EKEE | 2 | rooted | False |  | 19 | EEKE:18,EEEE:1 | 2 |
| nav | DEKA | 104 | rooted | False |  | 23 | DEEA:12,DEKG:3,DKEA:2,DENS:2,DEET:1,DENA:1 | 17 |
| nav | DEEA | 12 | rooted | False |  | 115 | DEKA:104,DEKG:3,DENS:2,DKEA:2,DETE:1,DEET:1 | 7 |
| nav | DEKG | 3 | rooted | False |  | 1 | DEKA:1 | 2 |
| nav | DENS | 2 | rooted | True | 100 | 0 |  | 1 |
| nav | DKEA | 2 | rooted | True | 100 | 0 |  | 1 |
| tpc | AN | 69 | rooted | False |  | 13 | SN:6,SG:4,VN:3 | 27 |
| tpc | DA | 6 | rooted | False |  | 89 | AN:69,SN:6,SG:4,ES:3,DS:3,VN:3 | 4 |
| tpc | SN | 6 | rooted | False |  | 76 | AN:69,SG:4,VN:3 | 5 |
| tpc | SG | 4 | rooted | True | 100 | 0 |  | 1 |
| tpc | DS | 3 | rooted | False |  | 92 | AN:69,DA:6,SN:6,SG:4,ES:3,VN:3 | 2 |
| tpc | ES | 3 | rooted | True | 100 | 0 |  | 1 |
| tpc | VN | 3 | rooted | True | 100 | 0 |  | 1 |

`origins` counts maximal carrier-only clades and is fragmented by any one nested change; read it with the core clade below.

### The Q5 test case: Cav3's EEDD

EEDD has 49 carriers and is not one clade under the strict test (146 intruders, because two carriers sit far away). **47 of the 49 form one clade at UFBoot 100**, a 49-chain clade that holds the three human T-type channels (O43497, O95180, Q9P0X4), invertebrate, cnidarian and placozoan chains, and nothing else but one EQDD and one chain with an unread repeat. The two EEDD chains outside it are both *Hydra* (A0ABM4BFT6, A0ABM4BFU2), read the same by both instruments. **EEDD sits where the tree says: one origin, in the lineage that carries it, with at most one independent EEDD (in *Hydra*).** Both of the positions where EEDD differs from EEEE (repeats III, IV) rank at the 95th–96th percentile of the Cav alignment's columns for tree congruence (§ 4). The 11 DDDD chains are all *Paramecium*; 8 form one clade (UFBoot 100), and with NDDN (also *Paramecium*) a 14-chain clade at UFBoot 100.

Other non-canonical strings, read descriptively (best-fitting clade = the clade maximising carriers minus others): Nav DEEA — 12 carriers, 6 in one clade (UFBoot 78), the other 6 scattered (*Hydra*, *Nematostella*, *Trichoplax*, amphioxus, and an elephant-shark and an opossum chain both instruments read DEEA, emergent); Nav DKEA — 2 carriers, one clade (UFBoot 100); NALCN EKEE — 2 carriers, not one clade; CatSper E — 2 carriers, scattered (RI 0).

### Tier 2 (P-loop pore modules, unrooted, 15 % of edges at UFBoot ≥ 95)

The locus column is congruent with the module tree beyond chance (table above), but the tree's deep order is unsupported (S8b), so no statement is made about where across families a filter type arose. Per-residue clade tests are in `string_clades.tsv`; no residue is one clade across families.

## 4. Post-hoc control: filter positions against every other column

Added after § 3 was read, and labelled so: the filter is inside the alignment each tree was inferred from, so § 3's permutation test shows only that the filter carries phylogenetic signal, as every column does. Here each single filter position's RI is ranked among every parsimony-informative column (≤ 50 % gaps) of the same tree's input alignment.

| tree | family | position | n_states | ri | n_columns | column_ri_median | percentile |
|---|---|---|---|---|---|---|---|
| tier1 | catsper | repeat 1 | 2 | 0.0 | 411 | 0.6364 | 1.2 |
| tier1 | cav | repeat 1 | 6 | 0.7647 | 1763 | 0.6867 | 66.2 |
| tier1 | cav | repeat 2 | 6 | 0.8667 | 1763 | 0.6867 | 84.9 |
| tier1 | cav | repeat 3 | 4 | 0.9315 | 1763 | 0.6867 | 95.1 |
| tier1 | cav | repeat 4 | 4 | 0.9394 | 1763 | 0.6867 | 96.1 |
| tier1 | cng | window 1 | 12 | 0.6949 | 767 | 0.4737 | 90.7 |
| tier1 | cng | window 2 | 12 | 0.7836 | 767 | 0.4737 | 96.7 |
| tier1 | cng | window 3 | 13 | 0.7234 | 767 | 0.4737 | 93.2 |
| tier1 | cng | window 4 | 11 | 0.75 | 767 | 0.4737 | 94.7 |
| tier1 | cng | window 5 | 13 | 0.7808 | 767 | 0.4737 | 96.6 |
| tier1 | k2p | window 1 | 11 | 0.7143 | 299 | 0.5543 | 87.6 |
| tier1 | k2p | window 2 | 13 | 0.6381 | 299 | 0.5543 | 71.9 |
| tier1 | k2p | window 4 | 9 | 0.7941 | 299 | 0.5543 | 95.7 |
| tier1 | kca_sk | window 1 | 4 | 0.7778 | 519 | 0.6667 | 72.1 |
| tier1 | kca_sk | window 2 | 4 | 0.5789 | 519 | 0.6667 | 33.1 |
| tier1 | kca_slo | window 2 | 5 | 0.0 | 1009 | 0.7115 | 2.6 |
| tier1 | kca_slo | window 4 | 2 | 0.7353 | 1009 | 0.7115 | 54.0 |
| tier1 | kir | window 1 | 5 | 0.9167 | 371 | 0.6588 | 96.5 |
| tier1 | kir | window 2 | 4 | 0.4 | 371 | 0.6588 | 5.7 |
| tier1 | kir | window 4 | 5 | 0.875 | 371 | 0.6588 | 93.0 |
| tier1 | kir | window 5 | 4 | 0.5 | 371 | 0.6588 | 15.9 |
| tier1 | kv_eag | window 1 | 3 | 0.8529 | 887 | 0.6941 | 81.2 |
| tier1 | kv_eag | window 2 | 4 | 0.7037 | 887 | 0.6941 | 52.1 |
| tier1 | kv_eag | window 4 | 2 | 1.0 | 887 | 0.6941 | 100.0 |
| tier1 | kv_kcnq | window 2 | 3 | 0.5714 | 602 | 0.6667 | 33.2 |
| tier1 | kv_kcnq | window 4 | 2 | 0.0 | 602 | 0.6667 | 3.7 |
| tier1 | kv_shaker | window 1 | 6 | 0.3333 | 444 | 0.6667 | 8.3 |
| tier1 | kv_shaker | window 2 | 4 | 0.9083 | 444 | 0.6667 | 96.2 |
| tier1 | nalcn | repeat 2 | 3 | 0.0 | 1103 | 0.6667 | 16.3 |
| tier1 | nalcn | repeat 3 | 2 | 0.8 | 1103 | 0.6667 | 67.1 |
| tier1 | nav | repeat 2 | 3 | 1.0 | 1806 | 0.5 | 100.0 |
| tier1 | nav | repeat 3 | 6 | 0.9444 | 1806 | 0.5 | 96.8 |
| tier1 | nav | repeat 4 | 5 | 0.6 | 1806 | 0.5 | 70.6 |
| tier1 | tpc | repeat 1 | 6 | 0.7619 | 727 | 0.6129 | 81.2 |
| tier1 | tpc | repeat 2 | 4 | 0.8462 | 727 | 0.6129 | 91.9 |
| tier1 | trpa | window 1 | 2 | 0.6667 | 1079 | 0.5294 | 77.6 |
| tier1 | trpa | window 2 | 5 | 0.7857 | 1079 | 0.5294 | 91.0 |
| tier1 | trpa | window 3 | 7 | 0.6957 | 1079 | 0.5294 | 81.2 |
| tier1 | trpa | window 4 | 8 | 0.6875 | 1079 | 0.5294 | 79.2 |
| tier1 | trpa | window 5 | 13 | 0.6364 | 1079 | 0.5294 | 71.3 |
| tier1 | trpc | window 1 | 8 | 0.9048 | 764 | 0.6667 | 98.2 |
| tier1 | trpc | window 2 | 12 | 0.575 | 764 | 0.6667 | 28.9 |
| tier1 | trpc | window 3 | 11 | 0.6562 | 764 | 0.6667 | 47.8 |
| tier1 | trpc | window 4 | 18 | 0.8421 | 764 | 0.6667 | 93.8 |
| tier1 | trpc | window 5 | 14 | 0.6721 | 764 | 0.6667 | 52.6 |
| tier1 | trpm | window 2 | 5 | 0.8 | 1324 | 0.6207 | 87.2 |
| tier1 | trpm | window 3 | 8 | 0.9306 | 1324 | 0.6207 | 98.0 |
| tier1 | trpm | window 4 | 8 | 0.8933 | 1324 | 0.6207 | 96.7 |
| tier1 | trpm | window 5 | 12 | 0.9038 | 1324 | 0.6207 | 97.1 |
| tier1 | trpml | window 3 | 3 | 0.0 | 519 | 0.5 | 4.6 |
| tier1 | trpml | window 4 | 4 | 0.5 | 519 | 0.5 | 52.6 |
| tier1 | trpml | window 5 | 4 | 0.2222 | 519 | 0.5 | 9.6 |
| tier1 | trpn | window 2 | 4 | 0.6667 | 1290 | 0.5 | 75.0 |
| tier1 | trpn | window 3 | 3 | 0.5 | 1290 | 0.5 | 58.5 |
| tier1 | trpn | window 4 | 6 | 1.0 | 1290 | 0.5 | 100.0 |
| tier1 | trpn | window 5 | 7 | 0.2 | 1290 | 0.5 | 23.9 |
| tier1 | trpp | window 2 | 8 | 1.0 | 720 | 0.5 | 100.0 |
| tier1 | trpp | window 3 | 9 | 0.6667 | 720 | 0.5 | 81.4 |
| tier1 | trpp | window 4 | 7 | 0.7143 | 720 | 0.5 | 86.7 |
| tier1 | trpp | window 5 | 14 | 0.6 | 720 | 0.5 | 69.7 |
| tier1 | trpv | window 1 | 5 | 0.875 | 663 | 0.6111 | 95.5 |
| tier1 | trpv | window 2 | 9 | 0.8507 | 663 | 0.6111 | 93.4 |
| tier1 | trpv | window 3 | 7 | 0.878 | 663 | 0.6111 | 95.5 |
| tier1 | trpv | window 4 | 12 | 0.8462 | 663 | 0.6111 | 93.4 |
| tier1 | trpv | window 5 | 12 | 0.4561 | 663 | 0.6111 | 16.6 |
| tier2 | ploop | locus | 16 | 0.6711 | 90 | 0.4248 | 95.6 |
| tier2 | ploop | window 1 | 17 | 0.6862 | 90 | 0.4248 | 98.9 |
| tier2 | ploop | window 2 | 16 | 0.6711 | 90 | 0.4248 | 95.6 |
| tier2 | ploop | window 3 | 17 | 0.6087 | 90 | 0.4248 | 92.2 |
| tier2 | ploop | window 4 | 17 | 0.7443 | 90 | 0.4248 | 100.0 |
| tier2 | ploop | window 5 | 19 | 0.4944 | 90 | 0.4248 | 70.0 |

**Reading.** In the four-repeat channels the positions that vary are more tree-congruent than the typical site — Cav repeats III/IV (EEDD) 95th/96th, Nav II/III 100th/97th, TPC 81st/92nd: selectivity at the four-repeat locus follows descent. The exceptions are positions whose minority state is carried by one or two chains (CatSper, NALCN repeat II). In the K channels the Y/F of GYG tracks descent (K2P, Kir, EAG, CNG ≥ 93rd) while the x of TxGYG is the homoplastic position (Kir 6th, Slo 3rd percentile; Shaker's T/S position 8th), i.e. re-tuned repeatedly within families.

## 5. Limits

* Members are S6's D39 sets (high-confidence calls, intact genome loci); medium calls and the r4 families are not in them.
* Six P-loop roots are undefined (S7d); their clade reads are unrooted.
* The tier-2 P-loop tree is a 90-column module tree; only its supported edges are read.
* Literature checks of the non-canonical strings (DEEA, DKEA, EKEE, DDDD) are *(pending: S14 references)*; the elephant-shark and opossum DEEA chains are an emergent row.
