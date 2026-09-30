# S3b — profile sweep of the panel proteomes + jackhmmer → census v3

_Rendered by `scripts/s3b_report.py` from the tables in this directory (D13). Do not edit by hand._

## 1. The sweep

S3a's **103 profiles**, unchanged (SHA-256 per profile in `sweep_runs.tsv`), searched against S4's declared denominator: **822,499 canonical entries from 50 reference proteomes** (release 2026_03, D35; DB SHA-256 `366dd895ad6b9a0c…`). `-Z` fixed to the DB size; 0.8 profile-hours. Assignment is S3a's D32 rule with its constants imported, not restated: ≥ 30 bits, ≥ 30% of the profile's match states, ≥ 10% relative margin over the best other profile.

**Instrument check.** 12,402 panel entries are also census v2 records, so S3a already scored the same sequence with the same profiles; only `-Z` differs, which moves domain i-E-values and so coverage. **Same verdict and family on 12,378 / 12,402 (99.8 %).** The differences:

| s3a_p_call | s3b_p_call | records |
|---|---|---|
| no_hit | module | 9 |
| no_hit | low_score | 9 |
| module | family | 4 |
| no_hit | family | 1 |
| low_score | module | 1 |

## 2. Census v3 on the panel — what domain search missed

The merge rule (fixed before the tables were read): an entry in census v2 keeps its v3a call; an entry outside it takes the S3b profile call (`panel_profile`); an entry only a jackhmmer run reaches is a `candidate`, never a family call (D14, D33).

**28,891 of 822,499 panel entries** carry evidence; 12,402 are census v2 records.

| v3_basis | v3_status | records |
|---|---|---|
| jackhmmer_only | candidate | 7980 |
| v3a:profile_only | channel | 4786 |
| panel_profile | non_channel_homolog | 4473 |
| v3a:both | channel | 3038 |
| panel_profile | channel_associated | 2960 |
| v3a:unassigned | unassigned | 1419 |
| v3a:superfamily_only | superfamily_only | 1234 |
| panel_profile | out_of_scope | 442 |
| v3a:profile_only | channel_contested | 423 |
| v3a:profile_only | non_channel_homolog | 382 |
| v3a:both | channel_contested | 355 |
| panel_profile | transporter | 308 |
| v3a:both | non_channel_homolog | 241 |
| panel_profile | channel | 190 |
| v3a:both | channel_associated | 175 |
| v3a:profile_only | transporter | 167 |
| v3a:s2_only | channel | 91 |
| panel_profile | superfamily_only | 73 |
| v3a:conflict | conflict | 72 |
| panel_profile | channel_contested | 63 |
| v3a:s2_only | channel_contested | 11 |
| v3a:s2_only | channel_associated | 8 |

**Channel-family calls: 8,957, of which 253 (2.8 %) are outside census v2** — members of a census family that carry none of the 67 enumerated pore signatures, so no domain search over those signatures could have found them. **Only 110 of them are high-confidence profile calls** (margin ≥ 30 % and ≥ half the profile). The other 143 are medium, and the largest medium blocks are **not channels**: the TRPN profile is mostly its ankyrin-repeat array and the LRRC8 profile half leucine-rich repeat, so an ankyrin- or LRR-repeat protein (ANKRD52, ankyrin, IκB, titin…) can cover 30 % of the profile with the repeat module alone and pass D32's coverage gate. D32 assumed a shared module is a small part of the profile; for repeat-dominated profiles it is not (emergent). Read the high column as the finding and the rest as an upper bound. Per family, largest first:

| family | superfamily | called | in_census_v2 | missed_by_domain_search | missed_high | missed_frac |
|---|---|---|---|---|---|---|
| trpn | ploop | 52 | 14 | 38 | 1 | 0.7308 |
| clic | clic | 123 | 91 | 32 | 27 | 0.2602 |
| hv1 | hv | 33 | 8 | 25 | 15 | 0.7576 |
| mitok | mitok | 23 | 0 | 23 | 14 | 1.0 |
| cng | ploop | 550 | 528 | 22 | 10 | 0.04 |
| lrrc8 | innexin_like | 128 | 108 | 20 | 0 | 0.1562 |
| trpv | ploop | 125 | 112 | 13 | 8 | 0.104 |
| trpa | ploop | 79 | 69 | 10 | 3 | 0.1266 |
| trpp | ploop | 108 | 99 | 9 | 0 | 0.0833 |
| pannexin | innexin_like | 43 | 35 | 8 | 8 | 0.186 |
| trpml | ploop | 66 | 58 | 8 | 3 | 0.1212 |
| plgic_prok | cysloop | 20 | 14 | 6 | 1 | 0.3 |
| tmc | tmem16_like | 155 | 149 | 6 | 4 | 0.0387 |
| osca_tmem63 | tmem16_like | 160 | 156 | 4 | 1 | 0.025 |
| kir | ploop | 292 | 289 | 3 | 2 | 0.0103 |
| trpc | ploop | 178 | 175 | 3 | 2 | 0.0169 |
| viroporin | viroporin | 3 | 0 | 3 | 3 | 1.0 |
| bestrophin | bestrophin | 118 | 116 | 2 | 2 | 0.0169 |
| nachr | cysloop | 825 | 823 | 2 | 0 | 0.0024 |
| iglur_nonvertebrate | iglur | 51 | 49 | 2 | 0 | 0.0392 |
| kainate | iglur | 168 | 166 | 2 | 0 | 0.0119 |
| hcn | ploop | 79 | 77 | 2 | 0 | 0.0253 |
| k2p | ploop | 408 | 406 | 2 | 1 | 0.0049 |
| kcsa_prok | ploop | 8 | 6 | 2 | 1 | 0.25 |
| kv_eag | ploop | 243 | 241 | 2 | 0 | 0.0082 |

Control families are reported apart: 8,183 of their 9,156 panel calls are outside census v2 — expected, since auxiliary subunits and non-channel homologues carry no pore signature by definition (D23).

By panel group (census families only):

| group | in_census_v2 | missed_by_domain_search | missed_high | missed |
|---|---|---|---|---|
| vertebrate | 4621 | 58 | 31 | 1.2 % |
| invertebrate | 1393 | 51 | 23 | 3.5 % |
| ciliate | 685 | 36 | 13 | 5.0 % |
| plant | 222 | 30 | 10 | 11.9 % |
| deuterostome | 848 | 25 | 9 | 2.9 % |
| basal_metazoan | 162 | 17 | 7 | 9.5 % |
| cnidarian | 520 | 16 | 5 | 3.0 % |
| holozoa | 72 | 6 | 2 | 7.7 % |
| prokaryote | 45 | 4 | 3 | 8.2 % |
| algae | 71 | 4 | 1 | 5.3 % |
| virus | 0 | 3 | 3 | 100.0 % |
| amoebozoa | 25 | 2 | 2 | 7.4 % |
| excavate | 10 | 1 | 1 | 9.1 % |
| fungi | 22 | 0 | 0 | 0.0 % |
| apicomplexa | 8 | 0 | 0 | 0.0 % |

Every missed record is listed in `missed_records.tsv`; the family × species counts S10 starts from are `family_by_species.tsv`.

## 3. jackhmmer to convergence (D10)

One run per census family from a **derived seed**: the panel entry with the highest-scoring high-confidence profile call to that family (`jackhmmer_seeds.tsv`, written before any run). `-N 10 -E 1e-5 --incE 1e-5` over the panel DB; 26.9 run-hours.

The kill criterion (`s3b_kill.py`, ported from the IP3R project with its constants): **K1** the share of included targets the profiles call to *another* family rises > 10% points above round 1; **K2** the included set grows > 10× (from ≥ 20); **K3** 10 rounds without converging. Rounds from the first firing on are excluded.

**68 runs: 24 clean, 44 killed (K1 23, K2 2, K3 19), 0 not run; 26 converged.**

| family | rule | killed_at | accepted_rounds | n_accepted | reason |
|---|---|---|---|---|---|
| kv_shaker | K1 | 2 | 1 | 820 | round 2: other-family share 73.0% (2128/2916), +36.1% on the round-1 baseline 36.8% (limit 10% points) |
| kv_modifier | K3 | 10 | 9 | 8603 | 10-round ceiling without converging; no completeness claim may rest on this run |
| kv_kcnq | K3 | 10 | 9 | 4411 | 10-round ceiling without converging; no completeness claim may rest on this run |
| kv_eag | K3 | 10 | 9 | 5181 | 10-round ceiling without converging; no completeness claim may rest on this run |
| kca_slo | K1 | 2 | 1 | 300 | round 2: other-family share 83.8% (1932/2307), +29.8% on the round-1 baseline 54.0% (limit 10% points) |
| kca_sk | K1 | 2 | 1 | 112 | round 2: other-family share 83.3% (1087/1305), +70.8% on the round-1 baseline 12.5% (limit 10% points) |
| k2p | K1 | 2 | 1 | 403 | round 2: other-family share 50.1% (525/1047), +50.1% on the round-1 baseline 0.0% (limit 10% points) |
| kcsa_prok | K2 | 2 | 1 | 128 | round 2 grew the included set to 2193, 17.13× the previous 128 (limit 10.0×) |
| nav | K1 | 2 | 1 | 581 | round 2: other-family share 80.5% (1164/1445), +14.3% on the round-1 baseline 66.3% (limit 10% points) |
| cav | K1 | 2 | 1 | 574 | round 2: other-family share 73.1% (1089/1490), +27.8% on the round-1 baseline 45.3% (limit 10% points) |
| nalcn | K3 | 10 | 9 | 8841 | 10-round ceiling without converging; no completeness claim may rest on this run |
| catsper | K3 | 10 | 9 | 3977 | 10-round ceiling without converging; no completeness claim may rest on this run |
| tpc | K1 | 2 | 1 | 441 | round 2: other-family share 84.0% (1240/1477), +11.6% on the round-1 baseline 72.3% (limit 10% points) |
| cng | K1 | 2 | 1 | 1028 | round 2: other-family share 40.0% (1118/2793), +11.3% on the round-1 baseline 28.7% (limit 10% points) |
| hcn | K3 | 10 | 9 | 4705 | 10-round ceiling without converging; no completeness claim may rest on this run |
| trpc | K3 | 10 | 9 | 10557 | 10-round ceiling without converging; no completeness claim may rest on this run |
| trpv | K2 | 2 | 1 | 267 | round 2 grew the included set to 5852, 21.92× the previous 267 (limit 10.0×) |
| trpm | K1 | 2 | 1 | 376 | round 2: other-family share 43.5% (358/823), +29.4% on the round-1 baseline 14.1% (limit 10% points) |
| trpa | K1 | 3 | 2 | 7919 | round 3: other-family share 16.0% (1436/8962), +13.6% on the round-1 baseline 2.4% (limit 10% points) |
| trpml | K1 | 2 | 1 | 84 | round 2: other-family share 75.4% (402/533), +57.6% on the round-1 baseline 17.9% (limit 10% points) |
| trpp | K1 | 2 | 1 | 416 | round 2: other-family share 75.2% (955/1270), +22.3% on the round-1 baseline 52.9% (limit 10% points) |
| trpn | K1 | 4 | 3 | 8503 | round 4: other-family share 20.8% (1965/9424), +12.1% on the round-1 baseline 8.8% (limit 10% points) |
| glyr | K3 | 10 | 9 | 1623 | 10-round ceiling without converging; no completeness claim may rest on this run |
| plgic_prok | K3 | 10 | 9 | 1613 | 10-round ceiling without converging; no completeness claim may rest on this run |
| ampa | K1 | 6 | 5 | 2734 | round 6: other-family share 71.7% (2326/3244), +13.9% on the round-1 baseline 57.8% (limit 10% points) |
| kainate | K1 | 2 | 1 | 714 | round 2: other-family share 65.5% (1407/2147), +19.3% on the round-1 baseline 46.2% (limit 10% points) |
| nmda | K1 | 3 | 2 | 1192 | round 3: other-family share 70.7% (1604/2269), +16.5% on the round-1 baseline 54.1% (limit 10% points) |
| delta_glur | K1 | 6 | 5 | 2693 | round 6: other-family share 71.6% (2245/3137), +12.4% on the round-1 baseline 59.2% (limit 10% points) |
| iglur_nonvertebrate | K3 | 10 | 9 | 4146 | 10-round ceiling without converging; no completeness claim may rest on this run |
| iglur_prok | K3 | 10 | 9 | 3670 | 10-round ceiling without converging; no completeness claim may rest on this run |
| clc_channel | K3 | 10 | 9 | 579 | 10-round ceiling without converging; no completeness claim may rest on this run |
| cftr | K3 | 10 | 9 | 3504 | 10-round ceiling without converging; no completeness claim may rest on this run |
| tweety | K3 | 10 | 9 | 179 | 10-round ceiling without converging; no completeness claim may rest on this run |
| ano_channel | K3 | 10 | 9 | 484 | 10-round ceiling without converging; no completeness claim may rest on this run |
| ano_scramblase | K1 | 8 | 7 | 414 | round 8: other-family share 22.7% (102/449), +12.0% on the round-1 baseline 10.7% (limit 10% points) |
| osca_tmem63 | K1 | 4 | 3 | 268 | round 4: other-family share 50.6% (278/549), +50.6% on the round-1 baseline 0.0% (limit 10% points) |
| piezo | K3 | 10 | 9 | 760 | 10-round ceiling without converging; no completeness claim may rest on this run |
| pannexin | K1 | 2 | 1 | 42 | round 2: other-family share 30.4% (21/69), +30.4% on the round-1 baseline 0.0% (limit 10% points) |
| innexin | K1 | 2 | 1 | 80 | round 2: other-family share 27.3% (39/143), +26.0% on the round-1 baseline 1.2% (limit 10% points) |
| lrrc8 | K3 | 10 | 9 | 9025 | 10-round ceiling without converging; no completeness claim may rest on this run |
| itpr | K1 | 3 | 2 | 341 | round 3: other-family share 55.9% (430/769), +38.8% on the round-1 baseline 17.1% (limit 10% points) |
| ryr | K3 | 10 | 9 | 6378 | 10-round ceiling without converging; no completeness claim may rest on this run |
| hv1 | K1 | 2 | 1 | 80 | round 2: other-family share 83.8% (656/783), +60.0% on the round-1 baseline 23.8% (limit 10% points) |
| clic | K3 | 10 | 9 | 1565 | 10-round ceiling without converging; no completeness claim may rest on this run |

## 4. The completeness argument

Two directions. **Does a single-sequence iterated search reach what the profile reaches?** Recall of each family's panel profile calls by its own accepted jackhmmer set: **8,610 / 8,736 (98.6 %)**. **Does it reach what the profile cannot call?** Accepted targets the profile library calls to no family are the upper bound on what the sweep misses; those outside census v2 are the `candidate` rows, by profile verdict: module 11,452, no_hit 2,805, low_score 1,481.

**The claim rests on the 24 clean runs**, which recover **3,352 / 3,356 (99.9 %)** of their families' profile calls and include only **1,460** targets the profile library calls to no family (at most, what the sweep misses in those families). Killed runs are reported with their pre-kill rounds; a K3 run contributes no candidates (D36).

Candidate rows by the run that found them (a target in several runs counts once per run): trpn 7,076, trpa 7,052, ampa 483, delta_glur 476, trpm 122, nmda 92. The TRPA and TRPN runs dominate: their pre-K1 rounds already hold thousands of ankyrin-repeat proteins the profiles call `module` or not at all. **K1 cannot see this drift** — uncalled targets are exempt by design, as in the parent project, whose S19 measured the same blind spot. The candidates are counted, never called.

Families whose jackhmmer run recovers < 90 % of their profile calls (3):

| family | verdict | profile_calls | recovered | recall | jh_accepted | jh_other_family | jh_own_superfamily_only | jh_module | jh_low_score | jh_no_hit | top_other_families |
|---|---|---|---|---|---|---|---|---|---|---|---|
| viroporin | clean | 3 | 1 | 0.3333 | 1 | 0 | 0 | 0 | 0 | 0 |  |
| kcsa_prok | killed | 8 | 4 | 0.5 | 128 | 107 | 2 | 7 | 2 | 6 | k2p:73,kca_slo:29,cng:5 |
| cng | killed | 553 | 468 | 0.8463 | 1028 | 295 | 117 | 144 | 0 | 4 | kv_eag:219,hcn:76 |

All families: `jackhmmer_completeness.tsv`.

## 5. Human positive control

The 328 human census genes against the human reference proteome (gene symbol used to *score* the call, never to make it — H15): **right_family 327**, **not_in_proteome 1**.

| gene | family | target | v3_family | v3_status | v3_basis | p_call | p_family | verdict |
|---|---|---|---|---|---|---|---|---|
| GLRA4 | glyr |  |  |  |  |  |  | not_in_proteome |

