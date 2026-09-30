# S2 — census v2: every UniProtKB protein carrying a pore signature

*Rendered from the tables in this directory by `scripts/s2_report.py`. UniProtKB release 2026_03. Revision r3 (S2c, 2026-09-28: H2/H4/H11/H12/H13 positive tests) — earlier revisions are in git history.*

## Headline

- **1,275,990 records** carry at least one of the catalogue's 79 pore signatures (union count 1,275,990; fetched 1,275,990).
- **314,205 (24.6 %) get a family call**, 268,583 of them to a family catalogued as a channel; 437,608 (34.3 %) reach a superfamily only; 524,177 (41.1 %) are `unassigned` and stay in the census with that label.
- **Every family call was made without the reference tier**: 282,837 by architecture or hazard rules and 31,368 by the selectivity-filter motif (37,393 four-repeat records projected onto Nav1.5).
- **Human census genes: 327/329 enumerated, 175/329 called to the right family** (reviewed human entries; the gene symbol is used here to *score*, never to classify — H15).
- **S1 panel: 75/97 enumerated; 45 correct without the reference tier vs 57 with it in S1.** Panel members from channel families enumerated: **71/71**. Non-channel members: 22/26 carry no pore signature and are excluded at enumeration; of the 4 that are enumerated, **0 are called to a channel family** (CLCN7 → superfamily_only, PKD1 → assoc_polycystin1, KCTD1 → nonchannel_kctd, TPTE → nonchannel_vsp).

## Completeness

The union query was cut into taxonomic shards that partition it; the shard counts must sum to the union count, and each shard's fetched records must equal UniProt's own count for it.

| shard | uniprot_count | fetched | status |
|---|---|---|---|
| bact_pseudomonadota | 213842 | 213842 | ok |
| bact_other | 229741 | 229741 | ok |
| archaea | 10281 | 10281 | ok |
| virus | 177 | 177 | ok |
| mammalia | 106776 | 106776 | ok |
| actinopterygii | 164799 | 164799 | ok |
| chordata_other | 165576 | 165576 | ok |
| metazoa_other | 194499 | 194499 | ok |
| viridiplantae | 74728 | 74728 | ok |
| fungi | 48186 | 48186 | ok |
| euk_other | 36307 | 36307 | ok |
| unplaced | 288 | 288 | ok |
| r4_bact_pseudomonadota | 1 | 1 | ok |
| r4_bact_other | 4 | 4 | ok |
| r4_archaea | 626 | 626 | ok |
| r4_virus | 4 | 4 | ok |
| r4_mammalia | 3249 | 3249 | ok |
| r4_actinopterygii | 3426 | 3426 | ok |
| r4_chordata_other | 4094 | 4094 | ok |
| r4_metazoa_other | 3548 | 3548 | ok |
| r4_viridiplantae | 6376 | 6376 | ok |
| r4_fungi | 4398 | 4398 | ok |
| r4_euk_other | 1039 | 1039 | ok |
| r4_unplaced | 18 | 18 | ok |
| r5_bact_pseudomonadota | 1 | 1 | ok |
| r5_bact_other | 9 | 9 | ok |
| r5_archaea | 0 | 0 | ok |
| r5_virus | 3994 | 3994 | ok |
| r5_mammalia | 0 | 0 | ok |
| r5_actinopterygii | 0 | 0 | ok |
| r5_chordata_other | 0 | 0 | ok |
| r5_metazoa_other | 1 | 1 | ok |
| r5_viridiplantae | 0 | 0 | ok |
| r5_fungi | 1 | 1 | ok |
| r5_euk_other | 1 | 1 | ok |
| r5_unplaced | 0 | 0 | ok |
| UNION | 1275990 | 1275990 | ok |

Shards or union short: **0**. Per-signature check (records fetched carrying the signature == UniProt's count for it alone): **79/79 ok**.

UniProt and InterPro run on different release cycles, so their counts for the same Pfam entry differ. The census is UniProt's; the largest differences are listed so the gap is visible rather than assumed away.

| pfam | uniprot_count | interpro_count | uniprot_minus_interpro |
|---|---|---|---|
| PF00520 | 207571 | 206358 | 1213 |
| PF02931 | 83437 | 84444 | -1007 |
| PF00654 | 47664 | 48464 | -800 |
| PF00497 | 118379 | 117979 | 400 |
| PF00027 | 185439 | 185042 | 397 |
| PF07885 | 68992 | 68610 | 382 |
| PF21082 | 37900 | 37581 | 319 |
| PF00060 | 56382 | 56688 | -306 |
| PF16178 | 13901 | 14184 | -283 |
| PF00924 | 69768 | 69514 | 254 |

## Calls by family — channels

| family | superfamily | records | reviewed | human | eukaryota | bacteria | archaea | viruses | gold | silver | bronze |
|---|---|---|---|---|---|---|---|---|---|---|---|
| k2p | ploop | 26007 | 61 | 38 | 26005 | 2 | 0 | 0 | 0 | 26007 | 0 |
| cav | ploop | 20080 | 43 | 77 | 20080 | 0 | 0 | 0 | 5180 | 14085 | 815 |
| connexin | connexin | 19855 | 111 | 52 | 19855 | 0 | 0 | 0 | 0 | 0 | 19855 |
| osca_tmem63 | tmem16_like | 18520 | 35 | 16 | 18520 | 0 | 0 | 0 | 0 | 0 | 18520 |
| kir | ploop | 17515 | 68 | 50 | 17464 | 51 | 0 | 0 | 0 | 17515 | 0 |
| vdac | porin | 17059 | 68 | 53 | 17059 | 0 | 0 | 0 | 0 | 0 | 17059 |
| mscl | msc | 16591 | 385 | 0 | 1049 | 15494 | 46 | 2 | 0 | 0 | 16591 |
| kca_slo | ploop | 13003 | 24 | 79 | 13003 | 0 | 0 | 0 | 0 | 0 | 13003 |
| nav | ploop | 12204 | 39 | 35 | 12204 | 0 | 0 | 0 | 6569 | 4175 | 1460 |
| kv_kcnq | ploop | 8842 | 21 | 45 | 8842 | 0 | 0 | 0 | 0 | 0 | 8842 |
| bestrophin | bestrophin | 8233 | 25 | 29 | 8233 | 0 | 0 | 0 | 0 | 0 | 8233 |
| kv_eag | ploop | 7492 | 27 | 20 | 7485 | 7 | 0 | 0 | 0 | 0 | 7492 |
| tmem175 | tmem175 | 6987 | 10 | 15 | 1131 | 5672 | 184 | 0 | 0 | 0 | 6987 |
| mcu | mcu | 6471 | 20 | 7 | 6468 | 3 | 0 | 0 | 0 | 0 | 6471 |
| trpc | ploop | 6221 | 17 | 20 | 6221 | 0 | 0 | 0 | 0 | 0 | 6221 |
| otop | otopetrin | 6107 | 10 | 7 | 6107 | 0 | 0 | 0 | 0 | 0 | 6107 |
| mscs | msc | 5992 | 10 | 0 | 8 | 5950 | 34 | 0 | 0 | 0 | 5992 |
| ryr | ca_release | 5618 | 13 | 20 | 5618 | 0 | 0 | 0 | 0 | 5297 | 321 |
| calhm | calhm | 4946 | 19 | 8 | 4946 | 0 | 0 | 0 | 0 | 0 | 4946 |
| cng | ploop | 4713 | 22 | 10 | 4713 | 0 | 0 | 0 | 0 | 0 | 4713 |
| trpml | ploop | 4692 | 8 | 32 | 4692 | 0 | 0 | 0 | 0 | 3865 | 827 |
| kca_sk | ploop | 4553 | 15 | 25 | 4553 | 0 | 0 | 0 | 0 | 0 | 4553 |
| lrrc8 | innexin_like | 4357 | 17 | 29 | 4356 | 1 | 0 | 0 | 0 | 0 | 4357 |
| orai | orai | 3586 | 16 | 11 | 3586 | 0 | 0 | 0 | 0 | 0 | 3586 |
| hcn | ploop | 3540 | 14 | 13 | 3540 | 0 | 0 | 0 | 0 | 0 | 3540 |
| piezo | piezo | 3513 | 7 | 7 | 3513 | 0 | 0 | 0 | 0 | 0 | 3513 |
| tric | tric | 2673 | 19 | 16 | 2672 | 1 | 0 | 0 | 0 | 0 | 2673 |
| trpm | ploop | 2280 | 6 | 3 | 2280 | 0 | 0 | 0 | 0 | 0 | 2280 |
| cftr | abc_channel | 1469 | 42 | 16 | 1469 | 0 | 0 | 0 | 0 | 1445 | 24 |
| nalcn | ploop | 1359 | 3 | 8 | 1359 | 0 | 0 | 0 | 0 | 1359 | 0 |
| trpp | ploop | 1358 | 7 | 8 | 1358 | 0 | 0 | 0 | 0 | 1358 | 0 |
| pacc | pac | 992 | 8 | 13 | 992 | 0 | 0 | 0 | 0 | 0 | 992 |
| kv_shaker | ploop | 932 | 11 | 3 | 932 | 0 | 0 | 0 | 0 | 0 | 932 |
| hv1 | hv | 685 | 3 | 5 | 685 | 0 | 0 | 0 | 0 | 685 | 0 |
| iglur_prok | iglur | 138 | 0 | 0 | 2 | 136 | 0 | 0 | 0 | 0 | 138 |

## Calls by family — catalogued non-channels (controls, auxiliaries, transporters, out of scope)

These records carry a pore signature and were positively called to a family the catalogue holds *in order to exclude it*. They are in the census so the exclusion is visible.

| family | superfamily | catalogue_status | records | reviewed | eukaryota | bacteria | archaea | viruses |
|---|---|---|---|---|---|---|---|---|
| nonchannel_kctd | ploop | non_channel_homolog | 14354 | 45 | 14354 | 0 | 0 | 0 |
| tmc | tmem16_like | channel_contested | 8494 | 27 | 8494 | 0 | 0 | 0 |
| clic | clic | channel_contested | 6354 | 25 | 6353 | 1 | 0 | 0 |
| tweety | tweety | channel_contested | 4678 | 20 | 4677 | 1 | 0 | 0 |
| assoc_polycystin1 | ploop | channel_associated | 4137 | 11 | 4137 | 0 | 0 | 0 |
| gphr | gphr | channel_contested | 3614 | 14 | 3614 | 0 | 0 | 0 |
| tmem87 | tmem87 | channel_contested | 1660 | 5 | 1660 | 0 | 0 | 0 |
| clcc1 | mclc | channel_contested | 1562 | 7 | 1558 | 0 | 0 | 4 |
| nonchannel_vsp | hv | non_channel_homolog | 767 | 3 | 767 | 0 | 0 | 0 |
| nonchannel_class_c_gpcr | iglur | non_channel_homolog | 2 | 0 | 2 | 0 | 0 | 0 |

## Superfamily only

`ambiguous` counts the records whose architecture named two or more candidate families and could not choose (D7: a tie is not a call); the rest matched superfamily-level evidence only.

| superfamily | records | ambiguous | reviewed | eukaryota | bacteria | archaea | viruses |
|---|---|---|---|---|---|---|---|
| ploop | 183801 | 28028 | 327 | 151339 | 31126 | 1232 | 103 |
| cysloop | 83413 | 0 | 262 | 83160 | 245 | 8 | 0 |
| iglur | 56369 | 0 | 127 | 55605 | 763 | 0 | 1 |
| clc | 47655 | 0 | 192 | 27507 | 19886 | 261 | 0 |
| deg_enac | 25927 | 25927 | 78 | 25927 | 0 | 0 | 0 |
| tmem16_like | 13287 | 13287 | 17 | 13287 | 0 | 0 | 0 |
| innexin_like | 11653 | 0 | 40 | 11651 | 0 | 0 | 2 |
| p2x | 8295 | 8295 | 29 | 8295 | 0 | 0 | 0 |
| ca_release | 7208 | 0 | 15 | 7208 | 0 | 0 | 0 |

## Which tier decided

| decisive_tier | status | records |
|---|---|---|
| architecture | channel | 181043 |
| hazard | channel | 56172 |
| motif | channel | 31368 |
| architecture | channel_contested | 26362 |
| hazard | non_channel_homolog | 15123 |
| hazard | channel_associated | 4137 |

## Four-repeat selectivity filters

| filter_string | family | records |
|---|---|---|
| EEEE | cav | 13966 |
| DEKA | nav | 10744 |
| EEDD | cav | 5299 |
| EEKE | nalcn | 1359 |
| NEEE | - | 1267 |
| EEE- | - | 438 |
| DEEA | - | 435 |
| DEKG | nav | 253 |
| DE-A | nav | 231 |
| QEEE | - | 221 |
| DEK- | nav | 178 |
| -EEE | - | 166 |
| DENA | - | 164 |
| DDDD | - | 129 |
| QEEE | cav | 94 |
| EED- | - | 93 |
| -EDD | - | 88 |
| EKEE | - | 73 |
| DE-A | - | 72 |
| EDEE | - | 54 |
| EDDD | - | 51 |
| -EEE | cav | 43 |
| E-EE | - | 40 |
| EE-E | - | 40 |
| DEKG | - | 37 |

## What brought the unassigned records in

`records` counts every unassigned record carrying the signature; `solo_records` those for which it is the only pore signature.

| pfam | records | solo_records | eukaryota | bacteria | archaea | viruses |
|---|---|---|---|---|---|---|
| PF00027 | 150616 | 149130 | 44306 | 106294 | 16 | 0 |
| PF00497 | 117066 | 115927 | 982 | 115648 | 434 | 2 |
| PF13426 | 112924 | 111788 | 18186 | 89759 | 4977 | 1 |
| PF00924 | 63770 | 29372 | 9342 | 51623 | 2804 | 1 |
| PF21082 | 31906 | 143 | 92 | 30000 | 1814 | 0 |
| PF23317 | 12289 | 10388 | 12289 | 0 | 0 | 0 |
| PF04547 | 11053 | 10965 | 11053 | 0 | 0 | 0 |
| PF06814 | 9565 | 9563 | 9565 | 0 | 0 | 0 |
| PF02932 | 7307 | 7294 | 7282 | 21 | 4 | 0 |
| PF01956 | 7252 | 7250 | 6625 | 1 | 626 | 0 |
| PF18139 | 4720 | 1954 | 4720 | 0 | 0 | 0 |
| PF05552 | 4569 | 3083 | 4 | 4042 | 523 | 0 |
| PF10613 | 3852 | 3838 | 3847 | 5 | 0 | 0 |
| PF25508 | 3682 | 643 | 3682 | 0 | 0 | 0 |
| PF02714 | 3046 | 1293 | 3046 | 0 | 0 | 0 |
| PF14703 | 2301 | 598 | 2301 | 0 | 0 | 0 |
| PF12166 | 2257 | 329 | 2257 | 0 | 0 | 0 |
| PF00599 | 2223 | 2223 | 0 | 4 | 0 | 2219 |
| PF24874 | 2201 | 91 | 2201 | 0 | 0 | 0 |
| PF22614 | 2058 | 2056 | 1911 | 147 | 0 | 0 |
| PF23188 | 1802 | 107 | 1802 | 0 | 0 | 0 |
| PF14965 | 1691 | 1691 | 1691 | 0 | 0 | 0 |
| PF01365 | 1459 | 1156 | 1458 | 1 | 0 | 0 |
| PF00558 | 1386 | 1386 | 3 | 6 | 0 | 1377 |
| PF13967 | 1331 | 528 | 1331 | 0 | 0 | 0 |
| PF02026 | 1326 | 1074 | 678 | 570 | 13 | 65 |
| PF24871 | 1125 | 766 | 1125 | 0 | 0 | 0 |
| PF17655 | 1074 | 1069 | 1046 | 28 | 0 | 0 |
| PF08344 | 933 | 933 | 933 | 0 | 0 | 0 |
| PF15917 | 788 | 252 | 788 | 0 | 0 | 0 |

## Unassigned records carrying part of one family's architecture

A derived family rule requires *every* `FAMILY`-level signature the catalogue declares for the family (`rules.derived_rules`). A record carrying some but not all of one family's signatures therefore goes unassigned. This table counts them: a measure of how far the declared architectures are from what the family's members actually carry, and the input to any rule revision — which must be re-benchmarked on S1 before it is adopted, not tuned here.

| family | declared_family_signatures | unassigned_records |
|---|---|---|
| (several) |  | 175715 |
| iglur_prok | PF00497,PF07885 | 115928 |
| kv_eag | PF00027,PF13426 | 111789 |
| mscs | PF00924,PF05552,PF21082 | 65508 |
| (none) |  | 35201 |
| osca_tmem63 | PF02714,PF13967,PF14703 | 4505 |
| trpm | PF16519,PF18139,PF23317,PF25508 | 3960 |
| piezo | PF12166,PF15917,PF23188,PF24871,PF24874 | 3939 |
| kca_slo | PF03493,PF22614 | 2243 |
| ryr | PF02026,PF06459 | 1497 |
| kir | PF01007,PF08466,PF17655 | 1076 |
| kca_sk | PF02888,PF03530 | 981 |
| trpc | PF08344,PF23317 | 933 |
| gphr | PF12430,PF12537 | 451 |
| hcn | PF00027,PF08412 | 291 |
| nav | PF06512,PF11933 | 99 |
| cng | PF00027,PF16526 | 61 |

## Human census genes

Not enumerated (2): the reviewed human entry carries none of the pore signatures, or has no reviewed entry under that symbol.

| gene | expected_family |
|---|---|
| GLRA4 | glyr |
| CCDC51 | mitok |

Enumerated but not called to the expected family (152):

| gene | expected_family | accession | s2_family | s2_status | confidence |
|---|---|---|---|---|---|
| GRIA1 | ampa | P42261 |  | superfamily_only | bronze |
| GRIA2 | ampa | P42262 |  | superfamily_only | bronze |
| GRIA3 | ampa | P42263 |  | superfamily_only | bronze |
| GRIA4 | ampa | P48058 |  | superfamily_only | bronze |
| ANO1 | ano_channel | Q5XXA6 |  | superfamily_only | bronze |
| ANO2 | ano_channel | Q9NQ90 |  | superfamily_only | bronze |
| ANO10 | ano_scramblase | Q9NW15 |  | unassigned | unassigned |
| ANO3 | ano_scramblase | Q9BYT9 |  | superfamily_only | bronze |
| ANO4 | ano_scramblase | Q32M45 |  | superfamily_only | bronze |
| ANO5 | ano_scramblase | Q75V66 |  | superfamily_only | bronze |
| ANO6 | ano_scramblase | Q4KMQ2 |  | superfamily_only | bronze |
| ANO7 | ano_scramblase | Q6IWH7 |  | superfamily_only | bronze |
| ANO8 | ano_scramblase | Q9HCE9 |  | unassigned | unassigned |
| ANO9 | ano_scramblase | A1A5B4 |  | superfamily_only | bronze |
| ASIC1 | asic | P78348 |  | superfamily_only | bronze |
| ASIC2 | asic | Q16515 |  | superfamily_only | bronze |
| ASIC3 | asic | Q9UHC3 |  | superfamily_only | bronze |
| ASIC4 | asic | Q96FT7 |  | superfamily_only | bronze |
| ASIC5 | asic | Q9NY37 |  | superfamily_only | bronze |
| CATSPER1 | catsper | Q8NEC5 |  | superfamily_only | bronze |
| CATSPER2 | catsper | Q96P56 |  | superfamily_only | bronze |
| CATSPER3 | catsper | Q86XQ3 |  | superfamily_only | bronze |
| CATSPER4 | catsper | Q7RTX7 |  | superfamily_only | bronze |
| CLCN1 | clc_channel | P35523 |  | superfamily_only | bronze |
| CLCN2 | clc_channel | P51788 |  | superfamily_only | bronze |
| CLCNKA | clc_channel | P51800 |  | superfamily_only | bronze |
| CLCNKB | clc_channel | P51801 |  | superfamily_only | bronze |
| CNGB1 | cng | Q14028 |  | superfamily_only | bronze |
| CNGB3 | cng | Q9NQW8 |  | superfamily_only | bronze |
| GRID1 | delta_glur | Q9ULK0 |  | superfamily_only | bronze |
| GRID2 | delta_glur | O43424 |  | superfamily_only | bronze |
| SCNN1A | enac | P37088 |  | superfamily_only | bronze |
| SCNN1B | enac | P51168 |  | superfamily_only | bronze |
| SCNN1D | enac | P51172 |  | superfamily_only | bronze |
| SCNN1G | enac | P51170 |  | superfamily_only | bronze |
| GABRA1 | gabaa | P14867 |  | superfamily_only | bronze |
| GABRA2 | gabaa | P47869 |  | superfamily_only | bronze |
| GABRA3 | gabaa | P34903 |  | superfamily_only | bronze |
| GABRA4 | gabaa | P48169 |  | superfamily_only | bronze |
| GABRA5 | gabaa | P31644 |  | superfamily_only | bronze |
| GABRA6 | gabaa | Q16445 |  | superfamily_only | bronze |
| GABRB1 | gabaa | P18505 |  | superfamily_only | bronze |
| GABRB2 | gabaa | P47870 |  | superfamily_only | bronze |
| GABRB3 | gabaa | P28472 |  | superfamily_only | bronze |
| GABRD | gabaa | O14764 |  | superfamily_only | bronze |
| GABRE | gabaa | P78334 |  | superfamily_only | bronze |
| GABRG1 | gabaa | Q8N1C3 |  | superfamily_only | bronze |
| GABRG2 | gabaa | P18507 |  | superfamily_only | bronze |
| GABRG3 | gabaa | Q99928 |  | superfamily_only | bronze |
| GABRP | gabaa | O00591 |  | superfamily_only | bronze |
| GABRQ | gabaa | Q9UN88 |  | superfamily_only | bronze |
| GABRR1 | gabaa | P24046 |  | superfamily_only | bronze |
| GABRR2 | gabaa | P28476 |  | superfamily_only | bronze |
| GABRR3 | gabaa | A8MPY1 |  | superfamily_only | bronze |
| GLRA1 | glyr | P23415 |  | superfamily_only | bronze |
| GLRA2 | glyr | P23416 |  | superfamily_only | bronze |
| GLRA3 | glyr | O75311 |  | superfamily_only | bronze |
| GLRB | glyr | P48167 |  | superfamily_only | bronze |
| HTR3A | ht3 | P46098 |  | superfamily_only | bronze |
| HTR3B | ht3 | O95264 |  | superfamily_only | bronze |
| HTR3C | ht3 | Q8WXA8 |  | superfamily_only | bronze |
| HTR3D | ht3 | Q70Z44 |  | superfamily_only | bronze |
| HTR3E | ht3 | A5X5Y0 |  | superfamily_only | bronze |
| ITPR1 | itpr | Q14643 |  | superfamily_only | bronze |
| ITPR2 | itpr | Q14571 |  | superfamily_only | bronze |
| ITPR3 | itpr | Q14573 |  | superfamily_only | bronze |
| GRIK1 | kainate | P39086 |  | superfamily_only | bronze |
| GRIK2 | kainate | Q13002 |  | superfamily_only | bronze |
| GRIK3 | kainate | Q13003 |  | superfamily_only | bronze |
| GRIK4 | kainate | Q16099 |  | superfamily_only | bronze |
| GRIK5 | kainate | Q16478 |  | superfamily_only | bronze |
| KCNH4 | kv_eag | Q9UQ05 |  | superfamily_only | bronze |
| KCNF1 | kv_modifier | Q9H3M0 |  | superfamily_only | bronze |
| KCNG1 | kv_modifier | Q9UIX4 |  | superfamily_only | bronze |
| KCNG2 | kv_modifier | Q9UJ96 |  | superfamily_only | bronze |
| KCNG3 | kv_modifier | Q8TAE7 |  | superfamily_only | bronze |
| KCNG4 | kv_modifier | Q8TDN1 |  | superfamily_only | bronze |
| KCNS1 | kv_modifier | Q96KK3 |  | superfamily_only | bronze |
| KCNS2 | kv_modifier | Q9ULS6 |  | superfamily_only | bronze |
| KCNS3 | kv_modifier | Q9BQ31 |  | superfamily_only | bronze |
| KCNV1 | kv_modifier | Q6PIU1 |  | superfamily_only | bronze |
| KCNV2 | kv_modifier | Q8TDN2 |  | superfamily_only | bronze |
| KCNA1 | kv_shaker | Q09470 |  | superfamily_only | bronze |
| KCNA10 | kv_shaker | Q16322 |  | superfamily_only | bronze |
| KCNA2 | kv_shaker | P16389 |  | superfamily_only | bronze |
| KCNA3 | kv_shaker | P22001 |  | superfamily_only | bronze |
| KCNA4 | kv_shaker | P22459 |  | superfamily_only | bronze |
| KCNA5 | kv_shaker | P22460 |  | superfamily_only | bronze |
| KCNA6 | kv_shaker | P17658 |  | superfamily_only | bronze |
| KCNA7 | kv_shaker | Q96RP8 |  | superfamily_only | bronze |
| KCNC1 | kv_shaker | P48547 |  | superfamily_only | bronze |
| KCNC2 | kv_shaker | Q96PR1 |  | superfamily_only | bronze |
| KCNC3 | kv_shaker | Q14003 |  | superfamily_only | bronze |
| KCNC4 | kv_shaker | Q03721 |  | superfamily_only | bronze |
| KCND1 | kv_shaker | Q9NSA2 |  | superfamily_only | bronze |
| KCND2 | kv_shaker | Q9NZV8 |  | superfamily_only | bronze |
| KCND3 | kv_shaker | Q9UK17 |  | superfamily_only | bronze |
| CHRNA1 | nachr | P02708 |  | superfamily_only | bronze |
| CHRNA10 | nachr | Q9GZZ6 |  | superfamily_only | bronze |
| CHRNA2 | nachr | Q15822 |  | superfamily_only | bronze |
| CHRNA3 | nachr | P32297 |  | superfamily_only | bronze |
| CHRNA4 | nachr | P43681 |  | superfamily_only | bronze |
| CHRNA5 | nachr | P30532 |  | superfamily_only | bronze |
| CHRNA6 | nachr | Q15825 |  | superfamily_only | bronze |
| CHRNA7 | nachr | P36544 |  | superfamily_only | bronze |
| CHRNA9 | nachr | Q9UGM1 |  | superfamily_only | bronze |
| CHRNB1 | nachr | P11230 |  | superfamily_only | bronze |
| CHRNB2 | nachr | P17787 |  | superfamily_only | bronze |
| CHRNB3 | nachr | Q05901 |  | superfamily_only | bronze |
| CHRNB4 | nachr | P30926 |  | superfamily_only | bronze |
| CHRND | nachr | Q07001 |  | superfamily_only | bronze |
| CHRNE | nachr | Q04844 |  | superfamily_only | bronze |
| CHRNG | nachr | P07510 |  | superfamily_only | bronze |
| SCN7A | nav | Q01118 |  | superfamily_only | bronze |
| GRIN1 | nmda | Q05586 |  | superfamily_only | bronze |
| GRIN2A | nmda | Q12879 |  | superfamily_only | bronze |
| GRIN2B | nmda | Q13224 |  | superfamily_only | bronze |
| GRIN2C | nmda | Q14957 |  | superfamily_only | bronze |
| GRIN2D | nmda | O15399 |  | superfamily_only | bronze |
| GRIN3A | nmda | Q8TCU5 |  | superfamily_only | bronze |
| GRIN3B | nmda | O60391 |  | superfamily_only | bronze |
| P2RX1 | p2x | P51575 |  | superfamily_only | bronze |
| P2RX2 | p2x | Q9UBL9 |  | superfamily_only | bronze |
| P2RX3 | p2x | P56373 |  | superfamily_only | bronze |
| P2RX4 | p2x | Q99571 |  | superfamily_only | bronze |
| P2RX5 | p2x | Q93086 |  | superfamily_only | bronze |
| P2RX6 | p2x | O15547 |  | superfamily_only | bronze |
| P2RX7 | p2x | Q99572 |  | superfamily_only | bronze |
| PANX1 | pannexin | Q96RD7 |  | superfamily_only | bronze |
| PANX2 | pannexin | Q96RD6 |  | superfamily_only | bronze |
| PANX3 | pannexin | Q96QZ0 |  | superfamily_only | bronze |
| TMCO1 | tmco1 | Q9UM00 |  | unassigned | unassigned |
| TMEM109 | tmem109 | Q9BVC6 |  | unassigned | unassigned |
| TPCN1 | tpc | Q9ULQ1 |  | superfamily_only | bronze |
| TPCN2 | tpc | Q8NHX9 |  | superfamily_only | bronze |
| TRPA1 | trpa | O75762 |  | superfamily_only | bronze |
| TRPC4 | trpc | Q9UBN4 |  | superfamily_only | bronze |
| TRPC5 | trpc | Q9UL62 |  | superfamily_only | bronze |
| TRPM1 | trpm | Q7Z4N2 |  | superfamily_only | bronze |
| TRPM2 | trpm | O94759 |  | superfamily_only | bronze |
| TRPM3 | trpm | Q9HCF6 |  | superfamily_only | bronze |
| TRPM4 | trpm | Q8TD43 |  | superfamily_only | bronze |
| TRPM5 | trpm | Q9NZQ8 |  | superfamily_only | bronze |
| TRPM8 | trpm | Q7Z2W7 |  | unassigned | unassigned |
| PKD2L2 | trpp | Q9NZM6 |  | superfamily_only | bronze |
| TRPV1 | trpv | Q8NER1 |  | superfamily_only | bronze |
| TRPV2 | trpv | Q9Y5S1 |  | superfamily_only | bronze |
| TRPV3 | trpv | Q8NET8 |  | unassigned | unassigned |
| TRPV4 | trpv | Q9HBA0 |  | superfamily_only | bronze |
| TRPV5 | trpv | Q9NQA5 |  | superfamily_only | bronze |
| TRPV6 | trpv | Q9H1D0 |  | superfamily_only | bronze |
| ZACN | zac | Q401N2 |  | superfamily_only | bronze |

## S1 panel: S1 call vs S2 call

| accession | gene | expected_family | s1_family | s2_family | s2_status | s1_correct | s2_correct |
|---|---|---|---|---|---|---|---|
| Q09470 | KCNA1 | kv_shaker | kv_shaker |  | superfamily_only | yes | no |
| Q14721 | KCNB1 | kv_shaker | kv_modifier | kv_shaker | channel | no | yes |
| Q7Z2W7 | TRPM8 | trpm | trpm |  | unassigned | yes | no |
| O94759 | TRPM2 | trpm | trpm |  | superfamily_only | yes | no |
| P02708 | CHRNA1 | nachr | nachr |  | superfamily_only | yes | no |
| P42261 | GRIA1 | ampa | ampa |  | superfamily_only | yes | no |
| Q9UBL9 | P2RX2 | p2x | p2x |  | superfamily_only | yes | no |
| P35523 | CLCN1 | clc_channel | clc_channel |  | superfamily_only | yes | no |
| Q5XXA6 | ANO1 | ano_channel | ano_channel |  | superfamily_only | yes | no |
| Q14643 | ITPR1 | itpr | itpr |  | superfamily_only | yes | no |
| Q13255 | GRM1 | nonchannel_class_c_gpcr | nonchannel_class_c_gpcr |  | not_enumerated | yes | no |
| Q9UBS5 | GABBR1 | nonchannel_class_c_gpcr | nonchannel_class_c_gpcr |  | not_enumerated | yes | no |
| Q9Y6A1 | POMT1 | nonchannel_pomt | nonchannel_pomt |  | not_enumerated | yes | no |
| P55087 | AQP4 | outofscope_aquaporin | outofscope_aquaporin |  | not_enumerated | yes | no |

(Only rows where S1 and S2 disagree on correctness are shown; all rows are in `s1_panel.tsv`.)

## Bulk files (data root)

| file | bytes | records | sha256 |
|---|---|---|---|
| raw_api/s2/census_v2.tsv.gz | 29273333 | 1275990 | 1d9b3a46c2cc2f2514f34ccd6e72157cb53b1925dda9f524f20a6edbca41b278 |
| raw_api/s2/census_v2.fasta.gz | 291875954 | 1275990 | 88e1d1fe9d8c4a8f767d5a9fa2650e8bd13fd2d9b79d351e0e908fe7a781b692 |


## Revision r2 — hazard rules H2, H4, H13 as positive tests (S2b)

S3a's calibration and its check against the IP3R project's census found three hazard rules that made a family call from what a protein *lacks* — which this project's conventions forbid. Each was rewritten so a family is called only on a domain it carries and its rivals do not; where no such domain exists (ITPR against RyR; AChBP against a receptor fragment) the architecture tier stops at the superfamily and a sequence-level tier makes the call (D33). **108,407 records** carry an accession the rewritten rules consult (`PF02931`, `PF08016`, `PF08709`, `PF20519`) and were re-classified; **27,028 calls changed**, and 0 changed outside that set (checked, not assumed). The previous call files are archived on the data root under `calls_r1/`.

| before | after | hazards_after | records |
|---|---|---|---|
| nonchannel_achbp | [cysloop] | H2 | 12567 |
| itpr | [ca_release] | H4 | 7208 |
| trpp | assoc_polycystin1 | H13 | 3125 |
| trpp | [ploop] | H13 | 3093 |
| [ploop] | assoc_polycystin1 | H13 | 1011 |
| unassigned | [ploop] | H13 | 13 |
| nav | [ploop] | H13 | 3 |
| nonchannel_achbp | unassigned | H2,H3 | 3 |
| unassigned | nonchannel_kctd | H12,H15,H2 | 1 |
| unassigned | trpml | H13,H2,H7 | 1 |
| unassigned | kir | H16,H2,H8 | 1 |
| unassigned | k2p | H16,H2,H8 | 1 |


## Revision r3 — hazard rules H11, H12 as positive tests (S2c)

The two remaining absence rules. KCTD is now called on a KCTD C-terminal domain, and the T1 domain alone is superfamily-only; CFTR is still called on its R domain, and the SUR rule is gone — no domain identifies SUR, and none of the 980 calls the absence rule made was one — 918 bacterial ABC transporters and 62 eukaryotic fused gene models. **63,411 records** carry an accession the rewritten rules consult (`PF00664`, `PF02214`) and were re-classified; **18,412 calls changed**, and 0 changed outside that set (checked, not assumed). The previous call files are archived on the data root under `calls_r2/`.

| before | after | hazards_after | records |
|---|---|---|---|
| nonchannel_kctd | [ploop] | H12 | 17423 |
| assoc_sur | unassigned | - | 944 |
| assoc_sur | [ploop] | H1 | 27 |
| unassigned | kir | H16,H8 | 4 |
| assoc_sur | [deg_enac] | - | 2 |
| assoc_sur | [ploop] | - | 2 |
| [ploop] | kir | H12,H16,H8 | 1 |
| unassigned | nonchannel_kctd | H12,H15 | 1 |
| assoc_sur | cav | H1 | 1 |
| assoc_sur | trpc | H7 | 1 |
| assoc_sur | kca_slo | H8 | 1 |
| nonchannel_kctd | unassigned | H12,H2 | 1 |


## Revision r4 — eight signatures added after S20 (D43)

S20 found eight proposed channels no catalogue family named; they were added (`src/catalogue/proposed.py`) and their eight Pfam signatures brought into the search space. r4 walked only the delta — records with a new signature and none of r3's 67 — in the same taxonomic shards: 26,783 records, r3's 1,245,200 + delta = UniProt's 1,271,983 for the 75-signature union, every shard and signature exact (release 2026_03). The 11 r3 records that carry a new signature were re-classified; every other call is copied. TMCO1, TMEM87A and TMEM109 carry only domains shared with a non-channel (H17–H19), so domain rules leave them unassigned by design; MITOK (CCDC51) carries no Pfam domain and is not enumerable at all. **11 records** carry an accession the rewritten rules consult (`PF01956`, `PF05934`, `PF06814`, `PF12430`, `PF12537`, `PF14965`, `PF15122`, `PF21901`) and were re-classified; **2 calls changed**, and 0 changed outside that set (checked, not assumed). The previous call files are archived on the data root under `calls_r3/`.

| before | after | hazards_after | records |
|---|---|---|---|
| connexin | gphr | - | 1 |
| unassigned | gphr | - | 1 |


## Revision r5 — the viroporin signatures (S4's row, D43)

The viroporin family's four signatures (Flu M2, Vpu, CoV E, bCoV viroporin) are SUBFAMILY-level — each proves one virus's channel — and so were never enumerated (D34 couples nothing any more, but the declaration had to change). r5 walked the delta the same way as r4: 4,007 records (3,994 viral), r4's 1,271,983 + delta = UniProt's 1,275,990 for the 79-signature union, every shard and signature exact. No earlier record carries a viroporin signature, so no call was re-classified. **0 records** carry an accession the rewritten rules consult (`PF00558`, `PF00599`, `PF02723`, `PF11289`) and were re-classified; **0 calls changed**, and 0 changed outside that set (checked, not assumed). The previous call files are archived on the data root under `calls_r4/`.

| before | after | hazards_after | records |
|---|---|---|---|
