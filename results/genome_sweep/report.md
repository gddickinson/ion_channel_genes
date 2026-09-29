# S5a — the genomic sweep's instrument, and a seven-genome pilot

Rendered by `scripts/s5_report.py` from the tables in this directory (D13). Bulk output (genomes, miniprot GFFs, per-genome loci, tblastn) is under `<data root>/genomes/`.

## What was measured

- **Bait panel**: 1808 baits for 91 catalogue families (controls included) from 50 panel species — rules B1–B3 in `s5_baits.py`, one bait per species (D37).
- **Pilot**: 7 genomes, 6.02 Gbp; 1881 loci, 1020 called to a family by the S3a profiles (D32).
- **Positive control, the whole instrument**: of 175 control cells (families the genome's own proteome carries at high confidence, or a genome-only species' group-core families), 166 are *found* (a profile-called locus) with the genome's own species' baits excluded.
- **Matched detection** — the number an absence inherits: **172 / 173** control cells with a bait from another species of the same group are detected (found, partial, gap or trace); **2 / 2** without one. Detection is a property of the nearest bait.

## Per-genome control

| species | group | control | found | detected | matched | m. detected | m. rate | m. undetected | unmatched | u. detected |
|---|---|---|---|---|---|---|---|---|---|---|
| Arabidopsis thaliana | plant | 9 | 8 | 9 | 8 | 8 | 1.0 |  | 1 | 1 |
| Caenorhabditis elegans | invertebrate | 39 | 36 | 38 | 38 | 37 | 0.9737 | tweety | 1 | 1 |
| Cornu aspersum | invertebrate | 27 | 27 | 27 | 27 | 27 | 1.0 |  | 0 | 0 |
| Drosophila melanogaster | invertebrate | 40 | 38 | 40 | 40 | 40 | 1.0 |  | 0 | 0 |
| Escherichia coli | prokaryote | 3 | 2 | 3 | 3 | 3 | 1.0 |  | 0 | 0 |
| Mus musculus | vertebrate | 55 | 55 | 55 | 55 | 55 | 1.0 |  | 0 | 0 |
| Saccharomyces cerevisiae | fungi | 2 | 0 | 2 | 2 | 2 | 1.0 |  | 0 | 0 |

The control floor is 90% matched detection; a genome below it has every absence read as `uncontrolled`.

## Zero cells: what the genomes say about S3b's absences

312 zero cells in the pilot genomes; 71 informative (the family is present at high confidence in another species of the group).

| verdict | all | informative |
|---|---|---|
| no_locus_unrescued | 184 | 0 |
| partial | 50 | 5 |
| genome_present | 48 | 46 |
| absent | 10 | 10 |
| gap | 10 | 1 |
| trace | 5 | 5 |
| absent_bar_unmeasured | 3 | 3 |
| genome_weak | 2 | 1 |

### Zero cells with genomic evidence

| species | family | informative | proteome_band | genome | n_found | best_confidence | n_strong | n_partial | n_traces | trace_best_bits | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | kv_modifier | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Escherichia coli | kcsa_prok | 1 | 0 | no_locus | 0 |  | 0 | 0 | 1 | 54.3 | trace |
| Saccharomyces cerevisiae | kv_modifier | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | kca_sk | 0 | 0 | partial | 0 |  | 0 | 7 | 0 |  | partial |
| Saccharomyces cerevisiae | catsper | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Saccharomyces cerevisiae | tpc | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | hcn | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Saccharomyces cerevisiae | trpc | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | trpm | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Saccharomyces cerevisiae | trpp | 0 | 0 | gap | 0 |  | 0 | 1 | 0 |  | gap |
| Saccharomyces cerevisiae | ano_channel | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | ano_scramblase | 1 | 0 | no_locus | 0 |  | 0 | 0 | 1 | 107.0 | trace |
| Saccharomyces cerevisiae | tmc | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | piezo | 0 | 0 | partial | 0 |  | 0 | 3 | 0 |  | partial |
| Saccharomyces cerevisiae | mscs | 1 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Saccharomyces cerevisiae | clic | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Arabidopsis thaliana | kca_sk | 0 | 0 | partial | 0 |  | 0 | 5 | 0 |  | partial |
| Arabidopsis thaliana | kir | 1 | 0 | partial | 0 |  | 0 | 5 | 0 |  | partial |
| Arabidopsis thaliana | nalcn | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Arabidopsis thaliana | hcn | 0 | 0 | partial | 0 |  | 0 | 4 | 0 |  | partial |
| Arabidopsis thaliana | trpc | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Arabidopsis thaliana | trpa | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Arabidopsis thaliana | trpp | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Arabidopsis thaliana | nmda | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Arabidopsis thaliana | tweety | 1 | 0 | no_locus | 0 |  | 0 | 0 | 1 | 52.8 | trace |
| Arabidopsis thaliana | tmc | 0 | 0 | partial | 0 |  | 0 | 3 | 0 |  | partial |
| Arabidopsis thaliana | pannexin | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Arabidopsis thaliana | orai | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Caenorhabditis elegans | kv_modifier | 0 | 6 | partial | 0 |  | 0 | 4 | 0 |  | partial |
| Caenorhabditis elegans | hcn | 1 | 0 | partial | 0 |  | 0 | 4 | 0 |  | partial |
| Caenorhabditis elegans | plgic_prok | 1 | 0 | no_locus | 0 |  | 0 | 0 | 1 | 52.0 | trace |
| Caenorhabditis elegans | ampa | 0 | 3 | found | 1 | medium | 0 | 1 | 0 |  | genome_weak |
| Caenorhabditis elegans | ano_channel | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Caenorhabditis elegans | mscs | 0 | 0 | partial | 0 |  | 0 | 3 | 0 |  | partial |
| Caenorhabditis elegans | lrrc8 | 0 | 2 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Caenorhabditis elegans | tmem175 | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Caenorhabditis elegans | viroporin | 0 | 0 | partial | 0 |  | 0 | 7 | 0 |  | partial |
| Drosophila melanogaster | kv_modifier | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Drosophila melanogaster | catsper | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Drosophila melanogaster | glyr | 0 | 5 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Drosophila melanogaster | plgic_prok | 1 | 0 | no_locus | 0 |  | 0 | 0 | 3 | 72.8 | trace |
| Drosophila melanogaster | mscl | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Drosophila melanogaster | mscs | 0 | 0 | partial | 0 |  | 0 | 5 | 0 |  | partial |
| Drosophila melanogaster | connexin | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Drosophila melanogaster | calhm | 1 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Drosophila melanogaster | viroporin | 0 | 0 | partial | 0 |  | 0 | 4 | 0 |  | partial |
| Mus musculus | zac | 1 | 0 | found | 1 | medium | 0 | 0 | 0 |  | genome_weak |
| Mus musculus | p2x_nonmetazoan | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Mus musculus | deg_invertebrate | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Mus musculus | mscl | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |
| Mus musculus | mscs | 0 | 0 | gap | 0 |  | 0 | 15 | 0 |  | gap |
| Mus musculus | innexin | 0 | 0 | partial | 0 |  | 0 | 1 | 0 |  | partial |
| Mus musculus | viroporin | 0 | 0 | partial | 0 |  | 0 | 2 | 0 |  | partial |

Genome-only species (no proteome to compare): 48 families present by a profile-called locus — Cornu ampa, Cornu ano_scramblase, Cornu asic, Cornu bestrophin, Cornu cav, Cornu clc_channel, Cornu clic, Cornu cng, Cornu deg_invertebrate, Cornu delta_glur, Cornu gabaa, Cornu hcn, Cornu hv1, Cornu iglur_nonvertebrate, Cornu innexin, Cornu itpr, Cornu k2p, Cornu kainate, Cornu kca_sk, Cornu kca_slo, Cornu kir, Cornu kv_eag, Cornu kv_kcnq, Cornu kv_shaker, Cornu mcu, Cornu nachr, Cornu nalcn, Cornu nav, Cornu nmda, Cornu orai, Cornu osca_tmem63, Cornu otop, Cornu p2x, Cornu piezo, Cornu plgic_prok, Cornu ryr, Cornu tmc, Cornu tpc, Cornu tric, Cornu trpa, Cornu trpc, Cornu trpm, Cornu trpml, Cornu trpn, Cornu trpp, Cornu trpv, Cornu tweety, Cornu vdac.

### Zero cells read as absent (or blocked)

| species | family | matched | verdict |
|---|---|---|---|
| Escherichia coli | kir | 1 | absent |
| Escherichia coli | plgic_prok | 1 | absent_bar_unmeasured |
| Escherichia coli | iglur_prok | 1 | absent_bar_unmeasured |
| Arabidopsis thaliana | kca_slo | 1 | absent |
| Caenorhabditis elegans | nav | 1 | absent |
| Caenorhabditis elegans | tpc | 1 | absent |
| Caenorhabditis elegans | p2x | 1 | absent |
| Caenorhabditis elegans | osca_tmem63 | 1 | absent |
| Caenorhabditis elegans | hv1 | 1 | absent |
| Drosophila melanogaster | tpc | 1 | absent |
| Drosophila melanogaster | p2x | 1 | absent |
| Drosophila melanogaster | deg_invertebrate | 1 | absent_bar_unmeasured |
| Drosophila melanogaster | hv1 | 1 | absent |

### Literature-expected absences (external check, references pending)

| species | family | why | proteome | genome | traces | verdict |
|---|---|---|---|---|---|---|
| Mus musculus | zac | ZAC has no rodent orthologue | 0 | found | 0 | genome_weak |
| Caenorhabditis elegans | nav | nematodes lack voltage-gated Na+ channels | 0 | no_locus | 0 | absent |
| Caenorhabditis elegans | p2x | no P2X receptor in C. elegans | 0 | no_locus | 0 | absent |
| Drosophila melanogaster | p2x | no P2X receptor in Drosophila | 0 | no_locus | 0 | absent |

## Calibration from annotation (not from miniprot)

Widest annotated intron at a high-confidence family-called locus, against the `-G` used:

| species | annotation_accession | called_loci | annotated_loci | widest_annotated_intron | widest_gene | max_intron_used | over_G |
|---|---|---|---|---|---|---|---|
| Escherichia coli | GCA_000005845.2 | 2 | 2 | 0 | mscL (mscl) | 2000 | 0 |
| Saccharomyces cerevisiae | GCA_000146045.2 | 0 | 0 |  |  | 200000 | 0 |
| Arabidopsis thaliana | GCA_000001735.2 | 33 | 32 | 766 | KCO2 (k2p) | 200000 | 0 |
| Caenorhabditis elegans | GCA_000002985.3 | 101 | 98 | 20089 | unc-7 (innexin) | 200000 | 0 |
| Drosophila melanogaster | GCA_000001215.4 | 68 | 68 | 87537 | shakB (innexin) | 200000 | 0 |
| Cornu aspersum |  | 152 | 0 |  |  | 1000000 |  |
| Mus musculus | GCF_000001635.27 | 231 | 220 | 996015 | Asic2 (asic) | 1000000 | 0 |

Gene spans (420 annotated loci) → D4's bar per family (median span; 61 families measured):

| family | n_genes | n_species | median_span | min_span | max_span | by_group |
|---|---|---|---|---|---|---|
| ampa | 4 | 1 | 298428 | 122778 | 378358 | vertebrate:298428 |
| ano_channel | 2 | 1 | 256605 | 163477 | 349734 | vertebrate:256605 |
| ano_scramblase | 9 | 4 | 22140 | 4165 | 396266 | invertebrate:14638;plant:4165;vertebrate:107569 |
| asic | 5 | 2 | 24518 | 4388 | 1088234 | invertebrate:4470;vertebrate:27720 |
| bestrophin | 21 | 3 | 2725 | 476 | 43470 | invertebrate:2700;vertebrate:29966 |
| calhm | 3 | 1 | 5168 | 3140 | 6063 | vertebrate:5168 |
| catsper | 2 | 1 | 19221 | 15468 | 22974 | vertebrate:19221 |
| cav | 16 | 3 | 66365 | 13159 | 610145 | invertebrate:25051;vertebrate:135197 |
| cftr | 1 | 1 | 152084 | 152084 | 152084 | vertebrate:152084 |
| clc_channel | 6 | 3 | 12959 | 5707 | 29142 | invertebrate:9955;vertebrate:29142 |
| clic | 8 | 3 | 8907 | 1222 | 145389 | invertebrate:9337;plant:1395;vertebrate:50943 |
| cng | 8 | 3 | 28713 | 4343 | 57291 | invertebrate:35938;vertebrate:28713 |
| connexin | 10 | 1 | 12637 | 2620 | 23571 | vertebrate:12637 |
| delta_glur | 2 | 1 | 1107391 | 763594 | 1451188 | vertebrate:1107391 |
| enac | 3 | 1 | 35250 | 34063 | 53691 | vertebrate:35250 |
| gabaa | 13 | 2 | 88575 | 13092 | 491652 | invertebrate:16435;vertebrate:89621 |
| glyr | 4 | 1 | 141809 | 25829 | 198851 | vertebrate:141809 |
| hcn | 5 | 2 | 23966 | 13457 | 378709 | invertebrate:23966;vertebrate:29851 |
| ht3 | 1 | 1 | 11939 | 11939 | 11939 | vertebrate:11939 |
| hv1 | 1 | 1 | 35489 | 35489 | 35489 | vertebrate:35489 |
| iglur_nonvertebrate | 7 | 2 | 4155 | 3647 | 4455 | invertebrate:4455;plant:4078 |
| innexin | 12 | 2 | 5658 | 2088 | 165980 | invertebrate:5658 |
| itpr | 5 | 3 | 65222 | 21757 | 393961 | invertebrate:22026;vertebrate:338034 |
| k2p | 25 | 4 | 7497 | 1765 | 194825 | invertebrate:6063;plant:2218;vertebrate:38733 |
| kainate | 5 | 1 | 394963 | 66148 | 696748 | vertebrate:394963 |
| kca_sk | 5 | 3 | 64903 | 22534 | 417201 | invertebrate:49380;vertebrate:152331 |
| kca_slo | 8 | 3 | 71897 | 10340 | 712818 | invertebrate:32528;vertebrate:228863 |
| kir | 14 | 3 | 17203 | 5829 | 252881 | invertebrate:11528;vertebrate:29639 |
| kv_eag | 9 | 3 | 52757 | 5625 | 490553 | invertebrate:31358;vertebrate:280635 |
| kv_kcnq | 8 | 3 | 56238 | 8984 | 564252 | invertebrate:9023;vertebrate:300263 |
| kv_modifier | 5 | 1 | 21417 | 6053 | 69739 | vertebrate:21417 |
| kv_shaker | 14 | 3 | 15148 | 3005 | 517704 | invertebrate:13886;vertebrate:18349 |
| lrrc8 | 3 | 1 | 26022 | 13202 | 94035 | vertebrate:26022 |
| mcu | 4 | 3 | 49222 | 1784 | 170426 | invertebrate:22493;vertebrate:112834 |
| mscl | 1 | 1 | 411 | 411 | 411 | prokaryote:411 |
| mscs | 8 | 2 | 3519 | 861 | 3898 | plant:3591;prokaryote:861 |
| nachr | 29 | 3 | 6760 | 2259 | 113835 | invertebrate:5239;vertebrate:18385 |
| nalcn | 4 | 3 | 12609 | 2835 | 350801 | invertebrate:11050;vertebrate:350801 |
| nav | 7 | 2 | 95627 | 30064 | 177220 | invertebrate:54067;vertebrate:146684 |
| nmda | 7 | 3 | 18129 | 5320 | 460440 | invertebrate:10428;vertebrate:427980 |
| orai | 4 | 3 | 12395 | 5336 | 25109 | invertebrate:12395;vertebrate:15222 |
| osca_tmem63 | 12 | 3 | 4222 | 3536 | 70485 | invertebrate:3623;plant:4176;vertebrate:32773 |
| otop | 7 | 3 | 4735 | 3176 | 26813 | invertebrate:3989;vertebrate:25978 |
| p2x | 3 | 1 | 22222 | 3509 | 40554 | vertebrate:22222 |
| pannexin | 3 | 1 | 15973 | 9321 | 39694 | vertebrate:15973 |
| piezo | 4 | 3 | 48348 | 18613 | 377504 | invertebrate:22839;vertebrate:223568 |
| plgic_invertebrate | 5 | 2 | 6080 | 2284 | 44517 | invertebrate:6080 |
| ryr | 5 | 3 | 121835 | 15515 | 586057 | invertebrate:21624;vertebrate:553847 |
| tmc | 7 | 3 | 41949 | 6717 | 170747 | invertebrate:14014;vertebrate:59120 |
| tmem175 | 1 | 1 | 17988 | 17988 | 17988 | vertebrate:17988 |
| tpc | 3 | 2 | 54510 | 6052 | 102374 | plant:6052;vertebrate:78442 |
| tric | 4 | 3 | 12745 | 3425 | 36032 | invertebrate:6409;vertebrate:26064 |
| trpa | 3 | 3 | 10797 | 2482 | 46219 | invertebrate:6639;vertebrate:46219 |
| trpc | 11 | 3 | 45314 | 5769 | 310854 | invertebrate:10051;vertebrate:129877 |
| trpm | 10 | 3 | 73322 | 5113 | 857986 | invertebrate:14139;vertebrate:102376 |
| trpml | 4 | 3 | 11088 | 3911 | 48532 | invertebrate:5656;vertebrate:31653 |
| trpn | 2 | 2 | 17242 | 14924 | 19561 | invertebrate:17242 |
| trpp | 4 | 2 | 44298 | 3750 | 99280 | invertebrate:3750;vertebrate:54307 |
| trpv | 12 | 3 | 6321 | 3185 | 36285 | invertebrate:4225;vertebrate:27174 |
| tweety | 3 | 1 | 28495 | 16789 | 45517 | vertebrate:28495 |
| vdac | 9 | 4 | 3086 | 1251 | 28320 | invertebrate:2278;plant:2580;vertebrate:16739 |

## Runs

| species | assembly | total_bp | n50 | max_intron | chunks | alignments | self_dropped | loci | family_called_loci | miniprot_s | call_s | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | GCA_000005845.2 | 4641652 | 4641652 | 2000 | 1 | 62 | 7 | 13 | 8 | 1.7 | 0.7 |  |
| Saccharomyces cerevisiae | GCA_000146045.2 | 12071326 | 924431 | 200000 | 1 | 338 | 8 | 54 | 19 | 20.8 | 4.7 |  |
| Arabidopsis thaliana | GCA_000001735.2 | 119668634 | 23459830 | 200000 | 1 | 985 | 27 | 166 | 105 | 45.7 | 9.5 |  |
| Caenorhabditis elegans | GCA_000002985.3 | 100272607 | 17493829 | 200000 | 1 | 1705 | 47 | 196 | 156 | 44.6 | 18.5 |  |
| Drosophila melanogaster | GCA_000001215.4 | 143726002 | 25286936 | 200000 | 1 | 1930 | 51 | 182 | 124 | 48.3 | 19.8 |  |
| Cornu aspersum | GCA_964187895.1 | 2908462832 | 110296109 | 1000000 | 2 | 3828 | 0 | 692 | 276 | 746.9 | 60.9 |  |
| Mus musculus | GCA_000001635.9 | 2728222451 | 130530862 | 1000000 | 2 | 3402 | 117 | 578 | 332 | 614.1 | 42.2 |  |
