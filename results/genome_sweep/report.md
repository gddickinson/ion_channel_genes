# S5 — the genomic sweep: what 52 genomes say about the proteome census

Rendered by `scripts/s5_report.py` from the tables in this directory (D13). S5a built and measured the instrument on seven genomes; S5b ran it over the whole S4 panel. Bulk output (genomes, miniprot GFFs, per-genome loci, tblastn, census v4) is under `<data root>/genomes/`.

## What was measured

- **Bait panel**: 1808 baits for 91 catalogue families (controls included) from 50 panel species, one per species (B1–B3, D37).
- **Sweep**: 52 genomes, 42.5 Gbp, 0 failures; 11,878 loci, **8,048 called to a family by the S3a profiles** (D32); miniprot 1.94 h. `-G` 1 Mb for genomes ≥ 1 Gbp, measured not assumed (D38).
- **Positive control**: 1,432 of 1,558 control cells *found* (a profile-called locus), 1,531 detected by the whole instrument, with the genome's own species' baits excluded.
- **Matched detection — the number every absence inherits: 1,441 / 1,451** control cells with an in-group non-self bait; unmatched 90 / 107. Below the 90% floor: Human immunodeficiency virus type 1, Influenza A virus, Severe acute respiratory syndrome coronavirus 2. No matched control at all (single-species groups — no absence readable): Chlamydomonas reinhardtii, Dictyostelium discoideum, Plasmodium falciparum, Trypanosoma brucei.
- **Census v4**: census v3's 28,891 proteome rows unchanged + **449 genome loci** (432 from the genome-only species, 17 proteome misses) in 128 cells (`s5_census_v4.py`; SHA-256 in `census_v4.json`). A genome row is a locus, never merged into a proteome call.

## Per-genome control

| species | group | control | found | detected | matched | m. detected | m. rate | m. undetected | unmatched | u. detected |
|---|---|---|---|---|---|---|---|---|---|---|
| Aliarcobacter butzleri | prokaryote | 2 | 1 | 2 | 2 | 2 | 1.0 |  | 0 | 0 |
| Amphimedon queenslandica | basal_metazoan | 29 | 25 | 26 | 17 | 17 | 1.0 |  | 12 | 9 |
| Anolis carolinensis | vertebrate | 62 | 62 | 62 | 62 | 62 | 1.0 |  | 0 | 0 |
| Aplysia californica | invertebrate | 46 | 44 | 46 | 46 | 46 | 1.0 |  | 0 | 0 |
| Arabidopsis thaliana | plant | 11 | 10 | 11 | 10 | 10 | 1.0 |  | 1 | 1 |
| Bacillus subtilis | prokaryote | 2 | 1 | 2 | 2 | 2 | 1.0 |  | 0 | 0 |
| Branchiostoma floridae | deuterostome | 51 | 45 | 49 | 46 | 45 | 0.9783 | tweety | 5 | 4 |
| Caenorhabditis elegans | invertebrate | 41 | 38 | 40 | 40 | 39 | 0.975 | tweety | 1 | 1 |
| Callorhinchus milii | vertebrate | 61 | 61 | 61 | 61 | 61 | 1.0 |  | 0 | 0 |
| Capsaspora owczarzaki | holozoa | 19 | 9 | 17 | 10 | 10 | 1.0 |  | 9 | 7 |
| Chlamydomonas reinhardtii | algae | 13 | 3 | 9 | 0 | 0 |  |  | 13 | 9 |
| Ciona intestinalis | deuterostome | 37 | 32 | 35 | 35 | 33 | 0.9429 | tmem175,tric | 2 | 2 |
| Cornu aspersum | invertebrate | 28 | 28 | 28 | 28 | 28 | 1.0 |  | 0 | 0 |
| Danio rerio | vertebrate | 60 | 60 | 60 | 60 | 60 | 1.0 |  | 0 | 0 |
| Daphnia pulex | invertebrate | 39 | 37 | 39 | 39 | 39 | 1.0 |  | 0 | 0 |
| Dictyostelium discoideum | amoebozoa | 12 | 6 | 10 | 0 | 0 |  |  | 12 | 10 |
| Drosophila melanogaster | invertebrate | 42 | 40 | 42 | 42 | 42 | 1.0 |  | 0 | 0 |
| Escherichia coli | prokaryote | 3 | 2 | 3 | 3 | 3 | 1.0 |  | 0 | 0 |
| Gallus gallus | vertebrate | 61 | 60 | 61 | 61 | 61 | 1.0 |  | 0 | 0 |
| Gloeobacter violaceus | prokaryote | 4 | 1 | 3 | 3 | 3 | 1.0 |  | 1 | 0 |
| Homo sapiens | vertebrate | 63 | 63 | 63 | 63 | 63 | 1.0 |  | 0 | 0 |
| Human immunodeficiency virus type 1 | virus | 1 | 0 | 0 | 1 | 0 | 0.0 | viroporin | 0 | 0 |
| Hydra vulgaris | cnidarian | 40 | 34 | 40 | 29 | 29 | 1.0 |  | 11 | 11 |
| Influenza A virus | virus | 1 | 0 | 0 | 1 | 0 | 0.0 | viroporin | 0 | 0 |
| Latimeria chalumnae | vertebrate | 62 | 61 | 62 | 62 | 62 | 1.0 |  | 0 | 0 |
| Lottia gigantea | invertebrate | 42 | 41 | 41 | 41 | 40 | 0.9756 | clcc1 | 1 | 1 |
| Lymnaea stagnalis | invertebrate | 43 | 42 | 43 | 42 | 42 | 1.0 |  | 1 | 1 |
| Methanothermobacter thermautotrophicus | prokaryote | 2 | 0 | 2 | 2 | 2 | 1.0 |  | 0 | 0 |
| Monodelphis domestica | vertebrate | 63 | 63 | 63 | 63 | 63 | 1.0 |  | 0 | 0 |
| Monosiga brevicollis | holozoa | 18 | 11 | 17 | 10 | 10 | 1.0 |  | 8 | 7 |
| Mus musculus | vertebrate | 62 | 62 | 62 | 62 | 62 | 1.0 |  | 0 | 0 |
| Nematostella vectensis | cnidarian | 38 | 34 | 37 | 34 | 33 | 0.9706 | tweety | 4 | 4 |
| Ornithorhynchus anatinus | vertebrate | 62 | 62 | 62 | 62 | 62 | 1.0 |  | 0 | 0 |
| Oryza sativa | plant | 10 | 10 | 10 | 10 | 10 | 1.0 |  | 0 | 0 |
| Paramecium tetraurelia | ciliate | 14 | 8 | 14 | 12 | 12 | 1.0 |  | 2 | 2 |
| Petromyzon marinus | vertebrate | 59 | 59 | 59 | 59 | 59 | 1.0 |  | 0 | 0 |
| Physcomitrium patens | plant | 12 | 10 | 12 | 9 | 9 | 1.0 |  | 3 | 3 |
| Plasmodium falciparum | apicomplexa | 5 | 0 | 3 | 0 | 0 |  |  | 5 | 3 |
| Rattus norvegicus | vertebrate | 62 | 61 | 62 | 62 | 62 | 1.0 |  | 0 | 0 |
| Saccharomyces cerevisiae | fungi | 3 | 0 | 3 | 2 | 2 | 1.0 |  | 1 | 1 |
| Schizosaccharomyces pombe | fungi | 4 | 1 | 4 | 2 | 2 | 1.0 |  | 2 | 2 |
| Severe acute respiratory syndrome coronavirus 2 | virus | 1 | 0 | 0 | 1 | 0 | 0.0 | viroporin | 0 | 0 |
| Streptomyces lividans | prokaryote | 3 | 1 | 3 | 3 | 3 | 1.0 |  | 0 | 0 |
| Strongylocentrotus purpuratus | deuterostome | 46 | 42 | 45 | 46 | 45 | 0.9783 | tmem175 | 0 | 0 |
| Synechocystis sp. PCC 6803 | prokaryote | 5 | 3 | 5 | 3 | 3 | 1.0 |  | 2 | 2 |
| Takifugu rubripes | vertebrate | 56 | 55 | 56 | 56 | 56 | 1.0 |  | 0 | 0 |
| Tetrahymena thermophila | ciliate | 12 | 4 | 12 | 12 | 12 | 1.0 |  | 0 | 0 |
| Thermus thermophilus | prokaryote | 2 | 1 | 2 | 2 | 2 | 1.0 |  | 0 | 0 |
| Torpedo marmorata | vertebrate | 53 | 53 | 53 | 53 | 53 | 1.0 |  | 0 | 0 |
| Trichoplax adhaerens | basal_metazoan | 31 | 28 | 31 | 22 | 22 | 1.0 |  | 9 | 9 |
| Trypanosoma brucei | excavate | 2 | 0 | 1 | 0 | 0 |  |  | 2 | 1 |
| Xenopus tropicalis | vertebrate | 58 | 58 | 58 | 58 | 58 | 1.0 |  | 0 | 0 |

## Zero cells: what the genomes say about S3b's absences

2,274 zero cells (census family × species with no proteome call); 279 informative (the family is present at high confidence in another species of the group). Every cell is in `cells.tsv`.

| verdict | all | informative |
|---|---|---|
| no_locus_unrescued | 1583 | 0 |
| partial | 328 | 29 |
| genome_present | 114 | 112 |
| gap | 111 | 11 |
| absent | 83 | 83 |
| genome_weak | 19 | 10 |
| trace | 18 | 18 |
| genome_found | 14 | 12 |
| unmatched | 4 | 4 |

Informative zero cells by group:

| group | genome_present | absent | partial | trace | genome_found | gap | genome_weak | unmatched | total |
|---|---|---|---|---|---|---|---|---|---|
| prokaryote |  | 29 |  | 3 |  |  |  |  | 32 |
| fungi |  | 1 | 1 | 1 |  |  |  |  | 3 |
| plant |  | 1 | 6 | 1 |  |  |  |  | 8 |
| ciliate |  |  |  | 1 |  |  |  |  | 1 |
| holozoa |  | 2 | 2 | 1 | 1 | 4 | 2 | 2 | 14 |
| basal_metazoan |  | 6 | 2 | 1 |  | 3 |  | 2 | 14 |
| cnidarian |  | 1 | 1 | 2 | 1 |  | 1 |  | 6 |
| invertebrate | 50 | 22 | 4 | 5 | 5 | 1 |  |  | 87 |
| deuterostome |  | 6 | 3 | 3 | 1 |  | 4 |  | 17 |
| vertebrate | 62 | 15 | 10 |  | 4 | 3 | 3 |  | 97 |

### Proteome misses: a high-confidence intact gene the proteome has no call for

`genome_found` (D37 (4)). The check beside each is reported, not applied (`genome_found_check.tsv`, added after the cells were read): the genome's loci of the family, and proteome entries whose runner-up profile is that family — a gene the proteome might hold under a sister family's call.

| species | family | genome_loci | genome_loci_called | proteome_near | near_min_margin |
|---|---|---|---|---|---|
| Monosiga brevicollis | p2x_nonmetazoan | 1 | 1 | 1 | 0.7876 |
| Trichoplax adhaerens | trpn | 1 | 1 | 1 | 0.7958 |
| Nematostella vectensis | trpn | 2 | 2 | 6 | 0.4029 |
| Nematostella vectensis | mscs | 6 | 2 | 0 |  |
| Daphnia pulex | hcn | 1 | 1 | 1 | 0.7983 |
| Daphnia pulex | trpn | 1 | 1 | 5 | 0.4192 |
| Daphnia pulex | gphr | 1 | 1 | 0 |  |
| Lottia gigantea | trpa | 2 | 2 | 0 |  |
| Lottia gigantea | trpn | 1 | 1 | 8 | 0.2036 |
| Ciona intestinalis | nalcn | 1 | 1 | 0 |  |
| Takifugu rubripes | ryr | 6 | 6 | 3 | 0.9355 |
| Takifugu rubripes | tmem175 | 1 | 1 | 0 |  |
| Xenopus tropicalis | tmem87 | 2 | 1 | 1 | 0.2736 |
| Xenopus tropicalis | mitok | 1 | 1 | 0 |  |

### Weak loci (not claims)

A profile-called locus without a high-confidence intact frame: retrocopies, frameshifted or chained models. Never read as a proteome miss.

| species | family | informative | n_found | best_confidence | n_strong | found_loci |
|---|---|---|---|---|---|---|
| Schizosaccharomyces pombe | hcn | 0 | 1 | medium | 0 | CU329670.1:3580313-3648173+ |
| Oryza sativa | orai | 0 | 1 | medium | 0 | AP014964.1:14295917-14353313+ |
| Chlamydomonas reinhardtii | tmem109 | 0 | 1 | medium | 0 | r4_CM008964.1:2476484-2543755- |
| Paramecium tetraurelia | kcsa_prok | 0 | 1 | medium | 0 | CT868665.1:186255-266224- |
| Monosiga brevicollis | hcn | 0 | 1 | medium | 0 | CH991580.1:26721-72525- |
| Monosiga brevicollis | trpa | 0 | 2 | medium | 0 | CH991544.1:112966-215252+;CH991570.1:73607-167889+ |
| Monosiga brevicollis | piezo | 1 | 2 | medium | 0 | CH991545.1:322709-350950-;CH991545.1:1054969-1074783- |
| Monosiga brevicollis | orai | 1 | 1 | high | 0 | CH991543.1:36116-88299+ |
| Nematostella vectensis | piezo | 1 | 1 | high | 0 | DS469520.1:1126602-1171342- |
| Caenorhabditis elegans | ampa | 0 | 1 | medium | 0 | BX284603.4:8584641-8588089- |
| Aplysia californica | glyr | 0 | 1 | medium | 0 | NW_004797315.1:100710-117076- |
| Strongylocentrotus purpuratus | tmem87 | 0 | 1 | medium | 0 | r4_AAGJ06000005.1:43944861-43968352+ |
| Ciona intestinalis | trpn | 1 | 1 | high | 0 | HT000012.2:4294893-4318182- |
| Ciona intestinalis | nmda | 1 | 1 | medium | 0 | HT000014.2:4000191-4004598- |
| Ciona intestinalis | piezo | 1 | 1 | high | 0 | HT000012.2:5256279-5295291- |
| Ciona intestinalis | gphr | 1 | 1 | high | 0 | r4_HT000005.2:1226492-1227891+ |
| Petromyzon marinus | tmem87 | 1 | 1 | medium | 0 | r4_NC_133697.1:6720396-7201240+ |
| Xenopus tropicalis | tmem175 | 1 | 1 | medium | 0 | NC_030684.2:112668717-113117206- |
| Mus musculus | zac | 1 | 1 | medium | 0 | CM000997.3:41266303-41666273+ |

### Absences (and blocked absences), by family

`absent` = no locus, no trace, a matched bait, a controlled genome and N50 ≥ the D4 bar; the bar is `cds` (prokaryote/virus: 3 × band), `group` (the group's own median span, ≥ 3 genes annotated) or `pooled` (D38).

| family | verdict | n | species | bar |
|---|---|---|---|---|
| bestrophin | absent | 6 | B. subtilis, S. lividans, G. violaceus, A. butzleri, M. thermautotrophicus, T. thermophilus | cds |
| calhm | absent | 5 | H. vulgaris, D. pulex, A. californica, L. stagnalis, S. purpuratus | pooled |
| catsper | absent | 1 | P. marinus | group |
| clcc1 | absent | 3 | C. elegans, D. melanogaster, D. pulex | pooled |
| deg_invertebrate | absent | 2 | D. melanogaster, D. pulex | group |
| enac | absent | 2 | D. rerio, T. rubripes | group |
| gphr | absent | 1 | S. pombe | pooled |
| hv1 | absent | 2 | C. elegans, D. melanogaster | group |
| iglur_prok | absent | 6 | E. coli, S. lividans, G. violaceus, A. butzleri, M. thermautotrophicus, T. thermophilus | cds |
| kca_slo | absent | 1 | A. thaliana | pooled |
| kcsa_prok | absent | 1 | T. thermophilus | cds |
| kir | absent | 7 | E. coli, B. subtilis, S. lividans, G. violaceus, A. butzleri, M. thermautotrophicus, T. thermophilus | cds |
| mcu | absent | 1 | T. adhaerens | pooled |
| mitok | absent | 2 | D. melanogaster, L. gigantea | pooled |
| mscl | absent | 2 | A. butzleri, M. thermautotrophicus | cds |
| nav | absent | 2 | A. queenslandica, C. elegans | group, pooled |
| osca_tmem63 | absent | 2 | C. elegans, C. intestinalis | group, pooled |
| p2x | absent | 3 | C. elegans, D. melanogaster, C. intestinalis | group |
| p2x_nonmetazoan | absent | 1 | T. adhaerens | pooled |
| pacc | absent | 2 | T. adhaerens, T. rubripes | group, pooled |
| pannexin | absent | 2 | S. purpuratus, C. intestinalis | pooled |
| piezo | absent | 1 | D. pulex | group |
| plgic_prok | absent | 9 | E. coli, B. subtilis, S. lividans, S. sp. PCC 6803, A. butzleri, M. thermautotrophicus, T. thermophilus, D. pulex, C. intestinalis | cds, pooled |
| tmco1 | absent | 1 | M. brevicollis | pooled |
| tmem109 | absent | 2 | P. marinus, C. milii | group |
| tmem175 | absent | 1 | T. adhaerens | pooled |
| tpc | absent | 2 | C. elegans, D. melanogaster | group |
| tric | absent | 1 | M. brevicollis | pooled |
| trpp | absent | 1 | D. pulex | group |
| tweety | absent | 2 | T. adhaerens, L. gigantea | pooled |
| zac | absent | 9 | C. milii, T. marmorata, D. rerio, T. rubripes, L. chalumnae, X. tropicalis, G. gallus, O. anatinus, R. norvegicus | pooled |

### Literature-expected absences (external check, references pending)

| species | family | why | proteome_records | genome | n_traces | verdict |
|---|---|---|---|---|---|---|
| Mus musculus | zac | ZAC has no rodent orthologue | 0 | found | 0 | genome_weak |
| Rattus norvegicus | zac | ZAC has no rodent orthologue | 0 | no_locus | 0 | absent |
| Caenorhabditis elegans | nav | nematodes lack voltage-gated Na+ channels | 0 | no_locus | 0 | absent |
| Caenorhabditis elegans | p2x | no P2X receptor in C. elegans | 0 | no_locus | 0 | absent |
| Drosophila melanogaster | p2x | no P2X receptor in Drosophila | 0 | no_locus | 0 | absent |

### Genome-only species

- *Cornu aspersum*: **52 census families present** by a profile-called locus — ampa, ano_scramblase, asic, bestrophin, cav, clc_channel, clcc1, clic, cng, deg_invertebrate, delta_glur, gabaa, gphr, hcn, hv1, iglur_nonvertebrate, innexin, itpr, k2p, kainate, kca_sk, kca_slo, kir, kv_eag, kv_kcnq, kv_shaker, mcu, mitok, nachr, nalcn, nav, nmda, orai, osca_tmem63, otop, p2x, piezo, plgic_prok, ryr, tmc, tmco1, tpc, tric, trpa, trpc, trpm, trpml, trpn, trpp, trpv, tweety, vdac.
- *Torpedo marmorata*: **62 census families present** by a profile-called locus — ampa, ano_channel, ano_scramblase, asic, bestrophin, calhm, catsper, cav, cftr, clc_channel, clcc1, clic, cng, connexin, delta_glur, enac, gabaa, glyr, gphr, hcn, ht3, hv1, itpr, k2p, kainate, kca_sk, kca_slo, kir, kv_eag, kv_kcnq, kv_modifier, kv_shaker, lrrc8, mcu, mitok, nachr, nalcn, nav, nmda, orai, osca_tmem63, otop, p2x, pacc, pannexin, piezo, ryr, tmc, tmco1, tmem175, tmem87, tpc, tric, trpa, trpc, trpm, trpml, trpn, trpp, trpv, tweety, vdac.

## Calibration from annotation (not from miniprot)

Widest annotated intron at a high-confidence family-called locus, against the `-G` used (`over_G` = annotated genes with an intron wider than it):

| species | annotation_accession | called_loci | annotated_loci | widest_annotated_intron | widest_gene | max_intron_used | over_G |
|---|---|---|---|---|---|---|---|
| Escherichia coli | GCA_000005845.2 | 2 | 2 | 0 | mscL (mscl) | 2000 | 0 |
| Bacillus subtilis | GCA_000009045.1 | 1 | 1 | 0 | mscL (mscl) | 2000 | 0 |
| Streptomyces lividans | GCA_000739105.1 | 1 | 1 | 0 | SLIV_21740 (mscl) | 2000 | 0 |
| Synechocystis sp. PCC 6803 | GCA_000009725.1 | 2 | 2 | 0 | mscL (mscl) | 2000 | 0 |
| Gloeobacter violaceus | GCA_000011385.1 | 1 | 1 | 0 | glr1073 (mscl) | 2000 | 0 |
| Aliarcobacter butzleri | GCA_000014025.1 | 1 | 1 | 0 | mscS (mscs) | 2000 | 0 |
| Methanothermobacter thermautotrophicus | GCA_000008645.1 | 0 | 0 |  |  | 2000 | 0 |
| Thermus thermophilus | GCA_000091545.1 | 1 | 1 | 0 | TTHA0627 (mscl) | 2000 | 0 |
| Saccharomyces cerevisiae | GCA_000146045.2 | 0 | 0 |  |  | 200000 | 0 |
| Schizosaccharomyces pombe | GCA_000002945.3 | 1 | 1 | 0 | SPOM_SPAC24H6.13 (osca_tmem63) | 200000 | 0 |
| Arabidopsis thaliana | GCA_000001735.2 | 37 | 36 | 798 | GTG2 (gphr) | 200000 | 0 |
| Oryza sativa | GCA_001433935.1 | 25 | 24 | 2026 | Os03g0202200 (vdac) | 200000 | 0 |
| Physcomitrium patens | GCF_000002425.5 | 22 | 22 | 1257 | LOC112275199 (clic) | 200000 | 0 |
| Chlamydomonas reinhardtii | GCA_000002595.3 | 5 | 4 | 614 | CHLRE_17g720600v5 (cav) | 200000 | 0 |
| Dictyostelium discoideum | GCA_000004695.1 | 8 | 8 | 227 | orfR1062 (osca_tmem63) | 200000 | 0 |
| Paramecium tetraurelia | GCA_000165425.1 | 12 | 11 | 37 | GSPATT00036989001 (cng) | 200000 | 0 |
| Tetrahymena thermophila | GCA_000189635.1 | 4 | 4 | 648 | TTHERM_00378480 (cav) | 200000 | 0 |
| Trypanosoma brucei | GCA_000002445.1 | 0 | 0 |  |  | 200000 | 0 |
| Plasmodium falciparum | GCA_000002765.3 | 0 | 0 |  |  | 200000 | 0 |
| Monosiga brevicollis | GCA_000002865.1 | 17 | 15 | 426 | MONBRDRAFT_33881 (kca_slo) | 200000 | 0 |
| Capsaspora owczarzaki | GCA_000151315.2 | 11 | 11 | 1886 | CAOG_009776 (itpr) | 200000 | 0 |
| Amphimedon queenslandica | GCF_000090795.2 | 32 | 31 | 4833 | LOC100636084 (gphr) | 200000 | 0 |
| Trichoplax adhaerens | GCA_000150275.1 | 61 | 45 | 8354 | TRIADDRAFT_51003 (k2p) | 200000 | 0 |
| Nematostella vectensis | GCA_000209225.1 | 123 | 88 | 12060 | NEMVEDRAFT_v1g118813 (kv_eag) | 200000 | 0 |
| Hydra vulgaris | GCF_038396675.1 | 74 | 74 | 35258 | LOC100204665 (cav) | 200000 | 0 |
| Caenorhabditis elegans | GCA_000002985.3 | 106 | 103 | 20089 | unc-7 (innexin) | 200000 | 0 |
| Drosophila melanogaster | GCA_000001215.4 | 71 | 71 | 87537 | shakB (innexin) | 200000 | 0 |
| Daphnia pulex | GCA_000187875.1 | 90 | 84 | 35180 | DAPPUDRAFT_321681 (nachr) | 200000 | 0 |
| Lottia gigantea | GCA_000327385.1 | 108 | 95 | 38063 | LOTGIDRAFT_105753 (kv_shaker) | 200000 | 0 |
| Aplysia californica | GCF_000002075.1 | 72 | 66 | 204025 | LOC101855406 (trpc) | 200000 | 1 |
| Lymnaea stagnalis | GCA_964033795.1 | 98 | 91 | 76732 | GSLYS_00019449001 (kca_sk) | 200000 | 0 |
| Cornu aspersum |  | 156 | 0 |  |  | 1000000 |  |
| Strongylocentrotus purpuratus | GCF_000002235.5 | 109 | 107 | 64435 | LOC589158 (trpc) | 200000 | 0 |
| Branchiostoma floridae | GCF_000003815.2 | 125 | 119 | 256433 | LOC118408435 (asic) | 200000 | 1 |
| Ciona intestinalis | GCF_000224145.3 | 82 | 81 | 28393 | LOC104265353 (nav) | 200000 | 0 |
| Petromyzon marinus | GCF_048934315.1 | 340 | 277 | 364747 | KCND1 (kv_shaker) | 1000000 | 0 |
| Callorhinchus milii | GCF_000165045.1 | 218 | 214 | 223342 | asic2 (asic) | 200000 | 2 |
| Torpedo marmorata |  | 236 | 0 |  |  | 1000000 |  |
| Danio rerio | GCF_049306965.2 | 340 | 338 | 728175 | asic4b (asic) | 1000000 | 0 |
| Takifugu rubripes | GCF_901000725.3 | 358 | 344 | 149854 | asic2 (asic) | 200000 | 0 |
| Latimeria chalumnae | GCF_000225785.1 | 244 | 243 | 637219 | GRID1 (delta_glur) | 1000000 | 0 |
| Xenopus tropicalis | GCF_000004195.4 | 240 | 238 | 217542 | asic2 (asic) | 1000000 | 0 |
| Anolis carolinensis | GCF_000090745.2 | 217 | 215 | 471577 | asic2 (asic) | 1000000 | 0 |
| Gallus gallus | GCF_016699485.2 | 210 | 209 | 294756 | TRPM3 (trpm) | 1000000 | 0 |
| Ornithorhynchus anatinus | GCF_004115215.2 | 230 | 229 | 410530 | KCNQ5 (kv_kcnq) | 1000000 | 0 |
| Monodelphis domestica | GCF_000002295.2 | 247 | 240 | 801144 | TRPM3 (trpm) | 1000000 | 0 |
| Mus musculus | GCF_000001635.27 | 238 | 227 | 996015 | Asic2 (asic) | 1000000 | 0 |
| Rattus norvegicus | GCF_036323735.1 | 248 | 233 | 975271 | Asic2 (asic) | 1000000 | 0 |
| Homo sapiens | GCF_000001405.40 | 314 | 277 | 1043910 | ASIC2 (asic) | 1000000 | 1 |
| Influenza A virus | GCF_000865725.1 | 0 | 0 |  |  | 2000 | 0 |
| Human immunodeficiency virus type 1 |  | 0 | 0 |  |  | 2000 |  |
| Severe acute respiratory syndrome coronavirus 2 | GCA_009858895.3 | 0 | 0 |  |  | 2000 | 0 |

Annotated genes with an intron beyond `-G` (5) — each still called, which is what D38 predicts (a split gene, not a lost one):

| species | family | gene_name | gene_span | annotated_max_intron | max_intron_used | locus |
|---|---|---|---|---|---|---|
| Aplysia californica | trpc | LOC101855406 | 381687 | 204025 | 200000 | NW_004797324.1:1395991-1513896+ |
| Branchiostoma floridae | asic | LOC118408435 | 273733 | 256433 | 200000 | NW_023365779.1:38557-149133- |
| Callorhinchus milii | asic | asic2 | 262382 | 223342 | 200000 | KI635942.1:2416117-2677415+ |
| Callorhinchus milii | trpm | trpm3 | 292682 | 208680 | 200000 | KI636021.1:118940-278611- |
| Homo sapiens | asic | ASIC2 | 1143682 | 1043910 | 1000000 | CM000679.2:33013965-34156532- |

Gene spans (4,484 annotated loci) → D4's bar (D38; 72 families measured):

| family | n_genes | n_species | median_span | min_span | max_span | by_group | by_group_n |
|---|---|---|---|---|---|---|---|
| ampa | 59 | 14 | 173907 | 25556 | 612057 | deuterostome:65612;vertebrate:175367 | deuterostome:1;vertebrate:58 |
| ano_channel | 32 | 13 | 138440 | 13313 | 568734 | vertebrate:138440 | vertebrate:32 |
| ano_scramblase | 118 | 28 | 36758 | 3353 | 618554 | basal_metazoan:5429;cnidarian:13225;deuterostome:18033;invertebrate:19750;plant:4165;vertebrate:57663 | basal_metazoan:3;cnidarian:7;deuterostome:9;invertebrate:14;plant:3;vertebrate:82 |
| asic | 74 | 22 | 28862 | 750 | 1143682 | basal_metazoan:4032;cnidarian:6034;deuterostome:15754;invertebrate:3807;vertebrate:120632 | basal_metazoan:4;cnidarian:5;deuterostome:13;invertebrate:5;vertebrate:47 |
| bestrophin | 69 | 21 | 7890 | 476 | 55796 | deuterostome:17890;invertebrate:2725;vertebrate:14788 | deuterostome:4;invertebrate:23;vertebrate:42 |
| calhm | 43 | 13 | 5585 | 693 | 32085 | vertebrate:5585 | vertebrate:43 |
| catsper | 27 | 13 | 12363 | 1315 | 26139 | basal_metazoan:1991;cnidarian:7750;deuterostome:8545;vertebrate:16405 | basal_metazoan:2;cnidarian:2;deuterostome:5;vertebrate:18 |
| cav | 185 | 29 | 111694 | 5255 | 996004 | algae:14442;basal_metazoan:12789;ciliate:5290;cnidarian:114056;deuterostome:97047;invertebrate:30011;vertebrate:158445 | algae:3;basal_metazoan:4;ciliate:3;cnidarian:12;deuterostome:10;invertebrate:14;vertebrate:139 |
| cftr | 13 | 13 | 114587 | 18949 | 204606 | vertebrate:114587 | vertebrate:13 |
| clc_channel | 46 | 25 | 27541 | 3353 | 214269 | basal_metazoan:7736;deuterostome:33992;holozoa:4435;invertebrate:9955;vertebrate:34098 | basal_metazoan:2;deuterostome:3;holozoa:3;invertebrate:9;vertebrate:29 |
| clcc1 | 13 | 13 | 24963 | 4401 | 50085 | vertebrate:24963 | vertebrate:13 |
| clic | 81 | 25 | 16767 | 600 | 248993 | basal_metazoan:1601;cnidarian:5831;deuterostome:16767;invertebrate:4552;plant:2367;vertebrate:29323 | basal_metazoan:1;cnidarian:1;deuterostome:3;invertebrate:9;plant:6;vertebrate:61 |
| cng | 85 | 27 | 13070 | 936 | 200372 | basal_metazoan:3533;ciliate:3239;cnidarian:5992;deuterostome:21941;invertebrate:12978;vertebrate:14114 | basal_metazoan:2;ciliate:3;cnidarian:5;deuterostome:7;invertebrate:16;vertebrate:52 |
| connexin | 220 | 14 | 7447 | 741 | 72007 | deuterostome:5816;vertebrate:7591 | deuterostome:11;vertebrate:209 |
| deg_invertebrate | 4 | 3 | 3758 | 1683 | 28686 | invertebrate:3758 | invertebrate:4 |
| delta_glur | 29 | 13 | 696754 | 69352 | 1915375 | vertebrate:696754 | vertebrate:29 |
| enac | 27 | 11 | 29250 | 7756 | 127086 | vertebrate:29250 | vertebrate:27 |
| gabaa | 156 | 18 | 62765 | 4748 | 491652 | deuterostome:47828;invertebrate:30264;vertebrate:65716 | deuterostome:7;invertebrate:2;vertebrate:147 |
| glyr | 47 | 13 | 81300 | 7109 | 646471 | vertebrate:81300 | vertebrate:47 |
| gphr | 36 | 32 | 8880 | 1432 | 97515 | algae:5352;amoebozoa:1834;basal_metazoan:6222;ciliate:1432;cnidarian:14048;deuterostome:8096;holozoa:2801;invertebrate:7283;plant:5664;vertebrate:28499 | algae:1;amoebozoa:1;basal_metazoan:3;ciliate:1;cnidarian:2;deuterostome:3;holozoa:2;invertebrate:6;plant:3;vertebrate:14 |
| hcn | 61 | 20 | 45326 | 4690 | 592243 | deuterostome:31861;invertebrate:45125;vertebrate:46649 | deuterostome:3;invertebrate:5;vertebrate:53 |
| ht3 | 29 | 13 | 9932 | 2883 | 24160 | vertebrate:9932 | vertebrate:29 |
| hv1 | 27 | 22 | 10311 | 769 | 96364 | basal_metazoan:1888;cnidarian:10434;deuterostome:7998;invertebrate:5196;vertebrate:12537 | basal_metazoan:2;cnidarian:2;deuterostome:4;invertebrate:7;vertebrate:12 |
| iglur_nonvertebrate | 15 | 5 | 4212 | 2427 | 6109 | invertebrate:4394;plant:4169 | invertebrate:2;plant:13 |
| innexin | 30 | 6 | 4016 | 1191 | 165980 | invertebrate:4016 | invertebrate:30 |
| itpr | 62 | 28 | 101695 | 9649 | 541015 | basal_metazoan:18715;cnidarian:93849;deuterostome:77899;holozoa:11281;invertebrate:44291;vertebrate:147444 | basal_metazoan:3;cnidarian:2;deuterostome:4;holozoa:4;invertebrate:8;vertebrate:41 |
| k2p | 231 | 30 | 11851 | 804 | 283451 | basal_metazoan:8930;cnidarian:3396;deuterostome:5302;holozoa:4004;invertebrate:5826;plant:2248;vertebrate:22353 | basal_metazoan:2;cnidarian:8;deuterostome:23;holozoa:1;invertebrate:31;plant:9;vertebrate:157 |
| kainate | 75 | 17 | 161469 | 4366 | 823858 | deuterostome:37748;invertebrate:24491;vertebrate:198374 | deuterostome:3;invertebrate:3;vertebrate:69 |
| kca_sk | 52 | 23 | 62515 | 1260 | 440519 | ciliate:1470;deuterostome:40970;invertebrate:7108;vertebrate:74463 | ciliate:1;deuterostome:3;invertebrate:9;vertebrate:39 |
| kca_slo | 82 | 30 | 76544 | 3107 | 1078759 | amoebozoa:4068;basal_metazoan:7064;ciliate:3139;cnidarian:108370;deuterostome:58137;holozoa:7394;invertebrate:29775;vertebrate:151752 | amoebozoa:2;basal_metazoan:3;ciliate:5;cnidarian:3;deuterostome:6;holozoa:1;invertebrate:13;vertebrate:49 |
| kir | 193 | 27 | 14751 | 861 | 309085 | basal_metazoan:1050;cnidarian:3526;deuterostome:12199;holozoa:1278;invertebrate:5829;vertebrate:17603 | basal_metazoan:3;cnidarian:5;deuterostome:8;holozoa:1;invertebrate:13;vertebrate:163 |
| kv_eag | 103 | 24 | 82813 | 2287 | 627641 | basal_metazoan:8744;cnidarian:14086;deuterostome:25276;invertebrate:27456;vertebrate:127008 | basal_metazoan:3;cnidarian:3;deuterostome:5;invertebrate:11;vertebrate:81 |
| kv_kcnq | 89 | 23 | 81838 | 1870 | 704504 | cnidarian:60647;deuterostome:42575;invertebrate:25067;vertebrate:124066 | cnidarian:4;deuterostome:7;invertebrate:9;vertebrate:69 |
| kv_modifier | 62 | 13 | 19092 | 1440 | 170021 | vertebrate:19092 | vertebrate:62 |
| kv_shaker | 186 | 25 | 8376 | 1071 | 603805 | basal_metazoan:1107;cnidarian:3184;deuterostome:14137;invertebrate:11198;vertebrate:8739 | basal_metazoan:1;cnidarian:17;deuterostome:9;invertebrate:12;vertebrate:147 |
| lrrc8 | 44 | 14 | 17509 | 1554 | 161222 | cnidarian:1554;vertebrate:17527 | cnidarian:1;vertebrate:43 |
| mcu | 38 | 27 | 42106 | 1784 | 195552 | basal_metazoan:1998;cnidarian:19175;deuterostome:39007;holozoa:4656;invertebrate:16207;plant:2461;vertebrate:70290 | basal_metazoan:1;cnidarian:2;deuterostome:4;holozoa:1;invertebrate:6;plant:1;vertebrate:23 |
| mitok | 16 | 16 | 7089 | 1547 | 30179 | deuterostome:3824;invertebrate:3159;vertebrate:10868 | deuterostome:1;invertebrate:2;vertebrate:13 |
| mscl | 6 | 6 | 424 | 375 | 545 | prokaryote:424 | prokaryote:6 |
| mscs | 15 | 6 | 3591 | 861 | 7207 | plant:3612;prokaryote:891 | plant:12;prokaryote:3 |
| nachr | 245 | 24 | 14153 | 1056 | 209967 | cnidarian:6371;deuterostome:11173;invertebrate:7440;vertebrate:19262 | cnidarian:22;deuterostome:37;invertebrate:48;vertebrate:138 |
| nalcn | 27 | 26 | 63865 | 2835 | 519213 | basal_metazoan:12785;cnidarian:31071;deuterostome:45603;invertebrate:14168;vertebrate:214682 | basal_metazoan:2;cnidarian:2;deuterostome:3;invertebrate:7;vertebrate:13 |
| nav | 100 | 24 | 73360 | 10281 | 339076 | basal_metazoan:12348;cnidarian:33439;deuterostome:64778;invertebrate:42020;vertebrate:97876 | basal_metazoan:2;cnidarian:5;deuterostome:10;invertebrate:10;vertebrate:73 |
| nmda | 86 | 22 | 51135 | 1950 | 555853 | cnidarian:12133;deuterostome:24985;invertebrate:18251;vertebrate:84997 | cnidarian:3;deuterostome:4;invertebrate:13;vertebrate:66 |
| orai | 43 | 25 | 8972 | 549 | 77135 | basal_metazoan:549;cnidarian:10320;deuterostome:3032;invertebrate:3798;vertebrate:12469 | basal_metazoan:1;cnidarian:2;deuterostome:3;invertebrate:7;vertebrate:30 |
| osca_tmem63 | 73 | 29 | 27853 | 2549 | 166446 | amoebozoa:2549;basal_metazoan:6209;cnidarian:11603;deuterostome:58197;fungi:2966;holozoa:3048;invertebrate:17488;plant:4970;vertebrate:42379 | amoebozoa:1;basal_metazoan:2;cnidarian:2;deuterostome:2;fungi:1;holozoa:1;invertebrate:5;plant:17;vertebrate:42 |
| otop | 50 | 22 | 9707 | 1128 | 74430 | cnidarian:2060;deuterostome:2527;invertebrate:3315;vertebrate:12653 | cnidarian:2;deuterostome:1;invertebrate:18;vertebrate:29 |
| p2x | 72 | 25 | 17913 | 1146 | 89348 | basal_metazoan:4685;cnidarian:27704;deuterostome:44420;holozoa:2549;invertebrate:12950;vertebrate:18385 | basal_metazoan:3;cnidarian:2;deuterostome:4;holozoa:2;invertebrate:4;vertebrate:57 |
| p2x_nonmetazoan | 2 | 2 | 7022 | 6156 | 7888 | basal_metazoan:6156;cnidarian:7888 | basal_metazoan:1;cnidarian:1 |
| pacc | 12 | 12 | 25219 | 8542 | 118533 | vertebrate:25219 | vertebrate:12 |
| pannexin | 39 | 13 | 16888 | 2997 | 69966 | vertebrate:16888 | vertebrate:39 |
| piezo | 48 | 25 | 91415 | 9690 | 592590 | amoebozoa:9690;basal_metazoan:13188;cnidarian:169344;deuterostome:88642;holozoa:9889;invertebrate:67905;vertebrate:100151 | amoebozoa:1;basal_metazoan:2;cnidarian:1;deuterostome:3;holozoa:1;invertebrate:4;vertebrate:36 |
| plgic_invertebrate | 7 | 3 | 7124 | 2284 | 44517 | invertebrate:7124 | invertebrate:7 |
| plgic_prok | 1 | 1 | 6302 | 6302 | 6302 | invertebrate:6302 | invertebrate:1 |
| ryr | 54 | 24 | 145702 | 15515 | 791805 | basal_metazoan:22879;deuterostome:129206;holozoa:20522;invertebrate:27734;vertebrate:184748 | basal_metazoan:2;deuterostome:3;holozoa:1;invertebrate:5;vertebrate:43 |
| tmc | 93 | 25 | 23313 | 1839 | 316690 | basal_metazoan:5830;cnidarian:11408;deuterostome:17773;invertebrate:15126;vertebrate:36081 | basal_metazoan:5;cnidarian:4;deuterostome:9;invertebrate:12;vertebrate:63 |
| tmco1 | 35 | 32 | 3835 | 652 | 54808 | amoebozoa:686;basal_metazoan:1652;ciliate:656;cnidarian:2692;deuterostome:12291;holozoa:948;invertebrate:1918;plant:1096;vertebrate:17476 | amoebozoa:1;basal_metazoan:2;ciliate:2;cnidarian:2;deuterostome:4;holozoa:1;invertebrate:6;plant:4;vertebrate:13 |
| tmem109 | 7 | 7 | 11087 | 7404 | 64975 | vertebrate:11087 | vertebrate:7 |
| tmem175 | 15 | 15 | 16107 | 2193 | 81178 | basal_metazoan:3112;cnidarian:2239;deuterostome:9929;vertebrate:17988 | basal_metazoan:1;cnidarian:2;deuterostome:1;vertebrate:11 |
| tmem87 | 37 | 28 | 22473 | 2255 | 96840 | amoebozoa:2255;basal_metazoan:4794;cnidarian:32484;deuterostome:21101;holozoa:2980;invertebrate:7747;vertebrate:27123 | amoebozoa:1;basal_metazoan:1;cnidarian:2;deuterostome:3;holozoa:2;invertebrate:6;vertebrate:22 |
| tpc | 51 | 28 | 23273 | 2494 | 121915 | amoebozoa:6260;basal_metazoan:4322;cnidarian:33929;deuterostome:26164;holozoa:3334;invertebrate:16324;plant:6555;vertebrate:42416 | amoebozoa:1;basal_metazoan:3;cnidarian:3;deuterostome:6;holozoa:2;invertebrate:7;plant:5;vertebrate:24 |
| tric | 32 | 20 | 16800 | 1728 | 100927 | deuterostome:6948;invertebrate:9078;vertebrate:25435 | deuterostome:1;invertebrate:7;vertebrate:24 |
| trpa | 29 | 22 | 38785 | 2482 | 365399 | basal_metazoan:3800;cnidarian:100591;deuterostome:30185;invertebrate:23384;vertebrate:44375 | basal_metazoan:1;cnidarian:1;deuterostome:7;invertebrate:6;vertebrate:14 |
| trpc | 125 | 25 | 50905 | 3846 | 381687 | cnidarian:7247;deuterostome:17250;holozoa:3846;invertebrate:13511;vertebrate:65400 | cnidarian:3;deuterostome:18;holozoa:1;invertebrate:19;vertebrate:84 |
| trpm | 118 | 26 | 64222 | 3818 | 1168467 | basal_metazoan:10162;cnidarian:40070;deuterostome:35630;holozoa:7433;invertebrate:24134;vertebrate:80408 | basal_metazoan:2;cnidarian:7;deuterostome:10;holozoa:2;invertebrate:15;vertebrate:82 |
| trpml | 44 | 24 | 21335 | 3094 | 184439 | basal_metazoan:3384;cnidarian:30979;deuterostome:13061;invertebrate:7401;vertebrate:25101 | basal_metazoan:1;cnidarian:5;deuterostome:3;invertebrate:5;vertebrate:30 |
| trpn | 14 | 11 | 64962 | 10559 | 133213 | cnidarian:99887;deuterostome:37462;invertebrate:17242;vertebrate:56454 | cnidarian:4;deuterostome:3;invertebrate:2;vertebrate:5 |
| trpp | 47 | 23 | 26667 | 3732 | 99280 | basal_metazoan:3732;cnidarian:58500;deuterostome:24676;invertebrate:17641;vertebrate:31176 | basal_metazoan:1;cnidarian:3;deuterostome:3;invertebrate:4;vertebrate:36 |
| trpv | 78 | 24 | 22142 | 3053 | 207776 | basal_metazoan:6678;cnidarian:16531;deuterostome:18332;invertebrate:5199;vertebrate:26298 | basal_metazoan:2;cnidarian:1;deuterostome:8;invertebrate:14;vertebrate:53 |
| tweety | 43 | 15 | 32817 | 6910 | 197591 | invertebrate:17878;vertebrate:33717 | invertebrate:2;vertebrate:41 |
| vdac | 74 | 29 | 10088 | 1077 | 142670 | basal_metazoan:3037;cnidarian:4768;deuterostome:7147;invertebrate:3622;plant:2904;vertebrate:16445 | basal_metazoan:3;cnidarian:2;deuterostome:4;invertebrate:13;plant:9;vertebrate:43 |
| zac | 2 | 2 | 6445 | 3625 | 9265 | vertebrate:6445 | vertebrate:2 |

## Runs

| species | assembly | total_bp | n50 | max_intron | chunks | alignments | self_dropped | loci | family_called_loci | miniprot_s | call_s | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | GCA_000005845.2 | 4641652 | 4641652 | 2000 | 1 | 62 | 7 | 13 | 8 | 1.7 | 1.5 |  |
| Bacillus subtilis | GCA_000009045.1 | 4215606 | 4215606 | 2000 | 1 | 95 | 5 | 12 | 9 | 1.8 | 1.1 |  |
| Streptomyces lividans | GCA_000739105.1 | 8345283 | 8345283 | 2000 | 1 | 81 | 6 | 30 | 12 | 2.7 | 1.1 |  |
| Synechocystis sp. PCC 6803 | GCA_000009725.1 | 3947019 | 3573470 | 2000 | 1 | 90 | 9 | 26 | 13 | 1.4 | 0.9 |  |
| Gloeobacter violaceus | GCA_000011385.1 | 4659019 | 4659019 | 2000 | 1 | 53 | 5 | 19 | 10 | 1.8 | 0.7 |  |
| Aliarcobacter butzleri | GCA_000014025.1 | 2341251 | 2341251 | 2000 | 1 | 81 | 4 | 13 | 6 | 0 | 0.8 |  |
| Methanothermobacter thermautotrophicus | GCA_000008645.1 | 1751377 | 1751377 | 2000 | 1 | 14 | 3 | 9 | 1 | 0 | 0.4 |  |
| Thermus thermophilus | GCA_000091545.1 | 2116056 | 1849742 | 2000 | 1 | 31 | 3 | 11 | 3 | 0 | 0.8 |  |
| Saccharomyces cerevisiae | GCA_000146045.2 | 12071326 | 924431 | 200000 | 1 | 338 | 8 | 54 | 19 | 20.8 | 5.0 |  |
| Schizosaccharomyces pombe | GCA_000002945.3 | 12571820 | 4539804 | 200000 | 1 | 237 | 9 | 37 | 16 | 14.9 | 3.2 |  |
| Arabidopsis thaliana | GCA_000001735.2 | 119668634 | 23459830 | 200000 | 1 | 985 | 27 | 166 | 105 | 45.7 | 11.9 |  |
| Oryza sativa | GCA_001433935.1 | 373795655 | 29958434 | 200000 | 1 | 960 | 17 | 217 | 102 | 59.9 | 10.5 |  |
| Physcomitrium patens | GCA_000002425.3 | 471971527 | 17432721 | 200000 | 1 | 855 | 28 | 141 | 85 | 53.3 | 9.2 |  |
| Chlamydomonas reinhardtii | GCA_000002595.3 | 111098438 | 7783580 | 200000 | 1 | 657 | 10 | 216 | 28 | 52.1 | 13.6 |  |
| Dictyostelium discoideum | GCA_000004695.1 | 34134454 | 5450249 | 200000 | 1 | 650 | 16 | 113 | 33 | 41.0 | 10.5 |  |
| Paramecium tetraurelia | GCA_000165425.1 | 72094543 | 413286 | 200000 | 1 | 559 | 29 | 141 | 68 | 43.9 | 11.0 |  |
| Tetrahymena thermophila | GCA_000189635.1 | 103014375 | 520943 | 200000 | 1 | 400 | 17 | 85 | 39 | 30.0 | 6.3 |  |
| Trypanosoma brucei | GCA_000002445.1 | 26075494 | 2481190 | 200000 | 1 | 207 | 7 | 49 | 11 | 25.2 | 3.1 |  |
| Plasmodium falciparum | GCA_000002765.3 | 23332839 | 1687656 | 200000 | 1 | 158 | 5 | 67 | 5 | 20.3 | 2.6 |  |
| Monosiga brevicollis | GCA_000002865.1 | 41633360 | 1073601 | 200000 | 1 | 653 | 23 | 112 | 46 | 40.3 | 10.4 |  |
| Capsaspora owczarzaki | GCA_000151315.2 | 27967784 | 1617775 | 200000 | 1 | 640 | 25 | 138 | 31 | 42.1 | 10.3 |  |
| Amphimedon queenslandica | GCA_000090795.2 | 164242647 | 123180 | 200000 | 1 | 1060 | 35 | 145 | 79 | 34.7 | 14.5 |  |
| Trichoplax adhaerens | GCA_000150275.1 | 105631681 | 5978658 | 200000 | 1 | 1559 | 34 | 156 | 115 | 37.2 | 20.6 |  |
| Nematostella vectensis | GCA_000209225.1 | 356613585 | 472588 | 200000 | 1 | 2359 | 43 | 328 | 243 | 64.4 | 32.3 |  |
| Hydra vulgaris | GCF_038396675.1 | 911675201 | 61145917 | 200000 | 1 | 1518 | 51 | 230 | 146 | 84.4 | 21.9 |  |
| Caenorhabditis elegans | GCA_000002985.3 | 100272607 | 17493829 | 200000 | 1 | 1705 | 47 | 196 | 156 | 44.6 | 21.7 |  |
| Drosophila melanogaster | GCA_000001215.4 | 143726002 | 25286936 | 200000 | 1 | 1930 | 51 | 182 | 124 | 48.3 | 20.6 |  |
| Daphnia pulex | GCA_000187875.1 | 197206209 | 642089 | 200000 | 1 | 1823 | 55 | 188 | 145 | 42.2 | 16.8 |  |
| Lottia gigantea | GCA_000327385.1 | 359505668 | 1870055 | 200000 | 1 | 1996 | 49 | 238 | 191 | 45.1 | 19.9 |  |
| Aplysia californica | GCF_000002075.1 | 927310431 | 917541 | 200000 | 1 | 1871 | 64 | 289 | 156 | 79.1 | 25.3 |  |
| Lymnaea stagnalis | GCA_964033795.1 | 942996421 | 957215 | 200000 | 1 | 2162 | 50 | 293 | 197 | 64.1 | 32.1 |  |
| Cornu aspersum | GCA_964187895.1 | 2908462832 | 110296109 | 1000000 | 2 | 3828 | 0 | 692 | 276 | 746.9 | 75.4 |  |
| Strongylocentrotus purpuratus | GCA_000002235.4 | 921840143 | 37282239 | 200000 | 1 | 2084 | 64 | 298 | 200 | 77.9 | 31.6 |  |
| Branchiostoma floridae | GCF_000003815.2 | 513460931 | 25441410 | 200000 | 1 | 2328 | 66 | 287 | 239 | 61.0 | 32.8 |  |
| Ciona intestinalis | GCA_000224145.2 | 115226814 | 5152901 | 200000 | 1 | 1605 | 49 | 162 | 137 | 33.6 | 16.3 |  |
| Petromyzon marinus | GCF_048934315.1 | 1808348452 | 13495987 | 1000000 | 1 | 3559 | 75 | 519 | 423 | 270.6 | 54.3 |  |
| Callorhinchus milii | GCA_000165045.2 | 974481817 | 4521921 | 200000 | 1 | 2676 | 82 | 322 | 294 | 52.1 | 30.4 |  |
| Torpedo marmorata | GCA_977017645.1 | 6221104460 | 111696861 | 1000000 | 3 | 3836 | 0 | 541 | 326 | 897.0 | 48.4 |  |
| Danio rerio | GCF_049306965.2 | 1448791966 | 59305174 | 1000000 | 1 | 4356 | 92 | 504 | 453 | 279.7 | 67.3 |  |
| Takifugu rubripes | GCA_901000725.3 | 384126663 | 16705553 | 200000 | 1 | 4400 | 89 | 525 | 477 | 60.4 | 55.3 |  |
| Latimeria chalumnae | GCA_000225785.1 | 2860575514 | 924513 | 1000000 | 2 | 3327 | 121 | 427 | 365 | 283.3 | 42.7 |  |
| Xenopus tropicalis | GCF_000004195.4 | 1451301209 | 153961319 | 1000000 | 1 | 3128 | 84 | 373 | 332 | 208.2 | 40.5 |  |
| Anolis carolinensis | GCA_000090745.2 | 1799143587 | 150641573 | 1000000 | 1 | 2829 | 82 | 377 | 309 | 200.6 | 41.8 |  |
| Gallus gallus | GCA_016699485.1 | 1053332251 | 90861225 | 1000000 | 1 | 2890 | 85 | 334 | 286 | 186.7 | 38.3 |  |
| Ornithorhynchus anatinus | GCA_004115215.4 | 1859281749 | 83338043 | 1000000 | 1 | 2872 | 85 | 343 | 308 | 218.2 | 45.0 |  |
| Monodelphis domestica | GCA_000002295.1 | 3600487649 | 527952102 | 1000000 | 2 | 3674 | 149 | 521 | 342 | 619.6 | 50.4 |  |
| Mus musculus | GCA_000001635.9 | 2728222451 | 130530862 | 1000000 | 2 | 3402 | 117 | 578 | 332 | 614.1 | 55.3 |  |
| Rattus norvegicus | GCA_036323735.1 | 2849597507 | 137014596 | 1000000 | 2 | 3621 | 124 | 613 | 332 | 590.4 | 58.6 |  |
| Homo sapiens | GCA_000001405.29 | 3298912062 | 145138636 | 1000000 | 2 | 4201 | 177 | 548 | 415 | 530.0 | 67.1 |  |
| Influenza A virus | GCA_000865725.1 | 13588 | 2233 | 2000 | 1 | 0 | 1 | 0 | 0 | 0 | 0.1 |  |
| Human immunodeficiency virus type 1 | GCA_003102975.1 | 9719 | 9719 | 2000 | 1 | 0 | 1 | 0 | 0 | 0 | 0.0 |  |
| Severe acute respiratory syndrome coronavirus 2 | GCA_009858895.3 | 29903 | 29903 | 2000 | 1 | 0 | 1 | 0 | 0 | 0 | 0.0 |  |
