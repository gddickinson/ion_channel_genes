# S4 — Proteome scope (the denominator)

UniProt reference proteomes, release **2026_03** (README dated 2026-09-03) — the release census v2 was enumerated from (D31).

## Headline

- **52 panel species** (`src/utils/species.py`) across **15 panel groups**: **50 reference proteomes** (47 cellular, 3 viral), **2 genome-only**, **0 with neither**.
- **822,499 canonical entries**, one per gene: the entry count equals the proteome's declared gene count in **50/50** proteomes.
- **Verified: 49/50 exact** against the release README (canonical entries = #(1), gene2acc rows = #(3)) **and** the per-proteome `RELEASE.metalink` MD5; 1 reissued (Homo sapiens), 0 failed.
- **Positive control:** 319/319 enumerated human census genes (census v2 `human_recall.tsv`) are entries of the human reference proteome.
- **Census v2 overlap:** 12,096 of census v2's 1,245,200 records are panel-proteome entries — the starting point S3b's profile sweep is measured against.
- Sweep DB for S3b: 822,499 sequences, 507 MB, SHA-256 `366dd895ad6b9a0c…` (`/Volumes/FANTOM/ION_CHANNEL_DATA/proteomes/s4/panel_refprot.fasta`); downloads 268 MB.

## The three rules

1. **Candidates**: reference proteomes in the panel taxon's subtree (strains included) that the pinned release ships.
2. **Selection**: most Swiss-Prot entries → BUSCO complete % → gene count → UPID (`s4_proteome_lib.select_proteome`).
3. **No reference proteome** → the species' best current NCBI assembly (annotated > RefSeq > level > scaffold N50), status `genome_only`: its families can then be found only by S5's genomic sweep.

Species with more than one candidate, and what the rule chose:

| species | UPID | organism | selected | reviewed | BUSCO % | genes |
|---|---|---|---|---|---|---|
| Thermus thermophilus | UP000000532 | Thermus thermophilus (strain ATCC 27634 / DSM 579 / HB8) | **yes** | 441 | 87.9 | 2,227 |
|  | UP000000592 | Thermus thermophilus (strain ATCC BAA-163 / DSM 7039 / HB27) |  | 378 | 86.3 | 2,200 |
| Saccharomyces cerevisiae | UP000002311 | Saccharomyces cerevisiae (strain ATCC 204508 / S288c) | **yes** | 6,067 | 99.6 | 6,066 |
|  | UP000008335 | Saccharomyces cerevisiae (strain RM11-1a) |  | 199 | 95.0 | 5,357 |
| Oryza sativa | UP000059680 | Oryza sativa subsp. japonica | **yes** | 4,184 | 84.3 | 43,678 |
|  | UP000007015 | Oryza sativa subsp. indica |  | 751 | 93.4 | 37,431 |
| Plasmodium falciparum | UP000001450 | Plasmodium falciparum (isolate 3D7) | **yes** | 318 | 99.1 | 5,361 |
|  | UP000030673 | Plasmodium falciparum (isolate NF54) |  | 14 | 80.0 | 5,887 |
| Influenza A virus | UP000009255 | Influenza A virus (strain A/Puerto Rico/8/1934 H1N1) | **yes** | 13 | — | 13 |
|  | UP000131152 | Influenza A virus (strain A/Goose/Guangdong/1/1996 H5N1 genotype Gs/Gd) |  | 13 | — | 13 |
| Human immunodeficiency virus type 1 | UP000002241 | Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) | **yes** | 9 | — | 9 |
|  | UP000007420 | Human immunodeficiency virus type 1 group N (isolate YBF30) |  | 9 | — | 9 |
|  | UP000007689 | Human immunodeficiency virus type 1 group O (isolate ANT70) |  | 9 | — | 9 |
|  | UP000105453 | Human immunodeficiency virus type 1 group M subtype B (isolate HXB2) |  | 5 | — | 9 |
|  | UP000134285 | Human immunodeficiency virus type 1 group M subtype A (isolate U455) |  | 6 | — | 6 |
|  | UP000107373 | HIV-1 CRF03_AB |  | 0 | — | 9 |

## Per group

| group | species | proteomes | genome-only / none | genes | census v2 records |
|---|---|---|---|---|---|
| algae | 1 | 1 | 0 | 17,625 | 145 |
| amoebozoa | 1 | 1 | 0 | 12,718 | 49 |
| apicomplexa | 1 | 1 | 0 | 5,361 | 10 |
| basal_metazoan | 2 | 2 | 0 | 31,692 | 271 |
| ciliate | 2 | 2 | 0 | 66,325 | 1365 |
| cnidarian | 2 | 2 | 0 | 46,896 | 845 |
| deuterostome | 3 | 3 | 0 | 69,757 | 1290 |
| excavate | 1 | 1 | 0 | 8,561 | 27 |
| fungi | 2 | 2 | 0 | 11,195 | 30 |
| holozoa | 2 | 2 | 0 | 17,894 | 120 |
| invertebrate | 7 | 6 | 1 | 129,208 | 2081 |
| plant | 3 | 3 | 0 | 102,530 | 351 |
| prokaryote | 8 | 8 | 0 | 30,384 | 189 |
| vertebrate | 14 | 13 | 1 | 272,314 | 5323 |
| virus | 3 | 3 | 0 | 39 | 0 |

## Manifest

Annotation source is carried per row (D9): a RefSeq gene set and a submitter gene set are not comparable evidence of absence.

| species | group | proteome | genes | BUSCO % | CPD | annotation | assembly | level | in census v2 | verified |
|---|---|---|---|---|---|---|---|---|---|---|
| Escherichia coli | prokaryote | UP000000625 | 4,403 | 100.0 | Unknown | ENA/EMBL | GCA_000005845.2 | Complete Genome | 27 | exact |
| Bacillus subtilis | prokaryote | UP000001570 | 4,288 | 99.6 | Unknown | ENA/EMBL | GCA_000009045.1 | Complete Genome | 20 | exact |
| Streptomyces lividans | prokaryote | UP000028682 | 7,430 | 99.7 | Unknown | ENA/EMBL | GCA_000739105.1 | Complete Genome | 25 | exact |
| Synechocystis sp. PCC 6803 | prokaryote | UP000001425 | 3,508 | 96.2 | Unknown | ENA/EMBL | GCA_000009725.1 | Complete Genome | 49 | exact |
| Gloeobacter violaceus | prokaryote | UP000000557 | 4,406 | 91.8 | Unknown | ENA/EMBL | GCA_000011385.1 | Complete Genome | 26 | exact |
| Aliarcobacter butzleri | prokaryote | UP000001136 | 2,254 | 98.2 | Unknown | ENA/EMBL | GCA_000014025.1 | Complete Genome | 22 | exact |
| Methanothermobacter thermautotrophicus | prokaryote | UP000005223 | 1,868 | 96.5 | Unknown | ENA/EMBL | GCA_000008645.1 | Complete Genome | 9 | exact |
| Thermus thermophilus | prokaryote | UP000000532 | 2,227 | 87.9 | Unknown | ENA/EMBL | GCA_000091545.1 | Complete Genome | 11 | exact |
| Saccharomyces cerevisiae | fungi | UP000002311 | 6,066 | 99.6 | Unknown | ENA/EMBL | GCA_000146045.2 | Complete Genome | 15 | exact |
| Schizosaccharomyces pombe | fungi | UP000002485 | 5,129 | 81.8 | Unknown | ENA/EMBL | GCA_000002945.2 (previous; now GCA_000002945.3) | Chromosome | 15 | exact |
| Arabidopsis thaliana | plant | UP000006548 | 27,496 | 100.0 | Unknown | ENA/EMBL | GCA_000001735.1 (previous; now GCA_000001735.2) | Chromosome | 120 | exact |
| Oryza sativa | plant | UP000059680 | 43,678 | 84.3 | Unknown | ENA/EMBL | GCA_001433935.1 | Chromosome | 104 | exact |
| Physcomitrium patens | plant | UP000006727 | 31,356 | 88.7 | Unknown | EnsemblPlants | GCA_000002425.2 (previous; now GCA_000002425.3) | Chromosome | 127 | exact |
| Chlamydomonas reinhardtii | algae | UP000006906 | 17,625 | 98.9 | Unknown | ENA/EMBL | GCA_000002595.3 | Chromosome | 145 | exact |
| Dictyostelium discoideum | amoebozoa | UP000002195 | 12,718 | 93.7 | Unknown | ENA/EMBL | GCA_000004695.1 | Chromosome | 49 | exact |
| Paramecium tetraurelia | ciliate | UP000000600 | 39,350 | 98.8 | Outlier (high value) | ENA/EMBL | GCA_000165425.1 | Scaffold | 765 | exact |
| Tetrahymena thermophila | ciliate | UP000009168 | 26,975 | 98.2 | Unknown | ENA/EMBL | GCA_000189635.1 | Scaffold | 600 | exact |
| Trypanosoma brucei | excavate | UP000008524 | 8,561 | 97.7 | Unknown | ENA/EMBL | GCA_000002445.1 | Chromosome | 27 | exact |
| Plasmodium falciparum | apicomplexa | UP000001450 | 5,361 | 99.1 | Unknown | ENA/EMBL | GCA_000002765.3 | Complete Genome | 10 | exact |
| Monosiga brevicollis | holozoa | UP000001357 | 9,156 | 78.8 | Unknown | ENA/EMBL | GCA_000002865.1 | Scaffold | 73 | exact |
| Capsaspora owczarzaki | holozoa | UP000008743 | 8,738 | 93.7 | Unknown | ENA/EMBL | GCA_000151315.2 | Scaffold | 47 | exact |
| Amphimedon queenslandica | basal_metazoan | UP000007879 | 20,174 | 90.8 | Unknown | EnsemblMetazoa | GCA_000090795.2 | Scaffold | 136 | exact |
| Trichoplax adhaerens | basal_metazoan | UP000009022 | 11,518 | 89.3 | Unknown | ENA/EMBL | GCA_000150275.1 | Scaffold | 135 | exact |
| Nematostella vectensis | cnidarian | UP000001593 | 24,428 | 89.7 | Unknown | ENA/EMBL | GCA_000209225.1 | Scaffold | 510 | exact |
| Hydra vulgaris | cnidarian | UP001652625 | 22,468 | 94.7 | Unknown | Refseq | GCF_038396675.1 | Complete Genome | 335 | exact |
| Caenorhabditis elegans | invertebrate | UP000001940 | 19,792 | 100.0 | Unknown | ENA/EMBL | GCA_000002985.3 | Complete Genome | 369 | exact |
| Drosophila melanogaster | invertebrate | UP000000803 | 13,814 | 100.0 | Unknown | ENA/EMBL | GCA_000001215.4 | Chromosome | 199 | exact |
| Daphnia pulex | invertebrate | UP000000305 | 30,118 | 94.7 | Unknown | ENA/EMBL | GCA_000187875.1 | Scaffold | 381 | exact |
| Lottia gigantea | invertebrate | UP000030746 | 23,671 | 88.4 | Unknown | ENA/EMBL | GCA_000327385.1 | Scaffold | 380 | exact |
| Aplysia californica | invertebrate | UP000694888 | 19,368 | 97.3 | Unknown | Refseq | GCF_000002075.1 | Scaffold | 384 | exact |
| Lymnaea stagnalis | invertebrate | UP001497497 | 22,445 | 91.7 | Unknown | ENA/EMBL | GCA_964033795.1 | Scaffold | 368 | exact |
| Cornu aspersum | invertebrate | *genome_only* | — | — | — | — | GCA_964187895.1 | Chromosome | — | — |
| Strongylocentrotus purpuratus | deuterostome | UP000007110 | 26,461 | 99.4 | Unknown | EnsemblMetazoa | GCA_000002235.4 | Scaffold | 371 | exact |
| Branchiostoma floridae | deuterostome | UP000001554 | 26,616 | 97.5 | Unknown | Refseq | GCF_000003815.2 | Chromosome | 744 | exact |
| Ciona intestinalis | deuterostome | UP000008144 | 16,680 | 78.7 | Unknown | Ensembl | GCA_000224145.1 (previous; now GCA_000224145.2) | Chromosome | 175 | exact |
| Petromyzon marinus | vertebrate | UP001735300 | 21,315 | 81.1 | Unknown | Refseq | GCF_048934315.1 | Chromosome | 497 | exact |
| Callorhinchus milii | vertebrate | UP000314986 | 19,409 | 90.1 | Close to standard (high value) | Ensembl | GCA_000165045.2 | Scaffold | 380 | exact |
| Torpedo marmorata | vertebrate | *genome_only* | — | — | — | — | GCA_977017645.1 | Chromosome | — | — |
| Danio rerio | vertebrate | UP000000437 | 26,567 | 99.1 | Unknown | Refseq | GCF_049306965.1 (previous; now GCF_049306965.2) | Complete Genome | 550 | exact |
| Takifugu rubripes | vertebrate | UP000005226 | 21,257 | 92.8 | Unknown | Ensembl | GCA_901000725.2 (previous; now GCA_901000725.3) | Chromosome | 491 | exact |
| Latimeria chalumnae | vertebrate | UP000008672 | 19,568 | 89.1 | Unknown | Ensembl | GCA_000225785.1 | Scaffold | 420 | exact |
| Xenopus tropicalis | vertebrate | UP000008143 | 21,732 | 96.8 | Unknown | Refseq | GCF_000004195.4 | Chromosome | 405 | exact |
| Anolis carolinensis | vertebrate | UP000001646 | 21,494 | 86.3 | Unknown | Ensembl | GCA_000090745.2 | Chromosome | 363 | exact |
| Gallus gallus | vertebrate | UP000000539 | 18,373 | 98.2 | Unknown | Ensembl | GCA_016699485.1 | Chromosome | 373 | exact |
| Ornithorhynchus anatinus | vertebrate | UP000002279 | 17,391 | 85.6 | Unknown | Ensembl | GCA_004115215.2 (previous; now GCA_004115215.4) | Chromosome | 365 | exact |
| Monodelphis domestica | vertebrate | UP000002280 | 21,225 | 87.4 | Unknown | Ensembl | GCA_000002295.1 | Chromosome | 370 | exact |
| Mus musculus | vertebrate | UP000000589 | 21,860 | 99.8 | Unknown | Ensembl | GCA_000001635.9 | Chromosome | 366 | exact |
| Rattus norvegicus | vertebrate | UP000002494 | 21,471 | 95.9 | Unknown | Ensembl | GCA_036323735.1 | Chromosome | 368 | exact |
| Homo sapiens | vertebrate | UP000005640 | 20,652 | 99.5 | Unknown | Ensembl | GCA_000001405.29 | Chromosome | 375 | reissued |
| Influenza A virus | virus | UP000009255 | 13 | — | Unknown | ENA/EMBL | GCA_000865725.1 | Complete Genome | 0 | exact |
| Human immunodeficiency virus type 1 | virus | UP000002241 | 9 | — | Unknown | ENA/EMBL | GCA_003102975.1 | Complete Genome | 0 | exact |
| Severe acute respiratory syndrome coronavirus 2 | virus | UP000464024 | 17 | — | Unknown | ENA/EMBL | GCA_009858895.3 | Complete Genome | 0 | exact |

## What S4 cannot decide — handed on

- **14 proteomes below 90 % BUSCO complete**: Thermus thermophilus (87.9), Schizosaccharomyces pombe (81.8), Oryza sativa (84.3), Physcomitrium patens (88.7), Monosiga brevicollis (78.8), Trichoplax adhaerens (89.3), Nematostella vectensis (89.7), Lottia gigantea (88.4), Ciona intestinalis (78.7), Petromyzon marinus (81.1), Latimeria chalumnae (89.1), Anolis carolinensis (86.3), Ornithorhynchus anatinus (85.6), Monodelphis domestica (87.4). An absence in one of these is weaker evidence; S5's D4 bars decide, this table only reports.
- **7 proteomes are annotated on a superseded assembly version** (the manifest names the current one): Schizosaccharomyces pombe, Arabidopsis thaliana, Physcomitrium patens, Ciona intestinalis, Danio rerio, Takifugu rubripes, Ornithorhynchus anatinus. S5 searches the assembly the gene set was built on, or says which it used.
- **14 assemblies are scaffold-level or below**: Paramecium tetraurelia, Tetrahymena thermophila, Monosiga brevicollis, Capsaspora owczarzaki, Amphimedon queenslandica, Trichoplax adhaerens, Nematostella vectensis, Daphnia pulex, Lottia gigantea, Aplysia californica, Lymnaea stagnalis, Strongylocentrotus purpuratus, Callorhinchus milii, Latimeria chalumnae.
- **Genome-only species**: Cornu aspersum — GCA_964187895.1 (Chromosome, annotated=0); Torpedo marmorata — GCA_977017645.1 (Chromosome, annotated=0). Neither has a reference proteome; both are catalogue exemplar sources.
- **The viral proteomes hold no census v2 record**: the viroporin signatures are declared `SUBFAMILY` and so were never in the enumerated search space (D34). S3b's viroporin profile reaches these three proteomes; the census-space gap is an emergent row.
- **The human proteome was reissued mid-release** (served file dated 2026-09-15): the README's counts describe the withdrawn file (147,503 entries against 20,652 served — the whole human UniProtKB set, not one per gene). The served file matches the metalink MD5 and the declared gene count, and is used (D35).

