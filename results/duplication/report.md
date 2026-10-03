# S11a — Duplication history by reconciliation

Rendered by `scripts/s11_report.py` from the tables beside it (D13). Rules: **D52**, fixed before any tree was reconciled. Driver: `scripts/s11_duplications.py` (`reconcile`, `pairs`); logic `scripts/s11_recon.py`; self-test `selftest_s11.py`.

## 1. What was reconciled

S7b's **63 tier-1 gene trees** (6,646 ingroup tips, outgroups pruned) against the NCBI Taxonomy tree of the 52 S4 panel species (D15; polytomies kept — Chordata is a three-way star of amphioxus, *Ciona* and vertebrates, because NCBI has no Olfactores). Duplications by species overlap; losses (tie-break only) by LCA mapping. Families with no tier-1 tree (the five < 4-sequence families and the r4 families, which S6 predates) are outside this task.

## 2. Roots

* **28** trees keep their S7b root (outgroup one clade, UFBoot ≥ 95 or one sequence, not an S7d family). The declared root is itself a reconciliation optimum in **6 / 28** — the outgroup root and the fewest-duplications root mostly disagree (reported, never replaced).
* **35** trees are rooted by reconciliation (no declared outgroup, a split or weak outgroup, or S7d's six); the optimum is a single edge in **27**. A duplication that holds under only some tied optima is not stated: 56 such nodes in all.

The six roots S7d left undefined, as reconciliation places them (for S11b to compare with an outgroup-free root):

| family | optimal roots | basal split | smaller side (species) |
|---|---|---|---|
| cng | 1 | 142|171 | Danio rerio 11; Takifugu rubripes 9; Gallus gallus 7; Rattus norvegicus 6 |
| k2p | 1 | 98|282 | Caenorhabditis elegans 38; Strongylocentrotus purpuratus 8; Drosophila melanogaster 7; Lottia gigantea 6 |
| kca_slo | 1 | 37|69 | Hydra vulgaris 4; Gallus gallus 3; Homo sapiens 2; Rattus norvegicus 2 |
| kv_eag | 1 | 7|151 | Paramecium tetraurelia 5; Tetrahymena thermophila 2 |
| kv_kcnq | 1 | 29|68 | Anolis carolinensis 2; Xenopus tropicalis 2; Takifugu rubripes 2; Danio rerio 2 |
| kv_shaker | 7 | – | tied |

## 3. Duplications

**2,316 stated duplications, 1,374 supported** (each child a leaf or UFBoot ≥ 95). Of the supported, **1,117 are species-specific** (36 involve a genome-locus tip) and **257 are shared by more than one species**.

Supported duplications shared by more than one species, by mapped taxon:

| taxon | supported |
|---|---|
| Clupeocephala | 58 |
| Gnathostomata | 52 |
| Gastropoda | 32 |
| Vertebrata | 24 |
| Oligohymenophorea | 11 |
| Euthyneura | 10 |
| Allotriocarida | 7 |
| Eumetazoa | 6 |
| Embryophyta | 6 |
| Chordata | 6 |
| Metazoa | 5 |
| Opisthokonta | 5 |
| Bilateria | 4 |
| Cnidaria | 4 |
| Eukaryota | 3 |

Species-specific supported duplications, top species (a species-specific node can also be an allelic or redundant gene model; the source of every tip is in `duplications.tsv`):

| species | supported |
|---|---|
| *Petromyzon marinus* | 122 |
| *Branchiostoma floridae* | 120 |
| *Caenorhabditis elegans* | 106 |
| *Paramecium tetraurelia* | 96 |
| *Strongylocentrotus purpuratus* | 90 |
| *Hydra vulgaris* | 63 |
| *Tetrahymena thermophila* | 59 |
| *Nematostella vectensis* | 39 |
| *Physcomitrium patens* | 34 |
| *Daphnia pulex* | 33 |
| *Takifugu rubripes* | 31 |
| *Gallus gallus* | 30 |

Teleost 3R window (mapped to Clupeocephala, *Danio* + *Takifugu*): cav 7, connexin 5, kv_shaker 5, kv_eag 4, nachr 4, nmda 4, k2p 3, kir 3, ampa 2, cng 2, gabaa 2, nav 2, trpml 2, glyr 1, hcn 1, itpr 1, kainate 1, kca_sk 1, kca_slo 1, kv_modifier 1, lrrc8 1, piezo 1, ryr 1, trpc 1, trpm 1, trpv 1.

2R window (Vertebrata or Gnathostomata): connexin 6, gabaa 5, k2p 5, kv_shaker 5, kir 4, lrrc8 4, nachr 4, kca_slo 3, kv_modifier 3, nmda 3, p2x 3, trpm 3, calhm 2, clc_channel 2, kainate 2, kv_eag 2, kv_kcnq 2, piezo 2, ano_channel 1, asic 1, cav 1, cng 1, enac 1, glyr 1, hcn 1, itpr 1, kca_sk 1, osca_tmem63 1, otop 1, ryr 1, tric 1, trpc 1, trpp 1, trpv 1.

## 4. Human paralogue pairs and OHNOLOGS v2

1,279 human paralogue pairs inside one family, by the window of the stated duplication separating them: 2R 494, bony_vertebrate 4, none 99, older 657, younger 25.

**OHNOLOGS v2 strict 2R pairs inside our families: 303. The gene trees date 176 (58 %) to the 2R window and 124 (41 %) older than vertebrates**; the rest: bony_vertebrate 1, none 2. No OHNOLOGS pair joins two of our families (strict 0, relaxed 0).

The too-old pairs concentrate in a few families (kir 53, nav 22, gabaa 18, trpm 7, ampa 6, hcn 5, kainate 5, trpv 4, trpp 2, mcu 1, trpc 1) and map to Chordata 50, cellular organisms 28, Eumetazoa 22, Bilateria 11, Opisthokonta 8, Eukaryota 5. Their duplications have a median species-overlap score of 0.6 against 0.8901 for the pairs dated to 2R: a non-vertebrate sequence nested among vertebrate paralogues (gene-tree incongruence, or a real pre-vertebrate duplication) drags the node to a deeper taxon. D52 is not changed for it; a synteny-free tree reconciliation dates 2R reliably only where the vertebrate paralogue clades are clean.

Per gene (286 human genes on tier-1 trees): OHNOLOGS strict lists 183 as 2R ohnologues; our trees give 158 of them (86 %) a 2R-window paralogue. Our trees give 241 genes a 2R-window paralogue, 158 (66 %) of them on the strict list (OHNOLOGS demands synteny with outgroup genomes).

## 5. Limits

* Reconciliation trusts the gene tree. Supported (UFBoot ≥ 95) child edges do not protect against a misplaced long branch: § 4 measures how often that dates a synteny-supported 2R pair too deep.
* The species tree is NCBI's, with its polytomies (D15); a duplication mapped to a polytomy node is dated to that node, not resolved.
* Copy number is the D39 set: proteome genes (one per gene) and intact genome loci. Species-specific duplications include redundant models and genome-locus splits; they are counted apart from the headline.
* The repeat-duplication order of the 4×6TM channels and outgroup-free roots are S11b's.

## Files

`recon_trees.tsv` (per tree: root mode, optima, declared-root agreement, counts, reconciliation basal split) · `duplications.tsv` (every stated duplication) · `expansions.tsv` (supported per family × taxon) · `human_pairs.tsv` · `human_genes.tsv` · `ohnolog_compare.tsv` · `summary.json` · figure `figures/duplications.png`. Archives under `<data root>/raw_api/s11/` (NCBI Taxonomy, HGNC complete set, OHNOLOGS v2 human 2R pairs).
