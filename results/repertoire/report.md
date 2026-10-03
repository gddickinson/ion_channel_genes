# S10a — Repertoire reconstruction

*Rendered by `scripts/s10_report.py` from the tables beside it (D13). Rules: D50, fixed before any family was reconstructed; addenda (a)–(c) labelled.*

## Headline

Every census family's presence was reconstructed on **437 eukaryotic orders** (S4b, one reference proteome each) over the NCBI Taxonomy tree (235 internal nodes, 93 polytomies, a 10-way star at the eukaryotic root). Under the primary cost (gain = 2 losses) **19 of 75 families are reconstructed present in the last eukaryotic common ancestor (LECA)**, 13 of them under all three costs; 173 gains and 673 losses are stated, 12 gains and 420 losses under all three costs. **Of the 673 stated losses, 19 are controlled by S5's genome sweep, 8 contradicted by it (all S5 proteome misses), and 646 are proteome-only** — the reconstruction's absences are, with few exceptions, annotation-level (D46).

## Inputs

- **Species tree** (D15): NCBI Taxonomy efetch LineageEx, retrieved 2026-10-03, archived under `<data root>/raw_api/s10/taxonomy/`; `species_tree.nwk`, `species_tree.tsv`. Root children: Amoebozoa, Apusomonadida, Discoba, Haptophyta, Metamonada, Opisthokonta, Pyrenomonadales, Rhodophyta, Sar, Viridiplantae.
- **Off the tree** (addendum (c)): Erysiphales (taxid 546991); Stephanodiscales (taxid 29204) — NCBI places each in an order already on the tree.
- **Characters** (`characters.tsv`, D50 (2)): present 12,046, absent 16,846, missing 3,883 (absence in a proteome below the BUSCO floor: 2,087, medium-only call: 1,796; BUSCO floor 70 %).
- **S4 species on the tree** for the overlay: 37 of 41 eukaryotic panel species (off: Capsaspora owczarzaki, Lottia gigantea, Lymnaea stagnalis, Torpedo marmorata).

## Costs compared

| cost | LECA present | absent | ambiguous | stated gains | stated losses | ambiguous edges |
|---|---:|---:|---:|---:|---:|---:|
| g1 | 13 | 62 | 0 | 280 | 441 | 594 |
| g2 | 19 | 54 | 2 | 173 | 673 | 451 |
| dollo | 43 | 32 | 0 | 28 | 1369 | 355 |

**Present at LECA under all three costs (13)**: ano_scramblase, cav, cng, gphr, kca_slo, kir, mscs, osca_tmem63, tmco1, tmem87, tpc, trpp, vdac.
**Present at LECA under the primary cost only (6)**: clc_channel, clic, hv1, itpr, k2p, tmc.
**Absent at LECA under all three (32)** — lineage-specific origins.

*Reading*: a LECA presence here is a profile call on proteins outside animals. Where a deep-lineage member is a family member by descent or an architecture-matched relative (plant/fungal K2P-like two-pore channels, plant CNBD channels called CNG) is a tree question this reconstruction does not ask (emergent row).

## Ancestral nodes (primary cost)

| node | present | absent | ambiguous | present under all costs |
|---|---:|---:|---:|---:|
| Eukaryota | 19 | 54 | 2 | 13 |
| Opisthokonta | 20 | 51 | 4 | 15 |
| Fungi | 10 | 63 | 2 | 9 |
| Metazoa | 31 | 41 | 3 | 24 |
| Bilateria | 49 | 24 | 2 | 48 |
| Ecdysozoa | 51 | 22 | 2 | 48 |
| Lophotrochozoa | 51 | 23 | 1 | 49 |
| Deuterostomia | 50 | 23 | 2 | 50 |
| Vertebrata | 62 | 13 | 0 | 62 |
| Mammalia | 62 | 12 | 1 | 62 |
| Viridiplantae | 19 | 55 | 1 | 13 |
| Embryophyta | 16 | 59 | 0 | 15 |
| Sar | 19 | 56 | 0 | 13 |
| Amoebozoa | 20 | 54 | 1 | 16 |
| Discoba | 15 | 57 | 3 | 9 |

## Families with the most stated losses (primary)

| family | present orders | stated losses | robust | gains | LECA (g1 / g2 / Dollo) |
|---|---:|---:|---:|---:|---|
| hv1 | 269 | 58 | 32 | 4 | absent / present / present |
| mscl | 74 | 30 | 6 | 13 | absent / absent / present |
| tmem87 | 313 | 29 | 15 | 1 | present / present / present |
| gphr | 382 | 25 | 21 | 0 | present / present / present |
| tpc | 260 | 23 | 12 | 0 | present / present / present |
| catsper | 86 | 19 | 12 | 13 | absent / absent / present |
| clc_channel | 179 | 19 | 13 | 0 | absent / present / present |
| k2p | 284 | 19 | 7 | 0 | absent / present / present |

## Losses against S5's genome cells (D50 (4))

| family | loss on | orders below | strength | S5 verdicts in the absent region |
|---|---|---:|---|---|
| clcc1 | Stylommatophora | 1 | contradicted | Cornu aspersum=genome_present |
| gphr | Diplostraca | 1 | contradicted | Daphnia pulex=genome_found |
| hcn | Diplostraca | 1 | contradicted | Daphnia pulex=genome_found |
| mitok | Batrachia | 2 | contradicted | Xenopus tropicalis=genome_found |
| ryr | Tetraodontiformes | 1 | contradicted | Takifugu rubripes=genome_found |
| tmem175 | Tetraodontiformes | 1 | contradicted | Takifugu rubripes=genome_found |
| trpn | Trichoplacea | 1 | contradicted | Trichoplax adhaerens=genome_found |
| trpn | Actiniaria | 1 | contradicted | Nematostella vectensis=genome_found |
| calhm | Mandibulata | 25 | controlled | Daphnia pulex=absent; Drosophila melanogaster=partial |
| catsper | Petromyzontiformes | 1 | controlled | Petromyzon marinus=absent |
| clcc1 | Diplostraca | 1 | controlled | Daphnia pulex=absent |
| clcc1 | Rhabditida | 1 | controlled | Caenorhabditis elegans=absent |
| hv1 | Rhabditida | 1 | controlled | Caenorhabditis elegans=absent |
| kca_slo | Magnoliopsida | 38 | controlled | Arabidopsis thaliana=absent; Oryza sativa=partial |
| mcu | Trichoplacea | 1 | controlled | Trichoplax adhaerens=absent |
| p2x | Neoptera | 10 | controlled | Drosophila melanogaster=absent |
| p2x | Rhabditida | 1 | controlled | Caenorhabditis elegans=absent |
| p2x_nonmetazoan | Trichoplacea | 1 | controlled | Trichoplax adhaerens=absent |
| pacc | Tetraodontiformes | 1 | controlled | Takifugu rubripes=absent |
| pannexin | Ascidiacea | 2 | controlled | Ciona intestinalis=absent |
| piezo | Branchiopoda | 2 | controlled | Daphnia pulex=absent |
| tmco1 | Craspedida | 1 | controlled | Monosiga brevicollis=absent |
| tmem175 | Trichoplacea | 1 | controlled | Trichoplax adhaerens=absent |
| tpc | Diptera | 1 | controlled | Drosophila melanogaster=absent |
| trpp | Branchiopoda | 2 | controlled | Daphnia pulex=absent |
| tweety | Trichoplacea | 1 | controlled | Trichoplax adhaerens=absent |
| zac | Rodentia | 1 | controlled | Mus musculus=genome_weak; Rattus norvegicus=absent |

Every `contradicted` loss is a cell S5 already reported as a proteome miss (`genome_found`) or a genome-only presence (*Cornu*): the proteome-level reconstruction states a loss the genome refutes, which is the reason D46 keeps the two strengths apart.

## Not done here

- Q4 (animal MscS) and the emergent absence checks — **S10b**, under D50 (6).
- Family-specific readings of gains (horizontal transfer, contamination: e.g. the one sponge order with a KcsA-like call) — S10b / S21.

## Files

`species_tree.{nwk,tsv,json}`, `s4_species_orders.tsv`, `characters.tsv`, `families.tsv` (per family × cost), `nodes.tsv` (every internal node × family × cost), `events.tsv` (every edge where a change is possible; `stated`, `robust`), `losses.tsv` (primary stated losses with strength), `figures/repertoire.png`.
