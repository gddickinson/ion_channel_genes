# S10b — Q4 and the absence checks

*Rendered by `scripts/s10b_report.py` from the tables beside it (D13). Rules: D51, fixed before any check was read; Q4 under D50 (6).*

## Headline

**Q4 (animal MscS): unresolved** — embedded but medium-only: ['A0A813VR16', 'A0A813VVU4', 'A0A815G4N8', 'A0A815VFJ4', 'A0AA35WXR1']; animal S5 cells not `absent`: {'gap': 13, 'genome_found': 1, 'partial': 13, 'present': 1}.

## 1. Q4 — the embedded-contig test (D51 (1)–(2))

Each subject's contig was searched by `blastx` against the S4 panel (own species removed); a region's lineage is its best hit's only past a 10 % bitscore margin. *embedded* = ≥ 1 region with a metazoan best hit.

| subject | species | order | confidence | contig | contig_len | regions | metazoan_regions | foreign_regions | unclear_regions | locus_best_lineage | locus_best_species | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Nematostella vectensis:DS475034.1:6095-6706+ | Nematostella vectensis |  | high | DS475034.1 | 6737 | 1 | 0 | 1 | 0 | prokaryote | Escherichia coli | foreign |
| Nematostella vectensis:DS477768.1:2677-3063+ | Nematostella vectensis |  | medium | DS477768.1 | 3083 | 0 | 0 | 0 | 0 | Metazoa | Strongylocentrotus purpuratus | no_evidence |
| A0A813VR16 | Adineta steineri | Adinetida | medium | CAJNOM010000026 | 682658 | 70 | 67 | 1 | 2 | Metazoa | Lymnaea stagnalis | embedded |
| A0A813VVU4 | Adineta steineri | Adinetida | medium | CAJNOM010000026 | 682658 | 55 | 49 | 1 | 5 | Metazoa | Lymnaea stagnalis | embedded |
| A0A815G4N8 | Adineta steineri | Adinetida | medium | CAJNOM010000299 | 192564 | 31 | 27 | 1 | 3 |  |  | embedded |
| A0A815VFJ4 | Rotaria sordida | Philodinida | medium | CAJNOL010002875 | 30252 | 2 | 2 | 0 | 0 |  |  | embedded |
| A0AA35R2R3 | Geodia barretti | Tetractinellida | high | CASHTH010000421 | 101297 | 7 | 0 | 2 | 5 |  |  | foreign |
| A0AA35WXR1 | Geodia barretti | Tetractinellida | medium | CASHTH010002892 | 55697 | 5 | 1 | 2 | 2 |  |  | embedded |
| B3SFM1 | Trichoplax adhaerens | Trichoplacea | high | DS986442.1 | 1151 | 0 | 0 | 0 | 0 | prokaryote | Escherichia coli | no_evidence |

*Reading*: no high-confidence animal MscS sits on a contig with an animal gene — the *Nematostella* locus and the *Geodia* (sponge) protein are on contigs whose other genes are non-animal, the *Trichoplax* protein on a 1.2 kb contig with nothing else. The contigs that are animal carry medium-only calls (bdelloid rotifers, a second *Geodia* protein), which D50 (2) never counts as a presence. So the presence arm fails; the absence arm fails too, because S5 never returned `absent` for an animal MscS cell (bait alignments the profiles do not call: `partial` / `gap`).

S5's animal MscS cells (the absence arm): gap 13, present 1, genome_found 1, partial 13.

### Prokaryotic-family calls in animal proteomes (descriptive, S21)

| family | embedded | foreign | no_evidence | failed |
|---|---|---|---|---|
| iglur_prok | 6 |  |  |  |
| kcsa_prok |  | 1 |  |  |
| mscl | 13 | 2 | 7 |  |

| subject | family | species | order | contig_len | regions | metazoan_regions | foreign_regions | locus_best_lineage | locus_best_species | status |
|---|---|---|---|---|---|---|---|---|---|---|
| A0A813YWE9 | mscl | Rotaria sordida | Philodinida | 7246 | 0 | 0 | 0 | prokaryote | Escherichia coli | no_evidence |
| A0A813YXJ7 | mscl | Rotaria sordida | Philodinida | 7246 | 0 | 0 | 0 | prokaryote | Escherichia coli | no_evidence |
| A0A814EJH1 | mscl | Adineta steineri | Adinetida | 173460 | 2 | 0 | 0 |  |  | no_evidence |
| A0A814EWP7 | mscl | Adineta steineri | Adinetida | 173460 | 2 | 0 | 0 |  |  | no_evidence |
| A0A814GSX9 | mscl | Rotaria sordida | Philodinida | 64040 | 11 | 10 | 0 |  |  | embedded |
| A0A814HMK0 | mscl | Adineta steineri | Adinetida | 246426 | 46 | 44 | 1 |  |  | embedded |
| A0A814HP98 | mscl | Adineta steineri | Adinetida | 246426 | 46 | 44 | 1 |  |  | embedded |
| A0A814N963 | mscl | Adineta steineri | Adinetida | 98453 | 8 | 2 | 5 |  |  | embedded |
| A0A814PY56 | mscl | Rotaria sordida | Philodinida | 62909 | 1 | 1 | 0 | prokaryote | Escherichia coli | embedded |
| A0A814QPV4 | mscl | Adineta steineri | Adinetida | 363664 | 17 | 10 | 3 |  |  | embedded |
| A0A814U258 | mscl | Rotaria sordida | Philodinida | 53720 | 12 | 9 | 3 | prokaryote | Escherichia coli | embedded |
| A0A815AM50 | mscl | Rotaria sordida | Philodinida | 29146 | 1 | 0 | 0 |  |  | no_evidence |
| A0A815BEV6 | mscl | Adineta steineri | Adinetida | 56840 | 1 | 0 | 0 |  |  | no_evidence |
| A0A815G1N3 | mscl | Adineta steineri | Adinetida | 192624 | 45 | 42 | 0 |  |  | embedded |
| A0A815G8S3 | mscl | Rotaria sordida | Philodinida | 69120 | 12 | 10 | 0 | prokaryote | Escherichia coli | embedded |
| A0A815I888 | mscl | Rotaria sordida | Philodinida | 64040 | 11 | 10 | 0 |  |  | embedded |
| A0A815K7G0 | mscl | Rotaria sordida | Philodinida | 53712 | 13 | 10 | 3 | prokaryote | Escherichia coli | embedded |
| A0A815LSB8 | mscl | Adineta steineri | Adinetida | 152700 | 27 | 24 | 0 |  |  | embedded |
| A0A815X1P8 | mscl | Adineta steineri | Adinetida | 58399 | 19 | 9 | 8 |  |  | embedded |
| A0A815X5T0 | mscl | Adineta steineri | Adinetida | 56840 | 1 | 0 | 0 |  |  | no_evidence |
| A0A816E4M0 | mscl | Adineta steineri | Adinetida | 8218 | 1 | 0 | 1 | prokaryote | Escherichia coli | foreign |
| A0AA35QW52 | kcsa_prok | Geodia barretti | Tetractinellida | 63081 | 18 | 0 | 15 | prokaryote | Gloeobacter violaceus | foreign |
| A0AA88H0Y7 | mscl | Artemia franciscana | Anostraca | 24125 | 11 | 0 | 9 | prokaryote | Aliarcobacter butzleri | foreign |
| A0ABN8NQA6 | iglur_prok | Porites lobata | Scleractinia | 3890156 | 6 | 6 | 0 |  |  | embedded |
| A0ABN8NSW9 | iglur_prok | Porites lobata | Scleractinia | 3872037 | 4 | 4 | 0 |  |  | embedded |
| A0ABN8P8Z2 | iglur_prok | Porites lobata | Scleractinia | 3077898 | 10 | 10 | 0 |  |  | embedded |
| A0ABN8PB31 | iglur_prok | Porites lobata | Scleractinia | 3077898 | 10 | 10 | 0 |  |  | embedded |
| A0ABN8RQ03 | iglur_prok | Porites lobata | Scleractinia | 442248 | 5 | 5 | 0 |  |  | embedded |
| A0ABN8RSD9 | iglur_prok | Porites lobata | Scleractinia | 442248 | 6 | 6 | 0 | Metazoa | Danio rerio | embedded |

## 2. *Daphnia pulex* on its current assembly (D51 (3))

Assembly GCA_000187875.1 (2011, scaffolds) → **GCF_021134715.1** (chosen by S4's `assembly_rank`; N50 12,288,052 bp). Matched detection 1.0 → **1.0** (40 matched control cells), so the new verdicts replace the old for S10's reading. 7 cells change.

| family | old_verdict | new_verdict | new_best_confidence | new_found_loci |
|---|---|---|---|---|
| kv_modifier | gap | partial |  |  |
| kcsa_prok | no_locus_unrescued | partial |  |  |
| catsper | partial | no_locus_unrescued |  |  |
| hcn | genome_found | genome_weak | medium | NC_060019.1:7491016-7576124+ |
| trpp | absent | partial |  |  |
| plgic_prok | absent | absent |  |  |
| deg_invertebrate | absent | absent |  |  |
| piezo | absent | absent |  |  |
| calhm | absent | partial |  |  |
| tmem175 | partial | no_locus_unrescued |  |  |
| clcc1 | absent | absent |  |  |

*Reading*: of the 2011 assembly's absences, plgic_prok, deg_invertebrate, piezo, clcc1 stay `absent` on a chromosome-level genome with every matched control detected; trpp, calhm become `partial` (a single weak bait alignment no profile scores), so they are no longer controlled absences. The two S10a contradictions in *Daphnia* (HCN, GPHR) stay genome presences (HCN now medium).

## 3. ZAC, PACC1, CLCC1 against two ortholog databases (D51 (4))

Every S5 `absent` cell of the three families; NCBI Gene orthologs and Ensembl Compara queried with the human gene.

| family | species | ncbi | ncbi_ids | ensembl | ensembl_ids | status |
|---|---|---|---|---|---|---|
| pacc | Trichoplax adhaerens | none |  | not_in_ensembl |  | not_covered |
| clcc1 | Caenorhabditis elegans | none |  | not_in_ensembl |  | not_covered |
| clcc1 | Drosophila melanogaster | none |  | not_in_ensembl |  | not_covered |
| clcc1 | Daphnia pulex | none |  | not_in_ensembl |  | not_covered |
| zac | Callorhinchus milii | none |  | none |  | agrees |
| zac | Torpedo marmorata | none |  | not_in_ensembl |  | agrees |
| zac | Danio rerio | none |  | listed | ENSDARG00000111812,ENSDARG00000115573,ENSDARG00000115752 | disputed |
| zac | Takifugu rubripes | none |  | listed | ENSTRUG00000027812 | disputed |
| pacc | Takifugu rubripes | listed | 101061834 | none |  | disputed |
| zac | Latimeria chalumnae | none |  | none |  | agrees |
| zac | Xenopus tropicalis | none |  | none |  | agrees |
| zac | Gallus gallus | none |  | none |  | agrees |
| zac | Ornithorhynchus anatinus | none |  | none |  | agrees |
| zac | Rattus norvegicus | none |  | none |  | agrees |

Every protein a database lists where S5 found none, scored against the S3a profiles (D32) and placed on S5's genome by `tblastn`:

| family | species | protein | length | p_call | p_family | p_confidence | win_score | runner | genome_hit | genome_pident | s5_locus |
|---|---|---|---|---|---|---|---|---|---|---|---|
| zac | Danio rerio | ensembl:ENSDARP00000153861 | 469 | family | ht3 | medium | 139.2 | nachr | NC_133178.1:60936301-60936849 | 98.378 | none |
| zac | Danio rerio | ensembl:ENSDARP00000152326 | 442 | family | ht3 | medium | 119.5 | nachr | NC_133178.1:60746520-60746816 | 100.0 | none |
| zac | Danio rerio | ensembl:ENSDARP00000149349 | 418 | family | ht3 | medium | 124.1 | nachr | NC_133178.1:60862386-60862958 | 77.551 | none |
| zac | Takifugu rubripes | ensembl:ENSTRUP00000060956 | 421 | superfamily_only |  | low | 114.4 | nachr | LR584247.1:8914221-8914694 | 86.076 | none |
| pacc | Takifugu rubripes | ncbi:XP_011618135.2 | 353 | family | pacc | high | 480.1 |  | LR584245.2:14721205-14721870 | 73.423 | r4_LR584245.2:14636898-14754054+=family:tmem87:high |

*Reading*: the seven ZAC absences no database disputes stand as measured; the two fish ZAC disputes are proteins Ensembl's gene trees call ZACN orthologues and our profiles call 5-HT3 or Cys-loop superfamily-only — a disagreement between orthology inference and profile assignment that neither resolves. The *Takifugu* PACC1 dispute is ours: the NCBI protein is a high-confidence PACC1 by our own profile and lies inside an S5 locus whose majority bait was TMEM87, so S5's clustering absorbed it and rescue (which ignores hits inside a locus) never looked. CLCC1 in ecdysozoans: neither database covers it (Ensembl 116's main site carries no fly or worm; NCBI's CLCC1 set names no invertebrate of the panel).

### Post-hoc screen: absences under a runner-up bait

Found by the PACC1 check, so labelled post hoc: every S5 `absent` cell whose family is the runner-up bait of a locus in that genome. 10 loci; 2 across superfamilies. Within a superfamily this is ordinary bait overlap (Nav under Cav, ENaC under ASIC); across superfamilies it is the PACC1 artefact and a no-call rat locus.

| species | family | locus | span | bait_family | p_family | p_confidence | same_superfamily |
|---|---|---|---|---|---|---|---|
| Amphimedon queenslandica | nav | GL345140.1:147500-196583+ | 49084 | cav | cav | high | 1 |
| Caenorhabditis elegans | nav | BX284604.4:7407965-7497639+ | 89675 | cav | cav | high | 1 |
| Caenorhabditis elegans | nav | BX284606.5:2716921-2806415- | 89495 | cav | cav | high | 1 |
| Caenorhabditis elegans | nav | BX284606.5:7849573-7857807- | 8235 | cav | cav | high | 1 |
| Petromyzon marinus | catsper | NC_133688.1:20464494-21259300- | 794807 | cav | cav | high | 1 |
| Danio rerio | enac | NC_133184.1:10971552-11348021+ | 376470 | asic | asic | high | 1 |
| Takifugu rubripes | zac | LR584240.2:8339978-8348161+ | 8184 | ht3 | ht3 | high | 1 |
| Takifugu rubripes | enac | LR584235.1:6710313-6901041- | 190729 | asic | asic | high | 1 |
| Takifugu rubripes | pacc | r4_LR584245.2:14636898-14754054+ | 117157 | tmem87 | tmem87 | high | 0 |
| Rattus norvegicus | zac | CM070411.1:63140811-63775205- | 634395 | nalcn |  | none | 0 |

## 4. What the K2P profile calls outside animals (D51 (5))

12 non-metazoan high-confidence K2P calls on the S4 panel:

| species | accession | protein_name | pfam | tm_features | win_score | runner |
|---|---|---|---|---|---|---|
| Monosiga brevicollis | A9VAE4 | Potassium channel domain-containing protein | PF00406x1;PF07885x2 | 6 | 158.4 | kcsa_prok |
| Arabidopsis thaliana | Q8LBL1 | Two-pore potassium channel 1 | PF07885x2 | 4 | 219.7 | kcsa_prok |
| Arabidopsis thaliana | Q9FL25 | Two-pore potassium channel 2 | PF07885x2 | 4 | 150.3 | kcsa_prok |
| Arabidopsis thaliana | Q9S6Z8 | Two-pore potassium channel 5 | PF07885x2 | 4 | 162.4 | kcsa_prok |
| Arabidopsis thaliana | Q9SVV6 | Two-pore potassium channel 3 | PF07885x2 | 4 | 148.6 | kcsa_prok |
| Oryza sativa | A0A0P0V6X9 | Os01g0696100 protein | PF07885x2 | 4 | 130.0 | kcsa_prok |
| Oryza sativa | Q69TN4 | Two pore potassium channel c | PF07885x2 | 4 | 222.4 | kcsa_prok |
| Oryza sativa | Q850M0 | Two pore potassium channel a | PF07885x2 | 4 | 150.9 | kcsa_prok |
| Oryza sativa | Q8LIN5 | Two pore potassium channel b | PF07885x2 | 4 | 158.9 | kcsa_prok |
| Physcomitrium patens | A0A2K1ITF6 | Potassium channel domain-containing protein | PF07885x2 | 6 | 163.5 | kcsa_prok |
| Physcomitrium patens | A0A2K1IUA0 | Potassium channel domain-containing protein | PF07885x2 | 5 | 152.4 | kcsa_prok |
| Physcomitrium patens | A0A2K1K7C7 | Uncharacterized protein | PF07885x2;PF13202x2 | 5 | 161.3 | kcsa_prok |

S4b orders with a high-confidence K2P call, by kingdom: protists 14, Fungi 30, Metazoa 192, Viridiplantae 48.

*Reading*: on the S4 panel the K2P profile's non-animal calls are the plant TPK two-pore K⁺ channels (every one carries two `PF07885` pore domains and four TM helices) and one *Monosiga* protein; yeast TOK1 is not called. The K2P LECA presence of S10a therefore rests on TPK-type (and, in S4b, fungal and protist) two-pore channels — whether those are K2P by descent is S17's tree question.

## Files

`contig_subjects.tsv`, `contig_test.tsv`, `contig_regions.tsv`, `q4.json`; `daphnia_assembly.json`, `daphnia_cells.tsv`, `daphnia_compare.tsv`, `daphnia_summary.json`; `ortholog_check.tsv`, `ortholog_disputed.tsv`; `absent_runner_screen.tsv`; `k2p_nonmetazoan.tsv`, `k2p_by_kingdom.tsv`; `s10b_summary.json`; `figures/s10b_checks.png`. Raw fetches under `<data root>/raw_api/s10b/`.
