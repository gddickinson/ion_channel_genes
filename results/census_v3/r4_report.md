# Census v2/v3a revision r4 — the families added after S20 (D43)

Rendered by `scripts/s3r4_compare.py` from `r4_transitions.tsv`, `r4_new_records.tsv`, `r4_reviewed_calls.tsv` and `benchmark_summary.json` (D13).

**26,783 records added** to census v2 (delta walk, `s2r4_delta.py`), swept by all profiles; **12 profiles added** (7 new census families; KChIP; the EMC3, GOST, BRI3BP and neuronal-calcium-sensor decoys), each swept over r3's NR database with `-Z` held at r3's size, so no r3 hit moved. **11 of 1,245,200 r3 calls changed**, every one with a new family as the new call, the winner or the runner-up that pulled the winner inside the D7 margin (checked, hard failure otherwise):

| before | after | basis | profile winner | records |
|---|---|---|---|---|
| unassigned | [ploop] | superfamily_only | — | 3 |
| trpml | conflict | conflict | clcc1 | 1 |
| connexin | gphr | both | gphr | 1 |
| iglur_nonvertebrate | conflict | conflict | nonchannel_gost | 1 |
| nachr | [cysloop] | superfamily_only | — | 1 |
| asic | conflict | conflict | mitok | 1 |
| unassigned | nonchannel_ncs | profile_only | nonchannel_ncs | 1 |
| osca_tmem63 | nonchannel_ncs | profile_only | nonchannel_ncs | 1 |
| unassigned | gphr | both | gphr | 1 |

## The new records' calls

| v3 call | records |
|---|---|
| tmem87 | 5,484 |
| nonchannel_gost | 4,441 |
| nonchannel_emc3 | 4,280 |
| gphr | 3,903 |
| tmco1 | 2,237 |
| unassigned | 1,916 |
| clcc1 | 1,562 |
| pacc | 992 |
| nonchannel_bri3bp | 804 |
| tmem109 | 763 |
| [tmem87] | 398 |
| [tmem109] | 3 |

## The look-alike tests (H17–H20)

S3a's held-out orthologue benchmark (set B) now scores 542 channel orthologues (542 correct) and 57 decoys (0 called a channel). Set B skips seeds and 'Precursor'-flagged records, so every reviewed r4 record is listed here with its call (the gene symbol scores, never makes, the call — H15):

| gene | organism | seed | profile call | confidence | margin | runner-up |
|---|---|---|---|---|---|---|
| — | Methanocaldococcus jannaschi |  | — | none |  | — |
| — | Schizosaccharomyces pombe (s |  | — | none |  | — |
| — | Ostreid herpesvirus 1 (isola |  | clcc1 | medium | 1.0 | — |
| CLCC1 | Homo sapiens | yes | clcc1 | high | 1.0 | — |
| CLCC1 | Bos taurus |  | clcc1 | high | 1.0 | — |
| Clcc1 | Rattus norvegicus |  | clcc1 | high | 1.0 | — |
| Clcc1 | Mus musculus |  | clcc1 | high | 1.0 | — |
| clcc1 | Danio rerio |  | clcc1 | high | 1.0 | — |
| clcc1 | Xenopus laevis |  | clcc1 | high | 1.0 | — |
| — | Saccharomyces cerevisiae (st | yes | gphr | high | 1.0 | — |
| COLD1 | Oryza sativa subsp. japonica | yes | gphr | high | 1.0 | — |
| COLD1 | Oryza sativa subsp. indica | yes | gphr | high | 1.0 | — |
| GPHR | Bos taurus | yes | gphr | high | 1.0 | — |
| GPHR | Cricetulus griseus | yes | gphr | high | 1.0 | — |
| GPHRA | Homo sapiens | yes | gphr | high | 1.0 | — |
| GPHRB | Homo sapiens | yes | gphr | high | 1.0 | — |
| GPR89 | Gallus gallus | yes | gphr | high | 1.0 | — |
| GTG1 | Arabidopsis thaliana |  | gphr | high | 1.0 | — |
| GTG2 | Arabidopsis thaliana | yes | gphr | high | 1.0 | — |
| Gphr | Mus musculus | yes | gphr | high | 1.0 | — |
| gphr-L | Xenopus laevis | yes | gphr | high | 1.0 | — |
| gpr89 | Salmo salar | yes | gphr | high | 1.0 | — |
| gpr89 | Dictyostelium discoideum | yes | gphr | high | 1.0 | — |
| BRI3BP | Homo sapiens | yes | nonchannel_bri3bp | high | 0.95 | tmem109 |
| Bri3bp | Mus musculus |  | nonchannel_bri3bp | high | 0.9405 | tmem109 |
| — | Schizosaccharomyces pombe (s |  | nonchannel_emc3 | high | 1.0 | — |
| AIM27 | Saccharomyces cerevisiae (st |  | nonchannel_emc3 | high | 1.0 | — |
| AIM27 | Saccharomyces cerevisiae (st |  | nonchannel_emc3 | high | 1.0 | — |
| AIM27 | Saccharomyces cerevisiae (st |  | nonchannel_emc3 | high | 1.0 | — |
| EMC3 | Homo sapiens | yes | nonchannel_emc3 | high | 1.0 | — |
| EMC3 | Bos taurus |  | nonchannel_emc3 | high | 1.0 | — |
| EMC3 | Pongo abelii |  | nonchannel_emc3 | high | 1.0 | — |
| EMC3 | Saccharomyces cerevisiae (st |  | nonchannel_emc3 | high | 1.0 | — |
| Emc3 | Rattus norvegicus |  | nonchannel_emc3 | high | 1.0 | — |
| Emc3 | Mus musculus |  | nonchannel_emc3 | high | 1.0 | — |
| emc3 | Danio rerio |  | nonchannel_emc3 | high | 1.0 | — |
| emc3 | Dictyostelium discoideum |  | nonchannel_emc3 | high | 1.0 | — |
| CAND6 | Arabidopsis thaliana |  | nonchannel_gost | high | 1.0 | — |
| CAND7 | Arabidopsis thaliana |  | nonchannel_gost | high | 1.0 | — |
| GPR107 | Homo sapiens | yes | nonchannel_gost | high | 1.0 | — |
| GPR108 | Homo sapiens | yes | nonchannel_gost | high | 1.0 | — |
| GPR108 | Bos taurus |  | nonchannel_gost | high | 1.0 | — |
| Gpr107 | Rattus norvegicus |  | nonchannel_gost | high | 0.9697 | tmem87 |
| Gpr107 | Mus musculus |  | nonchannel_gost | high | 1.0 | — |
| Gpr108 | Mus musculus |  | nonchannel_gost | high | 1.0 | — |
| Gpr108 | Rattus norvegicus |  | nonchannel_gost | high | 1.0 | — |
| TMEM87B | Homo sapiens | yes | nonchannel_gost | medium | 0.2654 | tmem87 |
| Tmem87b | Mus musculus |  | nonchannel_gost | medium | 0.1176 | tmem87 |
| PACC1 | Homo sapiens | yes | pacc | high | 1.0 | — |
| PACC1 | Bos taurus | yes | pacc | high | 1.0 | — |
| PACC1 | Pongo abelii | yes | pacc | high | 1.0 | — |
| Pacc1 | Mus musculus | yes | pacc | high | 1.0 | — |
| Pacc1 | Rattus norvegicus | yes | pacc | high | 1.0 | — |
| pacc1 | Danio rerio | yes | pacc | high | 1.0 | — |
| pacc1 | Xenopus tropicalis | yes | pacc | high | 1.0 | — |
| pacc1 | Xenopus laevis | yes | pacc | high | 1.0 | — |
| — | Caenorhabditis elegans |  | tmco1 | high | 1.0 | — |
| TMCO1 | Homo sapiens | yes | tmco1 | high | 1.0 | — |
| TMCO1 | Canis lupus familiaris |  | tmco1 | high | 1.0 | — |
| TMCO1 | Bos taurus |  | tmco1 | high | 1.0 | — |
| TMCO1 | Sus scrofa |  | tmco1 | high | 1.0 | — |
| TMCO1 | Pongo abelii |  | tmco1 | high | 1.0 | — |
| Tmco1 | Mus musculus |  | tmco1 | high | 1.0 | — |
| Tmco1 | Rattus norvegicus |  | tmco1 | high | 1.0 | — |
| tmco1 | Danio rerio |  | tmco1 | high | 1.0 | — |
| tmco1 | Dictyostelium discoideum |  | tmco1 | high | 1.0 | — |
| TMEM109 | Homo sapiens | yes | tmem109 | high | 1.0 | — |
| TMEM109 | Oryctolagus cuniculus |  | tmem109 | high | 1.0 | — |
| Tmem109 | Mus musculus |  | tmem109 | high | 1.0 | — |
| Tmem109 | Rattus norvegicus |  | tmem109 | high | 1.0 | — |
| — | Saccharomyces cerevisiae (st |  | tmem87 | medium | 1.0 | — |
| — | Schizosaccharomyces pombe (s |  | tmem87 | high | 0.3408 | nonchannel_gost |
| PTM1 | Saccharomyces cerevisiae (st |  | tmem87 | medium | 0.5967 | nonchannel_gost |
| PTM1 | Saccharomyces cerevisiae (st |  | tmem87 | medium | 0.5967 | nonchannel_gost |
| TMEM87A | Homo sapiens | yes | tmem87 | high | 0.7904 | nonchannel_gost |
| Tmem87a | Mus musculus |  | tmem87 | high | 0.7758 | nonchannel_gost |
| tmem87a | Xenopus tropicalis |  | tmem87 | high | 0.3223 | nonchannel_gost |
