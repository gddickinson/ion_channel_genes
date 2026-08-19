# S1 — classifier control benchmark

*Generated 2026-08-19 16:30:28 by `scripts/s1_report.py` from the tables in `results/benchmark_controls` (D13). Panel: `controls_benchmark`. 97/97 targets classified in 866.0 s.*

## Headline

| measure | value |
|---|---|
| recall (positives called correctly) | **50/72** (69.4%) |
| specificity (decoys not called as channels) | **25/25** (100.0%) |
| hazards exercised | 16/16 |
| four-repeat anchor validates | yes |
| MAFFT | v7.526 (2024/Apr/26) |
| reference panel | 162 exemplar sequences |
| panel proteins that are themselves exemplars | 93 — classified leave-one-out |

**Leave-one-out matters here.** 93 of the 97 panel proteins are themselves catalogue exemplars — both are drawn from the same curated gene lists — so without excluding the query's own accession the reference tier would score almost the whole panel against itself at 100 % identity and the benchmark would measure nothing but that overlap (**D29**).

## Which tier made the call

This is the number that says whether the classifier is doing anything a nearest-neighbour lookup could not. A call made by the architecture or hazard tier rests on domain composition and works on a sequence with no close relative in the panel; a call made only by the reference tier does not.

| decisive tier | n calls |
|---|---|
| none | 37 |
| reference | 24 |
| hazard | 17 |
| architecture | 12 |
| motif | 7 |

| confidence | n |
|---|---|
| bronze | 42 |
| gold | 25 |
| unassigned | 17 |
| silver | 13 |

## Per-family recall

37/57 **positive** families were called correctly for every member. The 18 decoy families are scored by specificity, not recall — a decoy left `unassigned` is a success, and listing it as a recall failure would invert the result.

**Positive families the classifier did not call correctly:**

| family | n | correct | recall | called instead |
|---|---|---|---|---|
| ano_scramblase | 1 | 0 | 0.000 | ano_channel |
| asic | 1 | 0 | 0.000 | unassigned |
| catsper | 1 | 0 | 0.000 | unassigned |
| clc_transporter | 1 | 0 | 0.000 | unassigned |
| delta_glur | 1 | 0 | 0.000 | unassigned |
| enac | 1 | 0 | 0.000 | unassigned |
| gabaa | 2 | 0 | 0.000 | unassigned |
| glyr | 1 | 0 | 0.000 | unassigned |
| ht3 | 1 | 0 | 0.000 | unassigned |
| kainate | 1 | 0 | 0.000 | unassigned |
| kv_modifier | 1 | 0 | 0.000 | unassigned |
| kv_shaker | 2 | 1 | 0.500 | kv_modifier |
| nachr | 2 | 1 | 0.500 | unassigned |
| nmda | 1 | 0 | 0.000 | unassigned |
| p2x | 2 | 1 | 0.500 | unassigned |
| pannexin | 1 | 0 | 0.000 | unassigned |
| tpc | 1 | 0 | 0.000 | unassigned |
| trpa | 1 | 0 | 0.000 | unassigned |
| trpv | 2 | 0 | 0.000 | unassigned |
| zac | 1 | 0 | 0.000 | nonchannel_achbp |

**Decoys:** 18 families; 18 were left unassigned or called as themselves, 0 were called as something else.

Full table: `recall.tsv`.

## Per-hazard

A hazard is only closed by a test. `n touched` counts panel members whose classification involved that hazard's rule; a hazard with zero is not solved, it is **untested** by this panel.

| hazard | severity | n touched | correct | wrong | examples | what it is |
|---|---|---|---|---|---|---|
| H1 | high | 7 | 7 | 0 | — | Nav, Cav, NALCN and CatSper share one domain architectur |
| H2 | high | 7 | 1 | 6 | P36544;P14867;P28472;P23415;P46098;Q401N2 | The Cys-loop ligand-binding domain exists without a chan |
| H3 | high | 6 | 3 | 3 | P39086;Q05586;O43424 | The iGluR clamshell is the class C GPCR ligand-binding d |
| H4 | high | 2 | 2 | 0 | — | IP3 receptors and ryanodine receptors share every diagno |
| H5 | medium | 1 | 1 | 0 | — | Half the CLC family are antiporters, not channels |
| H6 | medium | 2 | 1 | 1 | Q4KMQ2 | Anoctamin channels and scramblases are architecturally i |
| H7 | high | 5 | 5 | 0 | — | Most TRP families carry no `PF00520` at all |
| H8 | medium | 7 | 7 | 0 | — | Pfam's pore model does not track the number of transmemb |
| H9 | high | 2 | 2 | 0 | — | A voltage-sensor domain is not evidence of a channel |
| H10 | low | 3 | 3 | 0 | — | The MIR domain is shared with the O-mannosyltransferases |
| H11 | medium | 2 | 2 | 0 | — | CFTR and the sulfonylurea receptors are the same ABC arc |
| H12 | medium | 3 | 2 | 1 | Q14721 | The Kv tetramerisation domain is a generic BTB/POZ domai |
| H13 | medium | 3 | 2 | 1 | P98161 | `PF08016` covers TRPML, TRPP *and* polycystin-1 |
| H14 | high | 2 | 2 | 0 | — | Ion selectivity is not recoverable from sequence family |
| H15 | high | 1 | 1 | 0 | — | Gene-symbol prefixes group pores with their accessory su |
| H16 | low | 4 | 4 | 0 | — | Kir and K2P look like one 2TM/1P group and are modelled  |

## Panel coverage

41/68 census families have at least one member in this panel. The families with none are mostly the non-human ones — prokaryotic channels, invertebrate degenerins, plant OSCA — and a benchmark that only measures human proteins cannot claim anything about a census that is not human-only.

Families with no panel member: `kcsa_prok`, `catsper`, `tpc`, `trpv`, `trpa`, `trpn`, `gabaa`, `glyr`, `ht3`, `zac`, `plgic_invertebrate`, `plgic_prok`, `kainate`, `nmda`, `delta_glur`, `iglur_nonvertebrate`, `iglur_prok`, `p2x_nonmetazoan`, `enac`, `asic`, `deg_invertebrate`, `ano_scramblase`, `mscl`, `mscs`, `pannexin`, `innexin`, `viroporin`

## Every call

`calls.tsv` carries one row per protein with the full evidence string from every tier, including the tiers that did not decide. The point of keeping the losing evidence is that a wrong call is diagnosable without re-running anything.

**2 call(s) had disagreeing tiers:**

| protein | gene | expected | called | confidence | conflict |
|---|---|---|---|---|---|
| Q14721 | KCNB1 | kv_shaker | kv_modifier | silver | architecture=kv_shaker;reference=kv_modifier |
| P35523 | CLCN1 | clc_channel | clc_channel | bronze | topology: 5 TM observed vs 18 expected (±4) for clc_channel |

