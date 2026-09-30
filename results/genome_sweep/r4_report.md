# S5 r4 — the genome sweep for the families added after S20 (D43)

Rendered by `scripts/s5r4_compare.py` from `baits_r4.tsv`, `genome_runs_r4.tsv`, `r4_cells.tsv`, `r4_cell_changes.tsv` and `r4_control_changes.tsv` (D13).

![S5 r4](figures/genome_r4.png)

**Baits**: 287 for 12 families (the 7 new census families, KChIP and the 4 decoys), drawn by S5b's rules B1–B3 restricted to them, in a panel of their own so S5b's panel, runs and loci stay byte-identical. **Runs**: 52/52 genomes, 1,267 loci, 774 called to an r4 family, miniprot 50 min. An r4 locus counts only toward an r4 family (called to one, or partial on its bait); an old family's trace is judged against S5b's loci only.

**Old cells: 3,536 compared, 0 changed, 0 verdicts changed.** 41 genomes' control counts grew (the new families' high-confidence cells join the controls); none crossed the 0.90 floor.

## The new census families, per species

| verdict | clcc1 | gphr | mitok | pacc | tmco1 | tmem109 | tmem87 |
|---|---|---|---|---|---|---|---|
| absent | 3 | 1 | 2 | 2 | 1 | 2 | 1 |
| gap | 0 | 0 | 1 | 2 | 0 | 2 | 1 |
| genome_found | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
| genome_present | 2 | 2 | 2 | 1 | 2 | 0 | 2 |
| genome_weak | 0 | 1 | 0 | 0 | 0 | 1 | 0 |
| no_locus_unrescued | 23 | 11 | 21 | 30 | 11 | 28 | 15 |
| partial | 3 | 1 | 4 | 4 | 3 | 9 | 0 |
| present | 20 | 35 | 21 | 13 | 35 | 10 | 33 |
| trace | 1 | 0 | 0 | 0 | 0 | 0 | 0 |

**Proteome misses (high-confidence intact loci)**: *Daphnia pulex* gphr; *Xenopus tropicalis* mitok.

**Controlled absences (matched bait, detection ≥ 0.90, D4 bar met)**: *Saccharomyces cerevisiae* tmem87 (bar pooled); *Schizosaccharomyces pombe* gphr (bar pooled); *Monosiga brevicollis* tmco1 (bar pooled); *Trichoplax adhaerens* pacc (bar pooled); *Caenorhabditis elegans* clcc1 (bar pooled); *Drosophila melanogaster* clcc1 (bar pooled); *Drosophila melanogaster* mitok (bar pooled); *Daphnia pulex* clcc1 (bar pooled); *Lottia gigantea* mitok (bar pooled); *Petromyzon marinus* tmem109 (bar group); *Callorhinchus milii* tmem109 (bar group); *Takifugu rubripes* pacc (bar group).

**Weak genome evidence**: *Chlamydomonas reinhardtii* tmem109; *Ciona intestinalis* gphr.

**Gaps**: *Chlamydomonas reinhardtii* tmem87; *Monosiga brevicollis* mitok; *Lottia gigantea* tmem109; *Cornu aspersum* pacc; *Cornu aspersum* tmem109; *Strongylocentrotus purpuratus* pacc.

These are measurements under the S5b instrument, not literature statements: each absence is a candidate for the S10 checks, and the single-seed profiles (TMCO1, TMEM87, TMEM109, CLCC1, MITOK) call loci at the same D32 gates as every other family.
