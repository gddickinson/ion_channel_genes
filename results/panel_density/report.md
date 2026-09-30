# Panel density for S10 — the options, measured

Rendered by `scripts/s10_panel_density.py` from `options.tsv` (D13).

UniProt holds **3,521 eukaryotic reference proteomes**. A denser S10 panel takes one per taxonomic rank, chosen by the best BUSCO completeness (then most proteins) — S4's rule, restated for many species. The census sweep scales with sequences; the S3b-rate estimate below is for the profile sweep alone.

| option | proteomes | sequences | profile sweep (h, S3b rate) | kingdoms |
|---|---|---|---|---|
| current panel (S4; eukaryotic proteomes) | 39 | — | — |  |
| one per phylum | 53 | 1,107,971 | 0.3 | Metazoa 18; (protists) 17; Fungi 16; Viridiplantae 2 |
| one per class | 141 | 2,730,111 | 0.7 | Metazoa 48; Fungi 45; (protists) 32; Viridiplantae 16 |
| one per order | 439 | 10,949,076 | 2.8 | Metazoa 194; Fungi 126; (protists) 62; Viridiplantae 57 |
| one per family | 1093 | 25,189,186 | 6.3 | Metazoa 591; Fungi 313; Viridiplantae 112; (protists) 77 |

**What does not scale is the genome check.** S5's D4 absence bars (matched in-group bait, measured detection, contiguity) needed each genome downloaded and swept: 42.5 Gbp for 52 genomes. At one proteome per order that is ~400 genomes and hundreds of GB. A denser panel's absences are therefore **proteome absences** — annotation-level, weaker than the 52-species panel's controlled absences, which stay the high-evidence subset.
