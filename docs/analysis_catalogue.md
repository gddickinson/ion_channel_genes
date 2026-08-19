# What the harvested data can answer

The analysis menu behind ledger rows S15–S22, each with its inputs, method,
deliverable and the caveat that would invalidate it. Written before the data
exists, so that the analyses are chosen by what the questions need rather
than by what the tables happen to contain.

| # | Analysis | Inputs | Method | Deliverable | Invalidated by |
|---|----------|--------|--------|-------------|----------------|
| A1 | **Method contribution** | census v2/v3/v4 | recall of the S1 panel as each method is added | a recall curve per superfamily, and the answer to Q3 | a panel that does not span the superfamilies — the S1 coverage table is the check |
| A2 | **Annotation quality** | census v4, RefSeq/Ensembl/UniProt records | per-locus status: complete / split / fragmentary / unnamed / wrong family | the correction list, and a per-database quality table | comparing gene sets of different provenance (D9) |
| A3 | **Filter-vs-tree congruence** | S9 filter atlas, S8 tier-2 trees | map the filter locus onto the pore-module tree | whether selectivity tracks phylogeny in the P-loop superfamily (Q5) | filter calls made without `verify_anchor()` passing |
| A4 | **Mechanism vs clade** | S7 tier-1 trees for CLC and TMEM16 | do the channel and transporter/scramblase members form clades? | Q7 answered, and D24 either vindicated or overturned | trees rooted on the wrong outgroup |
| A5 | **Convergence audit** | S8 network, S12 structures, profile-profile | HHsearch + TM-align between families the catalogue calls unrelated | Q6: is "unrelated" a fact or a detection limit | treating an unmeasured edge as an absent one |
| A6 | **Repertoire evolution** | S10 presence/absence | ancestral-state reconstruction on the species tree (D15) | when each family arose and where it was lost | absence claims that fail D4's two bars |
| A7 | **Duplication order** | S7/S8 trees, S11 | reconciliation against the species tree | the order of the internal duplications that made the 4×6TM channels | a gene tree with poor support at the informative nodes |
| A8 | **Auxiliary-subunit inflation** | S20 | count `channel_associated` genes; compare with published human totals | how much of the 240–400 spread is scope rather than biology | published totals whose own scope is undocumented |
| A9 | **Channelopathy mapping** | S19, S12 | map ClinVar variants onto the pore module and the constrained core | which families' disease variants concentrate in the pore | variant sets of unequal ascertainment across genes |
| A10 | **Prokaryotic sisters** | S21 | best-profile assignment of prokaryotic hits to eukaryotic families | which superfamilies predate eukaryotes | horizontal transfer mistaken for vertical descent |
| A11 | **Pharmacology overlay** | S22 | approved-drug targets by family, against family size | whether drug targets concentrate in well-annotated families | target lists that are themselves annotation-biased |
| A12 | **Pore-module conservation** | S6 alignments | per-column conservation of the pore module across each superfamily | the constrained core each family shares | alignments dominated by one over-sampled lineage |

## The report skeleton

Which analysis feeds which manuscript section, so that a finished analysis
has a home and an orphan analysis is visible.

| Manuscript section | Fed by |
|--------------------|--------|
| *A declared scope for the channelome* | D23, S0, A8 |
| *A benchmarked classifier without gene symbols* | S1, A1 |
| *The census* | S2, S3, S5, A2 |
| *The forest* | S7, S8 |
| *Selectivity and phylogeny* | S9, A3 |
| *Mechanism cuts across clades* | S17, A4 |
| *What is convergent and what is undetected* | S18, A5, A12 |
| *Repertoire evolution across the tree of life* | S10, S11, A6, A7 |
| *Clinical and pharmacological weight* | S13, S19, S22, A9, A11 |
| *Methods, limits, and what a census cannot say* | A1, A2, H14, D23–D28 |

## Analyses deliberately not attempted

- **Predicting selectivity from sequence.** Hazard H14: Cys-loop charge
  selectivity flips with three substitutions and AMPA Ca²⁺ permeability is
  set by RNA editing. A predictor would be reporting family membership with
  extra steps.
- **A single tree of all ion channels.** D27.
- **Gating-mechanism inference.** HCN opens on hyperpolarisation and is
  architecturally indistinguishable from CNG, which is voltage-insensitive.
  Nothing in a sequence carries that difference reliably.
- **Counting channels in a species with no reference proteome.** The
  denominator would be undefined, which is the failure D23 exists to prevent.
