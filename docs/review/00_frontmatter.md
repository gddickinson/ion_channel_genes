# The ion channels: one function, twenty-five origins

*A review of the ion-channel superfamilies, their folds, their evolution and
what the sequence databases do and do not record about them.*

**Generated file — do not edit.** The text lives in `docs/review/*.md` and is
assembled by `scripts/s0_review_build.py`. Citations are written as stable
keys and renumbered on build; every reference was resolved against Europe PMC
by `scripts/s0_review_refs.py`, and a cited key with no verified record is a
build error rather than a silent entry in the bibliography.

---

## Summary

An ion channel is a gated aqueous pore. That is a description of what a
protein does, not of where it came from, and the distinction organises
everything that follows. The proteins the term covers do not form a family,
a superfamily, or any other clade: they are on the order of twenty-five
independent inventions, built on at least ten unrelated folds, that arrived
at the same solution to the same physical problem. There is no alignment
containing a nicotinic acetylcholine receptor and a potassium channel, and
therefore no tree.

This review is organised around that fact. We take the folds one at a time —
the P-loop superfamily that supplies nearly half the human channel genes, the
four unrelated ligand-gated superfamilies, the anion channels, the large
pores, the mechanosensors and the channels of intracellular membranes — and
ask of each what is genuinely shared and what merely looks shared. We then
turn to the two questions a modern census has to answer before it can count
anything: what the annotation databases actually record about these
proteins, and where their coverage fails.

The second question turns out to be sharper than expected. Domain models
were built family by family from the best-studied members, and they behave
accordingly: the canonical pore model `PF00520` is carried by twenty
catalogued families including a phosphatase that is not a channel, while most
of the TRP families carry no copy of it at all. Three of the four
prokaryotic outgroups that root the major superfamily trees cannot be found
by the domain rules that define their own superfamilies. These are not
database errors. They are what a set of models built from well-studied
proteins looks like when it is asked about the proteins nobody studied.

## Figures

| | |
|---|---|
| **1** | the folds, one subunit each — *schematic* |
| **2** | the potassium filter across the branch — *measured, no alignment used* |
| **3** | the four-repeat locus, DEKA / EEEE / EEDD / EEKE — *measured* |
| **4** | length range by superfamily — *measured* |
| **5** | which families carry which shared domain — *measured* |
| **6** | a rooted ML tree of the Cys-loop receptors — *measured* |
| **7** | the forest, and the four superfamilies that get no tree — *schematic of a rule* |
| **8** | the architecture traps, to scale — *measured* |

Six of the eight are rendered from committed tables in
`results/s0_baseline/` and `results/phylogeny/`; the two schematics say so
in their legends and carry no data beyond the counts they display.

## Key points

- **"Ion channel" is a functional class, not a clade.** Any statement of the
  form "ion channels evolved…" is either about one superfamily or is wrong.
- **The human channelome is lopsided.** Of ~320 pore-forming genes, 143 are
  in the P-loop superfamily and 79 of those are potassium channels; the
  remaining 177 are spread across twenty-four superfamilies.
- **Convergence is the rule.** Potassium selectivity evolved at least twice;
  vertebrates run two unrelated gap-junction systems; mechanosensitivity has
  been invented at least six times.
- **Selectivity and gating are not recoverable from sequence family.**
  Cys-loop charge selectivity flips with three substitutions; AMPA-receptor
  calcium permeability is set by RNA editing.
- **Domain annotation does not partition the subject the way the textbooks
  do**, and a census built on the obvious accessions silently loses whole
  divisions.
