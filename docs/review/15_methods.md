## 15. How the claims in this review were checked

This review is generated, and its verification is mechanical rather than
editorial. Three checks run on every build.

**Every reference was resolved against a live literature service.** The
source list (`scripts/review_sources.py`) contains *titles*, not citations.
`scripts/s0_review_refs.py` queries Europe PMC for each, accepts a record
only when its title is at least 90 % similar to the one requested, and
writes the returned PMID, DOI, year, journal and authors to
`results/s0_baseline/references.tsv`. A title that does not resolve is
reported and not written.

That check has teeth. On the first run, seven titles failed: five were short
generic titles where the search ranked another paper, and two were
misremembered — the Cav1.1 structure paper is *Structure of the voltage-gated
calcium channel Cav1.1 complex* (Science, 2015), and the Moran review is
*Evolution of voltage-gated ion channels at the emergence of Metazoa*, not
"of Nervous Systems". Both were corrected against the returned records. The
final list is 148 of 148 resolved.

**A citation with no verified reference is a build failure.**
`scripts/s0_review_build.py` collects every `[key]` in the section files,
renumbers them into order of first appearance, and refuses to write the
document if any key is missing from the reference table. Prose therefore
cannot cite a paper that was never resolved.

**Every `[db]` statement is re-derivable.** Section 12 reports measurements
rather than literature, and the queries that produced them are in
`scripts/s0_catalogue_verify.py`; the tables are in
`results/s0_baseline/`. Re-running the script re-derives every number in
that section, and the S0 report is rendered from the tables rather than
written alongside them.

### The computational toolchain

Where this review reports its own measurements, the tools are: MAFFT for
alignment [katoh2013], trimAl for alignment trimming
[capellagutierrez2009], IQ-TREE 2 with ModelFinder and ultrafast bootstrap
for phylogenies [minh2020, kalyaanamoorthy2017], HMMER for profile searches
[eddy2011], Foldseek for structure comparison [vankempen2024], and the
AlphaFold database for predicted structures where experimental ones are
absent [jumper2021, varadi2022]. Domain and sequence annotation comes from
InterPro, Pfam and UniProt [paysanlafosse2023, mistry2021, uniprot2023].

### What this does not verify

Resolution confirms that a paper exists with the title cited. It does not
confirm that the paper says what the sentence citing it claims. That
remains an editorial responsibility, and the mechanical checks should not be
mistaken for having discharged it.

Nor does it verify the *selection*: which papers are cited, and which
landmark results are missing, is a judgement, and a review assembled by
searching for known titles will under-represent work the assembler did not
already know about. The bibliography here is weighted towards structures,
identifications and reviews in high-visibility journals, and that weighting
is a property of how it was built.
