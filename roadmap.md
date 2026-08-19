# Feature-development history of the app itself

Distinct from `PUBLICATION_ROADMAP.md`, which tracks the *science*. This file
tracks the tool.

## Inherited (v1.0–v1.7, from `../piezo_genes` → `../ip3r_genes`)

- **v1.0** Tkinter GUI, NCBI / Ensembl / UniProt clients, disk cache,
  results bundles.
- **v1.1** Sequence analysis: MSA, identity matrix, NJ tree, clustering,
  mutation calling, conservation.
- **v1.2** Discovery scorer with the composite criteria and the D3 evidence
  gate.
- **v1.3** AlphaFold and Foldseek clients; structural evidence in the
  investigation pipeline.
- **v1.4** Ensembl Compara paralogue mining; BLAST bait search.
- **v1.5** Deep-dive investigation pipeline (seven evidence lines) and the
  case-file renderer.
- **v1.6** Domain-bait scan and the exhaustive hunt mode.
- **v1.7** Publication tooling: figure style, dashboard, manuscript assembly,
  claim checker, deposit manifest.
- **v1.7.1** (IP3R project) Ensembl symbol resolution switched from
  `/xrefs/symbol/` to `/lookup/symbol/` after the former was measured to
  stall indefinitely for `homo_sapiens`. Kept.
- **v1.7.2** (IP3R project) Sister-family test in the discovery scorer —
  a labelled-bait margin, the origin of this project's D14.

## This project

- **v2.0** (2026-08-19) **The catalogue.** `src/catalogue/` — 90 families,
  25 superfamilies, 16 hazards, with signature levels, provenance tags and
  machine-checked internal consistency. Replaces the parent's single
  `utils/family.py` with a subject that has to be looked up rather than
  assumed. `src/utils/scope.py` derives a run's scope from it and keeps the
  ported machinery working unchanged.
- **v2.1** **The classifier.** `src/classify/` — three independent tiers
  (architecture rules derived from the catalogue plus hand-written hazard
  rules; the selectivity-filter motif tier; reference identity with a
  margin), combined into one `ChannelCall` with a full audit trail, a
  confidence tier, and the gene symbol carried but never consulted.
- **v2.2** **The phylogeny forest.** `src/phylo/` — pore-module extraction,
  tier-1 and tier-2 tree builders that record every command they ran, a
  `NotAlignable` refusal for non-homologous superfamilies, and a fold network
  that cannot be turned into a tree.
- **v2.3** **Stdlib MAFFT wrapper** (`src/utils/mafft.py`) so alignment does
  not depend on Biopython, which is absent from the base interpreter. Raises
  rather than degrading (D28).
- **v2.4** Species table rebuilt for the tree of life (51 species, 15 groups,
  named panels) — the parent's vertebrate-plus-margins panel could not
  support a claim about a family that is bacterial, plant or invertebrate.

- **v2.5** `scripts/selftest.py` — the offline invariants, wired into the
  session protocol as step 4: catalogue validation, every hazard rule fired
  against the synthetic architecture it was written for, the classifier's
  blindness to gene symbols, the motif regex, the alignment helpers, the D27
  refusal and the scope's derivation of its sister panel. Under a second, no
  network.
- **v2.5.1** `ChannelFamily.architecture()` now brackets signatures shared
  with a known decoy, so `[PF00664×2] + [PF00005×2] + PF14396` says at a
  glance that only one of CFTR's three accessions is diagnostic. A display
  that hid them made every family look more distinguishable than it is.
- **v2.6** `src/gui/tools.py` split at 549 lines into `tools.py` (menu, task
  runner, the pipeline the CLI also runs) and `tools_extra.py` (the one-off
  investigations), and the Analysis menu gained **Classify selected row**,
  which runs the three-tier classifier and shows the full audit trail.

- **v2.7** (S1) `alignment_stats()` returns coverage of both sequences and
  the reference tier enforces a floor on the longer one (**D30**). Without
  it, covered-only identity — inherited from a project whose sequences were
  one family and all the same size — scored a 226-residue connexin at 61 %
  to a 4,967-residue ryanodine receptor. The benchmark found it; the
  self-test now locks it.

- **v2.8** (S1) The reference tier now inherits the architecture tier's
  result: it searches within the named ambiguity, or within the named
  superfamily when only that was reachable — which is the normal outcome
  for the Cys-loop and iGluR receptors, where every family carries the same
  two accessions. Falls back to the full panel if the restriction matches
  no exemplar.

- **v2.9** (S1) The reference prefilter applies the length ratio the D30
  coverage floor implies, before aligning. `kmer_containment` divides by the
  smaller profile, so the shortlist was filling with the longest references
  in the panel — the ones the floor then discarded. Runtime per protein went
  ~90 s → ~9 s, and the shortlists became the right ones: connexin-26 now
  shortlists connexin-43, not RYR2.

## Known gaps

- `src/analysis/alignment.py` still imports Biopython at module level, so
  `src/analysis`, `src/gui` and `src/cli` only run in the `piezo1` env.
  Migrating it onto `src/utils/mafft.py` would make the whole app
  stdlib-plus-`requests`.
- `src/gui` has no *catalogue browser*: a family can only be inspected from
  the command line (`run.py --catalogue --scope X --detail`).
- The discovery scorer's evidence strings still read as though they name
  paralogues, although the sister panel now comes from the hazard registry.
- `run.py --phylo <superfamily> --tier 2` aligns full-length exemplars and
  says so; real tier-2 trees need the pore-module extraction that S6 runs.
