# SESSION_LOG.md

Running notes per session: what ran, what resulted, what's next.
The biology goes in `FINDINGS.md`; this file is the operational record.

---

## 2026-08-19 — setup session: project created

**What this session did.** Created the project as a port of
`../ip3r_genes`, retargeted from one gene family to every ion channel, and
ran the S0 catalogue verification against live databases.

### Ported unchanged
`src/core`, `src/databases`, `src/analysis`, `src/discovery`,
`src/investigation`, `src/gui`, most of `src/utils`, and the
`scripts/` publication toolchain (`dashboard.py`, `figstyle.py`, the
`s14_*` manuscript assembly and claim checker, `s1_toolchain.py`,
`build_findings_page.py`). Docstrings retargeted; no logic changes except
where noted below.

### Written new
- **`src/catalogue/`** (12 modules) — the subject. 90 families, 25
  superfamilies, 16 hazards, machine-validated by `python3 -m src.catalogue`.
- **`src/classify/`** (5 modules) — the three-tier classifier.
- **`src/phylo/`** (3 modules) — pore-module extraction, the tier-1/tier-2
  tree builders with the `NotAlignable` refusal, the tier-3 fold network.
- **`src/utils/scope.py`** — replaces the parent's `utils/family.py`; derives
  a run's scope from the catalogue and keeps the ported machinery working.
- **`src/utils/mafft.py`** — stdlib-only MAFFT wrapper, because
  `src/analysis/alignment.py` imports Biopython at module level and Biopython
  is not in the base interpreter.
- **`src/cli_channel.py`**, new `run.py` flags, eight presets,
  `scripts/s0_lib.py`, `scripts/s0_catalogue_verify.py`,
  `scripts/s0_report.py`, `scripts/s1_benchmark.py`, `scripts/s1_report.py`.
- Docs: `CLAUDE.md`, `INTERFACE.md`, `PUBLICATION_ROADMAP.md` (26 ledger
  rows, D1–D18 inherited + D23–D28 new), `docs/channel_background.md`,
  `docs/scope_and_boundaries.md`, `docs/classification_rules.md`,
  `docs/phylogeny_protocol.md`, `docs/session_briefs.md`,
  `docs/analysis_catalogue.md`, `roadmap.md`.

### What ran
- `python3 -m src.catalogue` — validation clean.
- `python3 scripts/s1_toolchain.py` — 12/12 external binaries resolve;
  `results/toolchain_manifest.txt` written. Biopython reports MISSING from
  the base interpreter and is present in the `piezo1` env (D18).
- `python3 scripts/s0_catalogue_verify.py` — three passes (751 s / 627 s /
  666 s; 746, 762 and 769 live requests; 0 request failures in any).
  **Third pass clean**: 115/115 Pfam accessions verified, 162/162 exemplars
  resolved, 52/52 taxon ids, 162 sequences cached to
  `results/s0_baseline/reference_panel.fasta`.
- The four-repeat selectivity-filter projection was validated by hand before
  being written into `src/classify/motifs.py`: 6/6 correct
  (CACNA1C→EEEE, CACNA1G→EEDD, NALCN→EEKE, SCN1A/SCN11A/SCN5A→DEKA)
  at ~1.1 s per query.

### What the first S0 pass found, and what was done about it
| Finding | Action |
|---|---|
| 113 Pfam accessions, 0 missing, 2 short-name mismatches (`PF13833`, `PF15149`) | catalogue corrected |
| `kir` declared `PF07885`; KCNJ2/KCNJ11 carry none | declaration corrected; **hazard H16's discriminator rewritten** — it is two different models, not one at two copy numbers |
| `kca_slo`: KCNMA1 carries `PF00520`, KCNT1 carries `PF07885` | recorded in the family note as H8 inside one family |
| `clc_prokaryotic` exemplar declared gene `eriC`; UniProt calls it `clcA` | corrected, old name kept in the note |
| 13 exemplars unresolvable by gene symbol (almost all prokaryotic/invertebrate) | 12 resolved by protein-name search and their accessions added; the Aplysia AChBP exemplar dropped (no reviewed entry) |
| **`organism_id:` is an exact taxonomy node and misses strain-level entries** — GLIC is filed under *Gloeobacter violaceus* strain PCC 7421, so `organism_id:33072` returned nothing | `resolve_gene()` and `resolve_symbol()` switched to `taxonomy_id:`; this is why every prokaryotic exemplar was invisible |
| Taxonomy: 51 ids checked, 1 name mismatch (HIV-1) | species table corrected; *Arcobacter* → *Aliarcobacter*, *Methanocaldococcus* → *Methanothermobacter*, *Helix aspersa* → *Cornu aspersum* added |
| 29 exemplar architecture mismatches | most are subfamily-level signatures absent from a family's other members (expected); the informative ones are listed above |

Exemplars carrying a verified accession: **94 → 112** of 162.

### Corrections made after the second S0 pass
| Finding | Action |
|---|---|
| 7 exemplars flagged `DECLARED_GENE_MISMATCH` turned out to be entries with **no gene name in UniProt at all** (NavAb, ELIC, GluR0, FaNaC, NaK, the Klebsiella bestrophin, AChBP) | the checker gained a distinct `NO_GENE_NAME_IN_UNIPROT` status, so real mismatches are not buried in noise |
| GLIC and ELIC carry `PF02931` and **no** `PF02932` — the same shape as ZACN | **hazard H2 rewritten**: three measured false negatives, including the Cys-loop superfamily's own rooting outgroup. Rule kept, exceptions recorded |
| GluR0's measured architecture is `PF07885` + `PF00497` — a K+ channel pore plus a bacterial binding domain, with no iGluR model | `iglur_prok` signatures corrected to what was observed; the "inverted P-loop" story is now a `[db]` claim rather than a `[lit]` one |
| `PF16799`, the positive Hv1 test, is absent from *Ciona* Hv1 | **hazard H9 amended**: the test is a mammal-only instrument |
| Prokaryotic NavAb carries neither `PF06512` nor `PF11933` | recorded in the `nav` family note: the Nav rule cannot reach its own outgroup |
| TRPM7 carries `PF16519` (TRPM_tetra) and `PF02816` (alpha-kinase) | both added to the `trpm` declaration |
| HIV-1 Vpu had no accession | `P05919` added |

### Code corrections found by running the benchmark
| Finding | Action |
|---|---|
| `derived_rules()` folded `SUBFAMILY` signatures into the family requirement, so `kv_shaker` required the Kv2-only `PF03521` and KCNA1 — a Kv1 — fell through to `kv_modifier` | generator rewritten to emit a family rule and a higher-priority subfamily rule. KCNA1 now returns *ambiguous between `kv_shaker` and `kv_modifier`*, which is the honest architecture-tier answer |
| Copy-number requirements were generated from the catalogue and are brittle at family scale (OTOP1 carries `PF03189`×3, OTOP2 ×2) | copy requirements removed from derived rules; the two places they are diagnostic are hazard rules with the measurement attached |
| The reference tier prefiltered globally, so an architecture-tier ambiguity could survive a tier asked to resolve it | prefilter now runs over only the ambiguous families' exemplars |
| **93 of the 97 S1 panel proteins are themselves catalogue exemplars** — both come from the same curated gene lists, so the reference tier would have scored nearly the whole panel against itself at 100 % | `classify(..., leave_one_out=True)`; the benchmark excludes the query's own accession and reports the overlap. Now **decision D29** |
| `src/gui/tools.py` was 549 lines, over the project's 500-line rule | split at the seam between the core pipeline and the one-off investigations → `tools_extra.py:ToolsExtraActions`, and a **Classify selected row** action added to the Analysis menu |
| Nothing checked the invariants offline | `scripts/selftest.py` — catalogue validation, every hazard rule fired against a synthetic architecture, symbol-blindness, the motif regex, the alignment helpers, the D27 refusal, scope derivation. Wired into the session protocol as step 4 |

### Next session
**S1** — the classifier control benchmark. Read
`results/s0_baseline/report.md` first. Do not skip: `verify_anchor()` before
reporting any filter result; the per-tier attribution table; the list of
hazards no panel member exercises.
