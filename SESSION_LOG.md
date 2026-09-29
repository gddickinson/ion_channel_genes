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

### What the S1 benchmark found — including a bug in the method it inherited
The first full run (97 proteins, 910 s) returned recall 48/72 and specificity
24/25, and three of its errors were not classifier mistakes but a broken
metric:

| call | what it looked like |
|---|---|
| GJB2 (connexin-26, 226 aa) → `ryr` | "nearest Hs_RYR2 at 61.3 %" |
| GRIN1 → `viroporin` | "nearest SARS2_E at 50.9 %" |
| GABRA1 → ambiguous with `ryr` | "nearest Hs_RYR2 at 46.8 %" |

A 226-residue connexin is not 61 % identical to a 4,967-residue ryanodine
receptor. **Covered-only identity had no coverage floor**: MAFFT places the
short sequence inside the long one, the 225 covered columns are by
construction the best-matching 225 of ~5,000, and the metric scores exactly
those. Measured: connexin-26 scored 61.3 % to RYR2 and **49.8 % to
connexin-43, its actual relative**.

Coverage of the shorter sequence does not catch this (0.996). Coverage of
the longer does — 0.045 against 0.589 — so `alignment_stats()` now returns
both, `reference.MIN_COVERAGE = 0.30` enforces the floor on the longer, the
number of dropped hits is reported with every call, and `selftest.py` locks
the asymmetry in. This is **decision D30**, and it is the clearest thing the
benchmark bought: the metric was inherited from a project whose sequences
were one family and all the same size.

Also fixed: the recall table mixed positives with decoys, so a decoy
correctly left `unassigned` read as a recall failure. `recall_table()` now
carries a `panel_role` column and the report scores the two separately.

The benchmark was then re-run with the floor in place.

### The literature review
Added after the benchmark, on the same rule as everything else here: no
claim without a check.

- `scripts/review_sources.py` — 148 `(key, title)` entries. **Titles, not
  citations.**
- `scripts/s0_review_refs.py` — resolves each against Europe PMC, accepting a
  record only when its title is ≥ 90 % similar to the one requested, and
  writes `results/s0_baseline/references.tsv` with the returned PMID, DOI,
  year, journal and authors. **148/148 resolved**, after the check caught
  seven failures: five short generic titles that needed an author/year
  constraint to rank, and two I had misremembered — the Cav1.1 paper is
  *Structure of the voltage-gated calcium channel Cav1.1 complex* (Science
  2015, not a 2016 Nature paper at 3.6 Å), and the Moran review is
  *…at the emergence of Metazoa*, not "of Nervous Systems".
- `scripts/s0_review_build.py` — assembles `docs/channel_review_2026.md`
  from `docs/review/*.md`, renumbers citations into order of first
  appearance, renders the bibliography, and **refuses to write the document
  if any cited key has no verified reference**. `--pdf` typesets via pandoc
  + xelatex.
- Result: 16 sections, 8,424 words, 148 references, 0 unused. §12 reports
  this project's own S0 measurements rather than literature, and §15 states
  what the mechanical checks do *not* verify — that a paper exists with the
  title cited is not that it says what the citing sentence claims.

### S1 results (final run)
Recall **50/72 (69.4 %)**, specificity **25/25 (100 %)**, 16/16 hazards
exercised, 866 s. Per-tier: 17 hazard + 12 architecture + 7 motif + 24
reference. Three fixes went in between the first run and this one — the D30
coverage floor, the superfamily-scoped reference search, and the prefilter
length rule — and the first of them turned a false positive into a correct
rejection.

### Review figures
Eight figures added, six rendered from committed tables and two labelled
schematics. New scripts: `s0_domain_map.py` (InterPro match coordinates for
the 22 proteins the architecture figures compare), `s0_filter_atlas.py` (the
two filter alignments, extracted from real sequences), `review_figlib.py`
(shared drawing primitives) and three `s0_review_fig_*.py` renderers.

Two things the figures taught, both fixed:

- **The first fold gallery was misleading.** It drew three subunits per
  panel and truncated the helix count to fit, so a panel labelled "6 TM"
  showed four. Rewritten to draw one subunit at its true count with a
  stoichiometry top-view beside it.
- **`s0_review_build.py --pdf` was silently producing a figure-less PDF.**
  Pandoc resolves image paths against the working directory, not the input
  file. Fixed with `--resource-path`, and the build now warns if a PDF with
  figures comes out implausibly small — the failure mode was a clean exit
  code and a wrong document.

The atlas also recorded a useful negative: *Bacillus* NaK is the one panel
member that fails the TxGYG test, and its filter reads **TVGDG** — the
single substitution that makes it non-selective. It is drawn as a labelled
negative control rather than dropped.

### Next session
**S2** — the uncapped enumeration. Read
`results/benchmark_controls/report.md` first; the three findings it hands S2
are in the roadmap's next-session note. Read
`results/s0_baseline/report.md` first. Do not skip: `verify_anchor()` before
reporting any filter result; the per-tier attribution table; the list of
hazards no panel member exercises.

---

## 2026-09-28 — S2: census v2 (uncapped enumeration + a call on every record)

**Ran.** `s2_enumerate.py counts` → `walk -j 6` → `verify`; `s2_classify.py`
(12 shards, up to 9 workers); `s2_census.py`; `s2_report.py`. External drive
attached throughout (671 GB free at start); raw pages 0.5 GB, census TSV
28.6 MB + FASTA 289 MB, all under `<data root>/raw_api/s2/`.

**Decided (D31).** Enumerate from UniProtKB rather than InterPro's
per-signature listing: the rules consult 19 accessions outside the pore union
and need copy numbers, and only UniProt's record carries the whole Pfam list
with `MatchStatus` counts. InterPro's count is fetched per signature as the
independent check. The union query is cut into taxonomic shards that must sum
to the union count. Reference tier not run at census scale (`not_run` on every
row, D28); four-repeat projection gated on `PF00520` ≥ 4 (the H1 trigger).

**Result.** 1,245,200 records; every completeness check exact (12/12 shards,
union, 67/67 signatures). Family calls 346,627 (27.8 %), superfamily-only
398,282 (32.0 %), unassigned 500,291 (40.2 %). Human census genes 319/320
enumerated, 173/320 right family, 1 wrong (ZACN). S1 panel: 71/71 channel
members enumerated, 0/26 non-channel members called channel.

**Operational notes.** The chordate shard was re-split mid-walk into
Mammalia / Actinopterygii / other (resume by cursor is safe: a page file is
written before its cursor). The four-repeat projection dominates runtime:
37,393 MAFFT alignments ≈ 3.5 CPU-hours; everything else classifies 1.2 M
records in about a minute. The call files lack the gene column (the
classifier never reads it); `s2_census.py` joins it back for scoring.

**Emergent.** Co-domain signatures (cNMP, SBP_bac_3, PAS_9) bring in ~376 k
non-channels; derived rules demand full architectures (MscS 65,506 partial);
uncalled four-repeat filters (NEEE, DEEA, …); ZACN → AChBP persists.

### Next session
**S3** — profile-HMM sweep. Read `results/census_v2/report.md` first; the
superfamily-only table is S3's worklist.

---

## 2026-09-28 — S3a: profile library + best-profile assignment → census v3a

**Split.** S3 as written could not run: its brief built profiles "from the
S6 alignment" (S6 depends on S3) and swept "the reference proteomes" (S4
declares them). Split into **S3a** (this session: profiles + assignment over
census v2) and **S3b** (S4 proteome sweep + jackhmmer, pending on S4).

**Ran.** `s3_build_profiles.py` (90 profiles, 28 s); `s3_benchmark.py`
(80 s); `s3_sweep.py prep` (1,187,702 unique sequences) and `search --jobs 5
--cpu 2` (90 profiles, 83 min wall, 6.9 profile-hours); `s3_census_v3.py`
(65 s, 1.9 GB RSS); `s3_external_check.py`; `s3_report.py`. Drive attached
throughout (670 GB free at start). Bulk under `<data root>/hmmer/s3/`.

**Caught on the way.** The first profile build aborted on its own
disjointness rule: `s0_lib.resolve_gene` matches synonyms, so TRPC7 resolved
to TRPM2 (and six more pairs). S3 now resolves by primary gene name, and a
within-family duplicate-accession check was added. The benchmark's first
decoy explanation ("too few seeds") was wrong — assoc_k_beta has 13 seeds;
the family pools unrelated proteins. Report corrected before commit.

**Result.** Benchmark LOO 69/71 on S1 (vs 50/72); orthologues 534/534
(vertebrate). Calibration vs S2 96.8 %. Census v3a family calls 58.7 %
(27.8 % in S2); conflicts 12,857; human 318/320. External check vs
`../ip3r_genes` census v6: 12,204/15,601 same call; the 502 swapped records
are all S2's `H4-itpr` absence rule. The user asked mid-session whether
PIEZO/IP3R/RyR could be taken from the parent projects instead of re-run:
no (profiles are needed as competitors, cost already paid, D28 and the
port-the-method rule) — they became the external check instead.

**Emergent.** H4, H2, H13 absence rules; PF00520 superfamily conflicts;
resolver synonyms; heterogeneous auxiliary families; no non-vertebrate test
for the Cys-loop/iGluR/DEG/P2X splits; 2,258 PIEZO-census records
unassigned. D32 recorded.

### Next session
**S4** — proteome scope. S3b waits on it. The S2 hazard-rule fixes (H2, H4,
H13) are emergent rows and can precede S4 if the user prefers.
