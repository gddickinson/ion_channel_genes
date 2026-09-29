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

---

## 2026-09-28 — S2b: hazard rules H2, H4, H13 rewritten as positive tests

**Why.** User-directed, straight after S3a, whose calibration and external
check traced most S2 errors to three rules that call a family from an
absence (CLAUDE.md forbids it).

**Measured first.** Which domains separate each pair in census v2 (with the
S3a profile call as the independent label): `PF02026` and `PF06459` each sit
only on profile-RyR records; PLAT/REJ/GPS only on polycystin-1-like;
`PF18109` only (1,358/1,358) on profile-TRPP. ITPR has no domain RyR lacks.

**Changed.** `rules.py`: H4-ryr + H4-ryr-tm positive, H4-core superfamily;
H13 PC1 by PLAT/REJ/GPS/PKD×5, TRPP by `PF18109`, core superfamily; H2
superfamily-only for any LBD without the TM region. Rule engine gained
whole-sequence bounds (TM, length, completeness). Catalogue: `trpp` PF20519
FAMILY → SHARED_WITH_DECOY. Hazard records rewritten. D33 recorded.

**H2 took three passes, and the record should show it.** (1) A positive
"soluble, ≤ 300 aa" AChBP test called nothing at all — which exposed a bug:
`tm_count` 0 was passed as None in `s2_classify.py` *and* `cli_channel.py`.
(2) With that fixed it called 6,644, mostly 106–150-aa scraps and flagged
fragments; adding a lower bound and a completeness test exposed a second
bug: UniProt's flag "Precursor" (the Lymnaea AChBP itself) was read as a
fragment. (3) The final shape test called 2,164, agreed with the profiles on
53/1,263, and ~75 % of its calls were in Ecdysozoa/Chordata, which have no
AChBP. Rejected: no positive AChBP test exists at the architecture tier.
Each r2 intermediate call set is archived (`calls_r2a/`, `r2b/`, `r2c/`).

**Result.** 108,407 records re-classified, 27,028 calls changed, 0 outside
the re-checked set. S1 benchmark unchanged (50/72, 25/25) — ITPR1 now by
reference, PKD1 now positively. S2 r2 human: 169 right, 0 wrong. Census v3a
re-merged: conflicts 12,857 → 3,288; ITPR/RYR swaps vs IP3R 502 → 0; human
319/320. The S1 re-run used an intermediate H2 rule; no panel protein can
tell the versions apart (ZACN has 4 TM helices).

**Emergent.** H11/H12 still absence rules; census `fragment` column
conflates Precursor; S3a R3 seeds drawn from r1 calls.

### Next session
**S4** — proteome scope (S3b waits on it).

---

## 2026-09-28 — S2c: hazard rules H11, H12 rewritten as positive tests

**Measured first.** H12: of PF02214 carriers without a pore module, 14,354
carry one of five KCTD C-terminal domains (KCTD5-like, KCTD8/12/16 H1,
KCTD10-like, KCTD1/15, SHKBP1/KCTD3) and the S3a profiles call 14,317 of
them KCTD, 3 otherwise; the shape alternative (T1, 0 TM, complete) was
contradicted on 160 — rejected for the same reason as AChBP. H11: every one
of the 980 "SUR" calls was a SUR: 918 bacterial cNMP + C39 peptidase ABC
transporters, 62 eukaryotic fused gene models, none SUR-shaped (checked:
1,300–1,800 aa with ABC/TMD0 domains only). ABCC8 carries nothing
SUR-specific.

**Changed.** H12-kctd-<marker> ×5 + H12-t1 (superfamily); H11-abcc removed;
catalogue gains the five KCTD markers (FAMILY on the control family — no
derived rule, no change to the pore union) and `assoc_sur` ABC domains
re-levelled. Self-test invariant now covers every hazard.

**Caught.** Checking the pore union before re-levelling the ABC domains
showed S2b had already dropped PF20519 from it (67 → 66): enumeration was
derived from evidence level. Decoupled (`Signature.enumerate`, D34), union
restored to 67 and pinned by a self-test against `signature_counts.tsv`.
Also: the revision-transition script was hard-wired to S2b and its own
guard fired when S2c changes appeared — parametrised (`--base/--head`).

**Result.** 63,411 records re-checked, 18,412 calls changed, none outside
the set. S2 r3; S1 unchanged; v3a conflicts 2,907; human 319/320.

**Emergent.** 525 records (508 bacterial) called SUR by the S3a profile
alone — no generic ABC-transporter decoy profile exists. Proposed fix adds a
catalogue family; left for the user.

### Next session
**S4** — proteome scope.

---

## 2026-09-28 — S3a2: an ABC-transporter decoy for the profile library

**Why.** S2c found the S3a SUR profile winning 525 records (508
bacterial) because no generic ABC transporter was on the menu. User asked
for the decoy family.

**Ran.** Catalogue: `nonchannel_abc_transporter` (controls.py; H11 now
lists it). Exemplars chosen from reviewed UniProt entries with their Pfam
architecture checked live (HlyB, SunT, MsbA, ABCC1). `s0_catalogue_verify.py`
(752 s, 793 requests, 0 failures) → reference panel +4, no existing sequence
changed (diffed against a snapshot); `s0_report.py`. `s3_build_profiles.py
--only` and `s3_sweep.py search --only` — both now merge into the frozen
records instead of rebuilding (new in this session); the 90 existing builds
and seed sets verified byte-identical. `s3_benchmark.py`, `s3_census_v3.py`,
`s3_external_check.py`, `s3_report.py`.

**Result.** Profile-called SUR 525 → 0; decoy 1,026 (997 bacterial);
benchmark unchanged; real SUR and CFTR keep wide margins over the decoy.
Conflicts 2,919 (+12 fusions). Family count 91 updated in CLAUDE.md,
README, INTERFACE, channel_background, findings_summary; the generated
literature review still says 90 (it describes the S0 catalogue).

### Next session
**S4** — proteome scope.

---

## 2026-09-28 — S4: proteome scope (the denominator)

**Ran.** Session protocol clean (drive attached, 669 GB free; catalogue and
self-test clean). New: `scripts/s4_proteome_lib.py`, `s4_proteomes.py`
(`manifest` / `download` / `verify` / `all`), `s4_report.py`; eight S4
self-test invariants. UniProt proteomes API + release README + per-proteome
`RELEASE.metalink` + NCBI Datasets, all archived under
`<data root>/raw_api/s4/`; 100 files (268 MB) to `<data root>/proteomes/s4/`;
sweep DB `panel_refprot.fasta` (822,499 seqs, 507 MB). Offline rerun
reproduced every table byte for byte.

**Result.** 52 species → 50 reference proteomes (release 2026_03) + 2
genome-only (*Cornu*, *Torpedo*). Selection rule picked the model strain in
all 6 multi-candidate species. 49/50 exact on README counts + MD5; human
reissued by UniProt on 2026-09-15 (README describes a withdrawn file of
147,503 entries; served file 20,652, one per gene, MD5 OK) → D35. 319/319
human census genes in the human proteome. 12,096 census v2 records are
panel-proteome entries.

**Surprises.** (1) NCBI's dataset_report omits superseded assembly versions
by default; 7 proteomes sit on one → `all_assemblies`. (2) The reference
proteome `.fasta` is one entry per gene in 50/50, gene2acc "groups" are not
a gene count. (3) **The viroporin family was never in the census search
space** (SUBFAMILY signatures, D34 default) — emergent. (4) INTERFACE said
51 species; the table has 52.

### Next session
**S3b** — profile sweep over the S4 panel DB + jackhmmer to convergence (D10).

---

## 2026-09-29 — S3b: panel profile sweep + jackhmmer → census v3

**Ran.** Protocol clean (drive attached, 668 GB free; catalogue and
self-test clean). New: `scripts/s3b_lib.py`, `s3b_sweep.py`, `s3b_kill.py`,
`s3b_jackhmmer.py`, `s3b_census_v3.py`, `s3b_report.py`; seven D10 self-test
invariants; `s3_sweep._search` gained an `out_dir` parameter. 91 profiles ×
822,499 panel entries (12.4 min wall); 68 jackhmmer runs (26.9 run-hours,
6.1 h wall at 5×2 threads). Every table re-derivable offline from the
archived domtbls and logs under `<data root>/hmmer/s3b/`.

**Result.** Census v3 on the panel: 49,595 entries with evidence; 8,755
census-family calls, 230 (2.6 %) outside census v2, 96 of them high
confidence (CLIC 27, Hv1 15, CNG 10, TRPV 8, pannexin 8, viroporins 3/3).
Instrument check vs S3a 99.8 %; human 319/320. jackhmmer: 24 clean, 44
killed (K1 23, K2 2, K3 19); clean-run completeness 3,352 / 3,356.

**Surprises.** (1) TRPN's "misses" are ankyrin-repeat proteins — D32's
coverage gate fails for repeat-dominated profiles (emergent). (2) The first
report counted auxiliary-subunit families as channel families
(`startswith("channel")` matched `channel_associated`) — caught, fixed.
(3) The ported K3 handling let non-converged runs contribute rounds 1–9:
30,639 candidates → 7,980 once removed (D36). (4) K1 fires at round 2 in
multi-family superfamilies; Cys-loop/DEG/P2X runs converge on one
superfamily set — jackhmmer is a superfamily instrument there.

### Next session
**S5** — genomic tblastn + miniprot for absent cells of
`family_by_species.tsv` and the two genome-only species (D4).

## 2026-09-29 — S5a: the genomic sweep's instrument, and a seven-genome pilot

**Ran.** Protocol clean (drive attached, 666 GB free; catalogue and
self-test clean). S5 split in the ledger (S5a instrument + pilot, S5b full
sweep), as the parent project did. New: `scripts/s5_genome_io.py`,
`s5_baits.py`, `s5_lib.py`, `s5_classify.py`, `s5_rescue.py`,
`s5_annotation.py`, `s5_calibrate.py`, `s5_ledger.py`, `s5_verdict.py`,
`s5_sweep.py`, `s5_report.py` (all < 360 lines); eight S5 self-test
invariants, the exon-confirmation one mutation-tested (the first version
passed under the mutation — its case never reached the threshold — and was
rewritten). Seven genomes downloaded and MD5-verified (6.0 Gbp); pilot
sweep 29.8 min wall on 1,808 baits.

**Result.** 1,881 loci, 1,020 profile-called. Matched detection 172 / 173
control cells; found 166 / 175. Mouse introns up to 996 kb at called loci
(13 / 220 over 200 kb). Zero cells: 10 absent among the 25 informative
non-*Cornu* cells, 0 proteome misses of a high-confidence intact gene;
literature checks 4 / 4 consistent. *Cornu*: 46 families present by genome
alone.

**Surprises.** (1) The first panel (8 baits/family, one per group) found
only 5/9 of *Arabidopsis*'s own families with self-baits excluded: TPC1 drew
nothing from 7 non-plant TPC baits in either miniprot or tblastn → one bait
per species, and a *matched* control (D37). (2) `spread()`'s default cap was
bound at definition, so the first comparison silently rebuilt the old panel
(664 baits) — caught by the count, fixed by explicit parameters. (3) Python
`hash()` salting would have made the tblastn cache never hit — fixed with a
digest. (4) The N-run scan was a per-character Python loop over 5 Mb loci —
regex. (5) "Overlaps an annotated gene" confirmed a 400 kb chained ZAC model;
reciprocal extent then rejected real long-UTR genes (mouse 223 → 130);
CDS-on-exons settled it (220). (6) Mouse VDAC: 12 high-confidence loci for
3 genes, retrocopies — the `intact` requirement.

### Next session
**S5b** — `s5_sweep.py run --all` over the 52 genomes (~42 Gbp, ~3 h
miniprot + downloads), after deciding the `-G` margin and the D4 fallback
bar (emergent rows).
