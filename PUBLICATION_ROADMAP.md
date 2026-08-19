# PUBLICATION_ROADMAP.md — a census, classification and phylogeny of the ion channels

**Goal.** Identify every ion channel in a declared search space, classify
each one by positive tests with a recorded audit trail, and reconstruct the
phylogeny that the data actually supports — a forest of within-family and
within-superfamily trees plus a structural network between them. To the
evidence standard of an MBE / GBE / *Genome Research* paper: an enumerated
search space with a completeness argument, a benchmarked classifier with
per-family recall and per-hazard specificity, ML phylogenies with support
values, and one-command reproducibility.

**How the claim is scoped.** "A systematic census across N reference
proteomes spanning M lineages" (N and M fixed in S4), never an unbounded "all
ion channels". Every absence claim is a claim about a *declared* search
space, and every count is a count under the scope decisions recorded in
`docs/scope_and_boundaries.md`.

**The biology this rests on** is `docs/channel_background.md`, which tags
every statement `[db]` / `[lit]` / `[open]`. Read it before S1.

**The three things that will bite you.**

1. **Domain architecture is not a family call.** `PF00520` is carried by 20
   catalogued families including a phosphatase. Nav, Cav, NALCN and CatSper
   are architecturally identical. See Decisions **D25**, **D26**.
2. **Most TRP families carry no `PF00520` at all** (measured). A census
   enumerating P-loop channels by that accession loses the TRP division
   silently. Hazard **H7**.
3. **There is no tree of ion channels.** `build_tier2()` refuses
   non-alignable superfamilies and there is no tier 3. Decision **D27**.

---

## Session protocol

Claude: follow this protocol in every session that touches this project.

### Start of session
0. Open the project dashboard: `python3 scripts/dashboard.py --open` and
   start its watcher in the background (`python3 scripts/dashboard.py
   --watch`) so the page stays live through the session. Stdlib-only.
1. Read this file top to bottom. Read `docs/session_briefs.md` for the brief
   of the task you are about to work.
2. `git pull`.
3. Run `python3 -m src.utils.data_root --require`. It exits non-zero if the
   bulk-storage drive is not attached. **If it fails, stop and tell the
   user** — do not download anything, and do not "temporarily" use the
   internal disk. Note free space (`df -h`) before any bulk task.
4. Run `python3 -m src.catalogue` and `python3 scripts/selftest.py`. Both
   must be clean before anything else runs: the catalogue is the project's
   definition of its own subject, and the self-test is what proves the
   hazard rules, the D27 refusal and the symbol-blindness still hold after
   whatever the last session changed.
5. Pick the task: the single ledger row whose Status is `in_progress`, else
   the topmost `pending` row whose dependencies are all `completed`.
   Announce it. Set its Status to `in_progress` (commit at session end).
6. Work **only that task** until its completion criteria pass. Do not start
   the next task even if time remains — spend surplus on tests, docs, or
   hardening of the current task.

### During the session
- Big/raw files → `get_data_root()` subdirs, never the repo. Committed
  artefacts: summaries, TSV/CSV ≤ ~20 MB, figures, manifests, code.
- Record every load-bearing number and file path in the task's Results
  column (link out to a file under `results/` if long).
- Write `results/session_live.json` from any long-running driver
  (`{"task": "S5", "workers": 2, "steps": [{"label": …, "done": …}]}`) so the
  dashboard shows real progress. `scripts/s0_lib.py:live_progress()` does it.
- A new confusable pair discovered mid-task is a **hazard row**
  (`src/catalogue/hazards.py`) with a discriminating test, not a note.
- New leads, surprises, or scope changes → add a row to **Emergent tasks**.
  Never silently expand the current task.
- A task too large for one session: split it in the ledger (S5 → S5a/S5b),
  complete the first part properly, leave the rest `pending` with a note.
- Anything needing an interactive login or sudo: ask the user to run it via
  the `!` prefix.

### End of session (checklist — do all of it)
1. Completion criteria all pass → Status `completed` + date. Not all passing
   → stays `in_progress` with a NEXT note giving the exact resume point.
2. Update this file (ledger, Results, Emergent, Decisions).
3. Append a dated entry to `SESSION_LOG.md`.
3a. **Append a biological-findings entry to `FINDINGS.md`** — what the task
   changed in the biological story, plain language, findings before methods,
   unconfirmed claims marked *(pending: which task confirms it)*.
3b. **Refresh `README.md`** — the state-of-the-project summary for someone
   arriving cold.
4. `git add -A && git commit` (message: `S<n>: <one-line outcome>`) and
   `git push` if a remote is configured.
5. Tell the user: task status, headline results, what the next session does.

---

## Storage

- Active bulk-data root: see `data_root.txt`, read via
  `src/utils/data_root.py:get_data_root()`. Currently
  `/Volumes/FANTOM/ION_CHANNEL_DATA` — **an external drive that must be
  attached before any bulk session**.
- `require_data_root()` raises rather than falling back, so an unplugged
  drive stops the session at the top instead of after 200 GB has landed on
  the laptop.
- Expected bulk footprint: reference proteomes across ~60 lineages ~25 GB;
  HMMER raw output ~20 GB; AlphaFold structures for the reference panel
  ~5 GB; Foldseek databases ~15 GB; per-genome sweep evidence ~30 GB.
- Budget rule: any single download > 20 GB needs a note in Decisions with the
  running disk total.

---

## Task ledger

Statuses: `pending` / `in_progress` / `completed YYYY-MM-DD`.
Full step-by-step briefs: `docs/session_briefs.md`.

> **Next session: S1.** The setup session built the catalogue, ran S0 to a
> clean third pass, and ran the S1 benchmark once; read
> `results/s0_baseline/report.md` and `results/benchmark_controls/report.md`
> before doing anything else. S0's report is the list of things the
> catalogue claims that a live database did not confirm — 34 of 162
> exemplar architectures still differ from their family's declaration, and
> each one is either an expected subfamily-level absence or a correction
> waiting to be made.
>
> **Three things S1 must not skip.** (1) `verify_anchor()` on the four-repeat
> filter reference before any filter result is reported — a UniProt
> sequence-version bump would otherwise shift the numbering and make every
> Nav/Cav call silently wrong. (2) The per-tier breakdown: if most correct
> calls are made by the reference tier, the classifier is a nearest-neighbour
> lookup and will fail on exactly the sequences a census exists for. (3) The
> hazard table's zero rows — a hazard no panel member exercises is untested,
> not solved.
>
> **The external drive was not attached during setup or S0** (neither needs
> bulk storage). S4/S5 cannot start until it is.
>
> **This project is the IP3R project's machinery pointed at a much larger
> subject.** The app (`src/core`, `src/databases`, `src/analysis`,
> `src/gui`, `src/investigation`, `src/discovery`), the figure style, the
> dashboard, the manuscript-assembly and claim-checking tooling and this
> protocol are ported from `../ip3r_genes`, itself ported from
> `../piezo_genes`. The **Decisions** section carries their methodological
> scar tissue forward as pre-agreed rules. What is *not* ported is any
> result. ITPR and RYR appear here as two of ninety families, and the parent
> project's census of them is an external check on this one — never an input.

| ID | Task (one session each) | Depends | Status | Results (headline) |
|----|-------------------------|---------|--------|--------------------|
| S0 | Catalogue verification + scope confirmation: re-derive every Pfam accession, exemplar, architecture and taxon id against live databases; build the reference panel | — | completed 2026-08-19 | Catalogue built and verified: **90 families / 25 superfamilies / 320 human census genes / 16 hazards**, `validate()` clean. **Three passes** (769 requests, 0 failures on the last); clean on the third: **115/115 Pfam accessions verified, 162/162 exemplars resolved, 52/52 taxon ids**, 162 sequences cached to `reference_panel.fasta`. **Seven catalogue corrections came out of it**, incl. hazard **H16** rewritten (Kir carries *no* `PF07885` — it is two models, not one at two copy numbers), **H2** amended (three measured false negatives: ZACN, GLIC, ELIC — the Cys-loop outgroup fails its own superfamily's rule), **H9** amended (`PF16799` is absent from *Ciona* Hv1, so the Hv1 test is mammal-only), and `iglur_prok` corrected — GluR0's measured architecture is `PF07885` + `PF00497`, a **potassium-channel pore plus a bacterial binding domain**. Measured: `PF00520` spans 20 families incl. a phosphatase (H9); TRPC3/TRPM8/MCOLN1/PKD2 carry **no** `PF00520` (H7); ANO1 ≡ ANO6 (H6); CFTR ≡ ABCC8 minus `PF14396` (H11); KCTD1 carries the Kv T1 domain (H12). Four-repeat filter projection **6/6** correct. → `results/s0_baseline/report.md` |
| S1 | Classifier control benchmark: positives from every census family, one decoy per hazard; per-family recall, per-hazard specificity, per-tier attribution | S0 | in_progress | Panel built and running: 97 proteins (72 positives across the census families, 25 decoys — one per hazard). Evidence fetched for all 97; **93 of them are themselves catalogue exemplars**, so the run classifies leave-one-out (**D29**). Anchor validates, MAFFT v7.526, 162-sequence reference panel. **NEXT: the run is in `scripts/s1_benchmark.py`; when it finishes, `python3 scripts/s1_report.py` and then read the per-tier attribution table and the untested-hazard list before anything else.** |
| S2 | Uncapped InterPro enumeration of every pore signature → census v2, with a positive family call on every record | S1 | pending | |
| S3 | Profile-HMM sweep (one HMM per family, best-profile assignment) over reference proteomes + jackhmmer-to-convergence completeness argument → census v3 | S1, S2 | pending | |
| S4 | Proteome scope: declared reference-proteome manifest across the lineage panel (the denominator) + download tooling | S1 | pending | |
| S5 | Genomic tblastn + miniprot sweep for families the proteomes miss → per-lineage ledger (found / lost / assembly-gap) → census v4 | S4 | pending | |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl per family, and pore-module extraction for every tier-2 unit | S2, S3, S5 | pending | |
| S7 | **Tier-1 phylogenies**: IQ-TREE 2 + ModelFinder + 1000 UFBoot for every census family, rooted on its declared outgroup | S6 | pending | |
| S8 | **Tier-2 pore-module phylogenies** for every alignable superfamily, plus the **tier-3 fold network** (Foldseek/TM-score) for the refused ones | S6, S7 | pending | |
| S9 | Selectivity-filter atlas: the filter locus for every P-loop family, and whether it is congruent with the tier-2 tree (Q5) | S7, S8 | pending | |
| S10 | Repertoire evolution: presence/absence of every family across the lineage panel; ancestral reconstruction; the MscS animal-absence question (Q4) | S5, S8 | pending | |
| S11 | Duplication history: which families expanded in which lineage, 2R/3R ohnologue status, and the Kv/Nav/Cav repeat-duplication order | S7, S10 | pending | |
| S12 | Structures: AFDB coverage per family, TM-align against the experimental reference set, Foldseek all-vs-all for the fold network | S2, S6 | pending | |
| S13 | ML selection tests on the families with clinical variant sets (CFTR, SCN1A, KCNQ1, RYR1) | S6, S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3, S5, S7, S8, S9, S10, S12, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Method contribution** — what would domain search alone have missed, per superfamily (Q3); the recall curve as methods are added | S2, S3, S5 | high | pending | |
| S16 | **Annotation-quality audit** — how often a real channel locus is missing, fragmentary, split, unnamed or filed under the wrong family across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5, S15 | high | pending | |
| S17 | **Mechanism vs clade** — do the CLC channel/transporter and anoctamin channel/scramblase splits survive the tier-1 trees (Q7) | S7 | high | pending | |
| S18 | **Convergence audit** — connexin vs pannexin, TMEM175 vs the GYG channels, the mechanosensitive families: is "unrelated" a fact or a detection limit (Q6) | S8, S12 | high | pending | |
| S19 | **Channelopathy map** — clinical variants of every family mapped onto the pore module and the constrained core | S6, S12 | medium | pending | |
| S20 | **Auxiliary-subunit census** — the excluded 15 %: how many there are, and how often they are counted as channels in published totals | S2 | medium | pending | |
| S21 | **Prokaryotic and viral channels** — the outgroup census, and whether the eukaryotic families have prokaryotic sisters | S3, S5 | medium | pending | |
| S22 | **Pharmacology overlay** — which families carry approved-drug targets, against family size and annotation quality | S2, S19 | low | pending | |
| S14c | **Manuscript rewrite pass** — one full pass once every analysis has landed | S14a, S24 | medium | pending | |

---

## Emergent tasks & new aims

Anything discovered mid-session that deserves its own work goes here rather
than expanding the task in progress.

| Added | From | Task | Status |
|-------|------|------|--------|
| 2026-08-19 | setup | Confirm the external drive is attached and `data_root.txt` points at it before S4/S5 | open |
| 2026-08-19 | setup | **The TRP pore-model finding (H7) may generalise.** Pfam's coverage of pore modules was measured only on the exemplars. Enumerate which *census families* have any pore model at all, and how many members of each carry it — the answer sets S2's real recall ceiling | open (S2) |
| 2026-08-19 | setup | **`ZACN` breaks the H2 rule.** A human Cys-loop channel with no annotated TM domain. Check whether this is a UniProt annotation gap or a real truncation, and whether other single-exemplar families have the same problem | open (S1) |
| 2026-08-19 | setup | **69 exemplars carry no UniProt accession** in the catalogue and are resolved live each run. Fold the S0-resolved accessions back into the catalogue files so the reference panel is reproducible offline | open (S1) |
| 2026-08-19 | setup | **The four-repeat anchor works; the equivalent for other superfamilies does not exist yet.** A Cys-loop charge-selectivity anchor and a CLC gating-glutamate anchor would close hazards H14 and H5 the same way | open (S9) |
| 2026-08-19 | setup | Biopython is present in the `piezo1` env but absent from the base interpreter, so `src/analysis` and the GUI only run there. Either pin the env in a wrapper script or drop the Biopython dependency from `alignment.py` in favour of `src/utils/mafft.py` | open |
| 2026-08-19 | S0 | **Three of the four superfamily outgroups are unreachable by their own superfamily's rules.** GLIC and ELIC carry the Cys-loop LBD and no TM model, so the H2 rule rejects them; GluR0's architecture is a *potassium-channel* pore plus a bacterial binding domain, with no iGluR model at all; prokaryotic NavAb carries neither Nav-specific domain. Rooting therefore depends on S3's profile methods, not on domain search — and that is a result about annotation coverage, not a bug | open (S3 → S7/S8) |
| 2026-08-19 | S0 | **`PF16799`, the positive Hv1 test, is on human HVCN1 and not on *Ciona* Hv1.** The H9 discriminator is currently a mammal-only instrument; check how many other family-level tests are human-only by re-running the S1 panel on non-human orthologues | open (S1) |
| 2026-08-19 | S1 | **The reference tier does not scale to a census.** Measured on the S1 run: ~2.5 s per MAFFT pairwise alignment × up to 12 prefiltered candidates ≈ 30 s per protein, so the 97-protein panel takes ~45 min and a 5,000-protein census would take ~40 hours. S2 needs either a calibrated fast identity estimate (DIAMOND or BLAST, checked against MAFFT on the S1 panel so the D7 margin keeps its meaning) or a profile-based assignment that skips pairwise scoring entirely. **This is a scaling result, not a bug** — the tier is correct and slow, and the fix must not change what the margin means | open (S2/S3) |
| 2026-08-19 | setup | **The phylogeny pipeline was smoke-tested on the Cys-loop superfamily** (8 exemplars, MAFFT → trimAl → IQ-TREE 2 with 1000 UFBoot, rooted on GLIC/ELIC, 52 s) and returned the expected topology: anion- and cation-selective receptors separate at 100 %, AChBP sisters the cationic clade. `results/phylogeny/tier2_cysloop/`. **This is not S8's tree** — eight sequences, exemplars only — and must not be cited as one | closed 2026-08-19 |
| 2026-08-19 | S0 | **`PF00005` (ABC_tran) has 1.66 million UniProt proteins and `PF00520` has 206,115.** The census search space is dominated by two accessions that are mostly not channels. S2 needs a per-signature triage rule before it enumerates, or it will fetch a million transporters to find one CFTR | open (S2) |

---

## Decisions log

**D0 — This project inherits the PIEZO and IP3R projects' methodological
decisions.** They were paid for over 25+ sessions and are not to be
re-derived. D1–D18 below are those rules, restated. D23–D28 are new and
specific to a multi-family subject. A session may overturn one, but must say
so explicitly here with its reason.

**D1 — Bulk storage lives on the external drive**, `data_root.txt`. Sessions
call `require_data_root()`; an unplugged drive stops the session rather than
filling the internal disk.

**D2 — The repository stays private until S14b**, which flips it public
alongside the Zenodo DOI.

**D3 — The discovery scorer needs an evidence gate.** A score ≥ 40 needs at
least one family-specific component; size and novelty alone cap at 39.

**D4 — An absence claim must pass two bars**: a genome-wide contiguity floor
*and* a local check that the locus's own neighbourhood is present.

**D5 — Bait panels are screened by label, not padded for breadth.**

**D6 — Never issue a database correction without the integrity veto.**

**D7 — Best-hit assignment needs a margin.** Two references within ~10 % of
each other is not a call. Implemented as `DEFAULT_MARGIN = 0.10` in
`src/classify/reference.py`.

**D8 — Representatives are chosen per clade × per kingdom** by a rule in a
script, never "the longest sequence per species".

**D9 — Quote an annotation claim with its annotation source.**

**D10 — Iterative searches need a coded kill criterion.**

**D11 — Look at the figure.** Any session that changes a figure looks at it;
any session that writes a legend looks at the figure it describes.

**D12 — Every load-bearing number in the manuscript needs a claim row** in
`scripts/s14_claims.py`.

**D13 — Reports are rendered from the committed tables**, never written by
hand alongside them.

**D14 — Family assignment is a positive test at every stage.** Inherited
from the IP3R project, where ITPR/RYR was the single instance; here it is the
general rule, and the sixteen instances are `src/catalogue/hazards.py`.

**D15 — The species tree is an input, not a result.**

**D16 — Cross-family comparisons are paired within genome.**

**D17 — Permutation nulls are drawn from real genomic windows.**

**D18 — The project runs in the reused `piezo1` conda env**
(`/opt/anaconda3/envs/piezo1`, python 3.11.15) for anything needing
Biopython, matplotlib, BLAST+ or Foldseek. MAFFT, HMMER, trimAl and IQ-TREE 2
come from Homebrew and work in any shell. The catalogue, classifier and
phylogeny drivers deliberately need neither — stdlib plus `requests` — so the
session-protocol checks run anywhere. Versions:
`results/toolchain_manifest.txt`.

**D23 — The scope of "ion channel" is declared, not assumed.** Six boundary
questions, six recorded answers, in `docs/scope_and_boundaries.md`. Anything
excluded stays in the catalogue with a status, so an exclusion is a decision
rather than an omission. Aquaporins are the control that proves the boundary
is enforced.

**D24 — Membership and mechanism are separate calls.** CLC-3…7 are in the
CLC family by descent and are antiporters by mechanism; ANO3–10 are
anoctamins that scramble lipids. The classifier reports the family; the
mechanism is a separate, literature-derived field, and the two are never
merged into one label.

**D25 — A signature is evidence only at the level where it is diagnostic.**
`Signature.level` is `SUPERFAMILY` / `FAMILY` / `SUBFAMILY` / `ACCESSORY` /
`SHARED_WITH_DECOY`, and the rule engine will not make a family call on
superfamily evidence. Measured: `PF00520` spans 20 families.

**D26 — The four-repeat families are separated by their selectivity filter.**
Projection from human Nav1.5 `Q14524` positions 372/898/1419/1711, verified
against the reference before every use. Measured 6/6 correct.
`src/classify/motifs.py`.

**D27 — The phylogeny is a forest, and non-alignable comparisons go to the
fold network.** `build_tier2()` raises `NotAlignable`; there is no
`build_tier3()`. A tier-2 tree is labelled a *pore-module tree* in every
output because that is what it is.

**D29 — A benchmark drawn from the catalogue classifies leave-one-out.**
Measured: 93 of the 97 proteins in the S1 control panel are themselves
catalogue exemplars, because both come from the same curated gene lists.
Scoring them against a reference set that contains them returns 100 %
identity to themselves and measures nothing. `classify(...,
leave_one_out=True)` drops the query's own accession, and the overlap is
reported in `summary.json` rather than quietly handled. Any future panel
built from the catalogue inherits this rule.

**D28 — A missing tool disables a test, loudly.** `MafftUnavailable` is
raised, not caught; a missing IQ-TREE writes the alignment and no tree, with
the reason in `TreeRun.note`. No silent fallback to a weaker method, ever —
two incomparable methods in one figure is the failure this prevents.
