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

> **Next session: S4** (S3b needs S4's proteome manifest). S3a completed 2026-09-28 — read `results/census_v3/report.md`; the emergent rows it added (H2/H4/H13 absence rules, the PF00520 superfamily conflicts) are S2 rule fixes that can be done before or alongside S4.
>
> **Previously: S3.** S0, S1 and S2 are complete. Read
> `results/census_v2/report.md` before S3. S2's census is exact on
> enumeration (1,245,200 records, every count reconciled) and deliberately
> honest on calls: 27.8 % family, 32.0 % superfamily-only, 40.2 %
> unassigned, with the reference tier not run (D31).
>
> **What S2 hands S3.** (1) **The superfamily-only calls are the families
> the catalogue declares with identical architecture** — Kv1 vs silent
> modifiers (28,028), ASIC/ENaC/DEG (25,925), ANO channel vs scramblase
> (13,287), P2X metazoan vs non-metazoan (8,295), and the whole Cys-loop
> (70,846) and iGluR (56,369) superfamilies. That is S3's work: one profile
> per family, best-profile assignment with the D7 margin. (2) **Two
> catalogue problems inflate the unassigned 40 %** and are emergent tasks,
> not S3's: co-domain signatures declared `FAMILY`-level, and derived rules
> that require every declared signature. (3) The four-repeat projection
> costs ~3.5 CPU-hours per census pass (37,393 MAFFT alignments); a
> profile-anchored filter read (hmmalign) would make it cheap.
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
| S1 | Classifier control benchmark: positives from every census family, one decoy per hazard; per-family recall, per-hazard specificity, per-tier attribution | S0 | completed 2026-08-19 | **Recall 50/72 (69.4 %); specificity 25/25 (100 %); 16/16 hazards exercised.** 97 proteins classified leave-one-out (**93 of them are catalogue exemplars** → D29) in 866 s. **Per-tier attribution — the headline: 29 calls from architecture/hazard rules and 7 from the filter motif, against 24 from reference identity**, so the classifier is not a nearest-neighbour lookup. **The benchmark found a bug in the inherited identity metric**: covered-only identity without a coverage floor scored connexin-26 at 61.3 % to RYR2 and 49.8 % to connexin-43, and called it a ryanodine receptor. Fixed by requiring 30 % coverage of the longer sequence (**D30**); specificity 24/25 → 25/25. Recall failures concentrate where domain annotation genuinely cannot separate families — the Cys-loop receptors (same two accessions, 35–40 % mutual identity) and ANO1 vs ANO6 (identical architecture) — which is S3's problem, not a tuning problem. → `results/benchmark_controls/report.md` |
| S2 | Uncapped InterPro enumeration of every pore signature → census v2, with a positive family call on every record | S1 | completed 2026-09-28 | **1,245,200 UniProtKB records** (2026_03) carry ≥ 1 of the 67 pore signatures; enumeration exact on every check — **12/12 taxonomic shards, the union and 67/67 signatures fetched == UniProt's own count**, shards sum to the union, 0 duplicates; InterPro counts reported alongside (largest gap `PF00520` +1,213, release cycles) (**D31**). Calls: **346,627 family (27.8 %)** — 281,007 to channel families — **398,282 superfamily-only (32.0 %)**, **500,291 unassigned (40.2 %)**. **Every family call made without the reference tier**: 315,261 architecture/hazard + 31,366 filter motif (37,393 four-repeat records projected; EEEE 13,966 · DEKA 10,744 · EEDD 5,298 · EEKE 1,358). **Human: 319/320 census genes enumerated** (only GLRA4 missing) — the H7 fix works — and **173/320 called to the right family**; 141 stop at the superfamily, 4 unassigned, **1 wrong family (ZACN → AChBP, the known H2 failure)**. **S1 panel: 71/71 channel-family members enumerated; 0/26 non-channel members called to a channel family.** The unassigned 40 % is mostly **three co-domain signatures carried by non-channels** — cNMP_binding 148,190 solo, SBP_bac_3 115,927, PAS_9 111,788 — and **derived rules requiring a family's full declared architecture** (MscS 65,506 partial). Bulk census (TSV 28.6 MB + FASTA 289 MB, SHA-256 in manifest) on the data root. → `results/census_v2/report.md` |
| S3a | **Profile library + best-profile assignment over census v2** → census v3a: one HMM per catalogue family (90) from rule-enforced seeds, benchmarked leave-one-out on the S1 panel and on held-out orthologues, calibrated against S2's calls, then every v2 record assigned with the D7 margin | S1, S2 | completed 2026-09-28 | **Census v3a: 730,790 of 1,245,200 records (58.7 %) carry a family call, against 346,627 (27.8 %) in S2** — both instruments 310,386, S2 only 26,115, profile only 394,289; superfamily-only 93,766, unassigned 407,787, **conflict 12,857 kept and counted**. **90 profiles** (one per catalogue family, controls included) from **834 rule-enforced seeds** (R1 human 455 / R2 exemplars 36 / R3 reviewed S2 calls 343), 71,421 match states; 18 profiles rest on ≤ 2 sequences and are reported as thin. **Benchmarked before use**: S1 panel leave-one-out **69/71** (vs 50/72 in S1) — both misses are sole-seed families (TRPA1 no call, ZAC → 5-HT3); 26/26 decoys not called channel, 22/26 called to their own control family; held-out non-seed orthologues **534/534** (but 576/583 chordate). **Calibrated against S2 (seeds excluded): 96.8 % agreement where both call** (Bacteria 99.2 %, Eukaryota 96.7 %). S2's superfamily-only calls resolved: Cys-loop 93.9 %, TMEM16 98.3 %, P2X 95.3 %, CLC 86.3 %, iGluR 82.3 %, DEG/ENaC 63.8 %, P-loop 61.6 %. **Human 318/320** right family (173 in S2; ZACN now `conflict`, GLRA4 not enumerated). **Conflicts concentrate on S2's absence-test hazard rules**: H2 (AChBP) 6,633 records, median 383 aa — receptor length; H13 (TRPP vs polycystin-1) 2,989, median 2,263 aa. **External check vs parent censuses** (compared, never imported): IP3R v6 same call on **12,204 / 15,601**; **502 swapped ITPR/RYR are one S2 error** — `H4-itpr` calls ITPR on `PF08709` (shared by both) *without* RyR domains, and 1,431 ITPR calls rest on it alone; PIEZO v5 5,055 called piezo, 2,258 unassigned here. Sweep: 1,187,702 unique sequences, 6.9 profile-hours (83 min wall). → `results/census_v3/report.md` |
| S3b | Profile sweep of the S4 reference proteomes (what domain search misses) + jackhmmer-to-convergence completeness argument (D10) → census v3 | S3a, S4 | pending | Split from S3 on 2026-09-28: the brief's sweep needs S4's declared proteome set, and "from its S6 alignment" inverted the dependency order (S6 depends on S3). S3a builds its own seed alignments; S6 upgrades them. |
| S4 | Proteome scope: declared reference-proteome manifest across the lineage panel (the denominator) + download tooling | S1 | pending | |
| S5 | Genomic tblastn + miniprot sweep for families the proteomes miss → per-lineage ledger (found / lost / assembly-gap) → census v4 | S4 | pending | |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl per family, and pore-module extraction for every tier-2 unit | S2, S3a, S3b, S5 | pending | |
| S7 | **Tier-1 phylogenies**: IQ-TREE 2 + ModelFinder + 1000 UFBoot for every census family, rooted on its declared outgroup | S6 | pending | |
| S8 | **Tier-2 pore-module phylogenies** for every alignable superfamily, plus the **tier-3 fold network** (Foldseek/TM-score) for the refused ones | S6, S7 | pending | |
| S9 | Selectivity-filter atlas: the filter locus for every P-loop family, and whether it is congruent with the tier-2 tree (Q5) | S7, S8 | pending | |
| S10 | Repertoire evolution: presence/absence of every family across the lineage panel; ancestral reconstruction; the MscS animal-absence question (Q4) | S5, S8 | pending | |
| S11 | Duplication history: which families expanded in which lineage, 2R/3R ohnologue status, and the Kv/Nav/Cav repeat-duplication order | S7, S10 | pending | |
| S12 | Structures: AFDB coverage per family, TM-align against the experimental reference set, Foldseek all-vs-all for the fold network | S2, S6 | pending | |
| S13 | ML selection tests on the families with clinical variant sets (CFTR, SCN1A, KCNQ1, RYR1) | S6, S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3a, S3b, S5, S7, S8, S9, S10, S12, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Method contribution** — what would domain search alone have missed, per superfamily (Q3); the recall curve as methods are added | S2, S3a, S3b, S5 | high | pending | |
| S16 | **Annotation-quality audit** — how often a real channel locus is missing, fragmentary, split, unnamed or filed under the wrong family across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5, S15 | high | pending | |
| S17 | **Mechanism vs clade** — do the CLC channel/transporter and anoctamin channel/scramblase splits survive the tier-1 trees (Q7) | S7 | high | pending | |
| S18 | **Convergence audit** — connexin vs pannexin, TMEM175 vs the GYG channels, the mechanosensitive families: is "unrelated" a fact or a detection limit (Q6) | S8, S12 | high | pending | |
| S19 | **Channelopathy map** — clinical variants of every family mapped onto the pore module and the constrained core | S6, S12 | medium | pending | |
| S20 | **Auxiliary-subunit census** — the excluded 15 %: how many there are, and how often they are counted as channels in published totals | S2 | medium | pending | |
| S21 | **Prokaryotic and viral channels** — the outgroup census, and whether the eukaryotic families have prokaryotic sisters | S3b, S5 | medium | pending | |
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
| 2026-08-19 | S1 | **The benchmark re-fetches its whole panel on every run.** 97 proteins × 3 API calls, several minutes, repeated for each of the four runs the setup session needed. The parent project cached its control panels to `panel_positives.json` / `panel_decoys.json` so the benchmark reruns offline; do the same here — the evidence (sequence, Pfam counts, TM count) is what the classifier consumes and it does not change between runs | open (S1) |
| 2026-08-19 | S1 | **The reference tier does not scale to a census.** Measured on the S1 run: ~2.5 s per MAFFT pairwise alignment × up to 12 prefiltered candidates ≈ 30 s per protein, so the 97-protein panel takes ~45 min and a 5,000-protein census would take ~40 hours. S2 needs either a calibrated fast identity estimate (DIAMOND or BLAST, checked against MAFFT on the S1 panel so the D7 margin keeps its meaning) or a profile-based assignment that skips pairwise scoring entirely. **This is a scaling result, not a bug** — the tier is correct and slow, and the fix must not change what the margin means | open (S2/S3) |
| 2026-08-19 | setup | **The phylogeny pipeline was smoke-tested on the Cys-loop superfamily** (8 exemplars, MAFFT → trimAl → IQ-TREE 2 with 1000 UFBoot, rooted on GLIC/ELIC, 52 s) and returned the expected topology: anion- and cation-selective receptors separate at 100 %, AChBP sisters the cationic clade. `results/phylogeny/tier2_cysloop/`. **This is not S8's tree** — eight sequences, exemplars only — and must not be cited as one | closed 2026-08-19 |
| 2026-08-19 | S0 | **`PF00005` (ABC_tran) has 1.66 million UniProt proteins and `PF00520` has 206,115.** The census search space is dominated by two accessions that are mostly not channels. S2 needs a per-signature triage rule before it enumerates, or it will fetch a million transporters to find one CFTR | open (S2) |
| 2026-09-28 | S2 | **Three pore-union signatures are co-domains that mostly sit on non-channels.** Measured as the only pore signature on unassigned records: `PF00027` cNMP_binding **148,190**, `PF00497` SBP_bac_3 **115,927**, `PF13426` PAS_9 **111,788** — CRP-type regulators, periplasmic binding proteins and PAS sensors. They are declared `FAMILY`-level for kv_eag/HCN/CNG and iglur_prok, so `pore_signatures()` enumerates from them. Re-level to `ACCESSORY` (or add a pore-module requirement to those families' rules), re-run S1, and the ~376,000 records these three bring in alone leave the census space: 869,000 records remain, of which ~124,000 (~14 %) would still be unassigned. Check against `human_recall.tsv` that no channel is lost before adopting it. Supersedes the PF00005 triage row for S2's purposes | open |
| 2026-09-28 | S2 | **Derived family rules require a family's *complete* declared architecture, and many real members carry part of it.** Unassigned records carrying part of exactly one family's architecture: MscS **65,506** (mostly `PF00924`±`PF21082` without `PF05552`), OSCA 4,505, TRPM 3,960, Piezo 3,939, Slo 2,243, RyR 1,503, Kir 1,081. A "core pore signature suffices" rule per family would call them; it must be benchmarked on S1 (and on a non-human panel) before adoption, not tuned on the census. `results/census_v2/partial_architectures.tsv` | open |
| 2026-09-28 | S2 | **Four-repeat filter strings with no call.** Projected and uncalled: NEEE 1,267, DEEA 435 (the invertebrate Nav2/BSC1 calcium-selective filter, if the literature confirms), QEEE 221, DENA 164, DDDD 129, plus gapped reads. Each needs a literature-backed `FILTER_CALLS` entry or a recorded "no call" — never an inferred one. S9's atlas | open (S9) |
| 2026-09-28 | S2 | **ZACN is still called AChBP**, now in the census (the only wrong-family call among 320 human genes). The H2 rule's premise — an LBD without an annotated TM domain is AChBP — is false for ZACN. Needs a positive AChBP test (secreted, no TM features, `tm_count` = 0) rather than absence of `PF02932` | superseded by the S3a H2 row (profile calls ZACN `zac`; v3a records a conflict) |
| 2026-09-28 | S2 | **`results/session_live.json` is untracked scratch** written by every driver. Gitignored this session | closed 2026-09-28 |
| 2026-09-28 | S3a | **Hazard rule `H4-itpr` is an absence test and calls RyR fragments ITPR.** It fires on `PF08709` *without* `PF02026`/`PF06459`; `PF08709` is shared by ITPR and RyR. The external check found **502** records the IP3R project calls RYR that census v3a calls ITPR — all this rule, profile abstaining, median 76 aa — and **1,431** v3a ITPR calls rest on it alone. Rewrite as a positive test (an ITPR-only domain, or the S3a profile margin), re-run S2 classification on the affected records, re-merge | open |
| 2026-09-28 | S3a | **Hazard rule H2 (AChBP = Cys-loop LBD without annotated TM) is an absence test and fails at scale** — ZACN was the first case. 6,633 S2 AChBP calls are profile-called receptors (nAChR 3,726, GABA-A 1,566, …), median length 383 aa against AChBP's 210–240. Needs a positive AChBP test (length band + no TM features + profile) | open (supersedes the ZACN row) |
| 2026-09-28 | S3a | **H13 (TRPP vs polycystin-1) miscalls polycystin-1-like proteins with few PKD repeats**: 2,989 S2 TRPP calls are profile-called `assoc_polycystin1`, median 2,263 aa (TRPP is 600–900). The `PF00801` copy-number discriminator needs a length/profile positive test | open |
| 2026-09-28 | S3a | **`PF00520` assigns the P-loop superfamily in S2, but Hv1, VSP, ITPR and RyR (other superfamilies) carry it too.** 2,512 v3a conflicts are S2-`ploop` records the profile calls hv1 (1,973, median 252 aa — Hv1 length, probably the non-mammalian Hv1s the mammal-only H9 test misses), itpr 375, ryr 83, vsp 81. Either the signature is declared shared across those superfamilies or S2's superfamily call from it is dropped | open |
| 2026-09-28 | S3a | **`s0_lib.resolve_gene` matches synonyms**: `gene_exact` hits synonyms and the longest entry wins, so TRPC7 → TRPM2, GRIK2 → GRIK5, KCNG3 → KCNG4, GJC3 → GJE1, AQP7 → AQP9, CACNG6 → CACNG8, KCNMB2 → KCNMB3. S3a uses `resolve_primary()`; S0's exemplars and S2's human recall were checked and are unaffected. Fix the shared helper before anything else calls it | open |
| 2026-09-28 | S3a | **Four auxiliary-subunit families pool unrelated proteins** (`assoc_k_beta`: Kvβ aldo-keto reductases + BK β + LRRC γ + DPP6/10; also `assoc_mcu_reg`, `assoc_catsper_aux`, `assoc_clc_aux`). One profile cannot represent them; LOO decoys from them get no hit. Split into homologous families (S20) — no census effect, none carries a pore signature | open (S20) |
| 2026-09-28 | S3a | **The profile-only superfamily splits (Cys-loop, iGluR, DEG/ENaC, P2X, CLC, TMEM16) have no non-vertebrate test.** S2 never called these families, so calibration cannot reach them, and the orthologue benchmark is 576/583 chordate. Build a non-vertebrate labelled panel (e.g. *Drosophila* / *C. elegans* receptors with literature subunit identities) before S7 uses these calls | open |
| 2026-09-28 | S3a | **2,258 PIEZO-census (v5) records are unassigned in v3a.** The parent list includes fragments and short-motif hits; check whether these fall below S3a's coverage gate or are PIEZOs missed | open |

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

**D30 — Covered-only identity needs a coverage floor, measured on the
longer sequence.** Inherited from the parent project, covered-only identity
scores the columns where both sequences have a residue, and on its own it is
unsafe across a catalogue whose members differ 40-fold in length: an aligner
placing 226 residues inside 4,967 selects the 226 best-matching positions,
and the metric then scores exactly those. Measured in S1: **connexin-26
scored 61.3 % identity to ryanodine receptor 2 — higher than the 49.8 % to
its own relative connexin-43** — and GluN1 scored 50.9 % to the SARS-CoV-2
envelope protein. Coverage of the *shorter* sequence does not catch it
(0.996); coverage of the longer does (0.045 against 0.589). The reference
tier therefore drops any hit covering less than
`reference.MIN_COVERAGE` = 30 % of the longer sequence, reports how many it
dropped, and the self-test locks the asymmetry in.

**D29 — A benchmark drawn from the catalogue classifies leave-one-out.**
Measured: 93 of the 97 proteins in the S1 control panel are themselves
catalogue exemplars, because both come from the same curated gene lists.
Scoring them against a reference set that contains them returns 100 %
identity to themselves and measures nothing. `classify(...,
leave_one_out=True)` drops the query's own accession, and the overlap is
reported in `summary.json` rather than quietly handled. Any future panel
built from the catalogue inherits this rule.

**D31 — Census enumeration walks UniProtKB, sharded by taxonomy; InterPro
is the independent count.** S2's brief named InterPro's per-signature
listing, but that listing returns only the queried entry, and the
architecture tier needs each protein's *complete* Pfam list with copy
numbers — the rules consult 19 accessions outside the pore union (forbid
lists, subfamily markers) and `min_copies` needs `PF00520`×4. UniProt's JSON
carries every Pfam cross-reference with its `MatchStatus` copy count, the
transmembrane features and the sequence in one record. The union query
(1,245,200 records, release 2026_03) is cut into taxonomic shards whose
counts are checked to sum to the union before any walk starts, each shard
is cursor-walked to UniProt's own count, and per-signature membership is
checked against UniProt's count for that signature alone; InterPro's count
is reported alongside with the difference (the two are on different release
cycles — `PF00520`: 207,571 vs 206,358). **At census scale the tiers run
selectively, and say so on every row**: architecture, hazard and K⁺-regex on
every record; the four-repeat MAFFT projection only where `PF00520` ≥ 4 (the
H1 trigger, the only architecture where the filter can change the call);
the reference tier **not at all** (~30 s/protein × 1.25 M ≈ 10,000
CPU-hours) — per D28 a missing result with `reference_tier = not_run`, not a
weaker substitute. S3's profile assignment is the replacement.

**D32 — Best-profile assignment: three gates, fixed before the census was
seen.** One profile per catalogue family, controls included (a decoy needs
a profile to win). A record is called to the best-scoring profile only if
(1) it scores ≥ 30 bits, (2) its significant domains cover ≥ 30 % of that
profile's match states — D30's number, so a protein sharing one module
(cNMP, T1, SPRY) is `module`, not a call — and (3) it beats the best other
profile, whatever that profile's coverage, by ≥ 10 % of its own score (D7).
Inside the margin: `superfamily_only` if the band is one superfamily, else
`ambiguous`. The merge with S2 is a stated rule — agree → `both`, one speaks
→ that one, different families (or a profile family outside S2's
superfamily) → `conflict`, kept. A consequence to remember: an abstaining
profile leaves S2's call standing, so an S2 rule error survives the merge
(the H4 finding). Seeds are resolved by *primary* gene name only.

**D28 — A missing tool disables a test, loudly.** `MafftUnavailable` is
raised, not caught; a missing IQ-TREE writes the alignment and no tree, with
the reason in `TreeRun.note`. No silent fallback to a weaker method, ever —
two incomparable methods in one figure is the failure this prevents.
