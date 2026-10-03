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
3b. **Refresh `README.md` and the dashboard with the task's results *and*
   figures** — the README is the state-of-the-project summary for someone
   arriving cold, the dashboard is the live view; both must show what the
   task found, not just that it finished:
   - **Figure.** (Rules: CLAUDE.md § *Figures* — a title on every panel or the
     figure, enforced by `figstyle.save()`.) Every completed task that produces results draws at least
     one headline figure from its committed tables, through
     `scripts/figstyle.py` (`save()` refuses an off-canvas bbox), into
     `results/<task dir>/figures/<name>.png` (+ `.pdf`), by a script
     (`scripts/s<n>_figures.py`, D13 — never hand-made). Look at it before
     committing (D11). A task with genuinely nothing to plot says so in its
     Results cell.
   - **Description.** Add the figure to `scripts/figure_notes.py`: a
     plain-English *what it shows*, *how it was made* and *how to read it*
     (every colour, line and mark explained, no unexplained abbreviations),
     written after looking at the figure. This one entry feeds both the
     README and the dashboard.
   - **README.** Update the task's Status-board row (headline numbers +
     report link) and run `python3 scripts/readme_figures.py`, which
     rewrites the **Results in figures** section from `figure_notes.py`.
   - **Dashboard.** Run `python3 scripts/dashboard.py` (it reads
     `figure_notes.py`) and confirm the task card, Results text, figure and
     description appear.
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

> **Next session: S9** — selectivity-filter atlas: the filter locus for every P-loop family and whether it is congruent with the tier-2 tree (Q5); read `docs/session_briefs.md` § S9 and `results/phylogeny/tier2_report.md` § 6 first. **Caution from S8b**: the P-loop tier-2 tree is a 90-column pore-module tree with 15 % of edges at UFBoot ≥ 95 and no root (the KcsA outgroup is split) — congruence can only be tested against its supported clades, and S9 must say so; S0 already holds `filter_k.tsv` / `filter_four_repeat.tsv` and S6 the extracted modules (K⁺ filter inside the module 1,772/1,834). **S8b done 2026-10-02** (7 tier-2 trees; 3/7 rooted; 34/54 family groups one clade; animal `plgic_prok` group with GLIC/ELIC, not in any Cys-loop family). Previously: **S8b** — the tier-2 trees were launched detached at the end of S8a (`nohup caffeinate -i python3 -u scripts/s8_tier2.py run --jobs 4 --threads 2`, log `<data root>/trees/s8/run.log`). Check they finished (`pgrep -fl iqtree2`; one `run.json` per unit under `<data root>/trees/s8/iqtree/`); if the run died, re-run the same command (skips finished units, resumes from checkpoints). Then write the parse step (outgroup monophyly + root UFBoot as S7, per-family monophyly, support; where the 4 animal `plgic_prok` tips fall), report § 6, tree panels. Do not change D48 after reading a tree. **S8a done 2026-10-02** (D48; fold network 3/5 literature edges supported). Previously: **S8** — tier-2 pore-module phylogenies for every alignable superfamily plus the tier-3 fold network for the refused ones. Read `docs/session_briefs.md` § S8 and `results/phylogeny/tier1_report.md` § 4 first. Choose representatives per clade × kingdom by a rule in a script **before any run** (D8; the S7a compute row: nAChR alone took 26.7 h at tier 1). Inputs: S6's module FASTAs (`<data root>/alignments/s6/`, `results/alignments/modules.tsv`); check the ITPR module-span row (structure-derived span) and the `plgic_prok` animal-members row before the Ca²⁺-release and Cys-loop trees are read. **S7d done 2026-10-02: 0/6 re-rooted P-loop families resolved under D47** — the six roots stay undefined and pass to S11. Previously: **S7d** — the re-rooted trees. **S7c done 2026-10-01** (design + inputs, D47). Previously: **S7c** — re-root the six weakly rooted P-loop families (CNG, K2P, KCNQ, Shaker, Slo, EAG): closer architecture-matched outgroups declared in the catalogue and S0-verified, denser basal sampling from S4b, root acceptance criteria fixed before any tree (outgroup one clade, UFBoot ≥ 95, stable under two outgroup choices). Full brief: `docs/session_briefs.md` § S7c. The re-tree is an overnight detached IQ-TREE run. Then **S8** (choose representatives per clade × kingdom, D8, before any run). **S7b done 2026-10-01: 63 tier-1 trees.** **S4b dense panel done 2026-09-30 (D46): S10 reads repertoire on 439 orders.** **TMEM87 by descent (S2f, D45) done 2026-09-30.** **Catalogue additions (PACC1 + 7 contested, KChIP, decoys), census revision r4 (S2d) and its genome sweep (S5c) done (user-directed, 2026-09-29/30).** **S20 and S15 done out of order (user-directed, 2026-09-29, while S7b's trees run).** **Next session: S7b** — the tier-1 IQ-TREE run was launched detached at the end of S7a (`nohup caffeinate -i python3 -u scripts/s7_trees.py run --jobs 5 --threads 2`, log `<data root>/trees/s7/run.log`, ~15–25 h). Check it has finished (`pgrep -fl iqtree2`; one `run.json` per family under `<data root>/trees/s7/iqtree/`); if it died, re-run the same command — it skips finished families and resumes the rest from checkpoints. Then `s7_trees.py parse`, `s7_report.py`, `s7_figures.py` (piezo1 env), read root support and outgroup monophyly per family. Do not change D41 after reading any tree. S7a completed 2026-09-29. Previously: **S7** — tier-1 phylogenies: IQ-TREE 2 + ModelFinder + 1000 UFBoot on each of S6's 63 trimmed family alignments (`<data root>/alignments/s6/family/<fam>.trimal.fasta`), rooted by adding outgroup sequences with `mafft --add --keeplength`. Read `results/alignments/report.md` first, and **decide the trimming question before building any tree** (S6 emergent row: automated1 keeps 3–4 % of columns on K2P/nAChR/CNG) — measured, not tuned on the trees. S6 completed 2026-09-29 (D39, D40). Previously: **S6** — alignment upgrade: MAFFT L-INS-i + trimAl per census family, and pore-module extraction for every tier-2 unit (`src/phylo/modules.py`), extraction method recorded per sequence. Inputs: census v3 (proteomes) and census v4's genome loci (`<data root>/genomes/s5/census_v4.tsv.gz`), with the S2b emergent row's R3 seed re-draw (deferred by D39: a census revision). Read `results/genome_sweep/report.md` first; before S7 uses any TRPN/TRPA calls, see the S3b ankyrin emergent row. S5b completed 2026-09-29: 52 genomes, matched detection 99.3 %, 72 controlled absences, 11 proteome misses, census v4 (D38). Previously: **S5b** — the full sweep. S5a completed 2026-09-29 (S5 split: S5a instrument + pilot, S5b full sweep). Previously: **S5** — genomic tblastn + miniprot for families S3b finds absent in a lineage, plus the two genome-only species (*Cornu*, *Torpedo*), so an absence passes D4's two bars. S3b completed 2026-09-29: read `results/panel_sweep/report.md`. **What S3b hands S5**: `family_by_species.tsv` (census families × 50 proteome species) is the presence matrix whose zero cells S5 must test; the high-confidence calls are the evidence, medium calls in repeat-dominated profiles (TRPN, LRRC8, TRPA/C/V) are not. **Previously S3b** — sweep the S3a profiles over the S4 panel DB (`<data root>/proteomes/s4/panel_refprot.fasta`, 822,499 sequences) + jackhmmer (D10). S4 completed 2026-09-28: read `results/proteome_scope/report.md` — the denominator is 52 species, 50 reference proteomes + 2 genome-only (D35). **Previously S4** (S3b needs S4's proteome manifest). S3a2 (ABC decoy family) completed 2026-09-28. S2b and S2c (user-directed) completed 2026-09-28: every hazard rule is a positive test (D33), census v2 is r3, v3a re-merged. An emergent row proposes a generic ABC-transporter decoy family for the S3a profile library. S3a completed 2026-09-28 — read `results/census_v3/report.md`; the emergent rows it added (H2/H4/H13 absence rules, the PF00520 superfamily conflicts) are S2 rule fixes that can be done before or alongside S4.
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
| S3a | **Profile library + best-profile assignment over census v2** → census v3a: one HMM per catalogue family (90) from rule-enforced seeds, benchmarked leave-one-out on the S1 panel and on held-out orthologues, calibrated against S2's calls, then every v2 record assigned with the D7 margin | S1, S2 | completed 2026-09-28 | **Census v3a: 730,790 of 1,245,200 records (58.7 %) carry a family call, against 346,627 (27.8 %) in S2** — both instruments 310,386, S2 only 26,115, profile only 394,289; superfamily-only 93,766, unassigned 407,787, **conflict 12,857 kept and counted**. **90 profiles** (one per catalogue family, controls included) from **834 rule-enforced seeds** (R1 human 455 / R2 exemplars 36 / R3 reviewed S2 calls 343), 71,421 match states; 18 profiles rest on ≤ 2 sequences and are reported as thin. **Benchmarked before use**: S1 panel leave-one-out **69/71** (vs 50/72 in S1) — both misses are sole-seed families (TRPA1 no call, ZAC → 5-HT3); 26/26 decoys not called channel, 22/26 called to their own control family; held-out non-seed orthologues **534/534** (but 576/583 chordate). **Calibrated against S2 (seeds excluded): 96.8 % agreement where both call** (Bacteria 99.2 %, Eukaryota 96.7 %). S2's superfamily-only calls resolved: Cys-loop 93.9 %, TMEM16 98.3 %, P2X 95.3 %, CLC 86.3 %, iGluR 82.3 %, DEG/ENaC 63.8 %, P-loop 61.6 %. **Human 318/320** right family (173 in S2; ZACN now `conflict`, GLRA4 not enumerated). **Conflicts concentrate on S2's absence-test hazard rules**: H2 (AChBP) 6,633 records, median 383 aa — receptor length; H13 (TRPP vs polycystin-1) 2,989, median 2,263 aa. **External check vs parent censuses** (compared, never imported): IP3R v6 same call on **12,204 / 15,601**; **502 swapped ITPR/RYR are one S2 error** — `H4-itpr` calls ITPR on `PF08709` (shared by both) *without* RyR domains, and 1,431 ITPR calls rest on it alone; PIEZO v5 5,055 called piezo, 2,258 unassigned here. Sweep: 1,187,702 unique sequences, 6.9 profile-hours (83 min wall). → `results/census_v3/report.md` **Re-merged after S2b** (census v2 r2): conflicts 3,288, human 319/320, 0 ITPR/RYR swaps; the numbers above are the first merge. Profiles unchanged (SHA-256 in `profile_build.tsv`): their R3 seeds were drawn from S2 r1 calls. |
| S3b | Profile sweep of the S4 reference proteomes (what domain search misses) + jackhmmer-to-convergence completeness argument (D10) → census v3 | S3a, S4 | completed 2026-09-29 | **Census v3 over the declared denominator: 26,936 of 822,499 panel entries carry evidence; 8,755 census-family calls, of which 230 (2.6 %) are outside census v2 — only 96 at high confidence.** (The entry count was corrected in S5b from a mis-transcribed 49,595; the committed tables always said 26,936.) S3a's 91 profiles unchanged, swept over S4's `panel_refprot.fasta` (50 proteomes, `-Z` 822,499; 12.4 min wall), D32 assignment imported: 16,183 family calls (11,798 high). **Instrument check: same verdict as S3a on 12,072 / 12,096 (99.8 %)** of the panel entries that are also census v2 records (only `-Z` differs). **Human 319/320 right family** (GLRA4 not in the proteome). Where domain search misses channels (high-confidence): **CLIC 27**, **Hv1 15** (molluscs, moss, fish — 25/33 panel Hv1 carry no enumerated signature), CNG 10, TRPV 8, pannexin 8, and **all 3 viroporins** (Vpu, M2, SARS-CoV-2 E — the S4 emergent row confirmed). By group the rate runs from 1.0 % (vertebrate) to 12.8 % (plant). **The medium 134 are an upper bound, not a finding**: TRPN's 37 medium calls are ankyrin-repeat proteins (ANKRD52, ankyrin, IκB) and LRRC8's 20 are LRR proteins — repeat-dominated profiles let the module alone pass D32's 30 % coverage gate (emergent). **jackhmmer: 68/68 families seeded by rule** (top high-confidence panel call), `-N 10 -E 1e-5`, 26.9 run-hours (6.1 h wall). **D10: 24 clean, 44 killed (K1 23, K2 2, K3 19).** K1 fires at round 2 across P-loop/iGluR/TRP: a family-seeded model takes in its sister families in one iteration, while every Cys-loop, DEG/ENaC and P2X family run converges clean on the *same* superfamily set (1,625 / 490 / 125) — at those scales jackhmmer is a superfamily instrument. **Completeness (clean runs): 3,352 / 3,356 (99.9 %) of profile calls recovered**, 1,460 accepted targets no profile calls; all runs 8,610 / 8,736 (98.6 %). K3 runs contribute no candidates (**D36**). 7,980 `jackhmmer_only` candidates, dominated by the TRPA/TRPN pre-K1 rounds (ankyrin drift K1 cannot see) — counted, never called. Bulk on the data root under `hmmer/s3b/`. → `results/panel_sweep/report.md` |
| S2b | **Hazard-rule rewrite (H2, H4, H13) as positive tests** — the three S2 rules that call a family from an *absence*, found by S3a's calibration and external check; re-benchmark S1, re-classify the affected S2 records, rebuild census v2 and re-merge v3a | S3a | completed 2026-09-28 | **All three rules now call a family only on a domain it carries and its rivals lack; where none exists the architecture tier stops at the superfamily (D33).** H4: `PF02026` or `PF06459` with the shared core ⇒ RYR; **ITPR has no positive architectural test** (no domain RyR lacks) → superfamily, family from reference/profile. H13: PLAT/REJ/GPS or `PF00801`×≥5 ⇒ polycystin-1; `PF18109` ⇒ TRPP (measured: 1,358 records, all profile-TRPP); channel domain alone ⇒ superfamily; catalogue fix `trpp` `PF20519` FAMILY → SHARED_WITH_DECOY. H2: **no positive AChBP test exists at this tier** — the best candidate (complete, soluble, AChBP-sized LBD) agreed with the profiles on 53/1,263 and put ~75 % of its 2,164 calls in Ecdysozoa/Chordata, which have no known AChBP → LBD without TM is superfamily-only. **Re-classified 108,407 records** (`s2_classify.py --recheck`), **27,028 calls changed, 0 outside the re-checked set**; census v2 is now **r2** (r1 archived). S2 r2: family 324,769 (26.1 %), human 169/320 right (ITPR1–3, PKD2L2 now superfamily-only) and **0 wrong** (ZACN no longer AChBP). **S1 benchmark unchanged**: recall 50/72, specificity 25/25 — ITPR1 now called by reference identity, PKD1 now positively by PLAT/REJ. **Census v3a re-merged: conflicts 12,857 → 3,288; ITPR/RYR swaps vs the IP3R census 502 → 0**; family calls 732,191 (58.8 %); human **319/320** (ZACN now `zac` by profile). **Two bugs found on the way**: `tm_count` 0 was passed as unknown (`x if x else None`) in both `s2_classify.py` and `cli_channel.py`; UniProt's `flag` "Precursor" was read as a fragment. Both fixed, both self-tested. → `results/census_v2/report.md` § Revision r2, `results/census_v2/s2b_transitions.tsv` |
| S2c | **Hazard-rule rewrite (H11, H12) as positive tests** (D33) — SUR vs CFTR, KCTD vs Kv; re-benchmark S1, re-classify affected S2 records (census v2 r3), re-merge v3a | S2b | completed 2026-09-28 | **No hazard rule now calls a family from an absence** (self-test invariant over every hazard). **H12**: KCTD called on one of five KCTD C-terminal domains (catalogued as `nonchannel_kctd` FAMILY signatures; measured on 14,354 records, 14,317 profile-KCTD, 3 otherwise); T1 without a pore module is superfamily-only — the shape test (T1, 0 TM, complete) was rejected, contradicted by the profiles on 160 records. **H11**: CFTR still called on its R domain; **the SUR rule is removed — no domain identifies SUR** (ABCC8 carries only generic ABC domains; TMD0 is on the MRPs too), and none of its 980 census calls was a SUR (918 bacterial cNMP + C39 peptidase ABC transporters, 62 eukaryotic fused gene models); `assoc_sur` ABC domains re-levelled SHARED_WITH_DECOY. **Re-classified 63,411 records, 18,412 calls changed, 0 outside the re-checked set** (KCTD → [ploop] 17,423; SUR → unassigned 944). Census v2 **r3**: family 306,378, human 169/320 right, 0 wrong. **S1 unchanged** (50/72, 25/25, 16/16) — ABCC8 now uncalled, still rejected; KCTD1 still KCTD, by its C-terminal domain. **Census v3a re-merged: conflicts 3,288 → 2,907** (none from any hazard rule's absence logic; what remains is the PF00520 superfamily problem and K2P copy-number calls); KCTD agreement 99.98 %, CFTR 100 %; human 319/320. **Found on the way**: S2b had silently dropped PF20519 from the census search space (enumeration coupled to evidence level) — fixed (D34) and pinned by a self-test; and the S3a profile library's missing ABC-transporter decoy lets the SUR profile call 525 records (508 bacterial) — emergent, a catalogue addition. → `results/census_v2/report.md` § Revision r3 |
| S3a2 | **ABC-transporter decoy family** for the S3a profile library (user-directed, from the S2c emergent row): catalogue `nonchannel_abc_transporter`, its profile built and swept alone (the other 90 frozen), S0 verification, benchmark re-run, v3a re-merged | S2c | completed 2026-09-28 | **The SUR profile no longer wins anything that is not a SUR: profile-called SUR records 525 → 0.** New control family `nonchannel_abc_transporter` (91 families, 23 controls; census families, human genes and the 67-signature search space unchanged): human ABCC1–6/10–12 + *E. coli* HlyB, *B. subtilis* SunT (C39-peptidase exporters), *E. coli* MsbA; ABC domains SHARED_WITH_DECOY, `PF03412` ACCESSORY. **S0 re-run clean** — 793 requests, 0 failures, 166/166 exemplars, 121/121 Pfam, 52/52 taxa; **no existing exemplar sequence changed**, so the 90 frozen profiles still match their seeds. Profile: 12 seeds → 1,304 states, built alone (`--only`); the other 90 builds and seed sets byte-identical. Swept in 21 s (2,919 targets). **S3a benchmark unchanged** (LOO 69/71, orthologues 534/534, 0 decoys called channel); ABCC8 still SUR by 2,478 vs 1,147 bits, CFTR by 3,237 vs 971. **Census v3a**: the decoy takes 1,026 records (997 bacterial); conflicts 2,919 (+12, all ABC-domain/channel-domain fusions, S2 calling the channel domain and the profile the ABC part); human 319/320. → `results/census_v3/report.md` |
| S2d | **Census revision r4 for the families added after S20** (user-directed): enumerate their signatures as a delta, build and sweep their profiles (the other 91 frozen), benchmark the look-alikes, re-merge census v2 → v3a → v3 → v4 | S20 + catalogue additions | completed 2026-09-29 | **Enumeration exact (D43)**: 8 signatures walked as a delta (new signature, none of r3's 67) in the 12 taxonomic shards — **26,783 records; r3's 1,245,200 + delta = UniProt's 1,271,983** for the 75-signature union (release 2026_03, unchanged); `s2_enumerate.py verify` re-counted all 24 shards and 75 signatures, 0 problems, 0 duplicates. 11 r3 records carry a new signature and were re-classified (`--recheck --archive r3`): 2 calls changed, 0 outside that set. Census v2 **r4**: human 328 genes, 326 enumerated, 173 right by domain rules (+ PACC1, CLCC1, GPHRA/B); TMCO1/TMEM87A/TMEM109 unassigned by design (shared domains, H17–H19), MITOK not enumerable (no Pfam). **12 profiles** built alone (7 census families + KChIP + EMC3/GOST/BRI3BP decoys + a **neuronal-calcium-sensor decoy added mid-revision**: the KChIP profile alone called recoverin, hippocalcin, NCS-1 and the GCAPs; with the decoy, human KChIP calls = exactly KCNIP1–4, H20); **the other 91 byte-identical**. New profiles swept over r3's NR and all 103 over the delta NR with **`-Z` held at r3's size** — no r3 hit moved. **11 of 1,245,200 r3 v3a calls changed**, each with a new family as call, winner or margin-closing runner-up (a hard-failure check that fired once, correctly, before the runner-up case was added). New records: TMEM87 5,484 · GOST 4,441 · EMC3 4,280 · GPHR 3,903 · TMCO1 2,237 · CLCC1 1,562 · PACC1 992 · BRI3BP 804 · TMEM109 763. **Benchmark**: S3a set A unchanged; set B 542/542 orthologues (+8 TMCO1), 57 decoys (+7 EMC3, +1 BRI3BP), 0 called a channel — **H17 and H19 pass; H18 only partly**: mammalian TMEM87A → tmem87 high, TMEM87B → GOST at medium (margin 0.12–0.27), yeast PTM1 / *S. pombe* GOST → tmem87 — the TMEM87 count is an upper bound. **Panel census v3**: human 327/328 right family (all 8 new genes, MITOK by profile alone); 12,402 panel entries in census v2, instrument check 12,378/12,402. **Census v4**: 28,891 proteome rows + 434 genome loci (genome sweep not run for the new families). S15 re-run: headline unchanged (94.1 % enumerated, 44.1 % named); MITOK profile-only. S0 clean. → `results/census_v3/r4_report.md`, `results/census_v2/report.md` § Revision r4 |
| S2e | **Census revision r5: the viroporin signatures** (user-directed, from the open items) | S2d | completed 2026-09-30 | The viroporins' four Pfam signatures (Flu M2, Vpu, CoV E, bCoV viroporin), SUBFAMILY-level and so never enumerated, brought into the search space by the r4 delta tooling generalised to `--rev`: **+4,007 records → census v2 1,275,990**, every count exact; 0 earlier calls changed (S2 and v3a); 3,560 profile-called viroporins; panel census and census v4 re-merged; S15's viroporin row now 3/3 enumerated. → `results/census_v2/report.md` § Revision r5, `results/census_v3/r5_*.tsv` |
| S2f | **TMEM87 by descent** (user decision on H18): TMEM87A + TMEM87B one family; census revision r6 | S2e | completed 2026-09-30 | `tmem87` now holds TMEM87A and TMEM87B (human census genes 328 → **329**; the TMEM87B decoy retired, its files kept aside). **H18 became a positive test**: the TMEM87 GOLD domain `PF21901` — all 1,660 UniProt carriers metazoan TMEM87s, none a GOST protein — is family-level, so S2 calls TMEM87 on it (r6: 1,660 records unassigned → `tmem87`, 0 changed outside, all 36 shards checked; S2 human right-family 173 → 175). Profiles rebuilt (TMEM87 5 seeds; GOST decoy 6, now with *Arabidopsis* CAND6 and yeast PTM1, which lack PF21901); v3a: 5,068 calls changed, all within TMEM87/GOST (gated); **benchmark set B 744/744, 0 decoys called a channel** (the *Xenopus* TMEM87 miss gone); panel human 328/329; S5c re-called, **0 of 3,536 old cells changed**; **TMEM87 now present in every animal genome** (was scattered). Outside animals the GOST family still splits between the TMEM87 and GOST profiles (plants, *S. pombe*, *Dictyostelium*, holozoans 'present'; yeast 'absent') — none carries PF21901, and whether any is a TMEM87 orthologue is a phylogenetic question (emergent row). → `results/census_v2/report.md` § Revision r6, `results/census_v3/r6_transitions.tsv` |
| S4 | Proteome scope: declared reference-proteome manifest across the lineage panel (the denominator) + download tooling | S1 | completed 2026-09-28 | **The denominator: 52 panel species in 15 groups → 50 UniProt reference proteomes (47 cellular, 3 viral) + 2 genome-only (*Cornu aspersum*, *Torpedo marmorata*: chromosome-level, unannotated NCBI assemblies), 0 with neither** — release **2026_03**, the one census v2 was enumerated from (D35). **822,499 canonical entries, one per gene** (entry count = declared gene count 50/50). Selection by a stated rule (most Swiss-Prot entries → BUSCO → genes → UPID), which picked the model strain in all 6 multi-candidate species (HB8, S288C, Nipponbare, 3D7, PR8, HXB2). **Verified: 49/50 exact** against the release README *and* the per-proteome metalink MD5; **human was reissued by UniProt mid-release** (served 2026-09-15; README #(1) 147,503 = the whole human UniProtKB set, served file 20,652 = one per gene, MD5 matches) — accepted under D35. **Positive control 319/319** enumerated human census genes are human-proteome entries. 12,096 census v2 records are panel-proteome entries. Per row (D9): assembly accession + level + N50, annotation source, BUSCO, CPD; **7 proteomes sit on a superseded assembly version**, 14 below 90 % BUSCO, 14 scaffold-level — reported, not applied (S5's D4 bars decide). **Found on the way: the viroporin family is outside census v2's search space** (its 4 signatures are `SUBFAMILY`, so never enumerated; 0 census records in the 3 viral proteomes) — emergent. Sweep DB 507 MB on the data root (the ~25 GB estimate assumed a much wider set). Offline rerun reproduces all tables byte-identically. → `results/proteome_scope/report.md` |
| S4b | **Dense panel for S10**: one reference proteome per eukaryotic order (user decision, D46), the S3a profiles swept over it, a proteome-level family × order matrix | S4, S3b, S2f | completed 2026-09-30 | **439 orders, one release-2026_03 reference proteome each** (Metazoa 194, Fungi 126, Viridiplantae 57, protists 62), 35 keeping their S4-panel proteome; 86 proteomes have no order rank (stated gap). 439/439 files MD5-verified, **8,759,024 sequences** (entry counts = README for 438/439; human = the D35 reissue). The 103 S3a profiles swept (D32, `-Z` = DB size): 281,479 entries with a hit. **Instrument check: 2,624 / 2,625 cells agree with S3b** on the 35 shared proteomes (the one: *S. pombe* anoctamin scramblase, just below high confidence at the larger `-Z`). **Matrix: 12,060 present cells over 439 orders × 75 census families; 74 families present somewhere** (viroporin is viral). Widest: VDAC 422 orders, OSCA 402, GPHR 382, ANO scramblase 324, TMEM87 313; narrowest: GluR0 4, invertebrate degenerins 15, ZAC 17. Every cell a proteome call; absences annotation-level only (D46). → `results/panel_density/report.md` |
| S5a | **Genomic sweep: bait panel, instrument, pilot** — the instrument for S5, measured rather than ported: bait rules, locus → profile call (D32), tblastn rescue, matched positive control, D4 bar from annotation, seven-genome pilot | S3b, S4 | completed 2026-09-29 | **The instrument is measured, and its sensitivity is a property of the nearest bait.** Bait panel: **1,808 baits for all 91 catalogue families** (controls included) from 50 panel species, 1.63 M residues, by rules B1–B3 (`s5_baits.py`). Pilot: **7 genomes, 6.02 Gbp, 0 failures** (*E. coli*, yeast, *Arabidopsis*, *C. elegans*, *Drosophila*, mouse, genome-only *Cornu*), 1,881 loci, **1,020 called to a family by the S3a profiles** (D32 on miniprot translations, `-Z` = panel size; the bait family is a second witness), 29.8 min wall. **The headline is D37.** With 8 baits per family spread one per group, *Arabidopsis* recovered **5 of 9** of its own proteome's high-confidence families once its own baits were excluded — **AtTPC1 drew 0 miniprot alignments and 0 tblastn HSPs (E ≤ 1e-5) from 7 animal/ciliate TPC baits**. Neither pairwise method reaches a relative from another group at 20–30 % identity. One bait per species: 8 of 9. So the control is **matched**: detected / control cells that have a non-self bait from the same group — the condition every informative zero cell is judged under. **Matched detection 172 / 173** (the miss: *C. elegans* tweety); found (a profile call) 166 / 175; every genome ≥ 0.974 against a 0.90 floor. *Cornu*'s control is its group's 27 core families: 27 / 27. **Calibrated from annotation, not miniprot:** 420 called loci confirmed on annotated exons. **Mouse census genes carry introns far past miniprot's 200 kb default — 13 / 220, up to 996 kb (Asic2), Trpm3 558 kb, Kcnd2 497 kb** — so `-G` 1 Mb for genomes ≥ 1 Gbp is necessary and Asic2 sits 4 kb under it (emergent). D4 span bar per family = median annotated span, **61 families measured**. **No identity floor**: annotated loci run down to 0.21 identity, so the parent's 0.40 does not transfer. **Three instrument corrections found by the pilot**: "overlaps an annotated gene" confirmed a 400 kb chained ZAC model in mouse and a TRPM model on Ugt1a10 → confirmation is now CDS-on-exons ≥ 50 % (a reciprocal-extent rule was tried and rejected: mouse 223 → 130, real genes with long UTR introns); mouse carries **12 high-confidence VDAC loci for 3 genes**, most with frameshifts/stops → a proteome-miss claim needs a high-confidence **intact** locus (`genome_weak` otherwise); genome-only presence is `genome_present`, not a miss. **Zero cells (312 in the pilot; 71 informative, 46 of them *Cornu*)**: of the other 25 informative, **10 absent**, 5 trace, 3 absent with no measured bar, 5 partial, 1 gap, 1 weak — and **0 `genome_found`: no pilot proteome misses a high-confidence intact channel gene**. **Literature checks 4 / 4 consistent** (refs pending): *C. elegans* Nav and P2X, *Drosophila* P2X `absent`; mouse ZAC `genome_weak` (the 400 kb chained model, 32 % identity, off the syntenic chromosome). Self-test +8 invariants, the exon rule mutation-tested. Budget for S5b: ~225–260 s/Gbp miniprot on 1,808 baits → ~3 h over 42 Gbp. → `results/genome_sweep/report.md` |
| S5b | **The full panel sweep** (52 genomes incl. *Cornu*, *Torpedo*) → per-cell ledger (found / partial / gap / trace / absent) → census v4 | S5a | completed 2026-09-29 | **52 / 52 genomes swept, 42.5 Gbp, 0 failures**; 11,878 loci, **8,048 called to a family by the S3a profiles**; miniprot 1.94 CPU-h. **Two thresholds fixed before any cell was read (D38)**: `-G` stays 1 Mb — a 1.5 Mb margin *measured* on mouse lost 9 high-confidence loci (Hvcn1 among them) to chaining and gained no call needing > 1 Mb; D4 bar per cell = group median span (≥ 3 genes) → pooled median → 3 × band for prokaryote/virus → none. **Matched detection 1,317 / 1,326 (99.3 %)**; found 1,311 / 1,420 control cells. Below the 0.90 floor: only the 3 viruses (0/1 each: the viroporins are mutually unrelated, so no bait reaches another). No matched control at all: the 4 single-species groups (*Chlamydomonas*, *Dictyostelium*, *Plasmodium*, *Trypanosoma*) — no absence readable there. **3,536 cells; 2,056 zero cells, 247 informative**: **72 `absent`**, 102 `genome_present` (genome-only *Cornu* 48 families, *Torpedo* 56), 17 trace, 24 partial, 11 gap, 8 weak, 4 unmatched, **9 informative `genome_found`** (11 in all). **Proteome misses (high-confidence intact loci)**: *Takifugu* **RyR — 6 loci, 61–80 % identity, full profile coverage, 0 RyR in its reference proteome**; *Ciona* NALCN; *Takifugu* TMEM175; *Daphnia* HCN; *Monosiga* a second, non-metazoan-type P2X beside the proteome's one; *Lottia* TRPA; TRPN (NOMPC) in *Trichoplax*, *Nematostella*, *Daphnia*, *Lottia* (0.79–0.82 profile coverage, 1,350–2,000 bits — not ankyrin-only); *Nematostella* MscS (on 6.7 kb / 3.1 kb contigs — contamination suspected, emergent). **Absences consistent with the literature (refs pending)**: ZAC in rat (and 8 other non-mammalian/non-therian vertebrates), mouse ZAC `genome_weak`; Nav *C. elegans* + *Amphimedon*; P2X *C. elegans*, *Drosophila*; ENaC zebrafish + *Takifugu*; Hv1 *C. elegans*, *Drosophila*. **Calibration from annotation**: 4,328 annotated loci, D4 spans for 65 families; **5 annotated genes have an intron beyond `-G`** (human ASIC2 1,043,910 bp; *Callorhinchus* asic2, trpm3; amphioxus ASIC; *Aplysia* TRPC) — **every one still called**, as D38 predicts. **Two calibration bugs found and fixed**: GCF assemblies' sequence names were mapped to GenBank names (6 genomes had 0 annotated loci), and CDS-only gene structures (*Lymnaea*, *Paramecium*) read as intronless. **Census v4** = census v3's 26,936 proteome rows unchanged + **434 genome loci** (420 genome-only, 14 proteome misses) in 115 cells; bulk + SHA-256 on the data root (`genomes/s5/census_v4.tsv.gz`). → `results/genome_sweep/report.md` |
| S5c | **Genome sweep for the r4 families** (user-directed): the 12 families census revision r4 added, over all 52 genomes, as a delta | S5b, S2d | completed 2026-09-30 | **287 baits** (rules B1–B3 restricted to the 12 r4 families, a panel of their own — S5b's panel, runs and loci byte-identical, SHA checked); miniprot over **52/52 genomes, 0 failures, 50 min**; **1,267 loci, 774 called to an r4 family** (profiles over all 103). An r4 locus counts only toward an r4 family, and an old family's tblastn trace is judged against S5b's loci only (a first run leaked one Monosiga VDAC trace count through the r4 loci — found by the cell-by-cell check, fixed). **All 3,536 S5b cells unchanged, 0 verdicts changed**; 41 genomes' control sets grew by the new families' cells, min matched detection (cellular) 0.944 (*Ciona*, was 0.941). Old families' span bars unchanged; new families' bars measured from annotation (7 families). **The 364 new-family cells**: 167 in the proteome, **3 proteome misses** (*Trichoplax* TMEM87, *Daphnia* GPHR, *Xenopus* MITOK — no near call in their proteomes), 11 genome-only presences (*Cornu*, *Torpedo*), **11 controlled absences** (CLCC1 in *C. elegans*, *Drosophila*, *Daphnia*; MITOK in *Drosophila*, *Lottia*; TMEM109 in lamprey and *Callorhinchus*; PACC1 in *Trichoplax* and *Takifugu*; TMCO1 in *Monosiga*; GPHR in *S. pombe*), 30 partial/gap, 2 weak, 1 trace, 139 not informative. TMEM87 'presence' in fungi, plants and *Dictyostelium* is H18's GOST drift, not an established range. Census v4: 452 genome loci (+18). → `results/genome_sweep/r4_report.md` |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl per family, and pore-module extraction for every tier-2 unit | S2, S3a, S3b, S5b | completed 2026-09-29 | **63 family alignments, one method (MAFFT L-INS-i → trimAl `-automated1`), and 5,228 pore modules cut by one method per tier-2 unit, validated against UniProt at median Jaccard 0.97.** **Sets (D39, fixed before any alignment):** high-confidence S3a profile calls + intact genome loci → **6,658 sequences in 68 families** (323 genome loci); excluded and counted per row: 1,697 medium-confidence, 91 S2-only, 71 broken genome frames, 9 identical. **Alignments:** 63 families aligned, 5 with < 4 sequences (delta-GluR, GluR0, tweety, viroporin, ZAC); every run checked for row count, raggedness and unchanged residues, SHA-256 recorded; L-INS-i 4.4 h summed (nAChR 726 seqs 1.8 h), ~2 h wall. **trimAl keeps a median 22 % of columns but 3–4 % on the largest divergent families** (K2P 145 of 5,149; nAChR 208 of 8,159; CNG 365 of 8,643) — emergent, S7. **Modules (D40):** `Superfamily.module_rule` declared for the five module units (P-loop, iGluR, Ca²⁺-release `pore_loop`; innexin clan, Hv `tm_span`). Every member `hmmalign`ed to its own family profile and cut at the family's span of match states (`profile_projection`). Spans from UniProt topology of annotated references in **29 families**; the **6 with no annotated pore loop** (CNG, TRPC, TRPN, ITPR, plant + prokaryotic iGluR) by a vote of a unit module HMM, **measured held-out on the annotated P-loop/iGluR families: ≤ 12 profile states from the annotated span (median 4)**; a helix-snap rule moved ITPR's RyR-seeded vote (which began in the luminal loop) out to TM5. **A single module HMM as the extractor was measured first and rejected** (held-out Kir 22/298, TRPM 0/145, innexin clan ≤ 5/82). **Validation on 441 reviewed members never used as references: median Jaccard 0.97, 436/441 ≥ 0.8.** Expected module counts by construction (4 per Nav/Cav/NALCN chain, 2 per K2P/TPC); **K⁺ filter inside the extracted module 1,772/1,834 (96.6 %)** — of 62 misses, 45 chains carry no canonical filter anywhere. Self-test +5 invariants. Bulk under `<data root>/alignments/s6/`. → `results/alignments/report.md` |
| S7a | **Tier-1 design + inputs**: trimming measured and fixed, catalogue rooting, outgroups added, IQ-TREE driver (S7 split: the trees need ~15–25 h wall, beyond one session) | S6 | completed 2026-09-29 | **Trimming fixed before any tree (D41): trimAl `-gt 0.5`** — median 504 informative sites/family against 324 for S6's automated1, ≥ automated1 in 62/63 (GlyR 281 vs 325), median 91 % of each member's residues kept (automated1 61 %), gap fraction 7.3 %; K2P 145 → 302 columns, nAChR 208 → 440, CNG 365 → 768 (`tier1_trim_compare.tsv`). **Rooting from the catalogue only**: 36 families rooted on their superfamily's `root_with` exemplars (added by L-INS-i `--add --keeplength`, ingroup rows checked unchanged, same column mask); 5 are the superfamily outgroup, 22 in superfamilies with no declared outgroup — built unrooted, never midpoint. **P-loop outgroups (KcsA, MthK, NaK) fill only 5–16 % of Nav/Cav/NALCN/TRPM/TRPN columns** (0.95 for ASIC on degenerins). **Model search measured and restricted**: full ModelFinder ~500 models at 20–40 s each on the largest families → `-mset LG,WAG,JTT,Q.pfam` for every family; `-B 1000 -bnni -seed 1`. The first 2 h run finished none of the five largest (connexin: JTT+F+R9, tree search at iteration 31) — the driver resumes from IQ-TREE checkpoints; full run launched detached. `s7_trim.py`, `s7_trees.py` (prep/run/parse), `s7_newick.py`, `s7_report.py`, `s7_figures.py`; self-test +4. → `results/phylogeny/tier1_report.md` |
| S7b | **Tier-1 trees**: finish the detached IQ-TREE run, `s7_trees.py parse`, report + figure (support, root support, outgroup monophyly, RyR example tree) | S7a | completed 2026-10-01 | **63 tier-1 trees, 0 failures** (detached IQ-TREE run, 147.8 CPU-h summed; nAChR alone 26.7 h; `-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -bnni`, D41); models Q.pfam 23, JTT 22, LG 18. **Rooting**: 36 rooted on the catalogue's outgroup — 6 on a single outgroup sequence (the iGluR families on GluR0, P2X on *Dictyostelium* P2X, CLC on ClC-ec1: terminal root edge, no support value); of the 30 with several, **the outgroup is one clade in 26, root UFBoot ≥ 95 in 22**. **The four failures are all P-loop families rooted on KcsA/MthK/NaK** — CNG, K2P, Kv (KCNQ), Kv (Shaker) — and Slo (32) and EAG (70) are weakly rooted: S7a's emergent row confirmed (two-helix outgroups filling 15–42 % of columns). 27 unrooted by declaration. **Support**: median 54 % of a family's internal edges at UFBoot ≥ 95; lowest the invertebrate degenerins (22 %). Treefiles `results/phylogeny/tier1/`. → `results/phylogeny/tier1_report.md` |
| S7c | **Re-rooting design + inputs** (S7c split: the 12 trees need ~12–15 h wall): two outgroups per family declared in the catalogue and S0-verified, denser basal sampling from S4b, re-alignment, one input per family × outgroup, detached IQ-TREE run launched | S7b, S4b | completed 2026-10-01 | **Root criterion fixed before any tree (D47)**: in both trees the outgroup is one clade with root UFBoot ≥ 95, **and** the ingroup's basal split is identical under the two outgroups — else `unresolved`. **Catalogue**: new `RootSet` / `ChannelFamily.root_with` (validate: exactly two sets, never the family itself, inside its superfamily); prokaryotic sets **KvAP (Q9YDF8) + MVP (Q57603) + *E. coli* Kch (P31069)** for Shaker, KCNQ, EAG, Slo and **MloK1 (Q98GN8) + SthK (E0RR11)** for CNG; K2P's prokaryotic choice is the superfamily's KcsA/MthK/NaK (no prokaryotic K2P exists); sister sets: Shaker ← KCNQ+EAG, KCNQ ← Shaker+EAG, EAG ← HCN+CNG, Slo ← Shaker+KCNQ, CNG ← HCN+EAG, K2P ← Kir+Shaker. `s0_outgroups.py`: **5/5 inline exemplars verified** (accession, gene, species, pore model), kept in `outgroup_panel.fasta`, apart from the classifier's reference panel; full S0 re-run clean (196/196 exemplars, 132/132 Pfam, 52/52 taxa, 904 requests, 0 failures; reference panel byte-identical). **Basal sampling** (`s7c_basal.py`, one per kingdom × phylum × class outside Bilateria, top-scoring high-confidence complete in-band S4b call): **107 added** (CNG 23, K2P 28, Slo 39, EAG 10, KCNQ 2, Shaker 5; 13 clades already in the D39 sets). L-INS-i re-alignment (CNG 54 min). **Outgroup column occupancy far above S7b's KcsA**: prokaryotic 0.20–0.67 (MloK1 on CNG 0.13), sister families 0.33–0.87. Self-test +3 invariants (root-set rule, root split, criterion). Inputs → `results/phylogeny/tier1_reroot_inputs.tsv`, `reroot_basal.tsv`; report § 4. No figure yet (panel C of the S7 figure gains the new roots in S7d). |
| S7d | **Re-rooted trees**: finish the detached run (`nohup caffeinate -i python3 -u scripts/s7c_reroot.py run --jobs 5 --threads 2`, log `<data root>/trees/s7c/run.log`), `s7c_reroot.py parse`, `s7_report.py`, `bin/envpy scripts/s7_figures.py`; name the families still unresolved | S7c | completed 2026-10-02 | **12 / 12 re-rooted trees built, 0 failures** (87.9 run-hours summed at 2 threads each; CNG 15.6 + 16.4 h). **D47 read as written: 0 / 6 families resolved — CNG, K2P, KCNQ, Shaker, Slo and EAG keep no root.** (a) The outgroup is one clade in only **4 of 12 trees** (Slo-sister, EAG-prokaryotic, EAG-sister, KCNQ-sister); in the other 8 it is split by the ingroup. EAG is the only family with one clade under both outgroups (root UFBoot 100 in each), and **its root moves with the outgroup** (basal-split Jaccard 0.0); KCNQ-sister's root edge is UFBoot 68. Monophyly re-checked independently with Bio.Phylo: same verdict on all 12. **Where a root exists it isolates 1–4 sequences, never two substantial clades**: Slo-sister an *Acanthamoeba* basal pick (1 vs 144), EAG-prokaryotic a *Tetrahymena* EAG (1 vs 167), EAG-sister *Guillardia* + one other protist (2 vs 166, UFBoot 53/44), KCNQ-sister 4 vs 95 (UFBoot 57/29) — the pattern long-branch attraction produces (a reading, not a test). **Near-miss measure** (new, descriptive only: fewest ingroup sequences sharing a side with the whole outgroup, `s7_newick.outgroup_intruders`): K2P-sister **1** (*Hydra*), Shaker-prokaryotic **4** (3 *Hydra* + *Nanorana*); the rest 24–157, i.e. outgroup sequences land on opposite sides of deep splits. The basal picks sit disproportionately on the outgroup side of those splits (CNG 18 of 23 picks among 157 of 336 sequences; Slo 27 of 39 among 80 of 145). Ingroup support unchanged in kind (median 36–62 % of edges at UFBoot ≥ 95). The rooting of these six families passes to S11 (non-reversible-model root, reconciliation), which must state it. Figure: panel C of `results/phylogeny/figures/tier1_trees.png` (three bars per re-rooted family). → `results/phylogeny/tier1_report.md` § 4, `tier1_reroot.tsv`, `tier1_reroot_trees.tsv`, treefiles `tier1/reroot/` |
| S8a | **Tier-2 design + inputs + detached run; tier-3 fold network measured** (S8 split: the tier-2 trees need hours of IQ-TREE) | S6, S7b, S7d | completed 2026-10-02 | **Rules fixed before any tree or structure comparison (D48).** **Units**: 7 alignable multi-family superfamilies get a tier-2 tree — P-loop, iGluR, Ca²⁺-release, innexin clan (pore modules, one tip per module/repeat) and Cys-loop, P2X, DEG/ENaC (full length); 21 single-family superfamilies are their tier-1 tree (8 have none); **4 refused by D27** (CFTR, TMEM16/OSCA/TMC, MscL/MscS, viroporin) → fold network. **Representatives (D8)**: per family × panel group × module, central-first greedy clustering at 0.5 within-family identity — except where that leaves > 4 tips per informative site: the P-loop pore module (~105 aa) trims to 85 columns, so 840 tips → **317 at 0.3**, fixed on alignment properties alone. Tips: P-loop 319, Cys-loop 282, DEG/ENaC 101, innexin 59, P2X 33, iGluR 30, Ca²⁺-release 30; informative sites 85–451. **Root**: the superfamily's `root_with` family, its tips restricted to the kingdom of its catalogue exemplars (+ exemplars added): the 4 animal `plgic_prok` members are Cys-loop *ingroup* tips, their placement the answer to the S7a row. **ITPR module span checked against structure** (S6 row closed): on all 3 human ITPRs it starts at UniProt TM5, ends 2–3 past TM6 and holds the GVGD filter; on cryo-EM 6DQJ (ITPR3) it spans TM5 helix → pore helix → filter → TM6. **Fold network**: 74 census families with an AFDB model (RyR none — too long) + a Kv VSD node, cut to comparison units, 2,775 pairs (TM-align + Foldseek). **Literature edges (D48 rule: median TM ≥ 0.5 and mutual best other superfamily): supported 3/5** — iGluR–P-loop pore 0.64 (138 pairs), innexin clan–connexin 0.57, Hv1–Kv VSD 0.59; **not distinguished 2** — TMEM16/OSCA/TMC 0.48 (just under the bar, though twice its best outside partner 0.24; Foldseek E 1e-2–1e-63), ITPR–P-loop 0.56 (ties ITPR–iGluR 0.56: a P-loop-like pore, not placed nearer P-loop). **Unasserted**: PACC1–DEG/ENaC 0.51 (Foldseek E to 1e-7) — emergent. Tier-2 IQ-TREE run detached (`nohup caffeinate -i python3 -u scripts/s8_tier2.py run --jobs 4 --threads 2`, log `<data root>/trees/s8/run.log`). → `results/phylogeny/tier2_report.md`, `results/phylogeny/fold_network/`, figure `results/phylogeny/figures/fold_network.png` |
| S8b | **Tier-2 trees parsed**: finish the detached run, parse (outgroup monophyly + root UFBoot, per-family monophyly, support), the `plgic_prok` placement, report § 6, tree panels in the S8 figure | S8a | completed 2026-10-02 | **7 / 7 tier-2 trees, 0 failures** (9.1 run-hours summed; Cys-loop 5.5 h). Read under D48 unchanged (`s8_tier2.py parse` → `s8_parse.py`). **Roots: 3 of 7 rooted on the declared outgroup** — P2X (non-metazoan P2X one clade, UFBoot 100), Ca²⁺-release (ITPR modules one clade, 100), iGluR (single GluR0 tip, no support value); **3 outgroups split by the ingroup** — P-loop (KcsA/MthK/NaK/*Aliivibrio*: 32 tips of Slo/Shaker/KCNQ/SK/K2P in the way), DEG/ENaC (ENaC nests inside the invertebrate degenerins), Cys-loop (see below) — read unrooted, never re-rooted; innexin clan none declared. Where rooted, P2X and Ca²⁺-release isolate a single tip at the base (the S7d pattern); iGluR splits 7 of 8 non-vertebrate iGluR modules from the rest (UFBoot 96/82). **Support**: share of ingroup edges at UFBoot ≥ 95 — DEG/ENaC 0.51, Cys-loop 0.37, innexin 0.38, iGluR 0.19, P-loop 0.15 (319 tips on 90 columns: its deep order is unresolved and is not read), P2X 0.13, Ca²⁺-release 0/3. **Families as clades** (per module for multi-module chains): **34 of 54 groups one clade** — every Nav and NALCN repeat one clade at UFBoot ≥ 95, TPC repeat I one clade, Cav repeats scattered (a 90-column resolution limit; S11 owns the repeat order); near misses are nested relatives (Kv modifier in Shaker, ENaC in degenerins, 5-HT3 + ZAC in nAChR, HCN beside one CNG tip). **The four animal `plgic_prok` proteins** (*Branchiostoma*, *Aplysia*, *Lottia* ×2) **are nested in no Cys-loop family**: they form one clade (UFBoot 100) that joins ELIC (69) then GLIC (87) — the six-tip set is one side of an edge, separated from all 276 eukaryotic Cys-loop tips, and it is because they sit between GLIC and ELIC that the declared outgroup is split. Horizontal transfer, an ancient retained lineage, or long-branch attraction — passes to S21. Self-test +5 (`selftest_s8.py`). Figure `results/phylogeny/figures/tier2_trees.png`. → `results/phylogeny/tier2_report.md` § 6, `tier2_trees.tsv`, `tier2_families.tsv`, `tier2_placement.tsv`, treefiles `tier2/` |
| S9 | Selectivity-filter atlas: the filter locus for every P-loop family, and whether it is congruent with the tier-2 tree (Q5) | S7b, S8b | pending | |
| S10 | Repertoire evolution: presence/absence of every family across the lineage panel; ancestral reconstruction; the MscS animal-absence question (Q4) | S5b, S8b | pending | |
| S11 | Duplication history: which families expanded in which lineage, 2R/3R ohnologue status, and the Kv/Nav/Cav repeat-duplication order | S7b, S10 | pending | |
| S12 | Structures: AFDB coverage per family, TM-align against the experimental reference set, Foldseek all-vs-all for the fold network | S2, S6 | pending | |
| S13 | ML selection tests on the families with clinical variant sets (CFTR, SCN1A, KCNQ1, RYR1) | S6, S7b | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3a, S3b, S5b, S7b, S8, S9, S10, S12, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Method contribution** — what would domain search alone have missed, per superfamily (Q3); the recall curve as methods are added | S2, S3a, S3b, S5b | high | completed 2026-09-29 | **Q3 answered: domain search finds the channels but cannot name most of them.** Final census frame (census v4, high-confidence census-family profile calls, intact genome loci): **7,196 members — 94.2 % carry an enumerated pore signature, but S2's domain rules call only 44.1 % to the right family**; profiles add 96 proteome members never enumerated (1.3 %), genomes 323 loci (4.5 %). **Independent human frame (320 curated genes): enumerated 319, right family by domain rules 169, by profiles 319.** Per superfamily (thresholds fixed first, ≥ 95 % / ≥ 50 %): enumeration reaches ≥ 95 % in 22 of 25 superfamilies — the exceptions **CLIC 77 %, Hv 35 %, viroporin 0 %** (search-space and pore-model gaps, S3b/S4 rows); **the domain call is 0 % in Cys-loop, DEG/ENaC, P2X, CLC and viroporin and 0.3 % in iGluR** — identical architectures inside each superfamily (D25), so only profiles separate the families. By lineage the enumeration misses concentrate in CLIC (plants, ciliates 0 %), Hv (invertebrates 0 %), basal-metazoan P-loop (89 %). jackhmmer completeness per superfamily reported beside the curve (clean runs). The two new-instrument designs (superfamily-seeded jackhmmer, six-frame profile scan) stay open. → `results/method_contribution/report.md` |
| S16 | **Annotation-quality audit** — how often a real channel locus is missing, fragmentary, split, unnamed or filed under the wrong family across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5b, S15 | high | pending | |
| S17 | **Mechanism vs clade** — do the CLC channel/transporter and anoctamin channel/scramblase splits survive the tier-1 trees (Q7) | S7b | high | pending | |
| S18 | **Convergence audit** — connexin vs pannexin, TMEM175 vs the GYG channels, the mechanosensitive families: is "unrelated" a fact or a detection limit (Q6) | S8a, S12 | high | pending | |
| S19 | **Channelopathy map** — clinical variants of every family mapped onto the pore module and the constrained core | S6, S12 | medium | pending | |
| S20 | **Auxiliary-subunit census** — the excluded 15 %: how many there are, and how often they are counted as channels in published totals | S2 | completed 2026-09-29 | **The 240–400 spread of published human channelome counts is scope, measured.** Three curated database lists, archived with releases: **GtoPdb 2026.3 285, HGNC group 177 331, UniProt KW-0407 (2026_03) 338 — 400 in their union**; the **238 genes on all three are all pore-forming census genes**. Of the union, 314 are census pore genes (6 on no list: SCN7A, TMC3/5–8) and 86 are not: **44 auxiliary subunits** (39 catalogued + 5 uncatalogued: KChIP1–4, TMEM37), 14 aquaporins, 7 CLC/SLC26 transporters, 8 enzymes/transporters with the 'Ion channel' keyword, 8 proposed pores the catalogue lacks (PACC1, TMCO1, TMEM87A, TMEM109, CLCC1, CCDC51, GPHRA/B), 2 claudins, 3 pseudogenes (TRPC2 among them) — all 27 uncatalogued genes classified by hand (CURATED). Lists miss census pore genes by scope: GtoPdb 57 (no scramblase, TMC, CLIC, tweety, OSCA, otopetrin, CALHM, LRRC8, bestrophin), HGNC 30, UniProt 39 — **UniProt KW-0407 is on none of the 21 connexins**. **Auxiliaries: 76 human genes would inflate the 320 by 24 %** (the scope doc's 'roughly 15 %' corrected); the lists count 1.4 % (GtoPdb), 6.0 % (HGNC), 10.7 % (UniProt) of their totals as auxiliaries, each a different set. **6 of 11 auxiliary families pool unrelated proteins** (all-vs-all phmmer, E ≤ 1e-3): 22 homology groups. **Panel**: of 3,053 S3b profile calls to auxiliary families, **1,103 pass a reciprocal-best-human-hit test**; 1,949 fail — LRR (SLIT, LRR proteins) and Ig (hemicentin, titin) proteins called by the LRRC-γ and Navβ parts of pooled profiles; 4 groups (CATSPERZ, BSND, CNIH, EMRE) are unreachable by their pooled profile. 13 groups found only in vertebrates (not an absence claim). → `results/auxiliary/report.md` |
| S21 | **Prokaryotic and viral channels** — the outgroup census, and whether the eukaryotic families have prokaryotic sisters | S3b, S5b | medium | pending | |
| S22 | **Pharmacology overlay** — which families carry approved-drug targets, against family size and annotation quality | S2, S19 | low | pending | |
| S14c | **Manuscript rewrite pass** — one full pass once every analysis has landed | S14a, S24 | medium | pending | |

---

## Emergent tasks & new aims

Anything discovered mid-session that deserves its own work goes here rather
than expanding the task in progress.

| Added | From | Task | Status |
|-------|------|------|--------|
| 2026-08-19 | setup | Confirm the external drive is attached and `data_root.txt` points at it before S4/S5 | closed 2026-09-28 (S4 ran on it) |
| 2026-08-19 | setup | **The TRP pore-model finding (H7) may generalise.** Pfam's coverage of pore modules was measured only on the exemplars. Enumerate which *census families* have any pore model at all, and how many members of each carry it — the answer sets S2's real recall ceiling | closed 2026-09-30 (S15 measured it per family and superfamily: `results/method_contribution/curve_by_family.tsv` — domain enumeration reaches ≥ 95 % in all but CLIC, Hv1, viroporin, MITOK) |
| 2026-08-19 | setup | **`ZACN` breaks the H2 rule.** A human Cys-loop channel with no annotated TM domain. Check whether this is a UniProt annotation gap or a real truncation, and whether other single-exemplar families have the same problem | closed 2026-09-30 — a **Pfam coverage gap**, not truncation or a UniProt gap: UniProt annotates 4 TM helices (234–389, protein-level evidence) but Pfam gives ZACN only the LBD (`PF02931`), not `PF02932`; S2 now stops at the superfamily (positive rule, S2b) and the profiles call it `zac` |
| 2026-08-19 | setup | **69 exemplars carry no UniProt accession** in the catalogue and are resolved live each run. Fold the S0-resolved accessions back into the catalogue files so the reference panel is reproducible offline | closed 2026-09-30 — all 50 still missing folded in from S0's resolution (guard: status ok and resolved *primary* gene = declared gene; 50/50 passed); 183/183 exemplars now carry an accession |
| 2026-08-19 | setup | **The four-repeat anchor works; the equivalent for other superfamilies does not exist yet.** A Cys-loop charge-selectivity anchor and a CLC gating-glutamate anchor would close hazards H14 and H5 the same way | open (S9) |
| 2026-08-19 | setup | Biopython is present in the `piezo1` env but absent from the base interpreter, so `src/analysis` and the GUI only run there. Either pin the env in a wrapper script or drop the Biopython dependency from `alignment.py` in favour of `src/utils/mafft.py` | closed 2026-09-30 — `bin/envpy <script>` runs the piezo1 interpreter (override `ION_CHANNEL_PY`); CLAUDE.md updated |
| 2026-08-19 | S0 | **Three of the four superfamily outgroups are unreachable by their own superfamily's rules.** GLIC and ELIC carry the Cys-loop LBD and no TM model, so the H2 rule rejects them; GluR0's architecture is a *potassium-channel* pore plus a bacterial binding domain, with no iGluR model at all; prokaryotic NavAb carries neither Nav-specific domain. Rooting therefore depends on S3's profile methods, not on domain search — and that is a result about annotation coverage, not a bug | open (S3 → S7/S8) |
| 2026-08-19 | S0 | **`PF16799`, the positive Hv1 test, is on human HVCN1 and not on *Ciona* Hv1.** The H9 discriminator is currently a mammal-only instrument; check how many other family-level tests are human-only by re-running the S1 panel on non-human orthologues | open (S1) |
| 2026-08-19 | S1 | **The benchmark re-fetches its whole panel on every run.** 97 proteins × 3 API calls, several minutes, repeated for each of the four runs the setup session needed. The parent project cached its control panels to `panel_positives.json` / `panel_decoys.json` so the benchmark reruns offline; do the same here — the evidence (sequence, Pfam counts, TM count) is what the classifier consumes and it does not change between runs | closed 2026-09-30 — `s1_benchmark.py` caches the evidence in `results/benchmark_controls/panel_evidence.json` (`--refresh` refetches) |
| 2026-08-19 | S1 | **The reference tier does not scale to a census.** Measured on the S1 run: ~2.5 s per MAFFT pairwise alignment × up to 12 prefiltered candidates ≈ 30 s per protein, so the 97-protein panel takes ~45 min and a 5,000-protein census would take ~40 hours. S2 needs either a calibrated fast identity estimate (DIAMOND or BLAST, checked against MAFFT on the S1 panel so the D7 margin keeps its meaning) or a profile-based assignment that skips pairwise scoring entirely. **This is a scaling result, not a bug** — the tier is correct and slow, and the fix must not change what the margin means | closed 2026-09-30 — superseded: the census never runs it (D31); S3a's profile assignment (D32) is the scalable sequence tier |
| 2026-08-19 | setup | **The phylogeny pipeline was smoke-tested on the Cys-loop superfamily** (8 exemplars, MAFFT → trimAl → IQ-TREE 2 with 1000 UFBoot, rooted on GLIC/ELIC, 52 s) and returned the expected topology: anion- and cation-selective receptors separate at 100 %, AChBP sisters the cationic clade. `results/phylogeny/tier2_cysloop/`. **This is not S8's tree** — eight sequences, exemplars only — and must not be cited as one | closed 2026-08-19 |
| 2026-08-19 | S0 | **`PF00005` (ABC_tran) has 1.66 million UniProt proteins and `PF00520` has 206,115.** The census search space is dominated by two accessions that are mostly not channels. S2 needs a per-signature triage rule before it enumerates, or it will fetch a million transporters to find one CFTR | closed 2026-09-30 — superseded: ABC domains are SHARED_WITH_DECOY and not enumerated; census v2 is 1.25 M (r3) / 1.27 M (r4) records |
| 2026-09-28 | S2 | **Three pore-union signatures are co-domains that mostly sit on non-channels.** Measured as the only pore signature on unassigned records: `PF00027` cNMP_binding **148,190**, `PF00497` SBP_bac_3 **115,927**, `PF13426` PAS_9 **111,788** — CRP-type regulators, periplasmic binding proteins and PAS sensors. They are declared `FAMILY`-level for kv_eag/HCN/CNG and iglur_prok, so `pore_signatures()` enumerates from them. Re-level to `ACCESSORY` (or add a pore-module requirement to those families' rules), re-run S1, and the ~376,000 records these three bring in alone leave the census space: 869,000 records remain, of which ~124,000 (~14 %) would still be unassigned. Check against `human_recall.tsv` that no channel is lost before adopting it. Supersedes the PF00005 triage row for S2's purposes | open |
| 2026-09-28 | S2 | **Derived family rules require a family's *complete* declared architecture, and many real members carry part of it.** Unassigned records carrying part of exactly one family's architecture: MscS **65,506** (mostly `PF00924`±`PF21082` without `PF05552`), OSCA 4,505, TRPM 3,960, Piezo 3,939, Slo 2,243, RyR 1,503, Kir 1,081. A "core pore signature suffices" rule per family would call them; it must be benchmarked on S1 (and on a non-human panel) before adoption, not tuned on the census. `results/census_v2/partial_architectures.tsv` | closed 2026-09-30 — the profile tier recovers them: of the records S2 left unassigned, census v3a calls ~61,900 MscS (51,184 high), 2,831 TRPM, 2,817 OSCA; a rule rewrite would duplicate what S3a does |
| 2026-09-28 | S2 | **Four-repeat filter strings with no call.** Projected and uncalled: NEEE 1,267, DEEA 435 (the invertebrate Nav2/BSC1 calcium-selective filter, if the literature confirms), QEEE 221, DENA 164, DDDD 129, plus gapped reads. Each needs a literature-backed `FILTER_CALLS` entry or a recorded "no call" — never an inferred one. S9's atlas | open (S9) |
| 2026-09-28 | S2 | **ZACN is still called AChBP**, now in the census (the only wrong-family call among 320 human genes). The H2 rule's premise — an LBD without an annotated TM domain is AChBP — is false for ZACN. Needs a positive AChBP test (secreted, no TM features, `tm_count` = 0) rather than absence of `PF02932` | closed 2026-09-28 (S2b; ZACN is `zac` by profile) |
| 2026-09-28 | S2 | **`results/session_live.json` is untracked scratch** written by every driver. Gitignored this session | closed 2026-09-28 |
| 2026-09-28 | S3a | **Hazard rule `H4-itpr` is an absence test and calls RyR fragments ITPR.** It fires on `PF08709` *without* `PF02026`/`PF06459`; `PF08709` is shared by ITPR and RyR. The external check found **502** records the IP3R project calls RYR that census v3a calls ITPR — all this rule, profile abstaining, median 76 aa — and **1,431** v3a ITPR calls rest on it alone. Rewrite as a positive test (an ITPR-only domain, or the S3a profile margin), re-run S2 classification on the affected records, re-merge | closed 2026-09-28 (S2b: ITPR has no positive architectural test; 502 swaps → 0) |
| 2026-09-28 | S3a | **Hazard rule H2 (AChBP = Cys-loop LBD without annotated TM) is an absence test and fails at scale** — ZACN was the first case. 6,633 S2 AChBP calls are profile-called receptors (nAChR 3,726, GABA-A 1,566, …), median length 383 aa against AChBP's 210–240. Needs a positive AChBP test (length band + no TM features + profile) | closed 2026-09-28 (S2b, D33) |
| 2026-09-28 | S3a | **H13 (TRPP vs polycystin-1) miscalls polycystin-1-like proteins with few PKD repeats**: 2,989 S2 TRPP calls are profile-called `assoc_polycystin1`, median 2,263 aa (TRPP is 600–900). The `PF00801` copy-number discriminator needs a length/profile positive test | closed 2026-09-28 (S2b: PLAT/REJ/GPS or PF00801×≥5 ⇒ polycystin-1; PF18109 ⇒ TRPP) |
| 2026-09-28 | S3a | **`PF00520` assigns the P-loop superfamily in S2, but Hv1, VSP, ITPR and RyR (other superfamilies) carry it too.** 2,512 v3a conflicts are S2-`ploop` records the profile calls hv1 (1,973, median 252 aa — Hv1 length, probably the non-mammalian Hv1s the mammal-only H9 test misses), itpr 375, ryr 83, vsp 81. Either the signature is declared shared across those superfamilies or S2's superfamily call from it is dropped | open |
| 2026-09-28 | S3a | **`s0_lib.resolve_gene` matches synonyms**: `gene_exact` hits synonyms and the longest entry wins, so TRPC7 → TRPM2, GRIK2 → GRIK5, KCNG3 → KCNG4, GJC3 → GJE1, AQP7 → AQP9, CACNG6 → CACNG8, KCNMB2 → KCNMB3. S3a uses `resolve_primary()`; S0's exemplars and S2's human recall were checked and are unaffected. Fix the shared helper before anything else calls it | closed 2026-09-30 — an entry whose *primary* gene name matches now wins; a synonym-only match is the flagged fallback (`synonym_match`); TRPC7, GRIK2, KCNG3, GJC3, AQP7 now resolve to themselves |
| 2026-09-28 | S3a | **Four auxiliary-subunit families pool unrelated proteins** (S20 measured six: see its row) (`assoc_k_beta`: Kvβ aldo-keto reductases + BK β + LRRC γ + DPP6/10; also `assoc_mcu_reg`, `assoc_catsper_aux`, `assoc_clc_aux`). One profile cannot represent them; LOO decoys from them get no hit. Split into homologous families (S20) — no census effect, none carries a pore signature | open (S20) |
| 2026-09-28 | S3a | **The profile-only superfamily splits (Cys-loop, iGluR, DEG/ENaC, P2X, CLC, TMEM16) have no non-vertebrate test.** S2 never called these families, so calibration cannot reach them, and the orthologue benchmark is 576/583 chordate. Build a non-vertebrate labelled panel (e.g. *Drosophila* / *C. elegans* receptors with literature subunit identities) before S7 uses these calls | open |
| 2026-09-28 | S3a | **2,258 PIEZO-census (v5) records are unassigned in v3a.** The parent list includes fragments and short-motif hits; check whether these fall below S3a's coverage gate or are PIEZOs missed | closed 2026-09-30 — `scripts/s3_piezo_unassigned.py` → `results/census_v3/piezo_unassigned.tsv`: 2,196 (97 %) fall below D32's 30 % coverage gate (partial Piezo matches: 1,409 fragment-flagged or < 500 aa, 543 of 500–2,000 aa), 62 no hit / < 30 bits. **244 long (≥ 2,000 aa) records below the coverage gate** are the residue worth a look (split or chimeric models, or Piezo-module proteins) |
| 2026-09-28 | S2b | **H11-abcc and H12-kctd are still absence rules** (ABC architecture *without* the R domain ⇒ SUR; T1 *without* a pore module ⇒ KCTD). H12 accounts for 358 of the 3,288 v3a conflicts (profile calls Kv). Rewrite as positive tests the way S2b did H2/H4/H13; the S2b self-test invariant already names them as exempt | closed 2026-09-28 (S2c) |
| 2026-09-28 | S2b | **Census v2's `fragment` column stores UniProt's raw `flag`, which includes "Precursor"** | closed 2026-09-30 — S3's R3 seed rule and set-B filter now test for "Fragment" (D44); the column itself keeps the raw flag, read correctly everywhere it is used |
| 2026-09-28 | S2b | **S3a's R3 seeds were drawn from S2 r1 calls**, six ITPR and five TRPP R3 seeds were called by the now-removed H4 / H13 absence rules (all reviewed, non-fragment, inside the family length band ±25 %, so not the fragments and polycystin-1-like proteins those rules miscalled). Profiles were kept, not rebuilt, so the S2b before/after is a comparison on one instrument. A rebuild at S6 should re-draw R3 from the current census | open (S6) |
| 2026-09-28 | S2c | **The S3a profile library has no decoy for generic ABC transporters, so 525 records (508 bacterial) are called SUR by the profile alone.** They are ~1,000-aa ABC transporters with cNMP + C39 peptidase domains; the SUR profile (2 human seeds) is the only ABC-transporter profile they can score against. Profile agreement is evidence only when every family that could score has a profile. Fix: a catalogued non-channel ABC-transporter decoy family (human ABCC1–6/10–12 + a bacterial peptidase-containing ABC transporter as exemplars), its profile swept (~one profile), re-merge. Adds a 91st family — the user's call | closed 2026-09-28 (S3a2) |
| 2026-09-28 | S2c | **S2b silently shrank the census search space**: re-levelling PF20519 to SHARED_WITH_DECOY dropped it from `pore_signatures()` (67 → 66), because enumeration was derived from evidence level. Fixed by `Signature.enumerate` (D34) and a self-test pinning the union to the 67 signatures census v2 enumerated. No census data were affected (nothing was re-enumerated) | closed 2026-09-28 |
| 2026-09-28 | S4 | **The viroporin family is outside the census search space.** Its four signatures (`PF00599` Flu_M2, `PF00558` Vpu, `PF02723` CoV_E, `PF11289`) are declared `SUBFAMILY`, and `Signature.enumerate` defaults to SUPERFAMILY/FAMILY (D34), so census v2 never enumerated them: 177 viral records in the census, **0** in the three panel viral proteomes, and M2 (P06821) absent. It is the only census family with no enumerated signature (checked over all 68). Fix: `enumerate=True` on those four — a change to the search space and so a new census revision (D34), not an edit. S3b's viroporin profile covers the panel meanwhile | closed 2026-09-30 — census v2 **r5** (D43 delta): the four signatures set `enumerate=True` at SUBFAMILY level; 4,007 records (3,994 viral), r4's 1,271,983 + delta = UniProt's 1,275,990, all 36 shards and 79 signatures exact; no earlier record carries one (0 re-classified, 0 v3a calls changed); **3,560 called viroporin by profile** (3,554 high), 447 unassigned |
| 2026-09-28 | S4 | **The panel is 52 named species, a thin denominator for repertoire evolution.** It suffices for S3b/S5's completeness argument on the declared panel, but S10's gain/loss reconstruction and any "absent from lineage X" claim beyond these species needs density — e.g. one reference proteome per order by the same selection rule (the parent project swept 6,928 eukaryotic proteomes). A scope decision for the user before S10, not S4's | closed 2026-09-30 — user chose one proteome per order (S4b, D46) |
| 2026-09-28 | S4 | `INTERFACE.md` said the species table had 51 species; it has 52. Corrected | closed 2026-09-28 |
| 2026-09-29 | S3b | **D32's coverage gate does not exclude a shared module when the module dominates the profile.** TRPN's profile is mostly ~29 ankyrin repeats; 37 of its 38 panel calls outside census v2 are ankyrin-repeat proteins (ANKRD52, ankyrin-2, IκB, titin, synphilin) at medium confidence, coverage 0.30–0.46; LRRC8 (LRR half) 20 medium, TRPA/TRPC/TRPV similar. Fix candidates: require coverage of the *pore* segment of the profile (profile coordinates of the TM/pore region from the seed alignment), or re-level those calls to `module` when the covered span lies inside the repeat array. Must be benchmarked on S1/S3a before adoption; v3a is affected wherever S2 abstained | open |
| 2026-09-29 | S3b | **K1 is blind to drift into proteins no profile calls** (ankyrin/LRR/GST-fold repeats): TRPA/TRPN pre-kill rounds hold ~6,400 `module` targets, LRRC8 and CLIC K3 runs balloon to 9,025 and 1,565 accepted with 83 % / 92 % uncalled. The parent project's S19 measured the same blind spot and proposed a drift rule on the uncalled share of the finished model. Evaluate that rule on `jackhmmer_rounds.tsv` (`uncalled_frac` is recorded per round) as a classifier before adopting it — never retune D10 on the run it will judge | open (S19-like methods task) |
| 2026-09-29 | S3b | **In multi-family superfamilies a family-seeded jackhmmer is a superfamily search.** K1 fired at round 2 for 15 runs, 12 of them P-loop or iGluR; the Cys-loop, DEG/ENaC and P2X runs all converge on one set per superfamily. A superfamily-seeded design (one run per alignable superfamily, K1 on the *other-superfamily* share, recorded per round as `other_sf_frac`) would give a cleaner completeness argument for S15. Design note, not adopted | open (S15) |
| 2026-09-29 | S3b | **Hv1 is the family domain search misses most among real channels**: 25 of 33 panel Hv1 (molluscs, moss, fish, cnidarian, placozoan) carry none of the 67 enumerated signatures, 15 at high confidence — consistent with the S0 finding that the Hv1 domain `PF16799` is mammal-only. CLIC likewise (27 high, plants/ciliates/invertebrates). Both argue for profile-based enumeration in the next census revision | closed 2026-09-30 — measured in S15 (Hv 35 % enumerated) and S3b; profile-based enumeration for the next revision stays a census-revision decision |
| 2026-09-29 | S5a | ~~**`-G` has no margin at the top.**~~ Closed by D38: a 1.5 Mb margin measured worse on mouse; `-G` stays 1 Mb. Mouse Asic2's widest annotated intron is 996,015 bp against the 1 Mb `-G` for ≥ 1 Gbp genomes; 13 / 220 confirmed mouse loci have an intron > 200 kb. Before S5b: set `-G` from the widest measured intron with a stated margin (e.g. 1.5×) and measure what a larger `-G` does to chaining (the ZAC and Ugt1a10/TRPM chained models appeared at 1 Mb) | closed 2026-09-29 (S5b, D38) |
| 2026-09-29 | S5a | **D4 bar missing for families the pilot never confirmed on annotation** (3 informative cells: `plgic_prok`, `iglur_prok`, `deg_invertebrate`); prokaryotic and viral genes are intronless, so their bar could be the family's length band × 3 bp. S5b measures more; state the fallback rule before reading any cell | closed 2026-09-29 (S5b, D38) |
| 2026-09-29 | S5a | **Genome copy number is not locus count.** Mouse: 12 high-confidence VDAC loci for 3 genes, most with frameshifts or in-frame stops (retrocopies). Any genome-derived copy number (S11) must count intact, exon-confirmed loci | open (S11) |
| 2026-09-29 | S5a | **Partial loci at non-informative cells are shared modules, not misses** (e.g. yeast: 13 families partial — Kv, SK, CatSper, TPC, HCN, TRPs — from VSD/cNMP/CaM-binding modules). S16 must not read `partial` as an annotation failure | open (S16) |
| 2026-09-29 | S5a | **Unmatched cells cannot be judged by this instrument at all**: a family present in one panel group only can never be declared absent in another (no in-group bait). A profile-HMM scan of six-frame-translated genomes would reach further than pairwise baits; a method design for S15, and another argument for the panel-density row (S4) | open (S15) |
| 2026-09-29 | S5b | **The *Takifugu* reference proteome carries no ryanodine receptor**, while its genome (GCA_901000725.3) has 6 intact, full-coverage RyR loci at 61–80 % identity. The cleanest proteome miss in the panel — an annotation gap for S16 (check the proteome's source assembly, FUGU5 vs fTakRub1, and whether UniProt holds RyR fragments outside the canonical set) | open (S16) |
| 2026-09-29 | S5b | ***Nematostella* MscS `genome_found` sits on 6.7 kb and 3.1 kb contigs** (bait *E. coli* MscS, 41–43 % identity). Tiny contigs carrying a bacterial-like gene are the signature of contamination, and MscS in animals is Q4 itself. Before S10 reads it: contig composition/coverage, neighbouring genes, and whether any other cnidarian carries it | open (S10) |
| 2026-09-29 | S5b | ***Daphnia pulex* carries 5 of the 72 absences** (Piezo, TRPP, DEG, CALHM, pLGIC_prok) on the 2011 assembly (N50 0.64 Mb, scaffold-level). Piezo absent from an arthropod is surprising. Check against a current *Daphnia* assembly before S10 treats any as a loss | open (S10) |
| 2026-09-29 | S5b | **ZAC `absent` in 9 vertebrates** (chicken, *Xenopus*, platypus, coelacanth, both chondrichthyans, zebrafish, *Takifugu*, rat) with a `pooled` D4 bar. Only rat is a literature expectation; the rest are a result to check against ZACN's published distribution before it is stated | open (S10) |
| 2026-09-29 | S5b | **Five annotated genes have an intron beyond the `-G` they were searched with** (4 in genomes < 1 Gbp at 200 kb: *Aplysia* TRPC, amphioxus and *Callorhinchus* ASIC, *Callorhinchus* TRPM3). All were called. A per-genome `-G` from measured introns would remove the edge; low priority while every such gene is still found | open |
| 2026-09-29 | S5b | **Viroporin and the four single-species groups cannot be judged by a bait instrument at all** (viruses 0/1 matched detection each; algae, amoebozoa, apicomplexa, excavate have no in-group bait). Presence/absence there rests on the proteome alone — state it in S10/S21 rather than fill it | open (S10, S21) |
| 2026-09-29 | S6 | **trimAl `-automated1` keeps 3–4 % of columns on the largest, most divergent families** (K2P 145 / 5,149, nAChR 208 / 8,159, CNG 365 / 8,643; median over families 22 %). A 145-column alignment is thin for a 380-sequence tree. Before S7: compare automated1 with a fixed gap threshold (e.g. `-gt 0.5`) or ClipKIT on a stated criterion (retained informative sites, UFBoot on a small family), decided before the census trees are read | closed 2026-09-29 (S7a, D41: `-gt 0.5`) |
| 2026-09-29 | S6 | **ITPR's pore module rests on a RyR-seeded vote plus the helix snap**, the only module span not checkable against an annotated relative (no ITPR in UniProt carries an annotated pore loop, and RyR held out leaves no seed). A structure-derived span (rat ITPR1 cryo-EM, 3JAV/6MU2) would close it before S8's `ca_release` tree | closed 2026-10-02 (S8a: consistent — UniProt TM5 start, TM6 end +2–3, GVGD inside, on all 3 human ITPRs; on 6DQJ TM5 helix → pore helix → filter → TM6; `results/phylogeny/tier2_itpr_span.tsv`) |
| 2026-09-29 | S6 | **Five census families have < 4 sequences in their D39 set** (delta-GluR 3, GluR0 1, tweety 3, viroporin 3, ZAC 2) and get no tier-1 alignment. Tweety and viroporin are whole superfamilies; their tier-1/2 trees need either medium calls (flagged) or the wider census v3a | open (S7) |
| 2026-09-29 | user | **Headline figures for S1–S5b.** The protocol now requires one figure per completed task in the README and dashboard (end-of-session step 3b); only S0 has one. Candidates, all from committed tables: S1 per-family recall + per-tier attribution; S2 call status by superfamily; S3a profile vs S2 calibration; S4 panel species × BUSCO/N50; S3b domain-search misses by family × group; S5b the family × species presence matrix with genome verdicts. A `scripts/s<n>_figures.py` each, through `figstyle.py` | closed 2026-09-30 — `s1_figures.py`, `s2_figures.py`, `s3_figures.py`, `s4_figures.py`, `s3b_figures.py`, `s5_figures.py`; in README and dashboard |
| 2026-09-29 | S7a | **The P-loop families are rooted on two-helix prokaryotic pores that fill 5–16 % of a four-repeat or TRP alignment's columns.** The root then rests on the pore columns alone. If S7b finds the outgroup non-monophyletic or the root edge weakly supported, the fix is a per-family declared sister in the catalogue (e.g. Nav↔Cav, TRPC↔TRPV) — a catalogue edit with provenance, decided before rebuilding, never chosen by looking at which root looks best | → task **S7c** (user, 2026-10-01) |
| 2026-09-29 | S7a | **The `plgic_prok` D39 set holds four animal proteins** (*Branchiostoma*, *Aplysia*, *Lottia* ×2) called prokaryotic pLGIC at high confidence beside GLIC. Either real bacterial-type pLGICs in animals (horizontal transfer or an ancient lineage) or a profile built on two prokaryotic seeds calling divergent animal subunits. Check against the Cys-loop tree and the literature before the Cys-loop outgroup is read (S8, S21) | S8b read it: nested in no Cys-loop family — one clade (UFBoot 100) joining ELIC (69) then GLIC (87), apart from every eukaryotic Cys-loop tip; HGT vs ancient lineage vs long-branch attraction (and contamination: check the genomic context) → **S21** |
| 2026-09-29 | S7a | **Tree compute is the bottleneck for S8 too.** At `-mset` 4 matrices, 726-sequence nAChR was 17 models into ModelFinder after 2 h (2 threads). S8's per-superfamily trees must pick representatives by D8 before any run, not after | closed 2026-10-02 (S8a, D48) |
| 2026-09-29 | S20 | **Split the six pooled auxiliary families into their 22 homology groups** (`results/auxiliary/aux_groups.tsv`): each group its own catalogue family with its own profile. No census effect (none carries a pore signature), but the profile library changes, so it is a library revision with a benchmark, not an edit. Until then the panel auxiliary counts are per homology group by reciprocal best hit, and four groups are unreachable | open (library revision) |
| 2026-09-29 | S20 | **Five auxiliary genes are missing from the catalogue** — KCNIP1–4 (KChIPs, Kv4) and TMEM37 (Cav γ-like), all on UniProt KW-0407. Add to `channel_associated` families with provenance; no census effect | closed 2026-09-29 (user-directed: `assoc_kchip` added; TMEM37 into `assoc_cav_aux` as its own PF15108 group) |
| 2026-09-29 | S20 | **Eight proposed pore-forming channels are not in the catalogue** — PACC1 (proton-activated Cl⁻ channel, solved structures) above all; also TMCO1, TMEM87A, TMEM109, CLCC1, CCDC51 (MITOK), GPHRA/B. Each is a candidate census family, i.e. a change to the search space and a new census revision (D34) — the user's call, with literature checked first | closed 2026-09-29 (user-directed: added in `src/catalogue/proposed.py` — PACC1 `channel`, seven `channel_contested`; decoys EMC3, GOST, BRI3BP; H17–H19; declared not searched, D42 — see the census-revision row) |
| 2026-09-29 | S20 | **The review still quotes the human channelome as '240–400' from the literature** (`docs/review/01_introduction.md`, `14_open.md`). S20 measured it on three databases (285 / 331 / 338 / 400, 238 shared, all pore). Update the review sources and rebuild at the next review pass | open (S14c) |
| 2026-09-29 | user | **Census revision for the catalogue additions.** The seven new census families (pacc, tmco1, tmem87, tmem109, clcc1, mitok, gphr; 8 human genes) and four new controls (assoc_kchip, nonchannel_emc3, nonchannel_gost, nonchannel_bri3bp) have **no census members and no S3a profile**: their signatures are `enumerate=False` (D42). To bring them in: enumerate their FAMILY/SUPERFAMILY signatures (a census v2 revision), build their profiles alone as S3a2 did (the other 91 frozen), sweep them over census v2 and the S4 panel, benchmark (hazards H17–H19 need the decoy profiles to be testable), re-merge v3a/v3 and extend census v4. MITOK carries no Pfam domain — profile-only. Until then every census count stays on the 68-family / 320-gene catalogue | closed 2026-09-29 (S2d: census v2 r4 → v3a → v3 → v4, D43) |
| 2026-09-29 | S2d | **The genome sweep has not searched the families added in r4.** S5b's bait panel predates them, so census v4 carries their proteome rows only, and their presence matrix has no genome verdicts (no `absent` can be stated). Extend the bait panel with the 12 new families' baits (B1–B3) and run miniprot + the verdict chain for them alone; the existing cells must come out unchanged | closed 2026-09-30 (S5c: 0 of 3,536 old cells changed) |
| 2026-09-29 | S2d | **H18 after the reseed (2026-09-30): mammals split, fungi stop at the superfamily, non-mammalian TMEM87 is co-orthologous.** TMEM87B moved to its own decoy (`nonchannel_tmem87b`), GPR107/108 kept as GOST; TMEM87A/B now separate at margin ~0.60 (was 0.23), GPR107/108 at 0.96; yeast PTM1 / *S. pombe* GOST no longer called TMEM87 (superfamily only); TMEM87 calls 5,484 → 995. **But held-out *Xenopus* tmem87a is called TMEM87B** (the one set-B miss, 744/745), and the genome sweep shows TMEM87 'present' only in scattered species — invertebrate and fish TMEM87 are co-orthologues of A and B, and no profile can call them one or the other. **A scope decision for the user**: make the census family TMEM87 by descent (A and B, as the anoctamins are), or keep TMEM87A alone and state that non-mammalian TMEM87 stops at the superfamily | closed 2026-09-30 (user chose TMEM87 by descent — S2f, D45) |
| 2026-09-29 | S2d | **Seven r4 profiles rested on one seed** | closed 2026-09-30 — R2 exemplars added (reviewed orthologues, others held out) and R3's fragment test fixed: TMCO1 4 seeds, EMC3 5, CLCC1 6, GOST 4, TMEM87 2, TMEM87B 2, TMEM109 2; MITOK and BRI3BP stay single-seed (no reviewed orthologue to spare / BRI3BP's only held-out decoy kept). Set B 542 → 745 orthologues (744 correct), 59 decoys, 0 called a channel |
| 2026-09-29 | S2d | **Downstream tasks predate r4**: S6's D39 alignment sets and the running S7 trees hold no member of the new families; S20's auxiliary homology groups and panel test predate KChIP/TMEM37; S3b's jackhmmer completeness has no run for the new families. Each is noted in its report; extend when those tasks are next touched | open |
| 2026-09-30 | S5c | **Eleven controlled absences of the r4 families are measurements, not literature.** Two stand out: PACC1 absent from *Takifugu* (zebrafish and every other vertebrate carry it) and CLCC1 absent from all three ecdysozoans. Each is a candidate for S10's checks against published distributions before any is stated as a loss; the single-seed CLCC1/MITOK/TMEM109 profiles bound how far the instrument reaches | open (S10) |
| 2026-09-30 | S2f | **Are any non-animal GOST proteins TMEM87 orthologues?** None carries TMEM87's GOLD domain (PF21901, metazoan-only), yet plant, *S. pombe*, *Dictyostelium* and holozoan GOST proteins are called TMEM87 at high confidence against a GOST decoy seeded from animal GPR107/108, *Arabidopsis* CAND6 and yeast PTM1. A GOST-superfamily tree (TMEM87A/B, GPR107/108, fungal/plant GOST) answers it — and dates the TMEM87A/B duplication (the D45 premise). **Until then S10 restricts TMEM87 presence/absence to Metazoa**, and yeast TMEM87 'absent' is not a finding | closed 2026-09-30 — **yes: TMEM87 is a pan-eukaryotic lineage and its GOLD domain an animal addition.** GOST superfamily tree (`scripts/s7_gost_tree.py`, 115 panel proteins, Q.pfam+R7, 1000 UFBoot): the largest clade holding every animal TMEM87 and no animal GPR107/108 has **UFBoot 100** and contains 17 of 36 non-animal GOST proteins (plants, moss, both yeasts, *Dictyostelium*, *Plasmodium*, *Monosiga*, *Capsaspora*); the profile calls agree with the tree on 33/36 (the exceptions: the two yeast proteins, called GOST only because I had seeded the decoy with PTM1 — reverted — and *Plasmodium*, superfamily-only). A first reading took the smallest clade and concluded the opposite; corrected before any change was kept. The TMEM87A/B duplication is **not** dated by this tree (neither paralogue forms a clade) → `results/phylogeny/gost/report.md` |

| 2026-10-01 | S7c | **MloK1 fills only 13 % of CNG's trimmed columns** (SthK 53 %): the bacterial CNG outgroup is effectively one sequence plus a fragment of another. If CNG's root is unresolved in S7d, a third prokaryotic CNG (e.g. LliK, *Leptospira licerasiae*, if a UniProt entry resolves) is a catalogue addition for a later revision — not a substitution made after reading the tree | open — S7d: CNG unresolved, but the sister-family outgroup (HCN + EAG) fails as badly (152 ingroup sequences inside its clade) as MloK1 + SthK (157), so a third prokaryotic CNG alone is unlikely to rescue it; any addition is a new declaration with its criterion fixed first (S11) |
| 2026-10-02 | S7d | **No re-rooted P-loop family passes D47, and every defined root isolates 1–4 divergent sequences** (*Acanthamoeba*, *Tetrahymena*, *Guillardia* among them), several of them S7c basal picks. A design question for S11, not a re-run of S7c: (1) rooting without an outgroup — a non-reversible model (IQ-TREE `-m NONREV`/`UNREST` root, rootstrap) or gene-tree/species-tree reconciliation; (2) whether long-branch basal picks should be screened (e.g. a branch-length or TreeShrink rule) *before* any tree they would root. Either must be declared before it sees these trees' roots | open (S11) |
| 2026-10-02 | S7d | **Near misses on two trees**: K2P under Kir + Shaker has one *Hydra* sequence (A0ABM4D6C1) inside the outgroup clade; Shaker under KvAP/MVP/Kch has four (3 *Hydra*, 1 *Nanorana*). Worth checking whether those sequences are misassigned members (another family's sequence called K2P/Shaker at high confidence) — a census question, not a reason to drop them from a tree after reading it | open (S16/S11) |
| 2026-10-01 | S7c | **The K2P basal picks include plant and fungal proteins the K2P profile calls at high confidence** (S4b: 244 in Magnoliopsida alone; presumably the plant TPK and fungal TOK two-pore K⁺ channels — not checked). Whether these are K2P by descent or convergent two-pore architectures is a question the K2P tree can answer (S17-like); until then S10 should read plant/fungal K2P presence with that caveat | open (S10) |
| 2026-10-02 | S8a | **PACC1 resembles the DEG/ENaC fold**: median TM-score 0.51 over the three family pairs, Foldseek E 1e-5–1e-7 (ASIC, ENaC, degenerin vs PACC1), the strongest cross-superfamily pair no literature edge declares. Both are trimeric two-TM channels with a large extracellular domain; whether published PAC structures already note it, and whether profile-profile comparison (HHpred) finds sequence signal, is S18's question — no edge is added to `LITERATURE_EDGES` after seeing the matrix | open (S18) |
| 2026-10-02 | S8a | **TMEM16/OSCA/TMC misses the 0.5 bar by 0.016** while its within-superfamily median is twice its best outside partner (0.48 vs 0.24) and Foldseek finds every pair (E ≤ 1e-2). Whole-model units carry the cytosolic domains, which dilute an average-length TM-score; a TM-region unit would likely pass. A unit change is a new declaration for S12/S18, made before re-measuring — never applied to this verdict | open (S12/S18) |
| 2026-10-02 | S8a | **RyR has no AlphaFold DB model** (every exemplar exceeds the AFDB length limit), so the Ca²⁺-release node rests on ITPR alone, and ITPR3's module has only 65 residues at pLDDT ≥ 70. An experimental-structure route (PDB chains cut by the same module rule) for S12 | open (S12) |
| 2026-10-02 | S8b | **DEG/ENaC: ENaC nests inside the invertebrate degenerins** in the full-length tier-2 tree (deg_invertebrate not one clade; 10 ENaC tips in the way), so the declared outgroup is split and the tree is unrooted. Either the invertebrate degenerins are paraphyletic to ENaC (a real result about the outgroup choice) or a 7-tip outgroup sampling artefact; D48 is not revisited — S11 / S10 should state which before reading DEG/ENaC gene gains | open — S10/S11 |
| 2026-10-02 | S8b | **`scripts/selftest.py` is over the 500-line budget** (539 before S8b). S8b put its checks in `selftest_s8.py`; the file should be split by topic (catalogue / classifier / census / phylogeny / figures) | open — housekeeping |

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

**D33 — A family is called on a domain it carries and its rivals lack, or
not at all at the architecture tier.** D25 generalised, after S2b measured
what absence rules cost: H2 (6,633 receptor-length proteins called AChBP),
H4 (N-terminal RyR fragments called ITPR, 502 contradicting the IP3R
census), H13 (2,989 polycystin-1-like proteins called TRPP). A `forbid` list
may guard a rule, never be its basis. Where no positive domain exists —
ITPR against RyR, AChBP against a receptor fragment — the architecture tier
stops at the superfamily and a sequence-level tier (reference margin D7,
profile margin D32) makes the call. A candidate positive test is adopted
only if it agrees with an independent instrument: the AChBP shape test
(complete, soluble, 157–300 aa) was rejected at 53/1,263 agreement. Whole-
sequence measurements a rule reads (TM count, length, fragment flag) must be
*known*: unknown never passes as zero.

**D34 — Where the census looks is declared separately from what a signature
proves.** `Signature.enumerate` (default: SUPERFAMILY/FAMILY levels) decides
membership of `pore_signatures()`; `Signature.level` decides what the
classifier may conclude. Coupling them let S2b's evidence fix (PF20519 →
SHARED_WITH_DECOY) remove a signature from the census search space without
any enumeration being re-run. The self-test pins the union to the signatures
census v2 actually enumerated; changing the search space is a new census,
never a side effect.

**D35 — The denominator is the panel's reference proteomes in census v2's
release, chosen by rule, and a mid-release reissue is accepted only under
two checks.** One row per species in `src/utils/species.py` (52): the UniProt
reference proteome under the panel taxon that release 2026_03 ships, selected
by most Swiss-Prot entries → BUSCO → genes → UPID (the model strain, without
naming it); a species with none is kept as `genome_only` with its best NCBI
assembly, never dropped. N = 50 proteomes + 2 genome-only across M = 15 panel
groups. Every file must match its `RELEASE.metalink` MD5 (never waived) and
the release README's counts. The README count alone may be waived for a
file served after the release date whose entry count equals the proteome's
declared gene count. Measured case: human, reissued 2026-09-15, README 147,503
entries (the whole human UniProtKB set) against 20,652 served (one per gene,
like the other 49). The canonical FASTA is the gene count; gene2acc "groups"
are not (mouse 22,353 against 21,860 genes).

**D36 — A jackhmmer run killed by K3 contributes nothing.** The ported rule
kept a killed run's rounds before the firing, which for K3 is rounds 1–9 of a
run that never converged — contradicting K3's own statement that no
completeness claim may rest on it. S3b measured the cost: with K3 runs
included, 30,639 `jackhmmer_only` candidates, most from runs that ballooned
into repeat-domain proteins; without, 7,980. K1/K2 runs keep their pre-kill
rounds (the drift is dated), clean runs keep everything, and the
completeness claim is stated on clean runs alone. Changed before any
candidate was used, and the candidates are never called (D33).

**D37 — A genome absence is judged only against an in-group bait, and its
control is measured under that condition.** S5a measured that neither
miniprot nor tblastn detects a family member from a bait in another panel
group at 20–30 % identity (7 animal/ciliate TPC baits: 0 alignments, 0 HSPs
on *Arabidopsis* TPC1). So: (1) the bait panel takes **one bait per species**
(B2, no cap in practice); (2) a genome is never read with its own species'
baits (genome-scale leave-one-out, D29); (3) a zero cell can be `absent`
only if it is *matched* — a non-self bait of the family comes from the same
group — and the genome's **matched detection** (control cells detected by
miniprot *or* tblastn, since an absence needs the whole instrument silent)
is ≥ 0.90, and its N50 meets the family's median annotated gene span (D4);
(4) a locus is called by the S3a profiles (D32), never by its bait; a claim
that a proteome *missed* a gene needs a high-confidence call on an intact
reading frame; (5) annotation confirms a locus only where ≥ 50 % of its CDS
lies on the gene's exons. The 0.90 floor was set before any zero cell was
read and sits below every pilot genome (min 0.974).

**D38 — S5b's two thresholds, fixed before any S5b cell was read.**
(1) **`-G` stays 1 Mb for genomes ≥ 1 Gbp.** A 1.5 Mb margin was measured
on mouse, not assumed: high-confidence loci 320 → 311, annotation-confirmed
304 → 295, the real Hvcn1 locus lost, an ASIC locus chained into a 643 kb
chimera — and **no** call used an intron > 1 Mb (Asic2's 996,015 bp intron
is recovered intact at 1 Mb). A wider `-G` buys chaining across genes; an
intron beyond it splits a gene into partial loci rather than losing it, and
`introns.tsv` `over_G` reports every such annotated case. (2) **D4's bar per
cell** (`s5_verdict.d4_bar`): in a prokaryotic or viral genome, 3 bp × the
family's upper length band (no spliceosomal introns); otherwise the
family's median annotated span *in the cell's own group* when ≥ 3 of its
genes are annotated there, else the pooled median, else no bar —
`absent_bar_unmeasured`, never a guessed one (D28).

**D39 — A family alignment takes high-confidence profile calls only, fixed
before any alignment ran.** A census v4 row enters its family's tier-1
alignment only if the S3a profiles call it to that family at high
confidence (D32 margin ≥ 0.30 and ≥ half the profile covered); a genome
locus additionally needs an intact reading frame (D37 (4)). Medium calls
(the S3b ankyrin/LRR upper bound; sister families inside the margin),
S2-only calls with no profile call and broken genome frames are excluded
and counted per row in `results/alignments/members.tsv`; identical
sequences collapse to one, the duplicates listed. Every family is aligned
by the same method, MAFFT L-INS-i, whatever its size — none is moved to a
faster algorithm (D28). The profiles are S3a's frozen library, the
instrument that made the calls; the S2b R3 seed re-draw would change that
instrument and is therefore a census revision, not an S6 step.

**D40 — A tier-2 unit's modules are cut by one method, and how each
family's span was fixed is recorded.** `Superfamily.module_rule` declares
the module (`pore_loop`: the TM helix before each re-entrant pore loop
through the helix after it; `tm_span`: first TM to last). Every member is
`hmmalign`ed to its own family's S3a profile and cut at that family's span
of match states (`profile_projection`). The span comes from the UniProt
topology of the family's annotated references (29 families) or, where no
reference has an annotated pore loop (CNG, TRPC, TRPN, ITPR, plant and
prokaryotic iGluR), from a vote of a module HMM built from the annotated
references' modules — accepted because, run held-out on the annotated
P-loop and iGluR families, it lands within 12 profile states of the
annotated span; a vote boundary > 12 residues off every predicted helix
moves outward to the bracketing helix (fired once: ITPR, whose RyR-seeded
vote began in the luminal loop). A single module HMM as the extractor was
measured first and rejected (held-out Kir 22/298, TRPM 0/145, innexin clan
≤ 5/82).

**D41 — Tier-1 trees: trimming measured, roots declared, one command,
all fixed before any tree was read.** (1) **Trimming is trimAl `-gt 0.5`**
(keep a column where ≥ half the family has a residue) for every family,
chosen on alignment properties alone over four candidates
(`tier1_trim_compare.tsv`): median 504 informative sites against 324 for
S6's `-automated1`, ≥ automated1 in 62/63 families, a median 91 % of each
member's own residues kept against 61 %, gap fraction 7 %. The mask is
computed on the family alignment and applied unchanged to outgroup rows.
(2) **A family is rooted on the catalogue exemplars of its superfamily's
`root_with` family**, added by `mafft --localpair --maxiterate 1000 --add
--keeplength` (ingroup rows checked unchanged). The outgroup family itself
and families in superfamilies with no declared outgroup are built
unrooted and reported so — never midpoint-rooted. Whether the outgroup
forms one clade, and the UFBoot on that edge, is reported per tree; a
non-clade outgroup leaves the root undefined, not repaired. (3) **IQ-TREE 2
`-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -bnni -seed 1`** for every family:
the full ModelFinder set (~500 models) was measured at 20–40 s per model on
the largest families and restricted to four general empirical matrices
before any tree finished.

**D42 — A family added after the census is declared, not searched, until a
census revision brings it in.** S20 found eight proposed channels and five
auxiliary genes the catalogue lacked; they were added (user-directed,
2026-09-29) with every new signature `enumerate=False`, so census v2's
search space — pinned by the self-test (D34) — and every downstream count
(census v2–v4, S6 alignments, S7 trees, S15, S20's panel frame) stay on the
68-family / 320-gene catalogue they were built on. The catalogue's headline
(75 census families, 328 human genes) and the census's (68 / 320) therefore
differ until the census-revision emergent row is done, and every report says
which it uses. Where a new family shares its only domain with a non-channel
(TMCO1/EMC3, TMEM87A/GOST, TMEM109/BRI3BP) the decoy is catalogued with it and
the domain is SHARED_WITH_DECOY, so the architecture tier can never call the
family (D33); hazards H17–H19 are open until the profiles exist.

**D43 — A census revision for new families is a delta, and every old
result it touches is checked unchanged.** Census v2 r4 (S2d) brought the
families added after S20 into the census without re-walking or re-sweeping
what already existed: (1) enumeration walks only records carrying a new
signature and none of the old ones, per taxonomic shard, and must satisfy
old union + delta = new union and every per-signature count exactly;
(2) old records carrying a new signature are re-classified, and the
revision fails if any other call changes; (3) new profiles are built alone
and swept over the old database, every profile over the delta, with `-Z`
held at the old database size so no existing hit's E-value moves; (4) the
merge fails if an old call changes without a new family as its call, its
profile winner or the runner-up that pulled the winner inside the D7
margin. A new family's profile is evidence only once every family that
could score has a profile: the KChIP profile needed a neuronal-calcium-sensor
decoy before its calls meant anything (H20). Supersedes D42's "declared,
not searched" for these families; D42 stands for any future addition.

**D44 — The seed rules, corrected where they misread their inputs; the
frozen profiles keep the seeds they were built with.** (1) R3 and S3a's
held-out set B now read UniProt's flag correctly: "Fragment" excludes,
"Precursor" does not (S2b fixed the classifier; the seed rule and the
benchmark still read any flag as a fragment, which left CLCC1 with one
seed and kept 203 reviewed orthologues out of set B). (2) A catalogue gene's
family is always the catalogue's current one, never the family cached
beside its accession (TMEM87B moved families; the cache would have kept
seeding the old one — the seed-disjointness rule caught it). Profiles built
from 2026-09-30 on use the corrected rules; the 91 frozen profiles keep
their recorded seeds, and a full library rebuild would apply (1) to them.
(3) A decoy mixing two clades cannot draw a paralogue line: TMEM87B has
its own decoy apart from GPR107/108 (H18, measured 0.23 → 0.60).

**D45 — A family whose paralogue split is lineage-specific is counted by
descent.** TMEM87A and TMEM87B are vertebrate paralogues; an invertebrate's
single TMEM87 is co-orthologous to both, so a family boundary at the split
(TMEM87A a channel family, TMEM87B a decoy) left every non-mammalian TMEM87
unclassifiable (measured: the *Xenopus* orthologue fell to the decoy; the
genome matrix showed TMEM87 scattered). The family is TMEM87 (A and B) by
descent, `channel_contested`, with the mechanism — channel activity reported
for TMEM87A only — a literature note (D24). This differs from the CLC and
anoctamin channel/non-channel splits, which are old enough that every
lineage's genes fall clearly on one side. Counting by descent also turned
H18 into a positive test (`PF21901`, D33). User decision, 2026-09-30.
*Amended the same day*: the GOST tree shows TMEM87 is a pan-eukaryotic
lineage (UFBoot 100) whose GOLD domain is an animal addition, so the family
spans eukaryotes (PF21901 marks it only in animals; the profile margin marks
it elsewhere). The tree does **not** resolve TMEM87A or TMEM87B as clades,
so "the split is a vertebrate duplication" is unverified; D45 rests on the
measured fact that no sequence method separated non-mammalian TMEM87 into A
and B.

**D46 — S10 reads repertoire on a dense proteome panel, and states
absences at two strengths.** The user chose one UniProt reference proteome
per eukaryotic order (439, release 2026_03, S4b) after the options were
measured (`results/panel_density/`). Presence/absence on it is a proteome
call (the S3a profiles, D32); only the 52 S4 species carry S5's controlled
genome absences (D37/D38). S10 reconstructs gains and losses on the dense
panel and states any absence as either *controlled* (S5) or *proteome-only*
(S4b), never merging the two. Orders are NCBI order ranks; lineages without
one (86 proteomes, mostly protists) are outside the rule, a stated gap.

**D47 — A weak root is re-tested under two declared outgroups, and kept
only if it does not move (S7c).** For the six P-loop families S7b rooted
badly on KcsA/MthK/NaK (CNG, K2P, Kv-KCNQ, Kv-Shaker, Slo, Kv-EAG), fixed
before any re-rooted tree was built: (1) **two outgroups per family, declared
in the catalogue** (`ChannelFamily.root_with`, overriding the superfamily's
D41 declaration): a prokaryotic relative with the family's architecture —
KvAP + MVP + *E. coli* Kch (reviewed six-helix K⁺ channels) for the Kv
families and Slo, MloK1 + SthK for CNG, the superfamily's KcsA/MthK/NaK for
K2P (no prokaryotic K2P exists) — and two sister eukaryotic families (their
catalogue exemplars); inline exemplars verified live (`s0_outgroups.py`) and
kept out of the classifier's reference panel. (2) **Denser basal
sampling**: one S4b member per early-diverging clade (kingdom × phylum ×
class, outside Bilateria), the clade's top-scoring high-confidence, complete,
in-band profile call — so these six families' sets are S6's D39 sets plus
the picks, a stated exception to D39. (3) S7's method unchanged after
re-alignment (L-INS-i → `-gt 0.5` mask on the ingroup → IQ-TREE `IQ_ARGS`);
a family's two trees share the ingroup alignment and mask. (4) **A root is
accepted only if** in both trees the outgroup is one clade with UFBoot ≥ 95
on its edge **and** the ingroup's basal split is identical under the two
outgroups; otherwise `unresolved`, never picked. The non-reversible-model
root and reconciliation rooting are S11's.

**D48 — Tier-2 trees and the tier-3 network: units, representatives, roots
and edge verdicts fixed before any tier-2 tree or structure comparison
(S8a).** (1) **Units**: alignable superfamilies with ≥ 2 census families;
module superfamilies take S6's `full` modules, one tip per module (repeat);
the rest full length; single-family superfamilies are their tier-1 tree;
D27's refusals go to the fold network. (2) **Representatives (D8)**: per
family × panel group × module, members ordered by centrality (mean
within-family identity: module match states of S6's `hmmalign`, or the
family L-INS-i columns; ≥ 50 % mutual coverage) and greedily clustered at
0.5 identity — unless that leaves > 4 tips per informative site of the
trimmed unit alignment, when the highest threshold meeting the cap is used
(P-loop: 0.3, 317 tips on 85–90 sites; decided on alignment properties, no
tree seen). (3) **Root**: the superfamily's `root_with` family, its tips
restricted to the kingdom of its catalogue exemplars, plus those exemplars
as tips (module units cut by S6's method); no declaration → unrooted.
(4) MAFFT L-INS-i, trimAl `-gt 0.5`, S7's `IQ_ARGS`. (5) **Fold network**:
node = census family, the first S0 exemplar with an AFDB model of that
exact accession (none → unmeasured, no substitute); unit = module 1 for
module superfamilies, whole model otherwise, pLDDT ≥ 70; extra node Kv
VSD (UniProt S1–S4); TM-align average-length TM-score, Foldseek beside it;
an edge is `supported` iff its median family-pair TM-score ≥ 0.5 **and**
each side's best other superfamily is the other side.

**D28 — A missing tool disables a test, loudly.** `MafftUnavailable` is
raised, not caught; a missing IQ-TREE writes the alignment and no tree, with
the reason in `TreeRun.note`. No silent fallback to a weaker method, ever —
two incomparable methods in one figure is the failure this prevents.
