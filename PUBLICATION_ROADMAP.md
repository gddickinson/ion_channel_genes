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

> **Next session: S5b** — `python3 scripts/s5_sweep.py run --all` (52 genomes, ~42 Gbp; *Torpedo* takes the chunked path), then `s5_calibrate.py`, `s5_verdict.py`, `s5_report.py`, and census v4. Read `results/genome_sweep/report.md` first. **Decide before the run** (emergent rows): `-G` for ≥ 1 Gbp genomes (Asic2 996 kb sits 4 kb under 1 Mb) and a D4 bar rule for the families the pilot could not measure. S5a completed 2026-09-29 (S5 split: S5a instrument + pilot, S5b full sweep). Previously: **S5** — genomic tblastn + miniprot for families S3b finds absent in a lineage, plus the two genome-only species (*Cornu*, *Torpedo*), so an absence passes D4's two bars. S3b completed 2026-09-29: read `results/panel_sweep/report.md`. **What S3b hands S5**: `family_by_species.tsv` (census families × 50 proteome species) is the presence matrix whose zero cells S5 must test; the high-confidence calls are the evidence, medium calls in repeat-dominated profiles (TRPN, LRRC8, TRPA/C/V) are not. **Previously S3b** — sweep the S3a profiles over the S4 panel DB (`<data root>/proteomes/s4/panel_refprot.fasta`, 822,499 sequences) + jackhmmer (D10). S4 completed 2026-09-28: read `results/proteome_scope/report.md` — the denominator is 52 species, 50 reference proteomes + 2 genome-only (D35). **Previously S4** (S3b needs S4's proteome manifest). S3a2 (ABC decoy family) completed 2026-09-28. S2b and S2c (user-directed) completed 2026-09-28: every hazard rule is a positive test (D33), census v2 is r3, v3a re-merged. An emergent row proposes a generic ABC-transporter decoy family for the S3a profile library. S3a completed 2026-09-28 — read `results/census_v3/report.md`; the emergent rows it added (H2/H4/H13 absence rules, the PF00520 superfamily conflicts) are S2 rule fixes that can be done before or alongside S4.
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
| S3b | Profile sweep of the S4 reference proteomes (what domain search misses) + jackhmmer-to-convergence completeness argument (D10) → census v3 | S3a, S4 | completed 2026-09-29 | **Census v3 over the declared denominator: 49,595 of 822,499 panel entries carry evidence; 8,755 census-family calls, of which 230 (2.6 %) are outside census v2 — only 96 at high confidence.** S3a's 91 profiles unchanged, swept over S4's `panel_refprot.fasta` (50 proteomes, `-Z` 822,499; 12.4 min wall), D32 assignment imported: 16,183 family calls (11,798 high). **Instrument check: same verdict as S3a on 12,072 / 12,096 (99.8 %)** of the panel entries that are also census v2 records (only `-Z` differs). **Human 319/320 right family** (GLRA4 not in the proteome). Where domain search misses channels (high-confidence): **CLIC 27**, **Hv1 15** (molluscs, moss, fish — 25/33 panel Hv1 carry no enumerated signature), CNG 10, TRPV 8, pannexin 8, and **all 3 viroporins** (Vpu, M2, SARS-CoV-2 E — the S4 emergent row confirmed). By group the rate runs from 1.0 % (vertebrate) to 12.8 % (plant). **The medium 134 are an upper bound, not a finding**: TRPN's 37 medium calls are ankyrin-repeat proteins (ANKRD52, ankyrin, IκB) and LRRC8's 20 are LRR proteins — repeat-dominated profiles let the module alone pass D32's 30 % coverage gate (emergent). **jackhmmer: 68/68 families seeded by rule** (top high-confidence panel call), `-N 10 -E 1e-5`, 26.9 run-hours (6.1 h wall). **D10: 24 clean, 44 killed (K1 23, K2 2, K3 19).** K1 fires at round 2 across P-loop/iGluR/TRP: a family-seeded model takes in its sister families in one iteration, while every Cys-loop, DEG/ENaC and P2X family run converges clean on the *same* superfamily set (1,625 / 490 / 125) — at those scales jackhmmer is a superfamily instrument. **Completeness (clean runs): 3,352 / 3,356 (99.9 %) of profile calls recovered**, 1,460 accepted targets no profile calls; all runs 8,610 / 8,736 (98.6 %). K3 runs contribute no candidates (**D36**). 7,980 `jackhmmer_only` candidates, dominated by the TRPA/TRPN pre-K1 rounds (ankyrin drift K1 cannot see) — counted, never called. Bulk on the data root under `hmmer/s3b/`. → `results/panel_sweep/report.md` |
| S2b | **Hazard-rule rewrite (H2, H4, H13) as positive tests** — the three S2 rules that call a family from an *absence*, found by S3a's calibration and external check; re-benchmark S1, re-classify the affected S2 records, rebuild census v2 and re-merge v3a | S3a | completed 2026-09-28 | **All three rules now call a family only on a domain it carries and its rivals lack; where none exists the architecture tier stops at the superfamily (D33).** H4: `PF02026` or `PF06459` with the shared core ⇒ RYR; **ITPR has no positive architectural test** (no domain RyR lacks) → superfamily, family from reference/profile. H13: PLAT/REJ/GPS or `PF00801`×≥5 ⇒ polycystin-1; `PF18109` ⇒ TRPP (measured: 1,358 records, all profile-TRPP); channel domain alone ⇒ superfamily; catalogue fix `trpp` `PF20519` FAMILY → SHARED_WITH_DECOY. H2: **no positive AChBP test exists at this tier** — the best candidate (complete, soluble, AChBP-sized LBD) agreed with the profiles on 53/1,263 and put ~75 % of its 2,164 calls in Ecdysozoa/Chordata, which have no known AChBP → LBD without TM is superfamily-only. **Re-classified 108,407 records** (`s2_classify.py --recheck`), **27,028 calls changed, 0 outside the re-checked set**; census v2 is now **r2** (r1 archived). S2 r2: family 324,769 (26.1 %), human 169/320 right (ITPR1–3, PKD2L2 now superfamily-only) and **0 wrong** (ZACN no longer AChBP). **S1 benchmark unchanged**: recall 50/72, specificity 25/25 — ITPR1 now called by reference identity, PKD1 now positively by PLAT/REJ. **Census v3a re-merged: conflicts 12,857 → 3,288; ITPR/RYR swaps vs the IP3R census 502 → 0**; family calls 732,191 (58.8 %); human **319/320** (ZACN now `zac` by profile). **Two bugs found on the way**: `tm_count` 0 was passed as unknown (`x if x else None`) in both `s2_classify.py` and `cli_channel.py`; UniProt's `flag` "Precursor" was read as a fragment. Both fixed, both self-tested. → `results/census_v2/report.md` § Revision r2, `results/census_v2/s2b_transitions.tsv` |
| S2c | **Hazard-rule rewrite (H11, H12) as positive tests** (D33) — SUR vs CFTR, KCTD vs Kv; re-benchmark S1, re-classify affected S2 records (census v2 r3), re-merge v3a | S2b | completed 2026-09-28 | **No hazard rule now calls a family from an absence** (self-test invariant over every hazard). **H12**: KCTD called on one of five KCTD C-terminal domains (catalogued as `nonchannel_kctd` FAMILY signatures; measured on 14,354 records, 14,317 profile-KCTD, 3 otherwise); T1 without a pore module is superfamily-only — the shape test (T1, 0 TM, complete) was rejected, contradicted by the profiles on 160 records. **H11**: CFTR still called on its R domain; **the SUR rule is removed — no domain identifies SUR** (ABCC8 carries only generic ABC domains; TMD0 is on the MRPs too), and none of its 980 census calls was a SUR (918 bacterial cNMP + C39 peptidase ABC transporters, 62 eukaryotic fused gene models); `assoc_sur` ABC domains re-levelled SHARED_WITH_DECOY. **Re-classified 63,411 records, 18,412 calls changed, 0 outside the re-checked set** (KCTD → [ploop] 17,423; SUR → unassigned 944). Census v2 **r3**: family 306,378, human 169/320 right, 0 wrong. **S1 unchanged** (50/72, 25/25, 16/16) — ABCC8 now uncalled, still rejected; KCTD1 still KCTD, by its C-terminal domain. **Census v3a re-merged: conflicts 3,288 → 2,907** (none from any hazard rule's absence logic; what remains is the PF00520 superfamily problem and K2P copy-number calls); KCTD agreement 99.98 %, CFTR 100 %; human 319/320. **Found on the way**: S2b had silently dropped PF20519 from the census search space (enumeration coupled to evidence level) — fixed (D34) and pinned by a self-test; and the S3a profile library's missing ABC-transporter decoy lets the SUR profile call 525 records (508 bacterial) — emergent, a catalogue addition. → `results/census_v2/report.md` § Revision r3 |
| S3a2 | **ABC-transporter decoy family** for the S3a profile library (user-directed, from the S2c emergent row): catalogue `nonchannel_abc_transporter`, its profile built and swept alone (the other 90 frozen), S0 verification, benchmark re-run, v3a re-merged | S2c | completed 2026-09-28 | **The SUR profile no longer wins anything that is not a SUR: profile-called SUR records 525 → 0.** New control family `nonchannel_abc_transporter` (91 families, 23 controls; census families, human genes and the 67-signature search space unchanged): human ABCC1–6/10–12 + *E. coli* HlyB, *B. subtilis* SunT (C39-peptidase exporters), *E. coli* MsbA; ABC domains SHARED_WITH_DECOY, `PF03412` ACCESSORY. **S0 re-run clean** — 793 requests, 0 failures, 166/166 exemplars, 121/121 Pfam, 52/52 taxa; **no existing exemplar sequence changed**, so the 90 frozen profiles still match their seeds. Profile: 12 seeds → 1,304 states, built alone (`--only`); the other 90 builds and seed sets byte-identical. Swept in 21 s (2,919 targets). **S3a benchmark unchanged** (LOO 69/71, orthologues 534/534, 0 decoys called channel); ABCC8 still SUR by 2,478 vs 1,147 bits, CFTR by 3,237 vs 971. **Census v3a**: the decoy takes 1,026 records (997 bacterial); conflicts 2,919 (+12, all ABC-domain/channel-domain fusions, S2 calling the channel domain and the profile the ABC part); human 319/320. → `results/census_v3/report.md` |
| S4 | Proteome scope: declared reference-proteome manifest across the lineage panel (the denominator) + download tooling | S1 | completed 2026-09-28 | **The denominator: 52 panel species in 15 groups → 50 UniProt reference proteomes (47 cellular, 3 viral) + 2 genome-only (*Cornu aspersum*, *Torpedo marmorata*: chromosome-level, unannotated NCBI assemblies), 0 with neither** — release **2026_03**, the one census v2 was enumerated from (D35). **822,499 canonical entries, one per gene** (entry count = declared gene count 50/50). Selection by a stated rule (most Swiss-Prot entries → BUSCO → genes → UPID), which picked the model strain in all 6 multi-candidate species (HB8, S288C, Nipponbare, 3D7, PR8, HXB2). **Verified: 49/50 exact** against the release README *and* the per-proteome metalink MD5; **human was reissued by UniProt mid-release** (served 2026-09-15; README #(1) 147,503 = the whole human UniProtKB set, served file 20,652 = one per gene, MD5 matches) — accepted under D35. **Positive control 319/319** enumerated human census genes are human-proteome entries. 12,096 census v2 records are panel-proteome entries. Per row (D9): assembly accession + level + N50, annotation source, BUSCO, CPD; **7 proteomes sit on a superseded assembly version**, 14 below 90 % BUSCO, 14 scaffold-level — reported, not applied (S5's D4 bars decide). **Found on the way: the viroporin family is outside census v2's search space** (its 4 signatures are `SUBFAMILY`, so never enumerated; 0 census records in the 3 viral proteomes) — emergent. Sweep DB 507 MB on the data root (the ~25 GB estimate assumed a much wider set). Offline rerun reproduces all tables byte-identically. → `results/proteome_scope/report.md` |
| S5a | **Genomic sweep: bait panel, instrument, pilot** — the instrument for S5, measured rather than ported: bait rules, locus → profile call (D32), tblastn rescue, matched positive control, D4 bar from annotation, seven-genome pilot | S3b, S4 | completed 2026-09-29 | **The instrument is measured, and its sensitivity is a property of the nearest bait.** Bait panel: **1,808 baits for all 91 catalogue families** (controls included) from 50 panel species, 1.63 M residues, by rules B1–B3 (`s5_baits.py`). Pilot: **7 genomes, 6.02 Gbp, 0 failures** (*E. coli*, yeast, *Arabidopsis*, *C. elegans*, *Drosophila*, mouse, genome-only *Cornu*), 1,881 loci, **1,020 called to a family by the S3a profiles** (D32 on miniprot translations, `-Z` = panel size; the bait family is a second witness), 29.8 min wall. **The headline is D37.** With 8 baits per family spread one per group, *Arabidopsis* recovered **5 of 9** of its own proteome's high-confidence families once its own baits were excluded — **AtTPC1 drew 0 miniprot alignments and 0 tblastn HSPs (E ≤ 1e-5) from 7 animal/ciliate TPC baits**. Neither pairwise method reaches a relative from another group at 20–30 % identity. One bait per species: 8 of 9. So the control is **matched**: detected / control cells that have a non-self bait from the same group — the condition every informative zero cell is judged under. **Matched detection 172 / 173** (the miss: *C. elegans* tweety); found (a profile call) 166 / 175; every genome ≥ 0.974 against a 0.90 floor. *Cornu*'s control is its group's 27 core families: 27 / 27. **Calibrated from annotation, not miniprot:** 420 called loci confirmed on annotated exons. **Mouse census genes carry introns far past miniprot's 200 kb default — 13 / 220, up to 996 kb (Asic2), Trpm3 558 kb, Kcnd2 497 kb** — so `-G` 1 Mb for genomes ≥ 1 Gbp is necessary and Asic2 sits 4 kb under it (emergent). D4 span bar per family = median annotated span, **61 families measured**. **No identity floor**: annotated loci run down to 0.21 identity, so the parent's 0.40 does not transfer. **Three instrument corrections found by the pilot**: "overlaps an annotated gene" confirmed a 400 kb chained ZAC model in mouse and a TRPM model on Ugt1a10 → confirmation is now CDS-on-exons ≥ 50 % (a reciprocal-extent rule was tried and rejected: mouse 223 → 130, real genes with long UTR introns); mouse carries **12 high-confidence VDAC loci for 3 genes**, most with frameshifts/stops → a proteome-miss claim needs a high-confidence **intact** locus (`genome_weak` otherwise); genome-only presence is `genome_present`, not a miss. **Zero cells (312 in the pilot; 71 informative, 46 of them *Cornu*)**: of the other 25 informative, **10 absent**, 5 trace, 3 absent with no measured bar, 5 partial, 1 gap, 1 weak — and **0 `genome_found`: no pilot proteome misses a high-confidence intact channel gene**. **Literature checks 4 / 4 consistent** (refs pending): *C. elegans* Nav and P2X, *Drosophila* P2X `absent`; mouse ZAC `genome_weak` (the 400 kb chained model, 32 % identity, off the syntenic chromosome). Self-test +8 invariants, the exon rule mutation-tested. Budget for S5b: ~225–260 s/Gbp miniprot on 1,808 baits → ~3 h over 42 Gbp. → `results/genome_sweep/report.md` |
| S5b | **The full panel sweep** (52 genomes incl. *Cornu*, *Torpedo*) → per-cell ledger (found / partial / gap / trace / absent) → census v4 | S5a | pending | |
| S6 | Alignment upgrade: MAFFT L-INS-i + trimAl per family, and pore-module extraction for every tier-2 unit | S2, S3a, S3b, S5b | pending | |
| S7 | **Tier-1 phylogenies**: IQ-TREE 2 + ModelFinder + 1000 UFBoot for every census family, rooted on its declared outgroup | S6 | pending | |
| S8 | **Tier-2 pore-module phylogenies** for every alignable superfamily, plus the **tier-3 fold network** (Foldseek/TM-score) for the refused ones | S6, S7 | pending | |
| S9 | Selectivity-filter atlas: the filter locus for every P-loop family, and whether it is congruent with the tier-2 tree (Q5) | S7, S8 | pending | |
| S10 | Repertoire evolution: presence/absence of every family across the lineage panel; ancestral reconstruction; the MscS animal-absence question (Q4) | S5b, S8 | pending | |
| S11 | Duplication history: which families expanded in which lineage, 2R/3R ohnologue status, and the Kv/Nav/Cav repeat-duplication order | S7, S10 | pending | |
| S12 | Structures: AFDB coverage per family, TM-align against the experimental reference set, Foldseek all-vs-all for the fold network | S2, S6 | pending | |
| S13 | ML selection tests on the families with clinical variant sets (CFTR, SCN1A, KCNQ1, RYR1) | S6, S7 | pending | |
| S14a | **Manuscript assembly** — draft, figures, methods, deposit manifest, reviewer self-audit | S3a, S3b, S5b, S7, S8, S9, S10, S12, S15–S19 | pending | |
| S24 | Supplementary alignment + structure figures, and a figure-by-figure audit | S14a | pending | |
| S14b | **Deposit + release** — Zenodo DOI, repo public (D2 flip), reference verification, preprint upload | S14a | pending | **Human-gated; cannot be completed autonomously.** |

### Analysis & synthesis block (S15–S22)

These run on data the earlier tasks already produced. Priority orders them
when several are unblocked at once.

| ID | Task (one session each) | Depends | Priority | Status | Results |
|----|-------------------------|---------|----------|--------|---------|
| S15 | **Method contribution** — what would domain search alone have missed, per superfamily (Q3); the recall curve as methods are added | S2, S3a, S3b, S5b | high | pending | |
| S16 | **Annotation-quality audit** — how often a real channel locus is missing, fragmentary, split, unnamed or filed under the wrong family across RefSeq / Ensembl / UniProt / InterPro; the correction list | S5b, S15 | high | pending | |
| S17 | **Mechanism vs clade** — do the CLC channel/transporter and anoctamin channel/scramblase splits survive the tier-1 trees (Q7) | S7 | high | pending | |
| S18 | **Convergence audit** — connexin vs pannexin, TMEM175 vs the GYG channels, the mechanosensitive families: is "unrelated" a fact or a detection limit (Q6) | S8, S12 | high | pending | |
| S19 | **Channelopathy map** — clinical variants of every family mapped onto the pore module and the constrained core | S6, S12 | medium | pending | |
| S20 | **Auxiliary-subunit census** — the excluded 15 %: how many there are, and how often they are counted as channels in published totals | S2 | medium | pending | |
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
| 2026-09-28 | S2b | **H11-abcc and H12-kctd are still absence rules** (ABC architecture *without* the R domain ⇒ SUR; T1 *without* a pore module ⇒ KCTD). H12 accounts for 358 of the 3,288 v3a conflicts (profile calls Kv). Rewrite as positive tests the way S2b did H2/H4/H13; the S2b self-test invariant already names them as exempt | closed 2026-09-28 (S2c) |
| 2026-09-28 | S2b | **Census v2's `fragment` column stores UniProt's raw `flag`, which includes "Precursor" (2,203 records).** S3a's R3 seed rule and benchmark set B excluded those as fragments — conservative, not wrong, but mislabelled. S2b's classifier path now tests for "Fragment"; the column and the S3a filters should follow at the next rebuild | open |
| 2026-09-28 | S2b | **S3a's R3 seeds were drawn from S2 r1 calls**, six ITPR and five TRPP R3 seeds were called by the now-removed H4 / H13 absence rules (all reviewed, non-fragment, inside the family length band ±25 %, so not the fragments and polycystin-1-like proteins those rules miscalled). Profiles were kept, not rebuilt, so the S2b before/after is a comparison on one instrument. A rebuild at S6 should re-draw R3 from the current census | open (S6) |
| 2026-09-28 | S2c | **The S3a profile library has no decoy for generic ABC transporters, so 525 records (508 bacterial) are called SUR by the profile alone.** They are ~1,000-aa ABC transporters with cNMP + C39 peptidase domains; the SUR profile (2 human seeds) is the only ABC-transporter profile they can score against. Profile agreement is evidence only when every family that could score has a profile. Fix: a catalogued non-channel ABC-transporter decoy family (human ABCC1–6/10–12 + a bacterial peptidase-containing ABC transporter as exemplars), its profile swept (~one profile), re-merge. Adds a 91st family — the user's call | closed 2026-09-28 (S3a2) |
| 2026-09-28 | S2c | **S2b silently shrank the census search space**: re-levelling PF20519 to SHARED_WITH_DECOY dropped it from `pore_signatures()` (67 → 66), because enumeration was derived from evidence level. Fixed by `Signature.enumerate` (D34) and a self-test pinning the union to the 67 signatures census v2 enumerated. No census data were affected (nothing was re-enumerated) | closed 2026-09-28 |
| 2026-09-28 | S4 | **The viroporin family is outside the census search space.** Its four signatures (`PF00599` Flu_M2, `PF00558` Vpu, `PF02723` CoV_E, `PF11289`) are declared `SUBFAMILY`, and `Signature.enumerate` defaults to SUPERFAMILY/FAMILY (D34), so census v2 never enumerated them: 177 viral records in the census, **0** in the three panel viral proteomes, and M2 (P06821) absent. It is the only census family with no enumerated signature (checked over all 68). Fix: `enumerate=True` on those four — a change to the search space and so a new census revision (D34), not an edit. S3b's viroporin profile covers the panel meanwhile | open |
| 2026-09-28 | S4 | **The panel is 52 named species, a thin denominator for repertoire evolution.** It suffices for S3b/S5's completeness argument on the declared panel, but S10's gain/loss reconstruction and any "absent from lineage X" claim beyond these species needs density — e.g. one reference proteome per order by the same selection rule (the parent project swept 6,928 eukaryotic proteomes). A scope decision for the user before S10, not S4's | open (user, before S10) |
| 2026-09-28 | S4 | `INTERFACE.md` said the species table had 51 species; it has 52. Corrected | closed 2026-09-28 |
| 2026-09-29 | S3b | **D32's coverage gate does not exclude a shared module when the module dominates the profile.** TRPN's profile is mostly ~29 ankyrin repeats; 37 of its 38 panel calls outside census v2 are ankyrin-repeat proteins (ANKRD52, ankyrin-2, IκB, titin, synphilin) at medium confidence, coverage 0.30–0.46; LRRC8 (LRR half) 20 medium, TRPA/TRPC/TRPV similar. Fix candidates: require coverage of the *pore* segment of the profile (profile coordinates of the TM/pore region from the seed alignment), or re-level those calls to `module` when the covered span lies inside the repeat array. Must be benchmarked on S1/S3a before adoption; v3a is affected wherever S2 abstained | open |
| 2026-09-29 | S3b | **K1 is blind to drift into proteins no profile calls** (ankyrin/LRR/GST-fold repeats): TRPA/TRPN pre-kill rounds hold ~6,400 `module` targets, LRRC8 and CLIC K3 runs balloon to 9,025 and 1,565 accepted with 83 % / 92 % uncalled. The parent project's S19 measured the same blind spot and proposed a drift rule on the uncalled share of the finished model. Evaluate that rule on `jackhmmer_rounds.tsv` (`uncalled_frac` is recorded per round) as a classifier before adopting it — never retune D10 on the run it will judge | open (S19-like methods task) |
| 2026-09-29 | S3b | **In multi-family superfamilies a family-seeded jackhmmer is a superfamily search.** K1 fired at round 2 for 15 runs, 12 of them P-loop or iGluR; the Cys-loop, DEG/ENaC and P2X runs all converge on one set per superfamily. A superfamily-seeded design (one run per alignable superfamily, K1 on the *other-superfamily* share, recorded per round as `other_sf_frac`) would give a cleaner completeness argument for S15. Design note, not adopted | open (S15) |
| 2026-09-29 | S3b | **Hv1 is the family domain search misses most among real channels**: 25 of 33 panel Hv1 (molluscs, moss, fish, cnidarian, placozoan) carry none of the 67 enumerated signatures, 15 at high confidence — consistent with the S0 finding that the Hv1 domain `PF16799` is mammal-only. CLIC likewise (27 high, plants/ciliates/invertebrates). Both argue for profile-based enumeration in the next census revision | open |
| 2026-09-29 | S5a | **`-G` has no margin at the top.** Mouse Asic2's widest annotated intron is 996,015 bp against the 1 Mb `-G` for ≥ 1 Gbp genomes; 13 / 220 confirmed mouse loci have an intron > 200 kb. Before S5b: set `-G` from the widest measured intron with a stated margin (e.g. 1.5×) and measure what a larger `-G` does to chaining (the ZAC and Ugt1a10/TRPM chained models appeared at 1 Mb) | open (S5b) |
| 2026-09-29 | S5a | **D4 bar missing for families the pilot never confirmed on annotation** (3 informative cells: `plgic_prok`, `iglur_prok`, `deg_invertebrate`); prokaryotic and viral genes are intronless, so their bar could be the family's length band × 3 bp. S5b measures more; state the fallback rule before reading any cell | open (S5b) |
| 2026-09-29 | S5a | **Genome copy number is not locus count.** Mouse: 12 high-confidence VDAC loci for 3 genes, most with frameshifts or in-frame stops (retrocopies). Any genome-derived copy number (S11) must count intact, exon-confirmed loci | open (S11) |
| 2026-09-29 | S5a | **Partial loci at non-informative cells are shared modules, not misses** (e.g. yeast: 13 families partial — Kv, SK, CatSper, TPC, HCN, TRPs — from VSD/cNMP/CaM-binding modules). S16 must not read `partial` as an annotation failure | open (S16) |
| 2026-09-29 | S5a | **Unmatched cells cannot be judged by this instrument at all**: a family present in one panel group only can never be declared absent in another (no in-group bait). A profile-HMM scan of six-frame-translated genomes would reach further than pairwise baits; a method design for S15, and another argument for the panel-density row (S4) | open (S15) |

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

**D28 — A missing tool disables a test, loudly.** `MafftUnavailable` is
raised, not caught; a missing IQ-TREE writes the alignment and no tree, with
the reason in `TreeRun.note`. No silent fallback to a weaker method, ever —
two incomparable methods in one figure is the failure this prevents.
