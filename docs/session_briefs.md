# Session briefs — what each ledger task actually does

One brief per ledger row: goal, steps, completion criteria, outputs. Read the
brief for the task you are about to work, then work only that task.

The general shape of every task is the same: **produce tables, render the
report from the tables, and record what the task could not decide.** A task
that finishes with a conclusion and no table has not finished.

---

## S0 — Catalogue verification + scope confirmation

**Goal.** Turn `src/catalogue/` from a file full of assertions into a set of
checked claims, and build the reference panel the classifier needs.

**Steps.**
1. `python3 -m src.catalogue` and `python3 scripts/selftest.py` — both clean.
2. `python3 scripts/s0_catalogue_verify.py` (~40 min, ~800 live requests).
3. `python3 scripts/s0_report.py` and `python3 scripts/s0_figures.py`.
3a. The review and its figures:
   `python3 scripts/s0_review_refs.py` (resolve every citation),
   `python3 scripts/s0_domain_map.py`, `python3 scripts/s0_filter_atlas.py`,
   then `s0_review_fig_filters.py`, `s0_review_fig_domains.py`,
   `s0_review_fig_folds.py`, then
   `python3 scripts/s0_review_build.py --pdf`. **Look at every figure you
   regenerate** (D11).
4. Read `results/s0_baseline/report.md` end to end. For every flagged row
   decide: is the catalogue wrong, or is the database? Fix the catalogue,
   or add a note to the family's `notes` saying why the mismatch is
   expected. **Do not "fix" a mismatch by deleting the claim.**
5. Confirm the six scope answers in `docs/scope_and_boundaries.md` still
   describe what the catalogue does.

**Completion criteria.** `validate()` clean; every Pfam accession either
verified or removed with a reason; every exemplar resolved or marked
unresolvable in the table; `reference_panel.fasta` written; `report.md`
rendered from the tables only.

**A note on iterating.** S0 is expected to take more than one pass. The
setup session ran it three times: the first found two Pfam short-name
mismatches and thirteen unresolvable exemplars, the second found that a
"declared accession names the wrong gene" flag was really "UniProt has no
gene name for this entry", and only the third came back clean. Each pass
changed the catalogue, so each needed a fresh run — a report rendered from
tables that predate the edit is exactly the drift D13 exists to prevent.

**Outputs.** `results/s0_baseline/{catalogue_validation,pfam_verification,
exemplars_resolved,exemplar_architecture,signature_sharing,taxonomy_check}.tsv`,
`reference_panel.fasta`, `verification_summary.json`, `report.md`.

---

## S1 — Classifier control benchmark

**Goal.** Measure the classifier where it is meant to be hard, and find out
which tier is doing the work.

**Steps.**
1. `python3 scripts/s1_toolchain.py` — record exact tool versions first;
   a benchmark whose toolchain is unrecorded is not reproducible.
2. `python3 scripts/s1_benchmark.py` (~15 min).
3. `python3 scripts/s1_report.py`.
4. Read the per-tier table. **If most correct calls came from the reference
   tier, that is the headline finding**, not a footnote: the classifier is
   then a nearest-neighbour lookup and will fail on the sequences a census
   exists for.
4a. Check the leave-one-out line. Nearly the whole panel is drawn from the
   catalogue's own exemplars — 93 of 97 on the setup run — so the benchmark
   classifies with the query's own accession removed from the reference set.
   If that line ever says otherwise, the recall number is measuring the
   panel's overlap with itself and means nothing.
5. Read the hazard table. Every zero row is an untested hazard. Either
   extend the panel to exercise it, or record it in Emergent tasks.
6. Every wrong call gets diagnosed from its evidence string in `calls.tsv`
   and either fixed (a new hazard rule) or recorded.

**Completion criteria.** Recall and specificity reported with numerators and
denominators; per-tier attribution table present; every hazard either
exercised or listed as untested; `verify_anchor()` confirmed to pass, or the
motif tier reported as disabled.

**Outputs.** `results/benchmark_controls/{calls,confusion,recall,hazards,
coverage}.tsv`, `summary.json`, `report.md`, `results/toolchain_manifest.txt`.

---

## S2 — Uncapped enumeration → census v2

**Goal.** Every UniProt protein carrying any pore signature in the catalogue,
with a positive family call on every record.

**Steps.**
1. Enumerate from `registry.pore_signatures()` — **not** from `PF00520`
   (hazard H7).
2. `list_proteins_with_pfam(..., max_results=None, strict=True, dump_dir=…)`
   so pagination is complete and archived; record the API's own `count` per
   signature and the fetched total. A run that silently stopped early is the
   failure mode this guards.
3. Classify every record. Records the classifier cannot call are
   `unassigned` and stay in the census with that label.
4. Report per-family recall against the S1 panel, and the fraction of
   records the architecture tier could call without a reference.

**Completion criteria.** Fetched == API count for every signature, or the
shortfall is explained per signature; every record carries a call and a
confidence; the census is written as one FASTA + one CSV with provenance
columns.

---

## S3a — Profile library + best-profile assignment → census v3a

**Goal.** Resolve what S2's architecture tier could not (the 32 % that stop
at the superfamily) with one profile per catalogue family and a best-profile
call with a margin (D7) — and measure the instrument before using it.

**Steps.** (1) Seeds by stated, enforced rules (`s3_seed_spec.py`: R1
curated human genes, R2 exemplars, R3 reviewed S2 family calls one per
species); MAFFT L-INS-i single-threaded; hmmbuild — every catalogue family,
controls included, so decoys have profiles to win. (2) Benchmark: S1 panel
leave-one-out, plus held-out non-seed reviewed orthologues scored by gene
symbol (never classified by it). (3) Sweep every profile over the
non-redundant census v2. (4) Merge with S2 by a stated rule (agree / one
speaks / conflict) and calibrate against S2's calls with seeds excluded.

**Completion criteria.** Every family has a profile with a recorded seed
set and SHA-256s; benchmark reported per family with LOO; every v2 record
carries a profile verdict and a merged v3 call; conflicts counted, not
resolved by fiat; report rendered from tables.

---

## S3b — Proteome sweep + jackhmmer → census v3

**Goal.** Find what domain search misses: sweep S3a's profiles over S4's
declared reference proteomes and assign by *best profile* with a margin
(D7); jackhmmer runs to convergence with a coded kill criterion (D10),
including a sister-family contamination rule (the parent project's K1–K3).

**Notes.** This is the tier that reaches families whose Pfam coverage is
poor (H7) and lineages with no close reference. Needs S4's manifest.

---

## S4 — Proteome scope (the denominator)

**Goal.** A declared manifest of reference proteomes across the lineage
panel in `src/utils/species.py`, with assembly accession, source, date and
gene-set provenance per row (D9). Needs the external drive.

**Done 2026-09-28** (D35): `scripts/s4_proteomes.py all` → `results/proteome_scope/`.
Candidates are reference proteomes the pinned release ships; selection is a
rule; species without one are `genome_only`, never dropped; every file is
checked against the metalink MD5 and the README counts. Outputs S3b uses:
`<data root>/proteomes/s4/panel_refprot.fasta` and the manifest's `upid` per
species. Outputs S5 uses: the manifest's `assembly_id` (the assembly the gene
set was built on; `current_assembly` where it is superseded), level, N50 and
annotation source.

---

## S5 — Genomic sweep for what the proteomes miss

**Goal.** tblastn + miniprot for every family that S3 finds absent in a
lineage, so an absence claim passes both bars (D4).

**Split 2026-09-29** (as the parent project did): **S5a** builds and measures
the instrument on a pilot; **S5b** runs it over all 52 panel genomes.

### S5a — instrument + pilot

Scripts `s5_genome_io.py` (fetch, MD5, index), `s5_baits.py` (B1–B3),
`s5_lib.py` (miniprot, parse, overlap loci), `s5_classify.py` (locus
translations → S3a profiles → D32), `s5_rescue.py` (tblastn),
`s5_ledger.py` + `s5_verdict.py` (cells, controls, verdicts),
`s5_calibrate.py` (span / intron from annotation), `s5_report.py`;
driver `s5_sweep.py`. Pilot: *E. coli*, *S. cerevisiae*, *Arabidopsis*,
*C. elegans*, *Drosophila*, mouse, *Cornu* (genome-only).

**Completion criteria.** (1) Bait panel derived by rule, every family
baited. (2) Every pilot genome swept, loci called by the profiles, 0
failures. (3) Positive control measured per genome, with self-baits
excluded, *matched* (a same-group non-self bait exists) vs unmatched
reported separately. (4) `-G` and D4's span bar measured from annotation.
(5) The control floor set from the pilot, stated as a decision. (6) The
literature-expected absences (ZAC mouse; Nav, P2X *C. elegans*; P2X
*Drosophila*) read out. (7) Self-test covers the new invariants; report
rendered from tables.

### S5b — full sweep

`python3 scripts/s5_sweep.py run --all`, then `s5_calibrate.py`,
`s5_verdict.py`, `s5_report.py`. Budget from S5a's timings. Census v4 =
census v3 + genome-only loci (never merged into a proteome call).

---

## S6 — Alignments and pore modules

**Goal.** Per family: MAFFT L-INS-i + trimAl. Per alignable superfamily: pore
modules extracted by `src/phylo/modules.py`, with the extraction method
recorded per sequence — a tree built from modules found two different ways
has a confounder in it.

---

## S7 — Tier-1 phylogenies

**Goal.** One rooted ML tree per census family. IQ-TREE 2, ModelFinder, 1000
UFBoot, outgroup from the catalogue. Report the families that could not be
rooted because their outgroup was not in the sequence set.

---

## S7c — Re-rooting the weakly rooted P-loop families

**Why.** S7b rooted every P-loop family on the KcsA/MthK/NaK exemplars —
two-helix bacterial pores that fill only 15–42 % of a six-helix or 2×2P
alignment's columns. The root then rests on a few dozen pore positions at
the end of a long branch, the classic long-branch-attraction setting.
Measured: the outgroup is **not one clade** in CNG, K2P, Kv (KCNQ) and Kv
(Shaker), and the root is weak in Slo (UFBoot 32) and Kv (EAG) (70). S11's
duplication history reads these roots, so they are fixed first.

**The six families.** `cng`, `k2p`, `kv_kcnq`, `kv_shaker`, `kca_slo`,
`kv_eag`.

**Approach (the user agreed, 2026-10-01): closer outgroups + denser basal
sampling, with an outgroup-free check later.**

1. **Architecture-matched outgroups, declared in the catalogue first.**
   Prokaryotic relatives sharing the family's full architecture, verified
   live by S0 before use (accession, gene, Pfam architecture):
   * Kv families (Shaker, KCNQ, EAG) and Slo — **KvAP** (*Aeropyrum
     pernix* archaeal six-helix voltage-gated K⁺ channel) and any other
     verified prokaryotic 6TM VSD-K⁺ channel;
   * CNG — **MloK1** (*Mesorhizobium loti*) and **SthK** (*Spirochaeta
     thermophila*) cyclic-nucleotide-gated channels;
   * K2P — no prokaryotic K2P exists: use **two** sister eukaryotic K⁺
     families (candidates Kir, Kv) as a combined outgroup.
   The outgroups enter the catalogue as a per-family declaration (a new
   field, e.g. `ChannelFamily.root_with`, overriding the superfamily's) with
   provenance; `s7_trees.root_rule()` reads it. No outgroup is chosen by
   which root looks best.
2. **Denser basal sampling from S4b.** Add each family's high-confidence
   members from the early-diverging lineages of the dense panel
   (choanoflagellates, *Capsaspora*, sponges, placozoans, ctenophores if
   present, protists) — by rule, D8 (per clade × kingdom), stated before the
   run. This changes these six families' sequence sets away from S6's D39
   sets: record it as a stated exception, never silently.
3. **Re-align and re-tree** the six with the S7 method unchanged (L-INS-i →
   trimAl `-gt 0.5` → IQ-TREE 2 `-m MFP -mset LG,WAG,JTT,Q.pfam -B 1000
   -bnni -seed 1`). Run detached (the large three are 300–400 sequences,
   ~10–20 h overnight).
4. **A root is accepted only by a criterion fixed before any tree is read:**
   (a) the outgroup forms one clade; (b) root-edge UFBoot ≥ 95; (c) **the
   root is the same under two different outgroup choices** (e.g. the
   prokaryotic outgroup alone vs the sister-family outgroup). A root that
   moves with the outgroup is reported as **unresolved**, never picked.
5. **Optional cross-check now, required in S11:** a non-reversible-model
   root (IQ-TREE root inference with rootstrap) for each of the six;
   reconciliation rooting against the species tree is S11's job (minimise
   duplications + losses), and disagreements are reported, not resolved.

**Completion criteria.** Catalogue edit + S0 clean; six trees built; for
each, the three root criteria reported in a table
(`results/phylogeny/tier1_reroot.tsv`) with the old KcsA root alongside;
`tier1_report.md` and the S7 figure updated (panel C shows the new roots);
the families whose root is still unresolved named; self-test invariant for
the per-family outgroup rule.

**Outputs.** `results/phylogeny/tier1_reroot.tsv`, re-rooted treefiles under
`results/phylogeny/tier1/`, bulk under `<data root>/trees/s7c/`.

**Split (2026-10-01).** S7c delivered the design, catalogue edit, S0 check, basal sampling, alignments and the 12 inputs, and launched the trees. **S7d** finishes: wait for the detached run (`<data root>/trees/s7c/run.log`; re-run `scripts/s7c_reroot.py run --jobs 5 --threads 2` if it died — it resumes), then `s7c_reroot.py parse` → `tier1_reroot.tsv` / `tier1_reroot_trees.tsv`, `s7_report.py` (§ 4 fills itself), `bin/envpy scripts/s7_figures.py` (panel C shows the S7b root and both new roots with the D47 verdict), look at the figure (D11), and name the unresolved families in the ledger. The criterion and the outgroups are fixed (D47); nothing about them changes after a tree is read.

---

## S8 — Tier-2 pore-module trees and the tier-3 fold network

**Goal.** One pore-module tree per alignable superfamily, and the fold
network for everything else. `build_tier2()` will refuse the non-alignable
ones; **record the refusals in the report** — they are a result, not an
error. Then measure the `LITERATURE_EDGES` with Foldseek and flip each from
`measured=False` to a real value or leave it unmeasured with a reason.

---

## S9 — Selectivity-filter atlas

**Goal.** The filter locus for every P-loop family, and whether it is
congruent with the tier-2 tree (Q5). Cav3's `EEDD` is the test case: motif
and tree may disagree, and both get reported.

---

## S10 — Repertoire evolution

**Goal.** Presence/absence of every family across the lineage panel, with
ancestral reconstruction. The MscS animal-absence question (Q4) is the
sharpest claim and needs D4's two bars.

---

## S11 — Duplication history

**Goal.** Which families expanded in which lineage; 2R/3R ohnologue status;
the order of the internal repeat duplications that made the 4×6TM channels.

---

## S12 — Structures

**Goal.** AFDB coverage per family, TM-align against the experimental
reference set, Foldseek all-vs-all for the fold network. Needs the drive.

---

## S13 — Selection

**Goal.** ML selection tests on the families with real clinical variant sets.
Pairs with S19.

---

## S14a — Manuscript assembly

**Goal.** `scripts/s14_assemble.py` end to end: figures → claims → stitch →
deposit. Every load-bearing number needs a claim row (D12). Every figure gets
looked at (D11).

---

## S15–S22 — the analysis block

Each runs on data the earlier tasks produced; see the ledger for
dependencies and priority. The two with the most reach:

**S15 — Method contribution.** The recall curve as methods are added
(domain → profile → genomic). This is the answer to Q3 and the strongest
methods result the project can produce.

**S18 — Convergence audit.** Connexin vs pannexin, TMEM175 vs the GYG
channels, the mechanosensitive families. "Unrelated" is a statement about
detection methods, and this task is where it either survives profile-profile
comparison and structure or does not.
