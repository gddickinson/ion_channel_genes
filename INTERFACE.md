# INTERFACE.md — Module Navigation Map

> Read this **before** opening any source files. Each `.py` is < 500 lines
> and has one focused responsibility.

## Top-level

| File | Purpose |
|------|---------|
| `run.py` | CLI/GUI entry. `--catalogue`, `--classify`, `--phylo` need only stdlib + `requests`; `--headless` and the GUI need Biopython (D18). |
| `PUBLICATION_ROADMAP.md` | **The multi-session publication plan**: session protocol, task ledger (S0–S24; S3 split into S3a/S3b), emergent tasks, and the Decisions log — D1–D18 inherited from the PIEZO/IP3R projects, D23–D28 new to a multi-family subject. Read at the start of every session. |
| `docs/session_briefs.md` | Per-task instructions: goal, steps, completion criteria, outputs. |
| `docs/channel_review_2026.md` | **The literature review** — 16 sections, ~8,400 words, 148 references, covering the folds, selectivity, gating, each superfamily, evolution, what the annotation databases record (§12, this project's own measurements), disease and pharmacology, and the open questions. **Generated — never hand-edit it**; edit `docs/review/*.md` and re-run `scripts/s0_review_build.py`. |
| `docs/channel_review_2026.pdf` | The typeset review (A4, pandoc + xelatex), built by `s0_review_build.py --pdf`. |
| `docs/review/` | The review's **source**: numbered section files `00_frontmatter` … `15_methods`. Citations are stable keys (`[doyle1998]`, `[hilf2008, bocquet2009]`) resolved against `results/s0_baseline/references.tsv`; the build renumbers them into order of first appearance, and **a cited key with no verified reference is a build error**. |
| `docs/channel_background.md` | **The biology baseline.** Every statement tagged `[db]` / `[lit]` / `[open]`, including the measured annotation facts that decide how the census must run. |
| `docs/scope_and_boundaries.md` | **What counts as an ion channel** — the six boundary questions and this project's answers (D23). The reason the census has a denominator. |
| `docs/classification_rules.md` | How a protein gets classified: the three tiers, the rules they apply, the confidence scheme, and what the classifier deliberately does not do. |
| `docs/phylogeny_protocol.md` | **Why the phylogeny is a forest** (D27), the three tiers, the rooting table, and what replaces the tree nobody can build. |
| `docs/analysis_catalogue.md` | What the harvested data can answer — the analysis menu behind S15–S22, plus the analyses deliberately not attempted. |
| `CLAUDE.md` | Session-start instructions, conventions, the three things that will bite you. |
| `INTERFACE.md` | (this file) navigation map. |
| `README.md` | User-facing docs + the project status board. |
| `SESSION_LOG.md` | Running notes per session: what ran, what resulted, what's next. |
| `FINDINGS.md` | **The biological story, task by task** — plain-language, appended after every completed task. |
| `roadmap.md` | Feature-development history of the app itself (inherited v1.0–v1.7 plus this project's v2.x). |
| `data_root.txt` | Path to the bulk-data root, read by `src/utils/data_root.py`. |
| `requirements.txt` | `requests` (required), `biopython` + `matplotlib` (search / analysis / figures only). |
| `presets/` | `channelome_human` (all 320 census genes), `controls_benchmark` (the S1 panel: positives + one decoy per hazard), `ploop_domain_scan`, `channelome_all` (exhaustive mode), and the scoped surveys `kv_survey`, `cysloop_survey`, `intracellular_survey`, `channel_discovery`. |
| `cache/` | On-disk JSON cache (auto-created, gitignored). |
| `results/` | Committed analysis output, one directory per task. `s0_baseline/` holds the catalogue verification and `reference_panel.fasta`; `benchmark_controls/` holds the S1 benchmark; `census_v2/` holds S2's census tables (the census itself is on the data root); `proteome_scope/` holds S4's denominator — `proteome_manifest.tsv` (one row per panel species: proteome or genome-only, assembly, annotation source, BUSCO, CPD), `proteome_candidates.tsv`, `proteome_files.tsv` (MD5/README verification, census v2 overlap) and `summary.json` (the files and the sweep DB under `<data root>/proteomes/s4/`); `census_v3/` holds S3a's seed manifest, profile build record, benchmark, calibration and census v3a tables (profiles, alignments and hmmsearch output under `<data root>/hmmer/s3/`); `panel_sweep/` holds S3b's census v3 over the S4 panel — sweep and instrument-check records, `domain_search_missed.tsv` / `missed_by_group.tsv` / `missed_records.tsv` (what domain search missed, with a high-confidence column), `family_by_species.tsv` (the presence matrix S5/S10 start from), `human_recall.tsv`, the jackhmmer seeds / runs / rounds / completeness / candidates tables and `summary.json` (domtbls, jackhmmer logs and the bulk census v3 under `<data root>/hmmer/s3b/`); `genome_sweep/` holds S5's bait panel, `genome_runs.tsv`, `cells.tsv` (every species × census family with its genome verdict and D4 bar), `controls.tsv`, `traces.tsv`, the annotation calibration (`spans.tsv`, `introns.tsv`, `annotated_loci.tsv`), census v4's tables (`census_v4_status.tsv`, `family_by_species_v4.tsv`, `genome_found_check.tsv`, `census_v4.json`) and `report.md` (genomes, GFFs, loci, tblastn output and census v4 under `<data root>/genomes/`); `phylogeny/` holds one directory per tree with its alignment, trimmed alignment, treefile and `run.json` command record. |
| `manuscript/` | The submission package, built by `scripts/s14_assemble.py`. |

---

## `src/catalogue/` — the project's definition of its own subject

**The single place any channel family, signature, size band or exemplar is
declared.** Nothing else in `src/` hard-codes a gene name. 91 families in 25
superfamilies, of which 68 are census families covering 320 human genes; the
other 23 exist so they can be excluded.

| File | Contents |
|------|----------|
| `schema.py` | `ChannelFamily`, `Superfamily`, `Signature` (with `enumerate` — whether the census searches on it, declared separately from `level`, D34), `Exemplar`, `Hazard`, and the enums: `Level` (the level at which a signature is diagnostic — the idea D25 rests on), `Provenance`, `Gating`, `Selectivity`, `Fold`, `Status`. |
| `vgic_k.py` | P-loop superfamily, potassium branch: Kv1–4, the silent modifiers, KCNQ, EAG/ERG/ELK, Slo, SK/IK, Kir, K2P, and the prokaryotic channels that root them. `PLOOP_SUPERFAMILY` and the shared `ION_TRANS` / `ION_TRANS_2` signatures live here. |
| `vgic_cation.py` | P-loop, non-potassium: Nav, Cav, NALCN, CatSper, TPC, CNG, HCN, and TRPC/V/M/A/ML/P/N. Carries the measured finding that most TRP families have no `PF00520`. |
| `lgic.py` | Four unrelated ligand-gated superfamilies: Cys-loop (nAChR, GABA-A, GlyR, 5-HT3, ZAC, invertebrate GluCl, GLIC/ELIC), iGluR (AMPA, kainate, NMDA, delta, plant/insect, GluR0), P2X, DEG/ENaC. |
| `anion.py` | CLC channels *and* CLC transporters (D24), bestrophin, CFTR, tweety. |
| `tmem16_like.py` | Anoctamin channels, anoctamin scramblases, OSCA/TMEM63, TMC — one fold, `alignable=False`. |
| `mechano.py` | Piezo, MscL, MscS. |
| `largepore.py` | Connexins; the innexin/pannexin/LRRC8 clan; CALHM. |
| `intracellular.py` | ITPR, RYR (the parent project's whole subject, here two families of ninety-one), TRIC, MCU, VDAC, TMEM175. |
| `other.py` | ORAI, Hv1, otopetrins, CLIC, viroporins. |
| `controls.py` | The 23 non-census families: auxiliary subunits, domain-sharing non-channels (KCTD, class C GPCRs, AChBP, the generic ABC-transporter decoy (S3a2), POMT, VSP), transporters, and the out-of-scope channels (aquaporins, bacterial porins, gasdermins). |
| `hazards.py` | **The hazard registry** — 16 recorded ways to classify wrongly, each with the families involved, the shared evidence, the discriminating positive test, and the module that owns it. |
| `registry.py` | Assembly and lookup: `CATALOGUE`, `SUPERFAMILIES`, `census_families()`, `control_families()`, `human_genes()`, `signature_index()`, `shared_signatures()`, `pore_signatures()` (the H7 fix), `exemplars()`, `reference_panel()`, `hazards_for()`, `validate()`, `stats()`. |
| `__main__.py` | `python3 -m src.catalogue [--families] [--shared] [--census]` — validate and print. Exits non-zero on any validation problem, so it works as a pre-commit gate. |

## `src/classify/` — three tiers, one call, an audit trail

| File | Key types / functions |
|------|-----------------------|
| `rules.py` | `ArchitectureRule` (domain `require` / `min_copies` / `forbid`, plus whole-sequence bounds `max_tm`, `min_length`, `max_length`, `complete_only` — a rule that states one never fires on an unknown measurement), `derived_rules()` (generated from the catalogue), `HAZARD_RULES` (hand-written, priority 100, one per closable hazard), `match_architecture()` → `RuleMatch` (ties are *ambiguity*, never a coin toss), `topology_check()`. |
| `motifs.py` | `K_FILTER_RE` (TxGYG — the one test needing no reference), `FOUR_REPEAT_ANCHOR` (human Nav1.5 `Q14524` positions 372/898/1419/1711, verified 2026-08-19), `filter_signature()`, `verify_anchor()`, `scan()` → `MotifResult`. `FILTER_CALLS` maps DEKA/EEEE/EEDD/EEKE onto families. |
| `reference.py` | `ReferenceSet` (loads `reference_panel.fasta`, 3-mer prefilter, MAFFT identity over mutually covered columns), `best_match(...)` → `ReferenceResult` with the D7 margin, the D29 leave-one-out and the D30 coverage floor (`MIN_COVERAGE = 0.30` of the longer sequence; hits below it are dropped and counted). |
| `classifier.py` | `ChannelQuery` (Pfam counts, TM count, length, UniProt fragment flag), `ChannelCall`, `Evidence`, `classify()`, `classify_all()`. Tier precedence, confidence tiers, conflict recording. The gene symbol is carried and never consulted (H15). |
| `report.py` | `calls_table`, `confusion_table`, `recall_table`, `hazard_table`, `coverage_table`, `summary_counts`, `write_tsv` — the tables S1's report is rendered from (D13). |

## `src/phylo/` — a forest and a network

| File | Key types / functions |
|------|-----------------------|
| `modules.py` | `PoreModule`, `ModuleSet`, `module_from_envelope()`, `modules_from_envelopes()`, `module_by_projection()`, `whole_sequence_module()`. Extraction method is recorded per sequence. |
| `forest.py` | `build_tree()` (MAFFT → trimAl → IQ-TREE 2, every command recorded in `TreeRun`), `build_tier1()` (within family, full length), `build_tier2()` (within superfamily, pore module) which **raises `NotAlignable`** for a superfamily the catalogue marks non-alignable, `refused_superfamilies()`. There is no `build_tier3()` (D27). |
| `network.py` | `FoldNetwork`, `FoldEdge`, `seed_network()`, `LITERATURE_EDGES` — tier 3. Edges are structural similarity with no branch lengths and no support values; an unmeasured pair is recorded as unmeasured, never as unrelated. |

## `src/utils/`

| File | Functions |
|------|-----------|
| `scope.py` | **`Scope` and `scope_for(key)`** — narrows the catalogue to "all", a superfamily or a family, and derives the signatures, size band, reference panel, *sister panel* (from the hazard registry) and known genes for a run. The module-level constants at the bottom are the compatibility layer that keeps the ported machinery working. |
| `mafft.py` | Stdlib-only MAFFT wrapper: `align()`, `parse_fasta()`, `project_positions()`, `alignment_stats()` (identity **plus coverage of both sequences** — the D30 guard), `percent_identity()`, `MafftUnavailable`. Raises rather than degrading (D28). |
| `species.py` | 52 species across 15 groups, from *E. coli* to human plus three viruses; `PANELS` (human / vertebrate / metazoan / eukaryote / prokaryote / outgroup / virus / tree_of_life), `resolve_species()`, `ensembl_species_slug()`, `panel()`. |
| `data_root.py` | `get_data_root()`, **`require_data_root()`** (raises — use before anything bulk), `free_bytes()`. `python3 -m src.utils.data_root --require` is the session-protocol check. |
| `exporters.py` | `write_fasta()`, `write_csv()`, `write_json()`. |
| `results_writer.py` | `make_bundle_dir()`, `write_bundle()`. |
| `report.py` | `write_report()` — publication-style markdown + standalone HTML. |

## `src/core/` — data models, search orchestration, caching

| File | Key types / functions |
|------|-----------------------|
| `models.py` | `ProteinVariant`, `GeneRecord`, `SearchQuery`, `SearchResult`, `SearchStatus`. |
| `cache.py` | `DiskCache(root, ttl_s)` — JSON-file cache keyed by `(source, query)`. |
| `search.py` | `SearchOrchestrator` — fans queries across enabled DBs on worker threads, marshals results back via `queue.Queue`. Re-exported lazily from `src/core/__init__.py` (PEP 562) so importing `src.utils` does not drag in Biopython. |

## `src/databases/` — one client per source

| File | Class | Notes |
|------|-------|-------|
| `base.py` | `DatabaseClient` | Abstract; subclasses implement `search()`. |
| `ncbi.py` | `NCBIClient` | Biopython `Bio.Entrez` against the `protein` index. |
| `ensembl.py` | `EnsemblClient` | `rest.ensembl.org`. Uses `lookup/symbol` rather than `xrefs/symbol`, which the parent project measured to stall indefinitely for `homo_sapiens`. |
| `uniprot.py` | `UniProtClient` | `rest.uniprot.org/uniprotkb/search`. |
| `alphafold.py` | `AlphaFoldClient` | `alphafold.ebi.ac.uk` — pivots off UniProt; reports pLDDT. |
| `foldseek.py` | `FoldseekClient` | Structure-based remote homology; opt-in, async. The tier-3 network's measuring instrument. |
| `compara.py` | `ComparaClient` | Ensembl Compara gene trees — paralogues regardless of naming. |
| `blast.py` | `BlastClient` | Sequence-bait BLAST via NCBI's public queue. |
| `interpro.py` | (functions) | `list_proteins_with_pfam()` (census mode: paginates to the API's own count, `strict` failure, raw-page `dump_dir`), `batch_fetch_family_signatures()`, `fetch_uniprot_sequence()`. |

## `src/analysis/`, `src/discovery/`, `src/investigation/`, `src/gui/`

Ported from the parent project and retargeted at the catalogue; see
`roadmap.md` for what changed. The pieces that matter here:

| File | Key types / functions |
|------|-----------------------|
| `analysis/pipeline.py` | `analyse(variants, …)` → `AnalysisResult`; `write_analysis()`. |
| `analysis/distance.py` | `identity_matrix(covered_only=True)` — fragment-aware identity. |
| `analysis/motifs.py` | MSA-derived family-signature PSSMs — the domain-evidence fallback for candidates with no InterPro record. |
| `discovery/candidates.py` | `DiscoveryConfig` / `discover_novel_paralogs()` — the composite scorer with the D3 evidence gate and the **sister-family test** (D14), whose sister panel now comes from the hazard registry via `Scope`. |
| `discovery/domain_scan.py` | `run_domain_scan()` — enumerate by signature, then classify. |
| `discovery/exhaustive.py` | `run_exhaustive_hunt()` — the `"mode": "exhaustive"` pipeline. |
| `investigation/pipeline.py` | `Investigator(accession, options).run()` — seven evidence lines on one candidate; the default comparison panel is the scope's exemplars **plus its sister families**, so the hazard margin has something to measure against. |
| `gui/app.py` | `MainWindow`; `run(project_root, email)`. |
| `gui/tools.py` | The Analysis menu, the task runner, and the search → analyse → discover path the headless CLI also runs. |
| `gui/tools_extra.py` | `ToolsExtraActions`, mixed into `ToolsController`: the one-off investigations — selection test, fold check, presence matrix, domain scan, exhaustive hunt, deep dive, and **`classify_selection()`**, which runs the three-tier classifier on the selected row and shows the whole audit trail rather than just the answer. |

## `src/cli.py` and `src/cli_channel.py`

`cli.py` — `run_headless()` (preset search → optional analyse → optional
discovery → bundle) and `run_investigate()`, both scope-aware. Needs
Biopython.

`cli_channel.py` — `run_classify()`, `run_catalogue()`, `run_phylo()`, plus
the stdlib fetchers (`pfam_counts`, `resolve_symbol`, `entry_meta`,
`build_query`). **Deliberately Biopython-free**: classification is the core
operation and cannot be gated on an optional dependency.

---

## `scripts/` — roadmap-session tooling (not part of the app)

| File | Purpose |
|------|---------|
| `dashboard.py` | **Project dashboard** (stdlib-only) → `dashboard.html`: progress parsed from the roadmap ledger, a live panel driven by `results/session_live.json`, key figures base64-embedded, and `FINDINGS.md` rendered. Session protocol step 0. |
| `figstyle.py` | **The one figure style.** `SUPERFAMILY` (four safe categorical hues for twenty-five superfamilies — the constraint is stated, not hidden), `SELECTIVITY`, `CONFIDENCE` (gold→unassigned), `STATUS`/`QUALITY` (ordered evidence scales), `CLINICAL`. `save()` raises rather than writing a figure whose bbox runs off-canvas. |
| `selftest.py` | **The offline invariants**, under a second and no network: catalogue validation, every hazard rule fired against the synthetic architecture it was written for, the classifier's blindness to gene symbols, the motif regex, the alignment helpers, the D27 refusal, and the scope's derivation of its sister panel from the hazard registry. Session-protocol step 4. |
| `s0_lib.py` | S0 fetchers: `Fetcher` (polite, retrying, keeps a failure ledger), `pfam_entry`, `protein_pfams`, `pfam_protein_count`, `resolve_gene`, `entry_by_accession`, `sequence`, `taxon`, `write_tsv`, `read_tsv`, `live_progress`. |
| `s0_catalogue_verify.py` | **S0** — re-derives every claim the catalogue makes against live databases and writes six TSVs plus `reference_panel.fasta`. Reports mismatches; never edits the catalogue. |
| `s0_figures.py` | **S0's figure** — families and human genes per superfamily (the lopsidedness: 143 of 320 genes in one superfamily), and every domain signature carried by more than one family with the ones that cross the channel / non-channel boundary marked. Drawn from the committed S0 tables through `figstyle.py`. |
| `review_sources.py` | The review's source list: 148 `(key, title)` entries, optionally with a Europe PMC constraint for papers whose titles are too generic to rank. **Titles, not citations** — nothing here is a reference until it resolves. |
| `s0_review_refs.py` | Resolves every source title against Europe PMC and writes `results/s0_baseline/references.tsv` with the returned PMID, DOI, year, journal and authors. A record is accepted only if its title is ≥ 90 % similar to the one requested, so a misremembered paper fails loudly instead of becoming a plausible bibliography entry. Exits non-zero on any failure. |
| `s0_review_build.py` | Assembles `docs/channel_review_2026.md` from `docs/review/*.md`, renumbers citations, renders the bibliography; `--check` validates without writing, `--pdf` typesets via pandoc + xelatex. |
| `s0_domain_map.py` | Fetches InterPro match **coordinates** for the 22 proteins the review's architecture figures compare → `results/s0_baseline/domain_positions.tsv`. `exemplar_architecture.tsv` records *which* domains a protein carries; this records *where*. |
| `s0_filter_atlas.py` | Extracts the two selectivity-filter alignments from real sequences → `filter_k.tsv` (TxGYG located directly, no alignment used) and `filter_four_repeat.tsv` (the four-repeat locus projected from Nav1.5 by MAFFT). Records near misses explicitly, so *Bacillus* NaK's TVGDG appears as what it is rather than as a blank. |
| `review_figlib.py` | Drawing primitives shared by the review figures: coloured sequence rows, membrane cartoons, β-strands, stoichiometry top-views, scale domain bars, and the one superfamily short-name map every figure uses. |
| `s0_review_fig_filters.py` | Review figures 2 and 3 — the two filter alignments. |
| `s0_review_fig_domains.py` | Review figures 4, 5 and 8 — length range, the shared-signature matrix, and the architecture traps to scale. |
| `s0_review_fig_folds.py` | Review figures 1, 6 and 7 — the fold gallery (schematic), the real Cys-loop ML tree drawn from its Newick with bootstrap support, and the forest-not-a-tree diagram. |
| `s0_report.py` | **S0** — renders `results/s0_baseline/report.md` purely from those tables (D13). |
| `s1_toolchain.py` | **S1 step 1** — probes every external binary and Python package, writes `results/toolchain_manifest.txt`. |
| `s1_benchmark.py` | **S1 driver** — verifies the filter anchor, loads the reference panel, classifies the control panel, writes recall / confusion / hazard / coverage tables. |
| `s1_report.py` | **S1** — renders `results/benchmark_controls/report.md` from those tables, including the per-tier attribution and the untested-hazard list. |
| `s2_lib.py` | **S2** shared pieces: the taxonomic `SHARDS` that partition the census query, `union_clause()` (from `pore_signatures()`), strict backoff `get()`, `uniprot_count()` / `interpro_count()`, `parse_record()` (UniProt JSON → flat row with every Pfam and its `MatchStatus` copy count, TM-feature count, sequence), `pfam_dict()`, `iter_pages()`. Bulk output to `<data root>/raw_api/s2/` (D31). |
| `s2_enumerate.py` | **S2 step 1** — `counts` (partition check + per-signature UniProt and InterPro counts), `walk` (resumable parallel cursor walk, raw pages archived gzipped), `verify` (fetched == count per shard, union and signature). |
| `s2_classify.py` | **S2 step 2** — `--recheck PF…` re-classifies only records carrying the given accessions and copies every other call, archiving the previous files under `calls_<archive>/` (refuses to overwrite an archive; used by S2b). Otherwise: `classify()` on every record from its UniProt evidence; architecture + K⁺ regex everywhere, four-repeat projection only where `PF00520` ≥ 4, reference tier `not_run`. One gzipped call file per shard; resumes by shard. |
| `s2_census.py` | **S2 step 3** — assembles `census_v2.tsv.gz` + `census_v2.fasta.gz` on the data root (SHA-256 in `manifest.tsv`) and writes the committed tables: status, families, superfamily-only, tier attribution, four-repeat filters, unassigned signatures, partial architectures, human recall, S1-panel comparison. |
| `s2b_revision.py` | **S2b / S2c** — every call transition between two archived call revisions (`--tag`, `--base`, `--head`, `--accessions`) → `<tag>_transitions.tsv`, with a hard check that nothing changed outside the re-checked records. S2b = r1 → r2 (H2/H4/H13), S2c = r2 → r3 (H11/H12). |
| `s2_report.py` | **S2** — renders `results/census_v2/report.md` from those tables (D13). |
| `s3_hmm_lib.py` | **S3** shared pieces: `s3_dir()` (`<data root>/hmmer/s3/`, raises without the drive), streaming FASTA IO, `seq_id()` (the non-redundant key), `iter_domtblout()` / `hits_by_target()` (full-sequence score + merged profile coverage from domains with i-E ≤ 1e-3), `iter_census_v2()`, TSV IO. Ported from `../ip3r_genes`. |
| `s3_seed_spec.py` | **The seed rules**, enforced: R1 curated human genes (resolved by *primary* gene name — `resolve_primary()`, because `gene_exact` matches synonyms and resolved TRPC7 to TRPM2), R2 catalogue exemplars, R3 reviewed S2 family calls in band, one per species, round-robin over groups, ≤ 30. Seed sets disjoint across families or the build aborts. |
| `s3_build_profiles.py` | **S3a step 1** — one HMM per catalogue family (90): MAFFT L-INS-i single-threaded → hmmbuild; hard failure on MAFFT exit, row count or ragged output; SHA-256 of seeds / alignment / profile → `profile_build.tsv`; `profiles/all.hmm`. |
| `s3_assign.py` | **Best-profile assignment** (D32): score ≥ 30 bits, ≥ 30 % of the profile's match states covered (a shared module is not a call), ≥ 10 % relative margin over the best other profile (D7); inside the margin → `superfamily_only` if the band is one superfamily, else `ambiguous`. `collect_hits()` streams domtblouts keeping a per-target top list. |
| `s3_benchmark.py` | **S3a step 2** — the instrument measured before use: S1 panel with every seed query scored against its family profile **rebuilt without it** (LOO), plus held-out non-seed reviewed orthologues scored (never classified) by gene symbol. |
| `s3_sweep.py` | **S3a step 3** — `prep` collapses census v2 to unique sequences (+ accession map); `search` runs hmmsearch per profile with `-Z` fixed, resumable on the profile's SHA-256. |
| `s3_census_v3.py` | **S3a step 4** — `merge()`: both agree → `both`; one speaks → that one; different families, or a profile family outside S2's superfamily → `conflict` (kept, counted). Writes census v3a (bulk) and the calibration / conflict / resolution / human-recall tables. |
| `s3_external_check.py` | **S3a external check** — census v3a against the parent projects' censuses (`../ip3r_genes` census v6: call against call; `../piezo_genes` census v5: membership only, it has no per-record call). Read-only comparison on shared UniProt accessions, both directions, with each disagreement listed; parent files recorded by SHA-256. Never an input to a call. |
| `s3_report.py` | **S3a** — renders `results/census_v3/report.md` from those tables (D13). |
| `s4_proteome_lib.py` | **S4** shared pieces: the three scope rules as code — `select_proteome()` (most Swiss-Prot → BUSCO → genes → UPID), `assembly_rank()` (ported from `../ip3r_genes`), `file_verdict()` (metalink MD5 never waived; README counts waived only for a dated mid-release reissue whose count is the declared gene count, D35) — plus archived fetchers (UniProt proteomes API, release README, metalink, NCBI Datasets with `all_assemblies`, since proteomes often sit on a superseded assembly version) under `<data root>/raw_api/s4/`, `parse_readme()`, `flatten_proteome()`, `gene_groups()`. |
| `s4_proteomes.py` | **S4 driver** — `manifest` (candidates → selection → assembly per row; genome-only rows for species with no reference proteome), `download` (canonical FASTA + gene2acc → `<data root>/proteomes/s4/`), `verify` (MD5, README counts, human positive control, census v2 overlap, concatenated sweep DB `panel_refprot.fasta` for S3b). Reruns offline from the archive; `--refresh` re-queries. |
| `s4_report.py` | **S4** — renders `results/proteome_scope/report.md` from those tables (D13). |
| `s3b_lib.py` | **S3b** shared pieces: `s3b_dir()` (`<data root>/hmmer/s3b/`, raises without the drive), `panel_db()` (S4's sweep DB), `accession()`, `species_by_taxid()` (FASTA `OX=` → manifest row), and `parse_jackhmmer_log()` (per-round included target lists, ported from `../ip3r_genes`). |
| `s3b_sweep.py` | **S3b step 1** — `db` (panel universe table + DB SHA-256), `search` (S3a's 91 frozen profiles over the panel via `s3_sweep._search`, `-Z` = panel size, resumable on profile SHA-256), `assign` (S3a's D32 rule, imported). |
| `s3b_kill.py` | **D10's kill criterion for S3b**, ported from the IP3R project with its constants: K1 other-family share rises > 10 points over round 1, K2 > 10× growth, K3 10 rounds without converging. Uncalled targets never trip K1 (and so it cannot see repeat-domain drift — measured). |
| `s3b_jackhmmer.py` | **S3b step 2** — `seeds` (one per census family, derived: top high-confidence panel profile call), `run` (`jackhmmer -N 10 -E 1e-5`, resumable on seed + DB SHA-256), `parse` (offline from archived logs: D10 verdicts, per-round table with the reported-only other-superfamily and uncalled shares, completeness at family and superfamily level; K3 runs contribute no candidates, D36). |
| `s3b_census_v3.py` | **S3b step 3** — census v3 over the panel: in census v2 → its v3a call; outside it → the S3b profile call (`panel_profile`); jackhmmer-only → `candidate`, never a call. Also the S3a/S3b instrument check and the human-proteome positive control. |
| `s3b_report.py` | **S3b** — renders `results/panel_sweep/report.md` from those tables (D13). |
| `s5_genome_io.py` | **S5** genome plumbing: `fetch_one()` (NCBI Datasets, every file MD5-verified before install, `.done` resume, ported from `../ip3r_genes/scripts/fetch_genomes.py`), `verify_one()`, faidx-style `build_fai()` / `fetch_region()`, `longest_n_run()`, `assembly_stats()` (N50 of the file actually searched). Genomes under `<data root>/genomes/<assembly>/`. |
| `s5_baits.py` | **S5 bait panel**, rules B1–B3: high-confidence, complete, in-band S3b profile calls; **one per species** (D37) spread round-robin over groups; exemplar fallback. Every catalogue family, controls included. → `results/genome_sweep/baits.{faa,tsv}`. |
| `s5_lib.py` | **S5** alignment layer: `run_miniprot()` (atomic GFF), `run_miniprot_any()` (contig chunks above 2.5 Gbp), `parse_miniprot_gff()` (drops the genome's own species' baits — genome-scale LOO), `Aln`, `Locus`, `cluster_loci()` (**overlap**, not a 10 kb gap — tandem arrays hold different families), `max_intron_for()`. |
| `s5_classify.py` | **Locus → call**: each locus's top translations per bait family scored against all 91 S3a profiles (`-Z` = panel size) and assigned by D32's `assign_one`; the bait family is a second witness. → `<data root>/genomes/s5/<assembly>/loci.tsv.gz`. |
| `s5_rescue.py` | **tblastn rescue** for no-locus cells (control cells and informative zero cells): `run_tblastn()`, `outside_loci()` (an HSP inside any recorded locus is never a trace), `traces()`. |
| `s5_annotation.py` | The assembly's own (or RefSeq twin's) GFF3: `ensure_annotation()`, `seqid_map()` (any GenBank/RefSeq name → the name the searched FASTA carries), `read_genes()` (span + widest intron over transcripts; CDS blocks where exons carry no structure), `GeneIndex`. |
| `s5_calibrate.py` | **S5 calibration from annotation, never from miniprot**: gene span per family (D4's bar = median) and widest intron at called loci vs `-G` → `spans.tsv`, `introns.tsv`, `annotated_loci.tsv`. |
| `s5_ledger.py` | Cells (species × census family): proteome evidence, `informative`, `matched` (a same-group non-self bait exists), `control`, `genome_status()` (found / partial / gap / no_locus — the profile call decides). |
| `s5_verdict.py` | Rescue → per-genome **matched detection** control → D4 (`d4_bar()`, D38: group median → pooled → 3 × band in prokaryotes/viruses → none) → verdict per zero cell (`genome_found`, `partial`, `gap`, `trace`, `absent`, `absent_below_bar`, `uncontrolled`, `unmatched`, `no_locus_unrescued`). → `cells.tsv`, `controls.tsv`, `traces.tsv`. |
| `s5_sweep.py` | **S5 driver**: `baits`, `run --pilot / --species / --all` (fetch → miniprot → loci → profile calls, resumable on bait SHA-256) → `genome_runs.tsv`. |
| `s5_census_v4.py` | **S5b** — census v4: census v3's proteome rows unchanged + one row per genome locus the proteomes cannot supply (genome-only species: every profile-called locus; `genome_found` cells: high-confidence intact loci only) → `<data root>/genomes/s5/census_v4.tsv.gz`, `census_v4_status.tsv`, `family_by_species_v4.tsv`, `genome_found_check.tsv`, `census_v4.json`. Never merges a genome locus into a proteome call. |
| `s5_report.py` | **S5** — renders `results/genome_sweep/report.md` + `summary.json` from those tables (D13); panel-scale summaries, every cell in `cells.tsv`. |
| `s14_lib.py`, `s14_figures.py`, `s14_claims.py`, `s14_deposit.py`, `s14_pdf.py`, `s14_assemble.py` | Manuscript assembly: page geometry and figure maps; figure copying under publication numbers; **the claim checker** (every load-bearing number declared with the table and op that recovers it, D12); the deposit manifest with SHA-256s; the typeset PDF; the driver with ordered stages and non-zero exit on a failed claim. |
| `build_findings_page.py` + `findings_page.css` | Renders `docs/findings_summary.md` to one self-contained HTML page with figures inlined as WebP data URIs. |

Each ledger task adds its own `s<n>_*.py`. The parent projects'
equivalents (`../ip3r_genes/scripts/`, `../piezo_genes/scripts/`) are
reference implementations worth reading first.

---

## Data flow

### Classification (the core operation)

```
accession / gene symbol
        │  cli_channel.build_query()  → UniProt meta + InterPro pfam counts + sequence
        ▼
   ChannelQuery
        │
        ├─► rules.match_architecture()      hazard rules (100) > family rules (40-50)
        │                                   > superfamily-only (10); ties = ambiguous
        ├─► motifs.scan()                   TxGYG regex; four-repeat projection via MAFFT
        └─► reference.best_match()          3-mer prefilter → MAFFT → covered-only identity
                                            → D7 margin
        ▼
   combine by tier precedence (hazard > motif > reference > architecture)
        ▼
   ChannelCall  ── family · superfamily · status · confidence · margin
                   evidence[] from every tier · hazards[] · conflicts[]
```

### Phylogeny

```
catalogue exemplars / census sequences
        │
        ├─ tier 1 ─► build_tier1()  full length, rooted on the catalogue's outgroup
        ├─ tier 2 ─► modules.*  ─►  build_tier2()  pore module only; REFUSES if
        │                                          Superfamily.alignable is False
        └─ tier 3 ─► network.seed_network() + Foldseek  ── a network, never a tree
```

### Search (the ported path)

```
SearchPanel.on_search ─► MainWindow._on_search ─► SearchOrchestrator.start
                                                      │
                                    ┌─────────────────┼─────────────────┐
                                    ▼                 ▼                 ▼
                           NCBI / Ensembl / UniProt / Compara / AlphaFold
                                    └────────┬────────┴─────────────────┘
                                             ▼   (sequence DBs done)
                              gather longest sequences ─► BLAST + Foldseek
                                             ▼
                                        queue.Queue ─► _drain_queue (after())
                                             ▼
                                   ResultsView ─► DetailsView
```
