# CLAUDE.md — Project-specific instructions

@INTERFACE.md

## ⚠ Publication project — read this first every session

This project is executing a multi-session publication plan. At the start of
EVERY session:

1. Read `PUBLICATION_ROADMAP.md` and follow its **Session protocol**
   exactly: open the dashboard (`python3 scripts/dashboard.py --open` plus a
   background `--watch`), `git pull`, verify the data root
   (`python3 -m src.utils.data_root --require` — it exits non-zero if the
   external drive is not attached, and **that is a stop, not a warning**),
   pick the single next task from the ledger (first `in_progress`, else
   first unblocked `pending`), announce it, and work only that task to its
   completion criteria.
2. Detailed task instructions live in `docs/session_briefs.md`.
3. At session end, without being asked: update the roadmap ledger +
   Results (+ Emergent tasks for anything new), append to `SESSION_LOG.md`,
   **append a biological-findings entry to `FINDINGS.md`** (plain-language
   summary of what the finished task changed biologically; mark unconfirmed
   claims *(pending: task)*), **refresh `README.md` and the dashboard with
   the task's results and its headline figure** (drawn by
   `scripts/s<n>_figures.py` into `results/<dir>/figures/`, added to the
   README's *Results in figures* section and to `FIGURES` in
   `scripts/dashboard.py` — details in the roadmap's end-of-session step
   3b), then `git add -A`, commit
   (`S<n>: <one-line outcome>`), and push if a remote exists.
4. Bulk data (genomes, proteomes, HMM databases, structures) go under
   `get_data_root()` (`src/utils/data_root.py`) — never into the repo.

## Project summary

**Subject: every ion channel.** The project identifies, classifies and
reconstructs the phylogeny of the ion channels as a whole — not one family,
but the ~25 independent superfamilies and ~330 human pore-forming genes that
the term covers, plus their relatives across the tree of life.

Three deliverables, in order: a **census** (what exists, in a declared search
space), a **classification** (what each thing is, by positive tests with an
audit trail), and a **phylogeny** (a forest of trees plus a fold network,
because the superfamilies are not homologous to one another).

The subject is defined in `src/catalogue/` — 102 families, 32 superfamilies,
19 hazards — and nothing else in `src/` hard-codes a channel name.

## The three things that will bite you

**1. Domain architecture is not a family call.** `PF00520` (Ion_trans) is
carried by 20 catalogued families including a phosphatase that is not a
channel. Nav, Cav, NALCN and CatSper have *identical* architectures. Assign
by positive test — filter motif, reference margin, hazard rule — and record
which tier decided (**D25**, **D26**).

**2. Gene symbols group pores with their accessory subunits.** `KCNE*`,
`CACNB*`, `CACNG*`, `SCN*B`, `CATSPERB`, `KCTD*` all sort inside the channel
symbol space and none is a pore. The classifier never consults a symbol
(**H15**). Neither should any analysis.

**3. There is no tree of ion channels.** A nicotinic receptor and a Kv
channel share no alignable position. The phylogeny is tiered, and
`src/phylo/forest.py:build_tier2` *refuses* to build a tree across a
superfamily the catalogue marks `alignable=False` (**D27**). Do not work
around the refusal.

## Provenance

The app, the figure style, the dashboard, the manuscript-assembly and
claim-checking tooling and the session protocol are ported from the IP3
receptor project (`../ip3r_genes`), itself ported from the PIEZO project
(`../piezo_genes`). Their `INTERFACE.md` files map reference implementations
for most tasks here. **Port the method; never port a result.** The one
place the lineage shows in the science is `catalogue/intracellular.py`,
where ITPR and RYR are two of this catalogue's 102 families and the
parent project's census of them is an external check on ours.

## Conventions

- Every `.py` < 500 lines (split if it would exceed). Read `INTERFACE.md`
  before opening source files.
- The subject is defined in **one** place, `src/catalogue/`. A new family is
  an edit to one division file plus, at most, a hazard row. Never hard-code a
  gene name, Pfam id or size band anywhere else; read it from the catalogue
  through `src/utils/scope.py`.
- Every assertion in the catalogue carries a `Provenance` — `DB` (re-derived
  from a live database by an S0 script), `LIT`, `CURATED` (written by hand,
  not yet checked) or `OPEN`. Curated-but-unchecked is a legitimate state;
  silently wrong is not.
- A hazard is only closed by a **test**, never by a caveat, and the test must
  be a **positive** one: never "it is not a Nav because the symbol says
  CACNA1C", always "the filter reads EEEE".
- Each database client subclasses `src/databases/base.py:DatabaseClient` and
  normalises its records into `ProteinVariant`; new clients register in
  `src/databases/__init__.py:AVAILABLE_CLIENTS`.
- The GUI never makes network calls on the main thread — everything goes
  through `SearchOrchestrator` → worker threads → `queue.Queue` →
  `MainWindow._drain_queue` (polled with `root.after()`).
- Analysis scripts live in `scripts/s<n>_*.py`, one task per prefix, and
  render their `report.md` purely from their committed tables so report and
  data cannot drift (**D13**).
- A missing tool produces a missing result and a note, never a silently
  weaker method (**D28**). `MafftUnavailable` is raised, not caught.

## Running

```
python3 run.py --catalogue --scope ploop         # what the catalogue claims
python3 run.py --classify KCNQ1,KCTD1,TPTE       # classify, with audit trail
python3 run.py --classify-preset controls_benchmark --save-results
python3 run.py --phylo nav                       # tier 1 tree
python3 run.py --phylo tmem16_like --tier 2      # refused, with the reason
python3 scripts/selftest.py                      # offline invariants, <1 s
python3 scripts/s0_catalogue_verify.py           # re-derive every claim
python3 scripts/dashboard.py --open
```

The catalogue, classification and phylogeny commands need only the standard
library plus `requests`. The search / analysis / GUI paths need Biopython and
matplotlib, which live in the `piezo1` conda env — run those with
`/opt/anaconda3/envs/piezo1/bin/python` (**D18**).

## Adding a family to the catalogue

1. Add a `ChannelFamily` to the right `src/catalogue/<division>.py`, with
   `Provenance.CURATED` on anything you have not checked.
2. If it shares a signature with something it could be confused with, add a
   `Hazard` naming the discriminating test.
3. `python3 -m src.catalogue` and `python3 scripts/selftest.py` — both clean.
4. `python3 scripts/s0_catalogue_verify.py` — then read
   `results/s0_baseline/report.md` and fix what it flags.
5. Update `INTERFACE.md`.
