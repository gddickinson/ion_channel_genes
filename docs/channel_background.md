# The biology baseline

Every statement here is tagged by how far it can be trusted:

- **`[db]`** — re-derived from a live database by `scripts/s0_catalogue_verify.py`,
  with the query recorded. Re-runnable.
- **`[lit]`** — from the literature; audited in S0 against the references in
  `results/s0_baseline/`.
- **`[open]`** — a question this project answers. Not a fact yet.

Read this before S1. **Cite the review
(`docs/channel_review_2026.md`), never this file** — the review's references
are machine-resolved against Europe PMC and this file's are not.

---

## 1. What an ion channel is, and how many there are

An ion channel is a gated aqueous pore. The full scope decision — auxiliary
subunits, transporters with a channel fold, aquaporins, viroporins — is
`docs/scope_and_boundaries.md`.

**The catalogue holds 102 families in 32 superfamilies, of which 75 are
census families covering 328 human pore-forming genes.** `[db]` (Eight of those
genes — PACC1 and seven contested proposals — were added after S20 and are not
yet in any census; the census counts below are on the 320.) The remaining
27 families exist so they can be excluded: auxiliary subunits, domain-sharing
non-channels and channels that are not ion channels.

Published counts of the human channelome run from about 240 to about 400.
`[lit]` The spread is almost entirely the scope questions, not disagreement
about biology. Measured on three curated database lists (S20): GtoPdb 285,
HGNC 331, UniProt KW-0407 338, 400 in their union; the 238 genes on all three
are all pore-forming, and the differences are auxiliary subunits, aquaporins,
transporters and contested families (TMC, anoctamin scramblases, CLIC,
connexins in UniProt). `[db]`

The largest division is the **P-loop (voltage-gated-like) superfamily**: 23
catalogued families and 143 human genes — 45 % of the census on its own.
`[db]` Of those, **79 are potassium channels**, the largest single
ion-selectivity class in the genome. `[db]` The next four divisions are the
Cys-loop receptors (46), the anoctamin/OSCA/TMC clan (21), the connexins
(21) and the ionotropic glutamate receptors (18). `[db]`

## 2. The superfamilies are not one family

**There is no alignment containing a nicotinic receptor and a Kv channel.**
`[lit]` The ion channels are roughly twenty-five independent origins with at
least ten distinct pore folds: the P-loop module, the pentameric Cys-loop
bundle, the inverted P-loop of the iGluRs, the trimeric P2X dolphin fold, the
DEG/ENaC hand, the CLC double barrel, the TMEM16 groove, the Piezo propeller,
the β-barrel of VDAC, the connexin and innexin hexamers, CALHM, MscL, MscS,
otopetrin, TMEM175 and the ABC fold of CFTR.

Consequences the project is built around:

- The phylogeny is a forest plus a fold network (`docs/phylogeny_protocol.md`).
- "Ion channel" is a functional class, not a clade. Any statement of the form
  "ion channels evolved …" is either about one superfamily or is wrong. `[lit]`

## 3. Convergence is the rule, not the exception

- **Potassium selectivity evolved at least twice.** Every canonical K⁺
  channel uses the T-x-G-Y-G backbone carbonyl cage. TMEM175 is K⁺-selective
  and has no such filter — its pathway is lined by isoleucines. `[lit]`
- **Vertebrates run two unrelated gap-junction systems.** Connexins
  (`PF00029`) and pannexins (`PF00876`, the innexin model — `[db]`) build the
  same kind of large pore and are not detectably related; invertebrates have
  innexins and no connexins. `[lit]`
- **Mechanosensitivity is a mechanism, not a clade.** Piezo (38 TM per
  subunit), MscL (2 TM, 136 aa — `[db]`), MscS, TMEM63/OSCA, TMC and the
  mechanosensitive K2P channels are unrelated to one another. `[lit]`
- **The four-repeat architecture arose once and diversified.** Nav, Cav,
  NALCN, CatSper and the two-repeat TPCs are one clade, and their pore module
  is a duplicated Kv-like ancestor. `[lit]`

## 4. What the annotation databases actually say

These are the numbers that decide how a census must be run. All `[db]`,
measured 2026-08-19 unless noted; re-derive with
`scripts/s0_catalogue_verify.py`.

- **`PF00520` (Ion_trans) is carried by 20 catalogued families** — Kv, Nav,
  Cav, TRPV, TRPA, TRPM2, HCN, CNG, ITPR, RYR, and `nonchannel_vsp`, a
  phosphatase. A Pfam hit places a protein in a superfamily and never in a
  family.
- **Most TRP families carry no `PF00520` at all.** TRPV1, TRPV5, TRPA1 and
  TRPM2 do; TRPC3, TRPM8, MCOLN1, MCOLN2 and PKD2 do not — their pores are
  modelled by `PF08344`, `PF18139`/`PF25508` and `PF08016`. A census
  enumerating P-loop channels by `PF00520` alone loses most of the TRP
  division silently. This is hazard **H7** and it is the single most
  consequential annotation fact in the project.
- **Pfam's pore model does not track transmembrane count.** The 6TM SK
  channels (KCNN2) are annotated with `PF07885`, the 2TM model; so is KCNT1.
  KCNMA1 gets `PF00520`. Hazard **H8**.
- **`PF08016` covers three different things**: TRPML, TRPP and polycystin-1,
  which is not a pore. Hazard **H13**.
- **CFTR and ABCC8 have identical ABC architectures** (`PF00664`×2 +
  `PF00005`×2). One is a chloride channel, the other regulates a potassium
  channel. The only architectural difference is CFTR's R domain
  (`PF14396`). Hazard **H11**.
- **Human ZACN is annotated with `PF02931` and no `PF02932`** — a Cys-loop
  channel with no annotated transmembrane domain, which the H2 rule
  (LBD + TM required) therefore rejects. One rule, one true positive
  (AChBP) and one false negative (ZACN), both recorded.
- **113 declared Pfam accessions, 0 missing from InterPro, 2 short-name
  mismatches.** `[db]`

## 5. The hazards, in one line each

The full registry with tests is `src/catalogue/hazards.py`; nineteen entries (H17–H19 added with the S20 catalogue additions).
The five that will change a headline number:

- **H1** Nav / Cav / NALCN / CatSper share one architecture → the filter locus.
- **H3** the iGluR clamshell is the class C GPCR ligand-binding domain →
  require the pore region.
- **H7** most TRPs carry no `PF00520` → enumerate from the union of pore models.
- **H9** a voltage-sensor domain is not a channel → TPTE is a phosphatase.
- **H15** gene-symbol prefixes group pores with accessory subunits → symbols
  are never consulted; measured cost of getting this wrong, ~+15 % on the
  human count.

## 6. Clinical weight

The channelopathies are the largest Mendelian disease class attached to any
protein family. `[lit]` `CFTR` (cystic fibrosis), `SCN1A` (Dravet), `KCNQ1`
and `KCNH2` (long-QT), `RYR1` (malignant hyperthermia), `CLCN1` (myotonia),
`GJB2` (the commonest inherited deafness), `PKD1`/`PKD2` (polycystic kidney
disease), `KCNJ11` (neonatal diabetes), `SCNN1A` (Liddle syndrome). Roughly a
quarter of all approved drugs act on a channel or a channel complex. `[lit]`

## 7. The open questions this project answers

- **Q1** How many pore-forming ion-channel genes are there in a declared
  genome scope, and how does the number move when the scope questions are
  answered differently? `[open]`
- **Q2** Can a positive-test classifier reproduce the literature's family
  assignments without ever reading a gene symbol, and where does it fail? `[open]`
- **Q3** Which superfamilies are recoverable by domain search alone, and
  which are only reachable by reference or profile methods? `[open]`
- **Q4** Is the animal absence of MscS-family channels a genome fact or a
  database fact? `[open]`
- **Q5** Are the four-repeat channels' filter loci congruent with their
  pore-module phylogeny — does Cav3's EEDD sit where the tree says it should? `[open]`
- **Q6** Is the connexin/pannexin split a genuine convergence, or homology
  below the detection limit of sequence methods? `[open]`
- **Q7** Do the mechanism-based splits — CLC channel vs transporter,
  anoctamin channel vs scramblase — survive a phylogeny, or are they
  physiology cutting across clades? `[open]`
- **Q8** How complete is each superfamily's annotation across the tree of
  life, and how much of the apparent lineage-specific loss is assembly and
  gene-set quality? `[open]`
