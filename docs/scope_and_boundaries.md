# What counts as an ion channel — the scope decision

*Provenance tags follow the parent project's grammar: `[db]` verified against
a live database with the query recorded, `[lit]` literature, `[open]` a
question this project answers.*

Every census needs a denominator, and "all ion channels" does not have one
until somebody says what an ion channel is. Published counts of the human
channelome run from about 240 to about 400 — three curated databases hold
285 (GtoPdb), 331 (HGNC) and 338 (UniProt KW-0407) and 400 between them,
of which the 238 on all three are all pore-forming (**[db]**, S20) — and most
of the spread is not disagreement about biology — it is different answers to the six questions
below. This file records this project's answers so that every later number
can be read against them.

## The operational definition

> An **ion channel** is a protein assembly forming a gated aqueous pore that
> conducts ions across a membrane down their electrochemical gradient,
> without a coupled conformational transport cycle.

Four load-bearing words:

* **pore** — a continuous aqueous pathway, not a carrier site.
* **gated** — conductance changes in response to something. This is what
  excludes the ungated bacterial porins while keeping VDAC.
* **down their gradient** — no energy coupling. This is what excludes the
  pumps and the antiporters, including the CLC transporters.
* **assembly** — the unit counted is the *pore-forming subunit gene*, so a
  heterotetramer of four different genes counts four times and its auxiliary
  subunits count zero.

## The six boundary questions, and this project's answers

**1. Do auxiliary subunits count?** *No.* KCNE1, the Cav β and α2δ and γ
subunits, the Nav β subunits, SUR1/SUR2, STIM1, EMRE, barttin, the TARPs and
the CatSper auxiliaries are essential to the channels they serve and line no
pore. **Counting them all would inflate the human total by 24 %** — 71
catalogued auxiliary genes plus 5 the database lists carry that the catalogue
does not (KChIP1–4, TMEM37), against 320 pore-forming genes (**[db]**, S20;
this file said "roughly 15 %" before it was measured). The three database
channelomes count between 1.4 % (GtoPdb) and 10.7 % (UniProt KW-0407) of
their totals as auxiliaries, each a different set. They are catalogued with
status `channel_associated` so they are excluded by a recorded decision
rather than by omission. **[db]** — every one of them was resolved and its
domain architecture measured in S0.

**2. Do transporters with a channel fold count?** *No, but they stay in
their family.* CLC-3 to CLC-7 are 2Cl⁻/H⁺ antiporters in a family whose other
members are channels; the bacterial prototype ClC-ec1 is an antiporter too.
They are catalogued in the CLC superfamily with status `transporter`. This
splits membership from mechanism deliberately (**D24**): the family call and
the mechanism call are different questions with different evidence. **[lit]**

**3. Do scramblases with a channel fold count?** *Flagged, and counted.*
ANO1/ANO2 are chloride channels; ANO3–ANO10 are phospholipid scramblases
that also pass ions. Their Pfam architectures are identical (**[db]**: ANO1
and ANO6 both carry `PF04547` + `PF16178`), so no sequence-level test
separates them and the split is a literature assignment. Status
`channel_contested`; whether the split survives a phylogeny is question Q7.

**4. Do large-pore channels count?** *Yes.* Connexins, pannexins, innexins,
LRRC8 and CALHM conduct ions through a gated aqueous pore; that the pore also
passes ATP does not make it not a channel. Flagged `LARGE_PORE` so any
statement about ion selectivity can exclude them cleanly. **[lit]**

**5. Do non-ion channels count?** *No — and aquaporins are the control that
proves it.* Aquaporins are gated, selective, ancient and extremely well
characterised, and they conduct water. They are catalogued with status
`out_of_scope` and appear in the S1 decoy panel, so the boundary is enforced
by a test rather than assumed. The same rule excludes the gasdermin,
perforin and MLKL pores: those are lytic, ungated and not selective. **[db]**

**6. Do viral and prokaryotic channels count?** *Yes.* Influenza M2, HIV-1
Vpu and the coronavirus E protein are channels by every criterion above, and
KcsA, MscL, MscS, GLIC, ELIC, NavAb and ClC-ec1 are the structural prototypes
the whole field is built on — and the outgroups that make the trees rootable.
A census restricted to human genes cannot make an evolutionary claim. **[lit]**

## Contested memberships, kept and marked

The catalogue carries five families as `channel_contested`. Each is included
with the argument written down, because a census that resolves live disputes
by omission is a record of opinion.

| Family | Why contested |
|--------|---------------|
| `tweety` (TTYH1–3) | Proposed as VRAC and as the Ca²⁺-activated Cl⁻ channel; both roles went to other proteins, and recent structures show no obvious conduction path. |
| `clic` (CLIC1–6) | A soluble glutathione-S-transferase-fold enzyme that inserts into membranes and conducts; the physiological activity may be oxidoreductase. |
| `tmc` (TMC1–8) | TMC1/2 are the best candidate for the hair-cell mechanotransduction pore; no structure of a conducting state. |
| `ano_scramblase` (ANO3–10) | Scramblases with measurable ion flux. |
| `delta_glur` (GRID1/2) | Do not gate to glutamate; the Lurcher mutation shows the pore is real. |

## What the scope decision does not settle

- **Whether 320 is the right number.** The catalogue's human census count is
  a claim about a curated family list, not a measurement. S2 turns it into
  one, against a declared search space. `[open]`
- **Whether each family's boundary is right.** Splitting `ano_channel` from
  `ano_scramblase`, or `clc_channel` from `clc_transporter`, is a hypothesis
  the phylogeny can test. `[open]`
- **Selectivity.** Recorded per family from the literature, never predicted
  per sequence (hazard **H14**): Cys-loop charge selectivity flips with three
  substitutions, and AMPA-receptor Ca²⁺ permeability is set by RNA editing,
  which a genome does not show. `[lit]`
