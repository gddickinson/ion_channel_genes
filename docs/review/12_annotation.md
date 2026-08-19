## 12. What the databases record

Everything above is drawn from the primary literature. This section is
different: it reports what the annotation databases actually contain about
these proteins, measured directly against InterPro and UniProt on
2026-08-19 with the queries recorded in `results/s0_baseline/`
[paysanlafosse2023, mistry2021, uniprot2023]. It is included because every
computational statement about ion channels — every census, every "how many
are there", every claim of absence in a lineage — passes through these
resources, and their coverage is not uniform in ways that are rarely stated.

The measurements below come from a catalogue of 90 families in 25
superfamilies, of which 68 are ion-channel families covering 320 human
pore-forming genes, with 115 declared domain accessions and 162 reference
proteins. All 115 accessions and all 162 exemplars resolved.

### 12.1 A pore model is not a family, and not always a pore

`PF00520` (Ion_trans) is the canonical model of the P-loop pore module. In
this catalogue it is carried by **20 different families**, which is expected
— it is a superfamily-level model. Two of those families are the
informative ones:

- **HVCN1**, the voltage-gated proton channel, carries it and *has no pore
  domain at all*: its conduction pathway is the voltage sensor
  [ramsey2006, sasaki2006].
- **TPTE**, a voltage-sensing phosphatase, carries it and is not a channel
  [murata2005].

So the presence of the canonical pore model is not evidence that a protein
has a pore, let alone which family it belongs to.

### 12.2 Most TRP families carry no copy of it

Measured across exemplars: TRPV1, TRPV5, TRPA1 and TRPM2 carry `PF00520`;
**TRPC3, TRPM8, MCOLN1, MCOLN2 and PKD2 carry none of it.** Their pore
regions are covered by different models entirely — `PF08344` for TRPC,
`PF18139`/`PF25508` for TRPM, `PF08016` for the mucolipins and polycystins.

A census that enumerates P-loop channels from `PF00520` therefore loses most
of the TRP division without any error message. This is the single most
consequential annotation fact we measured, and it is not a database error:
the models were built family by family from the best-characterised members,
and behave accordingly.

![](figures/fig5_domain_matrix.png)

**Figure 5 | The shared signatures.** Presence of each widely shared domain
accession across the catalogued families, ordered by how many families carry
it. Blue marks an ion-channel family and red one the catalogue does not count
as a channel. The first column is `PF00520`, reaching twenty families and
including two red rows. Rendered from
`results/s0_baseline/exemplar_architecture.tsv`.

### 12.3 Domain models do not track topology

The small-conductance calcium-activated potassium channels are 6TM proteins,
and InterPro annotates KCNN2's pore with `PF07885` — the 2TM model. So does
KCNT1, while its close relative KCNMA1 gets `PF00520`. Any rule that infers
topology from which pore model matched is wrong for these families;
topology has to come from transmembrane predictions or experimental
annotation.

Similarly, the inward rectifiers carry **no** `PF07885` at all — they are
covered by `PF01007` and two cytoplasmic-domain models. The textbook
grouping of Kir with the two-pore-domain channels as "2TM/1P potassium
channels" has no counterpart in the annotation.

### 12.4 Architecturally identical families

Four cases where domain composition cannot distinguish families that are
functionally distinct:

| families | shared architecture | what separates them |
|---|---|---|
| Nav, Cav, NALCN, CatSper | `PF00520`×4 | the four filter residues (§3.2) |
| CFTR vs ABCC8/9 | `PF00664`×2 + `PF00005`×2 | `PF14396`, CFTR's R domain |
| TMEM16A/B vs the scramblases | `PF04547` + `PF16178` | nothing known |
| TRPML vs TRPP vs polycystin-1 | `PF08016` | `PF21381` vs `PF20519` vs PKD repeats |

The third row is the uncomfortable one: ANO1 and ANO6 have identical
architectures, ~40 % identity and different functions, and no
sequence-level test yet separates them.

![](figures/fig8_architecture_traps.png)

**Figure 8 | What domain composition cannot separate.** Scale diagrams from
measured InterPro coordinates. *Top:* four channels with the same four
copies of `PF00520` and three different permeant ions — separated only by
the residues in figure 3. *Second:* CFTR and the sulfonylurea receptor,
which differ by one accession, `PF14396`. *Third:* a chloride channel and a
lipid scramblase carrying the same two accessions in the same order. *Fourth:* `PF08016`
across TRPML, TRPP and polycystin-1, which is not a pore. *Bottom:* four
pairs in which the shared domain exists outside channels entirely — the Kv
T1 domain in a ubiquitin-ligase adaptor, the iGluR clamshell in a
metabotropic receptor, the Cys-loop ligand-binding domain in a soluble snail
protein, and the pore-module annotation on a phosphatase. Rendered from
`results/s0_baseline/domain_positions.tsv`.

### 12.5 Signatures shared with non-channels

Of 26 domain accessions carried by more than one family in this catalogue,
**15 cross the channel / non-channel boundary**. Beyond `PF00520`:

- `PF02214`, the Kv T1 tetramerisation domain, is a generic BTB/POZ domain
  carried by 25 human KCTD proteins that are mostly ubiquitin-ligase
  adaptors.
- `PF01094`, the iGluR clamshell, is the ligand-binding domain of every
  class C GPCR (§6.2).
- `PF02931`, the Cys-loop ligand-binding domain, is a complete soluble
  protein in AChBP [brejc2001].
- `PF02815`, the MIR domain of the intracellular release channels, is
  carried by the protein O-mannosyltransferases.

### 12.6 The outgroups are the worst-covered proteins in the subject

The finding with the most consequence for phylogenetics: **three of the four
prokaryotic outgroups that root the major superfamily trees cannot be found
by the domain rules that define their own superfamilies.**

- **GLIC and ELIC**, which root the Cys-loop tree, carry `PF02931` and no
  `PF02932` — so a rule requiring the transmembrane domain alongside the
  ligand-binding domain, which correctly rejects AChBP, also rejects them
  [hilf2008, bocquet2009].
- **GluR0**, the cyanobacterial glutamate receptor that roots the iGluR
  tree, is annotated as `PF07885` + `PF00497` — *a potassium-channel pore
  fused to a bacterial solute-binding protein* — with none of the eukaryotic
  iGluR models [chen1999]. That is the inverted-P-loop story visible
  directly in the annotation, and it means no iGluR-signature search will
  ever find it.
- **NavAb**, the prototype prokaryotic sodium channel, carries neither
  `PF06512` nor `PF11933`, the two Nav-specific models: it is the pore and
  nothing else [payandeh2011].

Seven of 162 reference proteins have **no gene name in UniProt at all**, and
they are almost all prokaryotic or invertebrate: NavAb, ELIC, GluR0, FaNaC,
the *Bacillus* NaK channel, the *Klebsiella* bestrophin and AChBP. One,
GLIC, was invisible to a symbol-based search because its entry is filed
under a bacterial strain rather than a species.

### 12.7 What follows

None of this argues against using these resources; there is no alternative,
and their coverage of well-studied families is excellent. It argues for
three specific practices. Enumerate from the union of pore models rather
than the canonical one. Treat a domain hit as evidence at the level where
the model is diagnostic and no further. And expect the proteins that root
your trees to be the ones your search cannot find.
