# FINDINGS.md — the biological story, task by task

A dated entry after every completed task: what it changed about the biology,
in plain language, findings before methods. Unconfirmed claims are marked
*(pending: which task confirms it)*. Paths and commands stay out of this
file — they live in `SESSION_LOG.md`.

---

## 2026-08-19 — setup: the catalogue, and what the databases say about it

**The subject has a size.** Ion channels are not one family and not one
number. Written out properly — every family that a search for ion channels
will return, with a recorded decision about whether it counts — the subject
comes to **90 families in 25 superfamilies**, of which 68 are ion channels
(**320 human pore-forming genes**) and 22 exist so that they can be
excluded: auxiliary subunits, proteins that carry a channel's diagnostic
domain without being channels, transporters built on a channel fold, and
channels that conduct something other than ions.

The subject is also lopsided in a way the textbooks flatten. The
voltage-gated-like superfamily accounts for 143 of the 320 genes — nearly
half — and 79 of those are potassium channels. The remaining 177 genes are
spread over twenty-four superfamilies that share no ancestor with it and,
mostly, none with each other.

**Twenty-five origins, not one.** There is no alignment that contains a
nicotinic acetylcholine receptor and a potassium channel. They share no
detectable ancestry, no fold, and no alignable position — and the same is
true of the P2X receptors, the chloride channels, the anoctamins, Piezo, the
connexins and about twenty other groups. "Ion channel" is a description of
what a protein does, not of where it came from. Everything else in this
project follows from that: there can be no single tree, and any statement
beginning "ion channels evolved…" is either about one superfamily or is
wrong.

**The domain annotations do not divide the subject the way the textbooks
do.** This was measured, not assumed, against InterPro on 2026-08-19:

- The pore model `PF00520` is carried by **20 of the catalogued families**,
  including a phosphatase (TPTE) that is not a channel at all, and including
  Hv1, which is a channel with no pore domain. A domain hit places a protein
  in a superfamily and never in a family.
- **Most TRP channels carry no `PF00520`.** TRPV1, TRPV5, TRPA1 and TRPM2
  do; TRPC3, TRPM8, the two mucolipins and polycystin-2 do not — their pores
  are modelled by three other accessions entirely. A census that enumerates
  voltage-gated-like channels by the obvious accession loses most of the TRP
  division without any error message.
- **Sodium, calcium, NALCN and CatSper channels are architecturally
  identical.** Four copies of one domain in one chain, nothing else to tell
  them apart. What separates them is four residues — one contributed by each
  repeat — and reading those four residues off a reference alignment
  recovers the right answer for all six four-repeat channels tested,
  including the T-type calcium channels, whose filter is genuinely different
  from the other calcium channels'.
- **CFTR and the sulfonylurea receptor have the same architecture.** One is
  a chloride channel; the other is an ABC transporter that regulates a
  potassium channel without conducting anything. One accession separates
  them.
- **A quarter of the shared domains cross the channel / non-channel
  boundary.** The potassium-channel tetramerisation domain is carried by 25
  proteins named after it that are ubiquitin-ligase adaptors; the
  glutamate-receptor clamshell is the ligand-binding domain of every
  metabotropic glutamate receptor; the calcium-release-channel MIR domain is
  carried by two glycosyltransferases.

**Two annotations are simply wrong for their families, and both were found
by asking.** Potassium channels of the SK family are 6-transmembrane
proteins annotated with the 2-transmembrane pore model; inward-rectifier
channels do not carry the 2-transmembrane model at all, although the
literature groups them with the two-pore channels that do. A rule written
from the textbook description would have failed on every member of both
families.

**One human channel breaks its own family's rule.** The zinc-activated
channel ZACN is annotated with the pentameric receptors' ligand-binding
domain and with no transmembrane domain — so the test that correctly rejects
the soluble snail acetylcholine-binding protein also rejects a real human
channel. The test keeps its priority and the exception is recorded, because
a rule quietly relaxed to fit one protein has stopped being a test.

**The structural prototypes are the hardest proteins in the catalogue to
look up.** Of 162 reference proteins, thirteen could not be found by gene
symbol on the first attempt, and they were almost all prokaryotic or
invertebrate: the archaeal MthK channel, the bacterial NaK and NavAb
channels, the *Gloeobacter* and *Dickeya* pentameric channels, the
cyanobacterial glutamate receptor, the snail acetylcholine-binding proteins.
Seven of them have **no gene name in UniProt at all**; one was invisible
because its entry is filed under a bacterial strain rather than a species.
All are now resolved and their sequences cached, but the lesson survives the
fix: the proteins whose structures the whole field is built on are the ones
a name-based search finds last.

**The outgroups are the hardest part.** Every superfamily tree needs a
prokaryotic or invertebrate outgroup to be rootable, and three of the four
main ones cannot be found by the rules that define their own superfamily.
The bacterial pentameric channels GLIC and ELIC carry the ligand-binding
domain of the nicotinic receptors and no transmembrane annotation, so the
test that correctly rejects a soluble snail protein rejects them too. The
cyanobacterial glutamate receptor GluR0 is annotated as a **potassium
channel pore fused to a bacterial nutrient-binding protein**, with none of
the eukaryotic glutamate-receptor domains — which is the "the glutamate
receptor pore is an upside-down potassium channel" story visible directly in
the annotation, and also the reason a glutamate-receptor search will never
find it. The prokaryotic sodium channel NavAb has neither of the two domains
that define the sodium-channel family: it is the pore and nothing else.

None of this is a mistake in the databases. It is what it looks like when a
family is defined from its best-studied members and then asked about its
own ancestors.

**The forest works, and it says something.** As an end-to-end check the
pipeline was pointed at the Cys-loop superfamily — the one that is alignable
end to end — and asked for a rooted maximum-likelihood tree of its eight
reference proteins. Fifty-two seconds later, rooted on the two bacterial
channels:

- the two nicotinic α subunits pair at 100 % support, the serotonin
  receptor joins them, and the soluble snail acetylcholine-binding protein
  sisters that whole cation-selective group;
- the GABA-A and glycine receptors form their own clade at 100 %.

The anion-selective and cation-selective receptors separate cleanly, and the
protein that is not a channel at all falls exactly where its ligand-binding
chemistry says it should. This is eight sequences and not the real tier-2
tree — that is S8's job, with the full census behind it — but it is the
pipeline producing a correct answer to a question with a known answer, which
is what a smoke test is for.

*(pending: S1 — whether the classifier's three tiers recover the
literature's family assignments, and which tier does the work.)*
