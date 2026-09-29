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

---

## 2026-08-19 — S1: what a classifier can and cannot tell you from sequence

**The classifier rejects every decoy and finds two thirds of the positives.**
Against a panel of 72 known channels and 25 proteins chosen because they look
like channels and are not, it called **50/72 (69 %)** of the positives
correctly and rejected **25/25 (100 %)** of the decoys. Every one of the
sixteen recorded confusions was exercised at least once. Nothing was called a
channel that is not one — including the potassium-channel-domain proteins
that are ubiquitin-ligase adaptors, the glutamate-receptor-clamshell proteins
that are metabotropic receptors, the voltage-sensor protein that is a
phosphatase, and the ABC transporter that regulates a channel without
conducting.

**The gene symbol was never read.** Classification used domain composition,
selectivity-filter residues and sequence identity only. The symbol is carried
through the audit trail so a reader can see it, and the answer is identical
when it is blanked or replaced with nonsense.

**Most of the work was not done by similarity.** Of the calls made, 29 came
from domain-composition rules and 7 from reading the selectivity filter,
against 24 from nearest-neighbour identity. That distinction matters: a call
from architecture or motif works on a sequence with no close relative in the
reference set, which is the case a census exists for.

**Where it fails is where the biology is genuinely ambiguous.** The 22 missed
positives are concentrated in the families that domain annotation cannot
separate and that sequence identity separates only weakly — the Cys-loop
receptors, whose families all carry the same two domains and are ~35–40 %
identical to one another; the anoctamin channels against the anoctamin
scramblases, which have identical architecture; and the families whose only
reference is the protein being tested. These are not tuning failures. They
are the point at which sequence stops carrying the answer and profile or
structural methods have to take over.

**The benchmark found a bug in a method inherited from the parent project,
and the bug had produced a false positive.** Scoring identity only over the
columns two sequences share is right for proteins of unequal length, and
badly wrong without a coverage floor: an aligner placing a 226-residue
connexin inside a 4,967-residue ryanodine receptor selects the 226
best-matching positions, and the score then reflects exactly those.
Connexin-26 measured **61 % identical to ryanodine receptor 2 and 50 % to
connexin-43, its actual relative** — and was duly called a ryanodine
receptor. Requiring the alignment to cover 30 % of the longer sequence fixed
it, and specificity went from 24/25 to 25/25.

**Two things became visible only once they were drawn.** The potassium
filter alignment (review fig. 2) puts *Streptomyces* KcsA, *Drosophila*
Shaker and human Kv1.1 in the same five columns, and the two rows that carry
the argument are the exceptions: *Bacillus* NaK, whose TVGDG differs from
the signature by one residue and is not potassium-selective, and GluR0 — a
cyanobacterial **glutamate receptor** — carrying a textbook TVGYG. The
second is the four-repeat locus (fig. 3), where DEKA, EEEE, EEDD and EEKE
sit in visibly similar sequence context: the residues that decide which ion
a channel conducts are four positions in a stretch that otherwise looks the
same.

*(pending: S2 — whether these recall figures hold when the panel is the whole
census rather than one member per family, and S3 — whether profile methods
recover the families sequence identity cannot.)*

---

## S2 — the census, uncapped (2026-09-28)

**Just over 1.2 million proteins in UniProt carry some recognisable piece of
an ion channel, and the number is exact.** Every one of the 1,245,200 was
fetched and counted. The query was split by taxonomy into pieces that must
add back up to the whole, and they do to the last protein. Each of the 67
domain signatures matches UniProt's own count for it.

**The signature union finds essentially every human channel.** 319 of the
320 human pore-forming genes are enumerated (GLRA4 is the exception). That
matters because the obvious shortcut would have silently dropped most of the
TRP channels, whose pores Pfam models under other names. Searching by the
Ion_trans pore model (`PF00520`) alone is that shortcut.

**Domain composition plus the selectivity filter names the family for about
a quarter of the census, and almost never names the wrong one.** 346,627
records get a family call without any nearest-neighbour comparison. 31,366
of those come from reading the four residues of the selectivity filter:
EEEE and EEDD for calcium channels, DEKA for sodium channels, EEKE for
NALCN. Among the 320 human genes the classifier is wrong exactly once
(ZACN, a known failure of the AChBP rule). Where it cannot decide, it stops
at the superfamily rather than guessing.

**Where it stops is where the biology says it should.** The ties are between
families that carry identical domains: the Kv1 channels and the silent Kv
modifiers, ASIC, ENaC and the invertebrate DEG channels, anoctamin channels
and anoctamin scramblases, and metazoan and non-metazoan P2X. Also every
Cys-loop and glutamate receptor, whose families all share the same two
domains. Telling these apart needs profile methods *(pending: S3)*.

**Much of the "unassigned" 40 % is not channels at all, and that is a
finding about domain databases.** Three signatures the catalogue uses to
recognise channels sit mostly on other proteins. The cyclic-nucleotide
binding domain of HCN and CNG channels is carried by 148,000 proteins with
no channel domain, mostly bacterial transcription regulators. The PAS domain
of hERG is on 112,000 sensor proteins. The binding-protein domain of the
bacterial glutamate receptor GluR0 is on 116,000 periplasmic transporters.
A channel's accessory domains are a poor handle on the channel.

**Some channel families are more varied than their textbook architecture.**
65,506 bacterial and archaeal proteins carry the core MscS mechanosensitive
channel domain but not the full three-domain set the catalogue declares. The
same partial architectures turn up in OSCA, TRPM, Piezo, Slo and RyR. How
many of these are working channels is open *(pending: rule revision
benchmarked on S1, then S3)*.

**Channels turn up where the catalogue expected them, plus a few surprises
to check.** MscL and MscS are overwhelmingly bacterial; TMEM175 is mostly
bacterial with a eukaryotic minority; K2P, Cav, Nav and Kir are eukaryotic.
39 viral records reach the P-loop superfamily; that they are Kcv-type K⁺ channels is expected but unchecked *(pending: S21)*. About a thousand eukaryotic
MscL records and 51 bacterial Kir records are unexpected and need checking
*(pending: S21)*.

---

## S3a — one profile per family resolves most of what domains could not (2026-09-28)

**Families that share every domain can still be told apart by sequence, and
mostly are.** A statistical profile of each of the 90 catalogued families
separates the receptors that S2's domain rules could only place in a
superfamily: 94 % of the Cys-loop receptors now fall into nicotinic, GABA-A,
glycine or 5-HT3 families; 98 % of the anoctamins into channels or
scramblases; 95 % of P2X; 82 % of glutamate receptors. Overall, 58.7 % of the
1.25 million candidate proteins now have a family, up from 27.8 %. How
reliable those splits are outside vertebrates is not yet measured *(pending:
non-vertebrate test panel, emergent)*.

**The profiles and the domain rules agree when both speak.** On 320,000
proteins both methods called, they name the same family 96.8 % of the time
(99.2 % in bacteria) — two different kinds of evidence converging. On a
held-out test where each protein was removed from its own profile, 69 of 71
known channels were placed correctly, and no decoy was called a channel.

**Where they disagree, the fault is mostly in rules that test for an
absence.** Three of S2's rules decide a family by what a protein *lacks*.
One calls any Cys-loop binding domain without an annotated membrane domain
an "AChBP" (the soluble snail protein) — but 6,600 of those proteins are
receptor-length, and the profiles call them receptors. Another calls any
fragment with the domain IP3 receptors and ryanodine receptors *share*, but
without the ryanodine-specific domains, an IP3 receptor; comparison with the
IP3R project's independent census found 502 such fragments that project
calls ryanodine receptors. A third splits the TRPP channels from
polycystin-1 by counting repeats and misplaces ~3,000 polycystin-1-like
proteins. None of these is corrected yet *(pending: rule rewrites,
emergent)* — the census records them as conflicts or flags them rather
than choosing.

**Hv1-like voltage-sensor proteins are more widespread than the domain rules
could see** — nearly 2,000 proteins, mostly Hv1-length (median 252 aa), that S2 only placed in the
P-loop superfamily are called Hv1 by the profile *(pending: the PF00520
superfamily-conflict fix; not yet verified as proton channels)*.

**Two independent projects agree on the IP3 and ryanodine receptors.** On
15,600 proteins both censuses contain, the same family is called for 12,200,
and the largest disagreement traces to one identifiable rule.

---

## S2b — rules that decide by absence, replaced (2026-09-28)

**Two pairs of channel families cannot be told apart by their domains at
all, and the census now says so instead of guessing.** IP3 receptors carry
no domain that ryanodine receptors lack; a protein fragment carrying only
their shared N-terminal domain could be either. Likewise, nothing in a
Cys-loop binding domain's annotation, size or lack of membrane helices
separates the soluble snail acetylcholine-binding protein from a receptor
gene model that stops before its membrane region — a test built on those
features put three quarters of its "AChBPs" in insects, worms and
vertebrates, which have none. These proteins now stop at the superfamily
unless a sequence-level method names them.

**Where a positive marker exists, it works.** Ryanodine receptors are
recognised by their own repeat or membrane region, polycystin-1 relatives by
their PLAT, REJ or GPS domains, and a large subset of TRPP channels by a
C-terminal domain found on no other family. Checked against the independent
profile method, those calls agree 99.9–100 %.

**The census got more honest, not smaller.** The profile method still names
most of these proteins, so 58.8 % of records carry a family. What changed
is that the calls resting on a missing domain are gone: conflicts between
the two methods fell from 12,857 to 3,288, the census no longer calls any
ryanodine-receptor fragment an IP3 receptor where the IP3R project's
independent census disagrees, and 319 of 320 human channel genes are placed
correctly.

---

## S2c — the last absence rules, and a lesson about what agreement means (2026-09-28)

**KCTD proteins are recognisable by their own tails.** The T1 domain that
KCTD proteins share with Kv channels says nothing on its own, but most KCTDs
also carry one of five clade-specific C-terminal domains, and on 14,354
proteins those domains and the independent profile method agree on all but
three. A T1 domain with neither a pore nor a KCTD tail could be either — a
KCTD from an uncatalogued clade or the front end of a Kv channel gene model
— and the census now leaves it at the superfamily.

**No domain marks a sulfonylurea receptor.** SUR1 and SUR2 carry only the
generic ABC-transporter domains shared with over a million other proteins.
The rule that called them "SUR" by the absence of CFTR's R domain had in
fact caught 918 bacterial peptide-exporting ABC transporters and 62
eukaryotic gene models fusing ABC domains to other proteins — not one
receptor. The rule is gone.

**Two methods agreeing is evidence only if the right answer was available
to both.** The profile method "confirmed" 528 of those bacterial
transporters as SUR, because the catalogue gives it no generic ABC
transporter to prefer. Agreement between classifiers measures something only
when every plausible family is on the menu *(pending: a decoy family,
emergent)*.

---

## S3a2 — giving the profiles the right wrong answer (2026-09-28)

**With a generic ABC transporter on the menu, nothing bacterial is called a
sulfonylurea receptor any more.** A decoy profile built from human ABCC
transporters and bacterial peptide exporters now claims the ~1,000
bacterial ABC transporters the receptor profile had been winning, while the
real receptor (SUR1) and CFTR still beat it by more than two to one in
score. The channel census does not change — these proteins were never
counted as channels — but a classification that names a family should name
the right one, and one that can only choose among the families it knows
needs the plausible wrong answers catalogued too.

---

## S4 — what "absent" will mean (2026-09-28)

**Every claim that a lineage lacks a channel family will now be a claim
about 52 named organisms and nothing wider.** The panel runs from *E. coli*
and an archaeon through yeasts, plants, ciliates, the malaria parasite, the
animals' single-celled relatives, a sponge, a placozoan and a sea anemone, to
fourteen vertebrates and three viruses. For 50 of them the census will
search UniProt's reference gene set: 822,499 genes in total, one protein per
gene, from the same database release the census was built from. The other
two are the marbled electric ray, whose electric organ gave biochemists
the nicotinic receptor, and the garden snail that carries the
peptide-gated FaNaC channel. Neither has a reference gene set, only a
genome, so a channel in either can only be found by searching its DNA
*(pending: S5)*.

**The gene sets are not equally good, and the table says which are weaker.**
Fourteen of them miss more than a tenth of the genes every organism in
their lineage is expected to have. They include the choanoflagellate
*Monosiga* and the sea squirt *Ciona* (both under 80 %), the lamprey, the
platypus and rice. An apparent loss in one of these is weaker evidence than
the same loss in mouse or fly *(pending: S5, which decides absences against
the genome)*.

**Viral channels were never searched, and this was found by accident.** The
influenza M2 proton channel is in the influenza gene set but not in the
census, because the census never searched the protein signatures that
define the viral channels. It is the only channel family with this gap. The
profile sweep will cover these three viruses *(pending: S3b)*. A fix to the
census itself is an open task.

## S3b — what searching by protein signature misses (2026-09-29)

**Searching by protein signature finds almost every ion channel in the 50
organisms with a reference gene set, but not all of them.** The census was
built by looking for the protein segments that databases annotate as
channel pores. Searching the same organisms' complete gene sets with
whole-family models finds 8,755 channels. 230 of them (2.6 %) carry no
annotated pore segment at all, and 96 of those are confident. The losses
are not random. **The voltage-gated proton channel Hv1** is the worst case
among real channels: 25 of the 33 Hv1 proteins in the panel, in snails, a
limpet, moss, fish, a sea anemone and *Trichoplax*, lack the annotation the
census searched for. That fits the earlier finding that the Hv1 domain
model was built from mammals. Chloride intracellular channels (CLICs) in
plants, ciliates and invertebrates follow, and so do all three viral
channels: HIV Vpu, influenza M2 and the SARS-CoV-2 envelope protein. The
miss rate is 1 % in vertebrates, 3 % in invertebrates and 13 % in plants.
Fungi and the malaria parasite show none, but they have few channels to
miss.

**One apparent discovery is an artefact, and the report says so.** Thirty-
seven "new" TRPN (NOMPC) channels turned up, including in human and mouse,
which have no TRPN. They are ankyrin-repeat proteins: the TRPN model is
mostly its long ankyrin-repeat spring, so any protein made of ankyrin
repeats matches enough of it to pass the coverage test. The same happens
with leucine-rich repeats for LRRC8. These are counted as an upper bound,
not as channels, and fixing the test is an open task.

**Iterative search from a single channel reaches its whole superfamily in
one step, not just its family.** Starting from one sodium channel, the
first round already includes calcium channels, TPCs and CatSper, and the
second round adds more of them. Starting
from any Cys-loop receptor subunit (nicotinic, GABA-A, glycine, 5-HT3,
ZAC, invertebrate or bacterial), the search converges on the same ~1,625
receptors. The same holds for the DEG/ENaC (490) and P2X (125) families.
In the 24 searches that converged cleanly, the iterative search recovers
3,352 of the 3,356 family members the models call (99.9 %), so the census
is not missing members that an unconstrained search would find. That
claim is limited to those 24 families *(pending: S5 for absences, S15 for
the per-superfamily recall curve)*.

## S5a — what it takes to say a gene is missing from a genome (2026-09-29)

**A gene can only be declared absent where a close relative was used to
look for it.** Searching genomes directly with protein "baits" finds a
channel reliably when the bait comes from the same broad lineage, and
barely at all when it doesn't. Seven animal and ciliate two-pore channels
put no alignment at all on *Arabidopsis*'s well-known TPC1 gene, by either
of the two search methods. With a bait from every panel species, and each
genome searched only with baits from *other* species, the search detects
172 of 173 channel families that a genome's own gene set says are there.
So every "absent" here is conditional on an in-lineage bait existing.
Where one doesn't, the result is recorded as unjudgeable, not as absent.

**In the six pilot genomes with gene sets, the gene sets are not hiding
channels.** No informative gap in the gene sets turned out to be an intact,
confidently identified channel gene sitting in the genome unannotated.
Ten absences survive every check. Among them are the long-reported ones:
no voltage-gated sodium channel in *C. elegans*, no P2X receptor in
*C. elegans* or *Drosophila*, and no two-pore channel in either. None of
the four literature expectations was contradicted *(pending: reference
verification before any is cited; S5b for the other 45 genomes)*.

**The apparent mouse ZAC is not a gene.** The zinc-activated channel is
reported absent from rodents. The genome search did return a mouse "ZAC"
model, but it is 32 % identical to the bait, stretched over 400 kb, and on
a different chromosome from the one that should carry the gene. It is
recorded as weak evidence, not a find. Mouse likewise carries 12
VDAC-like loci for three VDAC genes, most of them broken copies (retrocopies).
Counting genomic hits would overstate copy number.

**Mammalian channel genes are enormous.** At the loci the search called,
13 of 220 mouse genes have an intron longer than 200 kb (the aligner's
default limit), up to almost 1 Mb in the acid-sensing channel ASIC2.
At the default, a search would report these genes as fragments.

**The garden snail, which has no gene set, carries 46 channel families by
genome alone,** including Hv1 and CLIC
and the full molluscan complement of P-loop, Cys-loop and iGluR channels
*(pending: S5b for the snail's own ledger and for *Torpedo*)*.

## S5b — what 52 genomes add to the proteome census (2026-09-29)

**The gene sets are nearly complete for channels, and where they are not,
the gap can be large.** Across all 52 genomes the search recovered 1,317 of
1,326 channel families each genome's own gene set says are present, with
that species' own sequences excluded from the search. Against that control,
only 11 cells show an intact, confidently identified channel gene that the
reference proteome lacks. The most striking: **the pufferfish *Takifugu*
reference proteome contains no ryanodine receptor at all**, while its genome
carries six intact RyR genes at 61–80 % identity to other vertebrate RyRs.
Also found: a NALCN in *Ciona*, TRPN (NOMPC, the fly's mechanotransduction
channel) in *Trichoplax*, *Nematostella*, *Daphnia* and *Lottia*, and a
second P2X receptor of the non-animal type in the choanoflagellate
*Monosiga* *(pending: S16 for which are annotation gaps and which are
pseudogenes)*.

**An MscS in a sea anemone is probably not a sea-anemone gene.** The
*Nematostella* MscS loci sit on contigs of 6.7 and 3.1 kb, the usual mark
of bacterial contamination in an assembly. Whether any animal has MscS is
one of this project's open questions, and this is not yet an answer
*(pending: S10)*.

**72 absences survive every check** — an in-lineage bait, a genome whose
control passes, and contiguity sufficient for the family's typical gene
size. The long-reported absences hold: no voltage-gated sodium channel in
*C. elegans* or in the sponge *Amphimedon*, no P2X in *C. elegans* or
*Drosophila*, no ENaC in zebrafish or *Takifugu* (teleosts are reported to
have lost it), no Hv1 in *C. elegans* or *Drosophila*, and no ZAC in rat
(mouse has only a broken ZAC-like model) *(pending: reference verification)*.
ZAC also reads absent in eight non-mammalian or non-therian vertebrates, and
*Daphnia*'s old, fragmented assembly accounts for five absences including
Piezo — both need checking before they become claims *(pending: S10)*.

**Genome-only species.** The garden snail *Cornu* carries 48 channel
families by genome alone; the marbled electric ray *Torpedo*, whose electric
organ is the classical source of the nicotinic receptor, carries 56.

**Some channel genes are larger than the search's limits, and are still
found.** Human ASIC2 has an intron of 1.04 Mb, beyond the 1 Mb the aligner
was allowed. It was still recovered, as were four smaller-genome cases.
Allowing longer introns made results worse, not better: the aligner began
joining neighbouring genes into false ones.

## S6 — lining the channels up (2026-09-29)

**The part of every channel that forms the pore can be found the same way in
all of them, including genes nobody has annotated.** Within each of the five
groups of channels whose members can be compared by sequence at all —
the voltage-gated-like channels, the glutamate receptors, the
calcium-release channels, the gap-junction-like innexins/pannexins and the
proton channel — the pore was cut out of 5,228 chains by one method. Where
UniProt records the pore's position, the cut matches it almost exactly
(median overlap 97 %), and in potassium channels the signature filter
sequence (GYG) sits inside the cut piece 96.6 % of the time; most of the
exceptions have no canonical filter at all.

**For six families no database says where the pore is** — CNG channels,
TRPC, NOMPC, IP3 receptors and the plant and bacterial glutamate receptors.
Their pore was located by learning what pores look like from the annotated
families, a method measured to be within a dozen residues of the truth on
every annotated relative. For the IP3 receptor this needed a correction:
learned only from ryanodine receptors, the method found the IP3 receptor's
pore helix and filter but not the helix before them *(pending: S8, check
against the IP3 receptor structure)*.

**One pore detector does not fit a whole superfamily.** A single model of
"the P-loop pore" finds the pore in Kv, Nav and Cav channels it was not
trained on, but not in inward rectifiers, HCN or TRPM channels, and not in
the innexin clan — the pore modules within one superfamily have diverged
further than a shared model can bridge.

**Some families are too divergent to keep much alignment.** In the largest
families — two-pore-domain K⁺ channels, nicotinic receptors, CNG channels —
automatic trimming keeps only 3–4 % of aligned columns; whether that is
enough to resolve their trees is S7's first question *(pending: S7)*.
