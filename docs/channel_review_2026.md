# The ion channels: one function, twenty-five origins

*A review of the ion-channel superfamilies, their folds, their evolution and
what the sequence databases do and do not record about them.*

**Generated file — do not edit.** The text lives in `docs/review/*.md` and is
assembled by `scripts/s0_review_build.py`. Citations are written as stable
keys and renumbered on build; every reference was resolved against Europe PMC
by `scripts/s0_review_refs.py`, and a cited key with no verified record is a
build error rather than a silent entry in the bibliography.

---

## Summary

An ion channel is a gated aqueous pore. That is a description of what a
protein does, not of where it came from, and the distinction organises
everything that follows. The proteins the term covers do not form a family,
a superfamily, or any other clade: they are on the order of twenty-five
independent inventions, built on at least ten unrelated folds, that arrived
at the same solution to the same physical problem. There is no alignment
containing a nicotinic acetylcholine receptor and a potassium channel, and
therefore no tree.

This review is organised around that fact. We take the folds one at a time —
the P-loop superfamily that supplies nearly half the human channel genes, the
four unrelated ligand-gated superfamilies, the anion channels, the large
pores, the mechanosensors and the channels of intracellular membranes — and
ask of each what is genuinely shared and what merely looks shared. We then
turn to the two questions a modern census has to answer before it can count
anything: what the annotation databases actually record about these
proteins, and where their coverage fails.

The second question turns out to be sharper than expected. Domain models
were built family by family from the best-studied members, and they behave
accordingly: the canonical pore model `PF00520` is carried by twenty
catalogued families including a phosphatase that is not a channel, while most
of the TRP families carry no copy of it at all. Three of the four
prokaryotic outgroups that root the major superfamily trees cannot be found
by the domain rules that define their own superfamilies. These are not
database errors. They are what a set of models built from well-studied
proteins looks like when it is asked about the proteins nobody studied.

## Key points

- **"Ion channel" is a functional class, not a clade.** Any statement of the
  form "ion channels evolved…" is either about one superfamily or is wrong.
- **The human channelome is lopsided.** Of ~320 pore-forming genes, 143 are
  in the P-loop superfamily and 79 of those are potassium channels; the
  remaining 177 are spread across twenty-four superfamilies.
- **Convergence is the rule.** Potassium selectivity evolved at least twice;
  vertebrates run two unrelated gap-junction systems; mechanosensitivity has
  been invented at least six times.
- **Selectivity and gating are not recoverable from sequence family.**
  Cys-loop charge selectivity flips with three substitutions; AMPA-receptor
  calcium permeability is set by RNA editing.
- **Domain annotation does not partition the subject the way the textbooks
  do**, and a census built on the obvious accessions silently loses whole
  divisions.

## 1. What is being counted

The modern study of ion channels begins with a quantitative model that
required no molecules at all. Hodgkin and Huxley described the squid axon
action potential as the sum of independent, voltage-dependent sodium and
potassium conductances, and the description was so complete that the
proteins responsible could be inferred long before any was isolated
[1]. Single-channel recording made the discrete conducting unit
visible directly [2], and molecular cloning then produced the
sequences: the electric-ray acetylcholine receptor [3], the eel
sodium channel [4], the *Drosophila* Shaker potassium channel
[5,6]. Within a decade the field had a structure —
KcsA, at 3.2 Å, showing exactly how a filter of backbone carbonyls
discriminates potassium from sodium [7].

That history explains a persistent difficulty. Each family was cloned,
named and modelled on its own, and the resulting nomenclature and annotation
reflect the order of discovery rather than any underlying structure. The
term "ion channel" was never a phylogenetic claim; it was a description of
behaviour that happened to be shared.

### 1.1 A definition, and its six edges

Throughout this review an ion channel is:

> a protein assembly forming a **gated aqueous pore** that conducts ions
> across a membrane **down their electrochemical gradient**, without a
> coupled conformational transport cycle.

Every word does work, and each marks a boundary that published counts of the
human channelome — which range from roughly 240 to roughly 400 — draw
differently:

1. **Auxiliary subunits.** KCNE1 is 129 residues with one transmembrane
   helix and no pore, and it converts KCNQ1 into the cardiac IKs channel.
   The Cav β, α2δ and γ subunits, the Nav β subunits, the sulfonylurea
   receptors, STIM1 and barttin are likewise essential and likewise not
   pores. Counting them inflates the human total by roughly 15 %.
2. **Transporters on a channel fold.** Half the CLC family are 2Cl⁻/H⁺
   antiporters, including the bacterial prototype [8]. They are
   CLCs by descent and not channels by mechanism.
3. **Scramblases on a channel fold.** TMEM16A and TMEM16B are chloride
   channels; most of their relatives are phospholipid scramblases that also
   pass ions [9].
4. **Large pores.** Connexins, pannexins, LRRC8 and CALHM conduct ions
   through a gated aqueous pathway wide enough to pass ATP. That the pore is
   wide does not make it not a channel.
5. **Non-ion channels.** Aquaporins are gated, selective, ancient and
   exhaustively characterised, and they conduct water. They mark the outer
   edge of the definition and are excluded by it.
6. **Viral and prokaryotic channels.** Influenza M2, the bacterial
   mechanosensitive channels and KcsA are channels by every criterion, and
   they are the structural prototypes on which the rest of the field is
   built.

None of these six answers is obviously right. What matters is that a count
without them stated is not interpretable, and that the things excluded should
be excluded by a recorded decision rather than by never having been
considered.

### 1.2 What a census has to survive

Two further problems are specific to counting this subject rather than
describing it.

The first is that **gene symbols group pores with their accessory
subunits**. `KCNE*`, `KCNMB*`, `CACNB*`, `CACNG*`, `SCN*B` and the 25
`KCTD*` genes all sort inside the channel symbol space, and none is a pore;
the `CACNG` symbols are worse than uninformative, since CACNG2–8 are
AMPA-receptor auxiliary subunits with a claudin fold and no role in calcium
channels. A census that reads names counts the wrong things, and an unnamed
gene model in a newly assembled genome — the case a census exists for — has
no name to read.

The second is that **domain annotation is not a partition of the subject**.
Section 12 takes this up in detail with measurements; the short version is
that the canonical pore model is carried by families that are not channels,
and is absent from families that are.

## 2. The folds

A useful way to organise the subject is by the architecture of the pore
itself, because that is what is genuinely inherited. At least ten unrelated
folds solve the problem, and the boundaries between them are the boundaries
beyond which no alignment is meaningful.

**The P-loop (voltage-gated-like) fold.** Four subunits, or four repeats of
one chain, each contributing two transmembrane helices and a re-entrant loop
that dips into the membrane from the outside. The loop supplies the
selectivity filter; the last helix lines the inner cavity and carries the
gate. This is the fold of the potassium, sodium and calcium channels, the
TRP families, the cyclic-nucleotide-gated and hyperpolarisation-activated
channels, the intracellular calcium-release channels and the two-pore
channels [7,10]. It is also the only channel fold in this
review large enough to need a within-superfamily phylogeny of its own.

**The pentameric Cys-loop fold.** Five subunits around a central axis, each
with a β-sandwich extracellular ligand-binding domain over a four-helix
transmembrane bundle; the second helix of each subunit lines the pore. Known
from the *Torpedo* receptor [11], from the soluble acetylcholine
binding protein that turned out to be a complete ligand-binding domain
without a channel [12], and from bacterial members that lack the
disulfide the superfamily is named after [13,14].

**The iGluR fold.** Four subunits, a clamshell ligand-binding domain derived
from a bacterial periplasmic binding protein, and — pointing the other way —
a pore module that is an inverted P-loop [15]. The relationship
to potassium channels is structural and is confirmed by a cyanobacterial
glutamate receptor with an intact potassium filter [16].

**The trimeric folds.** P2X receptors are trimers with a large disulfide-rich
ectodomain and both termini inside [17]; the DEG/ENaC channels are
trimers of a different fold entirely [18,19].

**The CLC fold.** An antiparallel homodimer with two entirely independent
pores, one per subunit — inferred from noise analysis three decades before
the structure confirmed it [20].

**The TMEM16 fold.** Ten transmembrane helices per subunit with a
membrane-facing groove between TM4 and TM6; open to the lipid in the
scramblases, closed into a pore in the channels [21]. The OSCA
mechanosensitive channels and, on present evidence, the TMC proteins share
this fold without sharing detectable sequence [22].

**The remaining folds** are each a family of their own: the three-bladed
Piezo propeller with 38 transmembrane helices per subunit [23]; the
19-stranded β-barrel of VDAC [24]; the hexameric connexin
hemichannel [25] and the unrelated innexin/pannexin/LRRC8 architecture
[26,27]; the CALHM large pore [28]; MscL and MscS,
which are unrelated to each other [29,30]; the ABC fold of
CFTR [31]; the double-barrelled otopetrin fold [32]; the
voltage-sensor-only architecture of Hv1 [33,34]; and
TMEM175, whose potassium selectivity is achieved without a filter at all
[35].

### 2.1 What "unrelated" means here

Three of these relationships deserve care, because the word covers three
different situations.

*Undetectable but probably real.* TMEM16, OSCA and TMC share a fold and a
conduction pathway in the same place, and share no alignable sequence. The
common ancestry is likely; the data that would reconstruct it is not in the
sequences.

*Structurally equivalent, phylogenetically distant.* The iGluR pore is an
inverted potassium-channel pore. This is a real homology at the level of the
module, and it does not license putting an AMPA receptor and a Kir channel
in one alignment.

*Genuinely convergent, so far as anyone can tell.* Connexins and the
innexin/pannexin group build the same kind of large intercellular pore and
are not detectably related; vertebrates run both systems and invertebrates
only one [26,36]. "Not detectably related" is a statement
about detection methods, and section 11 treats it as a hypothesis rather
than a conclusion.

## 3. Selectivity

### 3.1 The potassium filter

The best-understood selectivity mechanism in biology is the potassium
filter: a short, absolutely conserved sequence — T-x-G-Y-G — whose backbone
carbonyls line a narrow pore and coordinate a dehydrated K⁺ ion in the same
geometry water would [7,37,38]. Sodium, which
is smaller, cannot be coordinated at that spacing and pays a dehydration
penalty it does not recover; selectivity is achieved by a structure that
fits the larger ion, which is the opposite of a sieve. The filter is
recognisable in *Streptomyces* KcsA, in *Drosophila* Shaker and in human
Kv1.1, and it is the one motif in this review that can be found in a raw
sequence with no reference and no annotation [39].

### 3.2 The four-repeat locus

The sodium and calcium channels arrived at selectivity differently. Both are
pseudo-tetramers — four homologous 6TM repeats in a single chain, the
product of two rounds of internal duplication from a Kv-like ancestor — and
their selectivity is set by one residue contributed by each repeat, at the
position corresponding to the potassium filter. In sodium channels those
four residues read D-E-K-A; in the high-voltage-activated calcium channels
E-E-E-E; and swapping them converts one into the other [40].

This locus is the only reliable way to tell the four-repeat families apart.
Their domain architectures are identical, their lengths overlap, and no
domain model distinguishes them. NALCN, a sodium leak channel, reads E-E-K-E
— a calcium-channel filter with a sodium-channel lysine in the third repeat
[41]. The T-type calcium channels read E-E-D-D rather than E-E-E-E,
which is a real difference within a family rather than an error.

Several of the families whose filters are least well characterised are also
the primary transducers of noxious stimuli, which is where much of their
pharmacological interest lies [42].

### 3.3 Selectivity that is not in the family

Two observations set a limit on how far family membership predicts what a
channel conducts.

The first is that charge selectivity in the Cys-loop receptors — the
difference between an excitatory nicotinic receptor and an inhibitory GABA
or glycine receptor — is set by a short ring of residues at the intracellular
end of the pore-lining helix, and can be inverted by three substitutions
[43]. Cation- and anion-selective members sit in the same superfamily
with the same fold.

The second is that AMPA-receptor calcium permeability is not encoded in the
genome at all: it is set by RNA editing at a single position in GluA2, and a
genomic census cannot see it [44].

The practical consequence is that selectivity is a property recorded per
family from experiment, never predicted per sequence. Any pipeline that
outputs a selectivity call from sequence alone is reporting family
membership with extra steps.

### 3.4 Potassium selectivity twice

TMEM175, the endolysosomal potassium channel, is selective for K⁺ and has no
T-x-G-Y-G filter and no structural relationship to the P-loop channels; its
conduction pathway is lined by isoleucines [35,45]. The
canonical filter is therefore one solution among at least two, and its near
universality reflects the dominance of one superfamily rather than a physical
necessity.

## 4. Gating

Gating is the property that most often gets mistaken for a clade. It is
worth stating plainly: **no gating mechanism in this review is confined to
one superfamily, and no superfamily is confined to one gating mechanism.**

**Voltage.** The canonical voltage sensor is a four-helix bundle whose
fourth helix carries a repeating pattern of arginines and moves outward on
depolarisation [46,47]. It is not exclusive to channels:
the same domain is coupled to a phosphatase in the voltage-sensing
phosphatases [48], and forms the entire conduction pathway of the
voltage-gated proton channel, which has no pore domain at all [33,34,49]. Conversely, several P-loop families retain a complete sensor and
are voltage-insensitive — the cyclic-nucleotide-gated channels among them
[50] — and the HCN channels use the same apparatus inverted, opening
on hyperpolarisation [51,52].

**Extracellular ligand.** Four unrelated superfamilies bind a transmitter
outside and open: the pentameric Cys-loop receptors, the tetrameric
glutamate receptors, the trimeric P2X receptors and, in a different sense,
the proton-gated ASICs [18]. Their ligand-binding domains are
unrelated clamshells, β-sandwiches and dolphin-shaped ectodomains
respectively.

**Intracellular ligand.** Cyclic nucleotides bind a conserved
cyclic-nucleotide-binding domain in CNG and HCN channels; the *same domain*,
in a form that binds nothing, sits in the EAG/ERG potassium channels
[50]. Calcium acts through an RCK ring in BK channels [53] and
through constitutively bound calmodulin in the SK channels; IP₃ and
ryanodine act on the huge cytosolic solenoids of the intracellular release
channels [54,55].

**Force.** Mechanosensitivity has been invented at least six times: Piezo
[56], the bacterial MscL and MscS — which are unrelated to each other
[30,57] — the plant and animal OSCA/TMEM63 family
[22,58], the TMC proteins of hair cells [59,60], the mechanosensitive two-pore-domain potassium channels
[61], and the *Drosophila* NOMPC channel with its 29-ankyrin
tether [62]. Some sense force from the lipid, some through a
tether; the mechanism is a physical strategy, not a lineage [63].

**Temperature, protons, lipids, volume.** Thermal gating in the TRP families
[64], proton gating in ASICs and Hv1, PIP₂ dependence across the
inward rectifiers and KCNQ channels [65], and cell-volume gating in
LRRC8 [66,67] complete the list without adding a clade.

### 4.1 The consequence for classification

None of this is recoverable from sequence family membership with any
reliability. HCN and CNG channels have the same domain architecture and
opposite voltage dependence. SK channels have an intact voltage sensor and
are voltage-insensitive. A classifier can say what family a protein belongs
to; it cannot say what opens it, and a review that blurs the two makes the
literature harder to use rather than easier.

## 5. The P-loop superfamily

Roughly 143 of the ~320 human pore-forming genes sit in this one
superfamily, and 79 of those are potassium channels [10,68]. It
is the only division of the subject large enough that its internal structure
is itself a research problem.

### 5.1 The potassium branch

Topology organises it [69]. The **6TM/1P** architecture — a voltage sensor plus a
pore module — covers the Shaker-related Kv1–Kv4 channels [5,6], the KCNQ channels, the EAG/ERG/ELK family, the
calcium-activated BK [53] and SK channels, and the electrically
silent modifier subunits that cannot conduct alone. The **2TM/1P** inward
rectifiers drop the sensor entirely and rectify through block by intracellular
magnesium and polyamines [65,70,71]. The **4TM/2P**
two-pore-domain channels are dimers carrying two pore loops per subunit, and
supply the leak conductances that set resting potential [61,72].

Prokaryotic members root all of it: KcsA [7], the calcium-gated
MthK [53], and KirBac [70] are not distant curiosities but the
structural basis of the whole branch [73].

### 5.2 The four-repeat channels

Two rounds of internal duplication produced a single chain with four 6TM
repeats: the voltage-gated sodium channels [4,74], the
calcium channels [75,76], the sodium leak channel NALCN
[41], and the sperm-specific CatSper complex, which is unusual in
assembling its pseudo-tetramer from four *different* genes [77,78]. The two-pore channels of endolysosomes, with two repeats
rather than one or four, are the intermediate the story predicts
[79].

### 5.3 The TRP families

TRPC, TRPV, TRPM, TRPA, TRPML and TRPP share the pore fold and little else:
their cytoplasmic domains range from a handful of ankyrin repeats to
fourteen, and TRPM2, TRPM6 and TRPM7 are chanzymes — a channel fused to an
enzyme [80,81]. The founding member came from a
*Drosophila* phototransduction mutant [82], and the division as a
whole is best understood as a set of cellular sensors rather than as a
functional class [83]; TRPV1, cloned as the
capsaicin receptor, became the first structure solved by single-particle
cryo-EM at near-atomic resolution and effectively started that era
[64,84]. TRPC2 is a pseudogene in humans and a functional
vomeronasal channel in mice, which makes it a useful test of whether a
census makes per-lineage presence calls correctly [85].

### 5.4 The intracellular calcium-release channels

The IP₃ and ryanodine receptors belong here structurally — their C-terminal
pore is a P-loop module — and nowhere near here in size: RyR1 is 5,038
residues, the largest ion channel known, against ~400 for a Kir channel
[54,55,86,87]. Almost all of the
difference is cytosolic solenoid with no counterpart elsewhere in the
superfamily, which is why any cross-family alignment has to be restricted to
the pore module.

## 6. The ligand-gated superfamilies

Four unrelated superfamilies open when a transmitter binds. They are worth
treating together only because their traps have the same shape: in each, the
ligand-binding domain has an independent existence outside channels, and a
search built on it returns proteins that do not conduct.

### 6.1 Pentameric (Cys-loop) receptors

Five subunits, each a β-sandwich ectodomain over four transmembrane helices,
with the second helix lining the pore [11,88]. The human
complement is large — nicotinic acetylcholine receptors, GABA_A receptors,
glycine receptors, the ionotropic serotonin receptor and the
zinc-activated channel — and the superfamily is alignable end to end, which
makes it the most tractable phylogenetic unit in this review.

Two members are informative precisely because they are not what the family
name implies. *AChBP*, secreted by snail glia, is a complete pentameric
ligand-binding domain with no transmembrane region at all, and its structure
made nicotinic pharmacology interpretable [12]. The bacterial
channels GLIC and ELIC lack the disulfide-bonded loop the superfamily is
named after [13,14] — the defining feature is not universal
within the group it defines. Both facts matter for classification: a rule
requiring the transmembrane domain alongside the ligand-binding domain
correctly rejects AChBP and, as measured in section 12, also rejects the
bacterial channels that root the tree.

Invertebrates add glutamate- and histamine-gated *chloride* channels absent
from vertebrates — the ivermectin target — and one of them provided the
first eukaryotic structure of the superfamily [89]. Human structures
followed for GABA_A and glycine receptors [90,91].

### 6.2 Ionotropic glutamate receptors

Four subunits, each with an amino-terminal clamshell, a second clamshell
that binds glutamate, and a pore module inserted the other way up
[15,92,93]. The AMPA, kainate, NMDA and
delta families differ in ligand requirements, kinetics and — for NMDA —
the coincidence detection that comes from needing glutamate, glycine and
depolarisation together [94,95].

The trap here is the clamshell. It is a bacterial periplasmic
binding-protein fold, and it is also the ligand-binding domain of every
class C G-protein-coupled receptor: the metabotropic glutamate receptors,
both GABA_B subunits, the calcium-sensing receptor and the sweet and umami
taste receptors. A search for glutamate receptors by clamshell returns a
receptor family that is not ionotropic at all, and the discriminating
evidence is the pore region.

The superfamily also extends well past animals with synapses. Plants carry
glutamate-receptor genes [96], and the insect ionotropic receptors are
divergent iGluRs that serve chemosensation and have mostly lost the
clamshell — so a clamshell-based search misses them entirely [97].
The delta receptors do not gate to glutamate at all, and are known to have a
functional pore only because a gain-of-function mutation holds it open
[98].

### 6.3 P2X receptors

Trimers gated by extracellular ATP, with two transmembrane helices per
subunit, both termini intracellular, and a large disulfide-rich ectodomain
shaped, in the original description, like a leaping dolphin [17,99]. Their distribution is patchy in an interesting way: present in
amoebozoa — where a *Dictyostelium* member works on intracellular vesicles
rather than the plasma membrane [100] — and absent from *Drosophila*
and *C. elegans*, which makes P2X a good test of whether a census can
distinguish gene loss from database absence.

### 6.4 DEG/ENaC channels

Trimers again, of an unrelated fold, and gated by almost everything except
voltage: protons in the acid-sensing channels [18], peptides in the
molluscan FaNaC, mechanical force in the *C. elegans* degenerins, and in the
epithelial sodium channel essentially nothing — ENaC is constitutively open
and regulated by proteolysis and trafficking [19,101].
The superfamily takes its name from gain-of-function alleles that kill the
neuron they are expressed in, found in the first genetic screens for touch
insensitivity [102]. Its size is strongly lineage-dependent: nine
members in humans against roughly thirty in *C. elegans* and thirty
pickpocket genes in *Drosophila*, so a human-only count badly
misrepresents it.

## 7. Anion channels

Anion channels were the last major group to yield their molecular identities,
and the delay shows in how many of them are still argued about. Chloride is
the abundant permeant anion, its equilibrium potential sits close to rest in
many cells, and small conductances are hard to attribute — so several
activities were characterised physiologically for decades before the gene
was found [103].

### 7.1 CLC: one family, two mechanisms

The CLC family is a homodimer in which each subunit contains a complete,
independent pore — a "double-barrelled" architecture inferred from noise
analysis long before the structure confirmed it [20]. What the
structure did not predict is that half the family are not channels at all.
The bacterial prototype ClC-ec1 is a 2Cl⁻/H⁺ antiporter [8], and so
are the mammalian endosomal and lysosomal members CLC-3 to CLC-7; only
CLC-1, CLC-2 and the CLC-K channels conduct passively [104].

The mechanistic difference turns on a conserved glutamate on the
intracellular side, and the transport phenotype is ancestral: channel
behaviour is derived within a family whose prototype transports. For a
census this is the cleanest possible demonstration that family membership
and mechanism are separate questions with separate evidence, and that
answering one does not answer the other.

### 7.2 CFTR: a transporter that is a channel

CFTR is an ABC transporter — two transmembrane domains, two
nucleotide-binding domains, the whole architecture — that conducts anions
down their gradient [31,105,106]. It is gated by ATP
binding at the nucleotide-binding domain interface and licensed by PKA
phosphorylation of a regulatory domain its siblings lack. Its nearest
relatives, including the sulfonylurea receptors that regulate the Kir6
potassium channels, are not channels.

CFTR is therefore the counter-example to every rule of the form "fold X
implies transporter". Any classifier that excludes proteins for carrying an
ABC domain excludes the most clinically important chloride channel in the
human genome; the discriminating evidence is one accession.

### 7.3 TMEM16: channels and scramblases on one fold

Three groups identified TMEM16A simultaneously as the long-sought
calcium-activated chloride channel [107,108,109].
The structure of a fungal relative then showed something unexpected: the same
fold, with the conduction pathway open to the lipid rather than closed into
a pore, is a phospholipid scramblase [21]. Most of the family
scrambles [9]; TMEM16A and TMEM16B conduct [110].

The two functions are not distinguishable by domain architecture, by length
or by any sequence feature yet identified, which makes the channel /
scramblase split a literature assignment rather than a measurement — and a
hypothesis a phylogeny can test.

### 7.4 Bestrophins, tweety, and the volume-regulated channel

Bestrophins are pentamers with a single central pore, unrelated to anything
else here, identified through a macular dystrophy gene [111,112].

The volume-regulated anion channel — a current known since the 1980s and
central to cell-volume regulation — was assigned to two different protein
families before two parallel genome-wide screens identified LRRC8 in 2014
[27,66,67]. The tweety proteins, one of the earlier
candidates, remain in catalogues as contested: recent structures show no
obvious conduction pathway. Keeping them, flagged, is the honest option;
dropping a family because the literature moved makes a census a record of
opinion.

## 8. Large pores

A pore wide enough to pass ATP is a different object from a filter that
discriminates K⁺ from Na⁺ by a fraction of an ångström, and the two are
easily conflated by a definition that only asks whether ions cross. The
large-pore channels are included here because they are gated aqueous
pathways, and flagged because no statement about ionic selectivity should
include them.

**Connexins** form hexameric hemichannels that dock across an extracellular
gap into a complete intercellular channel [25,113]. There
are 21 human genes, named twice over — by Greek-letter subfamily (`GJA1`)
and by predicted molecular weight (Cx43) — with the two systems not in
register, which is a pure nomenclature hazard. `GJB2` mutations are the
commonest cause of inherited deafness.

**Innexins** do the same job in invertebrates and are unrelated to connexins
[36]. **Pannexins** are the vertebrate members of the innexin family;
they are glycosylated, so they do not dock into junctions, and work instead
as regulated ATP-release conduits [26].

**LRRC8** completes the group: its transmembrane region is homologous to
pannexin — Pfam names the domain "Pannexin-like TM region of LRRC8" — with a
leucine-rich-repeat domain added [27]. Substrate selectivity,
including cisplatin uptake, is set by which of LRRC8B–E join the obligatory
LRRC8A.

**CALHM** is a fourth, unrelated solution: a large pore of variable
oligomeric state that releases ATP from taste cells, and is the transmitter
release mechanism for sweet, bitter and umami taste [28,114].

### 8.1 A convergence with a twist

Vertebrates run two independent intercellular-channel systems. One —
innexin/pannexin, with LRRC8 attached — they share with invertebrates. The
other, the connexins, they appear to have invented; no connexin has been
found outside chordates, and no sequence relationship to the innexin group
is detectable [26,36].

This is the best-supported case of fold-level convergence in the subject.
It is also, precisely because it rests on a negative, the one most worth
re-testing with profile-profile and structural methods rather than
sequence search: "not detectably related" is a statement about the method
used to look.

## 9. Mechanosensation

Turning force into a conductance is the clearest case in this review of one
problem solved repeatedly by unrelated proteins. At least six architectures
do it, and they share nothing but the physics [63].

**Piezo1 and Piezo2** are three-bladed propellers of 38 transmembrane
helices per subunit — the largest transmembrane count of any known channel —
that bow the surrounding membrane into a dome and flatten under tension
[23,56]. PIEZO2 carries touch and proprioception; PIEZO1
loss-of-function causes lymphatic dysplasia and gain-of-function a
haemolytic anaemia.

**MscL and MscS** are the bacterial safety valves, and they are unrelated to
each other as well as to everything else. MscL is a 136-residue pentamer
with the largest conductance known, opening under near-lytic tension to
release osmolytes [29,57]; MscS is a heptamer with a
lipid-facing gate and a large cytoplasmic vestibule [30]. The family
is present in bacteria, archaea, plants, fungi and protists and absent from
animals — one of the sharpest presence/absence claims in the subject, and
one that deserves testing rather than repeating.

**OSCA/TMEM63** were found in plants as hyperosmolality-gated calcium
channels [58] and only then recognised in animals as an evolutionarily
conserved family of mechanically activated channels [22]. They share
the TMEM16 fold without sharing its sequence — a relationship visible only
in structure. Their history is a caution: a family that had been present in
animal genomes all along was found in plants first, and a human-only survey
would have called it absent.

**TMC1 and TMC2** are the strongest candidates for the pore of the hair-cell
mechanotransduction channel, the fastest sensory transduction in vertebrates
[59,60]. The evidence is genetic, physiological and
structural-by-homology; there is still no structure of a conducting state,
which is why careful catalogues list them as contested. TMC6 and TMC8 do
something else entirely in keratinocytes.

**The mechanosensitive potassium channels** — TRAAK and TREK, in the
two-pore-domain family — are gated by membrane tension through a lipid-facing
opening in an otherwise ordinary P-loop architecture [61]. Force
gating here is a property of one branch of a superfamily whose other members
are gated by everything else.

**NOMPC**, a *Drosophila* TRP channel with 29 ankyrin repeats forming a
spring between the channel and the cytoskeleton, is the clearest tethered
mechanotransducer known [62]. It was lost in mammals.

### 9.1 Why this matters for a census

Mechanosensitivity is a mechanism, not a lineage. A catalogue organised by
what opens a channel would place Piezo, MscS, OSCA, TMC, TRAAK and NOMPC in
one group and would be wrong about every evolutionary question it was then
asked. A catalogue organised by fold keeps them apart and records the shared
physiology as an annotation — which is the only arrangement in which the
question "how many times has this been invented?" can even be posed.

## 10. Channels of intracellular membranes

Organelle channels are systematically under-represented in channel censuses,
for a reason worth naming: they were found by physiology on isolated
membranes rather than by patch clamp on cells, so they entered the
literature later, with less standardised nomenclature, and several still
have contested molecular identities. A census whose assumptions are shaped
by plasma-membrane channels — a signal peptide, an extracellular
ligand-binding domain, a voltage sensor — finds fewer of them than exist.

**IP₃ receptors and ryanodine receptors** release calcium from the
endoplasmic and sarcoplasmic reticulum. Both are homotetramers of enormous
subunits — 2,700 residues for an IP₃ receptor, 5,038 for RyR1, the largest
ion channel known — in which a P-loop pore module sits atop a cytosolic
solenoid that acts as a signal-integrating scaffold [54,55,86,87,115]. They share four Pfam
domains and are the textbook case of two families that domain annotation
cannot separate.

**TRIC channels** (TMEM38A/B) are trimeric monovalent cation channels of the
SR/ER membrane whose entire function is to provide the counter-ion flux that
lets the release channels work without building an opposing voltage
[116]. A channel whose job is to make another channel possible is
easy to leave out of a list of channels.

**The mitochondrial calcium uniporter** carried a measured current for fifty
years before its gene was found, by two groups simultaneously using
comparative genomics [117,118]. It is the most
calcium-selective channel known. Its regulators — MICU1–3, EMRE, MCUB — are
a good test of the pore/accessory distinction: MCUB is a pore subunit, EMRE
is required for conduction and lines no pore.

**VDAC** is the main conduit across the outer mitochondrial membrane and the
only all-β channel in the subject: a 19-stranded barrel, an odd number that
forces one parallel β-pair, unlike any bacterial porin [24]. It
is gated by voltage and regulated, which is what keeps it inside the
definition while ungated bacterial porins fall outside.

**TMEM175** is the endolysosomal potassium channel, and it is selective for
K⁺ without a T-x-G-Y-G filter — its pathway is lined by isoleucines
[35,45]. It is also a Parkinson's disease risk locus.
Together with the canonical filter it establishes that potassium selectivity
has evolved at least twice, which is the single most useful fact in this
review for calibrating how much weight a motif can carry.

**ORAI channels** are not intracellular, but belong beside them: the
store-operated calcium current is gated by physical contact with STIM
proteins in the ER membrane across a junction, so the gate for a
plasma-membrane channel sits in another organelle [119,120,121,122].

**Two-pore channels** (TPC1/2) sit in endolysosomal membranes, are gated by
PI(3,5)P₂ and NAADP, and have two 6TM repeats rather than one or four —
exactly the intermediate the duplication history of the four-repeat channels
predicts [79].

## 11. Evolution

### 11.1 Channels are older than the things they are famous for

Voltage-gated channels are usually introduced through the action potential,
which invites the assumption that they arose with nervous systems. They did
not. Potassium channels are present across bacteria and archaea, and the
prokaryotic structures are not distant analogues but the direct structural
basis of the eukaryotic families [73,123]. Sodium-channel
homologues predate the origin of nervous systems in animals
[124], and the voltage-gated set was already diversifying at the
emergence of the Metazoa [125]. Sponges and placozoans, which have no
neurons, carry substantial channel repertoires; the anemone *Nematostella*
carries a nearly complete one.

The ordering is therefore the reverse of the intuitive one. The channels
came first and the nervous system was assembled from parts that already
existed, which is why a census restricted to animals with synapses answers
the wrong question.

### 11.2 Repertoires are lineage-specific in both directions

Comparative work on complete genomes shows expansions and losses that
correlate with lifestyle rather than with taxonomic rank [68,126,127]. *C. elegans* has roughly ninety potassium
channel genes; humans have about seventy-nine. The first genome-scale survey
of a channel repertoire, in *Drosophila*, already made the point that the
set is not a scaled-down version of the vertebrate one [128]. The DEG/ENaC superfamily has
nine human members against roughly thirty in *C. elegans* and thirty
pickpocket genes in *Drosophila*. P2X receptors are present in amoebozoa and
absent from both major invertebrate models [100]. TRPC2 is a
functional vomeronasal channel in mice and a pseudogene in humans
[85]. NOMPC, the clearest tethered mechanotransducer known, was lost
in mammals [62].

Plants make the point from the other side. They have no voltage-gated
sodium channels and no nervous system, and they have glutamate receptors
[96], two-pore channels, mechanosensitive MscS-like channels and the
OSCA family that was found in plants before it was recognised in animals
[22,58,129].

### 11.3 Convergence, and the limits of the claim

Four convergences are well enough supported to state, with the caveat that
each rests on failure to detect homology:

- **Potassium selectivity**, achieved by the canonical filter and,
  independently, by TMEM175's isoleucine-lined pathway [35].
- **Intercellular channels**, built by connexins in chordates and by the
  unrelated innexin/pannexin family elsewhere [26,36].
- **Mechanosensitivity**, invented at least six times (section 9).
- **Large-pore ATP release**, by pannexin, LRRC8 and CALHM on at least two
  unrelated folds [27,114].

The honest formulation is that these are cases where sequence methods find
no relationship. Structure comparison and profile-profile search can push
the detection limit further, and the TMEM16/OSCA/TMC clan — three families
assigned separate domain models that turned out to share a fold — is the
warning that they sometimes do.

### 11.4 What a phylogeny of ion channels can and cannot be

There is no alignment containing a nicotinic receptor and a Kv channel, so
there is no tree of ion channels. What is well defined is a **forest**:
within-family trees, which are conventional and where most biologically
useful results live; within-superfamily trees, which for the P-loop
superfamily must be built on the pore module alone, because members range
from 400 to 5,038 residues and almost all of the difference is cytosolic
machinery with no counterpart across families.

Between superfamilies the correct object is a **network** of structural
similarity, with no branch lengths and no support values, in which an edge
records measured fold similarity and asserts no common ancestor. Every
published tree spanning non-homologous channel superfamilies is measuring
alignment artefacts.

## 12. What the databases record

Everything above is drawn from the primary literature. This section is
different: it reports what the annotation databases actually contain about
these proteins, measured directly against InterPro and UniProt on
2026-08-19 with the queries recorded in `results/s0_baseline/`
[130,131,132]. It is included because every
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
  [33,34].
- **TPTE**, a voltage-sensing phosphatase, carries it and is not a channel
  [48].

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

### 12.5 Signatures shared with non-channels

Of 26 domain accessions carried by more than one family in this catalogue,
**15 cross the channel / non-channel boundary**. Beyond `PF00520`:

- `PF02214`, the Kv T1 tetramerisation domain, is a generic BTB/POZ domain
  carried by 25 human KCTD proteins that are mostly ubiquitin-ligase
  adaptors.
- `PF01094`, the iGluR clamshell, is the ligand-binding domain of every
  class C GPCR (§6.2).
- `PF02931`, the Cys-loop ligand-binding domain, is a complete soluble
  protein in AChBP [12].
- `PF02815`, the MIR domain of the intracellular release channels, is
  carried by the protein O-mannosyltransferases.

### 12.6 The outgroups are the worst-covered proteins in the subject

The finding with the most consequence for phylogenetics: **three of the four
prokaryotic outgroups that root the major superfamily trees cannot be found
by the domain rules that define their own superfamilies.**

- **GLIC and ELIC**, which root the Cys-loop tree, carry `PF02931` and no
  `PF02932` — so a rule requiring the transmembrane domain alongside the
  ligand-binding domain, which correctly rejects AChBP, also rejects them
  [13,14].
- **GluR0**, the cyanobacterial glutamate receptor that roots the iGluR
  tree, is annotated as `PF07885` + `PF00497` — *a potassium-channel pore
  fused to a bacterial solute-binding protein* — with none of the eukaryotic
  iGluR models [16]. That is the inverted-P-loop story visible
  directly in the annotation, and it means no iGluR-signature search will
  ever find it.
- **NavAb**, the prototype prokaryotic sodium channel, carries neither
  `PF06512` nor `PF11933`, the two Nav-specific models: it is the pore and
  nothing else [74].

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

## 13. Disease and pharmacology

### 13.1 Channelopathies

The ion channels carry the largest set of Mendelian diseases attached to any
protein family, and the pattern is instructive: because a channel's output
is an electrical or chemical flux, small changes in gating kinetics produce
large, tissue-specific phenotypes, and both loss and gain of function are
pathogenic in ways that are often opposite [133,134].

The list spans nearly every family in this review. `CFTR` in cystic fibrosis
[31]; `SCN1A` in Dravet syndrome and `SCN5A` in Brugada and long-QT
syndromes; `KCNQ1` and `KCNH2` in the commonest inherited arrhythmias
[135]; `RYR1` in malignant hyperthermia and central core disease;
`CLCN1` in myotonia; `GJB2` as the commonest cause of inherited deafness;
`PKD1` and `PKD2` in autosomal dominant polycystic kidney disease; `KCNJ11`
in neonatal diabetes; `SCNN1A` in Liddle syndrome; `BEST1` in Best
vitelliform macular dystrophy [111]; `TMC1` in inherited deafness
[59]; `ORAI1` in a severe combined immunodeficiency [119];
`PIEZO1` in lymphatic dysplasia and hereditary xerocytosis [56];
`MCOLN1` in mucolipidosis IV; `TMEM38B` in a recessive osteogenesis
imperfecta [116]; and `TMEM175` as a Parkinson's disease risk locus
[35].

Two features of this list matter for a computational census. First, its
breadth is a reason to get family assignment right: a variant's
interpretation depends on which channel it is in. Second, the diseases
cluster in the *best-studied* families, which is a statement about
ascertainment as much as about biology — and any analysis correlating
disease burden with family properties has to say so.

### 13.2 Drug targets

Ion channels are among the most heavily drugged protein classes, second only
to GPCRs by most counts [136,137]. Local anaesthetics,
class I antiarrhythmics and many anticonvulsants act on voltage-gated sodium
channels; dihydropyridines and related agents on L-type calcium channels;
sulfonylureas on the Kir6/SUR complex; benzodiazepines, barbiturates and
most general anaesthetics on GABA_A receptors; ivermectin on invertebrate
glutamate-gated chloride channels; amantadine on the influenza M2 viroporin;
gabapentin and pregabalin on the calcium-channel α2δ subunit — which is an
auxiliary subunit and not a pore [138,139].

The last example is worth dwelling on. One of the most prescribed
"ion-channel drugs" targets a protein that this review's definition excludes
from the channel census. That is not an argument for changing the
definition; it is an argument for keeping the accessory subunits in the
catalogue, marked, so that the question can be asked at all.

Systematic target annotation is maintained by IUPHAR/BPS, whose concise
guide is the reference nomenclature for the field [140]. Where its
family boundaries differ from the ones used here, the difference is almost
always one of the six scope decisions in section 1.

## 14. Open questions

The questions below are the ones this review could not answer from the
literature, phrased so that a computational census could.

**Q1 — How many ion-channel genes are there, and how much of the published
range is scope?** Human counts run 240–400. Almost all of the spread is the
six decisions in section 1, not disagreement about biology, but nobody has
recomputed the total under each. A census that reports its answer once per
scope decision would replace an argument with a table.

**Q2 — Can family assignment be made without reading gene names, and where
does it fail?** Nomenclature encodes the order of discovery, and unnamed
gene models in new genomes have none. A classifier built only on positive
sequence and structural evidence would show which families are recoverable
and which depend on the annotation already being right.

**Q3 — Which superfamilies are reachable by domain search at all?** Section
12 shows that most TRP families and three of the four superfamily outgroups
are not. The recall of domain search, profile search and genomic search,
measured per superfamily, would set an upper bound on every existing census.

**Q4 — Is the absence of MscS-family channels from animals a genome fact or
a database fact?** It is repeated confidently and rests on annotation
coverage in exactly the lineages where coverage is weakest.

**Q5 — Do selectivity filters track phylogeny?** The T-type calcium
channels read E-E-D-D where their relatives read E-E-E-E. Mapping the filter
onto a pore-module tree would show whether selectivity follows descent or
has been re-tuned repeatedly.

**Q6 — Is connexin/pannexin convergence real, or homology below the
detection limit?** The claim rests on a negative from sequence search.
Profile-profile comparison and structure alignment can move that limit, and
the TMEM16/OSCA/TMC clan shows they sometimes do.

**Q7 — Do the mechanism-based splits survive a phylogeny?** CLC channels
against CLC transporters, TMEM16 channels against TMEM16 scramblases. If
each mechanism forms a clade, the split is evolutionary; if it cuts across
clades, the split is physiology imposed on a family that does not respect
it.

**Q8 — How much apparent lineage-specific loss is assembly and gene-set
quality?** Every claim in section 11.2 rests on absence in a particular
gene set, and gene sets differ in provenance and completeness far more than
the claims usually acknowledge.

## 15. How the claims in this review were checked

This review is generated, and its verification is mechanical rather than
editorial. Three checks run on every build.

**Every reference was resolved against a live literature service.** The
source list (`scripts/review_sources.py`) contains *titles*, not citations.
`scripts/s0_review_refs.py` queries Europe PMC for each, accepts a record
only when its title is at least 90 % similar to the one requested, and
writes the returned PMID, DOI, year, journal and authors to
`results/s0_baseline/references.tsv`. A title that does not resolve is
reported and not written.

That check has teeth. On the first run, seven titles failed: five were short
generic titles where the search ranked another paper, and two were
misremembered — the Cav1.1 structure paper is *Structure of the voltage-gated
calcium channel Cav1.1 complex* (Science, 2015), and the Moran review is
*Evolution of voltage-gated ion channels at the emergence of Metazoa*, not
"of Nervous Systems". Both were corrected against the returned records. The
final list is 148 of 148 resolved.

**A citation with no verified reference is a build failure.**
`scripts/s0_review_build.py` collects every `[key]` in the section files,
renumbers them into order of first appearance, and refuses to write the
document if any key is missing from the reference table. Prose therefore
cannot cite a paper that was never resolved.

**Every `[db]` statement is re-derivable.** Section 12 reports measurements
rather than literature, and the queries that produced them are in
`scripts/s0_catalogue_verify.py`; the tables are in
`results/s0_baseline/`. Re-running the script re-derives every number in
that section, and the S0 report is rendered from the tables rather than
written alongside them.

### The computational toolchain

Where this review reports its own measurements, the tools are: MAFFT for
alignment [141], trimAl for alignment trimming
[142], IQ-TREE 2 with ModelFinder and ultrafast bootstrap
for phylogenies [143,144], HMMER for profile searches
[145], Foldseek for structure comparison [146], and the
AlphaFold database for predicted structures where experimental ones are
absent [147,148]. Domain and sequence annotation comes from
InterPro, Pfam and UniProt [130,131,132].

### What this does not verify

Resolution confirms that a paper exists with the title cited. It does not
confirm that the paper says what the sentence citing it claims. That
remains an editorial responsibility, and the mechanical checks should not be
mistaken for having discharged it.

Nor does it verify the *selection*: which papers are cited, and which
landmark results are missing, is a judgement, and a review assembled by
searching for known titles will under-represent work the assembler did not
already know about. The bibliography here is weighted towards structures,
identifications and reviews in high-visibility journals, and that weighting
is a property of how it was built.

## References

1. HODGKIN AL, HUXLEY AF. A quantitative description of membrane current and its application to conduction and excitation in nerve. *J Physiol 1952*. doi:10.1113/jphysiol.1952.sp004764

2. Neher E, Sakmann B. Single-channel currents recorded from membrane of denervated frog muscle fibres. *Nature 1976*. doi:10.1038/260799a0

3. Noda M et al. Primary structure of alpha-subunit precursor of Torpedo californica acetylcholine receptor deduced from cDNA sequence. *Nature 1982*. doi:10.1038/299793a0

4. Noda M et al. Primary structure of Electrophorus electricus sodium channel deduced from cDNA sequence. *Nature 1984*. doi:10.1038/312121a0

5. Papazian DM et al. Cloning of genomic and complementary DNA from Shaker, a putative potassium channel gene from Drosophila. *Science 1987*. doi:10.1126/science.2441470

6. Tempel BL et al. Sequence of a probable potassium channel component encoded at Shaker locus of Drosophila. *Science 1987*. doi:10.1126/science.2441471

7. Doyle DA et al. The structure of the potassium channel: molecular basis of K+ conduction and selectivity. *Science 1998*. doi:10.1126/science.280.5360.69

8. Accardi A, Miller C. Secondary active transport mediated by a prokaryotic homologue of ClC Cl- channels. *Nature 2004*. doi:10.1038/nature02314

9. Suzuki J et al. Calcium-dependent phospholipid scrambling by TMEM16F. *Nature 2010*. doi:10.1038/nature09583

10. Yu FH, Catterall WA. The VGL-chanome: a protein superfamily specialized for electrical signaling and ionic homeostasis. *Sci STKE 2004*. doi:10.1126/stke.2532004re15

11. Unwin N. Refined structure of the nicotinic acetylcholine receptor at 4A resolution. *J Mol Biol 2005*. doi:10.1016/j.jmb.2004.12.031

12. Brejc K et al. Crystal structure of an ACh-binding protein reveals the ligand-binding domain of nicotinic receptors. *Nature 2001*. doi:10.1038/35077011

13. Hilf RJ, Dutzler R. X-ray structure of a prokaryotic pentameric ligand-gated ion channel. *Nature 2008*. doi:10.1038/nature06717

14. Bocquet N et al. X-ray structure of a pentameric ligand-gated ion channel in an apparently open conformation. *Nature 2009*. doi:10.1038/nature07462

15. Sobolevsky AI, Rosconi MP, Gouaux E. X-ray structure, symmetry and mechanism of an AMPA-subtype glutamate receptor. *Nature 2009*. doi:10.1038/nature08624

16. Chen GQ et al. Functional characterization of a potassium-selective prokaryotic glutamate receptor. *Nature 1999*. doi:10.1038/45568

17. Kawate T et al. Crystal structure of the ATP-gated P2X(4) ion channel in the closed state. *Nature 2009*. doi:10.1038/nature08198

18. Jasti J et al. Structure of acid-sensing ion channel 1 at 1.9 A resolution and low pH. *Nature 2007*. doi:10.1038/nature06163

19. Canessa CM et al. Amiloride-sensitive epithelial Na+ channel is made of three homologous subunits. *Nature 1994*. doi:10.1038/367463a0

20. Dutzler R et al. X-ray structure of a ClC chloride channel at 3.0 A reveals the molecular basis of anion selectivity. *Nature 2002*. doi:10.1038/415287a

21. Brunner JD et al. X-ray structure of a calcium-activated TMEM16 lipid scramblase. *Nature 2014*. doi:10.1038/nature13984

22. Murthy SE et al. OSCA/TMEM63 are an Evolutionarily Conserved Family of Mechanically Activated Ion Channels. *Elife 2018*. doi:10.7554/elife.41844

23. Ge J et al. Architecture of the mammalian mechanosensitive Piezo1 channel. *Nature 2015*. doi:10.1038/nature15247

24. Colombini M. A candidate for the permeability pathway of the outer mitochondrial membrane. *Nature 1979*. doi:10.1038/279643a0

25. Maeda S et al. Structure of the connexin 26 gap junction channel at 3.5 A resolution. *Nature 2009*. doi:10.1038/nature07869

26. Panchin Y et al. A ubiquitous family of putative gap junction molecules. *Curr Biol 2000*. doi:10.1016/s0960-9822(00)00576-5

27. Deneka D et al. Structure of a volume-regulated anion channel of the LRRC8 family. *Nature 2018*. doi:10.1038/s41586-018-0134-y

28. Ma Z et al. Calcium homeostasis modulator 1 (CALHM1) is the pore-forming subunit of an ion channel that mediates extracellular Ca2+ regulation of neuronal excitability. *Proc Natl Acad Sci U S A 2012*. doi:10.1073/pnas.1204023109

29. Chang G et al. Structure of the MscL homolog from Mycobacterium tuberculosis: a gated mechanosensitive ion channel. *Science 1998*. doi:10.1126/science.282.5397.2220

30. Bass RB et al. Crystal structure of Escherichia coli MscS, a voltage-modulated and mechanosensitive channel. *Science 2002*. doi:10.1126/science.1077945

31. Riordan JR et al. Identification of the cystic fibrosis gene: cloning and characterization of complementary DNA. *Science 1989*. doi:10.1126/science.2475911

32. Tu YH et al. An evolutionarily conserved gene family encodes proton-selective ion channels. *Science 2018*. doi:10.1126/science.aao3264

33. Ramsey IS et al. A voltage-gated proton-selective channel lacking the pore domain. *Nature 2006*. doi:10.1038/nature04700

34. Sasaki M, Takagi M, Okamura Y. A voltage sensor-domain protein is a voltage-gated proton channel. *Science 2006*. doi:10.1126/science.1122352

35. Cang C et al. TMEM175 Is an Organelle K(+) Channel Regulating Lysosomal Function. *Cell 2015*. doi:10.1016/j.cell.2015.08.002

36. Phelan P et al. Innexins: a family of invertebrate gap-junction proteins. *Trends Genet 1998*. doi:10.1016/s0168-9525(98)01547-9

37. Heginbotham L et al. Mutations in the K+ channel signature sequence. *Biophys J 1994*. doi:10.1016/s0006-3495(94)80887-2

38. Zhou Y et al. Chemistry of ion coordination and hydration revealed by a K+ channel-Fab complex at 2.0 A resolution. *Nature 2001*. doi:10.1038/35102009

39. MacKinnon R. Potassium channels and the atomic basis of selective ion conduction (Nobel Lecture). *Angew Chem Int Ed Engl 2004*. doi:10.1002/anie.200400662

40. Heinemann SH et al. Calcium channel characteristics conferred on the sodium channel by single mutations. *Nature 1992*. doi:10.1038/356441a0

41. Lu B et al. The neuronal channel NALCN contributes resting sodium permeability and is required for normal respiratory rhythm. *Cell 2007*. doi:10.1016/j.cell.2007.02.041

42. Giniatullin R. Ion Channels of Nociception. *Int J Mol Sci 2020*. doi:10.3390/ijms21103553

43. Galzi JL et al. Mutations in the channel domain of a neuronal nicotinic receptor convert ion selectivity from cationic to anionic. *Nature 1992*. doi:10.1038/359500a0

44. Sommer B et al. RNA editing in brain controls a determinant of ion flow in glutamate-gated channels. *Cell 1991*. doi:10.1016/0092-8674(91)90568-j

45. Lee C et al. The lysosomal potassium channel TMEM175 adopts a novel tetrameric architecture. *Nature 2017*. doi:10.1038/nature23269

46. Long SB, Campbell EB, Mackinnon R. Crystal structure of a mammalian voltage-dependent Shaker family K+ channel. *Science 2005*. doi:10.1126/science.1116269

47. Catterall WA. From ionic currents to molecular mechanisms: the structure and function of voltage-gated sodium channels. *Neuron 2000*. doi:10.1016/s0896-6273(00)81133-2

48. Murata Y et al. Phosphoinositide phosphatase activity coupled to an intrinsic voltage sensor. *Nature 2005*. doi:10.1038/nature03650

49. DeCoursey TE. Voltage-gated proton channels: molecular biology, physiology, and pathophysiology of the H(V) family. *Physiol Rev 2013*. doi:10.1152/physrev.00011.2012

50. Kaupp UB, Seifert R. Cyclic nucleotide-gated ion channels. *Physiol Rev 2002*. doi:10.1152/physrev.00008.2002

51. Robinson RB, Siegelbaum SA. Hyperpolarization-activated cation currents: from molecules to physiological function. *Annu Rev Physiol 2003*. doi:10.1146/annurev.physiol.65.092101.142734

52. Lee CH, MacKinnon R. Structures of the Human HCN1 Hyperpolarization-Activated Channel. *Cell 2017*. doi:10.1016/j.cell.2016.12.023

53. Jiang Y et al. Crystal structure and mechanism of a calcium-gated potassium channel. *Nature 2002*. doi:10.1038/417515a

54. Furuichi T et al. Primary structure and functional expression of the inositol 1,4,5-trisphosphate-binding protein P400. *Nature 1989*. doi:10.1038/342032a0

55. Takeshima H et al. Primary structure and expression from complementary DNA of skeletal muscle ryanodine receptor. *Nature 1989*. doi:10.1038/339439a0

56. Coste B et al. Piezo1 and Piezo2 are essential components of distinct mechanically activated cation channels. *Science 2010*. doi:10.1126/science.1193270

57. Sukharev SI et al. A large-conductance mechanosensitive channel in E. coli encoded by mscL alone. *Nature 1994*. doi:10.1038/368265a0

58. Yuan F et al. OSCA1 mediates osmotic-stress-evoked Ca2+ increases vital for osmosensing in Arabidopsis. *Nature 2014*. doi:10.1038/nature13593

59. Kawashima Y et al. Mechanotransduction in mouse inner ear hair cells requires transmembrane channel-like genes. *J Clin Invest 2011*. doi:10.1172/jci60405

60. Pan B et al. TMC1 Forms the Pore of Mechanosensory Transduction Channels in Vertebrate Inner Ear Hair Cells. *Neuron 2018*. doi:10.1016/j.neuron.2018.07.033

61. Brohawn SG, del Mármol J, MacKinnon R. Crystal structure of the human K2P TRAAK, a lipid- and mechano-sensitive K+ ion channel. *Science 2012*. doi:10.1126/science.1213808

62. Walker RG, Willingham AT, Zuker CS. A Drosophila mechanosensory transduction channel. *Science 2000*. doi:10.1126/science.287.5461.2229

63. Ranade SS, Syeda R, Patapoutian A. Mechanically Activated Ion Channels. *Neuron 2015*. doi:10.1016/j.neuron.2015.08.032

64. Caterina MJ et al. The capsaicin receptor: a heat-activated ion channel in the pain pathway. *Nature 1997*. doi:10.1038/39807

65. Hibino H et al. Inwardly rectifying potassium channels: their structure, function, and physiological roles. *Physiol Rev 2010*. doi:10.1152/physrev.00021.2009

66. Voss FK et al. Identification of LRRC8 heteromers as an essential component of the volume-regulated anion channel VRAC. *Science 2014*. doi:10.1126/science.1252826

67. Qiu Z et al. SWELL1, a plasma membrane protein, is an essential component of volume-regulated anion channel. *Cell 2014*. doi:10.1016/j.cell.2014.03.024

68. Jegla TJ et al. Evolution of the human ion channel set. *Comb Chem High Throughput Screen 2009*. doi:10.2174/138620709787047957

69. Yellen G. The voltage-gated potassium channels and their relatives. *Nature 2002*. doi:10.1038/nature00978

70. Kuo A et al. Crystal structure of the potassium channel KirBac1.1 in the closed state. *Science 2003*. doi:10.1126/science.1085028

71. Whorton MR, MacKinnon R. Crystal structure of the mammalian GIRK2 K+ channel and gating regulation by G proteins, PIP2, and sodium. *Cell 2011*. doi:10.1016/j.cell.2011.07.046

72. Enyedi P, Czirják G. Molecular background of leak K+ currents: two-pore domain potassium channels. *Physiol Rev 2010*. doi:10.1152/physrev.00029.2009

73. Derst C, Karschin A. Evolutionary link between prokaryotic and eukaryotic K+ channels. *J Exp Biol 1998*. PMID:9866872

74. Payandeh J et al. The crystal structure of a voltage-gated sodium channel. *Nature 2011*. doi:10.1038/nature10238

75. Wu J et al. Structure of the voltage-gated calcium channel Cav1.1 complex. *Science 2015*. doi:10.1126/science.aad2395

76. Catterall WA. Voltage-gated calcium channels. *Cold Spring Harb Perspect Biol 2011*. doi:10.1101/cshperspect.a003947

77. Jarow JP. A sperm ion channel required for sperm motility and male fertility. *J Urol 2002*. doi:10.1097/00005392-200211000-00119

78. Kirichok Y, Navarro B, Clapham DE. Whole-cell patch-clamp measurements of spermatozoa reveal an alkaline-activated Ca2+ channel. *Nature 2006*. doi:10.1038/nature04417

79. Calcraft PJ et al. NAADP mobilizes calcium from acidic organelles through two-pore channels. *Nature 2009*. doi:10.1038/nature08030

80. Venkatachalam K, Montell C. TRP channels. *Annu Rev Biochem 2007*. doi:10.1146/annurev.biochem.75.103004.142819

81. Nilius B, Owsianik G. The transient receptor potential family of ion channels. *Genome Biol 2011*. doi:10.1186/gb-2011-12-3-218

82. Montell C, Rubin GM. Molecular characterization of the Drosophila trp locus: a putative integral membrane protein required for phototransduction. *Neuron 1989*. doi:10.1016/0896-6273(89)90069-x

83. Clapham DE. TRP channels as cellular sensors. *Nature 2003*. doi:10.1038/nature02196

84. Liao M et al. Structure of the TRPV1 ion channel determined by electron cryo-microscopy. *Nature 2013*. doi:10.1038/nature12822

85. Liman ER, Corey DP, Dulac C. TRP2: a candidate transduction channel for mammalian pheromone sensory signaling. *Proc Natl Acad Sci U S A 1999*. doi:10.1073/pnas.96.10.5791

86. des Georges A et al. Structural Basis for Gating and Activation of RyR1. *Cell 2016*. doi:10.1016/j.cell.2016.08.075

87. Fan G et al. Gating machinery of InsP3R channels revealed by electron cryomicroscopy. *Nature 2015*. doi:10.1038/nature15249

88. Thompson AJ, Lester HA, Lummis SC. The structural basis of function in Cys-loop receptors. *Q Rev Biophys 2010*. doi:10.1017/s0033583510000168

89. Hibbs RE, Gouaux E. Principles of activation and permeation in an anion-selective Cys-loop receptor. *Nature 2011*. doi:10.1038/nature10139

90. Miller PS, Aricescu AR. Crystal structure of a human GABAA receptor. *Nature 2014*. doi:10.1038/nature13293

91. Du J et al. Glycine receptor mechanism elucidated by electron cryo-microscopy. *Nature 2015*. doi:10.1038/nature14853

92. Hollmann M, Heinemann S. Cloned glutamate receptors. *Annu Rev Neurosci 1994*. doi:10.1146/annurev.ne.17.030194.000335

93. Traynelis SF et al. Glutamate receptor ion channels: structure, regulation, and function. *Pharmacol Rev 2010*. doi:10.1124/pr.109.002451

94. Karakas E, Furukawa H. Crystal structure of a heterotetrameric NMDA receptor ion channel. *Science 2014*. doi:10.1126/science.1251915

95. Lee CH et al. NMDA receptor structures reveal subunit arrangement and pore architecture. *Nature 2014*. doi:10.1038/nature13548

96. Lam HM et al. Glutamate-receptor genes in plants. *Nature 1998*. doi:10.1038/24066

97. Benton R et al. Variant ionotropic glutamate receptors as chemosensory receptors in Drosophila. *Cell 2009*. doi:10.1016/j.cell.2008.12.001

98. Kohda K, Wang Y, Yuzaki M. Mutation of a glutamate receptor motif reveals its role in gating and delta2 receptor channel properties. *Nat Neurosci 2000*. doi:10.1038/73877

99. North RA. Molecular physiology of P2X receptors. *Physiol Rev 2002*. doi:10.1152/physrev.00015.2002

100. Fountain SJ et al. An intracellular P2X receptor required for osmoregulation in Dictyostelium discoideum. *Nature 2007*. doi:10.1038/nature05926

101. Kellenberger S, Schild L. Epithelial sodium channel/degenerin family of ion channels: a variety of functions for a shared structure. *Physiol Rev 2002*. doi:10.1152/physrev.00007.2002

102. Chalfie M, Sulston J. Developmental genetics of the mechanosensory neurons of Caenorhabditis elegans. *Dev Biol 1981*. doi:10.1016/0012-1606(81)90459-0

103. Jentsch TJ et al. Molecular structure and physiological function of chloride channels. *Physiol Rev 2002*. doi:10.1152/physrev.00029.2001

104. Miller C. ClC chloride channels viewed through a transporter lens. *Nature 2006*. doi:10.1038/nature04713

105. Gadsby DC, Vergani P, Csanády L. The ABC protein turned chloride channel whose failure causes cystic fibrosis. *Nature 2006*. doi:10.1038/nature04712

106. Negoda A et al. Conformational change of the extracellular parts of the CFTR protein during channel gating. *Cell Mol Life Sci 2018*. doi:10.1007/s00018-018-2777-0

107. Yang YD et al. TMEM16A confers receptor-activated calcium-dependent chloride conductance. *Nature 2008*. doi:10.1038/nature07313

108. Caputo A et al. TMEM16A, a membrane protein associated with calcium-dependent chloride channel activity. *Science 2008*. doi:10.1126/science.1163518

109. Schroeder BC et al. Expression cloning of TMEM16A as a calcium-activated chloride channel subunit. *Cell 2008*. doi:10.1016/j.cell.2008.09.003

110. Pedemonte N, Galietta LJ. Structure and function of TMEM16 proteins (anoctamins). *Physiol Rev 2014*. doi:10.1152/physrev.00039.2011

111. Sun H et al. The vitelliform macular dystrophy protein defines a new family of chloride channels. *Proc Natl Acad Sci U S A 2002*. doi:10.1073/pnas.052692999

112. Kane Dickson V, Pedi L, Long SB. Structure and insights into the function of a Ca(2+)-activated Cl(-) channel. *Nature 2014*. doi:10.1038/nature13913

113. Sosinsky GE, Nicholson BJ. Structural organization of gap junction channels. *Biochim Biophys Acta 2005*. doi:10.1016/j.bbamem.2005.04.001

114. Taruno A et al. CALHM1 ion channel mediates purinergic neurotransmission of sweet, bitter and umami tastes. *Nature 2013*. doi:10.1038/nature11906

115. Berridge MJ. The Inositol Trisphosphate/Calcium Signaling Pathway in Health and Disease. *Physiol Rev 2016*. doi:10.1152/physrev.00006.2016

116. Yazawa M et al. TRIC channels are essential for Ca2+ handling in intracellular stores. *Nature 2007*. doi:10.1038/nature05928

117. Baughman JM et al. Integrative genomics identifies MCU as an essential component of the mitochondrial calcium uniporter. *Nature 2011*. doi:10.1038/nature10234

118. De Stefani D et al. A forty-kilodalton protein of the inner membrane is the mitochondrial calcium uniporter. *Nature 2011*. doi:10.1038/nature10230

119. Feske S et al. A mutation in Orai1 causes immune deficiency by abrogating CRAC channel function. *Nature 2006*. doi:10.1038/nature04702

120. Vig M et al. CRACM1 is a plasma membrane protein essential for store-operated Ca2+ entry. *Science 2006*. doi:10.1126/science.1127883

121. Hou X et al. Crystal structure of the calcium release-activated calcium channel Orai. *Science 2012*. doi:10.1126/science.1228757

122. Prakriya M, Lewis RS. Store-Operated Calcium Channels. *Physiol Rev 2015*. doi:10.1152/physrev.00020.2014

123. Martinac B, Saimi Y, Kung C. Ion channels in microbes. *Physiol Rev 2008*. doi:10.1152/physrev.00005.2008

124. Liebeskind BJ, Hillis DM, Zakon HH. Evolution of sodium channels predates the origin of nervous systems in animals. *Proc Natl Acad Sci U S A 2011*. doi:10.1073/pnas.1106363108

125. Moran Y et al. Evolution of voltage-gated ion channels at the emergence of Metazoa. *J Exp Biol 2015*. doi:10.1242/jeb.110270

126. Anderson PA, Greenberg RM. Phylogeny of ion channels: clues to structure and function. *Comp Biochem Physiol B Biochem Mol Biol 2001*. doi:10.1016/s1096-4959(01)00376-1

127. Liebeskind BJ, Hillis DM, Zakon HH. Convergence of ion channel genome content in early animal evolution. *Proc Natl Acad Sci U S A 2015*. doi:10.1073/pnas.1501195112

128. Littleton JT, Ganetzky B. Ion channels and synaptic organization: analysis of the Drosophila genome. *Neuron 2000*. doi:10.1016/s0896-6273(00)81135-6

129. Hedrich R. Ion channels in plants. *Physiol Rev 2012*. doi:10.1152/physrev.00038.2011

130. Paysan-Lafosse T et al. InterPro in 2022. *Nucleic Acids Res 2023*. doi:10.1093/nar/gkac993

131. Mistry J et al. Pfam: The protein families database in 2021. *Nucleic Acids Res 2021*. doi:10.1093/nar/gkaa913

132. UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Res 2023*. doi:10.1093/nar/gkac1052

133. Ashcroft FM. From molecule to malady. *Nature 2006*. doi:10.1038/nature04707

134. Kullmann DM. Neurological channelopathies. *Annu Rev Neurosci 2010*. doi:10.1146/annurev-neuro-060909-153122

135. Abriel H, Rougier JS, Jalife J. Ion channel macromolecular complexes in cardiomyocytes: roles in sudden cardiac death. *Circ Res 2015*. doi:10.1161/circresaha.116.305017

136. Overington JP, Al-Lazikani B, Hopkins AL. How many drug targets are there?. *Nat Rev Drug Discov 2006*. doi:10.1038/nrd2199

137. Santos R et al. A comprehensive map of molecular drug targets. *Nat Rev Drug Discov 2017*. doi:10.1038/nrd.2016.230

138. Bagal SK et al. Ion channels as therapeutic targets: a drug discovery perspective. *J Med Chem 2013*. doi:10.1021/jm3011433

139. Wulff H, Castle NA, Pardo LA. Voltage-gated potassium channels as therapeutic targets. *Nat Rev Drug Discov 2009*. doi:10.1038/nrd2983

140. Alexander SPH et al. The Concise Guide to PHARMACOLOGY 2023/24: Ion channels. *Br J Pharmacol 2023*. doi:10.1111/bph.16178

141. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol 2013*. doi:10.1093/molbev/mst010

142. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics 2009*. doi:10.1093/bioinformatics/btp348

143. Minh BQ et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Mol Biol Evol 2020*. doi:10.1093/molbev/msaa015

144. Kalyaanamoorthy S et al. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nat Methods 2017*. doi:10.1038/nmeth.4285

145. Eddy SR. Accelerated Profile HMM Searches. *PLoS Comput Biol 2011*. doi:10.1371/journal.pcbi.1002195

146. van Kempen M et al. Fast and accurate protein structure search with Foldseek. *Nat Biotechnol 2024*. doi:10.1038/s41587-023-01773-0

147. Jumper J et al. Highly accurate protein structure prediction with AlphaFold. *Nature 2021*. doi:10.1038/s41586-021-03819-2

148. Varadi M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Res 2022*. doi:10.1093/nar/gkab1061
