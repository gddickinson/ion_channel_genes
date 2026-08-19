## 6. The ligand-gated superfamilies

Four unrelated superfamilies open when a transmitter binds. They are worth
treating together only because their traps have the same shape: in each, the
ligand-binding domain has an independent existence outside channels, and a
search built on it returns proteins that do not conduct.

### 6.1 Pentameric (Cys-loop) receptors

Five subunits, each a β-sandwich ectodomain over four transmembrane helices,
with the second helix lining the pore [unwin2005, thompson2010]. The human
complement is large — nicotinic acetylcholine receptors, GABA_A receptors,
glycine receptors, the ionotropic serotonin receptor and the
zinc-activated channel — and the superfamily is alignable end to end, which
makes it the most tractable phylogenetic unit in this review.

Two members are informative precisely because they are not what the family
name implies. *AChBP*, secreted by snail glia, is a complete pentameric
ligand-binding domain with no transmembrane region at all, and its structure
made nicotinic pharmacology interpretable [brejc2001]. The bacterial
channels GLIC and ELIC lack the disulfide-bonded loop the superfamily is
named after [hilf2008, bocquet2009] — the defining feature is not universal
within the group it defines. Both facts matter for classification: a rule
requiring the transmembrane domain alongside the ligand-binding domain
correctly rejects AChBP and, as measured in section 12, also rejects the
bacterial channels that root the tree.

Invertebrates add glutamate- and histamine-gated *chloride* channels absent
from vertebrates — the ivermectin target — and one of them provided the
first eukaryotic structure of the superfamily [hibbs2011]. Human structures
followed for GABA_A and glycine receptors [miller2014, du2015].

![](figures/fig6_cysloop_tree.png)

**Figure 6 | A tier-2 tree that behaves.** Maximum-likelihood tree of the
Cys-loop reference proteins — MAFFT, trimAl, IQ-TREE 2 with ModelFinder
(LG+G4) and 1000 ultrafast bootstrap replicates — rooted on the two
bacterial channels. Numbers are bootstrap support. The anion-selective GABA
and glycine receptors form a clade at 100 %, the cation-selective receptors
another, and **AChBP — which is not a channel — sisters the cationic group
at 100 %**, exactly where its cholinergic ligand-binding chemistry says it
should sit. This is the superfamily the review calls the most tractable
phylogenetic unit in the subject, and the tree is what that claim looks
like. Run recorded in `results/phylogeny/tier2_cysloop/cysloop.run.json`.

### 6.2 Ionotropic glutamate receptors

Four subunits, each with an amino-terminal clamshell, a second clamshell
that binds glutamate, and a pore module inserted the other way up
[hollmann1994, sobolevsky2009, traynelis2010]. The AMPA, kainate, NMDA and
delta families differ in ligand requirements, kinetics and — for NMDA —
the coincidence detection that comes from needing glutamate, glycine and
depolarisation together [karakas2014, lee2014].

The trap here is the clamshell. It is a bacterial periplasmic
binding-protein fold, and it is also the ligand-binding domain of every
class C G-protein-coupled receptor: the metabotropic glutamate receptors,
both GABA_B subunits, the calcium-sensing receptor and the sweet and umami
taste receptors. A search for glutamate receptors by clamshell returns a
receptor family that is not ionotropic at all, and the discriminating
evidence is the pore region.

The superfamily also extends well past animals with synapses. Plants carry
glutamate-receptor genes [lam1998], and the insect ionotropic receptors are
divergent iGluRs that serve chemosensation and have mostly lost the
clamshell — so a clamshell-based search misses them entirely [benton2009].
The delta receptors do not gate to glutamate at all, and are known to have a
functional pore only because a gain-of-function mutation holds it open
[kohda2000].

### 6.3 P2X receptors

Trimers gated by extracellular ATP, with two transmembrane helices per
subunit, both termini intracellular, and a large disulfide-rich ectodomain
shaped, in the original description, like a leaping dolphin [kawate2009,
north2002]. Their distribution is patchy in an interesting way: present in
amoebozoa — where a *Dictyostelium* member works on intracellular vesicles
rather than the plasma membrane [fountain2007] — and absent from *Drosophila*
and *C. elegans*, which makes P2X a good test of whether a census can
distinguish gene loss from database absence.

### 6.4 DEG/ENaC channels

Trimers again, of an unrelated fold, and gated by almost everything except
voltage: protons in the acid-sensing channels [jasti2007], peptides in the
molluscan FaNaC, mechanical force in the *C. elegans* degenerins, and in the
epithelial sodium channel essentially nothing — ENaC is constitutively open
and regulated by proteolysis and trafficking [canessa1994, kellenberger2002].
The superfamily takes its name from gain-of-function alleles that kill the
neuron they are expressed in, found in the first genetic screens for touch
insensitivity [chalfie1981]. Its size is strongly lineage-dependent: nine
members in humans against roughly thirty in *C. elegans* and thirty
pickpocket genes in *Drosophila*, so a human-only count badly
misrepresents it.
