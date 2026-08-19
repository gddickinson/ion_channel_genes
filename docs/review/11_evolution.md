## 11. Evolution

### 11.1 Channels are older than the things they are famous for

Voltage-gated channels are usually introduced through the action potential,
which invites the assumption that they arose with nervous systems. They did
not. Potassium channels are present across bacteria and archaea, and the
prokaryotic structures are not distant analogues but the direct structural
basis of the eukaryotic families [derst1998, martinac2008]. Sodium-channel
homologues predate the origin of nervous systems in animals
[liebeskind2011], and the voltage-gated set was already diversifying at the
emergence of the Metazoa [moran2015]. Sponges and placozoans, which have no
neurons, carry substantial channel repertoires; the anemone *Nematostella*
carries a nearly complete one.

The ordering is therefore the reverse of the intuitive one. The channels
came first and the nervous system was assembled from parts that already
existed, which is why a census restricted to animals with synapses answers
the wrong question.

### 11.2 Repertoires are lineage-specific in both directions

Comparative work on complete genomes shows expansions and losses that
correlate with lifestyle rather than with taxonomic rank [jegla2009,
anderson2001, liebeskind2015]. *C. elegans* has roughly ninety potassium
channel genes; humans have about seventy-nine. The first genome-scale survey
of a channel repertoire, in *Drosophila*, already made the point that the
set is not a scaled-down version of the vertebrate one [littleton2000]. The DEG/ENaC superfamily has
nine human members against roughly thirty in *C. elegans* and thirty
pickpocket genes in *Drosophila*. P2X receptors are present in amoebozoa and
absent from both major invertebrate models [fountain2007]. TRPC2 is a
functional vomeronasal channel in mice and a pseudogene in humans
[liman1999]. NOMPC, the clearest tethered mechanotransducer known, was lost
in mammals [walker2000].

Plants make the point from the other side. They have no voltage-gated
sodium channels and no nervous system, and they have glutamate receptors
[lam1998], two-pore channels, mechanosensitive MscS-like channels and the
OSCA family that was found in plants before it was recognised in animals
[hedrich2012, yuan2014, murthy2018].

### 11.3 Convergence, and the limits of the claim

Four convergences are well enough supported to state, with the caveat that
each rests on failure to detect homology:

- **Potassium selectivity**, achieved by the canonical filter and,
  independently, by TMEM175's isoleucine-lined pathway [cang2015].
- **Intercellular channels**, built by connexins in chordates and by the
  unrelated innexin/pannexin family elsewhere [panchin2000, phelan2001].
- **Mechanosensitivity**, invented at least six times (section 9).
- **Large-pore ATP release**, by pannexin, LRRC8 and CALHM on at least two
  unrelated folds [taruno2013, deneka2018].

The honest formulation is that these are cases where sequence methods find
no relationship. Structure comparison and profile-profile search can push
the detection limit further, and the TMEM16/OSCA/TMC clan — three families
assigned separate domain models that turned out to share a fold — is the
warning that they sometimes do.

![](figures/fig7_forest.png)

**Figure 7 | A forest, not a tree.** Every superfamily with a census family,
drawn as a unit that gets its own tree. Four are marked *no tree*: the
catalogue records them as non-alignable — their member families share a fold
with no detectable sequence homology — and `build_tier2()` raises rather
than producing an alignment artefact. Violet lines are the structural
relationships §11.3 discusses, drawn as network edges with no branch length
and no ancestor implied. There is no panel spanning two boxes, and no code
path that would draw one (**D27**). *Schematic of a rule, not a result.*

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
