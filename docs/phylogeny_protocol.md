# The phylogeny protocol — a forest, not a tree

## The problem, stated plainly

There is no alignment that contains a nicotinic acetylcholine receptor and a
Kv channel. They are both ion channels; they share no ancestor detectable in
sequence, no common fold, and no alignable position. The same is true of the
P2X receptors, the CLC family, the anoctamins, Piezo, the connexins, the
otopetrins and about twenty other groups.

"The phylogeny of ion channels" is therefore not one object. It is roughly
twenty-five objects plus a statement about how they relate, and the statement
is structural rather than phylogenetic.

Every published tree that spans non-homologous channel superfamilies is
measuring alignment artefacts. The only reliable way not to produce one by
accident is to make it impossible to ask for, which is what
`src/phylo/forest.py` does.

## Three tiers

### Tier 1 — within a family

Full-length alignment (MAFFT L-INS-i for ≤ 200 sequences), trimAl
`-automated1`, IQ-TREE 2 with ModelFinder and 1000 ultrafast bootstraps,
rooted on the sister family declared in the catalogue's `Superfamily.root_with`.

Conventional, uncontroversial, and where nearly every biologically useful
result lives: the Kv1–Kv4 duplications, the teleost 3R copies, the loss of
TRPC2 in humans, whether ITPR1 and ITPR3 are sisters.

### Tier 2 — within a superfamily

Full-length alignment is meaningless here: inside the P-loop superfamily
alone the members run from 400 aa (Kir) to 5,038 aa (RYR1), and almost all
of the difference is cytosolic machinery with no counterpart in the other
families. So tier 2 aligns the **pore module** only — the last two
transmembrane helices and the re-entrant loop between them, ~120 residues,
one per subunit (four per chain in Nav/Cav) — extracted either from the
InterPro domain envelope or by projection from a reference
(`src/phylo/modules.py`), with the method recorded per sequence.

**The result is labelled a pore-module tree in every output**, because that
is what it is. `TreeRun.label()` will not produce any other name for it.

A few superfamilies are the exception, and the catalogue says which: the
Cys-loop receptors are alignable end to end, and `Superfamily.anchor_module`
records that with the words *full-length*. `needs_modules(key)` reads it,
`build_tier2()` defaults to it, and a superfamily aligned full-length is
labelled a **protein** tree rather than a pore-module one. The label follows
the method, not the tier.

### Tier 3 — between superfamilies

**Refused.** `build_tier2()` raises `NotAlignable` for any superfamily the
catalogue marks `alignable=False`, and there is no `build_tier3()` at all.

What replaces it is `src/phylo/network.py`: a **fold-similarity network**
whose nodes are families and whose edges are structural similarity — Foldseek
TM-score or E-value between representative predicted structures. No branch
lengths. No support values. No edge implies a common ancestor. A missing edge
means "not measured", never "not related", and the manifest keeps those two
apart.

The relationships the network exists to hold are real and interesting:

- TMEM16, OSCA/TMEM63 and TMC share a ten-TM fold with the conduction groove
  in the same place, and share no sequence.
- The iGluR pore is an inverted Kir pore.
- Pannexin, innexin and LRRC8 are homologous (Pfam names LRRC8's TM region
  "Pannexin-like"); connexins build the same kind of pore and are unrelated
  to all of them.
- Hv1's voltage sensor is homologous to the P-loop channels' VSD, and Hv1 has
  no pore domain.

Each is entered as `measured=False` until S8 measures it, so the network
cannot ship with its conclusions pre-drawn.

## Rooting

| unit | outgroup | why |
|------|----------|-----|
| P-loop, K⁺ branch | `kcsa_prok` | prokaryotic 2TM channels predate the eukaryotic families |
| Cys-loop | `plgic_prok` (GLIC, ELIC) | bacterial pentamers; they lack the eponymous disulfide, which is itself a result |
| iGluR | `iglur_prok` (GluR0) | a K⁺-selective glutamate receptor: the missing link |
| P2X | `p2x_nonmetazoan` | *Dictyostelium* P2X |
| DEG/ENaC | `deg_invertebrate` | the *C. elegans* degenerins |
| CLC | `clc_prokaryotic` | ClC-ec1 — and the prototype is an antiporter, so transport is the ancestral state |
| ca_release | `itpr` | ITPR roots RYR within the superfamily |

Outgroups come from the catalogue, not from the analyst, and a run whose
sequence set does not include its declared outgroup is unrooted and reported
as unrooted.

## The rules

**D27 — the phylogeny is a forest, and non-alignable comparisons go to the
network.** Enforced in code, not by convention.

**D8 — representatives are chosen per clade × per kingdom** by a rule in a
script, never "the longest sequence per species" (inherited: longest-first
reliably selects chimeric gene models).

**D15 — the species tree is an input, not a result**, with a source on every
calibrated node.

**D5 — bait panels are screened by label, not padded for breadth**; a bait
wearing the wrong clade name does real damage.

**D28 — no tree without the tools.** Missing IQ-TREE means an alignment is
written and the tree is not, with the reason in `TreeRun.note`. The
alternative — quietly falling back to neighbour-joining — would put two
incomparable methods in one figure.
