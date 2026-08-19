## 3. Selectivity

### 3.1 The potassium filter

The best-understood selectivity mechanism in biology is the potassium
filter: a short, absolutely conserved sequence — T-x-G-Y-G — whose backbone
carbonyls line a narrow pore and coordinate a dehydrated K⁺ ion in the same
geometry water would [heginbotham1994, doyle1998, zhou2001]. Sodium, which
is smaller, cannot be coordinated at that spacing and pays a dehydration
penalty it does not recover; selectivity is achieved by a structure that
fits the larger ion, which is the opposite of a sieve. The filter is
recognisable in *Streptomyces* KcsA, in *Drosophila* Shaker and in human
Kv1.1, and it is the one motif in this review that can be found in a raw
sequence with no reference and no annotation [mackinnon2003].

![](figures/fig2_potassium_filter.png)

**Figure 2 | One motif, one billion years.** The T-x-G-Y-G filter located
directly in each sequence — no alignment was performed — with seven residues
of context on each side, from the *Streptomyces* prototype through
*Drosophila* Shaker to the human families. Residues are coloured by
chemistry. Two rows carry the argument. *Bacillus* **NaK** is the only panel
member that fails the test, and its filter reads TVGDG: a single Y→D
substitution, and it is not potassium-selective. **GluR0**, a cyanobacterial
*glutamate receptor*, carries a perfect TVGYG — which is the inverted-P-loop
relationship of §2 visible in a sequence rather than in a structure.
Rendered from `results/s0_baseline/filter_k.tsv`.

### 3.2 The four-repeat locus

The sodium and calcium channels arrived at selectivity differently. Both are
pseudo-tetramers — four homologous 6TM repeats in a single chain, the
product of two rounds of internal duplication from a Kv-like ancestor — and
their selectivity is set by one residue contributed by each repeat, at the
position corresponding to the potassium filter. In sodium channels those
four residues read D-E-K-A; in the high-voltage-activated calcium channels
E-E-E-E; and swapping them converts one into the other [heinemann1992].

This locus is the only reliable way to tell the four-repeat families apart.
Their domain architectures are identical, their lengths overlap, and no
domain model distinguishes them. NALCN, a sodium leak channel, reads E-E-K-E
— a calcium-channel filter with a sodium-channel lysine in the third repeat
[lu2007]. The T-type calcium channels read E-E-D-D rather than E-E-E-E,
which is a real difference within a family rather than an error.

Several of the families whose filters are least well characterised are also
the primary transducers of noxious stimuli, which is where much of their
pharmacological interest lies [mccleskey2005].

![](figures/fig3_four_repeat_filter.png)

**Figure 3 | Four residues, three ions.** The selectivity locus of the
four-repeat channels, projected from human Nav1.5 by MAFFT and shown with
its context in each query (decision **D26**). The projected residue is
boxed; the signature is given at the right. The two sodium channels read
DEKA, the L-type calcium channel EEEE, the T-type channel **EEDD** — a real
difference within a family, not an error — and the sodium leak channel
NALCN reads EEKE, a calcium-channel filter carrying the sodium channel's
lysine in repeat III. TPC1 and CatSper1 are included because they have two
repeats and one: the method declines rather than guessing, and a figure
showing only its successes would misrepresent it. Rendered from
`results/s0_baseline/filter_four_repeat.tsv`.

### 3.3 Selectivity that is not in the family

Two observations set a limit on how far family membership predicts what a
channel conducts.

The first is that charge selectivity in the Cys-loop receptors — the
difference between an excitatory nicotinic receptor and an inhibitory GABA
or glycine receptor — is set by a short ring of residues at the intracellular
end of the pore-lining helix, and can be inverted by three substitutions
[galzi1992]. Cation- and anion-selective members sit in the same superfamily
with the same fold.

The second is that AMPA-receptor calcium permeability is not encoded in the
genome at all: it is set by RNA editing at a single position in GluA2, and a
genomic census cannot see it [sommer1991].

The practical consequence is that selectivity is a property recorded per
family from experiment, never predicted per sequence. Any pipeline that
outputs a selectivity call from sequence alone is reporting family
membership with extra steps.

### 3.4 Potassium selectivity twice

TMEM175, the endolysosomal potassium channel, is selective for K⁺ and has no
T-x-G-Y-G filter and no structural relationship to the P-loop channels; its
conduction pathway is lined by isoleucines [cang2015, lee2017tmem175]. The
canonical filter is therefore one solution among at least two, and its near
universality reflects the dominance of one superfamily rather than a physical
necessity.
