"""The TMEM16/OSCA/TMC structural clan — one fold, three families, no
detectable sequence homology between them.

Anoctamins (TMEM16), the OSCA/TMEM63 mechanosensitive channels and the
transmembrane channel-like (TMC) proteins were each characterised
independently, and each was assigned its own Pfam model: `PF04547`,
`PF02714` and `PF07810` share no residues. When their structures arrived
they turned out to be the same ten-transmembrane fold, with the conduction
pathway in the same place — a groove between TM4 and TM6 that is open to
the lipid in the scramblases and closed into a pore in the channels.

This is the catalogue's clearest instance of the rule that decides the whole
phylogeny protocol: **a shared fold with no alignable sequence is a fold
edge, not a tree edge** (**D27**). The three families are recorded under one
superfamily with `alignable=False`, which routes them to the structural
network in tier 3 and forbids the pipeline from ever putting them in one
alignment. Their common ancestry is almost certainly real; the data that
would let us reconstruct it is not in the sequences.

The second lesson here is `ano_channel` versus `ano_scramblase`. ANO1 and
ANO6 have *identical* Pfam architectures (`PF04547` + `PF16178`, measured),
similar lengths and 40 % identity, and one is a calcium-activated chloride
channel while the other is a phospholipid scramblase that also passes ions.
No sequence-level classifier separates them. They are split here on
literature-assigned function, flagged as such, and the roadmap asks whether
the split survives a phylogeny (**Q7**).
"""

from __future__ import annotations

from .schema import ChannelFamily as CF
from .schema import Exemplar as Ex
from .schema import Fold
from .schema import Gating as G
from .schema import Level as L
from .schema import Provenance as P
from .schema import Selectivity as Sel
from .schema import Signature as Sig
from .schema import Status as St
from .schema import Superfamily as SF

ANOCTAMIN = Sig("pfam", "PF04547", "Anoctamin", L.FAMILY, 1, P.DB)
ANO_DIMER = Sig("pfam", "PF16178", "Anoct_dimer", L.FAMILY, 1, P.DB)
RSN1_7TM = Sig("pfam", "PF02714", "RSN1_7TM", L.FAMILY, 1, P.DB)
TMC_DOM = Sig("pfam", "PF07810", "TMC", L.FAMILY, 1, P.DB)

SUPERFAMILIES: list[SF] = [
    SF("tmem16_like", "TMEM16 / OSCA / TMC structural clan", Fold.TMEM16,
       alignable=False,
       anchor_module="TM4–TM6 conduction groove (structure only)",
       shared_signatures=(),
       notes=("Homology asserted from structure, undetectable in sequence. "
              "`alignable=False` is enforced by the phylogeny driver: a tier-2 "
              "tree across this superfamily is refused, not attempted and "
              "caveated.")),
]

FAMILIES: list[CF] = [
    CF(
        key="ano_channel",
        name="Anoctamin channels (TMEM16A, TMEM16B)",
        superfamily="tmem16_like", status=St.CHANNEL, fold=Fold.TMEM16,
        selectivity=Sel.ANION, gating=(G.CALCIUM, G.VOLTAGE),
        stoichiometry="dimer, one pore per subunit", tm_per_subunit=10,
        signatures=(ANOCTAMIN, ANO_DIMER),
        exemplars=(Ex("Hs_ANO1", "ANO1", "Homo sapiens", "Q5XXA6"),
                   Ex("Hs_ANO2", "ANO2", "Homo sapiens", "Q9NQ90")),
        human_genes=("ANO1", "ANO2"),
        length_band_aa=(950, 1010), confusable_with=("ano_scramblase", "H6"),
        notes=("The long-sought calcium-activated chloride channel of "
               "secretory epithelia and olfactory transduction."),
    ),
    CF(
        key="ano_scramblase",
        name="Anoctamin scramblases (TMEM16C–K)",
        superfamily="tmem16_like", status=St.CHANNEL_CONTESTED, fold=Fold.TMEM16,
        selectivity=Sel.ANION, gating=(G.CALCIUM,),
        stoichiometry="dimer", tm_per_subunit=10,
        signatures=(ANOCTAMIN, ANO_DIMER),
        exemplars=(Ex("Hs_ANO6", "ANO6", "Homo sapiens", "Q4KMQ2"),),
        human_genes=("ANO3", "ANO4", "ANO5", "ANO6", "ANO7",
                     "ANO8", "ANO9", "ANO10"),
        length_band_aa=(600, 1000), confusable_with=("ano_channel", "H6"),
        notes=("Phospholipid scramblases that also conduct ions; ANO6 loss "
               "causes Scott syndrome. Split from `ano_channel` on literature "
               "function, not on any sequence feature — the split is a "
               "hypothesis the phylogeny can test, and it is marked so nobody "
               "downstream mistakes it for a measurement."),
    ),
    CF(
        key="osca_tmem63",
        name="OSCA / TMEM63 mechanosensitive channels",
        superfamily="tmem16_like", status=St.CHANNEL, fold=Fold.TMEM16,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL,),
        stoichiometry="dimer (TMEM63 functions as a monomer)", tm_per_subunit=11,
        signatures=(RSN1_7TM,
                    Sig("pfam", "PF13967", "RSN1_TM", L.FAMILY, 1, P.DB),
                    Sig("pfam", "PF14703", "PHM7_cyt", L.FAMILY, 1, P.DB)),
        exemplars=(Ex("Hs_TMEM63A", "TMEM63A", "Homo sapiens", "O94886"),
                   Ex("Hs_TMEM63B", "TMEM63B", "Homo sapiens", "Q5T3F8"),
                   Ex("At_OSCA1.2", "OSCA1.2", "Arabidopsis thaliana", "Q5XEZ5")),
        human_genes=("TMEM63A", "TMEM63B", "TMEM63C"),
        other_genes=("OSCA1.1", "OSCA1.2"),
        length_band_aa=(730, 840),
        notes=("Discovered in plants as hyperosmolality-gated calcium channels "
               "and only then recognised in animals — the family exists across "
               "eukaryotes, so a human-only census would have called it absent. "
               "Pfam still names its core domain a 'putative phosphate "
               "transporter', a stale annotation the catalogue records rather "
               "than corrects silently."),
    ),
    CF(
        key="tmc",
        name="Transmembrane channel-like proteins (TMC1–TMC8)",
        superfamily="tmem16_like", status=St.CHANNEL_CONTESTED, fold=Fold.TMEM16,
        selectivity=Sel.CATION_NONSELECTIVE, gating=(G.MECHANICAL,),
        stoichiometry="dimer", tm_per_subunit=10,
        signatures=(TMC_DOM,),
        exemplars=(Ex("Hs_TMC1", "TMC1", "Homo sapiens", "Q8TDI8"),
                   Ex("Hs_TMC2", "TMC2", "Homo sapiens", "Q8TDI7")),
        human_genes=("TMC1", "TMC2", "TMC3", "TMC4", "TMC5",
                     "TMC6", "TMC7", "TMC8"),
        length_band_aa=(700, 1150),
        notes=("TMC1/TMC2 are the pore of the hair-cell mechanotransduction "
               "channel — the strongest claim in the family and still not a "
               "solved structure of a conducting state, which is why the family "
               "is flagged contested. TMC6/TMC8 (EVER1/EVER2) do something "
               "else entirely in keratinocytes: one Pfam domain, two "
               "biologies."),
    ),
]
