"""The types the ion-channel catalogue is written in.

One `ChannelFamily` record per family, grouped into `Superfamily` records.
Everything the project treats as an *assertion about biology* lives in these
records, so a claim can be traced to one line of one file and checked
against a live database by `scripts/s0_catalogue_verify.py`.

Three ideas carry most of the design:

**Signatures have a level.** A Pfam accession is only evidence for the level
at which it is diagnostic. `PF00520` (Ion_trans) is carried by Nav, Cav,
Kv, TRP, CNG, HCN and NALCN alike, so it is a `SUPERFAMILY` signature — it
says "P-loop channel", never "this is a Nav". Recording the level in the
data is what stops the classifier from making a family call on superfamily
evidence (see `docs/classification_rules.md`, decision **D25**).

**Status is not the same as membership.** A family is in the catalogue
because a search for ion channels will hit it, not because it is one. The
CLC transporters, the anoctamin scramblases, the auxiliary subunits, the
aquaporins and the periplasmic-binding-protein relatives of the iGluR
clamshell are all catalogued *and* marked, so an excluded thing is excluded
by a recorded decision rather than by never having been considered.

**Nothing here is trusted until it is verified.** Every accession carries a
`Provenance` tag; `results/s0_baseline/catalogue_verification.tsv` records
what a live check made of it. Curated-but-unverified is a legitimate state,
silently-wrong is not.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class Level(str, Enum):
    """The taxonomic level at which a signature is diagnostic."""
    SUPERFAMILY = "superfamily"   # shared across the superfamily; never a family call
    FAMILY = "family"             # diagnostic for this family
    SUBFAMILY = "subfamily"       # diagnostic for a subset of the family
    ACCESSORY = "accessory"       # present, but common to unrelated proteins
    SHARED_WITH_DECOY = "shared_with_decoy"  # also carried by a known non-channel


class Provenance(str, Enum):
    """How far an assertion can be trusted, in the parent project's grammar."""
    DB = "db"        # re-derived from a live database by an S0 script
    LIT = "lit"      # literature; audited in S0
    CURATED = "curated"   # hand-written here, not yet checked
    OPEN = "open"    # a question this project answers


class Gating(str, Enum):
    VOLTAGE = "voltage"
    LIGAND_EXTRACELLULAR = "ligand_extracellular"
    LIGAND_INTRACELLULAR = "ligand_intracellular"
    CYCLIC_NUCLEOTIDE = "cyclic_nucleotide"
    CALCIUM = "calcium"
    MECHANICAL = "mechanical"
    THERMAL = "thermal"
    PROTON = "proton"
    LIPID = "lipid"
    VOLUME = "volume"
    PHOSPHORYLATION = "phosphorylation"
    LIGHT = "light"
    LEAK = "leak"
    UNKNOWN = "unknown"


class Selectivity(str, Enum):
    K = "K+"
    NA = "Na+"
    CA = "Ca2+"
    ANION = "anion"
    CATION_NONSELECTIVE = "cation_nonselective"
    PROTON = "H+"
    LARGE_PORE = "large_pore"
    WATER = "water"
    UNKNOWN = "unknown"


class Fold(str, Enum):
    """Structural fold of the pore-forming module.

    The fold, not the sequence, is what makes two superfamilies comparable
    at all — and folds shared without detectable sequence homology are
    exactly the cases the phylogeny protocol forbids putting in one tree
    (decision **D27**).
    """
    P_LOOP = "p_loop"                 # VGIC superfamily: S5-P-S6 pore module
    CYS_LOOP = "cys_loop"             # pentameric LGIC
    IGLUR = "iglur"                   # tetrameric LGIC, clamshell + inverted pore
    P2X = "p2x"                       # trimeric, dolphin fold
    ENAC_DEG = "enac_deg"             # trimeric, hand fold
    CLC = "clc"                       # antiparallel double-barrelled
    TMEM16 = "tmem16"                 # scramblase/channel fold (also OSCA, TMC)
    PIEZO = "piezo"                   # three-bladed propeller
    BETA_BARREL = "beta_barrel"       # VDAC, bacterial porins
    CONNEXIN = "connexin"             # hexameric hemichannel
    INNEXIN = "innexin"               # octameric; pannexin fold
    CALHM = "calhm"                   # octa/undecameric large pore
    MSC = "msc"                       # MscL / MscS mechanosensitive
    AQUAPORIN = "aquaporin"           # hourglass; water, not ions
    ABC = "abc"                       # CFTR only
    TMEM175 = "tmem175"               # lysosomal K+, unrelated to P-loop
    GST = "gst"                       # CLIC soluble fold
    ORAI = "orai"                     # hexameric CRAC
    OTOPETRIN = "otopetrin"           # double-barrel proton channel
    VSD_ONLY = "vsd_only"             # Hv1: voltage sensor, no pore domain
    TRIC = "tric"                     # trimeric intracellular cation
    MCU = "mcu"                       # mitochondrial Ca2+ uniporter
    BESTROPHIN = "bestrophin"
    LRRC8 = "lrrc8"                   # VRAC; connexin-like fold
    VIROPORIN = "viroporin"
    UNKNOWN = "unknown"


class Status(str, Enum):
    """Why a record is in the catalogue."""
    CHANNEL = "channel"                       # pore-forming ion channel; in the census
    CHANNEL_CONTESTED = "channel_contested"   # channel activity disputed; in, flagged
    CHANNEL_ASSOCIATED = "channel_associated" # auxiliary subunit; not a pore
    TRANSPORTER = "transporter"               # same fold, coupled transport
    NON_CHANNEL_HOMOLOG = "non_channel_homolog"  # shares a signature, is not a channel
    OUT_OF_SCOPE = "out_of_scope"             # a channel, but not an ion channel


#: Statuses that count towards the census denominator.
CENSUS_STATUSES = (Status.CHANNEL, Status.CHANNEL_CONTESTED)


@dataclass(frozen=True)
class Signature:
    """One domain signature, with the level at which it means something."""
    db: str                 # "pfam" | "interpro"
    accession: str          # "PF00520"
    name: str               # "Ion_trans"
    level: Level = Level.FAMILY
    copies: int = 1         # expected copies per subunit (4 for Nav/Cav)
    provenance: Provenance = Provenance.CURATED
    note: str = ""
    # Whether the census enumerates from this signature. Separate from
    # `level` on purpose: level says what the signature *proves* (D25/D33),
    # enumeration says where the census *looks*. None = derive from level
    # (SUPERFAMILY/FAMILY enumerate). S2b re-levelled PF20519 to
    # SHARED_WITH_DECOY and, through the old coupling, silently removed it
    # from the census search space.
    enumerate: bool | None = None

    def key(self) -> str:
        return f"{self.db}:{self.accession}"

    @property
    def enumerates(self) -> bool:
        if self.enumerate is not None:
            return self.enumerate
        return self.level in (Level.SUPERFAMILY, Level.FAMILY)


@dataclass(frozen=True)
class Exemplar:
    """A reference protein for a family: the thing new sequences are scored against.

    `uniprot` may be empty when the record is curated from a gene symbol
    alone; `scripts/s0_catalogue_verify.py` resolves and fills the gap in
    `results/s0_baseline/exemplars_resolved.tsv` rather than editing code.
    """
    label: str              # "Hs_KCNA1"
    gene: str               # "KCNA1"
    species: str            # "Homo sapiens"
    uniprot: str = ""       # "Q09470"
    note: str = ""


@dataclass(frozen=True)
class ChannelFamily:
    """One family of the catalogue."""
    key: str                       # "kv_shaker"
    name: str                      # "Shaker-related K+ channels (Kv1-Kv4)"
    superfamily: str               # key into SUPERFAMILIES
    status: Status
    fold: Fold
    selectivity: Selectivity
    gating: tuple[Gating, ...]
    stoichiometry: str = "tetramer"
    tm_per_subunit: int = 0        # 0 = not applicable / unknown
    pore_loops_per_subunit: int = 0
    signatures: tuple[Signature, ...] = ()
    exemplars: tuple[Exemplar, ...] = ()
    human_genes: tuple[str, ...] = ()
    other_genes: tuple[str, ...] = ()      # non-human members worth naming
    length_band_aa: tuple[int, int] = (0, 0)
    filter_motif: str = ""         # selectivity-filter signature, if one exists
    iuphar_class: str = ""         # the IUPHAR/BPS name, where there is one
    confusable_with: tuple[str, ...] = ()  # family keys / hazard ids
    provenance: Provenance = Provenance.CURATED
    notes: str = ""

    def census_member(self) -> bool:
        return self.status in CENSUS_STATUSES

    def diagnostic_signatures(self) -> tuple[Signature, ...]:
        return tuple(s for s in self.signatures
                     if s.level in (Level.FAMILY, Level.SUBFAMILY))

    def superfamily_signatures(self) -> tuple[Signature, ...]:
        return tuple(s for s in self.signatures if s.level is Level.SUPERFAMILY)

    def architecture(self) -> str:
        """The domain pattern, with the non-diagnostic parts marked.

        Signatures shared with a known decoy are shown in brackets: they are
        part of the architecture and are not evidence for the family, and a
        display that omitted them would make the family look more
        distinguishable than it is. `PF08016` on TRPML is the case in point —
        it is also on TRPP and on polycystin-1, which is not a pore.
        """
        parts = []
        for s in self.signatures:
            if s.level is Level.ACCESSORY:
                continue
            token = f"{s.accession}×{s.copies}" if s.copies > 1 else s.accession
            parts.append(f"[{token}]" if s.level is Level.SHARED_WITH_DECOY
                         else token)
        return " + ".join(parts)


@dataclass(frozen=True)
class Superfamily:
    """A group of families with a common pore fold and a common ancestor.

    `alignable` is the load-bearing field: it says whether the member
    families can be put in one sequence alignment at all. Where it is False
    the relationship is structural only, and the phylogeny protocol routes
    the comparison to the fold network instead of to a tree (**D27**).
    """
    key: str
    name: str
    fold: Fold
    alignable: bool = True
    anchor_module: str = ""        # the region a cross-family alignment uses
    shared_signatures: tuple[Signature, ...] = ()
    root_with: tuple[str, ...] = ()   # outgroup family keys for rooting
    notes: str = ""
    #: How S6 extracts the tier-2 module from a member (D40), read from the
    #: UniProt topology of annotated references: "pore_loop" = the TM helix
    #: before each re-entrant (INTRAMEM) pore loop through the TM helix after
    #: it, one module per loop; "tm_span" = first TM start to last TM end;
    #: "" = the superfamily aligns full length and needs no module.
    module_rule: str = ""


@dataclass(frozen=True)
class Hazard:
    """A recorded way to get the classification wrong.

    The parent project (`../ip3r_genes`) carried exactly one of these —
    ITPR versus RYR — and it was load-bearing enough to become decision D14
    there. Across all ion channels there are dozens, so they are data, not
    prose: every hazard names the discriminating test the classifier must
    apply, and `scripts/s1_benchmark.py` measures whether that test works.
    """
    hid: str                       # "H1"
    title: str
    families: tuple[str, ...]      # catalogue keys involved
    shared_evidence: str           # what makes them look alike
    discriminator: str             # the positive test that separates them
    test_owner: str = ""           # which module implements the test
    severity: str = "high"         # high | medium | low
    provenance: Provenance = Provenance.CURATED
    notes: str = ""
