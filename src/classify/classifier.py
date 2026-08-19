"""The classifier — three tiers, one call, and a full audit trail.

    architecture  →  what the domain composition allows
    motif         →  what the selectivity filter says
    reference     →  what it is nearest to, with a margin

The tiers are deliberately independent: they use different evidence and fail
in different ways, so agreement between two of them is worth more than a
strong result from one. That is the whole confidence scheme:

    gold        two or more tiers name the same family
    silver      one tier names a family and nothing contradicts it
    bronze      superfamily only, or a family named by a derived rule alone
    unassigned  no tier could say anything

**The gene symbol is never consulted.** It is carried through into the audit
trail so a reader can see it, and the classification would be identical if
it were blank. This is hazard **H15**: symbols group pores with their
accessory subunits (`KCNE`, `CACNB`, `CATSPERB`, `KCTD`), and an unnamed
gene model from a new genome — the case the census actually cares about —
has no symbol at all.

**Disagreement is reported, not resolved by fiat.** When tiers conflict the
call follows a fixed precedence (hazard rule > motif > reference > derived
rule) and the conflict is recorded in `conflicts`, because a classifier that
hides its disagreements cannot be audited.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from ..catalogue import CATALOGUE, Status, hazards_for
from . import motifs as motif_tier
from . import reference as reference_tier
from .rules import RuleMatch, match_architecture, topology_check

GOLD, SILVER, BRONZE, UNASSIGNED = "gold", "silver", "bronze", "unassigned"

#: Which tier wins when two tiers name different families.
TIER_PRECEDENCE = ("hazard", "motif", "reference", "architecture")


@dataclass
class Evidence:
    tier: str
    verdict: str          # family key, superfamily key, or ""
    detail: str
    decisive: bool = False


@dataclass
class ChannelQuery:
    """One protein to classify. `gene_symbol` is recorded and never used."""
    accession: str
    sequence: str = ""
    gene_symbol: str = ""
    species: str = ""
    pfam_counts: dict[str, int] = field(default_factory=dict)
    tm_count: int | None = None
    pore_loops: int | None = None
    length_aa: int | None = None

    def __post_init__(self) -> None:
        if self.length_aa is None and self.sequence:
            self.length_aa = len(self.sequence)


@dataclass
class ChannelCall:
    query: str
    family: str = ""
    superfamily: str = ""
    status: str = ""
    confidence: str = UNASSIGNED
    margin: float = 0.0
    filter_string: str = ""
    evidence: list[Evidence] = field(default_factory=list)
    hazards: list[str] = field(default_factory=list)
    conflicts: list[str] = field(default_factory=list)
    gene_symbol: str = ""       # carried for the reader; never used to decide
    notes: str = ""

    def family_name(self) -> str:
        f = CATALOGUE.get(self.family)
        return f.name if f else ""

    def census_member(self) -> bool:
        f = CATALOGUE.get(self.family)
        return bool(f and f.census_member())

    def text(self) -> str:
        lines = [f"{self.query}  →  "
                 + (f"{self.family} ({self.family_name()})" if self.family
                    else (f"[{self.superfamily} superfamily only]"
                          if self.superfamily else "[unassigned]"))
                 + f"   {self.confidence}"]
        if self.status:
            lines.append(f"  status      : {self.status}")
        if self.margin:
            lines.append(f"  margin      : {self.margin:+.3f}")
        for e in self.evidence:
            mark = "*" if e.decisive else " "
            lines.append(f" {mark}{e.tier:13s} {e.detail}")
        if self.hazards:
            lines.append(f"  hazards     : {', '.join(self.hazards)}")
        if self.conflicts:
            lines.append(f"  CONFLICT    : {'; '.join(self.conflicts)}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["family_name"] = self.family_name()
        return d


def _tier_family(match: RuleMatch) -> tuple[str, str]:
    """(family, tier name) for the architecture result."""
    if match.family:
        hazard_fired = any(r.hazard and r.target == match.family
                           for r in match.fired)
        return match.family, ("hazard" if hazard_fired else "architecture")
    return "", "architecture"


def classify(query: ChannelQuery,
             refs: reference_tier.ReferenceSet | None = None,
             nav_reference: str = "",
             margin: float = reference_tier.DEFAULT_MARGIN,
             leave_one_out: bool = False) -> ChannelCall:
    """Classify one protein.

    `nav_reference` is the SCN5A sequence for the four-repeat filter
    projection; without it the motif tier runs the potassium-filter test
    only. `leave_one_out` drops the query's own accession from the reference
    panel — set it whenever the query might itself be an exemplar, which in
    the S1 benchmark is true of nearly half the panel.
    """
    call = ChannelCall(query=query.accession, gene_symbol=query.gene_symbol)
    votes: dict[str, str] = {}      # tier -> family key

    # -- tier 1: architecture ------------------------------------------
    arch = match_architecture(query.pfam_counts)
    fam, tier_name = _tier_family(arch)
    call.evidence.append(Evidence(tier_name, fam or arch.superfamily,
                                  arch.summary()))
    call.hazards = list(arch.hazards)
    if arch.superfamily:
        call.superfamily = arch.superfamily
    if fam:
        votes[tier_name] = fam

    # -- tier 2: motifs -------------------------------------------------
    mot = motif_tier.scan(query.sequence, nav_reference)
    call.filter_string = mot.filter_string
    if mot.filter_string or mot.k_filter_hits or not mot.available:
        call.evidence.append(Evidence("motif", mot.filter_family, mot.summary()))
    if mot.filter_family:
        votes["motif"] = mot.filter_family

    # -- tier 3: reference identity ------------------------------------
    if refs is not None and query.sequence:
        # Use whatever the architecture tier established. If it named the
        # candidates, break that tie; if it could only reach the superfamily
        # — which is the normal outcome for Cys-loop and iGluR, where every
        # family carries the same two accessions — search inside that
        # superfamily rather than across the whole catalogue. This is using
        # the evidence already gathered, not narrowing to fit an answer: a
        # restriction that matches no exemplar falls back to the full panel.
        restrict = tuple(arch.ambiguous) if arch.ambiguous else ()
        if not restrict and arch.superfamily and not arch.family:
            restrict = tuple(f.key for f in CATALOGUE.values()
                             if f.superfamily == arch.superfamily)
        ref = refs.best_match(
            query.sequence, margin=margin, restrict_to=restrict,
            exclude_accessions=(query.accession,) if leave_one_out else ())
        detail = ref.summary()
        best = ref.best()
        if best is not None and best.uniprot == query.accession:
            # Not an error — the catalogue's own exemplars are legitimate
            # classification targets — but a 100 % hit against yourself is
            # not evidence, and a reader must be able to see that it happened.
            detail += " — SELF-HIT: the query is itself a catalogue exemplar"
        call.evidence.append(Evidence("reference", ref.family, detail))
        call.margin = ref.margin
        if ref.family:
            votes["reference"] = ref.family

    # -- combine --------------------------------------------------------
    chosen = ""
    for tier in TIER_PRECEDENCE:
        if tier in votes:
            chosen = votes[tier]
            for e in call.evidence:
                if e.tier == tier:
                    e.decisive = True
            break
    distinct = set(votes.values())
    if len(distinct) > 1:
        call.conflicts = [f"{t}={f}" for t, f in sorted(votes.items())]

    if chosen:
        call.family = chosen
        fam_rec = CATALOGUE.get(chosen)
        if fam_rec:
            call.superfamily = fam_rec.superfamily or call.superfamily
            call.status = fam_rec.status.value
            call.hazards = sorted(set(call.hazards)
                                  | {h.hid for h in hazards_for(chosen)})
        agree = sum(1 for f in votes.values() if f == chosen)
        if agree >= 2 and not call.conflicts:
            call.confidence = GOLD
        elif agree >= 2:
            call.confidence = SILVER
        elif "hazard" in votes or "motif" in votes:
            call.confidence = SILVER
        else:
            call.confidence = BRONZE if len(distinct) == 1 else SILVER
        ok, detail = topology_check(chosen, query.tm_count, query.pore_loops)
        if query.tm_count is not None:
            call.evidence.append(Evidence("topology", "" if ok else "mismatch",
                                          detail))
            if not ok:
                call.conflicts.append(f"topology: {detail}")
                if call.confidence == GOLD:
                    call.confidence = SILVER
    elif arch.ambiguous:
        call.confidence = BRONZE
        call.notes = "ambiguous between " + ", ".join(arch.ambiguous)
    elif call.superfamily:
        call.confidence = BRONZE
        call.notes = "superfamily evidence only"
    else:
        call.confidence = UNASSIGNED
        call.notes = "no tier produced a call"

    if call.family and CATALOGUE[call.family].status is not Status.CHANNEL:
        call.notes = (call.notes + "; " if call.notes else "") + \
            f"catalogued as {CATALOGUE[call.family].status.value}"
    return call


def classify_all(queries: list[ChannelQuery],
                 refs: reference_tier.ReferenceSet | None = None,
                 nav_reference: str = "",
                 on_progress=None,
                 leave_one_out: bool = False) -> list[ChannelCall]:
    out: list[ChannelCall] = []
    for i, q in enumerate(queries, 1):
        if on_progress:
            on_progress(f"[{i}/{len(queries)}] {q.accession}")
        out.append(classify(q, refs, nav_reference, leave_one_out=leave_one_out))
    return out
