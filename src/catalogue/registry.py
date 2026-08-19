"""Assembly, lookup and self-validation of the ion-channel catalogue.

Everything downstream reads the catalogue through this module: the
classifier, the census, the phylogeny driver, the GUI and every `scripts/s*`
task. Nothing else imports the per-division files directly, so adding a
family is one edit in one division file plus, at most, a hazard row.

`validate()` is run by `python -m src.catalogue` and at the top of the S0
verification script. It checks the things that go wrong silently: a family
pointing at a superfamily that does not exist, the same human gene claimed
by two families, a duplicate key, an exemplar with no way to resolve it, a
`confusable_with` pointing at nothing. None of that is exotic — all of it is
the kind of drift a catalogue accumulates over sessions, and the point of a
machine-checked catalogue is that it cannot.
"""

from __future__ import annotations

from dataclasses import replace

from . import anion, controls, intracellular, largepore, lgic, mechano
from . import other, tmem16_like, vgic_cation, vgic_k
from .hazards import HAZARD_BY_ID, HAZARDS
from .schema import (CENSUS_STATUSES, ChannelFamily, Exemplar, Level,
                     Signature, Status, Superfamily)

_DIVISIONS = (vgic_k, vgic_cation, lgic, anion, tmem16_like, mechano,
              largepore, intracellular, other, controls)

#: Every superfamily, keyed by `Superfamily.key`.
SUPERFAMILIES: dict[str, Superfamily] = {}
for _mod in _DIVISIONS:
    for _sf in getattr(_mod, "SUPERFAMILIES", []):
        SUPERFAMILIES[_sf.key] = _sf
SUPERFAMILIES[vgic_k.PLOOP_SUPERFAMILY.key] = vgic_k.PLOOP_SUPERFAMILY

#: Every family, keyed by `ChannelFamily.key`, in division order.
CATALOGUE: dict[str, ChannelFamily] = {}
for _mod in _DIVISIONS:
    for _fam in getattr(_mod, "FAMILIES", []):
        CATALOGUE[_fam.key] = _fam


# ------------------------------------------------------------- lookups
def family(key: str) -> ChannelFamily:
    return CATALOGUE[key]


def superfamily(key: str) -> Superfamily:
    return SUPERFAMILIES[key]


def families(status: Status | None = None,
             superfamily_key: str = "") -> list[ChannelFamily]:
    out = list(CATALOGUE.values())
    if status is not None:
        out = [f for f in out if f.status is status]
    if superfamily_key:
        out = [f for f in out if f.superfamily == superfamily_key]
    return out


def census_families() -> list[ChannelFamily]:
    """Families that count towards the ion-channel census (D23)."""
    return [f for f in CATALOGUE.values() if f.status in CENSUS_STATUSES]


def control_families() -> list[ChannelFamily]:
    """Everything catalogued *so that it can be excluded* — the decoy panel."""
    return [f for f in CATALOGUE.values() if f.status not in CENSUS_STATUSES]


def human_genes(census_only: bool = True) -> dict[str, str]:
    """`{gene symbol: family key}` for every human gene in the catalogue."""
    src = census_families() if census_only else list(CATALOGUE.values())
    return {g: f.key for f in src for g in f.human_genes}


def family_for_gene(symbol: str) -> ChannelFamily | None:
    """Catalogue lookup by gene symbol.

    Provided for reporting and for building control panels — never for
    classification. Classifying by symbol is what hazard **H15** is about.
    """
    want = symbol.strip().upper()
    for f in CATALOGUE.values():
        if want in {g.upper() for g in f.human_genes + f.other_genes}:
            return f
    return None


def signature_index(census_only: bool = False) -> dict[str, list[str]]:
    """`{pfam accession: [family keys]}` across the whole catalogue."""
    src = census_families() if census_only else list(CATALOGUE.values())
    idx: dict[str, list[str]] = {}
    for f in src:
        for s in f.signatures:
            idx.setdefault(s.accession, []).append(f.key)
    return idx


def shared_signatures(min_families: int = 2) -> dict[str, list[str]]:
    """Signatures carried by several families — the raw material of hazards."""
    return {a: fams for a, fams in signature_index().items()
            if len(fams) >= min_families}


def pore_signatures() -> list[Signature]:
    """Every signature that models a pore-forming module, deduplicated.

    Hazard **H7**: enumerating the P-loop superfamily from `PF00520` alone
    loses TRPC, TRPM, TRPML and TRPP, whose pores Pfam models with different
    accessions. The census enumerates from this union instead.
    """
    seen: dict[str, Signature] = {}
    for f in census_families():
        for s in f.signatures:
            if s.level in (Level.SUPERFAMILY, Level.FAMILY):
                seen.setdefault(s.accession, s)
    return list(seen.values())


def exemplars(census_only: bool = True) -> list[tuple[str, Exemplar]]:
    """`[(family key, exemplar)]` — the reference panel the classifier scores against."""
    src = census_families() if census_only else list(CATALOGUE.values())
    return [(f.key, e) for f in src for e in f.exemplars]


def reference_panel(census_only: bool = True) -> list[tuple[str, str]]:
    """`[(label, uniprot)]` for exemplars that already carry an accession."""
    return [(e.label, e.uniprot) for _, e in exemplars(census_only) if e.uniprot]


def hazards_for(family_key: str) -> list:
    return [h for h in HAZARDS if family_key in h.families]


def with_resolved_exemplars(family_key: str,
                            resolved: dict[str, str]) -> ChannelFamily:
    """Return a copy of a family with `{gene: accession}` filled into exemplars.

    `scripts/s0_catalogue_verify.py` resolves missing accessions live and
    writes them to `results/s0_baseline/exemplars_resolved.tsv`; this applies
    that table in memory so the code file stays a curated claim and the
    resolved values stay a measurement.
    """
    f = CATALOGUE[family_key]
    new = tuple(e if e.uniprot else replace(e, uniprot=resolved.get(e.gene, ""))
                for e in f.exemplars)
    return replace(f, exemplars=new)


# ---------------------------------------------------------- validation
def validate() -> list[str]:
    """Internal consistency of the catalogue. Empty list == clean."""
    problems: list[str] = []
    seen_genes: dict[str, str] = {}
    seen_labels: dict[str, str] = {}

    for key, f in CATALOGUE.items():
        if f.key != key:
            problems.append(f"{key}: key mismatch ({f.key})")
        if f.superfamily and f.superfamily not in SUPERFAMILIES:
            problems.append(f"{key}: unknown superfamily {f.superfamily!r}")
        if f.census_member() and not f.exemplars:
            problems.append(f"{key}: census family with no exemplar")
        if f.census_member() and not (f.human_genes or f.other_genes):
            problems.append(f"{key}: census family with no named gene")
        lo, hi = f.length_band_aa
        if lo and hi and lo >= hi:
            problems.append(f"{key}: bad length band {f.length_band_aa}")
        for g in f.human_genes:
            if g in seen_genes and f.census_member():
                problems.append(
                    f"{key}: human gene {g} already claimed by {seen_genes[g]}")
            seen_genes.setdefault(g, key)
        for e in f.exemplars:
            if e.label in seen_labels:
                problems.append(
                    f"{key}: exemplar label {e.label} reused from {seen_labels[e.label]}")
            seen_labels[e.label] = key
            if not e.gene or not e.species:
                problems.append(f"{key}: exemplar {e.label} is unresolvable")
        for c in f.confusable_with:
            if c not in CATALOGUE and c not in HAZARD_BY_ID:
                problems.append(f"{key}: confusable_with {c!r} matches nothing")

    for h in HAZARDS:
        for fam in h.families:
            if fam not in CATALOGUE:
                problems.append(f"{h.hid}: unknown family {fam!r}")
    return problems


def stats() -> dict[str, int]:
    """Headline counts — what the census claims before any search runs."""
    cens = census_families()
    return {
        "superfamilies": len(SUPERFAMILIES),
        "families_total": len(CATALOGUE),
        "families_census": len(cens),
        "families_control": len(CATALOGUE) - len(cens),
        "human_genes_census": len(human_genes(True)),
        "human_genes_catalogued": len(human_genes(False)),
        "exemplars": len(exemplars(False)),
        "exemplars_with_accession": len(reference_panel(False)),
        "signatures": len(signature_index()),
        "shared_signatures": len(shared_signatures()),
        "hazards": len(HAZARDS),
    }
