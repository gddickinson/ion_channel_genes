"""The active analysis scope — one family, one superfamily, or every channel.

The parent project had one subject and declared it once in
`src/utils/scope.py`. This project has ninety families, so the equivalent
single place is `src/catalogue/`, and *scope* becomes a run-time choice: a
search, a discovery run or an investigation is scoped to the whole
channelome, to one superfamily, or to one family.

A `Scope` carries exactly what the ported machinery needs — the signatures
that define its search space, the size band, the reference panel, the sister
panel that the discovery scorer measures its margin against, and the
literature keyword. Everything is derived from the catalogue, so scoping to
`nav` and scoping to `all` differ only in which rows were selected.

The module-level constants at the bottom are the compatibility layer:
`src/databases/interpro.py`, `src/discovery/` and `src/investigation/` were
written against the parent project's flat names and keep working, now
defaulting to the whole-channelome scope. Passing a narrower `Scope`
explicitly is always better than relying on them.

**The sister panel is the important field.** In the parent project it held
one family (RYR, the thing that looks like an ITPR). Here it is computed
from the hazard registry: the sister panel of a scope is every family the
catalogue records as confusable with it. That is what turns the parent's
single hard-won D14 test into something that runs for all sixteen hazards.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..catalogue import (CATALOGUE, HAZARDS, SUPERFAMILIES, census_families,
                         exemplars, family, hazards_for, pore_signatures)


@dataclass(frozen=True)
class Scope:
    """What a run is about."""
    key: str                       # "all" | superfamily key | family key
    kind: str                      # "all" | "superfamily" | "family"
    label: str
    family_keys: tuple[str, ...]
    pfam_ids: tuple[str, ...]          # diagnostic signatures in scope
    census_pfam_ids: tuple[str, ...]   # what an enumerated search space uses
    min_length_aa: int
    max_length_aa: int
    reference_panel: tuple[tuple[str, str], ...]   # (label, accession)
    sister_panel: tuple[tuple[str, str], ...]
    sister_family_keys: tuple[str, ...]
    #: Gene symbols of the sister families. The ported discovery scorer
    #: matches *symbols*, not catalogue keys, so passing keys would silently
    #: disable the sister-family test — the one thing the parent project
    #: proved a scorer cannot run without.
    sister_genes: tuple[str, ...]
    known_genes: tuple[str, ...]
    literature_keyword: str
    hazard_ids: tuple[str, ...] = ()
    notes: str = ""

    def describe(self) -> str:
        return (f"scope {self.key} ({self.kind}): {len(self.family_keys)} "
                f"famil{'y' if len(self.family_keys) == 1 else 'ies'}, "
                f"{len(self.pfam_ids)} signature(s), "
                f"{len(self.known_genes)} named gene(s), "
                f"band {self.min_length_aa}-{self.max_length_aa} aa, "
                f"hazards {', '.join(self.hazard_ids) or 'none'}")


def _panel(keys: tuple[str, ...]) -> tuple[tuple[str, str], ...]:
    return tuple((e.label, e.uniprot) for k in keys
                 for e in CATALOGUE[k].exemplars if e.uniprot)


def _band(keys: tuple[str, ...]) -> tuple[int, int]:
    bands = [CATALOGUE[k].length_band_aa for k in keys
             if CATALOGUE[k].length_band_aa[1]]
    if not bands:
        return 0, 0
    return min(b[0] for b in bands), max(b[1] for b in bands)


def scope_for(key: str = "all") -> Scope:
    """Build a scope from a catalogue family key, superfamily key, or "all"."""
    if key in ("", "all", "channelome"):
        keys = tuple(f.key for f in census_families())
        kind, label = "all", "every ion-channel family in the catalogue"
        sisters: tuple[str, ...] = tuple(
            f.key for f in CATALOGUE.values() if not f.census_member())
        hz = tuple(h.hid for h in HAZARDS)
        sigs = tuple(s.accession for s in pore_signatures())
    elif key in SUPERFAMILIES:
        keys = tuple(f.key for f in CATALOGUE.values()
                     if f.superfamily == key and f.census_member())
        kind = "superfamily"
        label = SUPERFAMILIES[key].name
        sisters = tuple(f.key for f in CATALOGUE.values()
                        if f.superfamily == key and not f.census_member())
        hz = tuple(sorted({h.hid for k in keys for h in hazards_for(k)}))
        sigs = tuple(dict.fromkeys(s.accession for k in keys
                                   for s in CATALOGUE[k].signatures))
    elif key in CATALOGUE:
        keys = (key,)
        kind = "family"
        label = CATALOGUE[key].name
        confusable = set(CATALOGUE[key].confusable_with) & set(CATALOGUE)
        for h in hazards_for(key):
            confusable |= set(h.families)
        confusable.discard(key)
        sisters = tuple(sorted(confusable))
        hz = tuple(h.hid for h in hazards_for(key))
        sigs = tuple(s.accession for s in CATALOGUE[key].signatures)
    else:
        raise KeyError(f"unknown scope {key!r} — not a family, superfamily "
                       f"or 'all'")

    lo, hi = _band(keys)
    genes = tuple(dict.fromkeys(g for k in keys
                                for g in CATALOGUE[k].human_genes
                                + CATALOGUE[k].other_genes))
    kw = ("ion channel" if kind == "all"
          else (SUPERFAMILIES[key].name if kind == "superfamily"
                else CATALOGUE[key].name.split("(")[0].strip()))
    return Scope(
        key=key or "all", kind=kind, label=label, family_keys=keys,
        pfam_ids=sigs, census_pfam_ids=sigs,
        min_length_aa=lo or 50, max_length_aa=hi or 6000,
        reference_panel=_panel(keys), sister_panel=_panel(sisters),
        sister_family_keys=sisters,
        sister_genes=tuple(dict.fromkeys(
            g for k in sisters
            for g in CATALOGUE[k].human_genes + CATALOGUE[k].other_genes)),
        known_genes=genes,
        literature_keyword=kw, hazard_ids=hz,
    )


def all_scopes() -> list[str]:
    return ["all"] + sorted(SUPERFAMILIES) + sorted(CATALOGUE)


#: The default scope: every ion-channel family in the catalogue.
DEFAULT_SCOPE = scope_for("all")

# --------------------------------------------------------------------------
# Compatibility layer for machinery ported from the parent project. These are
# the whole-channelome defaults; pass a `Scope` explicitly where you can.
FAMILY_NAME = "ion channels"
FAMILY_LONG_NAME = "the ion-channel superfamilies"
FAMILY_PFAM_IDS = list(DEFAULT_SCOPE.pfam_ids)
PORE_PFAM_IDS = ["PF00520", "PF07885"]
CENSUS_PFAM_IDS = list(DEFAULT_SCOPE.census_pfam_ids)
MIN_LENGTH_AA = DEFAULT_SCOPE.min_length_aa
MAX_LENGTH_AA = DEFAULT_SCOPE.max_length_aa
REFERENCE_PANEL = list(DEFAULT_SCOPE.reference_panel)
SISTER_PANEL = list(DEFAULT_SCOPE.sister_panel)
SISTER_PARALOGS = list(DEFAULT_SCOPE.sister_genes)
KNOWN_PARALOGS = list(DEFAULT_SCOPE.known_genes)
KNOWN_NAME_SUBSTRINGS = [g.lower() for g in KNOWN_PARALOGS]
LITERATURE_KEYWORD = "ion channel"
POR_SEGMENT_FROM_END_AA = 400
