"""Selectivity-filter tests — the discriminators that domain architecture cannot make.

Two motifs carry most of the work.

**The potassium signature.** T-x-G-Y-G, the backbone carbonyl cage that
makes a K+ channel selective, is recognisable in *KcsA*, in Kv1.1 and in
GluR0 alike. It needs no reference sequence, so it is the one classification
test in this project that works on a bare sequence from an unannotated
genome.

**The four-repeat filter locus.** Nav, Cav, NALCN and CatSper have identical
domain architectures (hazard **H1**), and the thing that separates them is
one residue contributed by each of the four repeats. Those residues cannot
be found by a regular expression — the surrounding sequence is not conserved
enough — but they can be *projected from a reference*: align the query to
human Nav1.5 with MAFFT and read off whatever aligns to its D372, E898,
K1419 and A1711.

The anchor positions were verified against live UniProt on 2026-08-19, and
the projection was measured on six four-repeat channels:

    CACNA1C → EEEE    CACNA1G → EEDD    NALCN → EEKE
    CACNA1S → EEEE    SCN1A   → DEKA    SCN11A → DEKA

six of six correct, ~1.1 s per query. `CACNA1G` returning **EEDD** rather
than EEEE is not an error: the T-type channels really do differ from Cav1/2
at the filter, which is why `FILTER_CALLS` maps both strings to `cav` and
records the subfamily separately.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..utils.mafft import (MafftUnavailable, align, mafft_available,
                           project_positions)

#: The K+ selectivity filter. Matches TVGYG (Kv, KcsA), TIGYG (Kir, KCNQ),
#: SVGFG (EAG), CIGYG (HCN) and the K2P TxGFG.
K_FILTER_RE = re.compile(r"[TSAC][VILMTAF]G[YFL]G")

#: Motifs that are diagnostic on their own, with what they mean.
SIMPLE_MOTIFS: dict[str, tuple[re.Pattern, str]] = {
    "k_filter": (K_FILTER_RE, "potassium selectivity filter (TxGYG)"),
}


@dataclass(frozen=True)
class FilterAnchor:
    """A reference protein plus the positions its selectivity filter occupies."""
    key: str
    reference_label: str
    reference_uniprot: str
    positions: tuple[int, ...]
    expected: str            # what the reference itself must read
    applies_to: tuple[str, ...]   # catalogue family keys this anchor can call
    note: str = ""


#: The four-repeat locus. Verified against UniProt Q14524 on 2026-08-19:
#: positions 372/898/1419/1711 of human Nav1.5 read D/E/K/A.
FOUR_REPEAT_ANCHOR = FilterAnchor(
    key="four_repeat",
    reference_label="Hs_SCN5A",
    reference_uniprot="Q14524",
    positions=(372, 898, 1419, 1711),
    expected="DEKA",
    applies_to=("nav", "cav", "nalcn", "catsper", "tpc"),
    note="hazard H1; projection measured 6/6 correct",
)

ANCHORS: dict[str, FilterAnchor] = {FOUR_REPEAT_ANCHOR.key: FOUR_REPEAT_ANCHOR}

#: Filter string → (family key, subfamily note). Anything else is a non-call,
#: which is a legitimate answer: TPC1 projects to `--PN` because it has two
#: repeats, not four, and that is information rather than failure.
FILTER_CALLS: dict[str, tuple[str, str]] = {
    "DEKA": ("nav", "voltage-gated sodium channel locus"),
    "EEEE": ("cav", "Cav1/Cav2 high-voltage-activated locus"),
    "EEDD": ("cav", "Cav3 (T-type) locus"),
    "EEKE": ("nalcn", "NALCN leak-channel locus"),
}


@dataclass
class MotifResult:
    """What the motif tier found, and what it could not look for."""
    filter_string: str = ""
    filter_family: str = ""
    filter_note: str = ""
    k_filter_hits: list[tuple[int, str]] = field(default_factory=list)
    available: bool = True
    message: str = ""

    def summary(self) -> str:
        bits = []
        if self.filter_string:
            bits.append(f"four-repeat filter {self.filter_string}"
                        + (f" → {self.filter_family}" if self.filter_family else
                           " (no call)"))
        if self.k_filter_hits:
            pos = ", ".join(f"{p}:{m}" for p, m in self.k_filter_hits[:4])
            bits.append(f"{len(self.k_filter_hits)} K+ filter match(es) [{pos}]")
        if not self.available:
            bits.append(f"unavailable — {self.message}")
        return "; ".join(bits) or "no motif evidence"


def k_filter_hits(sequence: str) -> list[tuple[int, str]]:
    """1-based positions and matched text of every TxGYG-like motif."""
    return [(m.start() + 1, m.group(0))
            for m in K_FILTER_RE.finditer(sequence.upper())]


def filter_signature(query_sequence: str, reference_sequence: str,
                     anchor: FilterAnchor = FOUR_REPEAT_ANCHOR) -> str:
    """Project the anchor positions of the reference onto the query.

    Returns the four projected residues as a string (`-` where the query
    does not cover the position). Raises `MafftUnavailable` if MAFFT is not
    installed — the test does not have a weaker fallback and must not
    pretend otherwise (**D28**).
    """
    if not query_sequence or not reference_sequence:
        return ""
    rows = align([("ref", reference_sequence), ("qry", query_sequence)])
    projected = project_positions(rows["ref"], rows["qry"], list(anchor.positions))
    return "".join(projected[p] for p in anchor.positions)


def verify_anchor(reference_sequence: str,
                  anchor: FilterAnchor = FOUR_REPEAT_ANCHOR) -> bool:
    """Does the reference itself still read what the anchor claims?

    Run before trusting any projection: a UniProt sequence-version bump that
    shifts the numbering would otherwise turn every filter call silently
    wrong. `scripts/s1_benchmark.py` calls this first and refuses to report
    filter results if it fails.
    """
    if len(reference_sequence) < max(anchor.positions):
        return False
    got = "".join(reference_sequence[p - 1] for p in anchor.positions)
    return got == anchor.expected


def scan(query_sequence: str, reference_sequence: str = "",
         anchor: FilterAnchor = FOUR_REPEAT_ANCHOR) -> MotifResult:
    """Run every motif test available for this query."""
    res = MotifResult(k_filter_hits=k_filter_hits(query_sequence))
    if not reference_sequence:
        return res
    if not mafft_available():
        res.available = False
        res.message = "mafft not on PATH; four-repeat filter test skipped"
        return res
    if not verify_anchor(reference_sequence, anchor):
        res.available = False
        res.message = (f"anchor {anchor.key} does not validate against "
                       f"{anchor.reference_label} — refusing to project")
        return res
    try:
        sig = filter_signature(query_sequence, reference_sequence, anchor)
    except (MafftUnavailable, RuntimeError) as exc:
        res.available = False
        res.message = str(exc)[:200]
        return res
    res.filter_string = sig
    call = FILTER_CALLS.get(sig)
    if call:
        res.filter_family, res.filter_note = call
    return res
