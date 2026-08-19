"""Extracting the comparable part of a channel — the pore module.

A tier-2 tree cannot be built from full-length sequences. Inside the P-loop
superfamily alone the members run from 400 aa (Kir) to 5,038 aa (RYR1), and
almost all of that difference is cytosolic machinery that has no counterpart
in the other families: ankyrin repeats, RCK rings, IQ motifs, a 3,000-residue
solenoid. Aligning them end to end produces an alignment that is 90 % gap and
a tree that measures length.

What *is* comparable is the pore module: the last two transmembrane helices
and the re-entrant loop between them, ~120 residues in a P-loop channel.
Every family in the superfamily has exactly one per subunit (four per chain
in Nav/Cav), it carries the selectivity filter, and it is the part whose
homology is not in doubt.

Two ways to find it, in order of preference:

1. **From the domain envelope.** InterPro gives the `PF00520`/`PF07885`
   boundaries; the pore module is the C-terminal portion of that envelope.
   Exact, and available for anything with an InterPro record.
2. **By projection from a reference.** Align to a family exemplar whose
   module boundaries are known and take whatever aligns inside them. Works
   on an unannotated gene model, which is the case that matters for a census.

Both are recorded per sequence in the alignment manifest, because a tree
built from modules found two different ways is a tree with a confounder in
it, and the manifest is what lets S6 check.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..utils.mafft import align, project_positions

#: Fraction of a pore-domain envelope, measured from its C-terminal end,
#: that holds the pore helix + selectivity filter + inner helix.
MODULE_TAIL_FRACTION = 0.45

#: Minimum residues for an extracted module to be worth aligning.
MIN_MODULE_AA = 60


@dataclass
class PoreModule:
    label: str
    sequence: str
    start: int              # 1-based, in the parent sequence
    end: int
    method: str             # "envelope" | "projected" | "whole"
    source_accession: str = ""   # the Pfam envelope or reference used
    note: str = ""


@dataclass
class ModuleSet:
    """Modules for one alignment unit, plus the provenance of each."""
    unit: str
    modules: list[PoreModule] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)

    def labelled(self) -> list[tuple[str, str]]:
        return [(m.label, m.sequence) for m in self.modules]

    def manifest_rows(self) -> list[tuple]:
        return [(self.unit, m.label, m.method, m.source_accession,
                 m.start, m.end, len(m.sequence), m.note)
                for m in self.modules]

    def summary(self) -> str:
        methods: dict[str, int] = {}
        for m in self.modules:
            methods[m.method] = methods.get(m.method, 0) + 1
        parts = ", ".join(f"{n}×{k}" for k, n in sorted(methods.items()))
        s = f"{self.unit}: {len(self.modules)} module(s) [{parts}]"
        if self.skipped:
            s += f"; {len(self.skipped)} skipped"
        return s


def module_from_envelope(label: str, sequence: str, start: int, end: int,
                         accession: str = "",
                         tail_fraction: float = MODULE_TAIL_FRACTION
                         ) -> PoreModule | None:
    """Take the C-terminal `tail_fraction` of a pore-domain envelope."""
    if not sequence or end <= start:
        return None
    env_len = end - start + 1
    mod_start = start + int(env_len * (1 - tail_fraction))
    seg = sequence[mod_start - 1:end]
    if len(seg) < MIN_MODULE_AA:
        return None
    return PoreModule(label, seg, mod_start, end, "envelope", accession)


def modules_from_envelopes(label: str, sequence: str,
                           envelopes: list[tuple[int, int]],
                           accession: str = "") -> list[PoreModule]:
    """One module per envelope — four for a Nav or Cav chain, one for a Kv."""
    out = []
    for i, (a, b) in enumerate(sorted(envelopes), 1):
        m = module_from_envelope(
            f"{label}_r{i}" if len(envelopes) > 1 else label,
            sequence, a, b, accession)
        if m:
            out.append(m)
    return out


def module_by_projection(label: str, sequence: str,
                         reference_label: str, reference_sequence: str,
                         reference_span: tuple[int, int]) -> PoreModule | None:
    """Align to a reference and take whatever falls inside its module span."""
    if not sequence or not reference_sequence:
        return None
    rows = align([("ref", reference_sequence), ("qry", sequence)])
    a, b = reference_span
    mapped = project_positions(rows["ref"], rows["qry"], [a, b])
    # Walk the alignment to find the query residue indices at the span ends.
    ri = qi = 0
    q_start = q_end = 0
    for cr, cq in zip(rows["ref"], rows["qry"]):
        if cr != "-":
            ri += 1
        if cq != "-":
            qi += 1
        if ri == a and q_start == 0 and cq != "-":
            q_start = qi
        if ri == b and cq != "-":
            q_end = qi
    if not q_start or not q_end or q_end - q_start + 1 < MIN_MODULE_AA:
        return None
    seg = sequence[q_start - 1:q_end]
    return PoreModule(label, seg, q_start, q_end, "projected",
                      reference_label,
                      note=f"span {a}-{b} of {reference_label}; "
                           f"ends mapped to {mapped[a]}/{mapped[b]}")


def whole_sequence_module(label: str, sequence: str) -> PoreModule:
    """Tier-1 fallback: within a family, full-length alignment is fine."""
    return PoreModule(label, sequence, 1, len(sequence), "whole")
