"""The reference tier — best-hit identity against the catalogue's exemplars, with a margin.

Domain architecture answers "what kind of thing is this?" and stops. The
reference tier answers "which family is it closest to?", which is the only
route available for the families whose architectures are identical: Nav
versus Cav before the motif test runs, ANO1 versus ANO6 at all, the eight
Cys-loop families, the four iGluR families.

Two rules govern it, both inherited from the parent project and both paid
for by a benchmark failure there.

**D7 — a best hit is not a call.** Two families scoring within 10 % of each
other is an ambiguity, not a winner. The margin is reported with every call
so a downstream reader can see how close it was.

**Identity is scored over mutually covered columns.** Ion channels in one
superfamily differ in length by an order of magnitude — 136 aa for MscL,
5,038 aa for RYR1 — and full-alignment identity divides by the gaps, which
pushes real relationships below every threshold. The parent project measured
RyR-vs-ITPR identity at 0.105 full-alignment and 0.249 covered-only; the
first number is noise and the second is a signal.

MAFFT does the alignments, and a 3-mer containment prefilter keeps the
number of alignments down: scoring one query against 160-odd exemplars at a
second each would be three minutes per protein, which is not a census.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

from ..catalogue import CATALOGUE, exemplars
from ..utils.mafft import MafftUnavailable, align, alignment_stats

#: D7: the minimum gap between the best family and the next-best *different*
#: family for the reference tier to make a call.
DEFAULT_MARGIN = 0.10

#: How many prefiltered candidates get a real alignment.
DEFAULT_TOP_N = 12

#: Minimum fraction of the *longer* sequence the alignment must cover for the
#: identity to mean anything. Measured in S1: without it, connexin-26 (226 aa)
#: scored 61.3 % identity to ryanodine receptor 2 (4,967 aa) — higher than to
#: its own relative connexin-43 — because an aligner placing 226 residues
#: inside 4,967 picks the 226 best-matching positions and covered-only
#: identity then scores exactly those. Coverage of the shorter sequence does
#: not catch it (0.996); coverage of the longer does (0.045 vs 0.589).
MIN_COVERAGE = 0.30


@dataclass
class ReferenceHit:
    label: str
    family: str
    uniprot: str
    identity: float
    prefilter_score: float = 0.0
    coverage: float = 0.0        # of the longer sequence — the guard below
    covered_cols: int = 0


@dataclass
class ReferenceResult:
    family: str = ""
    margin: float = 0.0
    hits: list[ReferenceHit] = field(default_factory=list)
    ambiguous: list[str] = field(default_factory=list)
    rejected_low_coverage: list[str] = field(default_factory=list)
    available: bool = True
    message: str = ""

    def best(self) -> ReferenceHit | None:
        return self.hits[0] if self.hits else None

    def summary(self) -> str:
        if not self.available:
            return f"reference tier unavailable — {self.message}"
        if not self.hits:
            return "no reference hit"
        b = self.hits[0]
        head = (f"nearest {b.label} ({b.family}) at {b.identity:.1%} over "
                f"{b.covered_cols} cols ({b.coverage:.0%} of the longer), "
                f"margin {self.margin:+.3f}")
        if self.ambiguous:
            head += " — ambiguous between " + ", ".join(self.ambiguous)
        if self.rejected_low_coverage:
            head += (f"; {len(self.rejected_low_coverage)} hit(s) dropped "
                     f"below {MIN_COVERAGE:.0%} coverage")
        return head


def kmer_profile(seq: str, k: int = 3) -> set[str]:
    s = seq.upper()
    return {s[i:i + k] for i in range(max(0, len(s) - k + 1))}


def kmer_containment(query: set[str], ref: set[str]) -> float:
    """Fraction of the smaller profile shared — length-robust, unlike Jaccard."""
    if not query or not ref:
        return 0.0
    return len(query & ref) / min(len(query), len(ref))


class ReferenceSet:
    """The exemplar sequences the reference tier scores against.

    Sequences live in a FASTA cache (written by
    `scripts/s0_catalogue_verify.py`) so a classification run is
    reproducible offline and does not re-fetch UniProt once per query.
    """

    def __init__(self, sequences: dict[str, str],
                 family_of: dict[str, str],
                 accession_of: dict[str, str] | None = None) -> None:
        self.sequences = {k: v for k, v in sequences.items() if v}
        self.family_of = family_of
        self.accession_of = accession_of or {}
        self._profiles = {l: kmer_profile(s) for l, s in self.sequences.items()}

    # -- construction ----------------------------------------------------
    @classmethod
    def from_fasta(cls, path: Path) -> "ReferenceSet":
        """Load `>label|family|accession` records written by the S0 script."""
        seqs: dict[str, str] = {}
        fam: dict[str, str] = {}
        acc: dict[str, str] = {}
        label = ""
        chunk: list[str] = []
        for line in Path(path).read_text().splitlines():
            if line.startswith(">"):
                if label:
                    seqs[label] = "".join(chunk)
                parts = line[1:].strip().split("|")
                label = parts[0]
                fam[label] = parts[1] if len(parts) > 1 else ""
                acc[label] = parts[2] if len(parts) > 2 else ""
                chunk = []
            elif label:
                chunk.append(line.strip())
        if label:
            seqs[label] = "".join(chunk)
        return cls(seqs, fam, acc)

    def to_fasta(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(
            f">{lbl}|{self.family_of.get(lbl, '')}|{self.accession_of.get(lbl, '')}\n"
            f"{seq}\n" for lbl, seq in sorted(self.sequences.items())))
        return path

    @classmethod
    def fetch(cls, census_only: bool = True, sleep_s: float = 0.1,
              on_progress=None) -> "ReferenceSet":
        """Fetch every exemplar that already carries a UniProt accession."""
        import time
        seqs: dict[str, str] = {}
        fam: dict[str, str] = {}
        acc: dict[str, str] = {}
        for family_key, ex in exemplars(census_only):
            if not ex.uniprot:
                continue
            if on_progress:
                on_progress(f"{ex.label} ({ex.uniprot})")
            s = fetch_uniprot_sequence(ex.uniprot)
            if s:
                seqs[ex.label] = s
                fam[ex.label] = family_key
                acc[ex.label] = ex.uniprot
            time.sleep(sleep_s)
        return cls(seqs, fam, acc)

    # -- scoring ---------------------------------------------------------
    def labels_for(self, family_key: str) -> list[str]:
        return [l for l, f in self.family_of.items() if f == family_key]

    def prefilter(self, query_seq: str, top_n: int = DEFAULT_TOP_N,
                  labels: set[str] | None = None,
                  min_coverage: float = MIN_COVERAGE) -> list[tuple[str, float]]:
        """Cheap candidate shortlist, with the length filter D30 implies.

        `kmer_containment` divides by the smaller profile, so a 5,000-residue
        reference containing most of a 400-residue query's 3-mers scores near
        1.0 — the prefilter is biased towards the longest references in the
        panel. Those are exactly the hits the coverage floor then rejects, so
        without a length filter the tier spends most of its alignments on
        candidates it has already decided to discard. Measured: the S1 panel
        went from ~9 s to ~90 s per protein once the search was scoped to a
        superfamily containing the ryanodine receptors.

        The filter is not a heuristic — it is the coverage floor stated in
        advance. Coverage of the longer sequence can never exceed
        `len(shorter) / len(longer)`, so any reference outside that ratio
        cannot pass `MIN_COVERAGE` however it aligns.
        """
        qp = kmer_profile(query_seq)
        lq = max(1, len(query_seq))
        items = (self._profiles.items() if labels is None
                 else ((l, p) for l, p in self._profiles.items() if l in labels))
        scored = []
        for label, prof in items:
            lr = max(1, len(self.sequences[label]))
            if min(lq, lr) / max(lq, lr) < min_coverage:
                continue
            scored.append((label, kmer_containment(qp, prof)))
        scored.sort(key=lambda t: -t[1])
        return scored[:top_n]

    def best_match(self, query_seq: str, top_n: int = DEFAULT_TOP_N,
                   margin: float = DEFAULT_MARGIN,
                   restrict_to: tuple[str, ...] = (),
                   exclude_accessions: tuple[str, ...] = (),
                   min_coverage: float = MIN_COVERAGE) -> ReferenceResult:
        """Nearest exemplar, with the D7 margin.

        `exclude_accessions` is the benchmark's leave-one-out: 46 of the S1
        panel's 97 proteins **are** catalogue exemplars, so without it the
        reference tier would score each of them against itself at 100 %
        identity and the benchmark would measure nothing. Excluding the
        query's own accession makes the tier answer the question the census
        actually asks — *what is this nearest to, other than itself?*
        """
        if not query_seq or not self.sequences:
            return ReferenceResult(available=False, message="no reference sequences")
        drop = {l for l, a in self.accession_of.items()
                if a in set(exclude_accessions)} if exclude_accessions else set()
        # When the architecture tier came back ambiguous it names the families
        # in play, and the reference tier's whole job is to break that tie.
        # Prefiltering over *only* those families is what makes the tie-break
        # actually happen: a global top-N can easily contain none of them, and
        # then the ambiguity survives a tier that was asked to resolve it.
        keep: set[str] | None = None
        if restrict_to:
            keep = {l for l in self.sequences
                    if self.family_of.get(l, "") in restrict_to}
            if not keep:
                keep = None
        if drop:
            keep = (keep or set(self.sequences)) - drop
        candidates = self.prefilter(query_seq, top_n, keep, min_coverage)
        hits: list[ReferenceHit] = []
        dropped: list[str] = []
        for label, pre in candidates:
            try:
                rows = align([("q", query_seq), (label, self.sequences[label])])
            except (MafftUnavailable, RuntimeError) as exc:
                return ReferenceResult(available=False, message=str(exc)[:200])
            ident, _cov_short, cov_long, n_cols = alignment_stats(
                rows["q"], rows[label])
            if cov_long < min_coverage:
                dropped.append(f"{label}({cov_long:.0%})")
                continue
            hits.append(ReferenceHit(label, self.family_of.get(label, ""),
                                     self.accession_of.get(label, ""),
                                     ident, pre, cov_long, n_cols))
        hits.sort(key=lambda h: -h.identity)
        res = ReferenceResult(hits=hits, rejected_low_coverage=dropped)
        if not hits:
            return res
        best = hits[0]
        runner = next((h for h in hits[1:] if h.family != best.family), None)
        res.margin = best.identity - (runner.identity if runner else 0.0)
        if runner is None or res.margin >= margin:
            res.family = best.family
        else:
            res.ambiguous = sorted({best.family, runner.family})
        return res


def fetch_uniprot_sequence(accession: str, timeout_s: int = 30) -> str:
    """One sequence by accession, stdlib only (no `requests` dependency)."""
    if not accession:
        return ""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=timeout_s) as r:
            text = r.read().decode()
    except Exception:
        return ""
    return "".join(l.strip() for l in text.splitlines() if not l.startswith(">"))


def family_label_map() -> dict[str, str]:
    """`{exemplar label: family name}` — for report tables."""
    return {e.label: CATALOGUE[k].name for k, e in exemplars(False)}


def write_result(path: Path, result: ReferenceResult) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "family": result.family,
        "margin": result.margin,
        "ambiguous": result.ambiguous,
        "hits": [h.__dict__ for h in result.hits],
        "available": result.available,
        "message": result.message,
    }, indent=2))
    return path
