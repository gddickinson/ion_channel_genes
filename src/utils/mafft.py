"""Stdlib-only MAFFT wrapper.

`src/analysis/alignment.py` imports Biopython at module level, which makes it
unavailable in a bare interpreter — and biopython is currently *not*
installed in this project's environment (see `results/toolchain_manifest.txt`).
The classification and phylogeny code cannot depend on that, because
alignment is the discriminating test for several hazards and a test that
silently does not run is worse than no test.

So: this module shells out to MAFFT with nothing but the standard library,
and reports honestly when MAFFT is absent rather than falling back to a
weaker method that would change a classification without saying so
(decision **D28**).
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path


class MafftUnavailable(RuntimeError):
    """Raised when a test needs a real alignment and MAFFT is not installed."""


def mafft_available() -> bool:
    return shutil.which("mafft") is not None


def mafft_version() -> str:
    if not mafft_available():
        return ""
    p = subprocess.run(["mafft", "--version"], capture_output=True, text=True)
    return (p.stderr or p.stdout).strip()


def align(pairs: list[tuple[str, str]], *, auto: bool = True,
          timeout_s: int = 900, extra_args: tuple[str, ...] = ()) -> dict[str, str]:
    """Align `[(label, sequence)]` → `{label: aligned sequence}` (upper case).

    Labels are written to the FASTA verbatim, so keep them free of spaces.
    Raises `MafftUnavailable` rather than degrading silently.
    """
    pairs = [(l, s) for l, s in pairs if s]
    if len(pairs) < 2:
        return {l: s.upper() for l, s in pairs}
    if not mafft_available():
        raise MafftUnavailable(
            "mafft is not on PATH — install it (brew install mafft) or run "
            "this task in the environment recorded in "
            "results/toolchain_manifest.txt")

    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "in.fasta"
        fa.write_text("".join(f">{l}\n{s}\n" for l, s in pairs))
        cmd = ["mafft", "--quiet"]
        cmd += list(extra_args) or (["--auto"] if auto else [])
        cmd.append(str(fa))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout_s)
    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed (rc={proc.returncode}): "
                           f"{proc.stderr.strip()[:400]}")
    return parse_fasta(proc.stdout)


def parse_fasta(text: str) -> dict[str, str]:
    rows: dict[str, str] = {}
    label, chunk = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if label is not None:
                rows[label] = "".join(chunk).upper()
            label, chunk = line[1:].strip().split()[0], []
        elif label is not None:
            chunk.append(line.strip())
    if label is not None:
        rows[label] = "".join(chunk).upper()
    return rows


def project_positions(aligned_ref: str, aligned_query: str,
                      positions: list[int]) -> dict[int, str]:
    """Map 1-based reference positions onto the residue the query aligns there.

    A position the query does not cover comes back as `-`, which is a real
    answer (the query has no equivalent residue) and not a failure.
    """
    want = set(positions)
    out: dict[int, str] = {}
    ri = 0
    for cr, cq in zip(aligned_ref, aligned_query):
        if cr != "-":
            ri += 1
            if ri in want:
                out[ri] = cq
    for p in positions:
        out.setdefault(p, "-")
    return out


def alignment_stats(a: str, b: str) -> tuple[float, float, float, int]:
    """`(covered identity, coverage of the shorter seq, of the longer, columns)`.

    **Coverage is not optional.** Covered-only identity — scoring only the
    columns where both sequences have a residue — is the right metric for
    comparing proteins of unequal length, and on its own it is dangerous:
    the fewer columns two sequences share, the more the score is dominated
    by whichever residues the aligner happened to match.

    Measured in S1: human connexin-26 (226 aa) scored **61.3 % identity to
    ryanodine receptor 2 (4,967 aa)** — higher than its 49.8 % to its own
    relative connexin-43 — and GluN1 scored 50.9 % to the SARS-CoV-2
    envelope protein. The mechanism is cherry-picking: an aligner placing
    226 residues inside 4,967 chooses the 226 best-matching positions, and
    covered-only identity then scores exactly those. The parent project
    never met this because its sequences were one family and all the same
    size; a catalogue with a 40× length range meets it immediately.

    Coverage of the *shorter* sequence does not catch it (0.996 in the
    connexin case). **Coverage of the longer one does** — 225 covered
    columns out of 4,967 is 4.5 %, against 59 % for the true relative — so
    both are returned and `src/classify/reference.py:MIN_COVERAGE` enforces
    the floor on the longer.
    """
    same = covered = 0
    len_a = len_b = 0
    for ca, cb in zip(a, b):
        ga, gb = ca == "-", cb == "-"
        if not ga:
            len_a += 1
        if not gb:
            len_b += 1
        if ga or gb:
            continue
        covered += 1
        if ca == cb:
            same += 1
    shorter = min(len_a, len_b) or 1
    longer = max(len_a, len_b) or 1
    ident = (same / covered) if covered else 0.0
    return ident, covered / shorter, covered / longer, covered


def percent_identity(a: str, b: str, covered_only: bool = True) -> float:
    """Identity between two aligned strings.

    `covered_only=True` scores over mutually covered columns only. Prefer
    `alignment_stats()`, which also returns the coverage the number rests on
    — see the warning there.
    """
    same = cols = 0
    for ca, cb in zip(a, b):
        if covered_only and (ca == "-" or cb == "-"):
            continue
        if ca == "-" and cb == "-":
            continue
        cols += 1
        if ca == cb:
            same += 1
    return (same / cols) if cols else 0.0
