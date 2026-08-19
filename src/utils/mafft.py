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


def percent_identity(a: str, b: str, covered_only: bool = True) -> float:
    """Identity between two aligned strings.

    `covered_only=True` scores over mutually covered columns only, which is
    what the parent project measured to be the right default: full-alignment
    identity dilutes every comparison between proteins of unequal length, and
    ion channels differ in length by an order of magnitude within a single
    superfamily (`MscL` 136 aa, `RYR1` 5,038 aa).
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
