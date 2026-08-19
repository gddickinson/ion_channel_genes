"""The phylogeny is a forest, not a tree — and the code enforces that.

**Decision D27.** There is no alignment that contains a nicotinic receptor
and a Kv channel. The ion channels are not one homologous group: they are
roughly twenty-five independent origins, and the only honest global object
is a forest of trees plus a network of fold relationships between them.

Three tiers:

* **Tier 1 — within a family.** Full-length alignment, rooted on the sister
  family inside the same superfamily. Conventional, uncontroversial, and
  where nearly every biologically useful result lives.
* **Tier 2 — within a superfamily.** Pore-module alignment only
  (`src/phylo/modules.py`), rooted on the superfamily's declared outgroup.
  The tree is labelled a *pore-module tree* in every output, because that is
  what it is.
* **Tier 3 — between superfamilies.** Refused. `build_tier2` raises
  `NotAlignable` for any superfamily the catalogue marks
  `alignable=False`, and there is no `build_tier3`: cross-superfamily
  comparison goes to `src/phylo/network.py`, which produces a fold-similarity
  network and never a phylogram.

The refusal is the point. Every published "tree of ion channels" that spans
non-homologous superfamilies is measuring alignment artefacts, and the only
way not to produce one by accident is to make it impossible to ask for.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

from ..catalogue import CATALOGUE, SUPERFAMILIES
from ..utils.mafft import MafftUnavailable, align, mafft_version


class NotAlignable(RuntimeError):
    """Raised when a tree is requested across a non-homologous boundary."""


@dataclass
class TreeRun:
    """Everything needed to reproduce and audit one tree."""
    unit: str
    tier: int
    n_sequences: int
    alignment_path: str = ""
    trimmed_path: str = ""
    tree_path: str = ""
    newick: str = ""
    outgroup: tuple[str, ...] = ()
    module_based: bool = False
    tools: dict[str, str] = field(default_factory=dict)
    commands: list[str] = field(default_factory=list)
    elapsed_s: float = 0.0
    note: str = ""

    def label(self) -> str:
        kind = "pore-module tree" if self.module_based else "protein tree"
        return f"tier {self.tier} {kind} of {self.unit}"

    def to_dict(self) -> dict:
        return {**self.__dict__, "label": self.label()}


def _tool(name: str) -> str:
    return shutil.which(name) or ""


def toolchain() -> dict[str, str]:
    """What is actually available, recorded into every run."""
    out = {"mafft": mafft_version()}
    for t in ("trimal", "iqtree2", "iqtree", "fasttree"):
        p = _tool(t)
        if p:
            out[t] = p
    return out


def _run(cmd: list[str], run: TreeRun, timeout_s: int = 7200) -> subprocess.CompletedProcess:
    run.commands.append(" ".join(cmd))
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    if p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed (rc={p.returncode}): "
                           f"{p.stderr.strip()[:400]}")
    return p


def build_tree(unit: str, tier: int, sequences: list[tuple[str, str]],
               out_dir: Path, outgroup: tuple[str, ...] = (),
               module_based: bool = False, model: str = "MFP",
               bootstrap: int = 1000, trim: bool = True,
               linsi: bool = False) -> TreeRun:
    """MAFFT → trimAl → IQ-TREE 2, with every command recorded.

    Falls back to writing the alignment and stopping if IQ-TREE is absent —
    a missing tool produces a missing tree and a note, never a silently
    weaker method (**D28**).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run = TreeRun(unit=unit, tier=tier, n_sequences=len(sequences),
                  outgroup=outgroup, module_based=module_based,
                  tools=toolchain())
    t0 = time.time()

    if len(sequences) < 4:
        run.note = f"only {len(sequences)} sequences — too few for a tree"
        run.elapsed_s = time.time() - t0
        return run

    extra = ("--maxiterate", "1000", "--localpair") if linsi else ()
    try:
        aligned = align(sequences, auto=not linsi, extra_args=extra)
    except MafftUnavailable as exc:
        run.note = str(exc)
        run.elapsed_s = time.time() - t0
        return run
    run.commands.append(f"mafft {'--localpair --maxiterate 1000' if linsi else '--auto'}")
    aln = out_dir / f"{unit}.aln.fasta"
    aln.write_text("".join(f">{l}\n{s}\n" for l, s in aligned.items()))
    run.alignment_path = str(aln)

    current = aln
    if trim and _tool("trimal"):
        trimmed = out_dir / f"{unit}.trimal.fasta"
        _run([_tool("trimal"), "-in", str(aln), "-out", str(trimmed),
              "-automated1"], run)
        run.trimmed_path = str(trimmed)
        current = trimmed

    iq = _tool("iqtree2") or _tool("iqtree")
    if not iq:
        run.note = "iqtree2 not on PATH — alignment written, tree not built"
        run.elapsed_s = time.time() - t0
        return run

    prefix = out_dir / unit
    cmd = [iq, "-s", str(current), "-m", model, "-B", str(bootstrap),
           "--prefix", str(prefix), "-T", "AUTO", "--quiet", "-redo"]
    if outgroup:
        cmd += ["-o", ",".join(outgroup)]
    _run(cmd, run)
    tree_file = prefix.with_suffix(".treefile")
    if tree_file.exists():
        run.tree_path = str(tree_file)
        run.newick = tree_file.read_text().strip()
    run.elapsed_s = time.time() - t0
    return run


def build_tier1(family_key: str, sequences: list[tuple[str, str]],
                out_dir: Path, **kw) -> TreeRun:
    """A within-family tree, rooted on the family's declared sister group."""
    fam = CATALOGUE[family_key]
    sf = SUPERFAMILIES.get(fam.superfamily)
    outgroup = kw.pop("outgroup", ())
    if not outgroup and sf and sf.root_with:
        outgroup = tuple(l for l, _ in
                         [(e.label, e) for k in sf.root_with
                          for e in CATALOGUE[k].exemplars]
                         if any(l == s for s, _ in sequences))
    return build_tree(family_key, 1, sequences, out_dir,
                      outgroup=outgroup, module_based=False, **kw)


def needs_modules(superfamily_key: str) -> bool:
    """Does this superfamily require pore-module extraction for tier 2?

    A few superfamilies are alignable end to end — the Cys-loop receptors
    are the clearest case, and the catalogue records that in
    `Superfamily.anchor_module`. Everywhere else full-length alignment across
    the superfamily is meaningless (Kir 400 aa against RYR1 5,038 aa), and
    the tree must be built on the pore module and labelled as one.
    """
    sf = SUPERFAMILIES[superfamily_key]
    return "full-length" not in (sf.anchor_module or "").lower()


def build_tier2(superfamily_key: str, sequences: list[tuple[str, str]],
                out_dir: Path, module_based: bool | None = None,
                **kw) -> TreeRun:
    """A within-superfamily tree. Refuses non-alignable superfamilies (D27).

    `module_based` defaults to what the catalogue says the superfamily needs.
    It is a *label as well as a method*: a tier-2 run marked `module_based`
    is reported as a pore-module tree everywhere it appears, because calling
    a pore-module tree a protein tree is the misreading D27 exists to
    prevent.
    """
    sf = SUPERFAMILIES[superfamily_key]
    if not sf.alignable:
        raise NotAlignable(
            f"{superfamily_key} ({sf.name}) is marked alignable=False: its "
            f"member families share a fold with no detectable sequence "
            f"homology. A tier-2 tree here would be an alignment artefact — "
            f"use src/phylo/network.py instead (D27). {sf.notes}")
    outgroup = kw.pop("outgroup", ())
    if not outgroup and sf.root_with:
        labels = {s for s, _ in sequences}
        outgroup = tuple(e.label for k in sf.root_with
                         for e in CATALOGUE[k].exemplars if e.label in labels)
    if module_based is None:
        module_based = needs_modules(superfamily_key)
    return build_tree(superfamily_key, 2, sequences, out_dir,
                      outgroup=outgroup, module_based=module_based, **kw)


def write_run(out_dir: Path, run: TreeRun) -> Path:
    p = Path(out_dir) / f"{run.unit}.run.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(run.to_dict(), indent=2))
    return p


def alignable_superfamilies() -> list[str]:
    return [k for k, sf in SUPERFAMILIES.items() if sf.alignable]


def refused_superfamilies() -> list[tuple[str, str]]:
    """The ones D27 forbids treeing, with the reason — reported, not hidden."""
    return [(k, sf.notes or "member families share a fold, not a sequence")
            for k, sf in SUPERFAMILIES.items() if not sf.alignable]
