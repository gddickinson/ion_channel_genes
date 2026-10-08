"""s13_codeml.py — one codeml or HyPhy job: control file, resumable run, lnL.

Ported from `../ip3r_genes/scripts/s9_codeml_lib.py`. Each job runs in its
own directory under `<data root>/selection/s13/jobs/<name>/` because codeml
writes fixed filenames (`mlc`, `rst`, `2ML.dS`); a job directory is
**claimed** with a pid lock while it runs, since two drivers sharing one
interleave their output silently and `mlc` still parses. A job whose output
already parses is reused, never re-run, so an interrupted suite resumes.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s13_lib import tool_bin  # noqa: E402

CTL_DEFAULTS = {
    "noisy": "3", "verbose": "1", "runmode": "0", "seqtype": "1",
    "CodonFreq": "2", "clock": "0", "aaDist": "0", "icode": "0",
    "model": "0", "NSsites": "0", "fix_kappa": "0", "kappa": "2",
    "fix_omega": "0", "omega": "0.4", "fix_alpha": "1", "alpha": "0",
    "Malpha": "0", "ncatG": "10", "getSE": "0", "RateAncestor": "0",
    "Small_Diff": ".5e-6", "cleandata": "0", "method": "0",
}
LNL_RE = re.compile(r"lnL\(ntime:\s*(\d+)\s+np:\s*(\d+)\)[:\s]*(-?[\d.]+)")


@dataclass
class Job:
    name: str
    tool: str                    # "codeml" | "fel" | "relax"
    seqfile: Path
    treefile: Path
    workdir: Path
    settings: dict = field(default_factory=dict)
    group: str = ""              # the set or family the job belongs to
    cost: float = 1.0            # tips × codons, for scheduling the big ones first

    def ctl_text(self) -> str:
        ctl = dict(CTL_DEFAULTS)
        ctl.update({k: str(v) for k, v in self.settings.items()})
        head = [f"seqfile = {self.seqfile}", "treefile = tree.paml", "outfile = mlc"]
        return "\n".join(head + [f"{k} = {v}" for k, v in ctl.items()]) + "\n"

    def output(self) -> Path:
        return self.workdir / ("mlc" if self.tool == "codeml" else f"{self.tool}.json")

    def done(self) -> bool:
        out = self.output()
        if not out.exists() or not out.stat().st_size:
            return False
        if self.tool != "codeml":
            return True
        if str(self.settings.get("runmode", "0")) == "-2":
            return (self.workdir / "2ML.dS").exists()
        return parse_lnl(out) is not None


def parse_lnl(mlc: Path) -> tuple[float, int] | None:
    m = LNL_RE.search(mlc.read_text()) if mlc.exists() else None
    return (float(m.group(3)), int(m.group(2))) if m else None


def owner_alive(lock: Path) -> int:
    """The live pid holding `lock`, or 0 (missing, malformed or stale)."""
    try:
        pid = int(lock.read_text().split()[0])
        os.kill(pid, 0)
        return pid
    except (OSError, ValueError, IndexError):
        return 0


def run(job: Job, threads: int = 1, timeout_s: int = 7 * 86_400) -> str:
    """Run unless done or claimed; → 'done' | 'ran' | 'claimed' | 'failed'."""
    job.workdir.mkdir(parents=True, exist_ok=True)
    if job.done():
        return "done"
    lock = job.workdir / "RUNNING"
    if owner_alive(lock):
        return "claimed"
    lock.write_text(f"{os.getpid()}\n")
    try:
        if job.tool == "codeml":
            # A headed tree file ("ntaxa ntrees"): without it codeml finishes
            # the fit, then exits non-zero looking for a second tree.
            nwk = job.treefile.read_text().strip()
            (job.workdir / "tree.paml").write_text(f" {nwk.count(',') + 1} 1\n{nwk}\n")
            (job.workdir / "codeml.ctl").write_text(job.ctl_text())
            cmd = [tool_bin("codeml"), "codeml.ctl"]
        else:
            cmd = [tool_bin("hyphy"), f"CPU={threads}", job.tool,
                   "--alignment", str(job.seqfile), "--tree", str(job.treefile),
                   "--code", "Universal", "--output", str(job.output())]
            for k, v in job.settings.items():
                cmd += [f"--{k}", str(v)]
        (job.workdir / "cmd.txt").write_text(" ".join(cmd) + "\n")
        with open(job.workdir / "stdout.txt", "w") as out, \
                open(job.workdir / "stderr.txt", "w") as err:
            p = subprocess.run(cmd, cwd=job.workdir, stdout=out, stderr=err,
                               stdin=subprocess.DEVNULL, timeout=timeout_s)
        return "ran" if p.returncode == 0 and job.done() else "failed"
    except subprocess.TimeoutExpired:
        (job.workdir / "TIMEOUT").write_text(f"timeout after {timeout_s} s\n")
        return "failed"
    finally:
        lock.unlink(missing_ok=True)
