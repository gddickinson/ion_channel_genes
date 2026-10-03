"""s11b_lib.py — S11b shared pieces (D53).

* `s11b_dir()` — bulk under `<data root>/trees/s11b/` (raises without the drive).
* `REPEAT_FAMILIES`, `HYPOTHESES` — the repeat test's families and its three
  pairings of the four 4×6TM repeat classes, fixed in D53 (4).
* `ROOT_FAMILIES`, `CONTROL_FAMILIES` — the six S7d families and the five
  positive controls of D53 (5)–(7).
* `iqtree()` — one IQ-TREE 2 call, recorded in `run.json`, skipped when the
  same command already finished on the same input (resumable from `.ckp.gz`).
* `live()` — the dashboard's progress file.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import sha256  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "duplication"
LIVE = ROOT / "results" / "session_live.json"

REPEAT_FAMILIES = {"nav": 4, "cav": 4, "nalcn": 4, "tpc": 2}       # D53 (1)
CLASSES = ["I", "II", "III", "IV"]
HYPOTHESES = {"H13": (("I", "III"), ("II", "IV")),                  # D53 (4)
              "H12": (("I", "II"), ("III", "IV")),
              "H14": (("I", "IV"), ("II", "III"))}
MIN_REPEAT_COVER = 0.9                                              # D53 (2)
AU_ALPHA = 0.05

ROOT_FAMILIES = ["cng", "k2p", "kv_kcnq", "kv_shaker", "kca_slo", "kv_eag"]
CONTROL_FAMILIES = ["enac", "glyr", "ht3", "iglur_nonvertebrate", "nalcn"]
ROOTSTRAP_MIN = 95.0                                                # D53 (6)


def s11b_dir(*parts: str) -> Path:
    d = require_data_root() / "trees" / "s11b"
    for p in parts:
        d = d / p
    d.mkdir(parents=True, exist_ok=True)
    return d


def iqtree(src: Path, prefix: Path, args: list[str], threads: int,
           deps: tuple[Path, ...] = ()) -> dict:
    """Run IQ-TREE once per (input SHA, args, SHA of each file in `deps`);
    the record lives beside it."""
    rec = prefix.parent / f"{prefix.name}.run.json"
    key = {"input_sha256": sha256(src), "args": args}
    if deps:
        key["deps"] = [sha256(p) for p in deps]
    redo = []
    if rec.exists():
        r = json.loads(rec.read_text())
        if {k: r.get(k) for k in key} == key:
            return {**r, "skipped": True}
        redo = ["-redo"]           # a finished run with other inputs: replace it
    cmd = [shutil.which("iqtree2"), "-s", str(src), *args, *redo, "-T", str(threads),
           "--prefix", str(prefix), "--quiet"]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"iqtree2 failed ({prefix.name}): "
                           f"{(p.stderr or p.stdout).strip()[-400:]}")
    r = {**key, "cmd": " ".join(cmd), "seconds": round(time.time() - t0, 1)}
    rec.write_text(json.dumps(r, indent=1))
    return r


def live(steps: list[dict], workers: int = 1) -> None:
    LIVE.write_text(json.dumps({"task": "S11b", "workers": workers,
                                "steps": steps}, indent=1))
