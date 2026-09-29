"""s5_genome_io.py — genome download and FASTA plumbing for the S5 sweep.

Ported from `../ip3r_genes/scripts/fetch_genomes.py` and `s5_genome_io.py`
(method only). Resumable: a completed genome leaves a `.done` marker; nothing
is installed until every file in NCBI's `md5sum.txt` has verified, because a
truncated genome searched without complaint reports an absence.

Genomes live under `<data root>/genomes/<accession>/`; the sweep's own output
under `<data root>/genomes/s5/`. `require_data_root()` is called before
anything is written.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import sys
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.utils.data_root import require_data_root  # noqa: E402

DATASETS = "/opt/anaconda3/envs/piezo1/bin/datasets"     # D18
_N_RUN = re.compile(r"[Nn]+")
NGAP_RUN = 100            # an N-run this long counts as an assembly gap


def genomes_root() -> Path:
    p = require_data_root() / "genomes"
    p.mkdir(parents=True, exist_ok=True)
    return p


def genome_dir(accession: str) -> Path:
    return genomes_root() / accession


def is_done(accession: str) -> bool:
    return (genome_dir(accession) / ".done").exists()


def fna_path(accession: str) -> Path | None:
    hits = sorted(genome_dir(accession).glob("*_genomic.fna"))
    return hits[0] if hits else None


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_extract(zf: zipfile.ZipFile, dest: Path) -> None:
    for name in zf.namelist():
        if name.startswith(("/", "..")) or ".." in Path(name).parts:
            raise RuntimeError(f"unsafe path in zip: {name}")
    zf.extractall(dest)


def _verify_md5s(staging: Path) -> int:
    md5file = staging / "md5sum.txt"
    if not md5file.exists():
        raise RuntimeError("zip is missing md5sum.txt")
    n = 0
    for line in md5file.read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        want, rel = parts
        target = staging / rel
        if not target.exists():
            raise RuntimeError(f"md5sum.txt names a missing file: {rel}")
        if _md5(target) != want:
            raise RuntimeError(f"md5 mismatch for {rel}")
        n += 1
    if n == 0:
        raise RuntimeError("md5sum.txt contained no checksums")
    return n


def fetch_one(accession: str, quiet: bool = False, attempts: int = 4) -> Path:
    """Download, verify and install one genome; returns its directory.

    Download → extract → verify retried together: a dropped stream surfaces
    as an md5 failure, not as a download error.
    """
    dest = genome_dir(accession)
    if is_done(accession):
        return dest
    zips = genomes_root() / "_zips"
    staging = genomes_root() / "_staging" / accession
    zips.mkdir(parents=True, exist_ok=True)
    zip_path = zips / f"{accession}.zip"
    last_err, n_verified = "", 0
    for attempt in range(1, attempts + 1):
        if staging.exists():
            shutil.rmtree(staging)
        zip_path.unlink(missing_ok=True)
        proc = subprocess.run(
            [DATASETS, "download", "genome", "accession", accession,
             "--include", "genome,seq-report", "--filename", str(zip_path),
             "--no-progressbar"], capture_output=True, text=True)
        if proc.returncode != 0:
            last_err = f"datasets download failed: {proc.stderr.strip()[:200]}"
        else:
            try:
                with zipfile.ZipFile(zip_path) as zf:
                    _safe_extract(zf, staging)
                n_verified = _verify_md5s(staging)
                break
            except (zipfile.BadZipFile, RuntimeError) as exc:
                last_err = str(exc)[:200]
        if attempt == attempts:
            raise RuntimeError(f"{accession}: {last_err} (after {attempts} attempts)")
        if not quiet:
            print(f"  retry {attempt}/{attempts - 1} for {accession}: {last_err}")
        time.sleep(5 * attempt)

    payload = staging / "ncbi_dataset" / "data"
    src_dir = payload / accession
    if not src_dir.is_dir():
        raise RuntimeError(f"unexpected zip layout for {accession}")
    if dest.exists():
        shutil.rmtree(dest)
    shutil.move(str(src_dir), str(dest))
    for extra in payload.glob("*.jsonl"):
        shutil.move(str(extra), str(dest / extra.name))
    # Keep checksums of exactly the installed files, so --verify can treat a
    # missing name as a failure rather than skipping it.
    kept = []
    for line in (staging / "md5sum.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) == 2 and (dest / Path(parts[1]).name).exists():
            kept.append(f"{parts[0]}  {Path(parts[1]).name}")
    (dest / "md5sum.txt").write_text("\n".join(kept) + "\n")
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    (dest / ".done").write_text(f"accession={accession}\nfetched={stamp}\n"
                                f"md5_verified_files={n_verified}\n")
    zip_path.unlink(missing_ok=True)
    shutil.rmtree(staging, ignore_errors=True)
    if not quiet:
        gb = sum(p.stat().st_size for p in dest.rglob("*")) / 1e9
        print(f"  {accession}: {gb:.2f} GB, {n_verified} files md5-verified")
    return dest


def verify_one(accession: str) -> tuple[bool, str]:
    dest = genome_dir(accession)
    if not is_done(accession):
        return False, "not fetched"
    n = 0
    for line in (dest / "md5sum.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        target = dest / parts[1]
        if not target.exists():
            return False, f"missing file: {parts[1]}"
        if _md5(target) != parts[0]:
            return False, f"md5 mismatch: {parts[1]}"
        n += 1
    return n > 0, f"{n} files re-verified"


# ---------------------------------------------------------------- FASTA index

def build_fai(fna: Path) -> Path:
    """faidx-style index (name len offset linebases linewidth), reused if fresh."""
    fai = fna.with_suffix(fna.suffix + ".s5fai")
    if fai.exists() and fai.stat().st_mtime >= fna.stat().st_mtime:
        return fai
    rows = []
    with open(fna, "rb") as fh:
        name, seqlen, offset, lb, lw, pos = None, 0, 0, 0, 0, 0
        for raw in fh:
            n = len(raw)
            if raw.startswith(b">"):
                if name is not None:
                    rows.append((name, seqlen, offset, lb, lw))
                name = raw[1:].split()[0].decode()
                seqlen, offset, lb, lw = 0, pos + n, 0, 0
            else:
                stripped = len(raw.rstrip(b"\r\n"))
                if lb == 0:
                    lb, lw = stripped, n
                seqlen += stripped
            pos += n
        if name is not None:
            rows.append((name, seqlen, offset, lb, lw))
    with open(fai, "w") as out:
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    return fai


def read_fai(fai: Path) -> dict[str, tuple[int, int, int, int]]:
    idx = {}
    for line in fai.read_text().splitlines():
        name, ln, off, lb, lw = line.split("\t")
        idx[name] = (int(ln), int(off), int(lb), int(lw))
    return idx


def fetch_region(fna: Path, idx: dict, contig: str, start: int, end: int) -> str:
    """1-based inclusive region."""
    ln, off, lb, lw = idx[contig]
    start, end = max(1, start), min(ln, end)
    if end < start:
        return ""
    s0 = off + (start - 1) // lb * lw + (start - 1) % lb
    e0 = off + (end - 1) // lb * lw + (end - 1) % lb
    with open(fna, "rb") as fh:
        fh.seek(s0)
        raw = fh.read(e0 - s0 + 1)
    return raw.replace(b"\n", b"").replace(b"\r", b"").decode()


def longest_n_run(seq: str) -> int:
    return max((len(m) for m in _N_RUN.findall(seq)), default=0)


def assembly_stats(fna: Path) -> dict:
    """Total length, sequence count and N50 measured from the file itself.

    The manifest's N50 comes from NCBI's report of the assembly; this is the
    file that was searched, so D4's bar is judged against it.
    """
    lengths = sorted((v[0] for v in read_fai(build_fai(fna)).values()),
                     reverse=True)
    total = sum(lengths)
    acc, n50 = 0, 0
    for L in lengths:
        acc += L
        if acc >= total / 2:
            n50 = L
            break
    return {"total_bp": total, "n_sequences": len(lengths), "n50": n50,
            "longest": lengths[0] if lengths else 0}


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)
