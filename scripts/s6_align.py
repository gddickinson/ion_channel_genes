"""s6_align.py — S6 step 1: one MAFFT L-INS-i + trimAl alignment per family.

    python3 scripts/s6_align.py sets             # D39 sets → members.tsv + FASTA
    python3 scripts/s6_align.py align [--jobs 3 --threads 3] [--only nav,cav]
    python3 scripts/s6_align.py stats            # → alignments.tsv

Every census family with ≥ 4 sequences in its D39 set gets an L-INS-i
alignment (`mafft --localpair --maxiterate 1000`) — the one method for every
family, whatever its size; no family is quietly moved to a faster
algorithm (D28). trimAl `-automated1` then trims it. A run is skipped when
its input FASTA's SHA-256 matches the one recorded next to the finished
alignment, so an interrupted sweep resumes by family.

Rooting is S7's: its outgroup sequences are added to these alignments with
`mafft --add --keeplength`, so the family columns do not move.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s0_lib import live_progress  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_fasta, write_tsv  # noqa: E402
from s6_lib import LIVE, OUT_DIR, build_sets, s6_dir  # noqa: E402

MIN_SEQS = 4
MEMBER_FIELDS = ["family", "label", "target", "species", "group", "source",
                 "p_confidence", "win_coverage", "length", "verdict"]


def cmd_sets(_a) -> None:
    members, sets = build_sets()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT_DIR / "members.tsv", MEMBER_FIELDS, members)
    d = s6_dir("family")
    for fam, seqs in sorted(sets.items()):
        write_fasta(d / f"{fam}.fasta", seqs)
    v = Counter(m["verdict"].split(":")[0] for m in members)
    print(f"{len(members)} called rows; {dict(v)}")
    print(f"{len(sets)} families, {sum(map(len, sets.values()))} sequences "
          f"→ {d}")


def _mafft(src: Path, dst: Path, threads: int) -> tuple[list[str], float]:
    cmd = ["mafft", "--localpair", "--maxiterate", "1000", "--thread",
           str(threads), "--quiet", str(src)]
    t0 = time.time()
    with open(dst.with_suffix(".tmp"), "w") as out:
        p = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"mafft rc={p.returncode} on {src.name}: "
                           f"{p.stderr.strip()[:300]}")
    dst.with_suffix(".tmp").rename(dst)
    return cmd, time.time() - t0


def _check(src: Path, aln: Path) -> None:
    """Hard failure on row count or ragged output (ported from S3a)."""
    raw = {h.split()[0]: s for h, s in iter_fasta(src)}
    got = {h.split()[0]: s for h, s in iter_fasta(aln)}
    if set(raw) != set(got):
        raise RuntimeError(f"{aln.name}: {len(got)} rows for {len(raw)} inputs")
    if len({len(s) for s in got.values()}) != 1:
        raise RuntimeError(f"{aln.name}: ragged alignment")
    for k, s in got.items():
        if s.replace("-", "").upper() != raw[k].upper():
            raise RuntimeError(f"{aln.name}: {k} residues changed by MAFFT")


def align_family(fam: str, threads: int) -> dict:
    d = s6_dir("family")
    src, aln, trm = d / f"{fam}.fasta", d / f"{fam}.linsi.fasta", d / f"{fam}.trimal.fasta"
    rec = d / f"{fam}.run.json"
    digest = sha256(src)
    if rec.exists() and aln.exists() and trm.exists():
        old = json.loads(rec.read_text())
        if old.get("input_sha256") == digest:
            return {**old, "resumed": True}
    cmd, secs = _mafft(src, aln, threads)
    _check(src, aln)
    tcmd = ["trimal", "-in", str(aln), "-out", str(trm), "-automated1"]
    p = subprocess.run(tcmd, capture_output=True, text=True)
    if p.returncode != 0 or not trm.exists():
        raise RuntimeError(f"trimal failed on {fam}: {p.stderr.strip()[:300]}")
    out = {"family": fam, "input_sha256": digest,
           "alignment_sha256": sha256(aln), "trimmed_sha256": sha256(trm),
           "mafft_cmd": " ".join(cmd), "trimal_cmd": " ".join(tcmd),
           "mafft_seconds": round(secs, 1), "finished": time.strftime("%F %T")}
    rec.write_text(json.dumps(out, indent=2))
    return out


def cmd_align(a) -> None:
    d = s6_dir("family")
    fams = sorted(p.stem for p in d.glob("*.fasta")
                  if "." not in p.stem)
    sizes = {f: sum(1 for _ in iter_fasta(d / f"{f}.fasta")) for f in fams}
    if a.only:
        fams = [f for f in fams if f in set(a.only.split(","))]
    todo = [f for f in fams if sizes[f] >= MIN_SEQS]
    # Largest first, so the long runs start while the small ones fill in.
    todo.sort(key=lambda f: -sizes[f])
    done: list[str] = []
    steps = lambda: [(f"L-INS-i {f} (n={sizes[f]})", f in done)  # noqa: E731
                     for f in todo]
    live_progress(LIVE, "S6", steps(), a.jobs)
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(align_family, f, a.threads): f for f in todo}
        for fu in as_completed(futs):
            f = futs[fu]
            r = fu.result()
            done.append(f)
            tag = "resumed" if r.get("resumed") else f"{r['mafft_seconds']:.0f} s"
            print(f"[{len(done)}/{len(todo)}] {f} n={sizes[f]} {tag}", flush=True)
            live_progress(LIVE, "S6", steps(), a.jobs)


def _stats(path: Path) -> dict:
    rows = [s for _, s in iter_fasta(path)]
    if not rows:
        return {"n": 0, "cols": 0, "gap_frac": ""}
    cols = len(rows[0])
    gaps = sum(s.count("-") for s in rows)
    return {"n": len(rows), "cols": cols,
            "gap_frac": round(gaps / (cols * len(rows)), 4) if cols else ""}


def cmd_stats(_a) -> None:
    d = s6_dir("family")
    members = read_tsv(OUT_DIR / "members.tsv")
    inc = Counter(m["family"] for m in members if m["verdict"] == "include")
    exc = Counter(m["family"] for m in members if m["verdict"] != "include")
    out = []
    for fam in sorted(set(inc) | set(exc)):
        row = {"family": fam, "included": inc[fam], "excluded": exc[fam]}
        rec = d / f"{fam}.run.json"
        if inc[fam] < MIN_SEQS:
            row["status"] = f"too_few (<{MIN_SEQS})"
        elif not rec.exists():
            row["status"] = "not_run"
        else:
            r = json.loads(rec.read_text())
            a, t = _stats(d / f"{fam}.linsi.fasta"), _stats(d / f"{fam}.trimal.fasta")
            row.update(status="aligned", aln_cols=a["cols"], aln_gap=a["gap_frac"],
                       trim_cols=t["cols"], trim_gap=t["gap_frac"],
                       mafft_seconds=r["mafft_seconds"],
                       input_sha256=r["input_sha256"],
                       trimmed_sha256=r["trimmed_sha256"])
        out.append(row)
    fields = ["family", "included", "excluded", "status", "aln_cols", "aln_gap",
              "trim_cols", "trim_gap", "mafft_seconds", "input_sha256",
              "trimmed_sha256"]
    write_tsv(OUT_DIR / "alignments.tsv", fields, out)
    print(Counter(r["status"] for r in out))


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("sets")
    p = sub.add_parser("align")
    p.add_argument("--jobs", type=int, default=3)
    p.add_argument("--threads", type=int, default=3)
    p.add_argument("--only", default="")
    sub.add_parser("stats")
    a = ap.parse_args()
    {"sets": cmd_sets, "align": cmd_align, "stats": cmd_stats}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
