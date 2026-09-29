"""S3a step 1 — one profile HMM per catalogue family.

    python3 scripts/s3_build_profiles.py [--workers 10] [--refresh-human]

Seeds come from `s3_seed_spec.build_manifest()` (rules R1–R3, enforced).
Each family's seeds are aligned with MAFFT L-INS-i **single-threaded** —
the parent project measured `--thread -1` to be non-reproducible (the same
seeds gave 4,933 and 4,908 match states on consecutive runs) — and built
with `hmmbuild`. MAFFT's return code, the row count and a ragged-alignment
check are hard failures: a silent MAFFT failure degrades to a star
alignment with no other symptom (D28).

Writes (data root) `hmmer/s3/{seeds,aln,profiles}/<family>.*` and
`profiles/all.hmm`; (committed) `results/census_v3/seed_manifest.tsv`,
`seed_human_accessions.tsv`, `profile_build.tsv` with the SHA-256 of every
seed set, alignment and profile, so a rebuild that drifts is visible.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    LIVE, OUT_DIR, census_v2_fasta, iter_fasta, read_fasta, s3_dir, sha256,
    write_fasta, write_tsv,
)
from scripts.s3_seed_spec import (  # noqa: E402
    build_manifest, seed_sequence_cache, write_manifest,
)
from src.catalogue import registry  # noqa: E402

MAFFT = ["mafft", "--localpair", "--maxiterate", "1000", "--thread", "1",
         "--quiet"]
BUILD_FIELDS = ["family", "superfamily", "status", "n_seeds", "n_r1", "n_r2",
                "n_r3", "n_dup_dropped", "aln_columns", "match_states",
                "min_len", "max_len", "seconds", "seeds_sha256",
                "aln_sha256", "hmm_sha256"]


def gather_sequences(accs: set[str]) -> dict[str, str]:
    """Seed sequences: census v2 FASTA first, then a cached UniProt fetch."""
    seqs: dict[str, str] = {}
    for head, seq in iter_fasta(census_v2_fasta()):
        acc = head.split("|", 1)[0]
        if acc in accs:
            seqs[acc] = seq
    cache = seed_sequence_cache()
    cached = read_fasta(cache) if cache.exists() else {}
    need = sorted(accs - set(seqs) - set(cached))
    if need:
        f = s0_lib.Fetcher()
        for acc in need:
            s = s0_lib.sequence(f, acc)
            if not s:
                raise SystemExit(f"no sequence for seed {acc}")
            cached[acc] = s
        write_fasta(cache, sorted(cached.items()))
    for acc in accs - set(seqs):
        seqs[acc] = cached[acc]
    return seqs


def check_alignment(aln: dict[str, str], n: int, fam: str) -> int:
    if len(aln) != n:
        raise SystemExit(f"{fam}: MAFFT returned {len(aln)} rows for {n} seeds")
    widths = {len(s) for s in aln.values()}
    if len(widths) != 1:
        raise SystemExit(f"{fam}: ragged alignment {sorted(widths)[:5]}")
    return widths.pop()


def hmm_length(path: Path) -> int:
    for line in path.open():
        if line.startswith("LENG"):
            return int(line.split()[1])
    raise SystemExit(f"no LENG in {path}")


def build_one(job: tuple[str, list[tuple[str, str]], dict]) -> dict:
    fam, items, counts = job
    t0 = time.time()
    seeds_p = s3_dir("seeds") / f"{fam}.fa"
    aln_p = s3_dir("aln") / f"{fam}.afa"
    hmm_p = s3_dir("profiles") / f"{fam}.hmm"
    write_fasta(seeds_p, items)
    if len(items) == 1:
        write_fasta(aln_p, items)
    else:
        with aln_p.open("w") as out:
            r = subprocess.run(MAFFT + [str(seeds_p)], stdout=out,
                               stderr=subprocess.PIPE, text=True)
        if r.returncode:
            raise SystemExit(f"{fam}: mafft exit {r.returncode}: {r.stderr[:300]}")
    width = check_alignment(read_fasta(aln_p), len(items), fam)
    r = subprocess.run(["hmmbuild", "--cpu", "1", "-n", fam, "--informat",
                        "afa", str(hmm_p), str(aln_p)],
                       capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"{fam}: hmmbuild exit {r.returncode}: {r.stderr[:300]}")
    lens = [len(s) for _, s in items]
    return {"family": fam, **counts, "n_seeds": len(items),
            "aln_columns": width, "match_states": hmm_length(hmm_p),
            "min_len": min(lens), "max_len": max(lens),
            "seconds": round(time.time() - t0, 1),
            "seeds_sha256": sha256(seeds_p), "aln_sha256": sha256(aln_p),
            "hmm_sha256": sha256(hmm_p)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--refresh-human", action="store_true")
    ap.add_argument("--only", default="", help="comma-separated families")
    a = ap.parse_args()

    manifest = build_manifest(a.refresh_human)
    print(f"[seeds] {len(manifest)} seeds → {write_manifest(manifest)}")
    seqs = gather_sequences({r["accession"] for r in manifest})

    jobs = []
    only = set(filter(None, a.only.split(",")))
    for fam in registry.families():
        if only and fam.key not in only:
            continue
        mine = [r for r in manifest if r["family"] == fam.key]
        items, seen, dup = [], set(), 0
        for r in sorted(mine, key=lambda r: (r["rule"], r["accession"])):
            s = seqs[r["accession"]].upper()
            if s in seen:
                dup += 1
                continue
            seen.add(s)
            items.append((r["accession"], s))
        counts = {"superfamily": fam.superfamily, "status": fam.status.value,
                  "n_dup_dropped": dup,
                  **{f"n_{k.lower()}": sum(r["rule"] == k for r in mine)
                     for k in ("R1", "R2", "R3")}}
        jobs.append((fam.key, items, counts))
    # Longest total work first so the pool does not end on one RyR build.
    jobs.sort(key=lambda j: -len(j[1]) * max(len(s) for _, s in j[1]) ** 2)

    rows, t0 = [], time.time()
    with Pool(a.workers) as pool:
        for i, row in enumerate(pool.imap_unordered(build_one, jobs), 1):
            rows.append(row)
            print(f"[build] {i}/{len(jobs)} {row['family']:26s} "
                  f"{row['n_seeds']:3d} seeds → {row['match_states']:5d} "
                  f"states ({row['seconds']} s)", flush=True)
            s0_lib.live_progress(LIVE, "S3a", [
                (f"profiles built {i}/{len(jobs)}", i == len(jobs))])
    rows.sort(key=lambda r: r["family"])
    write_tsv(OUT_DIR / "profile_build.tsv", BUILD_FIELDS, rows)

    prof_dir = s3_dir("profiles")
    with (prof_dir / "all.hmm").open("w") as out:
        for fam in registry.families():
            p = prof_dir / f"{fam.key}.hmm"
            if p.exists():
                out.write(p.read_text())
    print(f"[done] {len(rows)} profiles in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
