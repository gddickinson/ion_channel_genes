"""s5_classify.py — call every genomic locus with the S3a profiles (D32).

A locus is found by baits and **called by profiles**: the top translations
of each bait family at the locus are scored against all 91 frozen S3a
profiles (`profiles/all.hmm`, `-Z` fixed at the panel size so the bit scores
are S3b's), and the locus takes the D32 call of its best-scoring
translation. The bait family and its margin are recorded beside it as a
second witness; disagreement is counted, never resolved by preference.

Per-genome output: `<data root>/genomes/s5/<acc>/loci.tsv.gz`.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import s3_dir, write_fasta, write_tsv  # noqa: E402
from s3_assign import assign_one, collect_hits, fam_superfamily  # noqa: E402
from s5_genome_io import fetch_region, longest_n_run, read_fai, build_fai  # noqa: E402
from s5_lib import EDGE_BP, Locus  # noqa: E402

PANEL_Z = 822_499            # S4's sweep DB size — S3b's -Z, so scores compare
LOCUS_FIELDS = [
    "locus", "contig", "strand", "start", "end", "span", "contig_len",
    "at_edge", "longest_n_run", "n_alns", "best_identity", "best_coverage",
    "best_bait", "bait_family", "bait_runner", "bait_margin",
    "p_call", "p_family", "p_superfamily", "p_confidence", "win_score",
    "win_coverage", "runner", "rel_margin", "band", "call_translation",
    "call_bait", "call_identity", "call_bait_coverage", "call_max_intron",
    "call_start", "call_end", "call_frameshifts", "call_stops", "call_intact", "call_cds",
    "agree"]


def _clean(seq: str) -> str:
    return seq.replace("*", "X").replace("-", "")


def score_translations(loci: list[Locus], work: Path, cpu: int = 8
                       ) -> dict[str, list[tuple]]:
    """hmmsearch all profiles over the loci's candidate translations."""
    items, seen = [], set()
    for i, L in enumerate(loci):
        for a in L.best_per_family(2):
            t = _clean(a.translation)
            if len(t) >= 30 and a.mp_id not in seen:
                seen.add(a.mp_id)
                items.append((f"L{i}|{a.mp_id}", t))
    faa, dom = work / "translations.faa", work / "translations.domtbl"
    write_fasta(faa, items)
    if not items:
        return {}
    cmd = ["hmmsearch", "--cpu", str(cpu), "--noali", "-E", "1e-3",
           "--domE", "1e-3", "-Z", str(PANEL_Z), "--domZ", str(PANEL_Z),
           "-o", "/dev/null", "--domtblout", str(dom),
           str(s3_dir("profiles") / "all.hmm"), str(faa)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"hmmsearch exit {r.returncode}: {r.stderr[:300]}")
    return collect_hits([dom])


def call_loci(loci: list[Locus], fna: Path, work: Path) -> list[dict]:
    """One row per locus: geometry, bait witness, profile call."""
    hits = score_translations(loci, work)
    sf_of = fam_superfamily()
    idx = read_fai(build_fai(fna))
    by_locus: dict[int, list[tuple[str, list]]] = {}
    for tid, h in hits.items():
        by_locus.setdefault(int(tid.split("|")[0][1:]), []).append((tid, h))
    alns = {a.mp_id: a for L in loci for a in L.alns}
    rows = []
    for i, L in enumerate(loci):
        best = L.best
        bf, br, bm = L.bait_family_margin()
        clen = idx[L.contig][0] if L.contig in idx else 0
        region = (fetch_region(fna, idx, L.contig, L.start, L.end)
                  if L.contig in idx and L.end - L.start < 5_000_000 else "")
        row = {"locus": f"{L.contig}:{L.start}-{L.end}{L.strand}",
               "contig": L.contig, "strand": L.strand, "start": L.start,
               "end": L.end, "span": L.end - L.start + 1, "contig_len": clen,
               "at_edge": int(L.start <= EDGE_BP or (clen and clen - L.end <= EDGE_BP)),
               "longest_n_run": longest_n_run(region), "n_alns": len(L.alns),
               "best_identity": round(best.identity, 4),
               "best_coverage": round(max(a.coverage for a in L.alns), 4),
               "best_bait": best.bait, "bait_family": bf, "bait_runner": br,
               "bait_margin": bm}
        cands = by_locus.get(i, [])
        if cands:
            tid, h = max(cands, key=lambda c: c[1][0][0])
            call = assign_one(tid, h, sf_of)
            a = alns[tid.split("|", 1)[1]]
            row.update({k: call.get(k, "") for k in (
                "p_call", "p_family", "p_superfamily", "p_confidence",
                "win_score", "win_coverage", "runner", "rel_margin", "band")})
            row.update(call_translation=tid, call_bait=a.bait,
                       call_identity=round(a.identity, 4),
                       call_bait_coverage=round(a.coverage, 4),
                       call_max_intron=a.max_intron,
                       call_start=a.start, call_end=a.end,
                       call_frameshifts=a.frameshifts, call_stops=a.stop_codons,
                       call_intact=int(a.frameshifts == 0 and a.stop_codons == 0),
                       call_cds=",".join(f"{x}-{y}" for x, y in sorted(a.cds_blocks)))
        else:
            row.update(p_call="no_hit", p_confidence="none")
        row["agree"] = int(row.get("p_family", "") == bf) if row.get("p_family") else ""
        rows.append(row)
    return rows


def write_loci(rows: list[dict], path: Path) -> int:
    return write_tsv(path, LOCUS_FIELDS, rows)
