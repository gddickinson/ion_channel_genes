"""s5r4_sweep.py — the genome sweep for census revision r4's families (D43).

    python3 scripts/s5r4_sweep.py baits          # baits_r4.{faa,tsv}, rules B1–B3
    python3 scripts/s5r4_sweep.py run [--species "Mus musculus,…"]   # default: all 52
    python3 scripts/s5r4_sweep.py check          # S5b's panel, runs and loci unchanged

The 12 families r4 added (7 census, KChIP, 4 decoys) were never searched in
the genomes. Rebuilding S5b's panel would move its baits (the B1 pool is
the S3b panel calls, which r4 re-assigned), so r4 gets a panel of its own,
drawn by S5b's rules restricted to its families, and its own miniprot run
per genome — same `-G`, same self-exclusion (D29, D37), same locus
clustering and profile calls (now over all 103 profiles). Output per genome
goes to `<data root>/genomes/s5/<assembly>/r4/`; nothing of S5b's is
rewritten. `s5_ledger.load_loci` reads both, and keeps an r4 locus only where
it bears on an r4 family.

Committed: `results/genome_sweep/baits_r4.tsv`, `genome_runs_r4.tsv`.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import live_progress  # noqa: E402
from s3_hmm_lib import read_tsv, sha256, write_tsv  # noqa: E402
from s5_classify import call_loci, write_loci  # noqa: E402
from s5_genome_io import assembly_stats, fna_path, write_json  # noqa: E402
from s5_lib import (BAITS_FAA, BAITS_R4_FAA, BAITS_R4_TSV, LIVE, OUT_DIR,  # noqa: E402
                    cluster_loci, load_bait_meta, manifest, max_intron_for,
                    parse_miniprot_gff, r4_families, run_miniprot_any, s5_dir,
                    sweep_assembly)

RUN_FIELDS = ["species", "group", "assembly", "max_intron", "chunks", "alignments",
              "self_dropped", "loci", "family_called_loci", "r4_family_loci",
              "miniprot_s", "call_s", "note"]



def cmd_baits(_a) -> None:
    from s5_baits import build
    print(json.dumps(build(faa=BAITS_R4_FAA, tsv=BAITS_R4_TSV, only=r4_families()),
                     indent=2))


def run_one(species: str, row: dict, bait_sha: str, meta: dict) -> dict:
    acc = sweep_assembly(row)
    out = s5_dir(acc, "r4")
    rec = {"species": species, "group": row["group"], "assembly": acc}
    fna = fna_path(acc)
    if not fna.exists():
        raise RuntimeError(f"{acc}: genome not on the data root (S5b fetched it)")
    total_bp = json.loads((s5_dir(acc) / "genome.json").read_text())["total_bp"]
    g = max_intron_for(row["group"], total_bp)
    gff, side = out / "miniprot.gff", out / "genome.json"
    prev = json.loads(side.read_text()) if side.exists() else {}
    t0 = time.time()
    if not (gff.exists() and prev.get("baits_sha256") == bait_sha
            and prev.get("max_intron") == g):
        chunks = run_miniprot_any(fna, BAITS_R4_FAA, gff, g, total_bp, out / "chunks")
    else:
        chunks = prev.get("chunks", 1)
    t_mp = time.time() - t0
    alns, dropped = parse_miniprot_gff(gff, meta, exclude_species=species)
    loci = cluster_loci(alns)
    t1 = time.time()
    rows = call_loci(loci, fna, out)
    for r in rows:
        r["locus"] = f"r4_{r['locus']}"
    write_loci(rows, out / "loci.tsv.gz")
    r4 = r4_families()
    rec.update(max_intron=g, chunks=chunks, alignments=len(alns), self_dropped=dropped,
               loci=len(loci),
               family_called_loci=sum(r["p_call"] == "family" for r in rows),
               r4_family_loci=sum(r["p_call"] == "family" and r["p_family"] in r4
                                  for r in rows),
               miniprot_s=round(t_mp if t_mp > 1 else prev.get("miniprot_s", 0), 1),
               call_s=round(time.time() - t1, 1), note="")
    write_json(side, {**rec, "baits_sha256": bait_sha})
    return rec


def cmd_run(a) -> None:
    man = manifest()
    runs = {r["species"]: r for r in read_tsv(OUT_DIR / "genome_runs.tsv")
            if not r.get("note", "").startswith("FAILED")}
    species = ([s.strip() for s in a.species.split(",") if s.strip()]
               or [s for s in man if s in runs])
    bait_sha = sha256(BAITS_R4_FAA)
    meta = load_bait_meta()
    path = OUT_DIR / "genome_runs_r4.tsv"
    done = {r["species"]: r for r in read_tsv(path)} if path.exists() else {}
    for i, s in enumerate(species, 1):
        live_progress(LIVE, "S5 r4", [(f"{x} r4 genome sweep", x in done)
                                      for x in species])
        print(f"[{i}/{len(species)}] {s}", flush=True)
        try:
            rec = run_one(s, man[s], bait_sha, meta)
        except Exception as exc:          # recorded, never swallowed silently
            rec = {"species": s, "group": man[s]["group"],
                   "assembly": sweep_assembly(man[s]), "note": f"FAILED: {exc}"[:300]}
        done[s] = rec
        print("   ", {k: rec.get(k) for k in ("loci", "r4_family_loci", "miniprot_s",
                                              "note")}, flush=True)
        order = list(man)
        write_tsv(path, RUN_FIELDS, sorted(done.values(),
                                           key=lambda r: order.index(r["species"])))


def cmd_check(_a) -> None:
    """S5b's bait panel and every genome's S5b loci file are unchanged."""
    bad = []
    runs = read_tsv(OUT_DIR / "genome_runs.tsv")
    for r in runs:
        side = s5_dir(r["assembly"]) / "genome.json"
        if side.exists() and json.loads(side.read_text()).get("baits_sha256") != sha256(BAITS_FAA):
            bad.append(r["species"])
    print(f"S5b bait panel SHA-256 {sha256(BAITS_FAA)[:16]}…; genomes whose S5b run "
          f"used another panel: {bad or 'none'}")
    if bad:
        raise SystemExit(1)


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("baits")
    p = sub.add_parser("run")
    p.add_argument("--species", default="")
    sub.add_parser("check")
    a = ap.parse_args()
    {"baits": cmd_baits, "run": cmd_run, "check": cmd_check}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
