"""s5_sweep.py — S5 driver: one genome at a time, fetch → align → call.

    python3 scripts/s5_sweep.py baits                  # build the bait panel
    python3 scripts/s5_sweep.py run --pilot            # the S5a pilot genomes
    python3 scripts/s5_sweep.py run --species "Mus musculus,Cornu aspersum"
    python3 scripts/s5_sweep.py run --all              # S5b

Per genome, under `<data root>/genomes/s5/<assembly>/`: `miniprot.gff`
(reused only if the bait panel's SHA-256 matches), `loci.tsv.gz` (one row
per locus with its profile call, `s5_classify.py`) and `genome.json`
(assembly stats measured from the searched file, `-G`, chunks, timings,
self-excluded alignments). Committed per-genome summary:
`results/genome_sweep/genome_runs.tsv`.
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
from s5_genome_io import assembly_stats, fetch_one, fna_path, write_json  # noqa: E402
from s5_lib import (BAITS_FAA, LIVE, OUT_DIR, cluster_loci, load_bait_meta,  # noqa: E402
                    manifest, max_intron_for, parse_miniprot_gff,
                    run_miniprot_any, s5_dir, sweep_assembly)

#: The S5a pilot: one genome per question the instrument must answer.
PILOT = ("Escherichia coli",          # prokaryote; Kir / TMEM175 zero cells
         "Saccharomyces cerevisiae",  # compact fungus
         "Arabidopsis thaliana",      # plant; Kir, ORAI, Slo, Hv1, tweety zeros
         "Caenorhabditis elegans",    # Nav, P2X, HCN, TPC zeros
         "Drosophila melanogaster",   # P2X, TPC, Hv1, LRRC8 zeros
         "Mus musculus",              # a mammalian genome; the ZAC zero
         "Cornu aspersum")            # genome-only: no proteome to compare
RUN_FIELDS = ["species", "group", "status", "assembly", "total_bp", "n_sequences",
              "n50", "max_intron", "chunks", "alignments", "self_dropped",
              "loci", "family_called_loci", "miniprot_s", "call_s", "note"]


def run_one(species: str, row: dict, bait_sha: str, meta: dict) -> dict:
    acc = sweep_assembly(row)
    out = s5_dir(acc)
    rec = {"species": species, "group": row["group"], "status": row["status"],
           "assembly": acc}
    fetch_one(acc)
    fna = fna_path(acc)
    stats = assembly_stats(fna)
    g = max_intron_for(row["group"], stats["total_bp"])
    gff, side = out / "miniprot.gff", out / "genome.json"
    prev = json.loads(side.read_text()) if side.exists() else {}
    t0 = time.time()
    if not (gff.exists() and prev.get("baits_sha256") == bait_sha
            and prev.get("max_intron") == g):
        chunks = run_miniprot_any(fna, BAITS_FAA, gff, g, stats["total_bp"],
                                  out / "chunks")
    else:
        chunks = prev.get("chunks", 1)
    t_mp = time.time() - t0
    alns, dropped = parse_miniprot_gff(gff, meta, exclude_species=species)
    loci = cluster_loci(alns)
    t1 = time.time()
    rows = call_loci(loci, fna, out)
    write_loci(rows, out / "loci.tsv.gz")
    rec.update(total_bp=stats["total_bp"], n_sequences=stats["n_sequences"],
               n50=stats["n50"], max_intron=g, chunks=chunks,
               alignments=len(alns), self_dropped=dropped, loci=len(loci),
               family_called_loci=sum(r["p_call"] == "family" for r in rows),
               miniprot_s=round(t_mp if t_mp > 1 else prev.get("miniprot_s", 0), 1),
               call_s=round(time.time() - t1, 1), note="")
    write_json(side, {**rec, **stats, "baits_sha256": bait_sha,
                      "fna": fna.name})
    return rec


def run(species_list: list[str]) -> None:
    man = manifest()
    bait_sha = sha256(BAITS_FAA)
    meta = load_bait_meta()
    runs_path = OUT_DIR / "genome_runs.tsv"
    done = {r["species"]: r for r in read_tsv(runs_path)} if runs_path.exists() else {}
    steps = [[s, False] for s in species_list]
    for i, s in enumerate(species_list):
        live_progress(LIVE, "S5a", [(f"{x} genome sweep", d) for x, d in steps])
        print(f"[{i + 1}/{len(species_list)}] {s}", flush=True)
        try:
            rec = run_one(s, man[s], bait_sha, meta)
        except Exception as exc:          # recorded, never swallowed silently
            rec = {"species": s, "group": man[s]["group"],
                   "status": man[s]["status"], "assembly": sweep_assembly(man[s]),
                   "note": f"FAILED: {exc}"[:300]}
            print("   ", rec["note"], flush=True)
        done[s] = rec
        steps[i][1] = True
        print("   ", {k: rec.get(k) for k in ("loci", "family_called_loci",
                                              "miniprot_s", "call_s")}, flush=True)
        order = list(man)
        write_tsv(runs_path, RUN_FIELDS,
                  sorted(done.values(), key=lambda r: order.index(r["species"])))
    live_progress(LIVE, "S5a", [(f"{x} genome sweep", True) for x, _ in steps])


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("baits")
    r = sub.add_parser("run")
    r.add_argument("--pilot", action="store_true")
    r.add_argument("--all", action="store_true")
    r.add_argument("--species", default="")
    a = ap.parse_args()
    if a.cmd == "baits":
        from s5_baits import build
        print(json.dumps(build(), indent=2))
        return 0
    species = (list(PILOT) if a.pilot else list(manifest()) if a.all
               else [s.strip() for s in a.species.split(",") if s.strip()])
    run(species)
    return 0


if __name__ == "__main__":
    sys.exit(main())
