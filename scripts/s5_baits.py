"""s5_baits.py — the S5 bait panel, derived from census v3 by stated rules.

One panel for every genome: every catalogue family, **controls included**.
A control family's baits are there for the same reason the IP3R project put
RyR baits beside its ITPR baits — without them a KCTD locus is claimed by a
Kv bait at partial coverage and reads out of the ledger as a Kv fragment. The
locus call itself is not made from the baits (the S3a profiles make it, D32);
the baits only find loci.

Rules (enforced here, recorded per bait in `baits.tsv`):

  B1  pool: S3b panel entries whose profile call is `family` F at **high**
      confidence, not a UniProt `(Fragment)`, length within F's catalogue
      band widened by BAND_SLACK either side (no length filter where the
      catalogue declares no band).
  B2  spread (D8): best entry (profile score) per species, then round-robin
      over panel groups — group order by name, species within a group by
      score — until MAX_BAITS. Breadth first, because a genome is found by
      the nearest bait, not the best one.
  B3  fallback: a family with no B1 entry takes its catalogue exemplars
      (`reference_panel.fasta`), up to MAX_BAITS.

Self-exclusion is applied at readout, not here: a genome is never read with
baits drawn from its own species (the genome-scale leave-one-out, D29).

    python3 scripts/s5_baits.py            # build baits.faa + baits.tsv
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from src.catalogue import registry  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_fasta, write_tsv  # noqa: E402
from s3b_lib import panel_db, s3b_dir  # noqa: E402
from s5_lib import BAITS_FAA, BAITS_TSV, MANIFEST  # noqa: E402

MAX_BAITS = 99            # one per species (D37): no cap in practice
BAND_SLACK = 0.25
REFERENCE_PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
BAIT_FIELDS = ["bait", "family", "superfamily", "family_status", "rule",
               "accession", "species", "group", "length", "win_score"]


def _in_band(length: int, band: tuple[int, int]) -> bool:
    lo, hi = band
    if not lo and not hi:
        return True
    return lo * (1 - BAND_SLACK) <= length <= hi * (1 + BAND_SLACK)


def _species_code_map() -> dict[str, str]:
    """Two-letter exemplar prefix ('Hs', 'Dm') → panel species, if unique."""
    by_code = defaultdict(set)
    for r in read_tsv(MANIFEST):
        g, s = (r["species"].split() + ["", ""])[:2]
        if g and s:
            by_code[(g[0] + s[0]).capitalize()].add(r["species"])
    return {c: next(iter(v)) for c, v in by_code.items() if len(v) == 1}


def pool_rows(families: dict) -> dict[str, list[dict]]:
    """B1 — high-confidence, complete, in-band profile calls per family."""
    frag = set()
    for head, _ in iter_fasta(panel_db()):
        if "(Fragment)" in head:
            frag.add(head.split()[0])
    pool = defaultdict(list)
    for r in read_tsv(s3b_dir() / "census_v3.tsv.gz"):
        f = r["p_family"]
        if (r["p_call"] != "family" or r["p_confidence"] != "high"
                or f not in families or r["target"] in frag):
            continue
        if _in_band(int(r["length"]), families[f].length_band_aa):
            pool[f].append(r)
    return pool


def spread(rows: list[dict], cap: int = MAX_BAITS) -> list[dict]:
    """B2 — best per species, then round-robin over groups."""
    best = {}
    for r in rows:
        s = r["species"]
        if s not in best or float(r["win_score"]) > float(best[s]["win_score"]):
            best[s] = r
    by_group = defaultdict(list)
    for r in best.values():
        by_group[r["group"]].append(r)
    for g in by_group:
        by_group[g].sort(key=lambda r: (-float(r["win_score"]), r["target"]))
    out, depth = [], 0
    while len(out) < cap and any(len(v) > depth for v in by_group.values()):
        for g in sorted(by_group):
            if depth < len(by_group[g]) and len(out) < cap:
                out.append(by_group[g][depth])
        depth += 1
    return out


def build(cap: int = MAX_BAITS, faa: Path = BAITS_FAA,
          tsv: Path = BAITS_TSV) -> dict:
    families = {f.key: f for f in registry.families()}
    pool = pool_rows(families)
    seqs = {h.split()[0]: s for h, s in iter_fasta(panel_db())}
    code_map = _species_code_map()
    group_of = {r["species"]: r["group"] for r in read_tsv(MANIFEST)}
    exemplars = defaultdict(list)
    for head, seq in iter_fasta(REFERENCE_PANEL):
        name, fam, acc = (head.split("|") + ["", ""])[:3]
        exemplars[fam].append((name, acc, seq))

    rows, fasta = [], []
    for key in sorted(families):
        f = families[key]
        chosen = spread(pool.get(key, []), cap)
        if chosen:
            for r in chosen:
                bait = f"{key}|{r['accession']}"
                rows.append({"bait": bait, "family": key,
                             "superfamily": f.superfamily,
                             "family_status": f.status.name.lower(),
                             "rule": "B1B2", "accession": r["accession"],
                             "species": r["species"], "group": r["group"],
                             "length": r["length"],
                             "win_score": r["win_score"]})
                fasta.append((bait, seqs[r["target"]]))
        else:
            for name, acc, seq in exemplars.get(key, [])[:cap]:
                sp = code_map.get(name.split("_")[0], "")
                bait = f"{key}|{acc or name}"
                rows.append({"bait": bait, "family": key,
                             "superfamily": f.superfamily,
                             "family_status": f.status.name.lower(),
                             "rule": "B3", "accession": acc, "species": sp,
                             "group": group_of.get(sp, ""),
                             "length": len(seq), "win_score": ""})
                fasta.append((bait, seq))
    faa.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(faa, fasta)
    write_tsv(tsv, BAIT_FIELDS, rows)
    missing = sorted(k for k in families if not any(r["family"] == k for r in rows))
    return {"baits": len(rows), "residues": sum(len(s) for _, s in fasta),
            "families": len({r["family"] for r in rows}),
            "b3_families": sorted({r["family"] for r in rows if r["rule"] == "B3"}),
            "families_without_bait": missing, "cap": cap, "baits_sha256": sha256(faa)}


if __name__ == "__main__":
    import json
    print(json.dumps(build(), indent=2))
