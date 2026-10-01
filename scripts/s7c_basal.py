"""s7c_basal.py — S7c step 1: denser basal sampling for the re-rooted families.

    python3 scripts/s7c_basal.py      # → results/phylogeny/reroot_basal.tsv

A long branch to a distant outgroup roots badly when nothing in the ingroup
breaks it. S7c therefore adds to each re-rooted family (every family with a
per-family `root_with`, read from the catalogue) members from the
**early-diverging lineages** of S4b's dense panel, by a rule fixed before any
re-rooted tree (D8, D47):

* **Evidence**: an S4b profile call to the family at high confidence (D32,
  the instrument D39 uses), sequence not flagged "(Fragment)", length inside
  the family's catalogue band widened by 25 % on each side (the R3 seed
  convention).
* **Lineage**: outside Bilateria — every non-metazoan eukaryote, plus
  Porifera, Placozoa, Ctenophora and Cnidaria.
* **Clade**: one representative per (kingdom, phylum, class), falling back to
  the order when the panel gives no class. Within a clade, the
  highest-scoring call (profile bits: the most typical member, never the
  longest); ties broken by accession.
* A pick whose sequence is already in the family's D39 set is recorded and
  not added twice.

This moves the six families' sequence sets away from S6's D39 sets; the
exception is stated here and in D47, never silent. Candidates' sequences
are cached in `<data root>/trees/s7c/basal_candidates.fasta`, so the 8.8 M-
sequence dense panel is streamed once.
"""

from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import iter_fasta, read_tsv, write_tsv  # noqa: E402
from s6_lib import s6_dir  # noqa: E402
from s7_lib import OUT_DIR, ROOT  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

NON_BILATERIAN = {"Porifera", "Placozoa", "Ctenophora", "Cnidaria"}
BAND_SLACK = 0.25
FIELDS = ["family", "clade", "kingdom", "phylum", "class", "order", "target",
          "label", "organism", "length", "win_score", "clade_candidates",
          "verdict"]


def reroot_families() -> list[str]:
    return sorted(k for k, f in CATALOGUE.items() if f.root_with)


def s7c_dir(*parts: str) -> Path:
    """`<data root>/trees/s7c/...` — raises without the drive (D1)."""
    p = require_data_root() / "trees" / "s7c"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def early(kingdom: str, phylum: str) -> bool:
    return kingdom != "Metazoa" or phylum in NON_BILATERIAN


def in_band(fam: str, n: int) -> bool:
    lo, hi = CATALOGUE[fam].length_band_aa
    return lo * (1 - BAND_SLACK) <= n <= hi * (1 + BAND_SLACK)


def label(target: str, organism: str) -> str:
    """The S6 label form: accession plus the first three letters of genus and species."""
    acc = target.split("|")[1] if "|" in target else target
    sp = "".join(w[:3] for w in organism.split()[:2])
    return f"{acc}__{sp}"


def candidates(fams: set[str]) -> list[dict]:
    bulk = require_data_root() / "proteomes" / "s4b"
    uni = {}
    with gzip.open(bulk / "universe.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if early(r["kingdom"], r["phylum"]):
                uni[r["target"]] = r
    out = []
    for r in read_tsv(bulk / "calls.tsv.gz"):
        if (r["p_family"] in fams and r["p_confidence"] == "high"
                and r["target"] in uni):
            out.append({**uni[r["target"]], "family": r["p_family"],
                        "target": r["target"], "win_score": float(r["win_score"])})
    return out


def sequences(targets: set[str]) -> dict[str, tuple[str, str]]:
    """{target: (header, seq)}, from the cache or one pass over the dense panel."""
    cache = s7c_dir() / "basal_candidates.fasta"
    have = {}
    if cache.exists():
        have = {h.split()[0]: (h, s) for h, s in iter_fasta(cache)}
    if not targets <= set(have):
        dense = require_data_root() / "proteomes" / "s4b" / "dense_panel.fasta"
        for h, s in iter_fasta(dense):
            t = h.split()[0]
            if t in targets:
                have[t] = (h, s)
        with open(cache, "w") as fh:
            for t in sorted(have):
                fh.write(f">{have[t][0]}\n{have[t][1]}\n")
    return have


def d39_seqs(fam: str) -> set[str]:
    return {s.upper() for _, s in iter_fasta(s6_dir("family") / f"{fam}.fasta")}


def main() -> int:
    fams = reroot_families()
    cands = candidates(set(fams))
    seqs = sequences({c["target"] for c in cands})
    by_clade: dict[tuple, list[dict]] = defaultdict(list)
    for c in cands:
        h, s = seqs[c["target"]]
        organism = h.split("OS=")[1].split(" OX=")[0] if "OS=" in h else ""
        c.update(organism=organism, length=len(s), seq=s,
                 fragment="(Fragment)" in h)
        if c["fragment"] or not in_band(c["family"], len(s)):
            continue
        clade = (c["kingdom"] or "-", c["phylum"] or "-",
                 c["class"] or f"order:{c['order']}")
        by_clade[(c["family"], clade)].append(c)
    rows, picks = [], defaultdict(list)
    have = {f: d39_seqs(f) for f in fams}
    for (fam, clade), cs in sorted(by_clade.items()):
        cs.sort(key=lambda c: (-c["win_score"], c["target"]))
        c = cs[0]
        dup = c["seq"].upper() in have[fam]
        rows.append({**c, "clade": " / ".join(clade),
                     "label": label(c["target"], c["organism"]),
                     "clade_candidates": len(cs),
                     "verdict": "in_d39_set" if dup else "add"})
        if not dup:
            picks[fam].append((rows[-1]["label"], c["seq"]))
            have[fam].add(c["seq"].upper())
    write_tsv(OUT_DIR / "reroot_basal.tsv", FIELDS, rows)
    for fam in fams:
        with open(s7c_dir() / f"{fam}.basal.fasta", "w") as fh:
            for lab, s in picks[fam]:
                fh.write(f">{lab}\n{s}\n")
        n_c = sum(1 for c in cands if c["family"] == fam)
        print(f"{fam:10s} {n_c:>4} high-confidence early-lineage calls → "
              f"{len(picks[fam]):>3} added "
              f"({sum(r['family'] == fam for r in rows)} clades)")
    print(f"→ {OUT_DIR / 'reroot_basal.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
