"""s6_lib.py — S6 shared pieces: the alignment sets and where their files live.

**The set rule (D39), fixed before any alignment was run.** A family's tier-1
alignment set is every census v4 row the S3a profiles call to that family at
**high** confidence (D32 margin ≥ 0.30 *and* ≥ half the profile covered):

* proteome rows with `v3_status == channel`, `p_family == v3_family`,
  `p_confidence == high`;
* genome loci (`v3_status == genome_locus`) with the same call on an
  **intact** reading frame (D37 (4)) — the proteome misses and the two
  genome-only species.

Everything else a family's census row carries is excluded *and counted* in
`members.tsv` with the reason: medium-confidence calls (the S3b ankyrin/LRR
upper bound, and sister-family calls inside 0.30 of the margin), S2-only
calls with no profile call, broken genome frames. Identical sequences are
collapsed to one representative and the duplicates listed, never dropped
silently. Bulk FASTA/alignments: `<data root>/alignments/s6/`.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import iter_fasta, seq_id  # noqa: E402
from s3b_lib import panel_db  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT_DIR = ROOT / "results" / "alignments"
LIVE = ROOT / "results" / "session_live.json"
CALLED_STATUS = ("channel", "genome_locus")


def s6_dir(*parts: str) -> Path:
    """`<data root>/alignments/s6/...` — raises without the drive (D1)."""
    p = require_data_root() / "alignments" / "s6"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def census_v4() -> Path:
    return require_data_root() / "genomes" / "s5" / "census_v4.tsv.gz"


def s3_profiles() -> Path:
    return require_data_root() / "hmmer" / "s3" / "profiles"


def iter_census_v4():
    with gzip.open(census_v4(), "rt") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def verdict(row: dict) -> str:
    """'include' or the exclusion reason — the D39 rule, one place."""
    if row["v3_status"] not in CALLED_STATUS:
        return "not_called"
    if not row["p_family"]:
        return "no_profile_call"
    if row["p_family"] != row["v3_family"]:
        return "profile_other_family"
    if row["p_confidence"] != "high":
        return f"profile_{row['p_confidence'] or 'none'}"
    if row["v3_status"] == "genome_locus" and row["call_intact"] != "1":
        return "genome_not_intact"
    return "include"


def genome_translations(needed: dict[str, set[str]]) -> dict[tuple[str, str], str]:
    """(assembly, locus) → called translation, for the genome rows needed."""
    s5 = require_data_root() / "genomes" / "s5"
    out: dict[tuple[str, str], str] = {}
    for asm, loci in needed.items():
        want: dict[str, str] = {}
        with gzip.open(s5 / asm / "loci.tsv.gz", "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["locus"] in loci:
                    want[r["call_translation"]] = r["locus"]
        for head, seq in iter_fasta(s5 / asm / "translations.faa"):
            tid = head.split()[0]
            if tid in want:
                out[(asm, want[tid])] = seq.rstrip("*").replace("*", "X")
    return out


def proteome_sequences(targets: set[str]) -> dict[str, str]:
    out = {}
    for head, seq in iter_fasta(panel_db()):
        t = head.split()[0]
        if t in targets:
            out[t] = seq
    return out


def label_for(row: dict) -> str:
    """A tree-safe, unique label: accession or assembly:locus, plus species."""
    sp = "".join(w[:3] for w in row["species"].split()[:2])
    if row["source"] == "genome":
        core = f"{row['assembly']}_{row['locus']}"
    else:
        core = row["accession"]
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in core)
    return f"{safe}__{sp}"


def build_sets():
    """→ (members rows, {family: [(label, seq)]}) under D39."""
    rows = list(iter_census_v4())
    inc = [r for r in rows if verdict(r) == "include"]
    genome_need: dict[str, set[str]] = defaultdict(set)
    for r in inc:
        if r["source"] == "genome":
            genome_need[r["assembly"]].add(r["locus"])
    gseq = genome_translations(genome_need)
    pseq = proteome_sequences({r["target"] for r in inc if r["source"] != "genome"})

    members, sets = [], defaultdict(list)
    seen: dict[tuple[str, str], str] = {}
    for r in rows:
        if r["v3_status"] not in CALLED_STATUS:
            continue
        v = verdict(r)
        fam, lab = r["v3_family"], label_for(r)
        seq = ""
        if v == "include":
            seq = (gseq.get((r["assembly"], r["locus"]))
                   if r["source"] == "genome" else pseq.get(r["target"])) or ""
            if not seq:
                v = "sequence_missing"
            else:
                key = (fam, seq_id(seq))
                if key in seen:
                    v = f"duplicate_of:{seen[key]}"
                else:
                    seen[key] = lab
                    sets[fam].append((lab, seq))
        members.append({
            "family": fam, "label": lab, "target": r["target"],
            "species": r["species"], "group": r["group"], "source": r["source"],
            "p_confidence": r["p_confidence"], "win_coverage": r["win_coverage"],
            "length": len(seq) if seq else "", "verdict": v})
    return members, dict(sets)
