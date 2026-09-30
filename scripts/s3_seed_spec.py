"""S3 — the rules that choose each family profile's seed sequences.

One profile per catalogue family (census *and* control families: a best-
profile assignment needs the decoys to have profiles too, or every AChBP
is called a nicotinic receptor by default). Seeds come from three sources,
in this order, and each seed carries the rule that admitted it:

R1  **Curated human genes.** The catalogue's `human_genes` for the family,
    resolved to the reviewed human UniProt entry. The catalogue is the
    project's definition of its subject; this is that definition, as
    sequences. (Resolving a symbol to an accession is how the catalogue's
    own curated list becomes sequences — it is not a classification, H15.)
R2  **Catalogue exemplars** — `reference_panel.fasta` (S0-verified).
R3  **Census-derived.** Reviewed, non-fragment census v2 records S2 called
    to the family by a positive test (architecture, hazard or filter motif —
    all evidence independent of profiles), inside the family's length band
    widened by `BAND_SLACK`, **one per species**, chosen round-robin across
    taxonomic groups so no clade dominates, at most `MAX_R3`.

Rules the build enforces (a violation aborts, it is not repaired):
* a seed accession belongs to exactly one family;
* R3 never re-admits an accession already an R1/R2 seed of any family;
* identical sequences within a family are kept once.

Families with no R1 genes and no S2 family calls get R2 only (one to three
sequences). Those thin profiles are reported as such, not padded.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import s0_lib  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    OUT_DIR, REFERENCE_PANEL, iter_census_v2, iter_fasta, read_tsv,
    s3_dir, write_tsv,
)
from src.catalogue import registry  # noqa: E402

MAX_R3 = 30
BAND_SLACK = 0.25          # R3 length must lie in [lo×0.75, hi×1.25]
HUMAN_TAXON = 9606

HUMAN_ACC_TSV = OUT_DIR / "seed_human_accessions.tsv"
SEED_FIELDS = ["family", "superfamily", "accession", "rule", "taxon_id",
               "organism", "group", "length", "v2_family", "v2_tier"]


# ------------------------------------------------------------------- R1
def resolve_primary(f: s0_lib.Fetcher, gene: str) -> dict | None:
    """Reviewed human entry whose **primary** gene name is `gene`.

    `s0_lib.resolve_gene` queries `gene_exact`, which UniProt matches against
    synonyms too, and keeps the longest hit. Measured on the first S3 build:
    TRPC7 resolved to TRPM2 (O94759, whose synonyms include TRPC7), GRIK2 to
    GRIK5, KCNG3 to KCNG4, GJC3 to GJE1, AQP7 to AQP9 — seven collisions,
    caught by the seed-disjointness rule. A seed must be the gene itself.
    """
    q = s0_lib.urllib.parse.quote(
        f"gene_exact:{gene} AND organism_id:{HUMAN_TAXON} AND reviewed:true")
    d = f.json(f"{s0_lib.UNIPROT}/uniprotkb/search?query={q}"
               "&fields=accession,gene_primary,length&size=25")
    for r in (d or {}).get("results", []):
        prim = ((r.get("genes") or [{}])[0].get("geneName") or {}).get("value", "")
        if prim.upper() == gene.upper():
            return {"accession": r["primaryAccession"],
                    "length": (r.get("sequence") or {}).get("length")}
    return None


def human_accessions(refresh: bool = False) -> dict[str, dict]:
    """Every catalogue human gene → reviewed human accession (cached TSV)."""
    genes = registry.human_genes(census_only=False)        # symbol → family
    cached = {} if refresh or not HUMAN_ACC_TSV.exists() else {
        r["gene"]: r for r in read_tsv(HUMAN_ACC_TSV)}
    missing = [g for g in genes if g not in cached]
    if missing:
        f = s0_lib.Fetcher()
        for g in missing:
            hit = resolve_primary(f, g)
            cached[g] = {"gene": g, "family": genes[g],
                         "accession": hit["accession"] if hit else "",
                         "length": hit["length"] if hit else "",
                         "status": "ok" if hit else "unresolved"}
        write_tsv(HUMAN_ACC_TSV, ["gene", "family", "accession", "length",
                                  "status"],
                  [cached[g] for g in sorted(cached)])
        if f.failures:
            raise SystemExit(f"{len(f.failures)} UniProt failures resolving "
                             f"human genes: {f.failures[:3]}")
    # the cache holds the accession; the *family* is always the catalogue's
    # current one (a gene moved between families — TMEM87B, r4 — must not
    # keep seeding its old family)
    out = {g: {**r, "family": genes[g]} for g, r in cached.items() if g in genes}
    by_acc: dict[str, list[str]] = defaultdict(list)
    for g, r in out.items():
        if r["accession"]:
            by_acc[r["accession"]].append(g)
    dup = {a: gs for a, gs in by_acc.items() if len(gs) > 1}
    if dup:
        raise SystemExit(f"human genes sharing one accession: {dup}")
    return out


# ------------------------------------------------------------------- R2
def exemplar_seeds() -> list[dict]:
    out = []
    for head, _seq in iter_fasta(REFERENCE_PANEL):
        label, fam, acc = head.split()[0].split("|")
        out.append({"family": fam, "accession": acc, "rule": "R2",
                    "label": label})
    return out


# ------------------------------------------------------------------- R3
def _band(fam) -> tuple[float, float]:
    lo, hi = fam.length_band_aa
    if not hi:
        return 0.0, float("inf")
    return lo * (1 - BAND_SLACK), hi * (1 + BAND_SLACK)


def r3_candidates(families: set[str]) -> dict[str, list[dict]]:
    """Reviewed, full-length, in-band census v2 family calls, per family."""
    cols = ("accession", "reviewed", "taxon_id", "organism", "group",
            "length", "fragment", "family", "decisive_tier")
    bands = {k: _band(registry.family(k)) for k in families}
    out: dict[str, list[dict]] = defaultdict(list)
    for r in iter_census_v2(cols):
        fam = r["family"]
        if (fam not in families or r["reviewed"] != "reviewed"
                # UniProt's flag: "Fragment(s)" is a fragment, "Precursor" is
                # not (S2b). The first S3a build read any flag as a fragment,
                # which left CLCC1 with one seed; the 91 frozen profiles keep
                # the seeds they were built with (seed_manifest.tsv).
                or "Fragment" in (r["fragment"] or "")
                or r["taxon_id"] == str(HUMAN_TAXON)):
            continue
        lo, hi = bands[fam]
        if not lo <= int(r["length"]) <= hi:
            continue
        out[fam].append(r)
    return out


def pick_r3(cands: list[dict], taken: set[str], cap: int = MAX_R3) -> list[dict]:
    """One per species, round-robin across groups, deterministic."""
    by_group: dict[str, list[dict]] = defaultdict(list)
    seen_taxa: set[str] = set()
    for r in sorted(cands, key=lambda r: r["accession"]):
        if r["accession"] in taken or r["taxon_id"] in seen_taxa:
            continue
        seen_taxa.add(r["taxon_id"])
        by_group[r["group"]].append(r)
    picked: list[dict] = []
    groups = sorted(by_group)
    while len(picked) < cap and any(by_group[g] for g in groups):
        for g in groups:
            if by_group[g] and len(picked) < cap:
                picked.append(by_group[g].pop(0))
    return picked


# ------------------------------------------------------------- manifest
def build_manifest(refresh_human: bool = False) -> list[dict]:
    """All seeds for all families, rule-tagged. Aborts on a rule violation."""
    fams = {f.key: f for f in registry.families()}
    rows: list[dict] = []
    owner: dict[str, str] = {}

    def admit(fam: str, acc: str, rule: str, **extra) -> None:
        if not acc:
            return
        if acc in owner and owner[acc] != fam:
            raise SystemExit(f"seed {acc} claimed by {owner[acc]} and {fam} "
                             f"({rule}) — the seed sets must be disjoint")
        if acc in owner:
            return
        owner[acc] = fam
        rows.append({"family": fam, "superfamily": fams[fam].superfamily,
                     "accession": acc, "rule": rule, **extra})

    for g, r in sorted(human_accessions(refresh_human).items()):
        admit(r["family"], r["accession"], "R1", taxon_id=HUMAN_TAXON,
              organism="Homo sapiens")
    for r in exemplar_seeds():
        admit(r["family"], r["accession"], "R2")

    taken = set(owner)
    cands = r3_candidates(set(fams))
    for fam in sorted(fams):
        for r in pick_r3(cands.get(fam, []), taken):
            admit(fam, r["accession"], "R3", taxon_id=r["taxon_id"],
                  organism=r["organism"], group=r["group"],
                  length=r["length"], v2_family=r["family"],
                  v2_tier=r["decisive_tier"])
    return rows


def write_manifest(rows: list[dict]) -> Path:
    p = OUT_DIR / "seed_manifest.tsv"
    write_tsv(p, SEED_FIELDS, sorted(rows, key=lambda r: (r["family"],
                                                          r["rule"],
                                                          r["accession"])))
    return p


def seed_sequence_cache() -> Path:
    return s3_dir() / "seed_sequences.fasta"
