"""Channel-specific headless commands: classify, census and phylogeny.

Kept separate from `src/cli.py` because these three commands must run in a
bare interpreter. `src/cli.py` reaches the sequence-analysis pipeline, which
imports Biopython; classification is the project's core operation and cannot
be gated on an optional dependency.

    python run.py --classify Q14524,O43497,Q719H9
    python run.py --classify-preset channelome_human --save-results
    python run.py --catalogue --scope ploop
    python run.py --phylo nav --tier 1
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .catalogue import CATALOGUE, SUPERFAMILIES, stats
from .classify import ChannelQuery, ReferenceSet, classify_all
from .classify.motifs import FOUR_REPEAT_ANCHOR
from .classify.reference import fetch_uniprot_sequence
from .classify.report import (CALLS_HEADER, calls_table, confusion_table,
                              recall_table, summary_counts, write_tsv)
from .phylo import (NotAlignable, build_tier1, build_tier2,
                    needs_modules, write_run)
from .utils.scope import scope_for

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb"


# ----------------------------------------------------------- fetching
def _json(url: str, timeout_s: int = 40):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout_s) as r:
        return json.load(r)


def pfam_counts(accession: str) -> dict[str, int]:
    """`{pfam accession: copies}` for one UniProt entry, via InterPro."""
    try:
        data = _json(f"{INTERPRO}/entry/pfam/protein/uniprot/{accession}/"
                     f"?page_size=100")
    except Exception:
        return {}
    out: dict[str, int] = {}
    for entry in data.get("results", []):
        acc = entry["metadata"]["accession"]
        n = 0
        for prot in entry.get("proteins", []):
            n += len(prot.get("entry_protein_locations", []) or [])
        out[acc] = max(1, n)
    return out


def resolve_symbol(symbol: str, taxon: int = 9606) -> str:
    """Gene symbol → reviewed UniProt accession (for building panels only)."""
    # taxonomy_id, not organism_id: the latter is the exact taxonomy node
    # and misses every strain-level prokaryotic entry (see scripts/s0_lib.py).
    q = f"gene_exact:{symbol} AND taxonomy_id:{taxon} AND reviewed:true"
    try:
        d = _json(f"{UNIPROT}/search?query={urllib.parse.quote(q)}"
                  f"&fields=accession&size=1")
    except Exception:
        return ""
    rs = d.get("results", [])
    return rs[0]["primaryAccession"] if rs else ""


def entry_meta(accession: str) -> dict:
    try:
        d = _json(f"{UNIPROT}/{accession}.json?fields=accession,gene_primary,"
                  f"organism_name,length,ft_transmem,fragment")
    except Exception:
        return {}
    genes = d.get("genes") or [{}]
    tm = sum(1 for f in d.get("features", [])
             if f.get("type") == "Transmembrane")
    return {
        "gene": (genes[0].get("geneName") or {}).get("value", ""),
        "species": (d.get("organism") or {}).get("scientificName", ""),
        "length": (d.get("sequence") or {}).get("length"),
        # 0 is a measurement (the entry was fetched with its TM features and
        # has none); only a failed fetch is unknown. `tm or None` turned every
        # soluble protein into "unknown" — found in S2b, where the positive
        # AChBP test (0 TM helices) then called nothing.
        "tm_count": tm,
        # UniProt's flag also carries "Precursor"; only "Fragment(s)" is one.
        "fragment": "Fragment" in ((d.get("proteinDescription") or {}).get("flag") or ""),
    }


def build_query(accession: str, sleep_s: float = 0.1) -> ChannelQuery:
    meta = entry_meta(accession)
    time.sleep(sleep_s)
    counts = pfam_counts(accession)
    time.sleep(sleep_s)
    seq = fetch_uniprot_sequence(accession)
    return ChannelQuery(
        accession=accession, sequence=seq, gene_symbol=meta.get("gene", ""),
        species=meta.get("species", ""), pfam_counts=counts,
        tm_count=meta.get("tm_count"), length_aa=meta.get("length"),
        fragment=meta.get("fragment"),
    )



#: UniProt accession pattern: 6 or 10 characters, letter-first, digit-bearing.
def _looks_like_accession(token: str) -> bool:
    t = token.strip().upper()
    if len(t) not in (6, 10) or not t[0].isalpha():
        return False
    return t[1].isdigit() and any(c.isdigit() for c in t[2:])


# ---------------------------------------------------------- commands
def run_classify(project_root: Path, targets: list[str],
                 save_results: bool = False, label: str = "",
                 use_reference: bool = True, expected: dict | None = None,
                 quiet: bool = False) -> dict:
    """Classify UniProt accessions (or human gene symbols) into the catalogue."""
    accs: list[str] = []
    alias: dict[str, str] = {}          # accession -> the token the caller used
    for t in targets:
        t = t.strip()
        if not t:
            continue
        if _looks_like_accession(t):
            accs.append(t)
            alias[t] = t
            continue
        symbol, _, organism = t.partition("@")
        a = resolve_symbol(symbol, int(organism) if organism else 9606)
        if not a:
            print(f"[classify] could not resolve symbol {t!r}")
            continue
        print(f"[classify] {t} → {a}")
        accs.append(a)
        alias[a] = t

    # `expected` may be keyed by symbol or by accession; normalise to accession
    # so a benchmark panel can be written in gene symbols and stay readable.
    if expected:
        expected = {acc: expected.get(acc, expected.get(alias.get(acc, ""), ""))
                    for acc in accs
                    if expected.get(acc) or expected.get(alias.get(acc, ""))}

    print(f"[classify] fetching evidence for {len(accs)} protein(s)…")
    queries = []
    for i, a in enumerate(accs, 1):
        q = build_query(a)
        print(f"[classify]   [{i}/{len(accs)}] {a} {q.gene_symbol or '?':10s} "
              f"{q.length_aa or '?'} aa  {len(q.pfam_counts)} pfam")
        queries.append(q)

    refs = None
    nav_ref = ""
    if use_reference:
        cache = project_root / "results" / "s0_baseline" / "reference_panel.fasta"
        if cache.exists():
            refs = ReferenceSet.from_fasta(cache)
            print(f"[classify] reference panel: {len(refs.sequences)} exemplars "
                  f"({cache})")
        else:
            print(f"[classify] no reference panel at {cache} — run "
                  f"scripts/s0_catalogue_verify.py first; reference tier off")
        nav_ref = fetch_uniprot_sequence(FOUR_REPEAT_ANCHOR.reference_uniprot)

    calls = classify_all(queries, refs, nav_ref,
                         on_progress=None if quiet else
                         (lambda m: print(f"[classify] {m}")))
    for c in calls:
        print()
        print(c.text())

    counts = summary_counts(calls, expected)
    print()
    print("[classify] " + ", ".join(f"{k}={v}" for k, v in counts.items()))

    out: dict = {"counts": counts, "n": len(calls)}
    if save_results:
        from .utils.results_writer import make_bundle_dir
        d = make_bundle_dir(project_root / "results", label or "classify")
        write_tsv(d / "calls.tsv", CALLS_HEADER, calls_table(calls, expected))
        (d / "calls.json").write_text(json.dumps(
            [c.to_dict() for c in calls], indent=2))
        if expected:
            write_tsv(d / "confusion.tsv", ("expected", "called", "n"),
                      confusion_table(calls, expected))
            write_tsv(d / "recall.tsv",
                      ("family", "name", "n", "correct", "recall",
                       "decisive_tiers", "confidence", "mis_calls"),
                      recall_table(calls, expected))
        print(f"[classify] saved → {d}")
        out["dir"] = str(d)
    return out


def run_catalogue(scope_key: str = "all", show_families: bool = False) -> dict:
    """Print the catalogue, or one scope of it."""
    sc = scope_for(scope_key)
    print(sc.describe())
    print()
    for k in sc.family_keys:
        f = CATALOGUE[k]
        print(f"  {f.key:24s} {f.superfamily:14s} {len(f.human_genes):3d} genes  "
              f"{f.status.value:18s} {f.name}")
        if show_families:
            print(f"      architecture: {f.architecture() or '(none declared)'}")
            print(f"      exemplars   : "
                  f"{', '.join(e.label for e in f.exemplars) or '(none)'}")
            if f.notes:
                print(f"      note        : {f.notes[:150]}")
    print()
    if sc.sister_family_keys:
        print("  sister / decoy families in scope: "
              + ", ".join(sc.sister_family_keys[:12])
              + (" …" if len(sc.sister_family_keys) > 12 else ""))
    print("  " + ", ".join(f"{k}={v}" for k, v in stats().items()))
    return {"scope": sc.key, "families": len(sc.family_keys)}


def run_phylo(project_root: Path, unit: str, tier: int = 0,
              bootstrap: int = 1000, linsi: bool = False) -> dict:
    """Build a tier-1 (family) or tier-2 (superfamily, pore-module) tree."""
    if tier == 0:
        tier = 2 if unit in SUPERFAMILIES else 1
    out_dir = project_root / "results" / "phylogeny" / f"tier{tier}_{unit}"
    keys = ([unit] if tier == 1 else
            [f.key for f in CATALOGUE.values() if f.superfamily == unit])
    seqs: list[tuple[str, str]] = []
    for k in keys:
        for e in CATALOGUE[k].exemplars:
            if not e.uniprot:
                continue
            s = fetch_uniprot_sequence(e.uniprot)
            if s:
                seqs.append((e.label, s))
            time.sleep(0.1)
    print(f"[phylo] tier {tier} {unit}: {len(seqs)} exemplar sequence(s)")
    if tier == 2 and unit in SUPERFAMILIES and needs_modules(unit):
        print(f"[phylo] NOTE: {unit} needs pore-module extraction for a valid "
              f"tier-2 tree, and this command aligns full-length exemplars. "
              f"Treat the result as a quick look, not as the S8 tree — "
              f"module extraction is S6's job (src/phylo/modules.py).")
    try:
        run = (build_tier1(unit, seqs, out_dir, bootstrap=bootstrap, linsi=linsi)
               if tier == 1 else
               build_tier2(unit, seqs, out_dir, module_based=False,
                           bootstrap=bootstrap, linsi=linsi))
    except NotAlignable as exc:
        print(f"[phylo] REFUSED — {exc}")
        return {"refused": str(exc)}
    write_run(out_dir, run)
    print(f"[phylo] {run.label()}: {run.n_sequences} seqs, "
          f"{run.elapsed_s:.1f}s{'; ' + run.note if run.note else ''}")
    if run.tree_path:
        print(f"[phylo] tree → {run.tree_path}")
    return run.to_dict()
