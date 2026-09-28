"""Offline self-test — the invariants that must hold before any session starts.

No network, no external binaries, under a second. It checks the things that
break silently when the catalogue is edited: that a hazard rule still fires
on the architecture it was written for, that the phylogeny driver still
refuses the comparisons decision **D27** forbids, that a scope still derives
its sister panel from the hazard registry, and that the classifier still
never consults a gene symbol.

Run it after any catalogue edit, and as step 4a of the session protocol:

    python3 scripts/selftest.py

Exit code 0 means every invariant holds. Each failure prints what was
expected and what happened, because a self-test whose output is just a count
tells you nothing at the moment you need it.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.catalogue import CATALOGUE, HAZARDS, validate
from src.classify import ChannelQuery, classify, k_filter_hits
from src.classify.motifs import FILTER_CALLS
from src.classify.rules import match_architecture
from src.phylo import NotAlignable, build_tier2, refused_superfamilies
from src.utils.mafft import (alignment_stats, percent_identity,
                             project_positions)
from src.utils.scope import scope_for

FAILURES: list[str] = []


def check(name: str, got, want) -> None:
    if got == want:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name}\n         want {want!r}\n         got  {got!r}")
        FAILURES.append(name)


def check_true(name: str, cond, detail: str = "") -> None:
    check(name, bool(cond), True) if cond else (
        print(f"  FAIL {name}  {detail}"), FAILURES.append(name))


# ------------------------------------------------------------ catalogue
print("catalogue")
check("validate() is clean", validate(), [])
check("every hazard names known families",
      sorted({f for h in HAZARDS for f in h.families if f not in CATALOGUE}), [])
check("every hazard has a discriminator",
      [h.hid for h in HAZARDS if not h.discriminator.strip()], [])

# ------------------------------------------- architecture / hazard rules
print("\narchitecture rules — the hazard cases, by synthetic architecture")
CASES = [
    # (name, pfam counts, expected family, expected ambiguity)
    ("Nav (H1)",        {"PF00520": 4, "PF06512": 1, "PF11933": 1}, "nav", []),
    ("Cav3 (H1)",       {"PF00520": 4}, "", []),          # filter tier's job
    ("KCTD (H12)",      {"PF02214": 1}, "nonchannel_kctd", []),
    ("Kv1 (H12)",       {"PF00520": 1, "PF02214": 1}, "",
     ["kv_modifier", "kv_shaker"]),                        # honest ambiguity
    ("K2P (H16)",       {"PF07885": 2}, "k2p", []),
    ("Kir (H16)",       {"PF01007": 1}, "kir", []),
    ("AChBP (H2)",      {"PF02931": 1}, "nonchannel_achbp", []),
    ("mGluR (H3)",      {"PF01094": 1, "PF00003": 1}, "nonchannel_class_c_gpcr", []),
    ("RyR (H4)",        {"PF08709": 1, "PF02026": 4, "PF02815": 1}, "ryr", []),
    ("ITPR (H4)",       {"PF08709": 1, "PF02815": 1, "PF01365": 2}, "itpr", []),
    ("POMT (H10)",      {"PF02815": 1, "PF02366": 1}, "nonchannel_pomt", []),
    ("CFTR (H11)",      {"PF00664": 2, "PF00005": 2, "PF14396": 1}, "cftr", []),
    ("SUR1 (H11)",      {"PF00664": 2, "PF00005": 2}, "assoc_sur", []),
    ("VSP (H9)",        {"PF00520": 1, "PF10409": 1}, "nonchannel_vsp", []),
    ("PKD1 (H13)",      {"PF08016": 1, "PF00801": 15}, "assoc_polycystin1", []),
    ("TRPML (H13)",     {"PF08016": 1, "PF21381": 1}, "trpml", []),
]
for name, counts, want_family, want_ambig in CASES:
    m = match_architecture(counts)
    check(f"{name} family", m.family, want_family)
    if want_ambig:
        check(f"{name} ambiguity", sorted(m.ambiguous), sorted(want_ambig))

# -------------------------------------------------------------- symbols
print("\nthe gene symbol is never consulted (H15)")
counts = {"PF02214": 1}
a = classify(ChannelQuery("X1", gene_symbol="KCNA1", pfam_counts=counts))
b = classify(ChannelQuery("X1", gene_symbol="", pfam_counts=counts))
c = classify(ChannelQuery("X1", gene_symbol="TOTAL_NONSENSE", pfam_counts=counts))
check("same call with, without and with a wrong symbol",
      {a.family, b.family, c.family}, {"nonchannel_kctd"})

# ---------------------------------------------------------------- motifs
print("\nmotifs")
check("TxGYG found in a KcsA-like filter",
      [m for _, m in k_filter_hits("AAATVGYGDLYP")], ["TVGYG"])
check("no K+ filter in a random stretch", k_filter_hits("MKKLLIVDDE"), [])
check("every filter string maps to a catalogued family",
      sorted({f for f, _ in FILTER_CALLS.values() if f not in CATALOGUE}), [])
check("EEDD is Cav, not a non-call", FILTER_CALLS["EEDD"][0], "cav")

# ---------------------------------------------------------- alignment ops
print("\nalignment helpers")
check("projection maps reference positions onto the query",
      project_positions("AB-CD", "AXYCD", [1, 2, 3, 4]),
      {1: "A", 2: "X", 3: "C", 4: "D"})
check("covered-only identity ignores overhang",
      round(percent_identity("AAAA----", "AAAAKKKK", covered_only=True), 3), 1.0)
check("full-alignment identity does not",
      round(percent_identity("AAAA----", "AAAAKKKK", covered_only=False), 3), 0.5)
# S1's measured failure: a short sequence scattered through a long one keeps
# a high covered identity and near-total coverage *of itself*, and only
# coverage of the longer sequence exposes it.
_id, _cs, _cl, _n = alignment_stats("A-A-A-A-", "AXAXAXAX")
check("coverage of the shorter sequence does not catch a cherry-picked hit",
      round(_cs, 3), 1.0)
check("coverage of the longer sequence does", round(_cl, 3), 0.5)

# -------------------------------------------------------------- phylogeny
print("\nphylogeny — D27")
refused = {k for k, _ in refused_superfamilies()}
check_true("some superfamilies are marked non-alignable", refused, "none marked")
for sf in sorted(refused):
    try:
        build_tier2(sf, [], ROOT / "results" / "_selftest")
        print(f"  FAIL tier-2 on {sf} was ACCEPTED")
        FAILURES.append(f"tier2-{sf}")
    except NotAlignable:
        print(f"  ok   tier-2 on {sf} refused")

# ------------------------------------------------------------------ scope
print("\nscope")
itpr = scope_for("itpr")
check("the ITPR scope derives its sister panel from the hazards",
      sorted(itpr.sister_family_keys), ["nonchannel_pomt", "ryr"])
check_true("the ITPR sister panel carries RyR gene symbols",
           "RYR1" in itpr.sister_genes, itpr.sister_genes)
allsc = scope_for("all")
check_true("the whole-channelome scope enumerates more than PF00520",
           len(allsc.pfam_ids) > 1 and "PF00520" in allsc.pfam_ids,
           allsc.pfam_ids[:5])

# ------------------------------------------------------- census parsing (S2)
print("\ncensus parsing — D31")
from scripts.s2_lib import SHARDS, parse_record, pfam_dict   # noqa: E402
_rec = parse_record({
    "entryType": "UniProtKB reviewed (Swiss-Prot)", "primaryAccession": "Q14524",
    "organism": {"taxonId": 9606, "scientificName": "Homo sapiens",
                 "lineage": ["Eukaryota", "Metazoa", "Chordata"]},
    "features": [{"type": "Transmembrane"}] * 24 + [{"type": "Domain"}],
    "uniProtKBCrossReferences": [
        {"database": "Pfam", "id": "PF00520",
         "properties": [{"key": "MatchStatus", "value": "4"}]},
        {"database": "PDB", "id": "6UZ3", "properties": []},
        {"database": "Pfam", "id": "PF06512",
         "properties": [{"key": "MatchStatus", "value": "1"}]}],
    "sequence": {"value": "MANF"}})
check("UniProt MatchStatus is read as a copy number",
      pfam_dict(_rec["pfam"]), {"PF00520": 4, "PF06512": 1})
check("only Transmembrane features are counted", _rec["tm_count"], 24)
check("census shard keys are unique", len({k for k, _ in SHARDS}), len(SHARDS))

print()
if FAILURES:
    print(f"{len(FAILURES)} invariant(s) FAILED: {', '.join(FAILURES)}")
    sys.exit(1)
print("all invariants hold")
