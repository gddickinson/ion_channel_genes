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
from src.classify.rules import HAZARD_RULES, match_architecture
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
import csv as _csv                                            # noqa: E402
from src.catalogue.registry import pore_signatures            # noqa: E402
_enum = ROOT / "results" / "census_v2" / "signature_counts.tsv"
if _enum.exists():
    check("the pore union is the search space census v2 enumerated (S2b drift)",
          sorted(s.accession for s in pore_signatures()),
          sorted(r["pfam"] for r in _csv.DictReader(_enum.open(), delimiter="\t")))
check("every hazard has a discriminator",
      [h.hid for h in HAZARDS if not h.discriminator.strip()], [])

# ------------------------------------------- architecture / hazard rules
print("\narchitecture rules — the hazard cases, by synthetic architecture")
CASES = [
    # (name, pfam counts, expected family, expected ambiguity[, tm, length])
    ("Nav (H1)",        {"PF00520": 4, "PF06512": 1, "PF11933": 1}, "nav", []),
    ("Cav3 (H1)",       {"PF00520": 4}, "", []),          # filter tier's job
    # S2c: KCTD on its C-terminal domain; T1 alone is superfamily-only.
    ("KCTD (H12)",      {"PF02214": 1, "PF31093": 1}, "nonchannel_kctd", []),
    ("T1 alone (H12)",  {"PF02214": 1}, "", []),
    ("Kv1 (H12)",       {"PF00520": 1, "PF02214": 1}, "",
     ["kv_modifier", "kv_shaker"]),                        # honest ambiguity
    ("K2P (H16)",       {"PF07885": 2}, "k2p", []),
    ("Kir (H16)",       {"PF01007": 1}, "kir", []),
    # S2b: no architecture-tier AChBP call exists (measured non-diagnostic).
    ("AChBP-shaped LBD (H2)",             {"PF02931": 1}, "", [], 0, 229, False),
    ("LBD + TM helices, no PF02932 (H2)", {"PF02931": 1}, "", [], 4, 383, False),
    ("mGluR (H3)",      {"PF01094": 1, "PF00003": 1}, "nonchannel_class_c_gpcr", []),
    ("RyR (H4)",        {"PF08709": 1, "PF02026": 4, "PF02815": 1}, "ryr", []),
    # S2b: ITPR has no domain RyR lacks — superfamily only at this tier.
    ("ITPR (H4)",       {"PF08709": 1, "PF02815": 1, "PF01365": 2}, "", []),
    ("RyR N-terminal fragment (H4)", {"PF08709": 1}, "", []),
    ("RyR TM4-6 (H4)",  {"PF08709": 1, "PF06459": 1}, "ryr", []),
    ("POMT (H10)",      {"PF02815": 1, "PF02366": 1}, "nonchannel_pomt", []),
    ("CFTR (H11)",      {"PF00664": 2, "PF00005": 2, "PF14396": 1}, "cftr", []),
    ("SUR1 (H11)",      {"PF00664": 2, "PF00005": 2}, "", []),   # S2c: no positive test
    ("VSP (H9)",        {"PF00520": 1, "PF10409": 1}, "nonchannel_vsp", []),
    ("PKD1 (H13)",      {"PF08016": 1, "PF00801": 15}, "assoc_polycystin1", []),
    ("TRPML (H13)",     {"PF08016": 1, "PF21381": 1}, "trpml", []),
    ("TRPP (H13)",      {"PF08016": 1, "PF20519": 1, "PF18109": 1}, "trpp", []),
    ("PKD1L-like, few PKD repeats (H13)",
     {"PF08016": 1, "PF20519": 1, "PF01477": 1, "PF02010": 1}, "assoc_polycystin1", []),
    ("polycystin channel domain alone (H13)", {"PF08016": 1, "PF20519": 1}, "", []),
]
for name, counts, want_family, want_ambig, *meas in CASES:
    tm, length, frag = (meas + [None, None, None])[:3]
    m = match_architecture(counts, tm_count=tm, length_aa=length, fragment=frag)
    check(f"{name} family", m.family, want_family)
    if want_ambig:
        check(f"{name} ambiguity", sorted(m.ambiguous), sorted(want_ambig))

# -------------------------------------------------------------- symbols
print("\nthe gene symbol is never consulted (H15)")
counts = {"PF02214": 1, "PF31093": 1}     # a KCTD by its C-terminal domain (S2c)
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
from scripts.s2_classify import classify_row                  # noqa: E402
# Whole-sequence measurements reach the rule engine intact (S2b bugs: a
# truthiness test turned 0 TM helices into "unknown", and UniProt's
# "Precursor" flag was read as "Fragment"). Probed with a synthetic rule.
from src.classify import rules as _rules                       # noqa: E402
_probe = _rules.ArchitectureRule("probe", "nonchannel_achbp", "cysloop", ("PF99999",),
                                 max_tm=0, min_length=150, max_length=300,
                                 complete_only=True, priority=_rules.P_HAZARD)
_saved = list(_rules.HAZARD_RULES)
_rules.HAZARD_RULES.append(_probe)
_row = {"accession": "X", "sequence": "M" * 229, "gene": "", "pfam": "PF99999:1",
        "tm_count": 0, "length": 229, "fragment": "Precursor"}
check("S2's row classifier passes 0 TM helices and 'Precursor' as measurements",
      classify_row(_row)["family"], "nonchannel_achbp")
check("… and 'Fragment' as a fragment",
      classify_row({**_row, "fragment": "Fragment"})["family"], "")
check("… and an unknown TM count as unknown, not zero",
      classify_row({**_row, "tm_count": ""})["family"], "")
_rules.HAZARD_RULES[:] = _saved

# ------------------------------------------------ profile assignment (S3a)
print("\nprofile assignment — D7 / D30 at HMM scale")
from scripts.s3_assign import assign_one, fam_superfamily   # noqa: E402
from scripts.s3_census_v3 import merge                      # noqa: E402
_sf = fam_superfamily()
# hit tuples: (score, profile, profile_coverage, target_coverage, evalue)
check("a clear winner over most of its profile is a family call",
      assign_one("t", [(900, "nav", .9, .9, 0), (500, "cav", .9, .9, 0)], _sf)["p_family"], "nav")
check("a high score over < 30 % of the profile is a shared module, not a call",
      assign_one("t", [(400, "hcn", .15, .9, 0)], _sf)["p_call"], "module")
check("two profiles within 10 % in one superfamily → superfamily_only",
      (lambda r: (r["p_call"], r["p_superfamily"]))(assign_one(
          "t", [(500, "nachr", .8, .9, 0), (470, "gabaa", .8, .9, 0)], _sf)),
      ("superfamily_only", "cysloop"))
check("two profiles within 10 % across superfamilies → ambiguous",
      assign_one("t", [(500, "nachr", .8, .9, 0), (470, "kir", .8, .9, 0)], _sf)["p_call"],
      "ambiguous")
check("the runner-up counts whatever its coverage (a close module is a competitor)",
      assign_one("t", [(500, "kv_shaker", .8, .9, 0), (480, "nonchannel_kctd", .2, .9, 0)],
                 _sf)["p_call"], "superfamily_only")
check("below the score floor nothing is called",
      assign_one("t", [(20, "mscl", .9, .9, 0)], _sf)["p_call"], "low_score")
_v2 = {"family": "cav", "superfamily": "ploop"}
check("S2 and profile agreeing → both",
      merge(_v2, {"p_call": "family", "p_family": "cav", "p_superfamily": "ploop"}, {})["v3_basis"], "both")
check("S2 and profile disagreeing → conflict, never a silent winner",
      merge(_v2, {"p_call": "family", "p_family": "nav", "p_superfamily": "ploop"}, {})["v3_basis"],
      "conflict")
check("a profile call outside S2's superfamily → conflict",
      merge({"family": "", "superfamily": "cysloop"},
            {"p_call": "family", "p_family": "kir", "p_superfamily": "ploop"}, {})["v3_basis"],
      "conflict")
check("no hazard rule makes a family call from an absence alone (S2b, S2c)",
      # a family-naming rule may forbid a domain only as a guard beside a
      # positive whole-sequence measurement (D33). Every hazard, since S2c.
      [r.rid for r in HAZARD_RULES if r.target and r.forbid and r.max_tm is None
       and not r.max_length], [])
check("a profile call inside S2's superfamily resolves it",
      merge({"family": "", "superfamily": "cysloop"},
            {"p_call": "family", "p_family": "nachr", "p_superfamily": "cysloop"}, {})["v3_family"],
      "nachr")

print()
if FAILURES:
    print(f"{len(FAILURES)} invariant(s) FAILED: {', '.join(FAILURES)}")
    sys.exit(1)
print("all invariants hold")
