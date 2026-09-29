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

# ------------------------------------------------------------ S4 proteome scope
print("proteome scope (S4)")
from scripts import s4_proteome_lib as _s4                    # noqa: E402

_c = lambda upid, rev, busco, genes: {"upid": upid, "reviewed": rev,  # noqa: E731
                                      "busco_complete_pct": busco, "gene_count": genes}
check("selection: most Swiss-Prot entries wins over a better BUSCO (the model strain)",
      _s4.select_proteome([_c("UP2", 199, 95.0, 5357), _c("UP1", 6067, 90.0, 6066)])["upid"], "UP1")
check("selection is total: an exact tie falls to the lower UPID",
      _s4.select_proteome([_c("UP9", 5, None, 9), _c("UP3", 5, None, 9)])["upid"], "UP3")
check("no candidate in the release → no proteome (the row becomes genome_only)",
      _s4.select_proteome([]), None)
check("a failed MD5 is never waived, even for a reissued file",
      _s4.file_verdict(False, True, True, 10, 10), "FAIL:md5")
check("counts off the README without a reissue date → FAIL",
      _s4.file_verdict(True, False, False, 10, 10), "FAIL:counts")
check("a reissued file passes only if its count is the declared gene count (D35)",
      (_s4.file_verdict(True, False, True, 10, 10), _s4.file_verdict(True, False, True, 11, 10)),
      ("reissued", "FAIL:counts"))
check("README statistics parse (upid → counts)",
      _s4.parse_readme("Release 2026_03, x\nProteome_ID\tTax_ID\tOSCODE\tSUPERREGNUM\t#(1)\t#(2)"
                       "\t#(3)\tSpecies Name\nUP000005640\t9606\tHUMAN\teukaryota\t1\t2\t3\tHomo"
                       )[1]["UP000005640"]["n_gene2acc"], 3)
_a = lambda acc, ann, lvl, n50: {"accession": acc, "annotated": ann, "level": lvl,  # noqa: E731
                                 "scaffold_n50": n50}
check("assembly rank: annotated beats a better unannotated assembly",
      max([_a("GCA_1", False, "Chromosome", 10**8), _a("GCA_2", True, "Scaffold", 10**5)],
          key=_s4.assembly_rank)["accession"], "GCA_2")

# ------------------------------------------------------------ S3b kill criterion
print("\njackhmmer kill criterion (S3b, D10)")
from scripts import s3b_kill as _k                            # noqa: E402
_rd = lambda i, inc: {"round": i, "new_targets": 0, "included": inc}  # noqa: E731
_fam = {f"o{i}": "nav" for i in range(50)} | {f"x{i}": "cav" for i in range(50)}
_own = [f"o{i}" for i in range(50)]
_x = [f"x{i}" for i in range(50)]
check("K1: a high but steady other-family share at round 1 is homology, not drift",
      _k.evaluate([_rd(1, _own[:20] + _x[:20]), _rd(2, _own[:25] + _x[:25])],
                  _fam, "nav", converged=True)["verdict"], "clean")
_ev = _k.evaluate([_rd(1, _own[:20]), _rd(2, _own[:20] + _x[:10])], _fam, "nav", True)
check("K1: a rise in other-family share kills the run and keeps the rounds before it",
      (_ev["rule"], _ev["accepted_rounds"]), ("K1", 1))
check("K2: a ten-fold jump in the included set kills the run",
      _k.evaluate([_rd(1, _own[:20]), _rd(2, _own[:20] + [f"u{i}" for i in range(300)])],
                  _fam, "nav", True)["rule"], "K2")
check("K3: converging on the last allowed round is not a ceiling kill",
      _k.evaluate([_rd(i, _own[:5]) for i in range(1, _k.MAX_ITER + 1)], _fam, "nav",
                  converged=True)["verdict"], "clean")
check("K3: the ceiling without convergence is a kill",
      _k.evaluate([_rd(i, _own[:5]) for i in range(1, _k.MAX_ITER + 1)], _fam, "nav",
                  converged=False)["rule"], "K3")
check("uncalled targets (module / no hit) never trip K1",
      _k.evaluate([_rd(1, _own[:5]), _rd(2, _own[:5] + [f"u{i}" for i in range(40)])],
                  _fam, "nav", True)["verdict"], "clean")
check("accepted targets are the last accepted round's set",
      _k.accepted_targets([_rd(1, ["a"]), _rd(2, ["a", "b"]), _rd(3, ["a", "b", "c"])], 2),
      {"a", "b"})

# ------------------------------------------------------------ S5 genomic sweep
print("\ngenomic sweep (S5)")
import tempfile as _tf                                        # noqa: E402
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import s5_lib as _l, s5_rescue as _r, s5_ledger as _g  # noqa: E402
from scripts.s5_verdict import d4_bar as _bar, verdict as _v  # noqa: E402
_gff = ("##PAF\tx\n##STA\tMKV\n"
        "c1\tminiprot\tmRNA\t100\t400\t50\t+\t.\tID=MP1;Identity=0.5;Target=kv|A 1 100\n"
        "c1\tminiprot\tCDS\t100\t200\t.\t+\t0\tParent=MP1;Target=kv|A 1 34\n"
        "c1\tminiprot\tCDS\t300\t400\t.\t+\t0\tParent=MP1;Target=kv|A 35 68\n"
        "##STA\tMSELF\n"
        "c1\tminiprot\tmRNA\t100\t400\t90\t+\t.\tID=MP2;Identity=0.9;Target=kv|B 1 100\n"
        "c1\tminiprot\tCDS\t100\t400\t.\t+\t0\tParent=MP2;Target=kv|B 1 100\n")
with _tf.NamedTemporaryFile("w", suffix=".gff", delete=False) as _fh:
    _fh.write(_gff)
_meta = {"kv|A": {"family": "kv", "species": "Other", "length": "100"},
         "kv|B": {"family": "kv", "species": "Self", "length": "100"}}
_al, _dr = _l.parse_miniprot_gff(Path(_fh.name), _meta, exclude_species="Self")
check("miniprot: a genome's own species' baits are dropped (genome-scale LOO)",
      ([a.bait for a in _al], _dr), (["kv|A"], 1))
check("miniprot: ##STA binds to the next mRNA; coverage is the CDS query union",
      (_al[0].translation, _al[0].aligned_aa, _al[0].max_intron), ("MKV", 68, 99))
_mk = lambda s, e, f: _l.Aln("c", s, e, "+", 10, 0.5, f"{f}|x", f, 1, 10, 10)  # noqa: E731
check("loci: overlapping alignments are one locus; a neighbour 1 kb away is not",
      len(_l.cluster_loci([_mk(1, 500, "a"), _mk(400, 900, "b"), _mk(1900, 2500, "c")])), 2)
check("rescue: an HSP inside any recorded locus is never a trace",
      [h["start"] for h in _r.outside_loci(
          [{"contig": "c", "start": 50, "end": 60}, {"contig": "c", "start": 5000, "end": 5100}],
          [{"contig": "c", "start": 1, "end": 100}])], [5000])
_L = lambda call, fam, bf, edge="0": {"p_call": call, "p_family": fam, "p_confidence":  # noqa: E731
                                     "high" if call == "family" else "none",
                                     "bait_family": bf, "at_edge": edge,
                                     "longest_n_run": "0", "locus": "l"}
check("genome status: the profile call decides, not the bait that found the locus",
      (_g.genome_status("kv", [_L("family", "kctd", "kv")])["genome"],
       _g.genome_status("kv", [_L("module", "", "kv", edge="1")])["genome"]),
      ("no_locus", "gap"))
_cell = lambda **k: {"proteome_records": "0", "genome": "no_locus", "n_traces": 0,  # noqa: E731
                     "informative": 1, "matched": 1, "species": "S", "family": "kv",
                     "bar_bp": 5000, **k}
_ctl = {"S": {"matched_detection": 0.95}}
check("verdict: absent needs a matched bait, a passing control and the D4 bar",
      (_v(_cell(), {"n50": "10000"}, _ctl),
       _v(_cell(), {"n50": "1000"}, _ctl),
       _v(_cell(matched=0), {"n50": "10000"}, _ctl),
       _v(_cell(), {"n50": "10000"}, {"S": {"matched_detection": 0.5}}),
       _v(_cell(n_traces=1), {"n50": "10000"}, _ctl)),
      ("absent", "absent_below_bar", "unmatched", "uncontrolled", "trace"))
_sp = {"kv": {"pooled": 50000, "groups": {"invertebrate": (8000, 3), "plant": (2000, 2)}}}
check("D4 bar (D38): prokaryote = CDS; own group with >= 3 genes; else pooled; else none",
      (_bar("kv", "prokaryote", _sp, 600), _bar("kv", "invertebrate", _sp, 600),
       _bar("kv", "plant", _sp, 600), _bar("nav", "plant", _sp, 600)),
      ((1800, "cds"), (8000, "group"), (50000, "pooled"), (None, "unmeasured")))
check("verdict: an unmeasured bar never reads as absent",
      _v(_cell(bar_bp=""), {"n50": "10000"}, _ctl), "absent_bar_unmeasured")

check("verdict: a proteome miss needs a high-confidence intact locus; genome-only is presence",
      (_v(_cell(genome="found", n_strong=1), {"n50": "1"}, _ctl),
       _v(_cell(genome="found", n_strong=0), {"n50": "1"}, _ctl),
       _v(_cell(genome="found", n_strong=0), {"n50": "1", "status": "genome_only"}, _ctl)),
      ("genome_found", "genome_weak", "genome_present"))
from scripts.s5_annotation import GeneIndex as _GI           # noqa: E402
_gi = _GI([{"contig": "c", "start": 1, "end": 100000, "strand": "+", "gene_id": "g",
            "name": "", "biotype": "protein_coding", "max_intron": 0, "n_exons": 2,
            "exons": [(1, 500), (99000, 100000)]}])
check("annotation confirms a model only where >= half its CDS lies on the gene's exons",
      (_gi.best_overlap("c", 100, 99500, "+", [(100, 400), (99100, 99500)]) is not None,
       _gi.best_overlap("c", 450, 60000, "+", [(450, 650), (59700, 60000)]) is None),
      (True, True))

# --- S6: alignment sets and pore modules (D39, D40) -------------------------
from src.catalogue import SUPERFAMILIES as _SFS                # noqa: E402
from src.phylo.forest import needs_modules as _nm              # noqa: E402
check("module_rule is declared exactly where a tier-2 unit needs modules (D40)",
      [k for k, sf in _SFS.items() if sf.alignable and _nm(k) != bool(sf.module_rule)],
      [])
from scripts.s6_lib import verdict as _v39                     # noqa: E402
_r = dict(v3_status="channel", p_family="nav", v3_family="nav",
          p_confidence="high", call_intact="", source="proteome")
check("D39: high profile call in; medium, other-family, broken genome frame out",
      (_v39(_r), _v39({**_r, "p_confidence": "medium"}),
       _v39({**_r, "p_family": "cav"}),
       _v39({**_r, "v3_status": "genome_locus", "call_intact": "0"}),
       _v39({**_r, "v3_status": "candidate"})),
      ("include", "profile_medium", "profile_other_family",
       "genome_not_intact", "not_called"))
from scripts.s6_module_refs import spans_for as _spans         # noqa: E402
_kv = [(165, 186), (221, 242), (254, 274), (288, 308), (324, 345), (387, 415)]
check("pore_loop: a loop annotated in two pieces is one module (S5 start → S6 end)",
      _spans("pore_loop", _kv, [(360, 371), (372, 379)]), [(324, 415)])
check("pore_loop: one module per loop; tm_span: first TM → last TM",
      (len(_spans("pore_loop", _kv * 1 + [(500, 520), (560, 580)],
                  [(360, 371), (530, 540)])),
       _spans("tm_span", [(30, 50), (80, 100), (150, 170), (200, 220)], [])),
      (2, [(30, 220)]))
from scripts.s6_project import a2m_map as _a2m, cut as _cut    # noqa: E402
_st, _res = _a2m("mkAC-Dx")     # states A C - D; inserts m k (before) and x (after)
check("A2M map: states, deletions and inserts; a cut keeps inserts inside the span",
      (_st[1:], _res, _cut(_st, _res, "MKACDX", 1, 4)),
      ([3, 4, None, 5], [0.5, 0.5, 1, 2, 4, 4.5], (3, 5, 0.75)))

# --- S7: tier-1 rooting and tree reading (D41) ------------------------------
from scripts.s7_newick import (ingroup_support as _isup, parse as _nwk,  # noqa: E402
                               split_support as _ssup)
from scripts.s7_lib import apply_mask as _mask, metrics as _met  # noqa: E402
from scripts.s7_trees import root_rule as _rr                  # noqa: E402
check("D41 rooting comes from the catalogue: sister family, else unrooted and said so",
      (_rr("kv_shaker"), _rr("kcsa_prok"), _rr("piezo"), _rr("ryr"), _rr("itpr")),
      (("rooted", "kcsa_prok"), ("unrooted:is_superfamily_outgroup", ""),
       ("unrooted:no_declared_outgroup", ""), ("rooted", "itpr"),
       ("unrooted:is_superfamily_outgroup", "")))
_t = _nwk("((a:1,b:1)90:1,(c:1,d:1)99:1,(OG_x:1,OG_y:1)80:1);")
check("split test: an outgroup clade is found with its UFBoot; a non-clade is not",
      (_ssup(_t, {"OG_x", "OG_y"}), _ssup(_t, {"a", "c"})[0],
       _ssup(_t, {"a", "b", "c", "d"})), ((True, 80.0), False, (True, 80.0)))
check("ingroup support counts only edges inside the family",
      _isup(_t, {"OG_x", "OG_y"})["internal_edges"], 2)
check("the trim mask is one column set for every row; informative = 2 states x 2",
      (_mask([("a", "AC-D"), ("b", "AG-E")], [0, 1, 3]),
       _met([("a", "AAC"), ("b", "AAC"), ("c", "GAT"), ("d", "GA-")], [0, 1, 2])["informative"]),
      ([("a", "ACD"), ("b", "AGE")], 1))

# --- S20: auxiliary subunits ------------------------------------------------
from scripts.s20_lib import category as _cat                    # noqa: E402
from scripts.s20_aux import _components as _comp                # noqa: E402
check("S20 category: census first; a channel_associated gene is auxiliary, never pore",
      (_cat("CHANNEL_CONTESTED", True), _cat("CHANNEL_ASSOCIATED", False),
       _cat("OUT_OF_SCOPE", False)), ("pore_census", "auxiliary", "out_of_scope"))
check("S20 homology groups are connected components (no edge = separate group)",
      _comp(["a", "b", "c", "d"], {("a", "b"), ("b", "c")}), [["a", "b", "c"], ["d"]])

print()
if FAILURES:
    print(f"{len(FAILURES)} invariant(s) FAILED: {', '.join(FAILURES)}")
    sys.exit(1)
print("all invariants hold")
