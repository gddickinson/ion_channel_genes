"""The architecture tier — domain-composition rules, and the hazard rules that override them.

A rule fires on the *multiset* of domain accessions a protein carries: which
signatures are present, and in how many copies. That is enough to make a
confident call for maybe two thirds of the catalogue and is provably not
enough for the rest, which is the point of keeping the rules in a table that
records their own priority.

Rules come from two places.

**Derived rules** are generated from the catalogue itself: a family's
`FAMILY`- and `SUBFAMILY`-level signatures become its requirement, and a
family with only `SUPERFAMILY` signatures generates a superfamily-level rule
that deliberately declines to name a family. Nothing is hand-maintained
twice — edit the catalogue and the rules follow.

**Hazard rules** are hand-written, sit at higher priority, and exist because
a hazard is only closed by a test. Each one names its `H` number, so
`scripts/s1_benchmark.py` can report per-hazard precision and recall, and so
that deleting a rule is visibly deleting a test.

Ties are not resolved. Two rules of equal priority naming different families
produce an `ambiguous` result carrying both, which the classifier hands to
the reference and motif tiers. Silently picking the first match would turn
the whole hazard registry into decoration.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..catalogue import CATALOGUE, Level, census_families

# Priorities. Higher wins; equal-and-different is ambiguity, not a coin toss.
P_HAZARD = 100
P_FAMILY_STRONG = 50    # two or more family-level signatures
P_FAMILY = 40           # one family-level signature
P_SUPERFAMILY = 10      # superfamily evidence only — no family call



@dataclass(frozen=True)
class ArchitectureRule:
    rid: str
    target: str                       # catalogue family key ("" = superfamily only)
    superfamily: str
    require: tuple[str, ...] = ()     # accessions that must be present
    min_copies: tuple[tuple[str, int], ...] = ()   # (accession, minimum copies)
    forbid: tuple[str, ...] = ()      # accessions that must be absent
    priority: int = P_FAMILY
    hazard: str = ""
    rationale: str = ""
    # Whole-sequence measurements a rule may require. A rule that states one
    # does not fire when the measurement is unknown: a positive test needs
    # the measurement, and "unknown" is not "zero".
    max_tm: int | None = None         # predicted TM helices ≤ this
    max_length: int = 0               # sequence length ≤ this (0 = no bound)
    min_length: int = 0               # sequence length ≥ this (0 = no bound)
    complete_only: bool = False       # UniProt does not flag it a fragment

    def matches(self, counts: dict[str, int], tm_count: int | None = None,
                length_aa: int | None = None,
                fragment: bool | None = None) -> bool:
        if any(a not in counts for a in self.require):
            return False
        if any(counts.get(a, 0) < n for a, n in self.min_copies):
            return False
        if any(a in counts for a in self.forbid):
            return False
        if self.max_tm is not None and (tm_count is None or tm_count > self.max_tm):
            return False
        if self.max_length and (not length_aa or length_aa > self.max_length):
            return False
        if self.min_length and (not length_aa or length_aa < self.min_length):
            return False
        if self.complete_only and fragment is not False:
            return False
        return bool(self.require or self.min_copies)

    def describe(self) -> str:
        bits = []
        if self.require:
            bits.append("+".join(self.require))
        for a, n in self.min_copies:
            bits.append(f"{a}×≥{n}")
        if self.forbid:
            bits.append("without " + "/".join(self.forbid))
        if self.max_tm is not None:
            bits.append(f"TM≤{self.max_tm}")
        if self.min_length or self.max_length:
            bits.append(f"{self.min_length or 0}–{self.max_length or '∞'} aa")
        if self.complete_only:
            bits.append("complete")
        return " ".join(bits)


# H12: the KCTD C-terminal domains, read from the catalogue (the one place a
# Pfam id is declared) — every FAMILY-level signature of nonchannel_kctd.
KCTD_MARKERS = tuple((sg.accession, sg.name)
                     for sg in CATALOGUE["nonchannel_kctd"].signatures
                     if sg.level is Level.FAMILY)


# --------------------------------------------------------- hazard rules
HAZARD_RULES: list[ArchitectureRule] = [
    # H2 — a pentameric channel needs the TM region, not just the clamshell LBD.
    ArchitectureRule("H2-cysloop", "", "cysloop", ("PF02931", "PF02932"),
                     priority=P_HAZARD, hazard="H2",
                     rationale="LBD + TM region ⇒ Cys-loop channel; family "
                               "within the superfamily needs the reference tier"),
    # H2 has no positive AChBP test at this tier (S2b, measured). The absence
    # rule (LBD without PF02932 ⇒ AChBP) called 6,633 receptor-length
    # proteins AChBP. Its best positive replacement — a complete, soluble
    # (0 TM helices), AChBP-sized LBD — agreed with the S3a profiles on 53 of
    # 1,263 records where both called, and ~75 % of its 2,164 calls were in
    # Ecdysozoa and Chordata, which have no known AChBP: it recognises
    # truncated receptor gene models as readily as AChBP. So an LBD without
    # the TM region is superfamily evidence, and AChBP is called by a
    # sequence-level tier (reference identity, or the S3a profile margin).
    ArchitectureRule("H2-lbd", "", "cysloop", ("PF02931",),
                     priority=P_HAZARD, hazard="H2",
                     rationale="a Cys-loop LBD without the TM region: AChBP, "
                               "a receptor lacking its TM annotation, or a "
                               "fragment — superfamily only"),
    # H3 — the iGluR clamshell is also the class C GPCR ligand-binding domain.
    ArchitectureRule("H3-gpcr", "nonchannel_class_c_gpcr", "iglur",
                     ("PF01094", "PF00003"), priority=P_HAZARD, hazard="H3",
                     rationale="clamshell + 7tm_3 ⇒ metabotropic receptor"),
    ArchitectureRule("H3-iglur", "", "iglur", ("PF00060",), forbid=("PF00003",),
                     priority=P_HAZARD, hazard="H3",
                     rationale="the pore region is what makes it ionotropic"),
    # H4 — ITPR versus RYR (inherited decision D14). ITPR has no domain RyR
    # lacks, so there is no positive ITPR test at this tier: the shared core
    # places a protein in the superfamily and the family call is left to a
    # sequence-level tier (reference identity, or S3a's profile margin).
    # S3a measured the old absence rule (PF08709 without RyR domains ⇒ ITPR)
    # calling N-terminal RyR fragments ITPR — 502 records the IP3R project
    # calls RYR. PF02026 and PF06459 are each positive RyR evidence: on
    # PF08709 carriers in census v2, 5,160 and 4,789 records respectively,
    # every one profile-called RYR.
    ArchitectureRule("H4-ryr", "ryr", "ca_release", ("PF08709", "PF02026"),
                     priority=P_HAZARD, hazard="H4",
                     rationale="the RyR repeat is present in RYR and absent "
                               "from ITPR"),
    ArchitectureRule("H4-ryr-tm", "ryr", "ca_release", ("PF08709", "PF06459"),
                     priority=P_HAZARD, hazard="H4",
                     rationale="RyR TM4-6 region with the shared N-terminal core"),
    ArchitectureRule("H4-core", "", "ca_release", ("PF08709",),
                     priority=P_HAZARD, hazard="H4",
                     rationale="the IP3R/RyR core both families carry: "
                               "superfamily only — ITPR has no positive "
                               "architectural test"),
    # H9 — a voltage-sensor domain is not evidence of a channel.
    ArchitectureRule("H9-vsp", "nonchannel_vsp", "hv", ("PF00520", "PF10409"),
                     priority=P_HAZARD, hazard="H9",
                     rationale="pore-module annotation + PTEN phosphatase ⇒ "
                               "voltage-sensing phosphatase, not a channel"),
    ArchitectureRule("H9-hv1", "hv1", "hv", ("PF16799",),
                     priority=P_HAZARD, hazard="H9",
                     rationale="the Hv1 C-terminal domain"),
    # H10 — MIR is shared with the O-mannosyltransferases.
    ArchitectureRule("H10-pomt", "nonchannel_pomt", "ca_release",
                     ("PF02815", "PF02366"), priority=P_HAZARD, hazard="H10",
                     rationale="MIR + mannosyltransferase ⇒ POMT"),
    # H11 — CFTR versus the sulfonylurea receptors. CFTR is called on its R
    # domain. There is no positive SUR test: ABCC8/9 carry only the generic
    # ABC domains, and TMD0 (PF24357) is shared with the long MRPs. The
    # absence rule it replaced (ABC architecture without the R domain ⇒ SUR,
    # removed in S2c) made 980 census calls and none was a SUR: 918 bacterial
    # ABC transporters with a cNMP and a C39 peptidase domain, and 62
    # eukaryotic gene models fusing ABC domains to other domains (none
    # SUR-shaped: 1,300–1,800 aa with ABC/TMD0 domains only).
    ArchitectureRule("H11-cftr", "cftr", "abc_channel", ("PF00664", "PF14396"),
                     priority=P_HAZARD, hazard="H11",
                     rationale="the R domain is CFTR's only architectural "
                               "difference from ABCC8"),
    # H12 — the Kv T1 domain is a generic BTB/POZ domain. KCTD is called on a
    # KCTD C-terminal domain (S2c, measured: 14,317 of the 14,320 carriers
    # the S3a profiles call are profile-KCTD). T1 without a pore module is
    # superfamily-only: a KCTD without a catalogued C-terminal domain and an
    # N-terminal Kv gene model look alike (the whole-protein shape test —
    # T1, 0 TM helices, complete — was contradicted by the profiles on 160
    # records). T1 *with* a pore module is left to the derived Kv rules.
    *[ArchitectureRule(f"H12-kctd-{mk}", "nonchannel_kctd", "ploop",
                       ("PF02214", mk), priority=P_HAZARD, hazard="H12",
                       rationale=f"T1/BTB + the KCTD C-terminal domain {name}")
      for mk, name in KCTD_MARKERS],
    ArchitectureRule("H12-t1", "", "ploop", ("PF02214",),
                     forbid=("PF00520", "PF07885"), priority=P_HAZARD, hazard="H12",
                     rationale="T1/BTB without a pore module or a KCTD "
                               "C-terminal domain: superfamily only"),
    # H13 — PF08016 / PF20519 cover TRPML, TRPP and polycystin-1. Each
    # family is called on a domain it carries and the others do not; the
    # shared channel domain alone is superfamily evidence. S3a measured the
    # old absence rule (polycystin domain without the mucolipin domain ⇒
    # TRPP) calling 2,989 polycystin-1-like proteins (median 2,263 aa) TRPP.
    *[ArchitectureRule(f"H13-pkd1-{mk}-{core}", "assoc_polycystin1", "ploop",
                       (core, mk), priority=P_HAZARD, hazard="H13",
                       rationale=f"polycystin-1 ectodomain marker {name}")
      for core in ("PF08016", "PF20519")
      for mk, name in (("PF01477", "PLAT"), ("PF02010", "REJ"), ("PF01825", "GPS"))],
    *[ArchitectureRule(f"H13-pkd1-repeats-{core}", "assoc_polycystin1", "ploop",
                       (core,), min_copies=(("PF00801", 5),), priority=P_HAZARD,
                       hazard="H13",
                       rationale="a long PKD-repeat ectodomain ⇒ polycystin-1")
      for core in ("PF08016", "PF20519")],
    ArchitectureRule("H13-trpml", "trpml", "ploop", ("PF08016", "PF21381"),
                     priority=P_HAZARD, hazard="H13",
                     rationale="the mucolipin extracytosolic domain"),
    ArchitectureRule("H13-trpp", "trpp", "ploop", ("PF20519", "PF18109"),
                     priority=P_HAZARD, hazard="H13",
                     rationale="the TRPP C-terminal domain PF18109 (measured: "
                               "1,358 census records, every one profile-called "
                               "TRPP)"),
    *[ArchitectureRule(f"H13-core-{core}", "", "ploop", (core,),
                       priority=P_HAZARD, hazard="H13",
                       rationale="polycystin channel domain shared by TRPP and "
                                 "polycystin-1: superfamily only")
      for core in ("PF08016", "PF20519")],
    # H16 — Kir and K2P differ by the copy number of one accession.
    ArchitectureRule("H16-k2p", "k2p", "ploop", min_copies=(("PF07885", 2),),
                     priority=P_HAZARD, hazard="H16",
                     rationale="two pore loops per subunit"),
    ArchitectureRule("H16-kir", "kir", "ploop", ("PF01007",),
                     priority=P_HAZARD, hazard="H16",
                     rationale="the inward-rectifier transmembrane domain"),
    # H1 — the four-repeat families are explicitly *not* callable here.
    ArchitectureRule("H1-fourrepeat", "", "ploop",
                     min_copies=(("PF00520", 4),), forbid=("PF06512", "PF08763"),
                     priority=P_HAZARD, hazard="H1",
                     rationale="four repeats with no family-specific domain: "
                               "Nav, Cav, NALCN and CatSper are "
                               "indistinguishable here — hand to the motif tier"),
]


# -------------------------------------------------------- derived rules
def derived_rules() -> list[ArchitectureRule]:
    """Rules generated from the catalogue's own signature declarations.

    Two rules per family where the catalogue supports them, and the split
    matters. A **family rule** requires only the `FAMILY`-level signatures,
    because those are the ones every member is expected to carry. A
    **subfamily rule** additionally requires the `SUBFAMILY`-level ones and
    sits at higher priority, so a protein that has them gets the more
    specific call and a protein that does not is still matched.

    Measured cost of getting this wrong: with `SUBFAMILY` signatures folded
    into the family requirement, `kv_shaker` required `PF02214` + `PF03521`
    (the Kv2-only model) and KCNA1 — a Kv1 — fell through to `kv_modifier`,
    whose rule needs `PF02214` alone. The correct architecture-tier answer
    for KCNA1 is *ambiguous between `kv_shaker` and `kv_modifier`*, which is
    what this generator now produces: the two families really are
    architecturally identical, and the reference tier is what separates them.

    Copy-number requirements are **not** generated here. They are brittle at
    family scale — measured: OTOP1 carries `PF03189`×3 and OTOP2 ×2 — so the
    two places where copy number is genuinely diagnostic (`PF07885`×2 for
    K2P, `PF00520`×4 for the four-repeat channels) are hazard rules instead,
    where the number is stated with the measurement behind it.
    """
    out: list[ArchitectureRule] = []
    for fam in census_families():
        family_sigs = tuple(s.accession for s in fam.signatures
                            if s.level is Level.FAMILY)
        subfamily_sigs = tuple(s.accession for s in fam.signatures
                               if s.level is Level.SUBFAMILY)
        if family_sigs:
            out.append(ArchitectureRule(
                f"cat-{fam.key}", fam.key, fam.superfamily, family_sigs,
                priority=(P_FAMILY_STRONG if len(family_sigs) > 1 else P_FAMILY),
                rationale=f"family-level signatures of {fam.name}"))
            if subfamily_sigs:
                out.append(ArchitectureRule(
                    f"cat-sub-{fam.key}", fam.key, fam.superfamily,
                    family_sigs + subfamily_sigs,
                    priority=P_FAMILY_STRONG + 5,
                    rationale=f"family + subfamily signatures of {fam.name}"))
            continue
        if subfamily_sigs:
            out.append(ArchitectureRule(
                f"cat-sub-{fam.key}", fam.key, fam.superfamily, subfamily_sigs,
                priority=P_FAMILY,
                rationale=f"subfamily signatures of {fam.name}"))
            continue
        shared = [s.accession for s in fam.signatures
                  if s.level is Level.SUPERFAMILY]
        if shared:
            out.append(ArchitectureRule(
                f"cat-sf-{fam.key}", "", fam.superfamily, tuple(shared[:1]),
                priority=P_SUPERFAMILY,
                rationale=(f"{fam.name} carries superfamily evidence only — "
                           f"no family call is possible from architecture")))
    return out


def all_rules() -> list[ArchitectureRule]:
    return HAZARD_RULES + derived_rules()


# ------------------------------------------------------------- matching
@dataclass
class RuleMatch:
    """The architecture tier's verdict."""
    family: str = ""
    superfamily: str = ""
    ambiguous: list[str] = field(default_factory=list)
    fired: list[ArchitectureRule] = field(default_factory=list)
    hazards: list[str] = field(default_factory=list)

    def called(self) -> bool:
        return bool(self.family)

    def summary(self) -> str:
        if self.family:
            head = f"family {self.family}"
        elif self.ambiguous:
            head = "ambiguous between " + ", ".join(sorted(self.ambiguous))
        elif self.superfamily:
            head = f"superfamily {self.superfamily} only"
        else:
            head = "no architecture match"
        rules = ", ".join(f"{r.rid}[{r.describe()}]" for r in self.fired[:4])
        return f"{head} — {rules}" if rules else head


def match_architecture(counts: dict[str, int],
                       rules: list[ArchitectureRule] | None = None,
                       tm_count: int | None = None,
                       length_aa: int | None = None,
                       fragment: bool | None = None) -> RuleMatch:
    """Match a `{accession: copies}` multiset (plus the whole-sequence TM
    count and length, for the rules that measure them) against the table."""
    rules = rules if rules is not None else all_rules()
    fired = [r for r in rules if r.matches(counts, tm_count, length_aa, fragment)]
    if not fired:
        return RuleMatch()
    best = max(r.priority for r in fired)
    top = [r for r in fired if r.priority == best]
    targets = {r.target for r in top if r.target}
    sfs = {r.superfamily for r in top if r.superfamily}
    res = RuleMatch(
        superfamily=(next(iter(sfs)) if len(sfs) == 1 else ""),
        fired=sorted(fired, key=lambda r: -r.priority),
        hazards=sorted({r.hazard for r in fired if r.hazard}),
    )
    if len(targets) == 1:
        res.family = next(iter(targets))
    elif len(targets) > 1:
        res.ambiguous = sorted(targets)
    return res


def topology_check(family_key: str, tm_count: int | None,
                   pore_loops: int | None = None) -> tuple[bool, str]:
    """Is an observed topology consistent with the family's declared one?

    Hazard **H8**: which Pfam pore model matched does *not* track the number
    of transmembrane helices — the 6TM SK channels are annotated with the 2TM
    model — so topology is taken from UniProt features or a prediction and
    compared here, never inferred from the domain hit.
    """
    fam = CATALOGUE.get(family_key)
    if fam is None or tm_count is None or not fam.tm_per_subunit:
        return True, "no topology expectation recorded"
    expected = fam.tm_per_subunit
    tol = max(2, expected // 4)
    ok = abs(tm_count - expected) <= tol
    detail = (f"{tm_count} TM observed vs {expected} expected "
              f"(±{tol}) for {fam.key}")
    if pore_loops is not None and fam.pore_loops_per_subunit:
        ok = ok and pore_loops == fam.pore_loops_per_subunit
        detail += f"; {pore_loops} pore loop(s) vs {fam.pore_loops_per_subunit}"
    return ok, detail
