# How a protein gets classified

The classifier answers one question — *which catalogued family is this?* —
through three independent tiers, and always says which one decided.

```
   architecture   what the domain composition allows        src/classify/rules.py
   motif          what the selectivity filter says          src/classify/motifs.py
   reference      what it is nearest to, with a margin      src/classify/reference.py
                                     ↓
                              ChannelCall
              family · superfamily · status · confidence · margin
              evidence[] from every tier · hazards[] · conflicts[]
```

## Why three tiers

Because each fails differently, and the families that defeat one are usually
caught by another.

*Architecture* is cheap, works from an InterPro record alone, and is decisive
for maybe two thirds of the catalogue. It cannot separate Nav from Cav from
NALCN from CatSper — they carry the same four copies of `PF00520` and nothing
else — and it cannot separate ANO1 from ANO6 at all.

*Motif* is the one test that works on a bare sequence from an unannotated
genome. The K⁺ signature T-x-G-Y-G needs no reference at all. The four-repeat
selectivity locus needs one alignment, and it resolves precisely the family
set architecture cannot.

*Reference identity* works whenever something related is in the panel and
fails exactly where a census needs it most: on the sequence with no close
relative. It is the last tier for that reason, not the first.

## The rules the tiers apply

**D25 — a signature is evidence only at the level where it is diagnostic.**
Every catalogued signature carries a `Level`. `PF00520` is `SUPERFAMILY`: it
says "P-loop channel" and never "this is a Nav". Measured: 20 catalogued
families carry it, including `nonchannel_vsp`, a phosphatase.

**D26 — the four-repeat families are separated by their filter, not their
domains.** Align the query to human Nav1.5 (`Q14524`) and read whatever falls
on its D372, E898, K1419, A1711:

| projection | call | note |
|-----------|------|------|
| `DEKA` | `nav` | the sodium locus |
| `EEEE` | `cav` | Cav1/Cav2 |
| `EEDD` | `cav` | Cav3 (T-type) — genuinely different, and correct |
| `EEKE` | `nalcn` | a lysine in repeat III, like Nav, in a leak channel |

Measured 2026-08-19 on six four-repeat channels: 6/6 correct at ~1.1 s per
query. `verify_anchor()` re-reads the reference's own filter before trusting
any projection, so a UniProt sequence-version bump cannot silently shift the
numbering.

**Ties are ambiguity, not a coin toss.** Rules are generated from the
catalogue at two levels — a family rule requiring the `FAMILY`-level
signatures, and a higher-priority rule that also requires the `SUBFAMILY`
ones — so a protein that lacks its family's subfamily marker still matches
the family rule. Where two families are genuinely indistinguishable by
architecture, the tier returns **ambiguous**, names both, and hands the
decision on. It never picks the first match. Measured: `kv_shaker` and
`kv_modifier` carry exactly the same two accessions, so KCNA1 comes back
ambiguous between them and is resolved by reference identity — which is the
truth about those two families, not a defect.

**D7 — a best hit is not a call** (inherited). The reference tier requires
the best family to beat the best *different* family by ≥ 0.10 identity.
Anything closer is an ambiguity, reported with the margin. When the
architecture tier has already named the candidates, the reference tier
prefilters over *only* those families — a global top-N can easily contain
none of them, and then the ambiguity survives the tier that was asked to
resolve it.

**Benchmarks classify leave-one-out (D29).** Ninety-three of the 97 proteins
in the S1 panel *are* catalogue exemplars — both come from the same curated
gene lists, so the overlap is near total. Scoring them against the panel that
contains them would return 100 % identity to themselves and measure
nothing, so `classify(..., leave_one_out=True)` drops the query's own
accession from the reference set. The overlap is reported in the S1
summary rather than quietly handled.

**Identity is scored over mutually covered columns.** Ion channels in one
superfamily differ in length by an order of magnitude (MscL 136 aa, RYR1
5,038 aa); full-alignment identity divides by the gaps. The parent project
measured RyR-vs-ITPR at 0.105 full-alignment and 0.249 covered-only — the
first is noise, the second is signal.

**H15 — the gene symbol is never consulted.** It is carried into the audit
trail for the reader and the call would be identical with it blank.
`KCNE*`, `CACNB*`, `CACNG*`, `SCN*B`, `CATSPERB` and `KCTD*` all sort inside
the channel symbol space and none is a pore; an unnamed gene model in a new
genome — the case a census exists for — has no symbol at all.

**D28 — a missing tool disables a tier, loudly.** No MAFFT means no motif
and no reference tier, and the call says so. There is no weaker fallback,
because a fallback that changes an answer without saying so is worse than a
missing answer.

## Confidence

| tier | meaning |
|------|---------|
| `gold` | two or more tiers name the same family, none disagrees |
| `silver` | one tier names it and nothing contradicts, or tiers agree with a conflict elsewhere |
| `bronze` | superfamily only, or a family named by a catalogue-derived rule alone |
| `unassigned` | no tier could say anything |

A topology mismatch (observed TM count against the family's declared one)
demotes gold to silver and is recorded as a conflict. It never overrides a
call, because hazard **H8** is precisely that Pfam's pore model does not
track TM count: the 6TM SK channels are annotated with the 2TM model.

## Conflicts are reported, not resolved by fiat

When tiers disagree the call follows a fixed precedence — hazard rule >
motif > reference > derived rule — and every disagreeing vote is written into
`conflicts`. A classifier that hides its disagreements cannot be audited, and
the disagreements are where the interesting biology is: a protein whose
architecture says one family and whose filter says another is either a
misannotation or something worth a paragraph.

## What the classifier does not do

- **It does not predict selectivity or gating.** Both are literature
  attributes of a family (**H14**).
- **It does not decide channel from non-channel by fold.** CFTR is an ABC
  transporter that is a channel; ABCC8 is the same architecture and is not.
  The discriminator is one accession (`PF14396`), and the general rule is
  that exclusion tests are positive tests on the pore.
- **It does not classify what it has never been shown.** A family with no
  exemplar and no diagnostic signature can only ever be reached by the
  reference tier, and the S1 coverage table reports which families those are.
