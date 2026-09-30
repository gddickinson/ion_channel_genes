"""s15_contribution.py — S15: what each census method contributes, per superfamily (Q3).

    python3 scripts/s15_contribution.py      # → results/method_contribution/*.tsv

Two frames, because no single one is both independent and complete:

* **Human frame (independent truth).** The catalogue's 320 human census
  genes — curated from the literature, never from any census method. For
  each: enumerated by domain search (census v2 carries it), called to the
  right family by the S2 domain rules, called to the right family by the
  S3a profiles over the panel (S3b).
* **Panel frame (the final census).** Every census v4 row the S3a profiles
  call to a census family at high confidence (D32), genome loci only with an
  intact frame (D37) — S6's D39 rule, with contested families kept (they
  are census families). The cumulative curve per superfamily:
  1. **domain call** — S2's domain rules made the same family call (`both`);
  2. **domain enumeration** — the record carries an enumerated pore
     signature at all (in census v2);
  3. **+ profile** — every proteome row (the profile call is what defines
     the frame);
  4. **+ genome** — loci the proteomes do not contain.
  The frame is defined by the best instrument, so step 3 is 100 % of the
  proteome rows by construction; the curve measures what domain search
  misses *relative to it*, the upper bound the review's Q3 asks for.
  Domain-only finds (S2 family calls the profiles do not make) are counted
  beside the curve, not in it.

**Q3 verdict per superfamily, thresholds fixed before any row was read:**
domain enumeration reaches ≥ 95 % of the frame → `domain_search`; 50–95 % →
`domain_partial`; < 50 % → `profile_or_genome_only`. The same thresholds on
the domain *call* give the second column.
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s6_lib import ROOT, iter_census_v4  # noqa: E402
from src.catalogue import CATALOGUE  # noqa: E402

OUT = ROOT / "results" / "method_contribution"
CENSUS_STATUS = ("channel", "channel_contested", "genome_locus")
HIGH, LOW = 0.95, 0.50


def q3(frac: float) -> str:
    return ("domain_search" if frac >= HIGH else
            "domain_partial" if frac >= LOW else "profile_or_genome_only")


def in_frame(r: dict) -> bool:
    """High-confidence census-family profile call; intact frame for a locus."""
    fam = r["v3_family"]
    return (r["v3_status"] in CENSUS_STATUS and fam in CATALOGUE
            and CATALOGUE[fam].census_member()
            and r["p_family"] == fam and r["p_confidence"] == "high"
            and (r["v3_status"] != "genome_locus" or r["call_intact"] == "1"))


def step(r: dict) -> str:
    """The first method on the curve that finds this row."""
    if r["source"] == "genome":
        return "genome"
    if r["v3a_basis"] == "both":
        return "domain_call"
    if r["in_v2"] == "1":
        return "domain_enumeration"
    return "profile"


def human_frame() -> list[dict]:
    s2 = {r["gene"]: r for r in read_tsv(ROOT / "results" / "census_v2" / "human_recall.tsv")}
    s3 = {r["gene"]: r for r in read_tsv(ROOT / "results" / "panel_sweep" / "human_recall.tsv")}
    rows = []
    for g, r in sorted(s2.items()):
        fam = r["expected_family"]
        rows.append({"gene": g, "family": fam, "superfamily": CATALOGUE[fam].superfamily,
                     "enumerated": int(r["s2_status"] != "absent"),
                     "domain_call": int(r["correct"] == "yes"),
                     "profile_call": int(s3.get(g, {}).get("verdict") == "right_family"),
                     "s2_status": r["s2_status"]})
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    hum = human_frame()
    write_tsv(OUT / "human_genes.tsv", list(hum[0]), hum)
    hsf = defaultdict(list)
    for r in hum:
        hsf[r["superfamily"]].append(r)

    frame, domain_only = [], Counter()
    for r in iter_census_v4():
        fam = r["v3_family"]
        if (r["v3a_basis"] == "s2_only" and fam in CATALOGUE
                and CATALOGUE[fam].census_member()):
            domain_only[CATALOGUE[fam].superfamily] += 1
        if in_frame(r):
            frame.append({"family": fam, "superfamily": CATALOGUE[fam].superfamily,
                          "group": r["group"], "species": r["species"],
                          "source": r["source"], "step": step(r)})
    steps = ("domain_call", "domain_enumeration", "profile", "genome")

    def curve(rows: list[dict]) -> dict:
        c = Counter(x["step"] for x in rows)
        n = len(rows)
        cum, out = 0, {"n": n}
        for s in steps:
            cum += c[s]
            out[f"n_{s}"] = c[s]
            out[f"cum_{s}"] = round(cum / n, 4) if n else ""
        return out

    by_sf = defaultdict(list)
    for x in frame:
        by_sf[x["superfamily"]].append(x)
    sf_rows = []
    for sf, rows in sorted(by_sf.items(), key=lambda kv: -len(kv[1])):
        c = curve(rows)
        prot = [x for x in rows if x["source"] != "genome"]
        enum_frac = (sum(x["step"] in ("domain_call", "domain_enumeration") for x in prot)
                     / len(prot)) if prot else 0.0
        call_frac = (sum(x["step"] == "domain_call" for x in prot) / len(prot)) if prot else 0.0
        h = hsf.get(sf, [])
        sf_rows.append({"superfamily": sf, **c,
                        "proteome_rows": len(prot),
                        "domain_enumeration_frac": round(enum_frac, 4),
                        "domain_call_frac": round(call_frac, 4),
                        "q3_enumeration": q3(enum_frac), "q3_call": q3(call_frac),
                        "domain_only_calls": domain_only.get(sf, 0),
                        "human_genes": len(h),
                        "human_enumerated": sum(r["enumerated"] for r in h),
                        "human_domain_call": sum(r["domain_call"] for r in h),
                        "human_profile_call": sum(r["profile_call"] for r in h)})
    write_tsv(OUT / "curve_by_superfamily.tsv", list(sf_rows[0]), sf_rows)

    fam_rows = []
    by_fam = defaultdict(list)
    for x in frame:
        by_fam[x["family"]].append(x)
    for fam, rows in sorted(by_fam.items()):
        fam_rows.append({"family": fam, "superfamily": CATALOGUE[fam].superfamily,
                         **curve(rows)})
    write_tsv(OUT / "curve_by_family.tsv", list(fam_rows[0]), fam_rows)

    grp = defaultdict(list)
    for x in frame:
        if x["source"] != "genome":
            grp[(x["superfamily"], x["group"])].append(x)
    g_rows = [{"superfamily": sf, "group": g, "proteome_rows": len(v),
               "domain_enumeration_frac": round(sum(x["step"] != "profile" for x in v) / len(v), 4),
               "domain_call_frac": round(sum(x["step"] == "domain_call" for x in v) / len(v), 4)}
              for (sf, g), v in sorted(grp.items())]
    write_tsv(OUT / "curve_by_group.tsv", list(g_rows[0]), g_rows)

    tot = curve(frame)
    hum_tot = {k: sum(r[k] for r in hum) for k in ("enumerated", "domain_call", "profile_call")}
    print(f"frame {tot}")
    print(f"human {len(hum)} {hum_tot}")
    for r in sf_rows:
        print(r["superfamily"], r["n"], r["domain_call_frac"], r["domain_enumeration_frac"],
              r["cum_genome"], r["q3_enumeration"], r["q3_call"],
              f"human {r['human_enumerated']}/{r['human_domain_call']}/{r['human_profile_call']} of {r['human_genes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
