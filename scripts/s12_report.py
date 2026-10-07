"""s12_report.py — S12's report, rendered from its committed tables (D13).

    python3 scripts/s12_report.py   # → results/structures/report.md + summary.json

Reads `coverage_members.tsv`, `coverage_families.tsv`, `experimental.tsv`,
`network_edges.tsv`, `network_units.tsv`, `dense_reps.tsv`, `dense_edges.tsv`,
`dense_recovery*.tsv`, `dense_detect.tsv` and S8a's `literature_edges.tsv`.
A missing table renders as "not yet run".
"""

from __future__ import annotations

import json
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv  # noqa: E402

OUT = ROOT / "results" / "structures"


def _t(name: str) -> list[dict]:
    p = OUT / name
    return read_tsv(p) if p.exists() else []


def _table(rows: list[dict], cols: list[str]) -> list[str]:
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return out


def _pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f} %" if b else "–"


def coverage(s: dict) -> list[str]:
    mem, fam = _t("coverage_members.tsv"), _t("coverage_families.tsv")
    L = ["## 1. AlphaFold DB and PDB coverage of the census (D54 (1))", ""]
    if not mem:
        return L + ["Not yet run.", ""]
    n = len(mem)
    k = {c: sum(int(m[c]) for m in mem) for c in ("model", "exact", "usable")}
    pdb = sum(1 for m in mem if m["n_pdb"] != "0")
    loci = sum(int(f["genome_loci"]) for f in fam)
    pl = [float(m["global_plddt"]) for m in mem if m["global_plddt"]]
    s["coverage"] = {"members": n, "genome_loci": loci, **k, "with_pdb": pdb,
                     "median_plddt": round(statistics.median(pl), 1)}
    L += [f"Frame: S15's final census — **{n:,} proteome members** in "
          f"{len(fam)} census families (plus {loci} genome loci, which have no UniProt "
          "accession and so no AlphaFold DB entry).", "",
          f"* **AlphaFold DB model of the exact accession: {k['model']:,} "
          f"({_pct(k['model'], n)})**; model sequence identical to the census sequence "
          f"{k['exact']:,} ({_pct(k['exact'], n)}).",
          f"* **Usable** (exact, global pLDDT ≥ 70): **{k['usable']:,} "
          f"({_pct(k['usable'], n)})**; median global pLDDT {s['coverage']['median_plddt']}.",
          f"* **Any experimental structure** (UniProt PDB cross-reference): **{pdb} "
          f"({_pct(pdb, n)})**.", ""]
    zero = [f for f in fam if f["usable"] == "0"]
    low = sorted((f for f in fam if f["frac_usable"] and f["usable"] != "0"),
                 key=lambda f: float(f["frac_usable"]))[:8]
    nomodel = sorted((f for f in fam if f["frac_model"] and int(f["members"]) >= 20),
                     key=lambda f: float(f["frac_model"]))[:5]
    L += ["**Families with no usable model**: " + ", ".join(
          f"{f['family']} ({f['members']} members, {f['model']} models, median pLDDT "
          f"{f['median_plddt'] or '–'})" for f in zero) + ".", "",
          "**Lowest usable share** (families with ≥ 1 usable model):", ""]
    L += _table(low, ["family", "superfamily", "members", "model", "usable", "frac_usable",
                      "median_plddt"])
    long_ = Counter(m["family"] for m in mem if m["model"] == "0" and int(m["length"]) > 2700)
    nm = n - k["model"]
    s["coverage"]["no_model_over_2700aa"] = sum(long_.values())
    for f in nomodel:
        f["no_model_over_2700aa"] = long_.get(f["family"], 0)
    L += ["", f"**Why members lack a model.** Of the {nm:,} members with no AlphaFold DB "
          f"entry, only {sum(long_.values())} are longer than 2,700 residues (AlphaFold "
          "DB's length limit); the rest are simply not in the database. Lowest model "
          "share (≥ 20 members):", ""]
    L += _table(nomodel, ["family", "members", "model", "frac_model",
                          "no_model_over_2700aa"])
    L += ["", "Per family: `coverage_families.tsv`; per member: `coverage_members.tsv`.", ""]
    return L


def experimental(s: dict) -> list[str]:
    rows = [r for r in _t("experimental.tsv") if r["node"] != "kv_shaker:VSD"]
    L = ["## 2. The experimental reference and the model check (D54 (2)–(3))", ""]
    if not rows:
        return L + ["Not yet run.", ""]
    have = [r for r in rows if r["pdb"]]
    meth = Counter(r["method"] for r in have)
    chk = Counter(r.get("model_check", "") for r in rows)
    tm = [float(r["tm_by_exp"]) for r in rows if r.get("tm_by_exp")]
    s["experimental"] = {"nodes": len(rows), "with_reference": len(have),
                         "methods": dict(meth), "checks": dict(chk),
                         "median_tm": round(statistics.median(tm), 3) if tm else None}
    L += [f"**{len(have)} of {len(rows)} census families have an experimental reference** "
          "among their catalogue exemplars (" + ", ".join(f"{v} {k}" for k, v in
                                                        meth.most_common()) + ").",
          "Chains are cut by SIFTS per-residue UniProt numbering; residues SIFTS maps to an "
          "isoform are renumbered to the canonical sequence by global alignment (addendum).",
          "", f"**Model check**: {sum(v for k, v in chk.items() if k.startswith('consistent'))}"
          f" of {len(tm)} comparable AFDB units are consistent with experiment "
          f"(TM-score ≥ 0.5, normalised by the experimental unit), "
          f"{chk.get('consistent_high', 0)} at ≥ 0.8; median {s['experimental']['median_tm']}."
          f" {chk.get('different_accession', 0)} families' reference is a different exemplar "
          f"from S8a's node, {chk.get('no_model', 0)} has no model (RyR).", ""]
    bad = [r for r in rows if r.get("model_check") == "inconsistent"]
    if bad:
        L += ["Inconsistent:", ""] + _table(bad, ["node", "pdb", "chain", "method",
                                                  "observed", "tm_by_exp", "rmsd"]) + [""]
    ryr = next((r for r in rows if r["node"] == "ryr"), None)
    if ryr and ryr["pdb"]:
        L += [f"**RyR** (no AFDB model) enters through {ryr['exemplar']} "
              f"{ryr['pdb']} chain {ryr['chain']} ({ryr['method']}, "
              f"{ryr['resolution']} Å), module residues {ryr['unit_start']}–"
              f"{ryr['unit_end']}.", ""]
    none = [r["node"] for r in rows if not r["pdb"]]
    L += [f"No experimental reference ({len(none)}): " + ", ".join(none) + ".", ""]
    return L


def network(s: dict) -> list[str]:
    edges, dense = _t("network_edges.tsv"), _t("dense_edges.tsv")
    L = ["## 3. The literature edges read four ways (D54 (4)–(6))", ""]
    if not edges:
        return L + ["Not yet run.", ""]
    L += ["S8a's AFDB reading is the primary verdict (D48); the other three are "
          "sensitivity readings. Each cell: median TM-score, rank b-for-a / a-for-b, "
          "verdict.", ""]
    look = {(e["variant"], e["a"], e["b"]): e for e in edges}
    look.update({("dense", e["a"], e["b"]): e for e in dense})
    variants = ["s8a_afdb", "experimental", "tm_region"] + (["dense"] if dense else [])
    keys = [(e["a"], e["b"]) for e in edges if e["variant"] == "s8a_afdb"]
    out, s["edges"] = [], {}
    for a, b in keys:
        row = {"edge": f"{a} – {b}"}
        for v in variants:
            e = look.get((v, a, b))
            row[v] = ("unmeasured" if not e or not e["median_tm"] else
                      f"{e['median_tm']} ({e['rank_b_for_a']}/{e['rank_a_for_b']}) "
                      f"**{e['verdict'].replace('_', ' ')}**")
            s["edges"].setdefault(f"{a}-{b}", {})[v] = e["verdict"] if e else "unmeasured"
        out.append(row)
    L += _table(out, ["edge"] + variants) + [""]
    units = Counter(r["source"].split(":")[0] for r in _t("network_units.tsv")
                    if r["variant"] == "experimental")
    L += [f"Experimental variant: {units.get('pdb', 0)} nodes on an experimental unit, "
          f"{units.get('afdb', 0)} on S8a's AFDB unit.", ""]
    return L


def dense_set(s: dict) -> list[str]:
    reps, rec = _t("dense_reps.tsv"), _t("dense_recovery_by_superfamily.tsv")
    L = ["## 4. The dense set: recovery and detectability (D54 (6))", ""]
    if not rec:
        return L + ["Not yet run.", ""]
    ok = [r for r in reps if r["residues"] and int(r["residues"]) >= 30]
    fams = {r["family"] for r in ok}
    tr = sum(int(r["tm_recovered"]) for r in rec)
    fr = sum(int(r["foldseek_recovered"]) for r in rec)
    n = sum(int(r["reps"]) for r in rec)
    s["dense"] = {"reps": len(ok), "families": len(fams), "recovery_reps": n,
                  "tm_recovered": tr, "foldseek_recovered": fr}
    L += [f"**{len(ok)} representatives** of {len(fams)} census families (one per family × "
          f"S4 group; {len(reps) - len(ok)} dropped, reasons in `dense_reps.tsv`), "
          "all-vs-all TM-align and Foldseek.", "",
          f"**Superfamily recovery** (superfamilies with ≥ 2 measured families, {n} "
          f"representatives): the closest structure outside a representative's own family "
          f"is in its own superfamily for **{tr} ({_pct(tr, n)}) by TM-align** and "
          f"**{fr} ({_pct(fr, n)}) by Foldseek E-value**.", ""]
    L += _table(rec, ["superfamily", "families", "reps", "tm_recovered",
                      "foldseek_recovered"]) + [""]
    miss = [r for r in _t("dense_recovery.tsv") if r["tm_recovered"] == "0"]
    if miss:
        c = Counter((r["superfamily"], r["best_tm_partner"].split("__")[0]) for r in miss)
        L += ["Not recovered by TM-align — the family the closest partner belongs to:", ""]
        L += _table([{"superfamily": k[0], "closest_family": k[1], "reps": v}
                     for k, v in c.most_common(12)], ["superfamily", "closest_family", "reps"])
        L += [""]
    alld = [d for d in _t("dense_detect.tsv") if d["a"] != d["b"]]
    det = sorted((d for d in alld if d["foldseek_detected"] != "0"),
                 key=lambda d: -float(d["frac_detected"]))
    s["dense"]["sf_pairs_detected"] = f"{len(det)}/{len(alld)}"
    L += [f"**Foldseek detectability** between superfamilies (share of representative "
          f"pairs with E ≤ 1e-3): **{len(det)} of {len(alld)} superfamily pairs** have any "
          "detected pair:", ""]
    L += _table(det, ["a", "b", "rep_pairs", "foldseek_detected", "frac_detected",
                      "median_tm"]) + [""]
    return L


def main() -> int:
    s: dict = {}
    L = ["# S12 — Structures", "",
         "AlphaFold DB coverage of the final census, experimental references and a "
         "model check, and S8a's fold network re-read with experimental units, "
         "TM-region units and a dense set of representatives. Rules: **D54**, fixed "
         "before any structure was fetched; S8a's verdicts (D48) stay primary.", ""]
    for part in (coverage, experimental, network, dense_set):
        L += part(s)
    L += ["Figure: `results/structures/figures/structures.png`.", ""]
    (OUT / "report.md").write_text("\n".join(L))
    (OUT / "summary.json").write_text(json.dumps(s, indent=2, ensure_ascii=False))
    print(f"wrote {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
