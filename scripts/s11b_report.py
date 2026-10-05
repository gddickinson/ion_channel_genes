"""s11b_report.py — S11b's report, rendered from its committed tables (D13).

    python3 scripts/s11b_report.py   # → results/duplication/s11b_report.md + s11b_summary.json

Reads `repeat_spans.tsv`, `repeat_reps.tsv`, `repeat_inputs.tsv`,
`repeat_au.tsv`, `repeat_classes.tsv`, `roots_inputs.tsv`, `roots_nonrev.tsv`
and S11a's `recon_trees.tsv`. A missing table renders as "not yet run".
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from s11b_lib import AU_ALPHA, CONTROL_FAMILIES, OUT, ROOT_FAMILIES, ROOTSTRAP_MIN  # noqa: E402
from s11b_repeats import UNITS, verdict  # noqa: E402
from s3_hmm_lib import read_tsv  # noqa: E402


def _t(name: str) -> list[dict]:
    p = OUT / name
    return read_tsv(p) if p.exists() else []


def _table(rows: list[dict], cols: list[str]) -> list[str]:
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(str(r.get(c, "")) for c in cols) + " |" for r in rows]
    return out


def repeats(summary: dict) -> list[str]:
    spans, reps, inputs = _t("repeat_spans.tsv"), _t("repeat_reps.tsv"), _t("repeat_inputs.tsv")
    au, cls = _t("repeat_au.tsv"), _t("repeat_classes.tsv")
    L = ["## 1. The order of the internal repeat duplications", ""]
    if not au:
        return L + ["Not yet run.", ""]
    picked = [r for r in reps if r["label"]]
    L += [f"**Rules (D53 (1)–(4)), fixed before any input was built.** Whole S1–S6 "
          f"repeats ({len(spans)} spans in profile states, from the S6 references' "
          "TRANSMEM helices) of one chain per family × S4 group — the most central chain "
          "whose every repeat covers ≥ 90 % of its profile states — aligned together "
          "(L-INS-i, trimAl `-gt 0.5`); the pore module of the same chains is the "
          "sensitivity unit. Three pairings of the four 4×6TM repeat classes are each an "
          "ML search constrained by their bipartition (TPC's two repeats free) and "
          f"compared with the unconstrained ML tree by the AU test (rejected iff p < {AU_ALPHA}).",
          "", f"**Representatives**: {len(picked)} chains "
          + ", ".join(f"{f} {sum(r['family'] == f for r in picked)}"
                      for f in dict.fromkeys(r["family"] for r in picked))
          + f"; cells with no eligible chain: {sum(not r['label'] for r in reps)}.", ""]
    L += _table(inputs, ["unit", "tips", "aligned_cols", "cols", "informative", "gap_frac"])
    L += ["", "**AU test** (ΔlnL from the best tree):", ""]
    L += _table(au, ["unit", "tree", "model", "logL", "deltaL", "p_SH", "p_AU", "satisfies", "rejected"])
    L += [""]
    for u in UNITS:
        v = verdict(au, u)
        summary[f"repeat_verdict_{u}"] = v
        L.append(f"* **{u}** unit: **{v}**")
    L += ["", "Pairings: **H13** {I,III}|{II,IV} (a tandem duplication of a two-repeat, "
          "TPC-like ancestor — stated in advance as the only pairing that fits it), "
          "**H12** {I,II}|{III,IV}, **H14** {I,IV}|{II,III}.", "",
          "**Repeat classes in the unconstrained ML tree** (descriptive): is each class "
          "one clade across Nav, Cav and NALCN, how many other tips stand in the way, and "
          "the class make-up of the smallest edge side that holds the whole class plus "
          "anything else.", ""]
    L += _table(cls, ["unit", "class", "tips", "one_clade", "ufboot", "intruders",
                      "intruder_classes", "nearest_side"])
    return L + [""]


def roots(summary: dict) -> list[str]:
    inp, res = _t("roots_inputs.tsv"), _t("roots_nonrev.tsv")
    recon = {r["family"]: r for r in _t("recon_trees.tsv")}
    L = ["## 2. Outgroup-free roots (non-reversible model)", "",
         "**Rules (D53 (5)–(7)).** S7b's ingroup alignment for each family (outgroup "
         "rows removed, D41 mask unchanged), NQ.pfam with S7b's rate-heterogeneity terms, "
         "1,000 UFBoot (→ rootstrap) and an AU test over every root branch. A root is "
         f"**resolved** iff rootstrap ≥ {ROOTSTRAP_MIN:g} and the non-reversible model fits "
         "better than Q.pfam on the same tree (ΔlnL > 0). It **agrees** with S11a iff its "
         "split is one of S11a's optimal reconciliation roots on S7b's tree. Controls: "
         "families whose declared outgroup root S11a accepted and also found optimal.", ""]
    L += _table(inp, ["family", "role", "n_tips", "cols", "s7b_model", "nq_model"])
    L += [""]
    if not res:
        return L + ["Runs not finished — results pending.", ""]
    done = {r["family"] for r in res}
    missing = [f for f in ROOT_FAMILIES + CONTROL_FAMILIES if f not in done]
    if missing:
        L += [f"**Not finished:** {', '.join(missing)}.", ""]
    for r in res:
        r["s11a_split"] = recon.get(r["family"], {}).get("root_split", "")
    ctl = [r for r in res if r["role"] == "control"]
    six = [r for r in res if r["role"] == "root"]
    L += ["**Controls** (declared root known):", ""]
    L += _table(ctl, ["family", "n_tips", "dlnl_nonrev", "root_rootstrap", "root_split",
                      "ml_root_tip", "au_set", "branches_tested", "au_set_frac",
                      "max_rootstrap", "control_recovered", "declared_in_au_set"])
    rec = sum(r["control_recovered"] == "True" for r in ctl)
    inau = sum(r.get("declared_in_au_set") == "True" for r in ctl)
    summary.update(controls=len(ctl), controls_recovered=rec, controls_declared_in_au=inau)
    L += ["", f"**{rec} / {len(ctl)} controls recover the declared root as the ML root; "
          f"the declared root is inside the AU confidence set in {inau} / {len(ctl)}.**", "",
          "**S7d's six:**", ""]
    L += _table(six, ["family", "n_tips", "dlnl_nonrev", "root_rootstrap", "root_split",
                      "small_side", "au_set", "au_set_frac", "max_rootstrap", "resolved",
                      "s11a_split",
                      "agrees_s11a", "recon_optimum_nq_topology", "s11a_root_in_au_set"])
    nres = sum(r["resolved"] == "True" for r in six)
    nag = sum(r["agrees_s11a"] == "True" for r in six)
    ins = sum(r["s11a_root_in_au_set"] == "True" for r in six)
    tip = sum(bool(r["ml_root_tip"]) for r in res)
    fr = sorted(float(r["au_set_frac"]) for r in res)
    summary.update(six_done=len(six), six_resolved=nres, six_agree_s11a=nag,
                   six_s11a_in_au=ins, ml_root_single_tip=tip,
                   au_set_frac_range=[fr[0], fr[-1]])
    L += ["", f"**{nres} / {len(six)} roots resolved under D53 (6); "
          f"{nag} / {len(six)} agree with S11a's reconciliation root**; S11a's root lies "
          f"inside the AU confidence set in {ins} / {len(six)}. Disagreements are reported, "
          "not resolved.", "",
          f"**What the roots look like** (descriptive): in {tip} of {len(res)} families the "
          "non-reversible ML root falls on a single terminal branch, and the AU confidence "
          f"set holds {fr[0]:.0%}–{fr[-1]:.0%} of all branches in every family — the "
          "non-reversible signal barely discriminates among roots at these depths. The "
          "controls, whose roots are known, measure it directly.", ""]
    return L


def main() -> int:
    summary: dict = {}
    L = ["# S11b — repeat-duplication order and outgroup-free roots", "",
         "Rendered by `scripts/s11b_report.py` from the tables beside it (D13). "
         "Rules: D53. Bulk under `<data root>/trees/s11b/`.", ""]
    L += repeats(summary)
    L += roots(summary)
    (OUT / "s11b_report.md").write_text("\n".join(L) + "\n")
    (OUT / "s11b_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
