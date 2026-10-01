"""s7c_report.py — the S7c section of results/phylogeny/tier1_report.md (D13).

Imported by `s7_report.py`; renders from `reroot_basal.tsv`,
`tier1_reroot_inputs.tsv`, `tier1_reroot_trees.tsv` and `tier1_reroot.tsv`.
"""

from __future__ import annotations

from collections import Counter

from s3_hmm_lib import read_tsv
from s7_lib import OUT_DIR
from src.catalogue import CATALOGUE


def _members(rs) -> str:
    return ", ".join(list(rs.families) + [e.label for e in rs.exemplars])


def reroot_section(L: list[str]) -> dict:
    ip = OUT_DIR / "tier1_reroot_inputs.tsv"
    if not ip.exists():
        return {}
    inputs = read_tsv(ip)
    basal = read_tsv(OUT_DIR / "reroot_basal.tsv")
    fams = sorted({r["family"] for r in inputs})
    L += ["## 4. Re-rooting the weakly rooted P-loop families (S7c, D47)", "",
          "S7b's KcsA/MthK/NaK outgroup was not one clade in four P-loop families and "
          "weakly supported in two more. Each of the six now declares **two outgroups "
          "in the catalogue** (`ChannelFamily.root_with`), fixed before any re-rooted "
          "tree: a prokaryotic relative with the family's architecture (inline "
          "exemplars verified live by `s0_outgroups.py` → `s0_baseline/outgroups.tsv`) "
          "and two sister eukaryotic families. K2P has no prokaryotic relative of its "
          "own, so its prokaryotic choice remains the superfamily outgroup.", "",
          "| family | outgroup 1 | outgroup 2 |", "|---|---|---|"]
    for f in fams:
        a, b = CATALOGUE[f].root_with
        L.append(f"| {f} | {a.name}: {_members(a)} | {b.name}: {_members(b)} |")
    added = Counter(r["family"] for r in basal if r["verdict"] == "add")
    dup = Counter(r["family"] for r in basal if r["verdict"] != "add")
    L += ["", "**Denser basal sampling (D47).** Each family's D39 set gains one member per "
          "early-diverging clade of S4b's dense panel (kingdom × phylum × class; outside "
          "Bilateria), the clade's highest-scoring high-confidence profile call, complete "
          "and inside the length band ±25 % (`reroot_basal.tsv`). The six sets therefore "
          "differ from S6's D39 sets — a stated exception. Re-aligned with L-INS-i; "
          "everything after is S7's method unchanged, and the ingroup alignment and "
          "trim mask are identical in a family's two trees.", "",
          "| family | D39 set | basal added | clades already in D39 | trimmed cols "
          "| outgroup occupancy (1 / 2) |", "|---|---|---|---|---|---|"]
    for f in fams:
        rs = [r for r in inputs if r["family"] == f]
        L.append(f"| {f} | {rs[0]['n_d39']} | {added[f]} | {dup[f]} | {rs[0]['cols']} | "
                 + " / ".join(r["og_occupancy"] for r in rs) + " |")
    L += ["", "**The root criterion, fixed before any re-rooted tree:** in both trees "
          "(a) the outgroup is one clade and (b) its edge has UFBoot ≥ 95, and (c) the "
          "ingroup's basal split is *identical* under the two outgroups. Anything else is "
          "`unresolved`; no root is picked by how it looks.", ""]
    tp, fp = OUT_DIR / "tier1_reroot_trees.tsv", OUT_DIR / "tier1_reroot.tsv"
    if not fp.exists():
        L += ["*Trees not yet built* (detached IQ-TREE run; `s7c_reroot.py parse` "
              "fills this section).", ""]
        return {"families": len(fams), "basal_added": sum(added.values())}
    trees, res = read_tsv(tp), read_tsv(fp)
    L += ["| family | tree | seqs | model | outgroup one clade | root UFBoot | "
          "basal split (sizes) | split UFBoot | ≥ 95 |", "|---|---|---|---|---|---|---|---|---|"]
    for t in trees:
        L.append(f"| {t['family']} | {t['root_set']} | {t['n_ingroup']} | {t['model']} | "
                 f"{t['outgroup_monophyletic']} | {t['root_ufboot'] or '—'} | "
                 f"{t['root_clade_sizes'] or '—'} | {t['root_clade_ufboot'] or '—'} | "
                 f"{t['frac_ge95']} |")
    L += ["", "| family | (a) one clade | (b) UFBoot ≥ 95 | (c) same root | smaller-clade "
          "Jaccard | verdict | S7b (KcsA): one clade / root UFBoot | S7b split vs new |",
          "|---|---|---|---|---|---|---|---|"]
    for r in res:
        L.append(f"| {r['family']} | {r['a_one_clade']} | {r['b_root_ufboot_ge95']} | "
                 f"{r['c_same_root']} | {r['smaller_clade_jaccard']} | "
                 f"**{r['verdict']}**{' — ' + r['reason'] if r['reason'] else ''} | "
                 f"{r['s7b_one_clade']} / {r['s7b_root_ufboot'] or '—'} | "
                 f"{r['s7b_split_on_shared'] or '—'} |")
    ok = [r["family"] for r in res if r["verdict"] == "resolved"]
    bad = [r["family"] for r in res if r["verdict"] != "resolved"]
    L += ["", f"**Resolved: {len(ok)}/{len(res)}** ({', '.join(ok) or 'none'}). "
          f"**Unresolved: {', '.join(bad) or 'none'}** — S11 must read these families' "
          "duplication order without a root, or from reconciliation with the species "
          "tree (its own job), and say so.", ""]
    return {"families": len(res), "resolved": ok, "unresolved": bad,
            "basal_added": sum(added.values())}
