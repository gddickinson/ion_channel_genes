"""s11b_repeats.py — S11b (1): the order of the internal repeat duplications (D53).

    python3 scripts/s11b_repeats.py spans    # repeat spans in profile states
    python3 scripts/s11b_repeats.py reps     # one chain per family × group
    python3 scripts/s11b_repeats.py prep     # cut, L-INS-i, trimAl -gt 0.5
    python3 scripts/s11b_repeats.py run [--threads 4]   # ML + 3 constrained + AU
    python3 scripts/s11b_repeats.py parse    # → repeat_au.tsv, repeat_classes.tsv

The repeat unit is the whole S1–S6 repeat (TRANSMEM helices 6r−5 … 6r of the
S6 references, projected into the family's S3a profile states); the
sensitivity unit is S6's pore module of the same chains. Each repeat is cut
from S6's `hmmalign` of the family members (`<fam>.members.a2m`). The test is
an AU test of the three pairings of the four 4×6TM repeat classes (H13
{I,III}|{II,IV}, H12, H14), each an ML search constrained by its bipartition;
TPC's two repeats are free. Rules fixed in D53 before any input was built.
"""

from __future__ import annotations

import argparse
import re
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import Fetcher  # noqa: E402
from s11b_lib import (AU_ALPHA, CLASSES, HYPOTHESES, MIN_REPEAT_COVER, OUT,  # noqa: E402
                      REPEAT_FAMILIES, iqtree, live, s11b_dir)
from s3_hmm_lib import iter_fasta, read_tsv, sha256, write_fasta, write_tsv  # noqa: E402
from s6_lib import s6_dir  # noqa: E402
from s6_module_refs import features, topology  # noqa: E402
from s6_modules import ref_sequence  # noqa: E402
from s6_project import a2m_map, cut, hmmalign, to_states  # noqa: E402
from s7_lib import apply_mask, metrics, read_aln, trim_columns, write_aln  # noqa: E402
from s7_newick import outgroup_intruders, parse, split_support, splits  # noqa: E402
from s7_trees import IQ_ARGS, TRIM  # noqa: E402
from s8_lib import base_label, centrality, family_alignment_rows, matrix, members  # noqa: E402

ALN = ROOT / "results" / "alignments"
UNITS = ("repeat", "pore")


def tip_class(fam: str, r: int) -> str:
    """Repeat r of a family → its class tag: I–IV for 4×6TM, tI/tII for TPC."""
    return ("t" if REPEAT_FAMILIES[fam] == 2 else "") + CLASSES[r - 1]


# ------------------------------------------------------------------ spans

def cmd_spans(_a) -> None:
    fetch, rows = Fetcher(), []
    refs = [r for r in read_tsv(ALN / "module_refs.tsv")
            if r["family"] in REPEAT_FAMILIES and r["status"] == "ok"]
    for fam, n in REPEAT_FAMILIES.items():
        seqs, want = [], []
        for r in (x for x in refs if x["family"] == fam):
            tm, _ = features(topology(fetch, r["reference"], False) or {})
            if len(tm) != 6 * n:
                print(f"  {fam} {r['reference']}: {len(tm)} TRANSMEM, not {6 * n} — skipped")
                continue
            lab = f"REF_{r['reference']}"
            seqs.append((lab, ref_sequence(fetch, r["reference"])))
            want.append((lab, [(tm[6 * i][0], tm[6 * i + 5][1]) for i in range(n)]))
        if not want:
            raise SystemExit(f"{fam}: no reference with {6 * n} TRANSMEM helices (D28)")
        maps = hmmalign(fam, seqs, "repeat_refs")
        votes = [[to_states(maps[lab][1], a, b) for a, b in sp] for lab, sp in want]
        votes = [v for v in votes if all(v)]
        for i in range(n):
            a = [v[i][0] for v in votes]
            b = [v[i][1] for v in votes]
            rows.append({"family": fam, "repeat": i + 1, "class": tip_class(fam, i + 1),
                         "k_start": int(statistics.median(a)),
                         "k_end": int(statistics.median(b)), "n_refs": len(votes),
                         "references": ",".join(lab[4:] for lab, _ in want),
                         "spread": max(max(a) - min(a), max(b) - min(b))})
            print(f"{fam:6s} R{i + 1} states {rows[-1]['k_start']}-{rows[-1]['k_end']} "
                  f"({len(votes)} refs, spread {rows[-1]['spread']})")
    write_tsv(OUT / "repeat_spans.tsv", list(rows[0]), rows)


def load_spans() -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list] = defaultdict(list)
    for r in read_tsv(OUT / "repeat_spans.tsv"):
        out[r["family"]].append((int(r["k_start"]), int(r["k_end"])))
    return out


# ------------------------------------------------------------------- reps

def member_cuts(fam: str, spans: list[tuple[int, int]]) -> dict[str, list[tuple]]:
    """label → [(start, end, state cover)] per repeat, from S6's members.a2m."""
    out = {}
    for h, s in iter_fasta(s6_dir("project") / f"{fam}.members.a2m"):
        st, res = a2m_map(s)
        out[h.split()[0]] = [cut(st, res, "", k0, k1) for k0, k1 in spans]
    return out


def cmd_reps(_a) -> None:
    spans, mem, rows = load_spans(), members(), []
    for fam in REPEAT_FAMILIES:
        cuts = member_cuts(fam, spans[fam])
        aln = family_alignment_rows(fam)
        cell: dict[str, list[str]] = defaultdict(list)
        for (f, lab), r in mem.items():
            if f == fam:
                cell[r["group"]].append(lab)
        for g, labs in sorted(cell.items()):
            ok = sorted(l for l in labs if l in cuts
                        and all(c[2] >= MIN_REPEAT_COVER for c in cuts[l]))
            if not ok:
                rows.append({"family": fam, "group": g, "label": "", "n_cell": len(labs),
                             "n_eligible": 0, "centrality": "", "covers": ""})
                continue
            cen = centrality(ok, matrix({l: aln[l] for l in ok}))
            pick = min(ok, key=lambda l: (-cen[l], l))
            rows.append({"family": fam, "group": g, "label": pick, "n_cell": len(labs),
                         "n_eligible": len(ok), "centrality": round(cen[pick], 3),
                         "covers": ",".join(f"{c[2]:.2f}" for c in cuts[pick])})
            print(f"{fam:6s} {g:15s} {len(ok):>3}/{len(labs):<3} → {pick}")
    write_tsv(OUT / "repeat_reps.tsv", list(rows[0]), rows)


# ------------------------------------------------------------------- prep

def segments(unit: str) -> list[tuple[str, str]]:
    spans = load_spans()
    seqs = {fam: {h.split()[0]: s for h, s in iter_fasta(s6_dir("family") / f"{fam}.fasta")}
            for fam in REPEAT_FAMILIES}
    pore = {(r["family"], base_label(r["label"]), int(r["module"])): (int(r["start"]), int(r["end"]))
            for r in read_tsv(ALN / "modules.tsv")
            if r["family"] in REPEAT_FAMILIES and r["quality"] == "full"}
    out = []
    for r in read_tsv(OUT / "repeat_reps.tsv"):
        fam, lab = r["family"], r["label"]
        if not lab:
            continue
        cuts = member_cuts(fam, spans[fam])[lab]
        for i, (a, b, _) in enumerate(cuts, 1):
            if unit == "pore":
                if (fam, lab, i) not in pore:
                    raise SystemExit(f"{fam} {lab} module {i}: no full S6 pore module")
                a, b = pore[(fam, lab, i)]
            out.append((f"{tip_class(fam, i)}__{fam}__{lab}", seqs[fam][lab][a - 1:b]))
    return out


def cmd_prep(_a) -> None:
    rows = []
    for unit in UNITS:
        d = s11b_dir("repeats", unit)
        raw = d / f"{unit}.fasta"
        write_fasta(raw, segments(unit))
        aln = d / f"{unit}.linsi.fasta"
        p = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "1",
                            "--quiet", str(raw)], capture_output=True, text=True)
        if p.returncode != 0 or not p.stdout.strip():
            raise RuntimeError(f"mafft {unit}: {p.stderr.strip()[:300]}")
        aln.write_text(p.stdout)
        rows_in = read_aln(aln)
        if len(rows_in) != len(list(iter_fasta(raw))):
            raise RuntimeError(f"{unit}: row count changed in alignment")
        cols = trim_columns(aln, TRIM)
        trimmed = apply_mask([(l, s.upper()) for l, s in rows_in], cols)
        out = d / f"{unit}.input.fasta"
        write_aln(out, trimmed)
        m = metrics(trimmed, list(range(len(cols))))
        rows.append({"unit": unit, "tips": len(trimmed), "aligned_cols": len(rows_in[0][1]),
                     "cols": m["cols"], "informative": m["informative"],
                     "gap_frac": m["gap_frac"], "input_sha256": sha256(out)})
        print(rows[-1])
    write_tsv(OUT / "repeat_inputs.tsv", list(rows[0]), rows)


# -------------------------------------------------------------------- run

def constraint(tips: list[str], hyp: str, free: str = "") -> str:
    """Newick of the hypothesis's bipartition over the 4×6TM tips only, less
    the `free` tip (see `free_tip`)."""
    side = lambda cls: ",".join(t for t in tips  # noqa: E731
                                if t.split("__")[0] in cls and t != free)
    (a, b) = HYPOTHESES[hyp]
    return f"(({side(a)}),({side(b)}));\n"


def pruned(node, keep: set[str]) -> str | None:
    """Topology of `node` restricted to `keep`, as Newick (no lengths)."""
    if not node.children:
        return node.name if node.name in keep else None
    kids = [k for k in (pruned(c, keep) for c in node.children) if k]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else f"({','.join(kids)})"


def free_tip(tips: list[str]) -> str:
    """The one 4×6TM tip left out of the constraint. IQ-TREE 2.3.6 rejects a
    constraint over exactly 64 taxa ("Initial tree is not compatible"; 63 and
    65 pass — measured, and it refuses an all-gap padding row), which 16
    chains × 4 repeats hit. The tree found is checked afterwards to hold the
    full bipartition (`satisfies`); if it does it is also the optimum under the
    full constraint, if not the hypothesis is flagged, never silently used."""
    return sorted(tips)[-1] if len(tips) == 64 else ""


def retry_tip(tips: list[str]) -> str:
    """The tip freed on the one retry when the first constrained tree missed
    its bipartition (added after H14's whole-repeat tree did, D53 note)."""
    return sorted(tips)[0] if len(tips) == 64 else ""


def constrained_tree(d: Path, hyp: str) -> Path:
    """The constrained tree used for a hypothesis: the retry's, if one ran."""
    b = d / f"{hyp}_retry.treefile"
    return b if b.exists() else d / f"{hyp}.treefile"


def satisfies(tree_text: str, hyp: str) -> bool:
    """Does the tree hold the hypothesis's bipartition over the 4×6TM tips?
    TPC tips are free, so they may sit on either side (or inside one): the
    test is on each edge's split restricted to the 4×6TM tips."""
    tree = parse(tree_text)
    fourx = frozenset(l for l in tree.leaves() if not l.startswith("t"))
    side = frozenset(l for l in fourx if l.split("__")[0] in HYPOTHESES[hyp][0])
    if not side or side == fourx:
        return False
    return any(s & fourx in (side, fourx - side) for s, _ in splits(tree))


def best_model(iq: Path) -> str:
    m = re.search(r"Best-fit model according to BIC:\s*(\S+)", iq.read_text())
    if not m:
        raise RuntimeError(f"no BIC model in {iq}")
    return m.group(1)


def cmd_run(a) -> None:
    steps = [{"label": f"{u} {s}", "done": False} for u in UNITS
             for s in ("ML", *HYPOTHESES, "AU")]
    for unit in UNITS:
        d = s11b_dir("repeats", unit)
        src = d / f"{unit}.input.fasta"
        tips = [l for l, _ in read_aln(src)]
        iqtree(src, d / "ml", IQ_ARGS, a.threads)
        model = best_model(d / "ml.iqtree")
        names = [l for l, _ in read_aln(src)]
        trees = [pruned(parse((d / "ml.treefile").read_text()), set(names)) + ";"]
        fourx = [t for t in tips if not t.startswith("t")]
        free = free_tip(fourx)
        for st in steps:
            st["done"] |= st["label"] == f"{unit} ML"
        live(steps)
        for h in HYPOTHESES:
            c = d / f"{h}.constraint.nwk"
            c.write_text(constraint(fourx, h, free))
            iqtree(src, d / h, ["-m", model, "-g", str(c), "-seed", "1"], a.threads)
            if free and not satisfies((d / f"{h}.treefile").read_text(), h):
                c2 = d / f"{h}_retry.constraint.nwk"
                c2.write_text(constraint(fourx, h, retry_tip(fourx)))
                iqtree(src, d / f"{h}_retry", ["-m", model, "-g", str(c2), "-seed", "1"],
                       a.threads)
            trees.append(pruned(parse(constrained_tree(d, h).read_text()), set(names)) + ";")
            for st in steps:
                st["done"] |= st["label"] == f"{unit} {h}"
            live(steps)
        cand = d / "candidates.trees"
        cand.write_text("\n".join(trees) + "\n")
        iqtree(src, d / "au", ["-m", model, "-z", str(cand), "-zb", "10000", "-au",
                               "-n", "0", "-seed", "1"], a.threads, deps=(cand,))
        for st in steps:
            st["done"] |= st["label"] == f"{unit} AU"
        live(steps)
        print(f"{unit}: done ({model})", flush=True)


# ------------------------------------------------------------------ parse

def au_table(iq: Path) -> list[dict]:
    text = iq.read_text()
    block = text[text.index("USER TREES"):]
    out = []
    for line in block.splitlines():
        tok = [t for t in line.split() if t not in "+-"]
        if len(tok) >= 8 and tok[0].isdigit():
            out.append({"tree": int(tok[0]), "logL": float(tok[1]), "deltaL": float(tok[2]),
                        "p_SH": float(tok[5]), "p_AU": float(tok[7])})
    return out


def class_rows(unit: str, tree) -> list[dict]:
    leaves = tree.leaves()
    rows = []
    for cls in sorted({l.split("__")[0] for l in leaves}, key=lambda c: (c[0] == "t", c)):
        grp = {l for l in leaves if l.split("__")[0] == cls}
        present, sup = split_support(tree, grp)
        intr = sorted(outgroup_intruders(tree, grp))
        # Smallest edge side holding the whole class plus any other-class tip.
        sides = [s for side, _ in splits(tree) for s in (side, frozenset(leaves) - side)]
        near = min((s for s in sides if grp <= s and len(s) > len(grp)),
                   key=len, default=frozenset())
        comp = defaultdict(int)
        for l in near - grp:
            comp[l.split("__")[0]] += 1
        rows.append({"unit": unit, "class": cls, "tips": len(grp), "one_clade": present,
                     "ufboot": "" if sup is None else sup, "intruders": len(intr),
                     "intruder_classes": ",".join(sorted({i.split("__")[0] for i in intr})),
                     "nearest_side": ";".join(f"{k}:{v}" for k, v in sorted(comp.items()))})
    return rows


def cmd_parse(_a) -> None:
    au, cls = [], []
    names = ["ML", *HYPOTHESES]
    for unit in UNITS:
        d = s11b_dir("repeats", unit)
        model = best_model(d / "ml.iqtree")
        for r in au_table(d / "au.iqtree"):
            name = names[r["tree"] - 1]
            au.append({"unit": unit, "tree": name, "model": model, **{k: r[k] for k in
                       ("logL", "deltaL", "p_SH", "p_AU")},
                       "satisfies": "" if name == "ML" else
                       satisfies(constrained_tree(d, name).read_text(), name),
                       "retry": "" if name == "ML" else
                       constrained_tree(d, name).name.endswith("_retry.treefile"),
                       "rejected": name != "ML" and r["p_AU"] < AU_ALPHA})
        tree = parse((d / "ml.treefile").read_text())
        cls.extend(class_rows(unit, tree))
        (OUT / "repeats").mkdir(exist_ok=True)
        (OUT / "repeats" / f"{unit}.treefile").write_text((d / "ml.treefile").read_text())
    write_tsv(OUT / "repeat_au.tsv", list(au[0]), au)
    write_tsv(OUT / "repeat_classes.tsv", list(cls[0]), cls)
    for r in au:
        print(f"{r['unit']:6s} {r['tree']:4s} ΔlnL {r['deltaL']:8.2f}  p-AU {r['p_AU']:.4f}"
              f"{'  rejected' if r['rejected'] else ''}")


def verdict(au: list[dict], unit: str) -> str:
    """D53 (4): one unrejected hypothesis names the pairing; more → unresolved."""
    rows = [r for r in au if r["unit"] == unit and r["tree"] != "ML"]
    bad = [r["tree"] for r in rows if str(r.get("satisfies", "True")) == "False"]
    if bad:                     # a constrained tree that misses its bipartition
        return "constraint_violated:" + ",".join(bad)
    alive = [r["tree"] for r in rows if str(r["rejected"]) == "False"]
    if len(alive) == 1:
        return alive[0]
    return "unresolved:" + ",".join(alive) if alive else "all_rejected"


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("spans", "reps", "prep", "parse"):
        sub.add_parser(name)
    r = sub.add_parser("run")
    r.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    {"spans": cmd_spans, "reps": cmd_reps, "prep": cmd_prep, "run": cmd_run,
     "parse": cmd_parse}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
