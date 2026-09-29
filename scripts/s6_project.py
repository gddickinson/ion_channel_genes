"""s6_project.py — S6 step 2c: pore modules by projection through the family profile.

    python3 scripts/s6_project.py spans      # module span per family, in profile states
    python3 scripts/s6_project.py extract    # every D39 member → modules.tsv + FASTA
    python3 scripts/s6_project.py validate   # extracted vs UniProt-annotated modules

**D40 (as adopted).** Every member of a tier-2 unit is cut the same way:
`hmmalign` to its own family's S3a profile (the instrument that called it),
then the residues aligned to that family's module span of match states —
flanking inserts inside the span included. Method `profile_projection`, one
method for every sequence in the unit.

What varies is how each *family's* span was fixed, recorded per family in
`module_spans.tsv`:

* `annotated` — the UniProt module spans of the family's annotated
  references (`s6_module_refs.py`), mapped to match states; median over
  references, spread reported;
* `unit_hmm_vote` — for the families with no annotated reference, the
  unit's module HMM (`s6_modules.py seeds`) is searched over the members;
  members carrying exactly the expected number of full-length module hits
  vote, median per module. The vote is itself measured: it is also run on
  every annotated family and its error against the annotated span reported.

A single module HMM per unit was measured first and rejected as the
extraction instrument (`module_loo.tsv`: held-out Kir 22/298, TRPM 0/145,
innexin clan ≤ 5/82 recovered).
"""

from __future__ import annotations

import argparse
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import Fetcher  # noqa: E402
from s3_hmm_lib import iter_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402
from s6_lib import OUT_DIR, s3_profiles, s6_dir  # noqa: E402
from s6_module_refs import evaluate, features, module_families, topology  # noqa: E402
from s6_modules import (MIN_HMM_COVER, expected, k_filter_family,  # noqa: E402
                        ref_sequence, search)
from src.classify.motifs import K_FILTER_RE  # noqa: E402

MIN_SPAN_COVER = 0.5     # member occupies ≥ half the span's states → `full`
MAX_VALIDATE = 25        # reviewed members per family checked against UniProt
SNAP_TOL = 12            # the measured vote error bound (held-out, module_spans)


def a2m_map(a2m: str) -> tuple[list, list]:
    """A2M row → (state→residue index or None, residue→state or insert-after)."""
    st, res = [None], []          # 1-based states
    k = ri = 0
    for c in a2m:
        if c == ".":
            continue
        if c.isupper() or c == "-":
            k += 1
            if c == "-":
                st.append(None)
                continue
            ri += 1
            st.append(ri)
            res.append(k)
        else:
            ri += 1
            res.append(k + 0.5)   # an insert after state k
    return st, res


def hmmalign(fam: str, seqs: list[tuple[str, str]], tag: str) -> dict[str, tuple]:
    d = s6_dir("project")
    fa, out = d / f"{fam}.{tag}.fasta", d / f"{fam}.{tag}.a2m"
    write_fasta(fa, seqs)
    p = subprocess.run(["hmmalign", "--amino", "--outformat", "A2M", "-o", str(out),
                        str(s3_profiles() / f"{fam}.hmm"), str(fa)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"hmmalign {fam}: {p.stderr.strip()[:300]}")
    return {h.split()[0]: a2m_map(s) for h, s in iter_fasta(out)}


def to_states(res: list, a: int, b: int) -> tuple[int, int] | None:
    """Residue span [a, b] → (first, last) match state inside it."""
    ks = [res[i - 1] for i in range(a, b + 1) if float(res[i - 1]).is_integer()]
    return (int(ks[0]), int(ks[-1])) if ks else None


def cut(st: list, res: list, seq: str, k0: int, k1: int) -> tuple[int, int, float]:
    """Residues aligned in [k0, k1] (inserts inside included) + state cover."""
    idx = [i for i, k in enumerate(res, 1) if k0 <= k <= k1]
    occ = sum(1 for k in range(k0, k1 + 1) if st[k] is not None)
    if not idx:
        return 0, 0, 0.0
    return idx[0], idx[-1], occ / (k1 - k0 + 1)


def family_members(fam: str) -> list[tuple[str, str]]:
    return [(h.split()[0], s) for h, s in iter_fasta(s6_dir("family") / f"{fam}.fasta")]


def _median_spans(votes: list[list[tuple[int, int]]], n: int) -> list[tuple]:
    out = []
    for m in range(n):
        a = [v[m][0] for v in votes]
        b = [v[m][1] for v in votes]
        out.append((int(statistics.median(a)), int(statistics.median(b)),
                    max(max(a) - min(a), max(b) - min(b))))
    return out


def annotated_spans(fam: str, refs: list[dict]) -> list[tuple] | None:
    fetch, seqs, want = Fetcher(), [], []
    for r in refs:
        seqs.append((f"REF_{r['reference']}", ref_sequence(fetch, r["reference"])))
        want.append((f"REF_{r['reference']}",
                     [tuple(map(int, s.split("-"))) for s in r["spans"].split(";")]))
    maps = hmmalign(fam, seqs, "refs")
    votes = []
    for lab, spans in want:
        if len(spans) != expected(fam):
            continue
        ks = [to_states(maps[lab][1], a, b) for a, b in spans]
        if all(ks):
            votes.append(ks)
    return _median_spans(votes, expected(fam)) + [len(votes)] if votes else None


def vote_spans(unit: str, fam: str, held_out: bool) -> list[tuple] | None:
    """The vote. For an annotated family it runs on the model built *without*
    that family's seeds (`s6_modules.py loo`), so its error is measured the
    way the unannotated families experience it."""
    members = family_members(fam)
    hmm = (s6_dir("modules", "loo") / f"{unit}_minus_{fam}.module.hmm"
           if held_out else s6_dir("modules") / f"{unit}.module.hmm")
    if not hmm.exists():
        return None
    hits = search(hmm,
                  s6_dir("family") / f"{fam}.fasta",
                  s6_dir("project") / f"{fam}.vote.domtbl")
    maps = hmmalign(fam, members, "members")
    votes = []
    for lab, _ in members:
        full = [r for r in hits.get(lab, [])
                if (r["hmm_to"] - r["hmm_from"] + 1) / r["qlen"] >= MIN_HMM_COVER]
        if len(full) != expected(fam):
            continue
        ks = [to_states(maps[lab][1], r["env_from"], r["env_to"]) for r in full]
        if all(ks):
            votes.append(ks)
    return _median_spans(votes, expected(fam)) + [len(votes)] if votes else None


def helix_snap(fam: str, refs: list[dict], span: tuple) -> tuple[tuple, str]:
    """A module boundary lies on a helix. A vote boundary more than SNAP_TOL
    residues outside every predicted helix of the family's first reference
    with a predicted topology moves outward to the bracketing helix."""
    fetch = Fetcher()
    for r in refs:
        d = topology(fetch, r["reference"], False) if r.get("reference") else None
        tm = features(d)[0] if d else []
        if tm:
            break
    else:
        return span, "no_predicted_topology"
    acc = r["reference"]
    lab = f"REF_{acc}"
    st, res = hmmalign(fam, [(lab, ref_sequence(fetch, acc))], "snap")[lab]
    a, b, _ = cut(st, res, "", span[0], span[1])
    if not a:
        return span, f"{acc}:not_covered"
    na, nb = a, b
    if not any(t0 - SNAP_TOL <= a <= t1 for t0, t1 in tm):
        prev = [t for t in tm if t[1] < a]
        na = prev[-1][0] if prev else a
    if not any(t0 <= b <= t1 + SNAP_TOL for t0, t1 in tm):
        nxt = [t for t in tm if t[0] > b]
        nb = nxt[0][1] if nxt else b
    if (na, nb) == (a, b):
        return span, f"{acc}:{a}-{b} on helices"
    ks = to_states(res, na, nb)
    return ((ks[0], ks[1]) + tuple(span[2:]) if ks else span,
            f"{acc}:{a}-{b} snapped to {na}-{nb}")


SPAN_FIELDS = ["unit", "family", "basis", "module", "k_start", "k_end", "spread",
               "n_votes", "vote_k_start", "vote_k_end", "vote_n", "vote_error",
               "helix_check"]


def cmd_spans(_a) -> None:
    refs, anyref = defaultdict(list), defaultdict(list)
    for r in read_tsv(OUT_DIR / "module_refs.tsv"):
        if r["status"] == "ok":
            refs[r["family"]].append(r)
        if r["reference"]:
            anyref[r["family"]].append(r)
    rows = []
    for fam, unit, _ in module_families():
        ann = annotated_spans(fam, refs[fam]) if refs[fam] else None
        vote = vote_spans(unit, fam, held_out=bool(ann))
        chosen, basis = (ann, "annotated") if ann else (vote, "unit_hmm_vote")
        if not chosen:
            rows.append({"unit": unit, "family": fam, "basis": "NONE"})
            print(f"{fam}: no span (D28)")
            continue
        checks = [""] * expected(fam)
        if basis == "unit_hmm_vote":
            snapped = [helix_snap(fam, anyref[fam], chosen[m])
                       for m in range(expected(fam))]
            chosen = [sp for sp, _ in snapped] + [chosen[-1]]
            checks = [c for _, c in snapped]
            if any("snapped" in c for c in checks):
                basis = "unit_hmm_vote+helix_snap"
        for m in range(expected(fam)):
            row = {"unit": unit, "family": fam, "basis": basis, "module": m + 1,
                   "helix_check": checks[m],
                   "k_start": chosen[m][0], "k_end": chosen[m][1],
                   "spread": chosen[m][2], "n_votes": chosen[-1]}
            if vote:
                row.update(vote_k_start=vote[m][0], vote_k_end=vote[m][1],
                           vote_n=vote[-1])
                if ann:
                    row["vote_error"] = max(abs(vote[m][0] - ann[m][0]),
                                            abs(vote[m][1] - ann[m][1]))
            rows.append(row)
        print(f"{fam}: {basis} {[c[:2] for c in chosen[:-1]]}", flush=True)
    write_tsv(OUT_DIR / "module_spans.tsv", SPAN_FIELDS, rows)


def load_spans() -> dict[str, list[tuple[int, int, str]]]:
    out = defaultdict(list)
    for r in read_tsv(OUT_DIR / "module_spans.tsv"):
        if r["basis"] != "NONE":
            out[r["family"]].append((int(r["k_start"]), int(r["k_end"]), r["basis"]))
    return out


MANIFEST_FIELDS = ["unit", "family", "label", "module", "method", "span_basis",
                   "k_start", "k_end", "start", "end", "length", "span_cover",
                   "quality", "k_filter", "k_filter_in_chain"]


def cmd_extract(_a) -> None:
    spans, manifest = load_spans(), []
    by_unit = defaultdict(list)
    for fam, unit, _ in module_families():
        if fam not in spans:
            manifest.append({"unit": unit, "family": fam, "quality": "no_span"})
            continue
        members = family_members(fam)
        maps = hmmalign(fam, members, "members")
        n = len(spans[fam])
        for lab, seq in members:
            st, res = maps[lab]
            for m, (k0, k1, basis) in enumerate(spans[fam], 1):
                a, b, cov = cut(st, res, seq, k0, k1)
                mod = seq[a - 1:b] if a else ""
                q = "full" if cov >= MIN_SPAN_COVER else ("partial" if a else "absent")
                ml = f"{lab}_m{m}" if n > 1 else lab
                kf = ("yes" if K_FILTER_RE.search(mod) else "no") \
                    if k_filter_family(fam) and mod else ""
                kc = ("yes" if K_FILTER_RE.search(seq) else "no") if kf else ""
                manifest.append({"unit": unit, "family": fam, "label": ml,
                                 "module": m, "method": "profile_projection",
                                 "span_basis": basis, "k_start": k0, "k_end": k1,
                                 "start": a or "", "end": b or "",
                                 "length": len(mod), "span_cover": round(cov, 3),
                                 "quality": q, "k_filter": kf,
                                 "k_filter_in_chain": kc})
                if q == "full":
                    by_unit[unit].append((f"{fam}|{ml}", mod))
    d = s6_dir("modules")
    for unit, rows in sorted(by_unit.items()):
        write_fasta(d / f"{unit}.modules.fasta", rows)
        print(f"{unit}: {len(rows)} full modules → {unit}.modules.fasta")
    write_tsv(OUT_DIR / "modules.tsv", MANIFEST_FIELDS, manifest)


def cmd_validate(_a) -> None:
    """Extracted module vs the member's own UniProt module (not a reference)."""
    spans, fetch, rows = load_spans(), Fetcher(), []
    used = {r["reference"] for r in read_tsv(OUT_DIR / "module_refs.tsv")}
    for fam, unit, rule in module_families():
        if fam not in spans:
            continue
        mem = [r for r in read_tsv(OUT_DIR / "members.tsv")
               if r["family"] == fam and r["verdict"] == "include"
               and r["target"].startswith("sp|")
               and r["target"].split("|")[1] not in used][:MAX_VALIDATE]
        seqs = dict(family_members(fam))
        members = [(m["label"], seqs[m["label"]]) for m in mem if m["label"] in seqs]
        if not members:
            continue
        maps = hmmalign(fam, members, "validate")
        for m in mem:
            acc = m["target"].split("|")[1]
            d = topology(fetch, acc, False)
            truth, _, _, st = evaluate(rule, d) if d else ([], 0, 0, "fetch_failed")
            if st != "ok" or len(truth) != len(spans[fam]):
                rows.append({"family": fam, "accession": acc, "status": st if st != "ok"
                             else "module_count_differs"})
                continue
            s, r = maps[m["label"]]
            for i, ((k0, k1, basis), (ta, tb)) in enumerate(zip(spans[fam], truth), 1):
                a, b, cov = cut(s, r, seqs[m["label"]], k0, k1)
                inter = max(0, min(b, tb) - max(a, ta) + 1) if a else 0
                union = max(b, tb) - min(a, ta) + 1 if a else tb - ta + 1
                rows.append({"family": fam, "accession": acc, "module": i,
                             "span_basis": basis, "extracted": f"{a}-{b}",
                             "annotated": f"{ta}-{tb}",
                             "jaccard": round(inter / union, 3),
                             "start_off": a - ta if a else "", "end_off": b - tb if a else "",
                             "status": "ok"})
        print(f"{fam}: {sum(1 for x in rows if x['family'] == fam and x['status'] == 'ok')}"
              f" modules checked", flush=True)
    write_tsv(OUT_DIR / "module_validation.tsv",
              ["family", "accession", "module", "span_basis", "extracted",
               "annotated", "jaccard", "start_off", "end_off", "status"], rows)
    print(f"{fetch.n_requests} requests, {len(fetch.failures)} failures")


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("spans", "extract", "validate"):
        sub.add_parser(c)
    a = ap.parse_args()
    {"spans": cmd_spans, "extract": cmd_extract, "validate": cmd_validate}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
