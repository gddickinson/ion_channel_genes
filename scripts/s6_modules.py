"""s6_modules.py — S6 step 2b: one module HMM per tier-2 unit, one method per unit.

    python3 scripts/s6_modules.py seeds     # annotated reference modules → HMM
    python3 scripts/s6_modules.py loo       # leave-one-family-out benchmark

**D40.** A tier-2 unit's modules are all cut by the same instrument, so a
pore-module tree has no extraction-method confounder in it:

1. *Seeds* — the module spans of every UniProt-annotated reference in the
   unit (`s6_module_refs.py`, the superfamily's declared `module_rule`),
   cut from the reference sequences themselves. No projection, no member.
2. *Model* — L-INS-i over the seeds → `hmmbuild` → `<unit>.module.hmm`.
3. *Its use* — measured first as the extraction instrument itself (`loo`)
   and **rejected** (below); it survives as the vote that fixes the module
   span of the six families with no annotated reference
   (`s6_project.py`, which does the extraction).

Measured before use (`loo`): for each family that contributed seeds, the
model is rebuilt without them and that family's members are searched — the
count of modules per chain must equal the catalogue's
`pore_loops_per_subunit` (1 for `tm_span`). That is the evidence the model
reaches the six families with no annotated reference at all.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import Fetcher, sequence  # noqa: E402
from s3_hmm_lib import (DOM_IVALUE_MAX, iter_domtblout, iter_fasta,  # noqa: E402
                        read_tsv, sha256, write_fasta, write_tsv)
from s6_lib import OUT_DIR, s6_dir  # noqa: E402
from s6_module_refs import module_families  # noqa: E402
from src.catalogue import CATALOGUE, SUPERFAMILIES  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

MIN_HMM_COVER = 0.5      # a domain covering less of the model is `partial`


def expected(fam: str) -> int:
    rule = SUPERFAMILIES[CATALOGUE[fam].superfamily].module_rule
    return 1 if rule == "tm_span" else CATALOGUE[fam].pore_loops_per_subunit


def k_filter_family(fam: str) -> bool:
    """A family whose declared filter is a K⁺ filter (…GYG / …GFG)."""
    return bool(re.search(r"G[YFL]G", CATALOGUE[fam].filter_motif or ""))


def units() -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for fam, sf, _ in module_families():
        out[sf].append(fam)
    return dict(out)


def ref_sequence(fetch: Fetcher, acc: str) -> str:
    p = require_data_root() / "raw_api" / "s6" / f"{acc}.fasta"
    if not p.exists():
        s = sequence(fetch, acc)
        if not s:
            raise RuntimeError(f"no sequence for reference {acc}")
        p.write_text(f">{acc}\n{s}\n")
    return next(iter_fasta(p))[1]


def seed_modules(exclude: str = "") -> dict[str, list[tuple[str, str, str]]]:
    """unit → [(label, family, module sequence)] from annotated references."""
    fetch, out = Fetcher(), defaultdict(list)
    for r in read_tsv(OUT_DIR / "module_refs.tsv"):
        if r["status"] != "ok" or r["family"] == exclude:
            continue
        seq = ref_sequence(fetch, r["reference"])
        for i, span in enumerate(r["spans"].split(";"), 1):
            a, b = map(int, span.split("-"))
            out[r["superfamily"]].append(
                (f"{r['reference']}_{r['family']}_m{i}", r["family"], seq[a - 1:b]))
    return dict(out)


def _run(cmd: list[str], stdout=None) -> None:
    p = subprocess.run(cmd, stdout=stdout or subprocess.PIPE,
                       stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{cmd[0]} rc={p.returncode}: {p.stderr.strip()[:300]}")


def build_model(unit: str, seeds: list, d: Path, tag: str = "") -> Path:
    name = f"{unit}{tag}"
    fa, aln, hmm = d / f"{name}.seeds.fasta", d / f"{name}.seeds.aln", d / f"{name}.module.hmm"
    write_fasta(fa, [(l, s) for l, _, s in seeds])
    with open(aln, "w") as out:
        _run(["mafft", "--localpair", "--maxiterate", "1000", "--quiet",
              "--thread", "4", str(fa)], stdout=out)
    _run(["hmmbuild", "--amino", "-n", name, str(hmm), str(aln)])
    return hmm


def unit_members(unit: str, fams: list[str], only: str = "") -> Path:
    fam_dir, d = s6_dir("family"), s6_dir("modules")
    out = d / f"{unit}{'_' + only if only else ''}.members.fasta"
    rows = []
    for f in fams:
        if only and f != only:
            continue
        src = fam_dir / f"{f}.fasta"
        if src.exists():
            rows += [(f"{f}|{h.split()[0]}", s) for h, s in iter_fasta(src)]
    write_fasta(out, rows)
    return out


def search(hmm: Path, db: Path, tbl: Path) -> dict[str, list[dict]]:
    """target → significant domains, N- to C-terminal."""
    _run(["hmmsearch", "--cpu", "4", "--noali", "-o", "/dev/null",
          "--domtblout", str(tbl), str(hmm), str(db)])
    hits: dict[str, list[dict]] = defaultdict(list)
    for r in iter_domtblout(tbl):
        if r["dom_ivalue"] <= DOM_IVALUE_MAX:
            hits[r["target_name"]].append(r)
    for t in hits:
        hits[t].sort(key=lambda r: r["env_from"])
    return hits


def cmd_seeds(_a) -> None:
    d, rows = s6_dir("modules"), []
    for unit, seeds in sorted(seed_modules().items()):
        hmm = build_model(unit, seeds, d)
        fams = Counter(f for _, f, _ in seeds)
        rows.append({"unit": unit, "seeds": len(seeds), "families": len(fams),
                     "seed_families": ",".join(f"{k}:{v}" for k, v in sorted(fams.items())),
                     "model_length": _hmm_len(hmm), "model_sha256": sha256(hmm)})
        print(f"{unit}: {len(seeds)} seeds from {len(fams)} families, "
              f"{rows[-1]['model_length']} match states")
    write_tsv(OUT_DIR / "module_models.tsv", list(rows[0]), rows)


def _hmm_len(p: Path) -> int:
    m = re.search(r"^LENG\s+(\d+)", p.read_text(), re.M)
    return int(m.group(1)) if m else 0


def _per_chain(hits: dict, fasta: Path, fam: str) -> Counter:
    """Distribution of full-coverage modules per chain, against expectation."""
    c = Counter()
    for h, _ in iter_fasta(fasta):
        t = h.split()[0]
        n = sum(1 for r in hits.get(t, [])
                if (r["hmm_to"] - r["hmm_from"] + 1) / r["qlen"] >= MIN_HMM_COVER)
        c["exact" if n == expected(fam) else ("under" if n < expected(fam) else "over")] += 1
    return c


def cmd_loo(_a) -> None:
    d, rows = s6_dir("modules", "loo"), []
    seeded = {r["family"] for r in read_tsv(OUT_DIR / "module_refs.tsv")
              if r["status"] == "ok"}
    for unit, fams in sorted(units().items()):
        for fam in fams:
            loo = fam in seeded
            seeds = seed_modules(exclude=fam if loo else "").get(unit, [])
            if not seeds:
                rows.append({"unit": unit, "family": fam, "held_out": loo,
                             "status": "no_seeds_left"})
                continue
            tag = f"_minus_{fam}" if loo else ""
            hmm = build_model(unit, seeds, d, tag) if loo else s6_dir("modules") / f"{unit}.module.hmm"
            db = unit_members(unit, [fam], only=fam)
            hits = search(hmm, db, d / f"{unit}_{fam}.domtbl")
            c = _per_chain(hits, db, fam)
            n = sum(c.values())
            rows.append({"unit": unit, "family": fam, "held_out": loo,
                         "expected": expected(fam), "members": n,
                         "exact": c["exact"], "under": c["under"], "over": c["over"],
                         "exact_frac": round(c["exact"] / n, 4) if n else "",
                         "status": "ok"})
            print(f"{unit}/{fam} {'LOO' if loo else 'no-ref'}: "
                  f"{c['exact']}/{n} exact ({dict(c)})", flush=True)
    write_tsv(OUT_DIR / "module_loo.tsv",
              ["unit", "family", "held_out", "expected", "members", "exact",
               "under", "over", "exact_frac", "status"], rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("seeds", "loo"):
        sub.add_parser(c)
    a = ap.parse_args()
    {"seeds": cmd_seeds, "loo": cmd_loo}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
