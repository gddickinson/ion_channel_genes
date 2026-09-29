"""s6_module_refs.py — S6 step 2a: where the module is, on annotated references.

    python3 scripts/s6_module_refs.py [--refresh]   # → module_refs.tsv

For every census family in a superfamily with a `module_rule` (D40), read
the UniProt topology (TRANSMEM / INTRAMEM features) of reference proteins
and derive their module spans by the declared rule:

* `pore_loop` — for each INTRAMEM pore loop, the nearest TRANSMEM helix
  ending before it through the nearest TRANSMEM helix starting after it
  (S5–P–S6, M1–P–M3, TM5–P–TM6; four per Nav/Cav chain);
* `tm_span` — first TRANSMEM start to last TRANSMEM end (TM1–TM4 of the
  innexin clan, S1–S4 of Hv1).

References, in order: the family's S0-resolved exemplars, then (only if no
exemplar is annotated) up to `MAX_EXTRA` reviewed members of its D39 set.
A family with no annotated reference gets **no** modules and a row saying
so (D28) — its members are never cut by another family's coordinates.
Raw UniProt JSON archived under `<data root>/raw_api/s6/`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s0_lib import UNIPROT, Fetcher  # noqa: E402
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s6_lib import OUT_DIR  # noqa: E402
from src.catalogue import SUPERFAMILIES  # noqa: E402
from src.catalogue.registry import census_families  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

MAX_EXTRA = 8
MIN_MODULE, MAX_MODULE = 60, 400   # a derived span outside this is rejected
FIELDS = ["family", "superfamily", "rule", "reference", "source", "n_tm",
          "n_intramem", "spans", "status"]


def module_families() -> list[tuple[str, str, str]]:
    out = []
    for f in census_families():
        rule = SUPERFAMILIES[f.superfamily].module_rule
        if rule:
            out.append((f.key, f.superfamily, rule))
    return out


def topology(fetch: Fetcher, acc: str, refresh: bool) -> dict | None:
    raw = require_data_root() / "raw_api" / "s6"
    raw.mkdir(parents=True, exist_ok=True)
    p = raw / f"{acc}.json"
    if p.exists() and not refresh:
        return json.loads(p.read_text())
    d = fetch.json(f"{UNIPROT}/uniprotkb/{acc}.json?fields=ft_transmem,"
                   f"ft_intramem,reviewed,length")
    if d is not None:
        p.write_text(json.dumps(d))
    return d


def features(d: dict) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    tm, im = [], []
    for f in d.get("features", []):
        loc = f.get("location", {})
        try:
            a, b = int(loc["start"]["value"]), int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue
        if f.get("type") == "Transmembrane":
            tm.append((a, b))
        elif f.get("type") == "Intramembrane":
            im.append((a, b))
    return sorted(tm), sorted(im)


def spans_for(rule: str, tm: list, im: list) -> list[tuple[int, int]]:
    """The declared rule, applied to one protein's annotated topology."""
    if rule == "tm_span":
        return [(tm[0][0], tm[-1][1])] if tm else []
    # A loop annotated in pieces (pore helix + filter) shares its bracketing
    # helices: one loop, one module.
    out: list[tuple[int, int]] = []
    for a, b in im:
        before = [t for t in tm if t[1] < a]
        after = [t for t in tm if t[0] > b]
        if before and after and (before[-1][0], after[0][1]) not in out:
            out.append((before[-1][0], after[0][1]))
    return out


def evaluate(rule: str, d: dict) -> tuple[list, int, int, str]:
    tm, im = features(d)
    sp = spans_for(rule, tm, im)
    if not sp:
        return [], len(tm), len(im), "no_topology"
    bad = [s for s in sp if not MIN_MODULE <= s[1] - s[0] + 1 <= MAX_MODULE]
    if bad:
        return sp, len(tm), len(im), "span_out_of_range"
    return sp, len(tm), len(im), "ok"


def candidates(fam: str, exemplars: dict, members: dict) -> list[tuple[str, str]]:
    out = [(a, "exemplar") for a in exemplars.get(fam, [])]
    seen = {a for a, _ in out}
    for m in members.get(fam, []):
        if m not in seen:
            out.append((m, "member"))
            seen.add(m)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    ex: dict[str, list[str]] = {}
    for r in read_tsv(ROOT / "results" / "s0_baseline" / "exemplars_resolved.tsv"):
        if r["resolved_accession"]:
            ex.setdefault(r["family"], []).append(r["resolved_accession"])
    mem: dict[str, list[str]] = {}
    for r in read_tsv(OUT_DIR / "members.tsv"):
        if r["verdict"] == "include" and r["target"].startswith("sp|"):
            mem.setdefault(r["family"], []).append(r["target"].split("|")[1])

    fetch, rows = Fetcher(), []
    for fam, sf, rule in module_families():
        n_ok, n_extra = 0, 0
        for acc, src in candidates(fam, ex, mem):
            if src == "member" and (n_ok or n_extra >= MAX_EXTRA):
                break
            n_extra += src == "member"
            d = topology(fetch, acc, a.refresh)
            if d is None:
                rows.append(dict(family=fam, superfamily=sf, rule=rule,
                                 reference=acc, source=src, status="fetch_failed"))
                continue
            sp, ntm, nim, st = evaluate(rule, d)
            n_ok += st == "ok"
            rows.append(dict(family=fam, superfamily=sf, rule=rule, reference=acc,
                             source=src, n_tm=ntm, n_intramem=nim,
                             spans=";".join(f"{x}-{y}" for x, y in sp), status=st))
        if not n_ok:
            rows.append(dict(family=fam, superfamily=sf, rule=rule,
                             status="NO_REFERENCE"))
    write_tsv(OUT_DIR / "module_refs.tsv", FIELDS, rows)
    fams = {r["family"] for r in rows}
    okf = {r["family"] for r in rows if r["status"] == "ok"}
    print(f"{len(okf)}/{len(fams)} families have an annotated reference; "
          f"missing: {sorted(fams - okf)}; {fetch.n_requests} requests, "
          f"{len(fetch.failures)} failures")
    return 1 if fetch.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
