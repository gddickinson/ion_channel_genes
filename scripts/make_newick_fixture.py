"""Write tests/fixtures/newick_fixture.json: what scripts/s7_newick.py reads from a few committed trees.

The web tree viewer's JS reader (web/static/newick.js) is tested against these values (web/static/newick.test.mjs),
and tests/test_web.py checks the fixture is still current.

Run: python scripts/make_newick_fixture.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import s7_newick  # noqa: E402

TREES = ["results/phylogeny/tier1/asic.treefile", "results/phylogeny/tier1/nachr.treefile",
         "results/phylogeny/tier1/reroot/cng.prok_cng.treefile", "results/phylogeny/tier2/ploop.treefile",
         "results/phylogeny/gost/gost.treefile", "results/duplication/repeats/repeat.treefile"]
OUT = ROOT / "tests" / "fixtures" / "newick_fixture.json"


def summarise(path: str) -> dict:
    root = s7_newick.parse((ROOT / path).read_text())
    internal = ge95 = 0
    total = 0.0
    stack = [(root, True)]
    while stack:
        n, is_root = stack.pop()
        total += n.length or 0.0
        if n.children and not is_root:
            internal += 1
            v = s7_newick.support(n.name)
            ge95 += v is not None and v >= 95
        stack.extend((c, False) for c in n.children)
    leaves = root.leaves()
    return {"path": path, "leaves": len(leaves), "internal": internal, "ge95": int(ge95),
            "total_length": total, "first_leaves": leaves[:5]}


def build() -> dict:
    return {"source": "scripts/s7_newick.py", "trees": [summarise(p) for p in TREES]}


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(build(), indent=1) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
