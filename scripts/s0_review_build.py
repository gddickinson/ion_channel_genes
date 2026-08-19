"""Assemble `docs/channel_review_2026.md` from `docs/review/*.md`.

The review is **generated**. Never hand-edit `docs/channel_review_2026.md` —
edit the numbered section files and re-run this, exactly as in the parent
project (`../ip3r_genes/scripts/s0_review_build.py`).

Citations are written in the sections as stable keys — `[hodgkin1952]`,
`[doyle1998, mackinnon2003]` — and resolved against
`results/s0_baseline/references.tsv`. The build renumbers them into order of
first appearance and renders the bibliography, so:

* adding a paper in the middle of section 4 does not renumber the source
  files, only the output;
* **a cited key with no verified reference row is a build error.** That is
  the whole point. The reference table is machine-generated from live Europe
  PMC records (`scripts/s0_review_refs.py`), so a citation that survives the
  build is one a service returned, not one anybody remembered.

    python3 scripts/s0_review_build.py            # write the markdown
    python3 scripts/s0_review_build.py --check    # validate, write nothing
    python3 scripts/s0_review_build.py --pdf      # typeset via pandoc+xelatex
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import read_tsv

SECTIONS = ROOT / "docs" / "review"
OUT_MD = ROOT / "docs" / "channel_review_2026.md"
OUT_PDF = ROOT / "docs" / "channel_review_2026.pdf"
REFS = ROOT / "results" / "s0_baseline" / "references.tsv"

CITE = re.compile(r"\[([a-z][a-z0-9]+(?:\d{4}[a-z0-9]*)(?:\s*,\s*[a-z][a-z0-9]+)*)\]")


def load_refs() -> dict[str, dict]:
    rows = read_tsv(REFS)
    if not rows:
        raise SystemExit(f"[review] no references at {REFS} — run "
                         f"scripts/s0_review_refs.py first")
    return {r["key"]: r for r in rows if r.get("status") == "ok"}


def section_files() -> list[Path]:
    files = sorted(SECTIONS.glob("*.md"))
    if not files:
        raise SystemExit(f"[review] no section files in {SECTIONS}")
    return files


def collect(text: str) -> list[str]:
    """Citation keys in order of first appearance."""
    out: list[str] = []
    for m in CITE.finditer(text):
        for key in (k.strip() for k in m.group(1).split(",")):
            if key not in out:
                out.append(key)
    return out


def renumber(text: str, order: dict[str, int]) -> str:
    def sub(m: re.Match) -> str:
        keys = [k.strip() for k in m.group(1).split(",")]
        nums = sorted(order[k] for k in keys if k in order)
        return "[" + ",".join(str(n) for n in nums) + "]" if nums else m.group(0)
    return CITE.sub(sub, text)


def bibliography(order: dict[str, int], refs: dict[str, dict]) -> str:
    lines = ["## References", ""]
    for key, n in sorted(order.items(), key=lambda kv: kv[1]):
        r = refs[key]
        bits = [f"{n}. {r['authors']}"]
        if r.get("title"):
            bits.append(f"{r['title']}.")
        tail = " ".join(x for x in (r.get("journal", ""), r.get("year", "")) if x)
        if tail:
            bits.append(f"*{tail}*.")
        if r.get("doi"):
            bits.append(f"doi:{r['doi']}")
        elif r.get("pmid"):
            bits.append(f"PMID:{r['pmid']}")
        lines.append(" ".join(bits))
        lines.append("")
    return "\n".join(lines)


def build(check_only: bool = False) -> int:
    refs = load_refs()
    files = section_files()
    body = "\n\n".join(f.read_text().rstrip() for f in files)

    cited = collect(body)
    missing = [k for k in cited if k not in refs]
    if missing:
        print(f"[review] {len(missing)} cited key(s) have no verified "
              f"reference — the build fails rather than dropping them:")
        for k in missing:
            print(f"   [{k}]")
        return 1

    unused = [k for k in refs if k not in cited]
    order = {k: i + 1 for i, k in enumerate(cited)}
    out = renumber(body, order) + "\n\n" + bibliography(order, refs)
    # pandoc/LaTeX housekeeping, same three traps as the parent project
    out = re.sub(r"<sub>(.*?)</sub>", r"~\1~", out)
    out = re.sub(r"<sup>(.*?)</sup>", r"^\1^", out)

    words = len(re.findall(r"\S+", body))
    print(f"[review] {len(files)} sections, {words:,} words, "
          f"{len(cited)} references cited, {len(unused)} resolved but unused")
    if unused:
        print(f"[review]   unused: {', '.join(sorted(unused)[:12])}"
              + (" …" if len(unused) > 12 else ""))
    if check_only:
        print("[review] --check: nothing written")
        return 0
    OUT_MD.write_text(out)
    print(f"[review] wrote {OUT_MD} ({OUT_MD.stat().st_size / 1000:.1f} kB)")
    return 0


def to_pdf() -> int:
    if not OUT_MD.exists():
        print("[review] build the markdown first")
        return 1
    cmd = ["pandoc", str(OUT_MD), "-o", str(OUT_PDF),
           "--pdf-engine=xelatex", "--toc", "--toc-depth=2",
           "-V", "geometry:a4paper,margin=2.2cm",
           "-V", "fontsize=10pt", "-V", "linkcolor=blue",
           "-V", "mainfont=Helvetica Neue"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"[review] pandoc failed: {p.stderr.strip()[:400]}")
        return 1
    print(f"[review] wrote {OUT_PDF} "
          f"({OUT_PDF.stat().st_size / 1e6:.2f} MB)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--pdf", action="store_true")
    args = ap.parse_args()
    rc = build(check_only=args.check)
    if rc or args.check:
        return rc
    return to_pdf() if args.pdf else 0


if __name__ == "__main__":
    sys.exit(main())
