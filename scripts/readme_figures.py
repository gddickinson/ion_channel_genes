"""readme_figures.py — rewrite README.md's "Results in figures" section from
`scripts/figure_notes.py`, so the README and the dashboard carry the same
plain-English descriptions.

    python3 scripts/readme_figures.py          # rewrite in place
    python3 scripts/readme_figures.py --check  # exit 1 if README is stale
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_notes import FIGURES, caption_markdown  # noqa: E402

README = Path(__file__).resolve().parents[1] / "README.md"
START = "## Results in figures"
END = "## What makes this hard"
INTRO = ("One headline figure per completed task, in task order, each drawn by a "
         "script from that task's committed tables. Every figure has a plain-"
         "English description: what it shows, how it was made, and how to read "
         "its colours and marks. The descriptions live in "
         "`scripts/figure_notes.py` (edit them there, then run "
         "`python3 scripts/readme_figures.py`); the dashboard shows the same text.")


def render() -> str:
    text = README.read_text()
    a, b = text.index(START), text.index(END)
    body = "\n".join(caption_markdown(f) for f in FIGURES)
    return text[:a] + f"{START}\n\n{INTRO}\n\n{body}\n---\n\n" + text[b:]


def main() -> int:
    new = render()
    if "--check" in sys.argv:
        return 0 if new == README.read_text() else 1
    README.write_text(new)
    print(f"{README}: {len(FIGURES)} figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
