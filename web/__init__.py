"""Read-only web dashboard of the committed census results (FastAPI app + static GitHub Pages export).

Contents:
    app.py         FastAPI routes under /api, STATIC_ENDPOINTS for the static export, /files for figure PNGs
    provenance.py  git provenance of every served file (last commit, date, GitHub link; refuses uncommitted files)
    stages.py      roadmap ledger (via scripts/dashboard.parse_ledger) and which ledger task each results file belongs to
    tables.py      the committed results tables (TSV, small JSON summaries), read at request/build time, row-capped
    headlines.py   headline numbers: each one a declared (file, key) read from a committed summary or table
    matrix.py      the S5 presence matrix (census family x panel species, genome verdicts)
    trees.py       committed Newick treefiles with their per-tree statistics rows
    docs.py        FINDINGS.md, the stage reports (results/**/report*.md) and the headline figures with captions
    static/        the plain ES-module front end (no build step)
"""
