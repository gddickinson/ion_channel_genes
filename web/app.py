"""FastAPI web app for the ion-channel census results dashboard (read-only).

Run with `python run_web.py` (http://127.0.0.1:8180, local only). Every response is read from files git tracks in
this repo, at request time; `scripts/build_static.py` exports every path in STATIC_ENDPOINTS to JSON for the GitHub
Pages site (and copies the figure PNGs served under /files), refusing files with uncommitted changes.
"""
from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from . import docs, headlines, matrix, provenance, stages, tables, trees

STATIC = Path(__file__).resolve().parent / "static"
TITLE = "Ion channel census"

app = FastAPI(title=TITLE)
app.mount("/static", StaticFiles(directory=STATIC), name="static")


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/files/{path:path}", include_in_schema=False)
def files(path: str):
    """Committed figure PNGs only."""
    if path not in docs.figure_paths():
        raise HTTPException(404, "not a committed figure")
    return FileResponse(provenance.ROOT / path)


def _found(v):
    if v is None:
        raise HTTPException(404, "not a committed results file")
    return v


@app.get("/api/info")
def info():
    return {"title": TITLE, "repo": provenance.REPO_URL, "head": provenance.head(),
            "counts": {"tables": len(tables.catalog()), "trees": len(trees.catalog()),
                       "reports": len(docs.report_paths()), "figures": len(docs.figures()),
                       "tasks": len(stages.ledger())}}


@app.get("/api/headlines")
def get_headlines():
    return headlines.payload()


@app.get("/api/ledger")
def ledger():
    return stages.ledger_payload()


@app.get("/api/tables")
def table_index():
    return tables.index_payload()


@app.get("/api/table/{path:path}")
def table(path: str):
    return _found(tables.load(path))


@app.get("/api/matrix")
def presence_matrix():
    return matrix.payload()


@app.get("/api/trees")
def tree_index():
    return trees.index_payload()


@app.get("/api/tree/{path:path}")
def tree(path: str):
    return _found(trees.load(path))


@app.get("/api/figures")
def figures():
    return docs.figures_payload()


@app.get("/api/findings")
def findings():
    return docs.findings_payload()


@app.get("/api/reports")
def reports():
    return docs.reports_payload()


@app.get("/api/report/{path:path}")
def report(path: str):
    return _found(docs.report_payload(path))


def static_endpoints() -> list[str]:
    fixed = ["/info", "/headlines", "/ledger", "/tables", "/matrix", "/trees", "/figures", "/findings", "/reports"]
    return (fixed + [f"/table/{p}" for p in tables.paths()] + [f"/tree/{p}" for p in trees.paths()]
            + [f"/report/{p}" for p in docs.report_paths()])


def served_files() -> list[str]:
    """Every repo file the static site's content is read from (checked clean before export)."""
    return sorted(set(tables.paths() + trees.paths() + docs.report_paths() + docs.figure_paths()
                      + headlines.files() + [matrix.SOURCE, docs.FINDINGS, stages.LEDGER]
                      + [t["stats_table"] for t in trees.catalog() if t["stats_table"]]))


# GET paths (without the /api prefix) written to the static site as api/<path>.json.
STATIC_ENDPOINTS = static_endpoints()
