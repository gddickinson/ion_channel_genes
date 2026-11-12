"""Web dashboard tests: every exported endpoint answers, numbers are read from the committed files, provenance and
ledger status are attached, the static build is complete, and the JS tests pass.

Run: python -m pytest -q tests/test_web.py
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web import headlines, matrix, provenance, stages, tables  # noqa: E402
from web.app import STATIC, STATIC_ENDPOINTS, app, served_files  # noqa: E402

client = TestClient(app)


def test_index_served():
    r = client.get("/")
    assert r.status_code == 200 and "static/main.js" in r.text


@pytest.mark.parametrize("path", STATIC_ENDPOINTS)
def test_endpoint(path):
    r = client.get("/api" + path)
    assert r.status_code == 200, r.text
    assert r.json() is not None


def test_only_committed_results_files_are_served():
    assert client.get("/api/table/results/does_not_exist.tsv").status_code == 404
    assert client.get("/api/table/README.md").status_code == 404
    assert client.get("/api/report/CLAUDE.md").status_code == 404
    assert client.get("/files/data_root.txt").status_code == 404
    for path in served_files():
        assert path in provenance.tracked(), path


def test_headlines_equal_the_committed_values():
    hl = {h["id"]: h for h in client.get("/api/headlines").json()["headlines"]}
    s0 = json.loads((ROOT / "results/s0_baseline/verification_summary.json").read_text())
    assert hl["families"]["value"] == s0["catalogue"]["families_census"]
    assert hl["families"]["of"] == s0["catalogue"]["families_total"]
    gs = json.loads((ROOT / "results/genome_sweep/summary.json").read_text())
    assert hl["genomes"]["value"] == gs["genomes"]
    _cols, rows = tables.read_tsv("results/panel_density/order_panel.tsv")
    assert hl["orders"]["value"] == len(rows)
    for h in hl.values():
        assert h["stage"]["id"] == h["task"]
        if not h.get("pending"):
            assert h["provenance"]["committed"] and h["provenance"]["file"].startswith("results/")


def test_pending_headline_follows_the_ledger():
    hl = {h["id"]: h for h in headlines.payload()["headlines"]}
    assert hl["selection"]["pending"] and hl["selection"]["value"] is None
    assert hl["selection"]["stage"]["state"] == stages.badge("S13b")["state"]


def test_ledger_counts_and_states():
    led = client.get("/api/ledger").json()
    assert led["total"] == len(led["tasks"]) == sum(led["counts"].values())
    by_id = {t["id"]: t for t in led["tasks"]}
    assert by_id["S11a"]["state"] == "completed"  # the '|' in its results cell no longer shifts the columns


def test_every_results_table_maps_to_ledger_tasks():
    ids = {t["id"] for t in stages.ledger()}
    for t in tables.catalog():
        assert t["tasks"], f"{t['path']} has no ledger task"
        assert set(t["tasks"]) <= ids, t


def test_table_cap_and_provenance():
    big = client.get("/api/table/results/repertoire/nodes.tsv").json()
    assert big["truncated"] and len(big["rows"]) == tables.MAX_ROWS and big["total_rows"] > tables.MAX_ROWS
    small = client.get("/api/table/results/census_v3/v3_status.tsv").json()
    assert not small["truncated"] and small["provenance"]["url"].startswith(provenance.REPO_URL)


def test_matrix_matches_table():
    m = matrix.payload()
    recs = tables.read_records(matrix.SOURCE)
    assert sum(c is not None for row in m["cells"] for c in row) == len(recs)
    r = recs[0]
    i = [s["name"] for s in m["species"]].index(r["species"])
    j = [f["id"] for f in m["families"]].index(r["family"])
    assert m["verdicts"][m["cells"][i][j][0]]["id"] == r["verdict"]


def test_verdict_scale_matches_repo_figure_scale():
    pytest.importorskip("matplotlib")
    sys.path.insert(0, str(ROOT / "scripts"))
    import s5r4_figures
    repo = [(v, label) for v, label, _c in s5r4_figures.VERDICT]
    ours = [(v, label) for v, label in matrix.VERDICTS]
    assert [v for v, _ in ours] == [v for v, _ in repo]
    assert all(lab == rl for (_, lab), (_, rl) in zip(ours, repo) if rl)


def test_newick_fixture_current():
    sys.path.insert(0, str(ROOT / "scripts"))
    import make_newick_fixture
    want = make_newick_fixture.build()
    have = json.loads((ROOT / "tests/fixtures/newick_fixture.json").read_text())
    assert [t["path"] for t in have["trees"]] == [t["path"] for t in want["trees"]]
    for a, b in zip(have["trees"], want["trees"]):
        assert (a["leaves"], a["internal"], a["ge95"], a["first_leaves"]) == \
            (b["leaves"], b["internal"], b["ge95"], b["first_leaves"])
        assert a["total_length"] == pytest.approx(b["total_length"], rel=1e-6, abs=1e-9)


def test_static_build(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_static
    build_static.build(tmp_path / "site", strict=False)
    site = tmp_path / "site"
    assert (site / "index.html").read_text().count("APP_STATIC") == 1
    for path in STATIC_ENDPOINTS:
        assert (site / "api" / (path.strip("/") + ".json")).exists()
    assert list((site / "files").rglob("*.png"))
    assert not list(site.rglob("*.test.mjs"))
    size = sum(f.stat().st_size for f in site.rglob("*") if f.is_file())
    assert size < 50e6


def test_require_clean_refuses_untracked():
    with pytest.raises(RuntimeError):
        provenance.require_clean(["results/not_a_committed_file.tsv"])


def test_static_assets_referenced_exist():
    html = (STATIC / "index.html").read_text()
    for name in ("style.css", "style_app.css", "main.js", "icon.svg"):
        assert f"static/{name}" in html and (STATIC / name).exists()


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_js_engine():
    for t in sorted(STATIC.glob("*.test.mjs")):
        r = subprocess.run(["node", str(t)], capture_output=True, text=True, env={**os.environ})
        assert r.returncode == 0, r.stdout + r.stderr
