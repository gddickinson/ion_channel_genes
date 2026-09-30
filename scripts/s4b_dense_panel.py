"""s4b_dense_panel.py — S4b: one reference proteome per eukaryotic order, for S10.

    python3 scripts/s4b_dense_panel.py select      # → results/panel_density/order_panel.tsv
    python3 scripts/s4b_dense_panel.py download [-j 4]
    python3 scripts/s4b_dense_panel.py db          # dense_panel.fasta + universe + SHA-256

**The decision (user, 2026-09-30; `results/panel_density/report.md`).** S10's
gain/loss reconstruction needs more lineages than S4's 52 species. The dense
panel takes one UniProt reference proteome per eukaryotic order from the
**census release (2026_03)** — the release README is the list, never the live
API — chosen by rule: the S4 panel's own proteome where the order already has
one (so the two panels nest), else the best BUSCO completeness, then the most
canonical proteins, then the lowest UPID. Every file is the release's
canonical FASTA, its MD5 checked against the proteome's `RELEASE.metalink`
(never waived, D35).

**What it can and cannot say.** Presence and absence here are *proteome*
calls (S3a profiles, D32). S5's controlled absences (in-group bait, measured
detection, D4 contiguity) exist only for the 52 S4 species; a dense-panel
absence is an annotation-level absence and is labelled so in S10.

Bulk: `<data root>/proteomes/s4b/`. Committed: `results/panel_density/`.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, sha256, write_tsv  # noqa: E402
from s4_proteome_lib import (FTP_DOMAIN, FTP_ROOT, download, fetch_metalink,  # noqa: E402
                             fetch_release_readme, md5, parse_readme)
from s10_panel_density import busco, fetch_proteomes, fetch_ranks  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "panel_density"
FIELDS = ["upid", "taxid", "organism", "kingdom", "phylum", "class", "order",
          "busco_c", "n_canonical", "in_s4_panel", "rule"]


def bulk(*parts: str) -> Path:
    p = require_data_root() / "proteomes" / "s4b"
    for part in parts:
        p = p / part
    p.mkdir(parents=True, exist_ok=True)
    return p


def cmd_select(_a) -> None:
    release, readme = parse_readme(fetch_release_readme())
    euk = {u: r for u, r in readme.items() if r["superregnum"] == "eukaryota"}
    api = {r["Proteome Id"]: r for r in fetch_proteomes()}
    ranks = fetch_ranks(sorted({str(r["taxid"]) for r in euk.values()}))
    s4 = {r["upid"] for r in read_tsv(ROOT / "results" / "proteome_scope" /
                                      "proteome_manifest.tsv") if r["status"] == "proteome"}
    by_order: dict[str, list[dict]] = {}
    for u, r in euk.items():
        rk = ranks.get(str(r["taxid"]), {})
        if not rk.get("order"):
            continue
        row = {"upid": u, "taxid": r["taxid"], "organism": r["name"],
               **{k: rk.get(k, "") for k in ("kingdom", "phylum", "class", "order")},
               "busco_c": busco(api.get(u, {}).get("BUSCO", "")),
               "n_canonical": r["n_canonical"], "in_s4_panel": int(u in s4)}
        by_order.setdefault(rk["order"], []).append(row)
    chosen = []
    for order, rows in sorted(by_order.items()):
        mine = [r for r in rows if r["in_s4_panel"]]
        if mine:
            pick, rule = mine[0], "S4 panel proteome"
        else:
            pick = sorted(rows, key=lambda r: (-r["busco_c"], -r["n_canonical"], r["upid"]))[0]
            rule = "best BUSCO, then proteins"
        chosen.append({**pick, "rule": rule})
    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(OUT / "order_panel.tsv", FIELDS, chosen)
    no_order = sum(1 for r in euk.values() if not ranks.get(str(r["taxid"]), {}).get("order"))
    stats = {"release": release, "eukaryotic_reference_proteomes": len(euk),
             "without_order_rank": no_order, "orders": len(chosen),
             "from_s4_panel": sum(r["in_s4_panel"] for r in chosen),
             "canonical_sequences": sum(r["n_canonical"] for r in chosen)}
    (OUT / "order_panel.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps(stats, indent=1))


def _one(r: dict) -> dict:
    upid, taxid = r["upid"], r["taxid"]
    url = f"{FTP_ROOT}/{FTP_DOMAIN['eukaryota']}/{upid}/{upid}_{taxid}.fasta.gz"
    dest = bulk("fasta") / f"{upid}_{taxid}.fasta.gz"
    status = download(url, dest)
    want = fetch_metalink(upid, "eukaryota").get(dest.name, "")
    got = md5(dest) if dest.exists() else ""
    return {"upid": upid, "file": dest.name, "download": status,
            "md5_expected": want, "md5_ok": int(bool(want) and want == got)}


def cmd_download(a) -> None:
    rows = read_tsv(OUT / "order_panel.tsv")
    out = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(_one, r): r["upid"] for r in rows}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                rec = f.result()
            except Exception as exc:                # recorded, never swallowed
                rec = {"upid": futs[f], "file": "", "download": f"failed:{exc}"[:120],
                       "md5_expected": "", "md5_ok": 0}
            out.append(rec)
            if i % 25 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
    write_tsv(OUT / "order_panel_files.tsv",
              ["upid", "file", "download", "md5_expected", "md5_ok"],
              sorted(out, key=lambda r: r["upid"]))
    bad = [r for r in out if not r["md5_ok"]]
    print(f"{len(out) - len(bad)}/{len(out)} files verified; failed: "
          f"{[r['upid'] for r in bad][:10]}")


def cmd_db(_a) -> None:
    rows = {r["upid"]: r for r in read_tsv(OUT / "order_panel.tsv")}
    files = [r for r in read_tsv(OUT / "order_panel_files.tsv") if r["md5_ok"] == "1"]
    db = bulk() / "dense_panel.fasta"
    n = 0
    with open(db, "w") as out, gzip.open(bulk() / "universe.tsv.gz", "wt") as uni:
        uni.write("target\tupid\torder\tclass\tphylum\tkingdom\n")
        for f in sorted(files, key=lambda r: r["upid"]):
            meta = rows[f["upid"]]
            with gzip.open(bulk("fasta") / f["file"], "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        t = line[1:].split()[0]
                        uni.write(f"{t}\t{f['upid']}\t{meta['order']}\t{meta['class']}\t"
                                  f"{meta['phylum']}\t{meta['kingdom']}\n")
                        n += 1
                    out.write(line)
    stats = {"proteomes": len(files), "sequences": n, "db_sha256": sha256(db)}
    (OUT / "dense_db.json").write_text(json.dumps(stats, indent=1))
    print(json.dumps(stats, indent=1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["select", "download", "db"])
    ap.add_argument("-j", "--jobs", type=int, default=4)
    a = ap.parse_args()
    {"select": cmd_select, "download": cmd_download, "db": cmd_db}[a.cmd](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
