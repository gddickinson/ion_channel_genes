"""s10_panel_density.py — what a denser panel would cost, for the S10 decision.

    python3 scripts/s10_panel_density.py

S4's emergent row: 52 named species suffice for the census's completeness
argument but are thin for S10's gain/loss reconstruction. This measures the
options before anyone chooses: every UniProt eukaryotic reference proteome
(the same release family S4 drew from), its rank-resolved lineage (UniProt
taxonomy), and — for a density rule of one proteome per phylum, class, order
or family (picked by S4's rule: most proteins as a proxy, then BUSCO) — how
many proteomes and sequences that is, how many of the current 52 species'
lineages it covers, and what the S3b-style sweep would cost at S3b's measured
rate. Raw pages archived under `<data root>/raw_api/s10/`.

→ `results/panel_density/options.tsv`, `proteomes.tsv`, `summary.json`,
`report.md` (rendered from them, D13).
"""

from __future__ import annotations

import csv
import io
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "panel_density"
PROT = "https://rest.uniprot.org/proteomes/stream"
TAX = "https://rest.uniprot.org/taxonomy/search"
RANKS = ("kingdom", "phylum", "class", "order", "family")
#: S3b measured 12.4 min wall for 91 profiles over 822,499 sequences (4 jobs × 2 cpu)
S3B_MIN_PER_SEQ = 12.4 / 822_499


def raw() -> Path:
    p = require_data_root() / "raw_api" / "s10"
    p.mkdir(parents=True, exist_ok=True)
    return p


def fetch_proteomes() -> list[dict]:
    path = raw() / "euk_reference_proteomes.tsv"
    if not path.exists():
        r = requests.get(PROT, params={"query": "(taxonomy_id:2759) AND (reference:true)",
                                       "format": "tsv",
                                       "fields": "upid,organism,organism_id,protein_count,busco"},
                         timeout=300)
        r.raise_for_status()
        path.write_text(r.text)
    return list(csv.DictReader(io.StringIO(path.read_text()), delimiter="\t"))


def fetch_ranks(taxids: list[str]) -> dict[str, dict]:
    path = raw() / "ranks.json"
    have = json.loads(path.read_text()) if path.exists() else {}
    todo = [t for t in taxids if t not in have]
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        q = " OR ".join(f"tax_id:{t}" for t in chunk)
        for attempt in range(5):
            r = requests.get(TAX, params={"query": q, "format": "json", "size": 100},
                             timeout=120)
            if r.status_code == 200:
                break
            time.sleep(2 ** attempt)
        r.raise_for_status()
        for x in r.json().get("results", []):
            lin = {e["rank"]: e["scientificName"] for e in x.get("lineage", [])}
            have[str(x["taxonId"])] = {k: lin.get(k, "") for k in RANKS}
        path.write_text(json.dumps(have))
        time.sleep(0.2)
    return have


def busco(s: str) -> float:
    m = re.match(r"C:([\d.]+)%", s or "")
    return float(m.group(1)) if m else 0.0


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = fetch_proteomes()
    ranks = fetch_ranks([r["Organism Id"] for r in rows])
    for r in rows:
        r.update(ranks.get(r["Organism Id"], {k: "" for k in RANKS}))
        r["busco_c"] = busco(r["BUSCO"])
        r["proteins"] = int(r["Protein count"] or 0)
    write_tsv(OUT / "proteomes.tsv",
              ["Proteome Id", "Organism", "Organism Id", *RANKS, "proteins", "busco_c"], rows)
    panel = [r for r in read_tsv(ROOT / "results" / "proteome_scope" / "proteome_manifest.tsv")
             if r["status"] == "proteome"]
    panel_ranks = fetch_ranks([r["proteome_taxid"] for r in panel])
    euk_panel = [panel_ranks.get(r["proteome_taxid"], {}) for r in panel
                 if r["superregnum"] == "eukaryota"]
    opts = [{"option": "current panel (S4; eukaryotic proteomes)", "rank": "", "proteomes": len(euk_panel),
             "sequences": "", "sweep_hours_91_profiles": ""}]
    for rank in ("phylum", "class", "order", "family"):
        best: dict[str, dict] = {}
        for r in rows:
            k = r[rank]
            if not k:
                continue
            if k not in best or (r["busco_c"], r["proteins"]) > (best[k]["busco_c"], best[k]["proteins"]):
                best[k] = r
        n_seq = sum(r["proteins"] for r in best.values())
        cover = sum(1 for p in euk_panel if p.get(rank) in best)
        by_kingdom = defaultdict(int)
        for r in best.values():
            by_kingdom[r["kingdom"] or "(protists)"] += 1
        opts.append({"option": f"one per {rank}", "rank": rank, "proteomes": len(best),
                     "sequences": n_seq,
                     "sweep_hours_91_profiles": round(n_seq * S3B_MIN_PER_SEQ / 60, 1),
                     "covers_current_species_lineages": f"{cover}/{len(euk_panel)}",
                     "by_kingdom": "; ".join(f"{k} {v}" for k, v in sorted(by_kingdom.items(),
                                                                          key=lambda x: -x[1]))})
    write_tsv(OUT / "options.tsv", ["option", "rank", "proteomes", "sequences",
                                    "sweep_hours_91_profiles",
                                    "covers_current_species_lineages", "by_kingdom"], opts)
    (OUT / "summary.json").write_text(json.dumps(
        {"eukaryotic_reference_proteomes": len(rows),
         "with_rank_resolved": sum(1 for r in rows if r["phylum"]),
         "options": opts}, indent=1))
    render(len(rows), opts)
    for o in opts:
        print(o)
    return 0


def render(n: int, opts: list[dict]) -> None:
    L = ["# Panel density for S10 — the options, measured", "",
         "Rendered by `scripts/s10_panel_density.py` from `options.tsv` (D13).", "",
         f"UniProt holds **{n:,} eukaryotic reference proteomes**. A denser S10 panel "
         "takes one per taxonomic rank, chosen by the best BUSCO completeness (then "
         "most proteins) — S4's rule, restated for many species. The census sweep "
         "scales with sequences; the S3b-rate estimate below is for the profile "
         "sweep alone.", "",
         "| option | proteomes | sequences | profile sweep (h, S3b rate) | kingdoms |",
         "|---|---|---|---|---|"]
    for o in opts:
        seqs = f"{o['sequences']:,}" if o["sequences"] != "" else "—"
        L.append(f"| {o['option']} | {o['proteomes']} | {seqs} | "
                 f"{o['sweep_hours_91_profiles'] or '—'} | {o.get('by_kingdom', '')} |")
    L += ["", "**What does not scale is the genome check.** S5's D4 absence bars "
          "(matched in-group bait, measured detection, contiguity) needed each "
          "genome downloaded and swept: 42.5 Gbp for 52 genomes. At one proteome "
          "per order that is ~400 genomes and hundreds of GB. A denser panel's "
          "absences are therefore **proteome absences** — annotation-level, weaker "
          "than the 52-species panel's controlled absences, which stay the "
          "high-evidence subset.", ""]
    (OUT / "report.md").write_text("\n".join(L))


if __name__ == "__main__":
    raise SystemExit(main())
