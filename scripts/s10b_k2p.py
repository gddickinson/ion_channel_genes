"""S10b — what the K2P profile calls outside animals (D51 (5)).

    bin/envpy scripts/s10b_k2p.py

Every non-metazoan high-confidence K2P call on the S4 panel, named by its
UniProt record (protein name, Pfam domains and copies, TM count), plus the
S4b high-confidence K2P orders per kingdom. Names are *read*, never used to
re-code a character: the descent question is S17's.
→ `results/repertoire/k2p_nonmetazoan.tsv`, `k2p_by_kingdom.tsv`.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import s10b_lib as L  # noqa: E402
from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402

FAMILY = "k2p"
OUT = ROOT / "results" / "repertoire"


def uniprot(acc: str) -> dict:
    return L.fetch_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                        L.raw_dir("uniprot") / f"{acc}.json")


def describe(acc: str) -> dict:
    r = uniprot(acc)
    desc = r.get("proteinDescription", {})
    name = (desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]
            ).get("fullName", {}).get("value", "")
    genes = ",".join(g.get("geneName", {}).get("value", "") for g in r.get("genes", [])
                     if g.get("geneName"))
    pf = Counter()
    for x in r.get("uniProtKBCrossReferences", []):
        if x["database"] == "Pfam":
            n = next((p["value"] for p in x.get("properties", [])
                      if p["key"] == "MatchStatus"), "1")
            pf[x["id"]] += int(n) if str(n).isdigit() else 1
    tm = sum(f["type"] == "Transmembrane" for f in r.get("features", []))
    return {"protein_name": name, "gene_names": genes, "reviewed":
            int(r.get("entryType", "").startswith("UniProtKB reviewed")),
            "pfam": ";".join(f"{k}x{v}" for k, v in sorted(pf.items())),
            "tm_features": tm}


def main() -> int:
    s3b = L.data_root() / "hmmer" / "s3b"
    uni = {r["target"]: r for r in read_tsv(s3b / "panel_universe.tsv.gz")}
    rows = []
    for c in read_tsv(s3b / "panel_profile_calls.tsv.gz"):
        if c["p_family"] != FAMILY or c["p_confidence"] != "high":
            continue
        u = uni[c["target"]]
        if L.lineage_of_group(u["group"]) == "Metazoa":
            continue
        rows.append({"accession": u["accession"], "species": u["species"],
                     "group": u["group"], "length": u["length"],
                     "win_score": c["win_score"], "runner": c["runner"],
                     "rel_margin": c["rel_margin"], **describe(u["accession"])})
    rows.sort(key=lambda r: (r["group"], r["species"], r["accession"]))
    write_tsv(OUT / "k2p_nonmetazoan.tsv", list(rows[0]), rows)
    ch = [r for r in read_tsv(OUT / "characters.tsv") if r["family"] == FAMILY]
    kin = Counter((r["kingdom"], r["state"]) for r in ch)
    krows = [{"kingdom": k, "present": kin[(k, "1")], "absent": kin[(k, "0")],
              "missing": kin[(k, "?")]} for k in sorted({r["kingdom"] for r in ch})]
    write_tsv(OUT / "k2p_by_kingdom.tsv", list(krows[0]), krows)
    print(len(rows), "non-metazoan high-confidence K2P calls on the S4 panel")
    for r in rows:
        print(r["species"], r["accession"], r["protein_name"], r["pfam"], r["tm_features"])
    for r in krows:
        print(r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
