"""s13_cds.py — a validated CDS for every S13 tip (D55 (4), IP3R D36).

    python3 scripts/s13_cds.py [--workers 4]

dN/dS is a statement about codons, so every tip needs a nucleotide sequence
that provably encodes the exact protein S6 aligned (`<data root>/alignments/
s6/family/<fam>.fasta`). Routes are **candidate generators**; each candidate
is translated and the first that validates (`s13_lib.validate_and_mask`) is
kept, the rejected ones counted:

* proteome tips — UniProt cross-references, in order: Ensembl translation
  and transcript ids (Ensembl REST CDS), EMBL `protein_id` (ENA browser CDS
  FASTA), RefSeq protein (Entrez GenPept `coded_by` → nuccore range);
* genome tips — the S5 call's own CDS blocks (`loci.tsv.gz` `call_cds`)
  spliced from the searched assembly, only where `call_frameshifts` = 0.

Raw fetches are cached under `<data root>/selection/s13/raw/`. Writes
`results/selection/cds_status.tsv` (one row per tip: route, candidates tried,
masked codons, reason on failure) and `<data root>/selection/s13/cds.fasta`
(masked, validated CDS). An anchor without a CDS is a hard failure.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s3_hmm_lib import read_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402
from s5_genome_io import build_fai, fna_path, fetch_region, read_fai  # noqa: E402
from s13_lib import (ANCHORS, OUT, anchor_label, revcomp, s13_dir,  # noqa: E402
                     validate_and_mask)
from src.analysis.selection import fetch_cds  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

FIELDS = ["label", "family", "species", "source", "target", "status", "route",
          "n_candidates", "cds_len", "prot_len", "n_masked", "n_mismatch",
          "n_comparable", "reason"]
PAUSE = 0.34


def _get(url: str, name: str, retries: int = 3) -> str | None:
    cache = s13_dir("raw") / name
    if cache.exists() and cache.stat().st_size:
        return cache.read_text()
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ion-channel-s13"})
            with urllib.request.urlopen(req, timeout=90) as r:
                text = r.read().decode()
            cache.write_text(text)
            time.sleep(PAUSE)
            return text
        except Exception:  # noqa: BLE001
            time.sleep(2.0 * (attempt + 1))
    return None


def _ensembl(eid: str) -> str | None:
    cache = s13_dir("raw") / f"ensembl_{eid}.json"
    if cache.exists():
        return json.loads(cache.read_text()).get("cds")
    for attempt in range(4):
        try:
            cds, _, _ = fetch_cds(eid)
            cache.write_text(json.dumps({"cds": cds}))
            return cds
        except ValueError:
            cache.write_text(json.dumps({"cds": None}))
            return None
        except Exception:  # noqa: BLE001  (HTTP 5xx bursts)
            time.sleep(3.0 * (attempt + 1))
    return None


def _refseq(pid: str) -> str | None:
    eu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    gp = _get(f"{eu}?db=protein&id={pid}&rettype=gp&retmode=text", f"gp_{pid}.txt")
    m = gp and re.search(r'/coded_by="(complement\()?(?:join\()?([A-Z_0-9.]+):[<>]?(\d+)'
                         r'\.\.[<>]?(\d+)', gp)
    if not m:
        return None
    comp, nuc, a, b = m.groups()
    fa = _get(f"{eu}?db=nuccore&id={nuc}&rettype=fasta&retmode=text&seq_start={a}"
              f"&seq_stop={b}" + ("&strand=2" if comp else ""),
              f"nuc_{nuc}_{a}_{b}{'_c' if comp else ''}.fasta")
    return "".join(fa.split("\n")[1:]).strip().upper() if fa and fa.startswith(">") else None


def uniprot_candidates(acc: str):
    """Every CDS the entry cross-references, in preference order."""
    text = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.json"
                "?fields=xref_ensembl,xref_embl,xref_refseq", f"uniprot_{acc}.json")
    xrefs = json.loads(text).get("uniProtKBCrossReferences", []) if text else []
    seen: set[str] = set()
    for x in xrefs:
        if x.get("database") != "Ensembl":
            continue
        ids = [p["value"].split(".")[0] for p in x.get("properties", [])
               if p.get("key") == "ProteinId" and p.get("value", "-") != "-"]
        ids.append(x.get("id", "").split(".")[0])
        for eid in ids:
            if eid and eid not in seen:
                seen.add(eid)
                cds = _ensembl(eid)
                if cds:
                    yield cds, f"ensembl:{eid}"
    for x in xrefs:
        if x.get("database") != "EMBL":
            continue
        pid = next((p["value"] for p in x.get("properties", [])
                    if p.get("key") == "ProteinId"), "-")
        if pid != "-" and pid not in seen:
            seen.add(pid)
            fa = _get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pid}", f"ena_{pid}.fasta")
            if fa and fa.startswith(">"):
                yield "".join(fa.split("\n")[1:]).strip().upper(), f"ena:{pid}"
    for x in xrefs:
        if x.get("database") == "RefSeq" and x.get("id") not in seen:
            seen.add(x["id"])
            cds = _refseq(x["id"])
            if cds:
                yield cds, f"refseq:{x['id']}"


def genome_candidates(target: str):
    """S5's call CDS blocks for `assembly:contig:start-end±`, if frame-intact."""
    asm, locus = target.split(":", 1)
    with gzip.open(require_data_root() / "genomes" / "s5" / asm / "loci.tsv.gz", "rt") as fh:
        row = next((r for r in csv.DictReader(fh, delimiter="\t")
                    if r["locus"] == locus), None)
    if row is None or row["call_frameshifts"] != "0":
        return
    fna = fna_path(asm)
    idx = read_fai(build_fai(fna))
    blocks = sorted((int(a), int(b)) for a, b in
                    (blk.split("-") for blk in row["call_cds"].split(",")))
    seq = "".join(fetch_region(fna, idx, row["contig"], a, b) for a, b in blocks).upper()
    yield (revcomp(seq) if row["strand"] == "-" else seq), f"s5_blocks:{asm}:{locus}"


def one_tip(row: dict, prot: str) -> tuple[dict, str | None]:
    out = {k: row.get(k, "") for k in ("label", "family", "species", "source", "target")}
    if row["source"] == "genome":
        gen = genome_candidates(row["target"])
    else:
        gen = uniprot_candidates(row["target"].split("|")[1])
    n, last = 0, {}
    for cds, route in gen:
        n += 1
        masked, st = validate_and_mask(cds, prot)
        last = st
        if masked:
            out.update(st, status="ok", route=route, n_candidates=n)
            return out, masked
    out.update(last, status="failed", route="", n_candidates=n,
               reason=((last.get("reason") or "no candidate CDS") if n
                       else ("frameshifted locus" if row["source"] == "genome"
                             else "no cross-referenced CDS")))
    return out, None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    tips: dict[str, dict] = {}
    for r in read_tsv(OUT / "sets.tsv"):
        tips.setdefault(r["label"], r)
    prot: dict[str, str] = {}
    for fam in ANCHORS:
        seqs = read_fasta(require_data_root() / "alignments" / "s6" / "family" / f"{fam}.fasta")
        prot.update({k: v.replace("*", "X") for k, v in seqs.items()})
    missing = [t for t in tips if t not in prot]
    if missing:
        raise SystemExit(f"{len(missing)} tips have no S6 sequence, e.g. {missing[:3]}")
    with ThreadPoolExecutor(a.workers) as ex:
        res = list(ex.map(lambda t: one_tip(tips[t], prot[t]), sorted(tips)))
    rows = [r for r, _ in res]
    write_tsv(OUT / "cds_status.tsv", FIELDS, rows)
    write_fasta(s13_dir() / "cds.fasta", [(r["label"], c) for r, c in res if c])
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"{ok}/{len(rows)} tips with a validated CDS; "
          f"{sum(int(r['n_masked'] or 0) for r in rows if r['status'] == 'ok')} codons masked")
    for r in rows:
        if r["status"] != "ok":
            print(f"  FAILED {r['label']}: {r['reason']} ({r['n_candidates']} candidates)")
    bad = [f for f in ANCHORS if anchor_label(f) not in {r["label"] for r in rows if r["status"] == "ok"}]
    if bad:
        raise SystemExit(f"anchor without a validated CDS: {bad} (D55 (4))")


if __name__ == "__main__":
    main()
