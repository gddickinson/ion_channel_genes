"""Fetch per-residue domain positions for the review's figure exemplars.

`exemplar_architecture.tsv` records *which* domains a protein carries.
Drawing an architecture — a scale diagram of where each domain sits along
the chain — needs *where*, so this fetches the InterPro match coordinates
for a named subset of proteins and writes them to a committed table.

The subset is deliberately small and explicit: these are the proteins the
review's figures compare, and each one is there to make a specific point
(the four-repeat families that domain composition cannot separate, CFTR
against the ABC transporter that regulates a channel, the anoctamin channel
against the anoctamin scramblase, the three families that share `PF08016`).

    python3 scripts/s0_domain_map.py [--out results/s0_baseline]

Output: `domain_positions.tsv` — accession, gene, length, pfam, name, start,
end, copy index. The figures are rendered from that table, never from a
live fetch, so a figure cannot disagree with the data behind it (D13).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import INTERPRO, Fetcher, write_tsv

#: (accession, gene, why it is in a figure)
FIGURE_PROTEINS: list[tuple[str, str, str]] = [
    # the four-repeat families — identical architecture, different ions
    ("Q14524", "SCN5A", "Nav: Ion_trans x4 + Na_trans_assoc"),
    ("Q13936", "CACNA1C", "Cav1: Ion_trans x4 + Ca_chan_IQ"),
    ("O43497", "CACNA1G", "Cav3: Ion_trans x4 and nothing else"),
    ("Q8IZF0", "NALCN", "leak: Ion_trans x4 and nothing else"),
    ("Q8NEC5", "CATSPER1", "one 6TM subunit of a heterotetramer"),
    # the ABC pair
    ("P13569", "CFTR", "ABC fold + the R domain"),
    ("Q09428", "ABCC8", "the same ABC fold, no R domain, not a channel"),
    # the anoctamin pair
    ("Q5XXA6", "ANO1", "channel"),
    ("Q4KMQ2", "ANO6", "scramblase, identical architecture"),
    # the PF08016 trio
    ("Q9GZU1", "MCOLN1", "TRPML"),
    ("Q13563", "PKD2", "TRPP"),
    ("P98161", "PKD1", "polycystin-1: same domain, not a pore"),
    # the shared-signature decoys
    ("Q09470", "KCNA1", "Kv: Ion_trans + the T1 domain"),
    ("Q719H9", "KCTD1", "the T1 domain without a channel"),
    ("P42261", "GRIA1", "iGluR: clamshell + pore"),
    ("Q13255", "GRM1", "the same clamshell, seven-TM, not ionotropic"),
    ("P02708", "CHRNA1", "Cys-loop: LBD + TM"),
    ("P58154", "AChBP", "the same LBD, soluble"),
    ("Q96D96", "HVCN1", "a channel with no pore domain"),
    ("P56180", "TPTE", "the same pore annotation, a phosphatase"),
    # scale
    ("P21817", "RYR1", "5,038 aa — the largest"),
    ("P0A742", "mscL", "136 aa — the smallest"),
]


def positions(f: Fetcher, accession: str) -> list[tuple]:
    d = f.json(f"{INTERPRO}/entry/pfam/protein/uniprot/{accession}/?page_size=100")
    if not d:
        return []
    out = []
    for entry in d.get("results", []):
        m = entry["metadata"]
        name = m.get("name")
        short = name.get("short") if isinstance(name, dict) else name
        for prot in entry.get("proteins", []):
            for i, loc in enumerate(prot.get("entry_protein_locations", []) or [], 1):
                frags = loc.get("fragments") or []
                if not frags:
                    continue
                out.append((m["accession"], short or "",
                            int(frags[0]["start"]), int(frags[-1]["end"]), i))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "s0_baseline")
    args = ap.parse_args()

    f = Fetcher()
    rows = []
    for i, (acc, gene, why) in enumerate(FIGURE_PROTEINS, 1):
        meta = f.json(f"https://rest.uniprot.org/uniprotkb/{acc}.json"
                      f"?fields=accession,length")
        length = ((meta or {}).get("sequence") or {}).get("length", 0)
        hits = positions(f, acc)
        for pfam, short, start, end, copy in sorted(hits, key=lambda h: h[2]):
            rows.append((acc, gene, length, pfam, short, start, end, copy, why))
        print(f"[domains] [{i}/{len(FIGURE_PROTEINS)}] {gene:9s} {acc} "
              f"{length or '?':>5} aa  {len(hits)} match(es)", flush=True)
        time.sleep(0.1)

    write_tsv(args.out / "domain_positions.tsv",
              ("accession", "gene", "length_aa", "pfam", "pfam_name",
               "start", "end", "copy", "role"), rows)
    print(f"[domains] {len(rows)} matches across {len(FIGURE_PROTEINS)} "
          f"proteins → {args.out / 'domain_positions.tsv'}")
    return 1 if f.failures else 0


if __name__ == "__main__":
    sys.exit(main())
