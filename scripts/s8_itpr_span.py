"""s8_itpr_span.py — the ITPR pore-module span checked against structure (S6 row).

ITPR's module span (S6, D40) came from a RyR-seeded vote plus the helix snap:
the only span with no annotated relative to check it against. Before S8's
Ca²⁺-release tree is read, this compares each human ITPR's extracted module
(`results/alignments/modules.tsv`) with

* UniProt TRANSMEM features (the convention every annotated family's span
  follows: TM5 start → TM6 end), and
* the helices of a cryo-EM structure (PDBe secondary structure, author
  numbering = UniProt numbering for the human chains): 6DQJ, human ITPR3,

and locates the GVGD selectivity-filter motif. Writes
`results/phylogeny/tier2_itpr_span.tsv`; raw JSON archived under
`<data root>/trees/s8/raw/`.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests  # noqa: E402

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s8_lib import ALN_DIR, OUT_DIR, s8_dir  # noqa: E402

STRUCTURES = {"Q14573": "6dqj"}          # human ITPR3, cryo-EM
FILTER = "GVGD"


def _get(name: str, url: str) -> dict:
    f = s8_dir("raw") / name
    if not f.exists():
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        f.write_text(r.text)
    return json.loads(f.read_text())


def helices(pdb: str, lo: int) -> list[tuple[int, int]]:
    d = _get(f"pdbe_ss_{pdb}.json",
             f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/{pdb}")
    ch = d[pdb]["molecules"][0]["chains"][0]
    return [(h["start"]["author_residue_number"], h["end"]["author_residue_number"])
            for h in ch["secondary_structure"].get("helices", [])
            if h["start"]["author_residue_number"] >= lo]


def main() -> int:
    mods = [r for r in read_tsv(ALN_DIR / "modules.tsv")
            if r["family"] == "itpr" and r["label"].endswith("__Homsap")]
    out = []
    for m in mods:
        acc = m["label"].split("__")[0]
        e = _get(f"entry_{acc}.json", f"https://rest.uniprot.org/uniprotkb/{acc}.json")
        seq = e["sequence"]["value"]
        tms = [(f["location"]["start"]["value"], f["location"]["end"]["value"])
               for f in e["features"] if f["type"] == "Transmembrane"]
        tm5, tm6 = tms[-2], tms[-1]
        s, t = int(m["start"]), int(m["end"])
        filt = [x.start() + 1 for x in re.finditer(FILTER, seq) if s <= x.start() + 1 <= t]
        row = {"accession": acc, "module_start": s, "module_end": t,
               "uniprot_tm5": f"{tm5[0]}-{tm5[1]}", "uniprot_tm6": f"{tm6[0]}-{tm6[1]}",
               "start_minus_tm5": s - tm5[0], "end_minus_tm6": t - tm6[1],
               "filter_GVGD": ",".join(map(str, filt)) or "absent",
               "structure": "", "tm5_helix": "", "pore_helix": "", "tm6_helix": "",
               "verdict": ""}
        pdb = STRUCTURES.get(acc)
        if pdb:
            hx = helices(pdb, tm5[0] - 20)
            h5 = next(h for h in hx if h[0] <= tm5[1] and h[1] >= tm5[0])
            hp = [h for h in hx if h5[1] < h[1] and filt and h[1] < filt[0]]
            h6 = [h for h in hx if filt and h[0] > filt[0] and h[0] <= tm6[1]]
            row.update(structure=pdb.upper(), tm5_helix=f"{h5[0]}-{h5[1]}",
                       pore_helix=f"{hp[-1][0]}-{hp[-1][1]}" if hp else "",
                       tm6_helix=f"{h6[0][0]}-{h6[-1][1]}" if h6 else "")
        ok = (abs(row["start_minus_tm5"]) <= 12 and abs(row["end_minus_tm6"]) <= 12
              and filt)
        row["verdict"] = "consistent" if ok else "INCONSISTENT"
        out.append(row)
        print(row, flush=True)
    write_tsv(OUT_DIR / "tier2_itpr_span.tsv", list(out[0]), out)
    return 0 if all(r["verdict"] == "consistent" for r in out) else 1


if __name__ == "__main__":
    raise SystemExit(main())
