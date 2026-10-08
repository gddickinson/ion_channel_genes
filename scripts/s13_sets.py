"""s13_sets.py — S13's tip sets, derived from S7b's trees under D55 (2)–(3).

    python3 scripts/s13_sets.py

Per family with an anchor: S7b's tier-1 tree restricted to the family's
vertebrate D39 members; every human tip's **orthogroup** (`s13_tree.
orthogroup`); the **common species** (a tip in every orthogroup). Writes

* `results/selection/orthogroups.tsv` — one row per human tip: size,
  species, bounding UFBoot, tips in the common species;
* `results/selection/sets.tsv` — one row per (set, tip). Set kinds:
  `site` (the anchor's orthogroup, all tips), `og` (one orthogroup at the
  common species), `family` (the union of a family's `og` sets);

Trees are pruned and written by `s13_codon.py`, after the CDS step has
decided which tips stay.

No gene symbol is read: human tips are named by accession, and the
display name of a non-anchor orthogroup is UniProt's entry name taken from
`members.tsv`'s `target`, used as a label only.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402
from s13_lib import (ANCHORS, HUMAN, MEMBERS, OUT, TREES, VERTEBRATE,  # noqa: E402
                     anchor_label, species_code)
from s13_tree import UTree, prune, orthogroup  # noqa: E402

OG_FIELDS = ["family", "human_tip", "entry_name", "is_anchor", "n_tips",
             "n_species", "species", "bounding_ufboot", "n_common_tips"]
SET_FIELDS = ["set", "kind", "family", "orthogroup", "label", "species",
              "source", "target"]


def members(fam: str) -> list[dict]:
    return [r for r in read_tsv(MEMBERS)
            if r["family"] == fam and r["verdict"] == "include"
            and r["group"] == VERTEBRATE]


def entry_name(row: dict) -> str:
    """'sp|P35498|SCN1A_HUMAN' → 'SCN1A' (a display label, never an input)."""
    parts = row["target"].split("|")
    return parts[2].rsplit("_", 1)[0] if len(parts) == 3 else row["label"]


def family_sets(fam: str):
    rows = members(fam)
    by_label = {r["label"]: r for r in rows}
    tree = UTree.from_newick((TREES / f"{fam}.treefile").read_text())
    vert = prune(tree, set(by_label) & tree.leaves())
    if set(by_label) - vert.leaves():
        raise SystemExit(f"{fam}: members missing from the S7b tree")
    humans = sorted(lab for lab in vert.leaves() if species_code(lab) == HUMAN)
    anchor = anchor_label(fam)
    if anchor not in humans:
        raise SystemExit(f"{fam}: anchor {anchor} is not a vertebrate member (D55 (1))")
    ogs = {h: orthogroup(vert, h, set(humans) - {h}) for h in humans}
    common = set.intersection(*[{species_code(x) for x in g} for g, _ in ogs.values()])
    og_rows, set_rows = [], []
    names = {h: entry_name(by_label[h]) for h in humans}

    def add(set_name: str, kind: str, og: str, tips) -> None:
        for lab in sorted(tips):  # a `family` row's orthogroup is its og's name
            r = by_label[lab]
            set_rows.append({"set": set_name, "kind": kind, "family": fam,
                             "orthogroup": og, "label": lab,
                             "species": r["species"], "source": r["source"],
                             "target": r["target"]})

    for h in humans:
        g, sup = ogs[h]
        in_common = sorted(x for x in g if species_code(x) in common)
        og_rows.append({"family": fam, "human_tip": h, "entry_name": names[h],
                        "is_anchor": int(h == anchor), "n_tips": len(g),
                        "n_species": len({species_code(x) for x in g}),
                        "species": ",".join(sorted({species_code(x) for x in g})),
                        "bounding_ufboot": "" if sup is None else sup,
                        "n_common_tips": len(in_common) if len(humans) > 1 else ""})
    site = ogs[anchor][0]
    add(f"site_{names[anchor]}", "site", names[anchor], site)
    if len(humans) > 1:
        for h in humans:
            tips = [x for x in ogs[h][0] if species_code(x) in common]
            add(f"og_{names[h]}", "og", names[h], tips)
            add(f"family_{fam}", "family", names[h], tips)
    return og_rows, set_rows, sorted(common)


def main() -> None:
    og_all, set_all = [], []
    for fam in ANCHORS:
        og_rows, set_rows, common = family_sets(fam)
        og_all += og_rows
        set_all += set_rows
        a = next(r for r in og_rows if r["is_anchor"])
        print(f"{fam:8s} anchor {a['entry_name']:6s} {a['n_tips']:>3} tips "
              f"(UFBoot {a['bounding_ufboot'] or '—'}); {len(og_rows)} orthogroups; "
              f"common species {len(common)}: {','.join(common)}")
    write_tsv(OUT / "orthogroups.tsv", OG_FIELDS, og_all)
    write_tsv(OUT / "sets.tsv", SET_FIELDS, set_all)
    kinds: dict[str, int] = {}
    for r in set_all:
        kinds[r["set"]] = kinds.get(r["set"], 0) + 1
    print(f"{len(kinds)} sets, {len({r['label'] for r in set_all})} distinct tips")


if __name__ == "__main__":
    main()
