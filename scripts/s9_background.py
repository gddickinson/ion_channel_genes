"""s9_background.py — S9 post-hoc control: is the filter more congruent with the
tree than any other column of the alignment the tree was built from?

    python3 scripts/s9_background.py

Added after the D49 tests were read, and labelled so: the filter residues are
inside the alignment each tree was inferred from, so a significant
permutation test (D49 (5)) shows only that the filter carries phylogenetic
signal — which every column does. The comparison that speaks to Q5 is the
filter position's retention index (RI) against the RI of every
parsimony-informative column of the same tree's input alignment (tier 1:
S7's `<fam>.tier1.fasta`, ingroup rows; tier 2: S8's `ploop.tier2.fasta`).
A filter position in a high percentile is conserved along the tree more than
the average site; a low one is re-tuned more often than the average site.

Single positions only, so filter and background are the same kind of
character: a four-repeat chain gives one character per repeat; a K-window
family one per window position (`/`-joined across modules).

Writes `results/filter_atlas/column_background.tsv`.
"""

from __future__ import annotations

import sys
from collections import Counter

from s3_hmm_lib import read_tsv, write_tsv
from s7_lib import s7_dir
from s7_newick import parse
from s8_lib import s8_dir
from s9_congruence import tier1_trees
from s9_lib import OUT_DIR, PHY_DIR, REPEAT_FAMILIES, UNIT, Fitch, is_missing, prune, read_aln

FIELDS = ["tree", "family", "position", "n_tips", "n_states", "ri", "n_columns",
          "column_ri_median", "percentile", "note"]


def ri(fitch: Fitch, states: dict[str, str]) -> float | None:
    known = {t: s for t, s in states.items() if not is_missing(s)}
    c = Counter(known.values())
    if sum(1 for v in c.values() if v >= 2) < 2:      # not parsimony-informative
        return None
    obs = fitch.length(known)
    m, g = len(c) - 1, len(known) - max(c.values())
    return (g - obs) / (g - m)


def column_ris(fitch: Fitch, rows: dict[str, str]) -> list[float]:
    width = len(next(iter(rows.values())))
    out = []
    for j in range(width):
        st = {t: (s[j] if s[j] not in "-.X" else None) for t, s in rows.items()}
        if sum(v is not None for v in st.values()) < 0.5 * len(st):
            continue
        r = ri(fitch, st)
        if r is not None:
            out.append(r)
    return out


def positions(fam: str, strings: dict[str, str]) -> dict[str, dict[str, str]]:
    """Single-position characters from the chain strings."""
    out = {}
    if fam in REPEAT_FAMILIES:
        n = max(len(s) for s in strings.values())
        for i in range(n):
            out[f"repeat {i + 1}"] = {t: (s[i] if len(s) > i else None)
                                      for t, s in strings.items()}
    else:
        for i in range(5):
            out[f"window {i + 1}"] = {t: "/".join(p[i] if len(p) > i else "?"
                                                  for p in s.split("/"))
                                      for t, s in strings.items()}
    return out


def row(tree, fam, pos, fitch, st, bg, note="") -> dict:
    r = ri(fitch, st)
    known = [v for v in st.values() if not is_missing(v)]
    med = sorted(bg)[len(bg) // 2] if bg else None
    return {"tree": tree, "family": fam, "position": pos, "n_tips": len(known),
            "n_states": len(set(known)), "ri": "" if r is None else round(r, 4),
            "n_columns": len(bg), "column_ri_median": "" if med is None else round(med, 4),
            "percentile": "" if r is None or not bg else
            round(100 * sum(x <= r for x in bg) / len(bg), 1), "note": note}


def main() -> None:
    chains = {(r["family"], r["chain"]): r["string"]
              for r in read_tsv(OUT_DIR / "filter_chains.tsv")}
    out = []
    for t in tier1_trees():
        fam = t["family"]
        tree = parse(t["path"].read_text())
        ingroup = set(tree.leaves()) - t["outgroup"]
        aln = {n: s for n, s in read_aln(s7_dir("family") / f"{fam}.tier1.fasta")
               if n in ingroup}
        fitch = Fitch(prune(tree, set(aln)))
        bg = column_ris(fitch, aln)
        strings = {lab: chains.get((fam, lab), "") for lab in aln}
        for pos, st in positions(fam, strings).items():
            out.append(row("tier1", fam, pos, fitch, st, bg))
        print(f"{fam}: {len(bg)} informative columns", flush=True)
    tree = parse((PHY_DIR / "tier2" / f"{UNIT}.treefile").read_text())
    aln = dict(read_aln(s8_dir("unit") / f"{UNIT}.tier2.fasta"))
    fitch = Fitch(prune(tree, set(aln)))
    bg = column_ris(fitch, aln)
    mods = {(r["family"], r["label"]): r for r in read_tsv(OUT_DIR / "filter_modules.tsv")}
    for char in ("locus",) + tuple(f"window {i + 1}" for i in range(5)):
        st = {}
        for tip in aln:
            fam, rest = tip.split("__", 1)
            m = mods.get((fam, rest))
            if not m:
                st[tip] = None
            elif char == "locus":
                st[tip] = m["locus"]
            else:
                w = m["window"]
                st[tip] = w[int(char[-1]) - 1] if len(w) == 5 else None
        out.append(row("tier2", UNIT, char, fitch, st, bg,
                       "locus column = KcsA V76 = window 2" if char == "locus" else ""))
    write_tsv(OUT_DIR / "column_background.tsv", FIELDS, out)


if __name__ == "__main__":
    sys.exit(main())
