"""Extract the two selectivity-filter alignments the review's figures show.

Both are measurements, not illustrations, and both are written to committed
tables so the figures are rendered from data rather than drawn by hand.

**The potassium filter.** `TxGYG` is found directly in each sequence by the
regex in `src/classify/motifs.py`, with flanking context. No alignment is
needed and none is used — that is the point of the motif, and a figure built
on it should not smuggle in an aligner.

**The four-repeat locus.** The four residues that separate Nav from Cav from
NALCN cannot be found by a regex; they are projected from human Nav1.5 by
MAFFT (decision **D26**). This writes the projected residue *and its
sequence context* for each of the four positions, so the figure can show the
alignment the call rests on rather than just the four letters.

    python3 scripts/s0_filter_atlas.py [--out results/s0_baseline]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.s0_lib import write_tsv
import re

from src.classify.motifs import FOUR_REPEAT_ANCHOR, K_FILTER_RE, verify_anchor
from src.classify.reference import ReferenceSet
from src.utils.mafft import align, project_positions

PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"
CONTEXT = 5          # residues either side of the filter

#: A deliberately looser pattern, used only when the strict K+ filter fails,
#: so a near miss is recorded as what it is rather than as a blank. The
#: informative case is *Bacillus* NaK, whose filter reads TVGDG: one
#: substitution from the potassium signature, and not potassium-selective.
NEAR_RE = re.compile(r"T[VILMA]G[A-Z]G")

#: Potassium-branch proteins, in the order the figure reads them: the
#: prokaryotic prototypes first, then the eukaryotic families, then the
#: cyanobacterial glutamate receptor whose pore is a K+ channel.
K_PANEL = ["Sl_KcsA", "Mt_MthK", "Bc_NaK", "Dm_Shaker", "Hs_KCNA1",
           "Hs_KCNB1", "Hs_KCND2", "Hs_KCNQ1", "Hs_KCNQ2", "Hs_KCNH2",
           "Hs_KCNH1", "Hs_KCNMA1", "Hs_KCNJ2", "Hs_KCNJ11", "Hs_KCNK2",
           "Hs_KCNK3", "Hs_HCN1", "Hs_HCN4", "Ss_GluR0"]

#: Four-repeat proteins. TPC and CatSper are included precisely because they
#: are expected NOT to give a four-residue call — a figure that showed only
#: the successes would misrepresent the method.
FOUR_PANEL = ["Hs_SCN5A", "Hs_SCN1A", "Hs_CACNA1C", "Hs_CACNA1G",
              "Hs_NALCN", "Hs_TPCN1", "Hs_CATSPER1"]


def k_filter_rows(refs: ReferenceSet) -> list[tuple]:
    rows = []
    for label in K_PANEL:
        seq = refs.sequences.get(label)
        if not seq:
            rows.append((label, 0, "", "", "", "NOT_IN_PANEL"))
            continue
        m = K_FILTER_RE.search(seq)
        if not m:
            near = NEAR_RE.search(seq)
            if near:
                a, b = near.start(), near.end()
                rows.append((label, len(seq), near.group(0),
                             seq[max(0, a - CONTEXT):a], seq[b:b + CONTEXT],
                             "NEAR_MISS"))
            else:
                rows.append((label, len(seq), "", "", "", "NO_MATCH"))
            continue
        a, b = m.start(), m.end()
        rows.append((label, len(seq), m.group(0),
                     seq[max(0, a - CONTEXT):a], seq[b:b + CONTEXT], "ok"))
    return rows


def four_repeat_rows(refs: ReferenceSet) -> list[tuple]:
    ref = refs.sequences.get(FOUR_REPEAT_ANCHOR.reference_label)
    if not ref or not verify_anchor(ref):
        return [("ANCHOR", 0, "", "", "", "", "ANCHOR_FAILS_TO_VALIDATE")]
    rows = []
    for label in FOUR_PANEL:
        seq = refs.sequences.get(label)
        if not seq:
            rows.append((label, 0, "", "", "", "", "NOT_IN_PANEL"))
            continue
        rows_aln = align([("ref", ref), ("qry", seq)])
        proj = project_positions(rows_aln["ref"], rows_aln["qry"],
                                 list(FOUR_REPEAT_ANCHOR.positions))
        # context in the *query*, around each projected residue
        contexts = []
        ri = qi = 0
        want = dict.fromkeys(FOUR_REPEAT_ANCHOR.positions)
        for cr, cq in zip(rows_aln["ref"], rows_aln["qry"]):
            if cr != "-":
                ri += 1
            if cq != "-":
                qi += 1
            if cr != "-" and ri in want and want[ri] is None:
                want[ri] = qi if cq != "-" else 0
        for p in FOUR_REPEAT_ANCHOR.positions:
            q = want.get(p) or 0
            contexts.append(seq[max(0, q - CONTEXT - 1):q + CONTEXT] if q else "")
        sig = "".join(proj[p] for p in FOUR_REPEAT_ANCHOR.positions)
        rows.append((label, len(seq), sig,
                     ";".join(str(want.get(p) or 0)
                              for p in FOUR_REPEAT_ANCHOR.positions),
                     "|".join(contexts), FOUR_REPEAT_ANCHOR.reference_label,
                     "ok" if "-" not in sig else "PARTIAL"))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=ROOT / "results" / "s0_baseline")
    args = ap.parse_args()
    if not PANEL.exists():
        print(f"[filters] no reference panel at {PANEL}")
        return 1
    refs = ReferenceSet.from_fasta(PANEL)

    k = k_filter_rows(refs)
    write_tsv(args.out / "filter_k.tsv",
              ("label", "length_aa", "motif", "left", "right", "status"), k)
    ok_k = sum(1 for r in k if r[5] == "ok")
    near = [r for r in k if r[5] == "NEAR_MISS"]
    print(f"[filters] potassium filter: {ok_k}/{len(k)} matched TxGYG"
          + (f"; {len(near)} near miss(es): "
             + ", ".join(f"{r[0]} {r[2]}" for r in near) if near else ""))

    fr = four_repeat_rows(refs)
    write_tsv(args.out / "filter_four_repeat.tsv",
              ("label", "length_aa", "signature", "positions", "contexts",
               "projected_from", "status"), fr)
    ok_f = sum(1 for r in fr if r[6] == "ok")
    print(f"[filters] four-repeat locus: {ok_f}/{len(fr)} gave a full "
          f"four-residue signature")
    for r in fr:
        print(f"[filters]   {r[0]:14s} {r[2] or '----':6s} {r[6]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
