"""s13_codon.py — codon alignments and trees for every S13 set (D55 (5)–(6)).

    python3 scripts/s13_codon.py

Per alignment unit — each `site` set, and each family's `family` set (its
`og` sets are row subsets of it, so every orthogroup's ω is read on the same
columns):

1. the S6 proteins of the tips with a validated CDS, every masked codon's
   residue written `X` so PAL2NAL is never asked to reconcile a pair that
   disagrees (it would drop the sequence, silently); MAFFT L-INS-i
   `--thread 1`;
2. codons placed by an in-house map (one codon per aligned residue) and by
   PAL2NAL; any difference aborts;
3. columns — `site`: exactly the codons where the anchor has a residue, in
   human UniProt numbering, occupancy recorded per site; `family`: trimAl
   `-gt 0.5` on the protein (D41), applied as whole codons;
4. the S7b tree pruned to the unit's tips and drawn so each orthogroup is a
   clade, tips renamed to short codes (`tip_codes.tsv`) for PAML.

Writes under `results/selection/`: `codon/<set>.fasta` + `.phy`,
`trees/<set>.nwk`, `tip_codes.tsv`, `site_map.tsv` (site set × codon →
human position, residue, occupancy), `codon_stats.tsv`.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s3_hmm_lib import read_fasta, read_tsv, write_fasta, write_tsv  # noqa: E402
from s7_lib import trim_columns  # noqa: E402
from s13_lib import (ANCHORS, OUT, TREES, anchor_label, s13_dir,  # noqa: E402
                     tool_bin, translate)
from s13_tree import UTree, drawn_for, paml_newick, prune  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402


def linsi(rows: list[tuple[str, str]], work: Path) -> dict[str, str]:
    src, out = work / "in.fasta", work / "aln.fasta"
    write_fasta(src, rows)
    p = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "1",
                        "--quiet", str(src)], capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"mafft failed: {p.stderr[:300]}")
    out.write_text(p.stdout)
    aln = read_fasta(out)
    if set(aln) != {k for k, _ in rows}:
        raise SystemExit("mafft lost or renamed rows")
    for k, s in rows:
        if aln[k].replace("-", "").upper() != s.upper():
            raise SystemExit(f"mafft changed the residues of {k}")
    return {k: v.upper() for k, v in aln.items()}


def codon_map(aln: dict[str, str], cds: dict[str, str]) -> dict[str, str]:
    """One codon per aligned residue, '---' per gap."""
    out = {}
    for k, row in aln.items():
        it, buf = iter(range(0, len(cds[k]), 3)), []
        for ch in row:
            buf.append("---" if ch == "-" else cds[k][next(it):][:3])
        if next(it, None) is not None:
            raise SystemExit(f"{k}: CDS longer than its aligned protein")
        out[k] = "".join(buf)
    return out


def pal2nal(aln: dict[str, str], cds: dict[str, str], work: Path) -> dict[str, str]:
    write_fasta(work / "p.fasta", aln.items())
    write_fasta(work / "n.fasta", [(k, cds[k]) for k in aln])
    p = subprocess.run([tool_bin("pal2nal.pl"), str(work / "p.fasta"), str(work / "n.fasta"),
                        "-output", "fasta", "-codontable", "1"],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"pal2nal failed: {p.stderr[:300]}")
    (work / "pal2nal.fasta").write_text(p.stdout)
    out = {k: v.upper() for k, v in read_fasta(work / "pal2nal.fasta").items()}
    if set(out) != set(aln):
        raise SystemExit(f"pal2nal dropped {sorted(set(aln) - set(out))[:5]}")
    return out


def write_phylip(path: Path, rows: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        fh.write(f" {len(rows)} {len(rows[0][1])}\n")
        for k, s in rows:
            fh.write(f"{k}  {s}\n")


def unit(name: str, tips: list[str], fam: str, cds: dict, prot: dict, codes: dict,
         groups: list[frozenset[str]] | None) -> tuple[dict, list[dict]]:
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        prot_x = {}
        for t in tips:   # masked codon → X in the protein given to the aligners
            p = list(prot[t])
            for i in range(len(p)):
                if cds[t][3 * i:3 * i + 3] == "NNN":
                    p[i] = "X"
            prot_x[t] = "".join(p)
        aln = linsi([(t, prot_x[t]) for t in tips], work)
        mine, theirs = codon_map(aln, cds), pal2nal(aln, cds, work)
        for t in tips:
            if mine[t] != theirs[t]:
                raise SystemExit(f"{name}: in-house codon map and PAL2NAL disagree on {t}")
        for t in tips:   # the codon alignment must translate back to the protein
            tr = translate(mine[t].replace("---", ""))
            if any(a != b and b != "X" for a, b in zip(tr, prot_x[t])):
                raise SystemExit(f"{name}: codon row of {t} does not encode its protein")
        sites: list[dict] = []
        if groups is None:   # site set: the anchor's residues, human numbering
            anc = anchor_label(fam)
            cols, pos = [], 0
            for c, ch in enumerate(aln[anc]):
                if ch != "-":
                    pos += 1
                    cols.append(c)
                    occ = sum(1 for t in tips if aln[t][c] not in "-X")
                    sites.append({"set": name, "codon": len(cols), "human_pos": pos,
                                  "human_aa": ch, "occupancy": occ, "n_tips": len(tips)})
        else:
            write_fasta(work / "aln.fasta", aln.items())
            cols = trim_columns(work / "aln.fasta", "gt0.5")
    codon_rows = [(t, "".join(mine[t][3 * c:3 * c + 3] for c in cols)) for t in tips]
    d = OUT / "codon"
    d.mkdir(parents=True, exist_ok=True)
    write_fasta(d / f"{name}.fasta", codon_rows)
    write_phylip(d / f"{name}.phy", [(codes[t], s) for t, s in codon_rows])
    stat = {"set": name, "family": fam, "n_tips": len(tips),
            "aln_cols": len(next(iter(aln.values()))), "codons_kept": len(cols),
            "masked_codons": sum(s.count("NNN") for _, s in codon_rows)}
    return stat, sites


def main() -> None:
    sets = read_tsv(OUT / "sets.tsv")
    ok = {r["label"] for r in read_tsv(OUT / "cds_status.tsv") if r["status"] == "ok"}
    cds = read_fasta(s13_dir() / "cds.fasta")
    prot = {}
    for fam in ANCHORS:
        prot.update({k: v.replace("*", "X") for k, v in read_fasta(
            require_data_root() / "alignments" / "s6" / "family" / f"{fam}.fasta").items()})
    labels = sorted({r["label"] for r in sets})
    codes = {lab: f"T{i + 1:03d}" for i, lab in enumerate(labels)}
    write_tsv(OUT / "tip_codes.tsv", ["code", "label", "has_cds"],
              [{"code": codes[lab], "label": lab, "has_cds": int(lab in ok)} for lab in labels])
    stats, site_rows = [], []
    (OUT / "trees").mkdir(parents=True, exist_ok=True)
    for name in sorted({r["set"] for r in sets}):
        rows = [r for r in sets if r["set"] == name and r["label"] in ok]
        fam, kind = rows[0]["family"], rows[0]["kind"]
        tips = sorted(r["label"] for r in rows)
        tree = prune(UTree.from_newick((TREES / f"{fam}.treefile").read_text()), set(tips))
        if kind == "family":
            ogs: dict[str, set[str]] = {}
            for r in rows:
                ogs.setdefault(r["orthogroup"], set()).add(r["label"])
            groups = [frozenset(g) for g in ogs.values() if len(g) > 0]
            drawn = drawn_for(tree, [g for g in groups if len(g) < len(tips)])
        else:
            groups = None
            drawn = drawn_for(tree, [frozenset(tips)])
        nwk = paml_newick(drawn, {})
        for lab in sorted(tips, key=len, reverse=True):
            nwk = nwk.replace(lab, codes[lab])
        (OUT / "trees" / f"{name}.nwk").write_text(nwk + "\n")
        if kind == "og":
            continue     # row subsets of the family unit, cut in s13_jobs.py
        st, sites = unit(name, tips, fam, cds, prot, codes, groups)
        stats.append(st)
        site_rows += sites
        print(f"{name:18s} {st['n_tips']:>3} tips  {st['aln_cols']:>5} cols → "
              f"{st['codons_kept']:>5} codons  ({st['masked_codons']} NNN)", flush=True)
    write_tsv(OUT / "codon_stats.tsv", list(stats[0]), stats)
    write_tsv(OUT / "site_map.tsv", list(site_rows[0]), site_rows)


if __name__ == "__main__":
    main()
