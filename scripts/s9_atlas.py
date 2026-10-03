"""s9_atlas.py — S9 driver, the instrument: every P-loop pore module read at the
selectivity filter in one coordinate system (D49).

    python3 scripts/s9_atlas.py add [--jobs 8]   # mafft --add --keeplength per family chunk
    python3 scripts/s9_atlas.py read             # anchors → per-module and per-chain reads
    python3 scripts/s9_atlas.py validate [--jobs 8]   # vs the classifier's pairwise projection

`add` aligns each family's S6 modules (chunks of ≤ CHUNK) to S8's untrimmed
P-loop tier-2 alignment with L-INS-i `--add --keeplength`; resumable on the
chunk's SHA-256. `read` finds the anchor columns (KcsA TVGYG; Nav1.5 DEKA,
after `verify_anchor()`), reads every module there, assembles the chain
strings, and checks every tier-2 tip's re-added copy against its own row.
`validate` projects every four-repeat-family chain onto Nav1.5 by the
classifier's own pairwise method and compares.

Writes `results/filter_atlas/anchors.tsv`, `filter_modules.tsv`,
`filter_chains.tsv`, `tip_check.tsv`, `projection_check.tsv`.
"""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from s3_hmm_lib import read_fasta, read_tsv, write_fasta, write_tsv
from s8_lib import base_label, s8_dir
from s9_lib import (KCSA, KCSA_ROW, NAV15_LABEL, OUT_DIR, REPEAT_FAMILIES, ROOT, UNIT,
                    col_of_residue, module_rows, read_aln, s9_dir)
from src.classify.motifs import FOUR_REPEAT_ANCHOR, K_FILTER_RE, filter_signature, verify_anchor

CHUNK = 60
MAFFT = ["mafft", "--localpair", "--maxiterate", "1000", "--thread", "1", "--quiet",
         "--mapout"]
REF_PANEL = ROOT / "results" / "s0_baseline" / "reference_panel.fasta"


def ref_aln():
    return s8_dir("unit") / f"{UNIT}.linsi.fasta"


def module_fasta() -> dict[str, str]:
    from s6_lib import s6_dir  # noqa: PLC0415
    return read_fasta(s6_dir("modules") / f"{UNIT}.modules.fasta")


def panel_seq(acc: str) -> str:
    for head, seq in read_fasta(REF_PANEL).items():
        if head.endswith("|" + acc):
            return seq
    raise SystemExit(f"{acc} not in {REF_PANEL}")


# ------------------------------------------------------------------ add

def _chunks() -> list[tuple[str, list[tuple[str, str]]]]:
    seqs = module_fasta()
    by_fam = defaultdict(list)
    for key, s in seqs.items():
        by_fam[key.split("|", 1)[0]].append((key, s))
    out = []
    for fam in sorted(by_fam):
        items = sorted(by_fam[fam])
        for i in range(0, len(items), CHUNK):
            out.append((f"{fam}.{i // CHUNK:02d}", items[i:i + CHUNK]))
    return out


def _run_chunk(name: str, items, ref) -> str:
    d = s9_dir("add")
    fa, out = d / f"{name}.fasta", d / f"{name}.aln"
    write_fasta(fa, items)
    sha = hashlib.sha256(fa.read_bytes() + ref.read_bytes()
                         + " ".join(MAFFT).encode()).hexdigest()
    done = d / f"{name}.sha256"
    if out.exists() and done.exists() and done.read_text() == sha and _map(fa).exists():
        return f"{name} cached"
    with open(out, "w") as fh:
        r = subprocess.run(MAFFT + ["--add", str(fa), "--keeplength", str(ref)],
                           stdout=fh, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"mafft failed on {name}: {r.stderr[-400:]}")
    rows = read_aln(out)
    width = {len(s) for _, s in rows}
    if len(rows) != len(read_aln(ref)) + len(items) or len(width) != 1:
        raise RuntimeError(f"{name}: {len(rows)} rows, widths {width}")
    for key, s in items:          # residues unchanged except --keeplength deletions
        got = iter(s)
        if not all(ch in got for ch in dict(rows)[key].replace("-", "")):
            raise RuntimeError(f"{name}: {key} residues changed")
    done.write_text(sha)
    return f"{name} aligned"


def _map(fa):
    return fa.with_name(fa.name + ".map")


def residue_columns(key: str) -> dict[int, int]:
    """MAFFT's `--mapout` for one added sequence: residue number (1-based) →
    reference-alignment column (0-based); residues `--keeplength` deleted are
    absent. Counting residues along the aligned row would be wrong exactly
    where an insertion was deleted before the position."""
    fam = key.split("|", 1)[0]
    for mp in sorted(s9_dir("add").glob(f"{fam}.*.fasta.map")):
        out, on = {}, False
        for line in mp.read_text().splitlines():
            if line.startswith(">"):
                if on:
                    return out
                on = line[1:].split()[0] == key
            elif on and not line.startswith("#"):
                _, pos, col = (x.strip() for x in line.split(","))
                if col != "-":
                    out[int(pos)] = int(col) - 1
        if on:
            return out
    raise KeyError(f"{key}: no --mapout record")


def cmd_add(a) -> None:
    ref = ref_aln()
    chunks = _chunks()
    print(f"{len(chunks)} chunks, {sum(len(c[1]) for c in chunks)} modules", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(_run_chunk, n, it, ref) for n, it in chunks]
        for i, f in enumerate(as_completed(futs), 1):
            print(f"[{i}/{len(futs)}] {f.result()}", flush=True)


# ------------------------------------------------------------------ read

def added_rows() -> dict[str, str]:
    rows = {}
    ref_names = {n for n, _ in read_aln(ref_aln())}
    for p in sorted(s9_dir("add").glob("*.aln")):
        for n, s in read_aln(p):
            if n not in ref_names:
                rows[n] = s
    return rows


def anchors(added: dict[str, str], ref: dict[str, str], mods: dict[str, str]) -> list[dict]:
    acc, first, motif = KCSA
    kcsa = panel_seq(acc)
    if kcsa[first - 1:first - 1 + len(motif)] != motif:
        raise SystemExit(f"KcsA anchor fails: {kcsa[first - 1:first + 4]} at {first}")
    row = ref[KCSA_ROW]
    i = row.replace("-", "").find(motif)
    if i < 0 or row.replace("-", "").count(motif) != 1:
        raise SystemExit("KcsA module row: TVGYG not found exactly once")
    out = [{"anchor": "k_window", "reference": f"KcsA {acc}", "position": first + k,
            "residue": motif[k], "repeat": "", "column": col_of_residue(row, i + 1 + k)}
           for k in range(len(motif))]
    nav = panel_seq(FOUR_REPEAT_ANCHOR.reference_uniprot)
    if not verify_anchor(nav):
        raise SystemExit("verify_anchor() fails on Nav1.5 — refusing to read (D49)")
    for rep, pos in enumerate(FOUR_REPEAT_ANCHOR.positions, 1):
        key = f"nav|{NAV15_LABEL}_m{rep}"
        mod = mods[key]
        start = nav.find(mod) + 1
        if start < 1:
            raise SystemExit(f"{key} is not a substring of Nav1.5")
        idx = pos - start + 1
        if mod[idx - 1] != FOUR_REPEAT_ANCHOR.expected[rep - 1]:
            raise SystemExit(f"{key}: {mod[idx - 1]} at {pos}")
        col = residue_columns(key).get(idx)
        if col is None or added[key][col] != FOUR_REPEAT_ANCHOR.expected[rep - 1]:
            raise SystemExit(f"{key}: anchor residue deleted by --keeplength")
        out.append({"anchor": "repeat_locus", "reference": "Nav1.5 Q14524", "position": pos,
                    "residue": FOUR_REPEAT_ANCHOR.expected[rep - 1], "repeat": rep,
                    "column": col})
    return out


def read_module(s: str, kcols: list[int], lcols: dict[int, int], module: int) -> dict:
    window = "".join(s[c] for c in kcols)
    lc = lcols.get(module, lcols[1])
    return {"window": window, "locus": s[lc], "locus_column": lc}


def cmd_read(_a) -> None:
    ref = dict(read_aln(ref_aln()))
    added = added_rows()
    mods = module_fasta()
    missing = set(mods) - set(added)
    if missing:
        raise SystemExit(f"{len(missing)} modules not aligned, e.g. {sorted(missing)[:3]}")
    anc = anchors(added, ref, mods)
    write_tsv(OUT_DIR / "anchors.tsv", list(anc[0]), anc)
    kcols = [r["column"] for r in anc if r["anchor"] == "k_window"]
    lcols = {r["repeat"]: r["column"] for r in anc if r["anchor"] == "repeat_locus"}
    print(f"K window columns {kcols}; repeat-locus columns {lcols}")
    rows = []
    for r in module_rows():
        key = f"{r['family']}|{r['label']}"
        s = added.get(key)
        hit = K_FILTER_RE.search(mods.get(key, ""))
        rec = {"family": r["family"], "label": r["label"], "chain": base_label(r["label"]),
               "module": r["module"], "quality": r["quality"],
               "k_regex": hit.group(0) if hit else ""}
        rec.update(read_module(s, kcols, lcols, int(r["module"])) if s else
                   {"window": "", "locus": "", "locus_column": ""})
        rec["window_is_k"] = "yes" if K_FILTER_RE.fullmatch(rec["window"] or "") else "no"
        rows.append(rec)
    write_tsv(OUT_DIR / "filter_modules.tsv", list(rows[0]), rows)
    chains = assemble(rows)
    write_tsv(OUT_DIR / "filter_chains.tsv", list(chains[0]), chains)
    tips = tip_check(ref, added, kcols, lcols)
    write_tsv(OUT_DIR / "tip_check.tsv", list(tips[0]), tips)
    agree = sum(t["agree"] == "yes" for t in tips)
    print(f"{len(rows)} modules, {len(chains)} chains; tier-2 tips re-read {agree}/{len(tips)}")


def assemble(rows: list[dict]) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        by[(r["family"], r["chain"])].append(r)
    out = []
    for (fam, chain), ms in sorted(by.items()):
        ms.sort(key=lambda r: int(r["module"]))
        if fam in REPEAT_FAMILIES:
            kind, string = "repeat_locus", "".join(m["locus"] or "?" for m in ms)
        else:
            kind, string = "k_window", "/".join(m["window"] or "?????" for m in ms)
        out.append({"family": fam, "chain": chain, "n_modules": len(ms), "character": kind,
                    "string": string,
                    "quality": ",".join(m["quality"] for m in ms)})
    return out


def tip_check(ref, added, kcols, lcols) -> list[dict]:
    out = []
    for tip, row in ref.items():
        key = tip.replace("__", "|", 1)
        if key not in added:
            continue
        m = int(key.rsplit("_m", 1)[1]) if "_m" in key.rsplit("__", 1)[-1] else 1
        a = read_module(row, kcols, lcols, m)
        b = read_module(added[key], kcols, lcols, m)
        out.append({"tip": tip, "ref_window": a["window"], "added_window": b["window"],
                    "ref_locus": a["locus"], "added_locus": b["locus"],
                    "agree": "yes" if (a["window"], a["locus"]) == (b["window"], b["locus"])
                    else "no"})
    return out


# -------------------------------------------------------------- validate

def cmd_validate(a) -> None:
    from s6_lib import s6_dir  # noqa: PLC0415
    nav = panel_seq(FOUR_REPEAT_ANCHOR.reference_uniprot)
    if not verify_anchor(nav):
        raise SystemExit("verify_anchor() fails on Nav1.5")
    chains = {(r["family"], r["chain"]): r for r in read_tsv(OUT_DIR / "filter_chains.tsv")
              if r["family"] in REPEAT_FAMILIES}
    seqs = {}
    for fam in REPEAT_FAMILIES:
        for lab, s in read_fasta(s6_dir("family") / f"{fam}.fasta").items():
            seqs[(fam, lab)] = s
    jobs = [k for k in chains if k in seqs]
    out = []
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(filter_signature, seqs[k], nav): k for k in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            k = futs[f]
            proj = f.result()
            mod = chains[k]["string"]
            out.append({"family": k[0], "chain": k[1], "module_string": mod,
                        "projection": proj, "agree": compare(mod, proj)})
            if i % 50 == 0:
                print(f"[{i}/{len(jobs)}]", flush=True)
    out.sort(key=lambda r: (r["family"], r["chain"]))
    write_tsv(OUT_DIR / "projection_check.tsv", list(out[0]), out)
    for fam in REPEAT_FAMILIES:
        c = Counter(r["agree"] for r in out if r["family"] == fam)
        print(fam, dict(c))


def compare(mod: str, proj: str) -> str:
    """Four-repeat chains position by position: `yes` (identical), `yes_read`
    (identical wherever both instruments read a residue, one missing a
    repeat), `no`. Chains with fewer than four modules (CatSper, TPC) are
    outside the projection's design (S0: TPC1 → `--PN`) and are recorded as
    `not_comparable` with both strings kept."""
    if len(mod) != 4 or len(proj) != 4:
        return "not_comparable"
    if mod == proj:
        return "yes"
    pairs = [(a, b) for a, b in zip(mod, proj) if a not in "-?X" and b not in "-?X"]
    return "yes_read" if all(a == b for a, b in pairs) else "no"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("add", cmd_add), ("read", cmd_read), ("validate", cmd_validate)):
        p = sub.add_parser(name)
        p.add_argument("--jobs", type=int, default=8)
        p.set_defaults(fn=fn)
    a = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
