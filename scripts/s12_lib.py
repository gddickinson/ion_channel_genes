"""s12_lib.py — S12 shared pieces: the coverage frame, archived structure
fetchers (AlphaFold DB, UniProt PDB cross-references, PDBe SIFTS), and the
D48 comparison unit computed on any sequence (D54).

Bulk: `<data root>/structures/s12/` (raw API pages, AFDB models, PDBe
updated mmCIF, cut units). Committed: `results/structures/`.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests  # noqa: E402

from s3_hmm_lib import read_tsv, write_tsv  # noqa: E402,F401
from s6_lib import ROOT, iter_census_v4, s3_profiles  # noqa: E402
from s6_project import a2m_map, cut, load_spans  # noqa: E402
from s15_contribution import in_frame  # noqa: E402
from src.catalogue import CATALOGUE, SUPERFAMILIES  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

OUT = ROOT / "results" / "structures"
EXEMPLARS = ROOT / "results" / "s0_baseline" / "exemplars_resolved.tsv"
S8_NET = ROOT / "results" / "phylogeny" / "fold_network"
MIN_PLDDT = 70.0
UNIT_COVER = 0.5
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}


def sdir(*p: str) -> Path:
    """`<data root>/structures/s12/...` — raises without the drive (D1)."""
    d = require_data_root().joinpath("structures", "s12", *p)
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------ frame

def frame_rows() -> list[dict]:
    """S15's final-census frame (D54 (1)), every row incl. genome loci."""
    return [r for r in iter_census_v4() if in_frame(r)]


# ------------------------------------------------------------ fetchers

def _get(url: str, tries: int = 5, **kw) -> requests.Response:
    for i in range(tries):
        try:
            r = requests.get(url, timeout=120, **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"{r.status_code}")
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError):
            if i == tries - 1:
                raise
            time.sleep(2 ** i)
    raise AssertionError


def afdb_record(acc: str) -> dict | None:
    """The AFDB prediction record of exactly this accession (archived)."""
    f = sdir("raw", "afdb") / f"{acc}.json"
    if not f.exists():
        r = _get(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}")
        if r.status_code in (400, 404, 422):
            f.write_text("[]")
        else:
            r.raise_for_status()
            f.write_text(r.text)
    hits = [e for e in json.loads(f.read_text()) if e.get("uniprotAccession") == acc]
    return hits[0] if hits else None


def afdb_model(rec: dict) -> Path:
    pdb = sdir("afdb") / Path(rec["pdbUrl"]).name
    if not pdb.exists():
        r = _get(rec["pdbUrl"])
        r.raise_for_status()
        pdb.write_bytes(r.content)
    return pdb


def uniprot_pdb(accs: list[str]) -> dict[str, list[str]]:
    """accession → PDB ids cross-referenced by UniProt (archived per batch)."""
    out: dict[str, list[str]] = {}
    for i in range(0, len(accs), 400):
        batch = accs[i:i + 400]
        f = sdir("raw", "uniprot_pdb") / f"batch_{i:06d}_{len(batch)}.tsv"
        if not f.exists():
            r = _get("https://rest.uniprot.org/uniprotkb/accessions",
                     params={"accessions": ",".join(batch), "format": "tsv",
                             "fields": "accession,xref_pdb"})
            r.raise_for_status()
            f.write_text(r.text)
        for line in f.read_text().splitlines()[1:]:
            a, _, x = line.partition("\t")
            out[a] = [p for p in x.strip().split(";") if p]
    return out


def sifts_best(acc: str) -> list[dict]:
    """PDBe SIFTS best_structures for an accession (ranked by PDBe; archived)."""
    f = sdir("raw", "sifts") / f"{acc}.json"
    if not f.exists():
        r = _get(f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{acc}")
        f.write_text(r.text if r.status_code == 200 else "{}")
    return json.loads(f.read_text()).get(acc, [])


def updated_cif(pdb_id: str) -> Path:
    f = sdir("pdbe") / f"{pdb_id}_updated.cif.gz"
    if not f.exists():
        r = _get(f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id}_updated.cif")
        r.raise_for_status()
        f.write_bytes(gzip.compress(r.content))
    return f


def experimental_ca(pdb_id: str, chain: str, acc: str) -> list[tuple[int, str, str]]:
    """Observed CA atoms of one chain → (UniProt position, aa, PDB line) by SIFTS.

    The chain is PDBe's author chain id; the label id is tried if it matches
    nothing. A residue SIFTS maps to an isoform is renumbered to the canonical
    sequence by `isoform_to_canonical` (dropped if it has no counterpart)."""
    out = _experimental_ca(pdb_id, chain, acc, "auth_asym_id")
    return out or _experimental_ca(pdb_id, chain, acc, "label_asym_id")


def _experimental_ca(pdb_id: str, chain: str, acc: str, key: str):
    lines = gzip.open(updated_cif(pdb_id), "rt").read().splitlines()
    cols, i = [], 0
    while i < len(lines) and not lines[i].startswith("_atom_site."):
        i += 1
    while i < len(lines) and lines[i].startswith("_atom_site."):
        cols.append(lines[i].split(".", 1)[1].strip())
        i += 1
    ix = {c: n for n, c in enumerate(cols)}
    if "pdbx_sifts_xref_db_acc" not in ix:
        return []                      # no SIFTS residue mapping in this file
    model0, out, seen = None, [], set()
    for line in lines[i:]:
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith(("#", "loop_")):
                break
            continue
        f = line.split()
        if len(f) != len(cols) or f[ix["label_atom_id"]] != "CA" or f[ix[key]] != chain:
            continue
        model = f[ix["pdbx_PDB_model_num"]]
        model0 = model0 or model
        if model != model0 or f[ix["pdbx_sifts_xref_db_acc"]].split("-")[0] != acc:
            continue
        num = f[ix["pdbx_sifts_xref_db_num"]]
        if not num.lstrip("-").isdigit():
            continue
        num = isoform_to_canonical(f[ix["pdbx_sifts_xref_db_acc"]], acc).get(int(num))
        if num is None or num in seen:
            continue
        seen.add(num)
        x, y, z = (float(f[ix[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))
        aa = f[ix["label_comp_id"]]
        pdbline = (f"ATOM  {len(out) + 1:5d}  CA  {aa:3s} A{int(num):4d}    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
        out.append((int(num), AA3.get(aa, "X"), pdbline))
    return out


def isoform_to_canonical(iso: str, acc: str) -> dict[int, int]:
    """Position map isoform → canonical by global alignment (identity for canonical)."""
    if iso in (acc, f"{acc}-1"):
        return _Identity()
    if iso not in _ISO:
        f = sdir("raw", "uniprot") / f"{iso}.fasta"
        if not f.exists():
            r = _get(f"https://rest.uniprot.org/uniprotkb/{iso}.fasta")
            r.raise_for_status()
            f.write_text(r.text)
        iseq = "".join(f.read_text().splitlines()[1:])
        cseq = uniprot_json(acc)["sequence"]["value"]
        from Bio import Align
        al = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=-1,
                                   open_gap_score=-5, extend_gap_score=-0.5)
        aln = al.align(iseq, cseq)[0]
        m = {}
        for (a0, a1), (b0, b1) in zip(*aln.aligned):
            for k in range(a1 - a0):
                if iseq[a0 + k] == cseq[b0 + k]:
                    m[a0 + k + 1] = b0 + k + 1
        _ISO[iso] = m
    return _ISO[iso]


class _Identity(dict):
    def get(self, k, default=None):
        return k


_ISO: dict[str, dict[int, int]] = {}


# ------------------------------------------------------------ units

def module_range(fam: str, seq: str, tag: str) -> tuple[int, int] | None:
    """D48 module 1 on a sequence: hmmalign to the family profile, S6's span."""
    d = sdir("hmm")
    fa, out = d / f"{tag}.fasta", d / f"{tag}.a2m"
    fa.write_text(f">q\n{seq}\n")
    p = subprocess.run(["hmmalign", "--amino", "--outformat", "A2M", "-o", str(out),
                        str(s3_profiles() / f"{fam}.hmm"), str(fa)],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"hmmalign {fam}: {p.stderr.strip()[:300]}")
    a2m = "".join(l.strip() for l in out.read_text().splitlines()[1:])
    st, res = a2m_map(a2m)
    k0, k1, _ = load_spans()[fam][0]
    a, b, cover = cut(st, res, seq, k0, k1)
    return (a, b) if cover >= 0.5 else None


def has_module(fam: str) -> bool:
    return bool(SUPERFAMILIES[CATALOGUE[fam].superfamily].module_rule)


def uniprot_json(acc: str) -> dict:
    f = sdir("raw", "uniprot") / f"{acc}.json"
    if not f.exists():
        r = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
        r.raise_for_status()
        f.write_text(r.text)
    return json.loads(f.read_text())


def tm_region(acc: str) -> tuple[int, int] | None:
    """First TRANSMEM/INTRAMEM start to last end (D54 (5))."""
    feats = [(x["location"]["start"]["value"], x["location"]["end"]["value"])
             for x in uniprot_json(acc).get("features", [])
             if x["type"] in ("Transmembrane", "Intramembrane")]
    feats = [(a, b) for a, b in feats if isinstance(a, int) and isinstance(b, int)]
    return (min(a for a, _ in feats), max(b for _, b in feats)) if feats else None


def read_af_ca(pdb: Path) -> list[tuple[int, str, float, str]]:
    """AFDB model CA atoms → (UniProt position, aa, pLDDT, line)."""
    out = []
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append((int(line[22:26]), AA3.get(line[17:20], "X"),
                        float(line[60:66]), line))
    return out


def write_ca(lines: list[str], out: Path) -> int:
    out.write_text("\n".join(lines) + "\nEND\n")
    return len(lines)
