"""s5_annotation.py — an assembly's own gene annotation, for calibration.

The pilot measures two thresholds that must not come from miniprot itself,
because `-G` shapes the loci miniprot reports (the parent's S23b lesson):
the **gene span** that sets D4's contiguity bar, and the **widest intron**
that sets `-G`. Both are read from NCBI's annotation of the genes that
family-called loci overlap.

The searched assembly is usually the GenBank (GCA) one, often carrying no
annotation; its RefSeq twin (`pairedAccession`) usually does, on different
sequence names. `sequence_report.jsonl` maps one onto the other.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import zipfile
from bisect import bisect_right
from pathlib import Path

from s5_genome_io import DATASETS, build_fai, fna_path, genome_dir, read_fai

MIN_EXON_FRAC = 0.50


def paired_accession(accession: str) -> str:
    p = genome_dir(accession) / "assembly_data_report.jsonl"
    if not p.exists():
        return ""
    row = json.loads(p.read_text().splitlines()[0])
    return row.get("pairedAccession", "") or ""


def _download_gff(acc: str, dest: Path) -> bool:
    zpath, tmp = dest.parent / "_gff.zip", dest.parent / "_gff_tmp"
    proc = subprocess.run([DATASETS, "download", "genome", "accession", acc,
                           "--include", "gff3", "--filename", str(zpath),
                           "--no-progressbar"], capture_output=True, text=True)
    ok = False
    if proc.returncode == 0:
        shutil.rmtree(tmp, ignore_errors=True)
        with zipfile.ZipFile(zpath) as zf:
            for n in zf.namelist():
                if n.startswith(("/", "..")) or ".." in Path(n).parts:
                    raise RuntimeError(f"unsafe path in zip: {n}")
            zf.extractall(tmp)
        src = tmp / "ncbi_dataset" / "data" / acc / "genomic.gff"
        if src.exists() and src.stat().st_size > 0:
            with open(src, "rb") as fi, gzip.open(dest, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            ok = True
    shutil.rmtree(tmp, ignore_errors=True)
    zpath.unlink(missing_ok=True)
    return ok


def ensure_annotation(accession: str) -> tuple[Path | None, str]:
    """GFF3 for the assembly, else for its paired twin. (path, source acc)."""
    d = genome_dir(accession)
    gz, marker = d / "genomic.gff.gz", d / ".gff.done"
    if marker.exists():
        src = marker.read_text().strip()
        return (gz if gz.exists() else None), (src if gz.exists() else "")
    for acc in (accession, paired_accession(accession)):
        if acc and _download_gff(acc, gz):
            marker.write_text(acc + "\n")
            return gz, acc
    marker.write_text("\n")
    return None, ""


def seqid_map(accession: str) -> dict[str, str]:
    """Any sequence accession (GenBank or RefSeq) → the name the searched
    FASTA actually carries.

    A GCA assembly's FASTA names its sequences by GenBank accession, a GCF
    assembly's by RefSeq accession; mapping onto GenBank names alone left
    every GCF genome with no annotated locus (found in S5b: six genomes).
    """
    p = genome_dir(accession) / "sequence_report.jsonl"
    fna = fna_path(accession)
    in_fasta = set(read_fai(build_fai(fna))) if fna else set()
    out = {}
    if p.exists():
        for line in p.read_text().splitlines():
            r = json.loads(line)
            names = [r.get(k) for k in ("genbankAccession", "refseqAccession")
                     if r.get(k) and r[k] != "na"]
            target = next((n for n in names if n in in_fasta), None)
            for n in names:
                if target:
                    out[n] = target
    return out


def _attr(attrs: str, key: str) -> str:
    for part in attrs.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


def read_genes(gff_gz: Path, names: dict[str, str]) -> list[dict]:
    """Protein-coding genes with span and widest intron (from mRNA exons).

    The widest intron of a gene is the widest over its transcripts: `-G`
    must admit every isoform's structure, not just the longest's.
    """
    genes, mrna_gene, exons, cds = {}, {}, {}, {}
    with gzip.open(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t, a = f[2], f[8]
            if t == "gene":
                gid = _attr(a, "ID")
                genes[gid] = {"contig": names.get(f[0], f[0]), "start": int(f[3]),
                              "end": int(f[4]), "strand": f[6], "gene_id": gid,
                              "name": _attr(a, "Name") or _attr(a, "gene"),
                              "biotype": _attr(a, "gene_biotype"),
                              "max_intron": 0, "n_exons": 0, "exons": []}
            elif t == "mRNA":
                mrna_gene[_attr(a, "ID")] = _attr(a, "Parent")
            elif t == "exon":
                exons.setdefault(_attr(a, "Parent"), []).append((int(f[3]), int(f[4])))
            elif t == "CDS":
                cds.setdefault(_attr(a, "Parent"), []).append((int(f[3]), int(f[4])))
    # Some annotations carry the structure only in CDS: one exon spanning
    # the transcript (Lymnaea, GCA_964033795.1), or no mRNA/exon at all with
    # CDS parented on the gene (Paramecium, GCA_000165425.1). A transcript
    # with at most one exon takes its CDS blocks instead.
    for tid, blocks in cds.items():
        if len(exons.get(tid, [])) <= 1 and len(blocks) > 1:
            exons[tid] = blocks
        if tid in genes and tid not in mrna_gene:
            mrna_gene[tid] = tid
    for tid, blocks in exons.items():
        g = genes.get(mrna_gene.get(tid, ""))
        if g is None:
            continue
        b = sorted(blocks)
        widest = max((b[i + 1][0] - b[i][1] - 1 for i in range(len(b) - 1)), default=0)
        g["max_intron"] = max(g["max_intron"], widest)
        g["n_exons"] = max(g["n_exons"], len(b))
        g["exons"] += b
    for g in genes.values():          # prokaryotic GFFs: gene + CDS, no exons
        if not g["exons"]:
            g["exons"] = [(g["start"], g["end"])]
    return [g for g in genes.values() if g["biotype"] in ("protein_coding", "")]


class GeneIndex:
    def __init__(self, genes: list[dict]):
        self.by_contig: dict[str, list[dict]] = {}
        for g in genes:
            self.by_contig.setdefault(g["contig"], []).append(g)
        for lst in self.by_contig.values():
            lst.sort(key=lambda g: g["start"])
        self.starts = {c: [g["start"] for g in v] for c, v in self.by_contig.items()}

    def overlapping(self, contig: str, start: int, end: int,
                    strand: str = "") -> list[dict]:
        lst = self.by_contig.get(contig, [])
        i = bisect_right(self.starts.get(contig, []), end)
        return [g for g in lst[:i] if g["end"] >= start
                and (not strand or g["strand"] == strand)]

    def best_overlap(self, contig: str, start: int, end: int, strand: str,
                     cds: list[tuple[int, int]] | None = None,
                     min_frac: float = MIN_EXON_FRAC):
        """The same-strand gene whose annotated exons carry ≥ `min_frac` of
        the model's CDS bases (the most, if several).

        Measured in S5a: extent overlap is not confirmation. A 400 kb chained
        ZAC model in mouse "confirmed" against the large gene whose intron it
        sat in, and a reciprocal-extent rule instead rejected real genes whose
        span is mostly 5' UTR introns (mouse 223 → 130 confirmed loci). Exons
        are what an annotation asserts.
        """
        cds = cds or [(start, end)]
        total = sum(e - s + 1 for s, e in cds)
        best, best_bp = None, 0
        for g in self.overlapping(contig, start, end, strand):
            bp = _shared(cds, g["exons"])
            if bp > best_bp:
                best, best_bp = g, bp
        if best is None or best_bp / total < min_frac:
            return None
        return best


def _shared(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    """Bases of `a` covered by the union of `b`."""
    merged = []
    for s, e in sorted(b):
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(max(0, min(e1, e2) - max(s1, s2) + 1)
               for s1, e1 in a for s2, e2 in merged)
