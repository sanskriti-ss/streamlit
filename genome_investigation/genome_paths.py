"""Locate genome files under NCBI Datasets download directories."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple


def find_input_genbank(genome_dir: Path) -> Tuple[Optional[Path], str]:
    """Prefer GBFF/GenBank over FASTA."""
    for pattern in ("*.gbff", "*.gbk", "*.gb", "*.gbff.gz", "*.gbk.gz"):
        hits = sorted(genome_dir.rglob(pattern))
        if hits:
            return hits[0], "genbank"
    for pattern in ("*.fna", "*.fasta", "*.fa"):
        hits = sorted(genome_dir.rglob(pattern))
        if hits:
            return hits[0], "fasta"
    return None, ""


def find_protein_fasta(genome_dir: Path) -> Optional[Path]:
    for pattern in ("*.faa", "*.protein.faa", "*_protein.faa", "*.faa.gz"):
        hits = sorted(genome_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None


def find_nucleotide_fasta(genome_dir: Path) -> Optional[Path]:
    for pattern in ("*_genomic.fna", "*.fna", "*.fasta", "*.fa"):
        hits = sorted(genome_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None


def find_gff3(genome_dir: Path) -> Optional[Path]:
    for pattern in ("*.gff", "*.gff3"):
        hits = sorted(genome_dir.rglob(pattern))
        if hits:
            return hits[0]
    return None
