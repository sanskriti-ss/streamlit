"""
Predict fungal carbon utilization from protein FASTA annotations.

Reuses keyword matching from genome_investigation.phenotype_confidence (CAZyme-style
headers). Scans genome protein files under data/genomes/ or a custom FASTA root.

For MycoCosm/dbCAN annotations, export protein FASTA or annotation tables and point
--fasta-root at them.

Usage:
  python -m fungi_investigation.fungi_utilization_predict \\
    --fasta-root data/genomes \\
    --output fungi_data/predicted/fungi_utilization_confidence.csv

  # Species list filter (one species per line)
  python -m fungi_investigation.fungi_utilization_predict \\
    --fasta-root data/genomes \\
    --species-list fungi_investigation/selected_fungi.txt
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set

import pandas as pd

from genome_investigation.genome_paths import find_protein_fasta
from genome_investigation.phenotype_confidence import UTILIZATION_CAP, utilization_confidence

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO_ROOT / "fungi_data" / "predicted" / "fungi_utilization_confidence.csv"
DEFAULT_GENOME_ROOT = REPO_ROOT / "data" / "genomes"


def _iter_fasta_jobs(
    root: Path,
    species_filter: Optional[Set[str]] = None,
) -> Iterator[tuple[str, str, Path]]:
    """
    Yield (species_slug, accession, fasta_path) from data/genomes layout or flat *.faa.
    """
    if not root.exists():
        return

    # Layout: root/{species_slug}/{accession}/*.faa
    for species_dir in sorted(root.iterdir()):
        if not species_dir.is_dir():
            continue
        species_slug = species_dir.name
        species_name = species_slug.replace("_", " ")
        if species_filter and species_name not in species_filter and species_slug not in species_filter:
            continue
        for acc_dir in sorted(species_dir.iterdir()):
            if not acc_dir.is_dir():
                continue
            faa = find_protein_fasta(acc_dir)
            if faa:
                yield species_name, acc_dir.name, faa

    # Flat *.faa in root
    for faa in sorted(root.glob("*.faa")) + sorted(root.glob("*.fasta")):
        stem = faa.stem
        species_name = stem.replace("_", " ")
        if species_filter and species_name not in species_filter:
            continue
        yield species_name, "", faa


def _entity_id(species: str, accession: str, idx: int) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", species.strip().lower())[:40]
    acc = re.sub(r"[^A-Za-z0-9_.-]+", "", accession)[:20]
    return f"FG_{slug}_{acc or idx}"


def predict_utilization_from_genomes(
    fasta_root: Path,
    *,
    species_filter: Optional[Set[str]] = None,
    kingdom: str = "fungi",
) -> pd.DataFrame:
    rows: List[dict] = []
    for idx, (species, accession, faa) in enumerate(_iter_fasta_jobs(fasta_root, species_filter)):
        util_map: Dict[str, float] = utilization_confidence(faa)
        if not util_map:
            continue
        genus = species.split()[0] if species else ""
        eid = _entity_id(species, accession, idx)
        for met, conf in util_map.items():
            rows.append(
                {
                    "entity_id": eid,
                    "entity_key": f"{eid} | {species}",
                    "kingdom": kingdom,
                    "species": species,
                    "genus": genus,
                    "activity": "utilization",
                    "metabolite": met,
                    "confidence": min(conf, UTILIZATION_CAP),
                    "evidence": "protein annotation (predicted)",
                    "observed": False,
                    "genome_accession": accession,
                    "fasta_path": str(faa),
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "entity_id",
                "entity_key",
                "kingdom",
                "species",
                "genus",
                "activity",
                "metabolite",
                "confidence",
                "observed",
            ]
        )
    return pd.DataFrame(rows)


def _load_species_list(path: Path) -> Set[str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    return {ln.strip() for ln in lines if ln.strip() and not ln.startswith("#")}


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict fungal utilization from protein FASTA")
    parser.add_argument("--fasta-root", type=Path, default=DEFAULT_GENOME_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--species-list", type=Path, default=None)
    parser.add_argument("--kingdom", default="fungi")
    args = parser.parse_args()

    species_filter = _load_species_list(args.species_list) if args.species_list else None
    df = predict_utilization_from_genomes(
        args.fasta_root,
        species_filter=species_filter,
        kingdom=args.kingdom,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"[ok] Wrote {len(df):,} utilization predictions → {args.output}")


if __name__ == "__main__":
    main()
