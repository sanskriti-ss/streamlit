"""
Ingest antiSMASH JSON results into fungal production phenotype predictions.

Works with:
  - Local antiSMASH output dirs (e.g. results/antismash/{species}/{accession}/)
  - Downloaded antiSMASH-DB JSON files (one JSON per genome)

antiSMASH-DB bulk: https://dl.secondarymetabolites.org/database/5.0/

Usage:
  # Parse local antiSMASH output tree
  python -m fungi_investigation.antismash_db_ingest \\
    --json-root results/antismash \\
    --output fungi_data/predicted/fungi_phenotype_confidence.csv

  # Single JSON files directory (flat or nested)
  python -m fungi_investigation.antismash_db_ingest \\
    --json-root /path/to/antismash_jsons \\
    --kingdom fungi
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

import pandas as pd

from genome_investigation.antismash_runner import parse_antismash_json
from genome_investigation.phenotype_confidence import (
    CONF_COLUMNS,
    PRODUCTION_PUTATIVE_CAP,
    production_confidence,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = REPO_ROOT / "fungi_data" / "predicted" / "fungi_phenotype_confidence.csv"


def _iter_antismash_jsons(root: Path) -> Iterator[Tuple[Path, str, str]]:
    """Yield (json_path, species_guess, accession_guess)."""
    patterns = ("genomic.json", "index.json", "*.antismash.json")
    seen: set = set()
    for pattern in patterns:
        for path in sorted(root.rglob(pattern)):
            if path in seen or not path.is_file():
                continue
            seen.add(path)
            parts = path.parts
            accession = ""
            species = ""
            if len(parts) >= 2:
                accession = parts[-2]
                species = parts[-3].replace("_", " ")
            yield path, species, accession


def _species_from_json(path: Path, data: dict) -> str:
    for key in ("species", "organism", "taxon"):
        if isinstance(data.get(key), str) and data[key].strip():
            return data[key].strip()
    records = data.get("records")
    if isinstance(records, list) and records:
        rec = records[0]
        if isinstance(rec, dict):
            name = rec.get("name") or rec.get("description") or ""
            if isinstance(name, str) and name.strip():
                return name.strip()
    stem = path.parent.name
    return stem.replace("_", " ")


def _entity_id(species: str, accession: str, idx: int) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", species.strip().lower())[:40]
    acc = re.sub(r"[^A-Za-z0-9_.-]+", "", accession)[:20]
    return f"FG_{slug}_{acc or idx}"


def ingest_antismash_jsons(
    json_root: Path,
    *,
    kingdom: str = "fungi",
    confidence_cap: float = PRODUCTION_PUTATIVE_CAP,
) -> pd.DataFrame:
    rows: List[dict] = []
    for idx, (json_path, species_guess, accession) in enumerate(_iter_antismash_jsons(json_root)):
        try:
            parsed = parse_antismash_json(json_path)
        except (json.JSONDecodeError, OSError):
            continue

        with json_path.open(encoding="utf-8") as f:
            raw = json.load(f)
        species = species_guess or _species_from_json(json_path, raw if isinstance(raw, dict) else {})
        if not species:
            continue

        bgc_types = parsed.get("bgc_types", "")
        prod_map = production_confidence(bgc_types)
        if not prod_map:
            continue

        genus = species.split()[0] if species else ""
        eid = _entity_id(species, accession, idx)
        for met, conf in prod_map.items():
            rows.append(
                {
                    "entity_id": eid,
                    "entity_key": f"{eid} | {species}",
                    "kingdom": kingdom,
                    "species": species,
                    "genus": genus,
                    "activity": "production",
                    "metabolite": met,
                    "confidence": min(conf, confidence_cap),
                    "evidence": "antiSMASH BGC (predicted)",
                    "observed": False,
                    "genome_accession": accession,
                    "source_json": str(json_path),
                }
            )

    if not rows:
        return pd.DataFrame(columns=CONF_COLUMNS + ["entity_id", "entity_key", "kingdom", "genus"])
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest antiSMASH JSON → fungal production predictions")
    parser.add_argument("--json-root", type=Path, required=True, help="Root directory of antiSMASH JSON output")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--kingdom", default="fungi")
    args = parser.parse_args()

    df = ingest_antismash_jsons(args.json_root, kingdom=args.kingdom)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"[ok] Wrote {len(df):,} production predictions → {args.output}")


if __name__ == "__main__":
    main()
