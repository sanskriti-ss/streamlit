"""
Ingest FUNG-GROWTH experimental carbon-source utilization data.

FUNG-GROWTH: https://www.fung-growth.org/
Paper: WUR e-depot 701979

Expected raw CSV columns (export from FUNG-GROWTH / BioloMICS or manual curation):
  species, genus, substrate, growth_positive [, entity_id, genome_accession, growth_score]

growth_positive: 1/0, yes/no, +/−, or growth_score > 0 treated as positive.

Usage:
  python -m fungi_investigation.fung_growth_ingest \\
    --input fungi_data/raw/fung_growth/export.csv \\
    --output fungi_data/experimental/fungi_phenotypes_long.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ALIASES = REPO_ROOT / "fungi_data" / "metabolite_aliases.csv"
DEFAULT_OUTPUT = REPO_ROOT / "fungi_data" / "experimental" / "fungi_phenotypes_long.csv"

POSITIVE = {"1", "yes", "y", "true", "+", "positive", "pos"}
NEGATIVE = {"0", "no", "n", "false", "-", "negative", "neg"}


def _is_positive(val) -> bool:
    if pd.isna(val):
        return False
    if isinstance(val, (int, float)):
        return float(val) > 0
    s = str(val).strip().lower()
    if s in POSITIVE:
        return True
    if s in NEGATIVE:
        return False
    try:
        return float(s) > 0
    except ValueError:
        return False


def _slug_entity_id(species: str, idx: int) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", species.strip().lower()).strip("_")
    return f"FG_{slug[:40]}_{idx}"


def load_aliases(path: Path) -> dict:
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    return {
        str(a).strip().lower(): str(c).strip()
        for a, c in zip(df["alias"], df["canonical"])
        if pd.notna(a) and pd.notna(c)
    }


def ingest_fung_growth(
    input_path: Path,
    *,
    aliases_path: Path = DEFAULT_ALIASES,
    include_negative: bool = False,
) -> pd.DataFrame:
    df = pd.read_csv(input_path)
    required = {"species", "substrate"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input missing columns: {sorted(missing)}")

    aliases = load_aliases(aliases_path)
    rows: List[dict] = []
    for idx, r in df.iterrows():
        species = str(r["species"]).strip()
        genus = str(r.get("genus", species.split()[0] if species else "")).strip()
        substrate = str(r["substrate"]).strip()
        met = aliases.get(substrate.lower(), substrate)
        entity_id = str(r.get("entity_id", "") or _slug_entity_id(species, idx))

        if "growth_positive" in df.columns:
            positive = _is_positive(r["growth_positive"])
        elif "growth_score" in df.columns:
            positive = _is_positive(r["growth_score"])
        else:
            raise ValueError("Input needs growth_positive or growth_score column")

        if not positive and not include_negative:
            continue

        rows.append(
            {
                "entity_id": entity_id,
                "entity_key": f"{entity_id} | {species}",
                "kingdom": "fungi",
                "species": species,
                "genus": genus,
                "activity": "utilization",
                "metabolite": met,
                "confidence": 1.0 if positive else 0.0,
                "observed": True,
                "source": "fung_growth",
                "genome_accession": str(r.get("genome_accession", "") or ""),
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest FUNG-GROWTH carbon utilization export")
    parser.add_argument("--input", type=Path, required=True, help="Raw FUNG-GROWTH CSV export")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--aliases", type=Path, default=DEFAULT_ALIASES)
    parser.add_argument(
        "--include-negative",
        action="store_true",
        help="Include negative growth rows (confidence 0.0)",
    )
    args = parser.parse_args()

    out_df = ingest_fung_growth(
        args.input,
        aliases_path=args.aliases,
        include_negative=args.include_negative,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, index=False)
    print(f"[ok] Wrote {len(out_df):,} rows → {args.output}")


if __name__ == "__main__":
    main()
