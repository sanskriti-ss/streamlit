"""
Phase 1: Genome metadata enrichment for BacDive strain tables.

Usage:
  python -m genome_investigation.genome_enrichment \\
    --input species_data/step3_met_util_exploded.csv.zip \\
    --output genome_investigation/results/Step2_5_genome_enriched.csv \\
    --limit 50
"""

from __future__ import annotations

import argparse
import csv
import time
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from genome_investigation.api_cache import DEFAULT_CACHE_DIR
from genome_investigation.bacdive_client import fetch_strain
from genome_investigation.ncbi_client import enrich_from_ncbi, fetch_assembly_details

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "results" / "Step2_5_genome_enriched.csv"
LOG_DIR = REPO_ROOT / "logs"
FAILURE_LOG = LOG_DIR / "genome_lookup_failures.csv"

OUTPUT_COLUMNS = [
    "BacID",
    "species",
    "strain",
    "genus",
    "order",
    "type_strain",
    "genome_accession",
    "assembly_level",
    "genome_size_bp",
    "gc_percent",
    "gene_count",
    "cds_count",
    "ncbi_taxid",
    "source_database",
    "match_confidence",
    "match_notes",
]

METADATA_COLS = {"BacID", "species", "genus", "order", "type_strain", "is_strain", "species_with_id", "strain"}


def read_input_table(path: Path, limit: Optional[int] = None) -> pd.DataFrame:
    """Load Step1/Step2/Step3 CSV or zip; keep unique strain rows."""
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            csvs = [n for n in zf.namelist() if n.lower().endswith(".csv")]
            if not csvs:
                raise ValueError(f"No CSV in zip: {path}")
            with zf.open(csvs[0]) as f:
                df = pd.read_csv(f, low_memory=False)
    else:
        df = pd.read_csv(path, low_memory=False)

    if "strain" not in df.columns:
        df["strain"] = ""

    keep_cols = [c for c in ["BacID", "species", "strain", "genus", "order", "type_strain"] if c in df.columns]
    if "BacID" not in keep_cols:
        raise ValueError("Input must include BacID column")

    out = df[keep_cols].copy()
    out = out.drop_duplicates(subset=["BacID"], keep="first")
    if limit is not None:
        out = out.head(int(limit))
    return out


def _merge_genome_fields(base: dict, found: dict) -> dict:
    row = {c: base.get(c, "") for c in OUTPUT_COLUMNS}
    for k in OUTPUT_COLUMNS:
        if k in found and found[k] is not None and found[k] != "":
            row[k] = found[k]
    if found.get("strain") and not row.get("strain"):
        row["strain"] = found["strain"]
    if found.get("order") and not row.get("order"):
        row["order"] = found["order"]
    return row


def enrich_row(
    row: dict,
    *,
    cache_dir: Path,
    force_refresh: bool,
    allow_species_only: bool,
    rate_limit_s: float,
) -> dict:
    base = {c: row.get(c, "") for c in OUTPUT_COLUMNS}
    bacid = str(row.get("BacID", "")).strip()

    time.sleep(rate_limit_s)
    bacdive = fetch_strain(bacid, cache_dir=cache_dir, force_refresh=force_refresh)
    if bacdive:
        merged = _merge_genome_fields(base, bacdive)
        acc = merged.get("genome_accession")
        if acc:
            time.sleep(rate_limit_s)
            ncbi_extra = fetch_assembly_details(acc, cache_dir=cache_dir, force_refresh=force_refresh)
            if ncbi_extra:
                for k in ("genome_size_bp", "gc_percent", "gene_count", "cds_count", "assembly_level"):
                    if ncbi_extra.get(k) is not None and merged.get(k) in (None, "", float("nan")):
                        merged[k] = ncbi_extra[k]
                merged["source_database"] = "bacdive+ncbi_datasets"
                merged["match_notes"] = str(merged.get("match_notes", "")) + "; NCBI stats merged by accession"
        return merged

    time.sleep(rate_limit_s)
    ncbi = enrich_from_ncbi(
        row,
        cache_dir=cache_dir,
        force_refresh=force_refresh,
        allow_species_only=allow_species_only,
    )
    if ncbi:
        return _merge_genome_fields(base, ncbi)

    base["match_confidence"] = 0.0
    base["match_notes"] = "no genome match in BacDive or NCBI"
    base["source_database"] = ""
    return base


def enrich_dataframe(
    df: pd.DataFrame,
    *,
    cache_dir: Path,
    force_refresh: bool,
    min_confidence: float,
    allow_species_only: bool,
    rate_limit_s: float = 0.2,
) -> tuple[pd.DataFrame, List[dict]]:
    rows: List[dict] = []
    failures: List[dict] = []

    for i, (_, r) in enumerate(df.iterrows(), start=1):
        row_in = r.to_dict()
        enriched = enrich_row(
            row_in,
            cache_dir=cache_dir,
            force_refresh=force_refresh,
            allow_species_only=allow_species_only,
            rate_limit_s=rate_limit_s,
        )
        conf = float(enriched.get("match_confidence") or 0)
        if conf < min_confidence or not enriched.get("genome_accession"):
            failures.append(
                {
                    "BacID": row_in.get("BacID"),
                    "species": row_in.get("species"),
                    "strain": row_in.get("strain"),
                    "match_confidence": conf,
                    "match_notes": enriched.get("match_notes"),
                }
            )
        rows.append(enriched)
        if i % 25 == 0:
            print(f"[info] enriched {i}/{len(df)} rows …")

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS), failures


def write_failure_log(failures: List[dict], path: Path = FAILURE_LOG) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["BacID", "species", "strain", "match_confidence", "match_notes"]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in failures:
            w.writerow({k: row.get(k, "") for k in fields})


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Phase 1: BacDive/NCBI genome metadata enrichment")
    parser.add_argument("--input", required=True, help="Step1/Step2/Step3 CSV or zip")
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--limit", type=int, default=None, help="Max rows to process")
    parser.add_argument("--force-refresh", action="store_true", help="Ignore API cache")
    parser.add_argument("--min-confidence", type=float, default=0.0, help="Log failures below this")
    parser.add_argument("--allow-species-only", action="store_true", help="Allow weaker NCBI species-only matches")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR))
    args = parser.parse_args(argv)

    in_path = Path(args.input)
    out_path = Path(args.output)
    cache_dir = Path(args.cache_dir)

    print(f"[info] loading {in_path} …")
    df = read_input_table(in_path, limit=args.limit)
    print(f"[info] enriching {len(df)} strains (BacDive → NCBI fallback) …")

    enriched, failures = enrich_dataframe(
        df,
        cache_dir=cache_dir,
        force_refresh=args.force_refresh,
        min_confidence=args.min_confidence,
        allow_species_only=args.allow_species_only,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out_path, index=False)
    write_failure_log(failures)
    n_ok = int((enriched["match_confidence"].astype(float) >= args.min_confidence).sum())
    print(f"[done] wrote {out_path} ({n_ok}/{len(enriched)} rows ≥ min confidence)")
    print(f"[done] failures logged to {FAILURE_LOG} ({len(failures)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
