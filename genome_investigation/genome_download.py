"""
Phase 1b: Download genomes only for selected species (high-confidence matches).

Usage:
  python -m genome_investigation.genome_download \\
    --selected-species genome_investigation/selected_species.yaml \\
    --genome-enriched genome_investigation/results/Step2_5_genome_enriched.csv \\
    --download-genomes \\
    --min-confidence 0.7
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import List, Optional

import pandas as pd

from genome_investigation.io_utils import load_selected_species, species_slug

DEFAULT_GENOME_DIR = Path(__file__).resolve().parents[1] / "data" / "genomes"
DEFAULT_ENRICHED = Path(__file__).resolve().parent / "results" / "Step2_5_genome_enriched.csv"
SELECTED_YAML = Path(__file__).resolve().parent / "selected_species.yaml"


def _datasets_available() -> bool:
    return shutil.which("datasets") is not None


def filter_selected_rows(
    enriched: pd.DataFrame,
    selected: dict,
    *,
    min_confidence: float,
    max_genomes: Optional[int],
) -> pd.DataFrame:
    species_set = {s.strip().lower() for s in selected.get("species", [])}
    id_set = {str(x).strip() for x in selected.get("bacdive_ids", [])}

    def _match(row) -> bool:
        sp = str(row.get("species", "")).strip().lower()
        bid = str(row.get("BacID", "")).strip()
        if species_set and sp in species_set:
            return True
        if id_set and bid in id_set:
            return True
        return False

    df = enriched[enriched.apply(_match, axis=1)].copy()
    df = df[df["genome_accession"].notna() & (df["genome_accession"].astype(str).str.len() > 3)]
    df["match_confidence"] = pd.to_numeric(df["match_confidence"], errors="coerce").fillna(0)
    df = df[df["match_confidence"] >= min_confidence]
    df = df.sort_values("match_confidence", ascending=False)
    if max_genomes is not None:
        df = df.head(int(max_genomes))
    return df


def download_genome_package(
    accession: str,
    dest_dir: Path,
    *,
    extra_args: Optional[List[str]] = None,
) -> tuple[bool, str]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    meta_path = dest_dir / "download_metadata.json"

    if not _datasets_available():
        return False, "NCBI 'datasets' CLI not found; install from https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"

    zip_path = dest_dir / "ncbi_dataset.zip"
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        str(accession),
        "--filename",
        str(zip_path),
        "--include",
        "genome,gff3,protein,gbff",
    ]
    if extra_args:
        cmd.extend(extra_args)

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as exc:
        err = getattr(exc, "stderr", None) or str(exc)
        return False, f"download failed: {err}"

    if zip_path.exists():
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(dest_dir)
        meta = {"accession": accession, "source": "ncbi_datasets_cli", "zip": str(zip_path)}
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        return True, "downloaded via NCBI datasets CLI"

    return False, "download completed but zip missing"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Phase 1b: selected genome download")
    parser.add_argument("--selected-species", default=str(SELECTED_YAML))
    parser.add_argument("--genome-enriched", default=str(DEFAULT_ENRICHED))
    parser.add_argument("--genome-output-dir", default=str(DEFAULT_GENOME_DIR))
    parser.add_argument("--download-genomes", action="store_true", help="Perform downloads")
    parser.add_argument("--max-genomes", type=int, default=None)
    parser.add_argument("--min-confidence", type=float, default=0.7)
    args = parser.parse_args(argv)

    selected = load_selected_species(Path(args.selected_species))
    enriched_path = Path(args.genome_enriched)
    if not enriched_path.exists():
        print(f"[error] genome enriched file not found: {enriched_path}")
        print("[hint] run genome_enrichment.py first")
        return 1

    enriched = pd.read_csv(enriched_path)
    targets = filter_selected_rows(
        enriched,
        selected,
        min_confidence=args.min_confidence,
        max_genomes=args.max_genomes,
    )
    print(f"[info] {len(targets)} selected strains pass filters")

    if not args.download_genomes:
        print("[info] dry run only (pass --download-genomes to fetch files)")
        for _, r in targets.iterrows():
            slug = species_slug(r["species"])
            print(f"  would download {r['genome_accession']} → data/genomes/{slug}/")
        return 0

    out_root = Path(args.genome_output_dir)
    ok, fail = 0, 0
    for _, r in targets.iterrows():
        acc = str(r["genome_accession"]).strip()
        slug = species_slug(r["species"])
        dest = out_root / slug / acc
        success, msg = download_genome_package(acc, dest)
        if success:
            ok += 1
            print(f"[ok] {r['species']} ({acc}) → {dest}")
        else:
            fail += 1
            print(f"[warn] {r['species']} ({acc}): {msg}")

    print(f"[done] downloaded {ok}, failed {fail}")
    return 0 if fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
