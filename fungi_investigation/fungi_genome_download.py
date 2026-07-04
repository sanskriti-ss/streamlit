"""
Download fungal reference genomes via the NCBI Datasets v2 API.

Usage:
  python -m fungi_investigation.fungi_genome_download
  python -m fungi_investigation.fungi_genome_download --species-list fungi_investigation/selected_fungi.yaml
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path
from typing import List, Optional
from urllib.parse import quote
from urllib.request import Request, urlopen

from genome_investigation.io_utils import load_yaml, species_slug

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIST = Path(__file__).resolve().parent / "selected_fungi.yaml"
DEFAULT_OUT = REPO_ROOT / "data" / "genomes_fungi"

NCBI_DOWNLOAD = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
    "?include_annotation_type=PROT_FASTA"
    "&include_annotation_type=GENOME_GBFF"
    "&include_annotation_type=GENOME_GFF"
    "&hydrated=FULLY_HYDRATED"
)


def load_fungi_list(path: Path) -> List[dict]:
    data = load_yaml(path)
    rows = data.get("species", data) if isinstance(data, dict) else data
    out: List[dict] = []
    for row in rows:
        if isinstance(row, str):
            out.append({"species": row, "accession": ""})
        else:
            out.append(
                {
                    "species": str(row["species"]).strip(),
                    "accession": str(row.get("accession", "")).strip(),
                }
            )
    return out


def resolve_accession(species: str) -> str:
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"
        f"{quote(species)}/dataset_report?limit=1"
    )
    req = Request(url, headers={"Accept": "application/json"})
    with urlopen(req, timeout=120) as resp:
        body = json.load(resp)
    reports = body.get("reports") or []
    if not reports:
        raise ValueError(f"No NCBI genome found for {species}")
    report = reports[0]
    return str(report.get("accession") or report.get("current_accession") or "").strip()


def download_genome(accession: str, dest_dir: Path) -> tuple[bool, str]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    meta_path = dest_dir / "download_metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("accession") == accession and any(dest_dir.rglob("*.faa")):
            return True, "already downloaded"

    zip_path = dest_dir / "ncbi_dataset.zip"
    url = NCBI_DOWNLOAD.format(acc=quote(accession, safe=""))
    req = Request(url, headers={"accept": "application/zip"})
    try:
        with urlopen(req, timeout=600) as resp:
            data = resp.read()
    except OSError as exc:
        return False, f"download failed: {exc}"

    if len(data) < 10_000:
        return False, f"download too small ({len(data)} bytes); accession may be invalid"

    zip_path.write_bytes(data)
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(dest_dir)

    meta = {
        "accession": accession,
        "source": "ncbi_datasets_v2_api",
        "zip": str(zip_path),
        "bytes": len(data),
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return True, "downloaded via NCBI Datasets API"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Download fungal reference genomes")
    parser.add_argument("--species-list", type=Path, default=DEFAULT_LIST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    targets = load_fungi_list(args.species_list)
    ok, fail = 0, 0
    for row in targets:
        species = row["species"]
        acc = row.get("accession") or ""
        if not acc:
            try:
                acc = resolve_accession(species)
            except ValueError as exc:
                print(f"[warn] {species}: {exc}")
                fail += 1
                continue

        slug = species_slug(species)
        dest = args.output_dir / slug / acc
        if args.dry_run:
            print(f"[dry] {species} ({acc}) → {dest}")
            continue

        success, msg = download_genome(acc, dest)
        if success:
            ok += 1
            print(f"[ok] {species} ({acc}) → {dest}: {msg}")
        else:
            fail += 1
            print(f"[warn] {species} ({acc}): {msg}")

    if args.dry_run:
        return 0
    print(f"[done] downloaded {ok}, failed {fail}")
    return 0 if fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
