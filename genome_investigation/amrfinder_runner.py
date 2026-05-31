"""
Run NCBI AMRFinderPlus on downloaded genomes and summarize hits.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

import pandas as pd

from genome_investigation.genome_paths import find_nucleotide_fasta, find_protein_fasta
from genome_investigation.io_utils import load_selected_species, species_slug

DEFAULT_GENOME_DIR = Path(__file__).resolve().parents[1] / "data" / "genomes"
DEFAULT_SUMMARY = Path(__file__).resolve().parent / "results" / "amrfinder_summary.csv"
SELECTED_YAML = Path(__file__).resolve().parent / "selected_species.yaml"

AMR_SUMMARY_COLUMNS = [
    "BacID",
    "species",
    "strain",
    "genome_accession",
    "amrfinder_status",
    "amrfinder_notes",
    "amr_gene_count",
    "amr_classes",
    "amr_genes",
    "amr_output_path",
]


def amrfinder_installed() -> bool:
    return shutil.which("amrfinder") is not None


def run_amrfinder(
    input_path: Path,
    output_path: Path,
    *,
    protein: bool = False,
    organism: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
) -> tuple[bool, str]:
    if not amrfinder_installed():
        return False, "amrfinder executable not found (install AMRFinderPlus)"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "amrfinder",
        "-o",
        str(output_path),
        "--report_all_equal",
        "--log",
        str(output_path.with_suffix(".log")),
    ]
    if protein:
        cmd.extend(["-p", str(input_path)])
    else:
        cmd.extend(["-n", str(input_path)])
    if organism:
        cmd.extend(["--organism", organism])
    if extra_args:
        cmd.extend(extra_args)

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
        return True, "completed"
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as exc:
        err = getattr(exc, "stderr", None) or str(exc)
        return False, f"AMRFinder failed: {err}"


def parse_amrfinder_tsv(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {
            "amr_gene_count": 0,
            "amr_classes": "",
            "amr_genes": "",
        }

    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.empty:
        return {"amr_gene_count": 0, "amr_classes": "", "amr_genes": ""}

    gene_col = "Element symbol" if "Element symbol" in df.columns else "Gene symbol"
    class_col = "Class" if "Class" in df.columns else "Drug class"
    genes = df[gene_col].astype(str).tolist() if gene_col in df.columns else []
    classes = df[class_col].astype(str).unique().tolist() if class_col in df.columns else []
    return {
        "amr_gene_count": len(df),
        "amr_classes": "; ".join(sorted(set(classes))[:30]),
        "amr_genes": "; ".join(sorted(set(genes))[:40]),
    }


def run_for_genome_dir(
    job: dict,
    out_root: Path,
    *,
    protein_preferred: bool = True,
) -> dict:
    sp = str(job.get("species", ""))
    acc = str(job.get("genome_accession", ""))
    gdir = Path(job["genome_dir"])
    row = {c: "" for c in AMR_SUMMARY_COLUMNS}
    row.update(
        {
            "BacID": job.get("BacID", ""),
            "species": sp,
            "strain": job.get("strain", ""),
            "genome_accession": acc,
        }
    )

    out_dir = out_root / species_slug(sp) / acc
    tsv_path = out_dir / "amrfinder.tsv"
    row["amr_output_path"] = str(tsv_path)

    prot = find_protein_fasta(gdir)
    nucl = find_nucleotide_fasta(gdir)

    if protein_preferred and prot is not None:
        ok, msg = run_amrfinder(prot, tsv_path, protein=True)
        row["amrfinder_notes"] = f"protein input: {prot.name}"
    elif nucl is not None:
        ok, msg = run_amrfinder(nucl, tsv_path, protein=False)
        row["amrfinder_notes"] = f"nucleotide input: {nucl.name}"
    else:
        row["amrfinder_status"] = "skipped"
        row["amrfinder_notes"] = "no FASTA in genome directory"
        return row

    row["amrfinder_status"] = "success" if ok else "failed"
    if not ok:
        row["amrfinder_notes"] = msg
        return row

    row.update(parse_amrfinder_tsv(tsv_path))
    return row


def main(argv: Optional[List[str]] = None) -> int:
    from genome_investigation.antismash_runner import collect_genome_jobs

    parser = argparse.ArgumentParser(description="AMRFinderPlus on downloaded genomes")
    parser.add_argument("--genome-enriched", required=True)
    parser.add_argument("--genome-dir", default=str(DEFAULT_GENOME_DIR))
    parser.add_argument("--selected-species", default=str(SELECTED_YAML))
    parser.add_argument("--species", action="append", default=None, help="Limit to species name(s)")
    parser.add_argument("--output", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--amrfinder-output-dir", default=str(Path(__file__).resolve().parents[1] / "results" / "amrfinder"))
    args = parser.parse_args(argv)

    if not amrfinder_installed():
        print("[error] amrfinder not on PATH")
        return 1

    enriched = pd.read_csv(args.genome_enriched)
    selected = load_selected_species(Path(args.selected_species))
    jobs = collect_genome_jobs(Path(args.genome_dir), enriched, selected)

    if args.species:
        want = {s.strip().lower() for s in args.species}
        jobs = [j for j in jobs if str(j.get("species", "")).strip().lower() in want]

    if not jobs:
        print("[warn] no genome jobs found")
        return 0

    rows = [run_for_genome_dir(j, Path(args.amrfinder_output_dir)) for j in jobs]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=AMR_SUMMARY_COLUMNS).to_csv(out, index=False)
    print(f"[done] wrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
