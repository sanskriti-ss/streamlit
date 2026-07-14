"""
Phase 1c: Optional antiSMASH integration (disabled by default).

Usage:
  python -m genome_investigation.antismash_runner --parse-antismash-only \\
    --antismash-output-dir results/antismash
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from genome_investigation.io_utils import load_selected_species, load_yaml, species_slug

DEFAULT_ANTISMASH_DIR = Path(__file__).resolve().parent.parent / "results" / "antismash"
DEFAULT_SUMMARY = Path(__file__).resolve().parent / "results" / "antismash_summary.csv"
CONFIG_PATH = Path(__file__).resolve().parent / "config" / "genome_config.yaml"
SELECTED_YAML = Path(__file__).resolve().parent / "selected_species.yaml"
DEFAULT_GENOME_DIR = Path(__file__).resolve().parent.parent / "data" / "genomes"

SUMMARY_COLUMNS = [
    "BacID",
    "species",
    "strain",
    "genome_accession",
    "antismash_version",
    "bgc_count_total",
    "bgc_types",
    "nrps_count",
    "pks_count",
    "terpene_count",
    "ribosomal_peptide_count",
    "saccharide_count",
    "siderophore_related_count",
    "other_bgc_count",
    "knownclusterblast_hits",
    "most_similar_known_clusters",
    "antismash_output_dir",
    "antismash_status",
    "antismash_notes",
]

BGC_TYPE_MAP = {
    "nrps": "nrps_count",
    "pks": "pks_count",
    "terpene": "terpene_count",
    "ripp": "ribosomal_peptide_count",
    "ribosomal": "ribosomal_peptide_count",
    "lassopeptide": "ribosomal_peptide_count",
    "proteusin": "ribosomal_peptide_count",
    "rre-containing": "ribosomal_peptide_count",
    "cyclic-lactone": "ribosomal_peptide_count",
    "saccharide": "saccharide_count",
    "siderophore": "siderophore_related_count",
    "ectoine": "siderophore_related_count",
    "betalactone": "other_bgc_count",
}


CONDA_ANTISMASH = Path.home() / "miniconda3" / "envs" / "antismash_env" / "bin" / "antismash"


def _antismash_binary() -> Optional[str]:
    """Return path to a working antismash binary or None."""
    for candidate in [shutil.which("antismash"), str(CONDA_ANTISMASH)]:
        if candidate and Path(candidate).exists():
            try:
                subprocess.run(
                    [candidate, "--version"],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                return candidate
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
                continue
    return None


def antismash_installed() -> bool:
    return _antismash_binary() is not None


def docker_installed() -> bool:
    if shutil.which("docker") is None:
        return False
    try:
        subprocess.run(["docker", "info"], check=True, capture_output=True, timeout=30)
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        return False


def antismash_runnable(*, allow_docker: bool = True) -> bool:
    if antismash_installed():
        return True
    return allow_docker and docker_installed()


def load_config(path: Path = CONFIG_PATH) -> dict:
    if not path.exists():
        return {"enable_antismash": False}
    return load_yaml(path)


def find_input_genbank(genome_dir: Path) -> tuple[Optional[Path], str]:
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


def parse_antismash_json(json_path: Path) -> dict:
    with json_path.open(encoding="utf-8") as f:
        data = json.load(f)

    version = ""
    if isinstance(data, dict):
        version = str(data.get("version") or data.get("antismash_version") or "")

    counts = {v: 0 for v in BGC_TYPE_MAP.values()}
    counts["other_bgc_count"] = 0
    types_seen: List[str] = []
    kcb_hits: List[str] = []
    similar_clusters: List[str] = []

    records = data.get("records") if isinstance(data, dict) else None
    if not isinstance(records, list):
        records = [data] if isinstance(data, dict) else []

    total = 0
    for rec in records:
        if not isinstance(rec, dict):
            continue
        areas = rec.get("areas") or rec.get("region_predictions") or []
        if isinstance(areas, dict):
            areas = list(areas.values())
        for area in areas:
            if not isinstance(area, dict):
                continue
            products = area.get("products") or []
            if isinstance(products, dict):
                products = list(products.values())
            for prod in products:
                # antiSMASH v8: products is a list of strings
                if isinstance(prod, str):
                    pclass = prod.lower()
                    total += 1
                    types_seen.append(pclass)
                    mapped = False
                    for key, col in BGC_TYPE_MAP.items():
                        if key in pclass:
                            counts[col] += 1
                            mapped = True
                            break
                    if not mapped:
                        counts["other_bgc_count"] += 1
                    continue
                # older format: products is a list of dicts
                if not isinstance(prod, dict):
                    continue
                total += 1
                pclass = str(prod.get("product_class") or prod.get("category") or "other").lower()
                types_seen.append(pclass)
                mapped = False
                for key, col in BGC_TYPE_MAP.items():
                    if key in pclass:
                        counts[col] += 1
                        mapped = True
                        break
                if not mapped:
                    counts["other_bgc_count"] += 1
            kcb = area.get("knownclusterblast") or area.get("knowncluster")
            if isinstance(kcb, dict):
                for cl_name, cl_data in kcb.items():
                    kcb_hits.append(str(cl_name))
                    if isinstance(cl_data, dict) and cl_data.get("similarity"):
                        similar_clusters.append(f"{cl_name}:{cl_data.get('similarity')}")

    return {
        "antismash_version": version,
        "bgc_count_total": total,
        "bgc_types": "; ".join(sorted(set(types_seen))),
        "knownclusterblast_hits": "; ".join(kcb_hits[:20]),
        "most_similar_known_clusters": "; ".join(similar_clusters[:20]),
        **counts,
    }


def parse_output_directory(out_dir: Path) -> dict:
    json_files = list(out_dir.glob("*.json"))
    # antiSMASH v8 produces genomic.json as the main summary
    summary = next((p for p in json_files if p.name == "genomic.json"), None)
    if summary is None:
        # fall back to any .json at top level, then recurse
        summary = json_files[0] if json_files else next(iter(out_dir.rglob("*.json")), None)
    if summary is None:
        return {"bgc_count_total": 0, "antismash_notes": "no JSON found in output directory"}
    parsed = parse_antismash_json(summary)
    parsed["antismash_notes"] = f"parsed {summary.name}"
    return parsed


def run_antismash(
    input_path: Path,
    out_dir: Path,
    *,
    input_type: str,
    extra_args: Optional[List[str]] = None,
) -> tuple[bool, str]:
    if not antismash_installed():
        return False, "antiSMASH executable not found (optional dependency)"

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["antismash", str(input_path), "--output-dir", str(out_dir)]
    extra = list(extra_args or [])
    if input_type == "genbank":
        cmd.extend(["--genefinding-tool", "none"])
    elif input_type == "fasta" and "--taxon" not in extra:
        cmd.extend(["--taxon", "bacteria"])
    if extra:
        cmd.extend(extra)

    binary = _antismash_binary() or "antismash"
    cmd[0] = binary
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)
        return True, "completed"
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        err = getattr(exc, "stderr", None) or str(exc)
        return False, f"antiSMASH failed: {err}"


def run_antismash_docker(
    input_path: Path,
    out_dir: Path,
    *,
    input_type: str,
    image: str = "antismash/standalone:latest",
    extra_args: Optional[List[str]] = None,
) -> tuple[bool, str]:
    """Run antiSMASH via Docker when the local conda binary is unavailable."""
    if not docker_installed():
        return False, "docker not available"

    out_dir.mkdir(parents=True, exist_ok=True)
    in_dir = input_path.resolve().parent
    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform", "linux/amd64",
        "-v",
        f"{in_dir}:/input:ro",
        "-v",
        f"{out_dir.resolve()}:/output",
        image,
        "antismash",
        f"/input/{input_path.name}",
        "--output-dir",
        "/output",
        "--genefinding-tool", "none",
    ]
    if input_type == "fasta":
        cmd.extend(["--taxon", "bacteria"])
    if extra_args:
        cmd.extend(extra_args)

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)
        return True, "completed (docker)"
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        err = getattr(exc, "stderr", None) or str(exc)
        return False, f"antiSMASH docker failed: {err}"


def run_antismash_auto(
    input_path: Path,
    out_dir: Path,
    *,
    input_type: str,
    use_docker: bool = True,
    docker_image: str = "antismash/standalone:latest",
    extra_args: Optional[List[str]] = None,
) -> tuple[bool, str]:
    if antismash_installed():
        ok, msg = run_antismash(input_path, out_dir, input_type=input_type, extra_args=extra_args)
        if ok:
            return ok, msg
    if use_docker and docker_installed():
        return run_antismash_docker(
            input_path,
            out_dir,
            input_type=input_type,
            image=docker_image,
            extra_args=extra_args,
        )
    return False, "antiSMASH not available (install locally or use Docker)"


def collect_genome_jobs(
    genome_root: Path,
    enriched: pd.DataFrame,
    selected: dict,
) -> List[dict]:
    species_set = {s.strip().lower() for s in selected.get("species", [])}
    jobs: List[dict] = []
    for _, row in enriched.iterrows():
        sp = str(row.get("species", "")).strip()
        if species_set and sp.lower() not in species_set:
            continue
        acc = str(row.get("genome_accession", "")).strip()
        if not acc:
            continue
        gdir = genome_root / species_slug(sp) / acc
        if not gdir.exists():
            continue
        jobs.append(row.to_dict() | {"genome_dir": str(gdir)})
    return jobs


def build_summary_rows(
    jobs: List[dict],
    antismash_root: Path,
    *,
    run: bool,
    parse_only: bool,
    use_docker: bool = False,
    docker_image: str = "antismash/standalone:latest",
    extra_args: Optional[List[str]] = None,
) -> List[dict]:
    rows: List[dict] = []
    cfg = load_config()
    if not cfg.get("enable_antismash", False) and run and not parse_only:
        print("[info] enable_antismash is false in config; use --run-antismash to override for this run")

    for job in jobs:
        sp = job.get("species", "")
        acc = job.get("genome_accession", "")
        gdir = Path(job["genome_dir"])
        out_dir = antismash_root / species_slug(sp) / acc
        row = {c: "" for c in SUMMARY_COLUMNS}
        row.update(
            {
                "BacID": job.get("BacID", ""),
                "species": sp,
                "strain": job.get("strain", ""),
                "genome_accession": acc,
                "antismash_output_dir": str(out_dir),
            }
        )

        inp, itype = find_input_genbank(gdir)
        if inp is None:
            row["antismash_status"] = "skipped"
            row["antismash_notes"] = "no genome file in download directory"
            rows.append(row)
            continue

        if run and not parse_only:
            if itype == "fasta":
                row["antismash_notes"] = "FASTA input; annotation quality may be lower than GenBank/GBFF"
            ok, msg = run_antismash_auto(
                inp,
                out_dir,
                input_type=itype,
                use_docker=use_docker,
                docker_image=docker_image,
                extra_args=extra_args,
            )
            row["antismash_status"] = "success" if ok else "failed"
            if not ok:
                row["antismash_notes"] = msg
                rows.append(row)
                continue

        # Also check parent dir (species-only, without accession subfolder)
        effective_dir = out_dir
        if not effective_dir.exists():
            parent_dir = antismash_root / species_slug(sp)
            if parent_dir.exists() and any(parent_dir.glob("*.json")):
                effective_dir = parent_dir

        if effective_dir.exists():
            parsed = parse_output_directory(effective_dir)
            row["antismash_output_dir"] = str(effective_dir)
            row.update(parsed)
            row["antismash_status"] = row.get("antismash_status") or "parsed"
        else:
            row["antismash_status"] = "missing_output"
            row["antismash_notes"] = "output directory not found"

        rows.append(row)
    return rows


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Phase 1c: optional antiSMASH runner")
    parser.add_argument("--selected-species", default=str(SELECTED_YAML))
    parser.add_argument("--genome-enriched", default=str(Path(__file__).resolve().parent / "results" / "Step2_5_genome_enriched.csv"))
    parser.add_argument("--genome-dir", default=str(DEFAULT_GENOME_DIR))
    parser.add_argument("--run-antismash", action="store_true")
    parser.add_argument("--parse-antismash-only", action="store_true")
    parser.add_argument("--antismash-output-dir", default=str(DEFAULT_ANTISMASH_DIR))
    parser.add_argument("--summary-output", default=str(DEFAULT_SUMMARY))
    parser.add_argument("--use-docker", action="store_true", help="Fall back to Docker when local antiSMASH is unavailable")
    parser.add_argument("--antismash-extra-args", nargs=argparse.REMAINDER, default=None)
    args = parser.parse_args(argv)

    if args.run_antismash and not antismash_runnable(allow_docker=args.use_docker):
        print(
            "[warn] antiSMASH is not runnable in this environment.\n"
            "       Install from https://docs.antismash.secondarymetabolites.org/latest/install/\n"
            "       or pass --use-docker if Docker is available,\n"
            "       or run with --parse-antismash-only to import existing results."
        )
        if not args.parse_antismash_only:
            return 1

    selected = load_selected_species(Path(args.selected_species))
    enriched_path = Path(args.genome_enriched)
    if not enriched_path.exists():
        print(f"[error] missing {enriched_path}")
        return 1

    enriched = pd.read_csv(enriched_path)
    jobs = collect_genome_jobs(Path(args.genome_dir), enriched, selected)
    if not jobs:
        print("[warn] no downloaded genomes found for selected species")
        return 0

    rows = build_summary_rows(
        jobs,
        Path(args.antismash_output_dir),
        run=args.run_antismash,
        parse_only=args.parse_antismash_only,
        use_docker=args.use_docker,
        extra_args=args.antismash_extra_args,
    )
    out = Path(args.summary_output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=SUMMARY_COLUMNS).to_csv(out, index=False)
    print(f"[done] wrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
