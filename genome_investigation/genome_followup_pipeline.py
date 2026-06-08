"""
End-to-end genomic follow-up for focal BacDive species.

Steps:
  1. Build / refresh focal genome-enriched table
  2. Download assemblies (NCBI Datasets)
  3. antiSMASH (high-confidence strains)
  4. AMRFinderPlus (resistance focal species)
  5. BacDive metabolite phenotype export (API cache)
  6. Targeted gene search (protein/GFF + AMR)
  7. Summary tables + figures

Usage:
  python -m genome_investigation.genome_followup_pipeline
  python -m genome_investigation.genome_followup_pipeline --skip-download --skip-antismash
"""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from genome_investigation.amrfinder_runner import (
    AMR_SUMMARY_COLUMNS,
    amrfinder_installed,
    run_for_genome_dir,
)
from genome_investigation.antismash_runner import (
    SUMMARY_COLUMNS as ANTISMASH_COLUMNS,
    antismash_runnable,
    build_summary_rows,
    collect_genome_jobs,
    docker_installed,
)
from genome_investigation.bacdive_phenotype_extract import build_phenotype_table
from genome_investigation.phenotype_confidence import build_confidence_table
from genome_investigation.genome_download import (
    download_genome_package,
    filter_selected_rows,
)
from genome_investigation.genome_enrichment import enrich_row, OUTPUT_COLUMNS
from genome_investigation.io_utils import load_selected_species, normalize_bacid, species_slug
from genome_investigation.targeted_gene_search import run_targeted_search
from genome_investigation.visualize_genomic_followup import generate_followup_figures

PKG_DIR = Path(__file__).resolve().parent
REPO_ROOT = PKG_DIR.parent
RESULTS_DIR = PKG_DIR / "results"
FOLLOWUP_DIR = RESULTS_DIR / "genomic_followup"
SELECTED_YAML = PKG_DIR / "selected_species.yaml"
DEFAULT_INPUT = RESULTS_DIR / "selected_test_input.csv"
DEFAULT_ENRICHED_ALL = RESULTS_DIR / "Step2_5_genome_enriched.csv"
FOCAL_ENRICHED = FOLLOWUP_DIR / "focal_genome_enriched.csv"

DEFAULT_GENOME_DIR = REPO_ROOT / "data" / "genomes"
DEFAULT_ANTISMASH_DIR = REPO_ROOT / "results" / "antismash"
DEFAULT_AMRFINDER_DIR = REPO_ROOT / "results" / "amrfinder"

AMR_FOCAL_SPECIES = ["Bacillus cereus", "Micrococcus luteus"]
# Primary BGC / production candidates (avoid running antiSMASH on all nine assemblies)
ANTISMASH_FOCAL_SPECIES = [
    "Cohnella algarum",
    "Clostridium cellulovorans",
    "Acetomicrobium thermoterrenum",
]


def _stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def build_focal_enriched(
    input_csv: Path,
    *,
    cache_dir: Path,
    allow_species_only: bool = True,
    existing_enriched: Optional[Path] = None,
) -> pd.DataFrame:
    """One row per input strain; enrich via BacDive/NCBI."""
    inp = pd.read_csv(input_csv)
    inp["BacID"] = inp["BacID"].map(normalize_bacid)

    rows: List[dict] = []
    for _, r in inp.iterrows():
        base = {c: r.get(c, "") for c in OUTPUT_COLUMNS if c in r.index}
        for c in OUTPUT_COLUMNS:
            base.setdefault(c, "")
        enriched = enrich_row(
            base,
            cache_dir=cache_dir,
            force_refresh=False,
            allow_species_only=allow_species_only,
            rate_limit_s=0.2,
        )
        rows.append(enriched)

    focal = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    # Prefer higher confidence if merging with catalog
    if existing_enriched and existing_enriched.exists():
        catalog = pd.read_csv(existing_enriched)
        catalog["BacID_norm"] = catalog["BacID"].map(normalize_bacid)
        focal["BacID_norm"] = focal["BacID"].map(normalize_bacid)
        for i, row in focal.iterrows():
            bid = row["BacID_norm"]
            sp = str(row.get("species", ""))
            alt = catalog[(catalog["BacID_norm"] == bid) | (catalog["species"] == sp)]
            if not alt.empty:
                best = alt.sort_values("match_confidence", ascending=False).iloc[0]
                if float(best["match_confidence"]) > float(row.get("match_confidence") or 0):
                    for c in OUTPUT_COLUMNS:
                        if pd.notna(best.get(c)) and str(best.get(c)) not in ("", "nan"):
                            focal.at[i, c] = best[c]
        focal = focal.drop(columns=["BacID_norm"], errors="ignore")

    return focal


def write_manifest(records: List[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(path, index=False)


def download_targets(
    targets: pd.DataFrame,
    genome_dir: Path,
) -> List[dict]:
    manifest: List[dict] = []
    for _, r in targets.iterrows():
        acc = str(r["genome_accession"]).strip()
        sp = str(r["species"])
        slug = species_slug(sp)
        dest = genome_dir / slug / acc
        t0 = time.time()
        ok, msg = download_genome_package(acc, dest)
        manifest.append(
            {
                "species": sp,
                "BacID": r.get("BacID", ""),
                "genome_accession": acc,
                "match_confidence": r.get("match_confidence"),
                "download_status": "success" if ok else "failed",
                "download_notes": msg,
                "genome_dir": str(dest),
                "elapsed_s": round(time.time() - t0, 1),
            }
        )
        status = "ok" if ok else "warn"
        print(f"[{status}] download {sp} ({acc}): {msg}")
    return manifest


def run_pipeline(
    *,
    input_csv: Path = DEFAULT_INPUT,
    selected_yaml: Path = SELECTED_YAML,
    out_dir: Path = FOLLOWUP_DIR,
    genome_dir: Path = DEFAULT_GENOME_DIR,
    antismash_dir: Path = DEFAULT_ANTISMASH_DIR,
    amrfinder_dir: Path = DEFAULT_AMRFINDER_DIR,
    min_confidence_download: float = 0.7,
    antismash_min_confidence: float = 1.0,
    skip_enrichment: bool = False,
    skip_download: bool = False,
    skip_antismash: bool = False,
    skip_amrfinder: bool = False,
    skip_phenotype: bool = False,
    skip_gene_search: bool = False,
    skip_figures: bool = False,
    run_antismash: bool = True,
    run_amrfinder: bool = True,
    antismash_docker: bool = True,
    force_refresh_phenotype: bool = False,
) -> dict:
    from genome_investigation.api_cache import DEFAULT_CACHE_DIR

    out_dir.mkdir(parents=True, exist_ok=True)
    run_id = _stamp()
    run_dir = out_dir / f"run_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)

    latest_link = {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "started": datetime.now().isoformat(),
        "steps": {},
    }

    selected = load_selected_species(selected_yaml)
    focal_path = run_dir / "focal_genome_enriched.csv"

    if skip_enrichment and FOCAL_ENRICHED.exists():
        focal = pd.read_csv(FOCAL_ENRICHED)
        focal.to_csv(focal_path, index=False)
        print(f"[info] using existing focal enrichment: {FOCAL_ENRICHED}")
    else:
        print("[step 1] enriching focal strains...")
        focal = build_focal_enriched(
            input_csv,
            cache_dir=DEFAULT_CACHE_DIR,
            existing_enriched=DEFAULT_ENRICHED_ALL,
        )
        focal.to_csv(focal_path, index=False)
        focal.to_csv(FOCAL_ENRICHED, index=False)
        print(f"[ok] wrote {focal_path}")

    latest_link["steps"]["enrichment"] = {"rows": len(focal), "path": str(focal_path)}

    # Download: all focal >= min_confidence_download
    dl_targets = filter_selected_rows(
        focal,
        selected,
        min_confidence=min_confidence_download,
        max_genomes=None,
    )
    dl_manifest_path = run_dir / "download_manifest.csv"
    if skip_download:
        print("[step 2] download skipped")
        if dl_manifest_path.exists():
            manifest = pd.read_csv(dl_manifest_path).to_dict("records")
        else:
            manifest = []
    else:
        print(f"[step 2] downloading {len(dl_targets)} genomes (confidence >= {min_confidence_download})...")
        manifest = download_targets(dl_targets, genome_dir)
        write_manifest(manifest, dl_manifest_path)
    latest_link["steps"]["download"] = {"n": len(manifest), "path": str(dl_manifest_path)}

    jobs = collect_genome_jobs(genome_dir, focal, selected)
    if not jobs:
        print("[warn] no downloaded genomes found — later steps may be empty")

    # antiSMASH
    asm_summary_path = run_dir / "antismash_summary.csv"
    high_conf = focal[pd.to_numeric(focal["match_confidence"], errors="coerce") >= antismash_min_confidence]
    asm_species = set(ANTISMASH_FOCAL_SPECIES) & set(high_conf["species"].astype(str))
    asm_jobs = [j for j in jobs if str(j.get("species")) in asm_species]

    if skip_antismash:
        print("[step 3] antiSMASH skipped")
    elif not asm_jobs:
        print("[step 3] no focal genomes on disk for antiSMASH")
    elif run_antismash and not antismash_runnable(allow_docker=antismash_docker):
        print("[step 3] antiSMASH not available (install or enable Docker) — skipping run")
        skip_antismash = True

    if not skip_antismash:
        do_run = run_antismash and antismash_runnable(allow_docker=antismash_docker)
        mode = "run"
        if do_run and not antismash_runnable(allow_docker=False) and antismash_docker and docker_installed():
            mode = "docker"
        print(f"[step 3] antiSMASH ({mode}) on {len(asm_jobs)} strain(s)...")
        asm_rows = build_summary_rows(
            asm_jobs,
            antismash_dir,
            run=do_run,
            parse_only=not do_run,
            use_docker=antismash_docker,
        )
        pd.DataFrame(asm_rows, columns=ANTISMASH_COLUMNS).to_csv(asm_summary_path, index=False)
        pd.DataFrame(asm_rows, columns=ANTISMASH_COLUMNS).to_csv(RESULTS_DIR / "antismash_summary.csv", index=False)
        print(f"[ok] antiSMASH summary → {asm_summary_path}")
    else:
        asm_summary_path = RESULTS_DIR / "antismash_summary.csv" if (RESULTS_DIR / "antismash_summary.csv").exists() else None

    latest_link["steps"]["antismash"] = {"path": str(asm_summary_path) if asm_summary_path else None}

    # AMRFinder
    amr_summary_path = run_dir / "amrfinder_summary.csv"
    amr_jobs = [j for j in jobs if str(j.get("species", "")) in AMR_FOCAL_SPECIES]

    if skip_amrfinder:
        print("[step 4] AMRFinder skipped")
    elif not amr_jobs:
        print("[step 4] resistance species genomes not downloaded — skipping AMRFinder")
    elif run_amrfinder and not amrfinder_installed():
        print("[step 4] amrfinder not on PATH — skipping")
        skip_amrfinder = True

    if not skip_amrfinder and amr_jobs:
        print(f"[step 4] AMRFinder on {len(amr_jobs)} genome(s)...")
        amr_rows = [run_for_genome_dir(j, amrfinder_dir) for j in amr_jobs]
        pd.DataFrame(amr_rows, columns=AMR_SUMMARY_COLUMNS).to_csv(amr_summary_path, index=False)
        pd.DataFrame(amr_rows, columns=AMR_SUMMARY_COLUMNS).to_csv(RESULTS_DIR / "amrfinder_summary.csv", index=False)
        print(f"[ok] AMRFinder summary → {amr_summary_path}")
    else:
        amr_summary_path = RESULTS_DIR / "amrfinder_summary.csv" if (RESULTS_DIR / "amrfinder_summary.csv").exists() else None

    latest_link["steps"]["amrfinder"] = {"path": str(amr_summary_path) if amr_summary_path else None}

    # BacDive phenotypes
    pheno_path = run_dir / "bacdive_phenotype_metabolites.csv"
    if skip_phenotype:
        print("[step 5] BacDive phenotype export skipped")
    else:
        print("[step 5] extracting BacDive metabolite lists from API cache...")
        pheno = build_phenotype_table(
            focal[["BacID", "species", "strain"]],
            force_refresh=force_refresh_phenotype,
        )
        pheno.to_csv(pheno_path, index=False)
        pheno.to_csv(FOLLOWUP_DIR / "bacdive_phenotype_metabolites.csv", index=False)
        print(f"[ok] {len(pheno)} phenotype rows → {pheno_path}")
    latest_link["steps"]["phenotype"] = {"path": str(pheno_path)}

    # Targeted gene search
    hits_path = run_dir / "targeted_gene_hits.csv"
    if skip_gene_search:
        print("[step 6] targeted gene search skipped")
    elif not jobs:
        print("[step 6] no genomes — skipping gene search")
    else:
        print("[step 6] targeted gene / AMR hit search...")
        pheno_df = pd.read_csv(pheno_path) if pheno_path.exists() else pd.DataFrame()
        amr_df = pd.read_csv(amr_summary_path) if amr_summary_path and Path(amr_summary_path).exists() else pd.DataFrame()
        hit_frames: List[pd.DataFrame] = []
        for job in jobs:
            amr_row = None
            if not amr_df.empty:
                m = amr_df[amr_df["species"].astype(str) == str(job.get("species"))]
                if not m.empty:
                    amr_row = m.iloc[0].to_dict()
            hit_frames.append(run_targeted_search(job, pheno_df, amr_row))
        hits = pd.concat([h for h in hit_frames if not h.empty], ignore_index=True) if hit_frames else pd.DataFrame()
        hits.to_csv(hits_path, index=False)
        hits.to_csv(FOLLOWUP_DIR / "targeted_gene_hits.csv", index=False)
        print(f"[ok] {len(hits)} gene hits → {hits_path}")
    latest_link["steps"]["gene_hits"] = {"path": str(hits_path)}

    # Per-metabolite phenotype confidence (resistance / production / utilization)
    conf_path = run_dir / "phenotype_confidence.csv"
    print("[step 6b] building per-metabolite phenotype confidence...")
    pheno_df = pd.read_csv(pheno_path) if pheno_path.exists() else pd.DataFrame()
    asm_df = pd.read_csv(asm_summary_path) if asm_summary_path and Path(asm_summary_path).exists() else pd.DataFrame()
    conf = build_confidence_table(
        focal,
        pheno_df,
        amr_root=amrfinder_dir,
        genome_root=genome_dir,
        antismash_summary=asm_df,
    )
    conf.to_csv(conf_path, index=False)
    conf.to_csv(FOLLOWUP_DIR / "phenotype_confidence.csv", index=False)
    print(f"[ok] {len(conf)} confidence rows → {conf_path}")
    latest_link["steps"]["confidence"] = {"path": str(conf_path)}

    # Integrated summary
    integrated_path = run_dir / "integrated_genomic_followup.csv"
    integrated = _build_integrated_summary(focal, run_dir)
    integrated.to_csv(integrated_path, index=False)
    integrated.to_csv(FOLLOWUP_DIR / "integrated_genomic_followup.csv", index=False)
    latest_link["steps"]["integrated"] = {"path": str(integrated_path)}

    # Figures
    fig_dir = run_dir / "figures"
    if not skip_figures:
        print("[step 7] generating figures...")
        paths = generate_followup_figures(
            focal,
            run_dir,
            fig_dir,
            phenotype_path=pheno_path if pheno_path.exists() else None,
            hits_path=hits_path if hits_path.exists() else None,
        )
        latest_link["steps"]["figures"] = paths
        print(f"[ok] figures in {fig_dir}")

    latest_link["finished"] = datetime.now().isoformat()
    meta_path = run_dir / "pipeline_run.json"
    meta_path.write_text(json.dumps(latest_link, indent=2), encoding="utf-8")
    (FOLLOWUP_DIR / "latest_run.json").write_text(json.dumps(latest_link, indent=2), encoding="utf-8")
    print(f"\n[done] pipeline run {run_id} → {run_dir}")
    return latest_link


def _build_integrated_summary(focal: pd.DataFrame, run_dir: Path) -> pd.DataFrame:
    rows = []
    asm = pd.read_csv(run_dir / "antismash_summary.csv") if (run_dir / "antismash_summary.csv").exists() else pd.DataFrame()
    amr = pd.read_csv(run_dir / "amrfinder_summary.csv") if (run_dir / "amrfinder_summary.csv").exists() else pd.DataFrame()
    pheno = pd.read_csv(run_dir / "bacdive_phenotype_metabolites.csv") if (run_dir / "bacdive_phenotype_metabolites.csv").exists() else pd.DataFrame()
    hits = pd.read_csv(run_dir / "targeted_gene_hits.csv") if (run_dir / "targeted_gene_hits.csv").exists() else pd.DataFrame()

    for _, r in focal.iterrows():
        sp = str(r["species"])
        row = r.to_dict()
        if not asm.empty:
            m = asm[asm["species"].astype(str) == sp]
            if not m.empty:
                a = m.iloc[0]
                row["bgc_count_total"] = a.get("bgc_count_total", 0)
                row["bgc_types"] = a.get("bgc_types", "")
                row["antismash_status"] = a.get("antismash_status", "")
        if not amr.empty:
            m = amr[amr["species"].astype(str) == sp]
            if not m.empty:
                a = m.iloc[0]
                row["amr_gene_count"] = a.get("amr_gene_count", 0)
                row["amr_genes"] = a.get("amr_genes", "")
                row["amrfinder_status"] = a.get("amrfinder_status", "")
        if not pheno.empty:
            p = pheno[pheno["species"].astype(str) == sp]
            row["bacdive_util_n"] = int((p["activity"] == "utilization").sum())
            row["bacdive_prod_n"] = int((p["activity"] == "production").sum())
            row["bacdive_resistance_n"] = int((p["activity"] == "resistance").sum())
        if not hits.empty:
            h = hits[hits["species"].astype(str) == sp]
            row["genomic_hit_n"] = len(h)
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Automated genome follow-up pipeline")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--selected-species", default=str(SELECTED_YAML))
    parser.add_argument("--output-dir", default=str(FOLLOWUP_DIR))
    parser.add_argument("--genome-dir", default=str(DEFAULT_GENOME_DIR))
    parser.add_argument("--min-confidence-download", type=float, default=0.7)
    parser.add_argument("--antismash-min-confidence", type=float, default=1.0)
    parser.add_argument("--skip-enrichment", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-antismash", action="store_true")
    parser.add_argument("--skip-amrfinder", action="store_true")
    parser.add_argument("--skip-phenotype", action="store_true")
    parser.add_argument("--skip-gene-search", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--no-run-antismash", action="store_true", help="Parse existing antiSMASH output only")
    parser.add_argument("--no-run-amrfinder", action="store_true", help="Skip AMRFinder execution")
    parser.add_argument("--no-antismash-docker", action="store_true", help="Do not fall back to Docker for antiSMASH")
    parser.add_argument("--force-phenotype-refresh", action="store_true")
    args = parser.parse_args(argv)

    run_pipeline(
        input_csv=Path(args.input),
        selected_yaml=Path(args.selected_species),
        out_dir=Path(args.output_dir),
        genome_dir=Path(args.genome_dir),
        min_confidence_download=args.min_confidence_download,
        antismash_min_confidence=args.antismash_min_confidence,
        skip_enrichment=args.skip_enrichment,
        skip_download=args.skip_download,
        skip_antismash=args.skip_antismash,
        skip_amrfinder=args.skip_amrfinder,
        skip_phenotype=args.skip_phenotype,
        skip_gene_search=args.skip_gene_search,
        skip_figures=args.skip_figures,
        run_antismash=not args.no_run_antismash,
        run_amrfinder=not args.no_run_amrfinder,
        antismash_docker=not args.no_antismash_docker,
        force_refresh_phenotype=args.force_phenotype_refresh,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
