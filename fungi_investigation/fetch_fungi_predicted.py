"""
Orchestrate predicted fungi layers: genome download → antiSMASH → utilization scan.

Skips FUNG-GROWTH experimental ingest (manual website export).

Usage:
  python -m fungi_investigation.fetch_fungi_predicted
  python -m fungi_investigation.fetch_fungi_predicted --skip-download --skip-antismash
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

ANTISMASH_ENV_BIN = Path.home() / "miniconda3" / "envs" / "antismash_env" / "bin"

from genome_investigation.antismash_runner import (
    antismash_installed,
    find_input_genbank,
    run_antismash,
)
from genome_investigation.io_utils import load_yaml, species_slug


def _find_genome_fasta(genome_dir: Path) -> Optional[Path]:
    for pattern in ("*.fna", "*.fasta", "*.fa"):
        hits = sorted(genome_dir.rglob(pattern))
        # Prefer genomic.fna over misc assembly files when present.
        preferred = [p for p in hits if "genomic" in p.name.lower()]
        if preferred:
            return preferred[0]
        if hits:
            return hits[0]
    return None


def _should_retry_with_fasta(message: str) -> bool:
    blob = (message or "").lower()
    markers = (
        "multiple cds features have the same location",
        "duplicate cds",
        "invalid location",
        "could not parse",
    )
    return any(m in blob for m in markers)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_LIST = Path(__file__).resolve().parent / "selected_fungi.yaml"
DEFAULT_GENOMES = REPO_ROOT / "data" / "genomes_fungi"
DEFAULT_ANTISMASH = REPO_ROOT / "results" / "antismash_fungi"
PROD_OUT = REPO_ROOT / "fungi_data" / "predicted" / "fungi_phenotype_confidence.csv"
UTIL_OUT = REPO_ROOT / "fungi_data" / "predicted" / "fungi_utilization_confidence.csv"


def _run(cmd: List[str]) -> int:
    print(f"[cmd] {' '.join(cmd)}")
    return subprocess.call(cmd)


def _load_targets(path: Path) -> List[dict]:
    data = load_yaml(path)
    rows = data.get("species", data) if isinstance(data, dict) else data
    return [r if isinstance(r, dict) else {"species": r, "accession": ""} for r in rows]


def run_antismash_on_fungi(
    targets: List[dict],
    genome_root: Path,
    antismash_root: Path,
) -> int:
    if not antismash_installed():
        print("[error] antiSMASH not installed; cannot build production predictions")
        return 1

    if ANTISMASH_ENV_BIN.is_dir():
        os.environ["PATH"] = f"{ANTISMASH_ENV_BIN}{os.pathsep}{os.environ.get('PATH', '')}"

    ok, fail = 0, 0
    for row in targets:
        species = str(row["species"]).strip()
        acc = str(row.get("accession", "")).strip()
        slug = species_slug(species)
        gdir = genome_root / slug / acc if acc else genome_root / slug
        if not gdir.exists():
            candidates = sorted((genome_root / slug).glob("*")) if (genome_root / slug).exists() else []
            gdir = candidates[0] if candidates else gdir
        if not gdir.exists():
            print(f"[warn] no genome dir for {species}")
            fail += 1
            continue

        inp, itype = find_input_genbank(gdir)
        if inp is None:
            print(f"[warn] no GBFF/FASTA for {species} in {gdir}")
            fail += 1
            continue

        out_dir = antismash_root / slug / (acc or inp.parent.name)
        if (out_dir / "genomic.json").exists():
            print(f"[skip] antiSMASH output exists for {species}")
            ok += 1
            continue

        # Fungal antiSMASH always declares taxon=fungi (extra_args override defaults
        # for FASTA; GBFF uses embedded annotations with --genefinding-tool none).
        extra = ["--taxon", "fungi"]
        success, msg = run_antismash(inp, out_dir, input_type=itype, extra_args=extra)
        if (not success) and itype == "genbank" and _should_retry_with_fasta(msg):
            fasta = _find_genome_fasta(gdir)
            if fasta is not None:
                print(
                    f"[retry] {species}: GBFF failed ({msg.splitlines()[0][:120]}); "
                    f"retrying with FASTA {fasta.name}"
                )
                # Clear partial failed output so antiSMASH can rewrite the directory.
                for leftover in out_dir.glob("*"):
                    if leftover.is_file():
                        leftover.unlink()
                success, msg = run_antismash(
                    fasta, out_dir, input_type="fasta", extra_args=extra
                )
            else:
                print(
                    f"[warn] {species}: GBFF failed and no genomic FASTA found "
                    "(re-run fungi_genome_download to fetch GENOME_FASTA)"
                )

        if success:
            ok += 1
            print(f"[ok] antiSMASH {species}: {msg}")
        else:
            fail += 1
            print(f"[warn] antiSMASH {species}: {msg}")

    print(f"[antiSMASH] completed {ok}, failed {fail}")
    # Partial success is OK — still ingest whatever genomic.json files exist.
    return 0 if ok > 0 else 2


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Build predicted fungi data layers")
    parser.add_argument("--species-list", type=Path, default=DEFAULT_LIST)
    parser.add_argument("--genome-dir", type=Path, default=DEFAULT_GENOMES)
    parser.add_argument("--antismash-dir", type=Path, default=DEFAULT_ANTISMASH)
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-antismash", action="store_true")
    parser.add_argument("--skip-utilization", action="store_true")
    args = parser.parse_args(argv)

    targets = _load_targets(args.species_list)

    if not args.skip_download:
        rc = _run(
            [
                sys.executable,
                "-m",
                "fungi_investigation.fungi_genome_download",
                "--species-list",
                str(args.species_list),
                "--output-dir",
                str(args.genome_dir),
            ]
        )
        if rc != 0:
            return rc

    if not args.skip_antismash:
        rc = run_antismash_on_fungi(targets, args.genome_dir, args.antismash_dir)
        if rc != 0:
            print(
                "[warn] antiSMASH finished with failures; continuing to ingest "
                "any successful outputs and build utilization layer"
            )
        ingest_rc = _run(
            [
                sys.executable,
                "-m",
                "fungi_investigation.antismash_db_ingest",
                "--json-root",
                str(args.antismash_dir),
                "--output",
                str(PROD_OUT),
            ]
        )
        if ingest_rc != 0:
            print("[warn] antismash_db_ingest failed; utilization will still run")

    if not args.skip_utilization:
        rc = _run(
            [
                sys.executable,
                "-m",
                "fungi_investigation.fungi_utilization_predict",
                "--fasta-root",
                str(args.genome_dir),
                "--output",
                str(UTIL_OUT),
            ]
        )
        if rc != 0:
            return rc

    print("[done] predicted fungi layers updated")
    print(f"  production:   {PROD_OUT}")
    print(f"  utilization:  {UTIL_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
