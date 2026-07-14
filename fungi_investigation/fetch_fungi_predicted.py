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

from fungi_investigation.gbff_dedupe_cds import write_deduped_gbff
from genome_investigation.antismash_runner import (
    antismash_installed,
    find_input_genbank,
    run_antismash,
)
from genome_investigation.io_utils import load_yaml, species_slug


def _should_retry_dedupe_gbff(message: str) -> bool:
    blob = (message or "").lower()
    markers = (
        "multiple cds features have the same location",
        "duplicate cds",
        "same name for mapping",
        "overlapping exons",
    )
    return any(m in blob for m in markers)


def _clear_outdir_files(out_dir: Path) -> None:
    if not out_dir.exists():
        return
    for leftover in out_dir.glob("*"):
        if leftover.is_file():
            leftover.unlink()


def _normalize_antismash_json(out_dir: Path) -> None:
    """Ensure ingest can find genomic.json when input was *.dedup.gbff."""
    canonical = out_dir / "genomic.json"
    if canonical.exists():
        return
    candidates = sorted(out_dir.glob("*.json"))
    # Prefer files that look like the primary antiSMASH result.
    preferred = [
        p
        for p in candidates
        if "dedup" in p.name.lower() or p.name.lower().startswith("genomic")
    ]
    src = (preferred or candidates)[0] if (preferred or candidates) else None
    if src is None:
        return
    src.rename(canonical)

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

        # Fungal antiSMASH requires annotated GenBank (no fungal gene finder for FASTA).
        # Proactively clean RefSeq GBFF quirks (dupe CDS loc/name, overlapping exons).
        extra = ["--taxon", "fungi"]
        run_input, run_type = inp, itype
        if itype == "genbank":
            dedup_path = inp.with_suffix(inp.suffix + ".dedup.gbff")
            try:
                nfix = write_deduped_gbff(inp, dedup_path)
                if nfix > 0:
                    print(f"[clean] {species}: applied {nfix} CDS fix(es) → {dedup_path.name}")
                    run_input, run_type = dedup_path, "genbank"
            except OSError as exc:
                print(f"[warn] {species}: GBFF clean skipped ({exc}); using original")

        success, msg = run_antismash(
            run_input, out_dir, input_type=run_type, extra_args=extra
        )
        if (not success) and run_type == "genbank" and _should_retry_dedupe_gbff(msg):
            # Second pass after clearing partial antiSMASH output.
            dedup_path = inp.with_suffix(inp.suffix + ".dedup.gbff")
            try:
                nfix = write_deduped_gbff(inp, dedup_path)
            except OSError as exc:
                print(f"[warn] {species}: could not write deduped GBFF: {exc}")
                nfix = -1
            if nfix >= 0:
                print(f"[retry] {species}: re-running cleaned GBFF ({nfix} fixes)")
                _clear_outdir_files(out_dir)
                success, msg = run_antismash(
                    dedup_path, out_dir, input_type="genbank", extra_args=extra
                )

        if success:
            _normalize_antismash_json(out_dir)
            if (out_dir / "genomic.json").exists():
                ok += 1
                print(f"[ok] antiSMASH {species}: {msg}")
            else:
                fail += 1
                print(
                    f"[warn] antiSMASH {species}: finished but genomic.json missing in {out_dir}"
                )
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
    parser.add_argument(
        "--only-species",
        nargs="+",
        default=None,
        help="Restrict to these species names (exact match as in selected_fungi.yaml)",
    )
    args = parser.parse_args(argv)

    targets = _load_targets(args.species_list)
    if args.only_species:
        wanted = {s.strip().lower() for s in args.only_species}
        targets = [t for t in targets if str(t.get("species", "")).strip().lower() in wanted]
        print(f"[filter] running {len(targets)} species: {[t['species'] for t in targets]}")
        if not targets:
            print("[error] --only-species matched nothing")
            return 2

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
