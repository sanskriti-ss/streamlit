"""Tests for antiSMASH runner helpers."""

from pathlib import Path

import pandas as pd

from genome_investigation.antismash_runner import (
    antismash_installed,
    build_summary_rows,
    parse_antismash_json,
    parse_output_directory,
)


FIXTURE = Path(__file__).parent / "fixtures" / "mock_antismash.json"


def test_parse_mock_antismash_json():
    parsed = parse_antismash_json(FIXTURE)
    assert parsed["bgc_count_total"] == 3
    assert parsed["nrps_count"] == 1
    assert parsed["pks_count"] == 1
    assert parsed["terpene_count"] == 1
    assert "cluster-1" in parsed["knownclusterblast_hits"]


def test_parse_output_directory(tmp_path: Path):
    import shutil

    shutil.copy(FIXTURE, tmp_path / "regions.json")
    parsed = parse_output_directory(tmp_path)
    assert parsed["bgc_count_total"] == 3


def test_missing_antismash_executable(monkeypatch):
    monkeypatch.setattr("genome_investigation.antismash_runner.shutil.which", lambda _: None)
    assert antismash_installed() is False


def test_skip_species_without_genome_dir(tmp_path: Path):
    jobs = [
        {
            "BacID": "1",
            "species": "Test sp.",
            "strain": "",
            "genome_accession": "GCA_TEST",
            "genome_dir": str(tmp_path / "missing"),
        }
    ]
    rows = build_summary_rows(jobs, tmp_path / "asm_out", run=False, parse_only=True)
    assert rows[0]["antismash_status"] == "skipped"


def test_summary_csv_columns(tmp_path: Path):
    out = tmp_path / "summary.csv"
    pd.DataFrame([{"BacID": "1", "species": "X", "bgc_count_total": 0}]).to_csv(out, index=False)
    df = pd.read_csv(out)
    assert "species" in df.columns
