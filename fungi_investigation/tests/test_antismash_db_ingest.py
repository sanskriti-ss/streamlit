"""Tests for antiSMASH DB ingest."""

import json
import shutil
from pathlib import Path

from fungi_investigation.antismash_db_ingest import ingest_antismash_jsons

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "genome_investigation"
    / "tests"
    / "fixtures"
    / "mock_antismash.json"
)


def test_ingest_mock_json(tmp_path):
    species_dir = tmp_path / "aspergillus_niger" / "GCA_000001"
    species_dir.mkdir(parents=True)
    shutil.copy(FIXTURE, species_dir / "genomic.json")
    df = ingest_antismash_jsons(tmp_path, kingdom="fungi")
    assert not df.empty
    assert (df["activity"] == "production").all()
    assert (df["kingdom"] == "fungi").all()
    assert (df["confidence"] > 0).all()
