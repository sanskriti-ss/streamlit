"""Tests for FUNG-GROWTH ingest."""

from pathlib import Path

import pandas as pd

from fungi_investigation.fung_growth_ingest import ingest_fung_growth

SAMPLE = Path(__file__).resolve().parents[2] / "fungi_data" / "raw" / "fung_growth" / "sample_export.csv"


def test_ingest_sample_positive_only():
    df = ingest_fung_growth(SAMPLE)
    assert not df.empty
    assert set(df["activity"]) == {"utilization"}
    assert (df["confidence"] > 0).all()
    assert "glucose" in df["metabolite"].values or "D-glucose" in SAMPLE.read_text()


def test_glucose_alias_applied():
    aliases = Path(__file__).resolve().parents[2] / "fungi_data" / "metabolite_aliases.csv"
    df = ingest_fung_growth(SAMPLE, aliases_path=aliases)
    assert "glucose" in set(df["metabolite"].str.lower())
