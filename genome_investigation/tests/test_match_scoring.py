"""Tests for genome match confidence scoring."""

from genome_investigation.match_scoring import pick_best_candidate, score_bacdive_direct, score_ncbi_match


def test_bacdive_direct_confidence():
    conf, note = score_bacdive_direct()
    assert conf == 1.0
    assert "BacDive" in note


def test_ncbi_species_and_strain_match():
    row = {"species": "Cohnella algarum", "genus": "Cohnella", "strain": "Pch-40", "type_strain": "yes"}
    cand = {
        "organism_name": "Cohnella algarum",
        "assembly_name": "ASM1693751v1 assembly for Cohnella algarum Pch-40",
        "assembly_accession": "GCA_016937515",
    }
    score, note = score_ncbi_match(row, cand, allow_species_only=False)
    assert score >= 0.9


def test_ncbi_species_only_requires_flag():
    row = {"species": "Unknown species xyz", "genus": "Escherichia", "strain": ""}
    cand = {"organism_name": "Escherichia coli K-12", "assembly_name": "reference"}
    score_no, _ = score_ncbi_match(row, cand, allow_species_only=False)
    score_yes, _ = score_ncbi_match(row, cand, allow_species_only=True)
    assert score_yes > score_no
    assert score_no == 0.0


def test_pick_best_candidate():
    row = {"species": "Test species", "genus": "Test", "strain": "A1"}
    cands = [
        {"organism_name": "Other", "assembly_name": "x"},
        {"organism_name": "Test species", "assembly_name": "Test species strain A1"},
    ]
    best, score, _ = pick_best_candidate(row, cands, allow_species_only=True)
    assert best is not None
    assert score > 0.5
