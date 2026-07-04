"""Tests for BacDive phenotype metabolite extraction."""

from genome_investigation.bacdive_phenotype_extract import _activity_rows
from genome_investigation.targeted_gene_search import metabolite_search_terms


SAMPLE_STRAIN = {
    "Name and taxonomic classification": {
        "species": "Testus fictus",
        "strain designation": "ABC",
    },
    "Physiology and metabolism": {
        "metabolite utilization": [
            {"metabolite": "glucose", "utilization activity": "+", "Chebi-ID": 17234},
            {"metabolite": "acetate", "utilization activity": "-"},
        ],
        "metabolite production": [{"metabolite": "lactate", "production": "yes"}],
        "antibiotic resistance": [
            {"metabolite": "tetracycline", "is resistant": "yes", "resistance conc.": "100 µg/mL"}
        ],
    },
}


def test_activity_rows_counts():
    rows = _activity_rows(SAMPLE_STRAIN, "99", "Testus fictus", "ABC")
    activities = {r["activity"] for r in rows}
    assert "utilization" in activities
    assert "production" in activities
    assert "resistance" in activities
    assert any(r["metabolite"] == "glucose" and r["result"] == "+" for r in rows)


def test_metabolite_search_terms_aliases():
    terms = metabolite_search_terms("tetracycline")
    assert "tet" in terms or "tetracycline" in terms
