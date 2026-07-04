"""Tests for symbiosis data loading."""

from utils.symbiosis_data import (
    list_species,
    load_layers,
    load_metabolite_aliases,
    normalize_metabolite,
    prettify_metabolite,
    species_display_name,
)
from pathlib import Path


def test_load_bacteria_experimental_genus():
    df, statuses = load_layers(
        {"bacteria_experimental"},
        max_entities=50,
        aggregation="genus",
        confidence_threshold=0.5,
    )
    bac = [s for s in statuses if s.name == "bacteria_experimental"][0]
    assert bac.available
    assert not df.empty
    assert (df["kingdom"] == "bacteria").all()


def test_missing_fungi_graceful():
    df, statuses = load_layers({"fungi_experimental"}, confidence_threshold=0.5)
    fung = [s for s in statuses if s.name == "fungi_experimental"][0]
    # May be available if sample was ingested
    if not fung.available:
        assert df.empty


def test_aliases_file_loads():
    path = Path(__file__).resolve().parents[2] / "fungi_data" / "metabolite_aliases.csv"
    aliases = load_metabolite_aliases(path)
    assert aliases.get("d-glucose") == "glucose"


def test_normalize_lowercases_without_alias():
    """Differently-cased names with no alias must still collapse to lowercase."""
    assert normalize_metabolite("Fructose", {}) == "fructose"
    assert normalize_metabolite("Arabinose", {}) == "arabinose"
    assert normalize_metabolite("Cellulose", {"cellulose": "cellulose"}) == "cellulose"


def test_prettify_preserves_bgc_display():
    assert prettify_metabolite("polyketide (type i)") == "Polyketide (type I)"
    assert prettify_metabolite("ripp") == "RiPP"
    assert prettify_metabolite("siderophore") == "Siderophore"
    # carbon sources stay lowercase (experimental convention)
    assert prettify_metabolite("glucose") == "glucose"


def test_species_display_name_prefers_proper_casing():
    assert species_display_name(["aspergillus niger", "Aspergillus niger"]) == "Aspergillus niger"
    # only lowercase available → capitalize the genus
    assert species_display_name(["aspergillus niger"]) == "Aspergillus niger"
    # strain designations are not mangled
    assert species_display_name(["Streptomyces sp. NPDC059396"]) == "Streptomyces sp. NPDC059396"


def test_species_dropdown_dedupes_across_layers():
    """A. niger has experimental (title case) + predicted (lowercase) rows → one entry."""
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        confidence_threshold=0.5,
    )
    species_items = list_species(df)
    aspergillus = [s for s in species_items if "aspergillus niger" in s.lower()]
    assert len(aspergillus) == 1
    assert aspergillus[0] == "Aspergillus niger (fungi)"
    # no duplicate that differs only by casing
    lowered = [s.lower() for s in species_items]
    assert len(lowered) == len(set(lowered))


def test_metabolite_casing_merges_across_fungi_layers():
    """Experimental (lowercase) and predicted (was capitalized) substrates merge."""
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        confidence_threshold=0.5,
    )
    mets = set(df["metabolite"])
    # no capitalized carbon-source duplicates remain
    assert "Cellulose" not in mets
    assert "Glucose" not in mets
    assert "Fructose" not in mets
    # canonical lowercase forms are present
    assert "cellulose" in mets
    assert "glucose" in mets
