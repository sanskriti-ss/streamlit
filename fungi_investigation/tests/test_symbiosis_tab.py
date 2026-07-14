"""Smoke tests for the integrated Symbiosis Network tab."""

from tabs.symbiosis_network import SCOPE_PRESETS
from utils.symbiosis_data import load_layers
from utils.symbiosis_graph import (
    build_activity_matrices,
    build_metabolite_focus_graph,
    compute_synergy_pairs,
)


def test_scope_presets_cover_kingdom_separations():
    assert set(SCOPE_PRESETS) == {"Both kingdoms", "Bacteria only", "Fungi only"}
    both = set(SCOPE_PRESETS["Both kingdoms"])
    bac = set(SCOPE_PRESETS["Bacteria only"])
    fung = set(SCOPE_PRESETS["Fungi only"])
    assert both & bac
    assert both & fung
    assert not bac & fung


def test_bacteria_only_scope_loads_bacteria():
    layers = set(SCOPE_PRESETS["Bacteria only"])
    df, statuses = load_layers(
        layers, max_entities=50, aggregation="genus", confidence_threshold=0.5
    )
    assert any(s.available for s in statuses)
    assert not df.empty
    assert set(df["kingdom"].unique()) == {"bacteria"}


def test_fungi_only_scope_loads_fungi():
    layers = set(SCOPE_PRESETS["Fungi only"])
    df, statuses = load_layers(layers, confidence_threshold=0.5)
    assert any(s.available for s in statuses)
    assert not df.empty
    assert set(df["kingdom"].unique()) == {"fungi"}


def test_both_kingdoms_scope_cross_synergy():
    layers = set(SCOPE_PRESETS["Both kingdoms"])
    df, statuses = load_layers(
        layers, max_entities=80, aggregation="genus", confidence_threshold=0.5
    )
    assert any(s.available for s in statuses if "bacteria" in s.name)
    assert any(s.available for s in statuses if "fungi" in s.name)
    assert {"bacteria", "fungi"} <= set(df["kingdom"].unique())

    prod, util, meta = build_activity_matrices(df)
    pairs, _ = compute_synergy_pairs(prod, util, meta, min_edge_weight=1)
    assert not pairs.empty
    cross = pairs[pairs["kingdom_a"] != pairs["kingdom_b"]]
    assert not cross.empty

    # metabolite focus graph builds for a shared carbon source
    nodes, edges = build_metabolite_focus_graph(df, "glucose", max_hops=1)
    assert nodes
    assert edges


def test_app_imports_symbiosis_tab():
    import app as main_app

    assert "Symbiosis Network" in main_app.available_tabs
    from tabs import symbiosis_network

    assert callable(symbiosis_network.display)
