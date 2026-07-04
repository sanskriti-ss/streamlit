"""Tests for synergy matrix alignment."""

from utils.symbiosis_data import load_layers
from utils.symbiosis_graph import (
    build_activity_matrices,
    build_species_focus_graph,
    compute_synergy_pairs,
    synergy_key,
)


def test_predicted_fungi_does_not_crash_synergy():
    """Mixed entity_ids across fungi layers must not break matmul."""
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        max_entities=100,
        aggregation="species",
        confidence_threshold=0.5,
    )
    prod, util, meta = build_activity_matrices(df)
    pairs, directed = compute_synergy_pairs(prod, util, meta)
    assert directed.shape == (0,) or directed.ndim == 2
    # aspergillus niger has both prod (predicted) and util (experimental) rows
    assert "aspergillus niger|fungi" in meta


def test_bacteria_fungi_combined_synergy():
    df, _ = load_layers(
        {"bacteria_experimental", "fungi_experimental", "fungi_predicted"},
        max_entities=50,
        aggregation="genus",
        confidence_threshold=0.5,
    )
    prod, util, meta = build_activity_matrices(df)
    pairs, _ = compute_synergy_pairs(prod, util, meta, min_edge_weight=1)
    assert prod.shape[0] == util.shape[0]
    assert prod.shape[1] == util.shape[1]
    cross = pairs[
        ((pairs.kingdom_a == "fungi") & (pairs.kingdom_b == "bacteria"))
        | ((pairs.kingdom_a == "bacteria") & (pairs.kingdom_b == "fungi"))
    ]
    assert not cross.empty


def test_cross_kingdom_synergy_without_degradation_bridges(tmp_path):
    """Experimental fungi (utilization-only) must still cross-feed with bacteria."""
    cfg = tmp_path / "paths.yaml"
    cfg.write_text(
        "bacteria_experimental:\n"
        "  production: species_data/step3_met_prod_exploded.csv.zip\n"
        "  utilization: species_data/step3_met_util_exploded.csv.zip\n"
        "fungi_experimental: fungi_data/experimental/fungi_phenotypes_long.csv\n"
        "metabolite_aliases: fungi_data/metabolite_aliases.csv\n",
        encoding="utf-8",
    )
    df, _ = load_layers(
        {"bacteria_experimental", "fungi_experimental"},
        config_path=cfg,
        max_entities=50,
        aggregation="genus",
        confidence_threshold=0.5,
    )
    # experimental fungi have utilization only — no degradation bridge rows
    assert not df["source_layer"].str.contains(r"\+degradation", regex=True).any()
    prod, util, meta = build_activity_matrices(df)
    pairs, _ = compute_synergy_pairs(prod, util, meta, min_edge_weight=1)
    cross = pairs[
        ((pairs.kingdom_a == "fungi") & (pairs.kingdom_b == "bacteria"))
        | ((pairs.kingdom_a == "bacteria") & (pairs.kingdom_b == "fungi"))
    ]
    assert not cross.empty
    # bacteria produce carbon sources that experimental fungi utilize
    met_blob = " ".join(cross["a_produces_b_utilizes_mets"].astype(str)).lower()
    met_blob += " " + " ".join(cross["b_produces_a_utilizes_mets"].astype(str)).lower()
    assert "glucose" in met_blob or "cellulose" in met_blob


def test_degradation_bridge_creates_fungi_synergy():
    """Polymer→monomer degradation bridge must yield non-empty fungi synergy."""
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        max_entities=100,
        aggregation="species",
        confidence_threshold=0.5,
    )
    # inferred degradation-derived production rows are present
    assert df["source_layer"].str.contains(r"\+degradation", regex=True).any()
    # they land in the shared carbon-source vocabulary (e.g. glucose)
    derived = df[df["source_layer"].str.contains(r"\+degradation", regex=True)]
    assert "production" in set(derived["activity"])
    assert "glucose" in set(derived["metabolite"].str.lower())

    prod, util, meta = build_activity_matrices(df)
    pairs, _ = compute_synergy_pairs(prod, util, meta, min_edge_weight=1)
    assert not pairs.empty


def test_species_focus_merges_layers_into_one_node():
    """Selecting a species by synergy key shows all its layers under one focal node."""
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        confidence_threshold=0.5,
    )
    prod, util, meta = build_activity_matrices(df)
    pairs, _ = compute_synergy_pairs(prod, util, meta, min_edge_weight=1)

    focal = synergy_key("Aspergillus niger", "fungi")
    nodes, edges = build_species_focus_graph(df, focal, pairs, max_hops=2)

    species_nodes = [n for n in nodes if n.node_type == "species"]
    focal_nodes = [n for n in species_nodes if n.node_id == focal]
    # exactly one focal node despite experimental + predicted rows
    assert len(focal_nodes) == 1
    assert focal_nodes[0].label == "Aspergillus niger"
    # it carries both produced (degradation) and utilized metabolite links
    assert any(e.edge_type == "prod_link" for e in edges)
    assert any(e.edge_type == "util_link" for e in edges)


def test_bridges_disabled_when_config_missing(tmp_path):
    """No metabolite_bridges config → no derived degradation rows."""
    cfg = tmp_path / "paths.yaml"
    cfg.write_text(
        "fungi_experimental: fungi_data/experimental/fungi_phenotypes_long.csv\n"
        "fungi_predicted_production: fungi_data/predicted/fungi_phenotype_confidence.csv\n"
        "fungi_predicted_utilization: fungi_data/predicted/fungi_utilization_confidence.csv\n"
        "metabolite_aliases: fungi_data/metabolite_aliases.csv\n",
        encoding="utf-8",
    )
    df, _ = load_layers(
        {"fungi_experimental", "fungi_predicted"},
        config_path=cfg,
        confidence_threshold=0.5,
    )
    assert not df["source_layer"].str.contains(r"\+degradation", regex=True).any()
