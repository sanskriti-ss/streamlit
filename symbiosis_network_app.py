"""
Cross-kingdom metabolite synergy network explorer.

Launch:
    streamlit run symbiosis_network_app.py
"""

from __future__ import annotations

import streamlit as st

from utils.symbiosis_data import (
    LAYER_LABELS,
    list_metabolites,
    list_species,
    load_layers,
)
from utils.symbiosis_graph import (
    build_activity_matrices,
    build_metabolite_focus_graph,
    build_species_focus_graph,
    compute_synergy_pairs,
    plot_network,
    synergy_key,
)

st.set_page_config(page_title="Symbiosis Network", layout="wide")

LAYER_KEYS = list(LAYER_LABELS.keys())

st.title("Symbiosis metabolite network")
st.caption(
    "Explore production ↔ utilization synergy between bacteria and fungi. "
    "Works with partial data — enable only the layers you have."
)

with st.sidebar:
    st.header("Focus")
    focus_mode = st.radio("View centered on", ["Metabolite", "Species"], horizontal=True)

    st.header("Data layers")
    default_layers = ["bacteria_experimental"]
    selected_labels = st.multiselect(
        "Include layers",
        options=[LAYER_LABELS[k] for k in LAYER_KEYS],
        default=[LAYER_LABELS[k] for k in default_layers],
    )
    label_to_key = {v: k for k, v in LAYER_LABELS.items()}
    selected_layers = {label_to_key[lbl] for lbl in selected_labels if lbl in label_to_key}

    st.header("Filters")
    confidence_threshold = st.slider("Confidence threshold", 0.0, 1.0, 0.5, 0.05)
    aggregation = st.selectbox(
        "Bacteria aggregation",
        ["genus", "species"],
        index=0,
        help="Genus mode aggregates BacDive strains; species mode keeps top strains by breadth.",
    )
    max_entities = st.slider("Max bacterial entities", 50, 2000, 400, 50)
    min_edge_weight = st.slider("Min synergy edge weight", 1, 10, 1)
    max_hops = st.slider("Max hops from focal node", 1, 3, 2)

    st.header("Layer status")


@st.cache_data(show_spinner="Loading phenotype layers...")
def cached_load_layers(
    selected_layers_tuple: tuple,
    max_entities: int,
    aggregation: str,
    confidence_threshold: float,
):
    return load_layers(
        set(selected_layers_tuple),
        max_entities=max_entities,
        aggregation=aggregation,
        confidence_threshold=confidence_threshold,
    )


df, statuses = cached_load_layers(
    tuple(sorted(selected_layers)),
    max_entities,
    aggregation,
    confidence_threshold,
)

for st_status in statuses:
    icon = "✅" if st_status.available else "⚠️"
    st.sidebar.markdown(f"{icon} **{st_status.label}**")
    st.sidebar.caption(st_status.message or f"{st_status.row_count:,} rows")

if df.empty:
    st.warning(
        "No phenotype data loaded. Enable at least one available layer "
        "(e.g. Experimental bacteria) or lower the confidence threshold."
    )
    st.stop()

prod_mat, util_mat, meta = build_activity_matrices(df)
pairs_df, _ = compute_synergy_pairs(prod_mat, util_mat, meta, min_edge_weight=min_edge_weight)

col_search, col_opts = st.columns([2, 1])

if focus_mode == "Metabolite":
    with col_search:
        query = st.text_input("Search metabolite", placeholder="e.g. glucose, cellulose")
    choices = list_metabolites(df, query)
    with col_opts:
        selected_met = st.selectbox(
            "Select metabolite",
            choices if choices else ["—"],
            index=0,
        )

    if not choices or selected_met == "—":
        st.info("Type a metabolite name to search, then pick from the dropdown.")
        st.stop()

    nodes, edges = build_metabolite_focus_graph(df, selected_met, max_hops=max_hops)
    fig = plot_network(
        nodes,
        edges,
        center_id=selected_met,
        title=f"Metabolite focus: {selected_met}",
    )
else:
    with col_search:
        query = st.text_input("Search species", placeholder="e.g. Pseudomonas, Aspergillus")
    choices = list_species(df, query)
    with col_opts:
        selected_display = st.selectbox(
            "Select species",
            choices if choices else ["—"],
            index=0,
        )

    if not choices or selected_display == "—":
        st.info("Type a species name to search, then pick from the dropdown.")
        st.stop()

    name_part, _, kingdom_part = selected_display.rpartition(" (")
    kingdom_part = kingdom_part.rstrip(")")
    focal_key = synergy_key(name_part.strip(), kingdom_part.strip())
    nodes, edges = build_species_focus_graph(
        df,
        focal_key,
        pairs_df,
        max_hops=max_hops,
    )
    label = meta.get(focal_key, {}).get("species", name_part.strip())
    fig = plot_network(
        nodes,
        edges,
        center_id=focal_key,
        title=f"Species focus: {label}",
    )

st.plotly_chart(fig, use_container_width=True)

tab_pairs, tab_pheno, tab_help = st.tabs(["Synergy pairs", "Loaded phenotypes", "Help"])

with tab_pairs:
    if pairs_df.empty:
        st.info("No synergy pairs at current filters. Try lowering min edge weight or adding layers.")
    else:
        show = pairs_df.head(50)
        st.dataframe(show, use_container_width=True, hide_index=True)
        st.caption(
            f"Showing top {len(show)} of {len(pairs_df):,} pairs "
            f"(species A produces metabolites species B utilizes, and vice versa)."
        )

with tab_pheno:
    if focus_mode == "Metabolite":
        view = df[df["metabolite"].str.lower() == selected_met.lower()]
    else:
        row_keys = (
            df["species"].astype(str).str.strip().str.lower()
            + "|"
            + df["kingdom"].astype(str).str.strip().str.lower()
        )
        view = df[row_keys == focal_key]
    st.dataframe(
        view.sort_values(["activity", "confidence"], ascending=[True, False]),
        use_container_width=True,
        hide_index=True,
    )

with tab_help:
    st.markdown(
        """
**Layers**
- *Experimental bacteria*: BacDive Step3 production/utilization (positive tests only).
- *Predicted bacteria*: `phenotype_confidence.csv` from the genome follow-up pipeline.
- *Experimental fungi*: FUNG-GROWTH carbon-source utilization (`fungi_data/experimental/`).
- *Predicted fungi*: antiSMASH production + CAZyme utilization predictions.

**Synergy** = bidirectional complementarity: metabolites species A produces that species B utilizes.

**Degradation cross-feeding**: because BGC-based *production* (siderophores, terpenes, …)
and carbon-source *utilization* (glucose, cellulose, …) use different vocabularies, we
bridge them via extracellular degradation — a species that degrades a polymer (e.g.
cellulose) is treated as *producing* the released monomers (glucose, cellobiose) that
partners can consume. These inferred rows are tagged with a `+degradation` source layer.
See `fungi_data/metabolite_bridges.csv`.

**Fungi data not ready?** The app still runs with bacteria only. See `fungi_investigation/README.md`.
        """
    )
