"""
Symbiosis network tab: bacteria ↔ fungi production/utilization synergy.

Usable from the main app (`app.py`) or the standalone launcher
`symbiosis_network_app.py`. Kingdom scope presets let you view both kingdoms,
bacteria only, or fungi only.
"""

from __future__ import annotations

import pandas as pd
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

LAYER_KEYS = list(LAYER_LABELS.keys())

# Preset layer keys for the kingdom-scope radio.
SCOPE_PRESETS = {
    "Both kingdoms": [
        "bacteria_experimental",
        "fungi_experimental",
        "fungi_predicted",
    ],
    "Bacteria only": [
        "bacteria_experimental",
        "bacteria_predicted",
    ],
    "Fungi only": [
        "fungi_experimental",
        "fungi_predicted",
    ],
}


@st.cache_data(show_spinner="Loading phenotype layers...")
def _cached_load_layers(
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


def _sync_layers_for_scope(scope: str) -> None:
    """When kingdom scope changes, reset the layer multiselect to that preset."""
    prev = st.session_state.get("_symbiosis_prev_scope")
    if prev == scope and "symbiosis_layers" in st.session_state:
        return
    labels = [LAYER_LABELS[k] for k in SCOPE_PRESETS[scope] if k in LAYER_LABELS]
    st.session_state.symbiosis_layers = labels
    st.session_state._symbiosis_prev_scope = scope


def _partner_rows_for_focal(pairs_df, focal_key: str, *, top_n: int = 15):
    """Top synergy partners for one organism (by synergy_key)."""
    if pairs_df is None or pairs_df.empty or not focal_key:
        return pairs_df.iloc[0:0] if pairs_df is not None else None
    a = pairs_df[pairs_df["species_a_key"].astype(str) == focal_key].copy()
    b = pairs_df[pairs_df["species_b_key"].astype(str) == focal_key].copy()
    rows = []
    for _, r in a.iterrows():
        rows.append(
            {
                "partner": r.get("species_b"),
                "partner_kingdom": r.get("kingdom_b"),
                "synergy_score": r.get("synergy_score"),
                "focal_produces_partner_utilizes": r.get("a_produces_b_utilizes_n"),
                "partner_produces_focal_utilizes": r.get("b_produces_a_utilizes_n"),
                "metabolites_out": r.get("a_produces_b_utilizes_mets"),
                "metabolites_in": r.get("b_produces_a_utilizes_mets"),
            }
        )
    for _, r in b.iterrows():
        rows.append(
            {
                "partner": r.get("species_a"),
                "partner_kingdom": r.get("kingdom_a"),
                "synergy_score": r.get("synergy_score"),
                "focal_produces_partner_utilizes": r.get("b_produces_a_utilizes_n"),
                "partner_produces_focal_utilizes": r.get("a_produces_b_utilizes_n"),
                "metabolites_out": r.get("b_produces_a_utilizes_mets"),
                "metabolites_in": r.get("a_produces_b_utilizes_mets"),
            }
        )
    if not rows:
        return pairs_df.iloc[0:0]

    out = pd.DataFrame(rows).sort_values("synergy_score", ascending=False)
    return out.head(top_n)


def display(_data_frames=None) -> None:
    st.title("Symbiosis metabolite network")
    st.caption(
        "Explore production ↔ utilization synergy between bacteria and fungi. "
        "Use kingdom scope to view both, bacteria only, or fungi only."
    )

    # Deep-link support: ?tab=Symbiosis+Network&symbiosis_species=Aspergillus+niger
    q_species = (st.query_params.get("symbiosis_species") or "").strip()
    if q_species and "symbiosis_focus_mode" not in st.session_state:
        st.session_state.symbiosis_focus_mode = "Species"
        st.session_state.symbiosis_sp_query = q_species

    with st.sidebar:
        st.header("Symbiosis")
        scope = st.radio(
            "Kingdom scope",
            list(SCOPE_PRESETS.keys()),
            index=0,
            horizontal=False,
            key="symbiosis_scope",
            help="Quick preset for which kingdoms to load. Refine with Data layers below.",
        )
        _sync_layers_for_scope(scope)

        st.header("Focus")
        focus_mode = st.radio(
            "View centered on",
            ["Metabolite", "Species"],
            horizontal=True,
            key="symbiosis_focus_mode",
        )

        st.header("Data layers")
        selected_labels = st.multiselect(
            "Include layers",
            options=[LAYER_LABELS[k] for k in LAYER_KEYS],
            key="symbiosis_layers",
            help="Preset comes from Kingdom scope; add/remove layers freely.",
        )
        label_to_key = {v: k for k, v in LAYER_LABELS.items()}
        selected_layers = {
            label_to_key[lbl] for lbl in selected_labels if lbl in label_to_key
        }

        st.header("Filters")
        confidence_threshold = st.slider(
            "Confidence threshold",
            0.0,
            1.0,
            0.5,
            0.05,
            key="symbiosis_confidence",
        )
        aggregation = st.selectbox(
            "Bacteria aggregation",
            ["genus", "species"],
            index=0,
            key="symbiosis_aggregation",
            help="Genus mode aggregates BacDive strains; species mode keeps top strains by breadth.",
        )
        max_entities = st.slider(
            "Max bacterial entities",
            50,
            2000,
            100,
            50,
            key="symbiosis_max_entities",
            help="Lower = faster loads. Genus aggregation still keeps all genera "
            "that pass the breadth filter when at genus level.",
        )
        min_edge_weight = st.slider(
            "Min synergy edge weight",
            1,
            10,
            1,
            key="symbiosis_min_edge",
        )
        max_hops = st.slider(
            "Max hops from focal node",
            1,
            3,
            2,
            key="symbiosis_max_hops",
        )
        cross_kingdom_only = st.checkbox(
            "Show only cross-kingdom pairs",
            value=True,
            key="symbiosis_cross_only",
            help="In the Synergy pairs table, hide same-kingdom links "
            "(bacteria–bacteria or fungi–fungi).",
        )

        st.header("Layer status")

    if not selected_layers:
        st.warning("Select at least one data layer in the sidebar.")
        return

    df, statuses = _cached_load_layers(
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
        return

    kingdoms = sorted(df["kingdom"].dropna().unique())
    org_keys = df[["species", "kingdom"]].drop_duplicates()
    n_species = len(org_keys)
    n_bac = int((org_keys["kingdom"].astype(str) == "bacteria").sum())
    n_fungi = int((org_keys["kingdom"].astype(str) == "fungi").sum())
    n_mets = df["metabolite"].nunique()
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Kingdoms", ", ".join(kingdoms) if kingdoms else "—")
    m2.metric("Organisms", f"{n_species:,}")
    m3.metric("Metabolites", f"{n_mets:,}")
    m4.metric("Phenotype rows", f"{len(df):,}")
    st.caption(
        f"Loaded **{n_bac:,}** bacterial and **{n_fungi:,}** fungal organisms "
        f"at current layers/filters"
        + (
            " (fungal production from antiSMASH may still be catching up — "
            "utilization predictions are already included)."
            if n_fungi > 0
            else ""
        )
    )

    prod_mat, util_mat, meta = build_activity_matrices(df)
    pairs_df, _ = compute_synergy_pairs(
        prod_mat, util_mat, meta, min_edge_weight=min_edge_weight
    )

    if not pairs_df.empty and len(kingdoms) >= 2:
        cross = pairs_df[
            pairs_df["kingdom_a"].astype(str) != pairs_df["kingdom_b"].astype(str)
        ]
        st.caption(
            f"Synergy pairs: {len(pairs_df):,} total · {len(cross):,} cross-kingdom"
        )
    elif not pairs_df.empty:
        st.caption(f"Synergy pairs: {len(pairs_df):,} (same-kingdom view)")

    col_search, col_opts = st.columns([2, 1])

    if focus_mode == "Metabolite":
        with col_search:
            query = st.text_input(
                "Search metabolite",
                placeholder="e.g. glucose, cellulose",
                key="symbiosis_met_query",
            )
        choices = list_metabolites(df, query)
        with col_opts:
            selected_met = st.selectbox(
                "Select metabolite",
                choices if choices else ["—"],
                index=0,
                key="symbiosis_met_select",
            )

        if not choices or selected_met == "—":
            st.info("Type a metabolite name to search, then pick from the dropdown.")
            return

        nodes, edges = build_metabolite_focus_graph(
            df, selected_met, max_hops=max_hops
        )
        fig = plot_network(
            nodes,
            edges,
            center_id=selected_met,
            title=f"Metabolite focus: {selected_met}",
        )
        focal_key = None
    else:
        with col_search:
            query = st.text_input(
                "Search species",
                placeholder="e.g. Pseudomonas, Aspergillus",
                key="symbiosis_sp_query",
            )
        choices = list_species(df, query)
        with col_opts:
            selected_display = st.selectbox(
                "Select species",
                choices if choices else ["—"],
                index=0,
                key="symbiosis_sp_select",
            )

        if not choices or selected_display == "—":
            st.info("Type a species name to search, then pick from the dropdown.")
            return

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
        selected_met = None

    st.plotly_chart(fig, use_container_width=True)

    if focus_mode == "Species" and focal_key:
        partners = _partner_rows_for_focal(pairs_df, focal_key, top_n=15)
        if partners is not None and not partners.empty:
            partner_label = meta.get(focal_key, {}).get("species", focal_key)
            st.subheader(f"Top partners for {partner_label}")
            cross = partners[
                partners["partner_kingdom"].astype(str).str.lower()
                != str(meta.get(focal_key, {}).get("kingdom", "")).lower()
            ]
            show_partners = cross if not cross.empty else partners
            if not cross.empty:
                st.caption(
                    f"Showing top {len(show_partners)} cross-kingdom partner(s) by synergy score."
                )
            else:
                st.caption("No cross-kingdom partners at current filters; showing same-kingdom.")
            st.dataframe(show_partners, use_container_width=True, hide_index=True)

    tab_pairs, tab_pheno, tab_help = st.tabs(
        ["Synergy pairs", "Loaded phenotypes", "Help"]
    )

    with tab_pairs:
        if pairs_df.empty:
            st.info(
                "No synergy pairs at current filters. Try lowering min edge weight "
                "or adding layers / switching kingdom scope."
            )
        else:
            table = pairs_df
            if cross_kingdom_only and {"kingdom_a", "kingdom_b"} <= set(table.columns):
                table = table[
                    table["kingdom_a"].astype(str) != table["kingdom_b"].astype(str)
                ]
            if table.empty:
                st.info(
                    "No cross-kingdom pairs at current filters. "
                    "Unset “Show only cross-kingdom pairs” or add fungi + bacteria layers."
                )
            else:
                show = table.head(50)
                st.dataframe(show, use_container_width=True, hide_index=True)
                scope_note = "cross-kingdom " if cross_kingdom_only else ""
                st.caption(
                    f"Showing top {len(show)} of {len(table):,} {scope_note}pairs "
                    f"(species A produces metabolites species B utilizes, and vice versa)."
                )

    with tab_pheno:
        if focus_mode == "Metabolite" and selected_met:
            view = df[df["metabolite"].str.lower() == selected_met.lower()]
        elif focal_key:
            row_keys = (
                df["species"].astype(str).str.strip().str.lower()
                + "|"
                + df["kingdom"].astype(str).str.strip().str.lower()
            )
            view = df[row_keys == focal_key]
        else:
            view = df.head(0)
        st.dataframe(
            view.sort_values(["activity", "confidence"], ascending=[True, False]),
            use_container_width=True,
            hide_index=True,
        )

    with tab_help:
        st.markdown(
            """
**Kingdom scope**
- *Both kingdoms*: bacteria + fungi layers (cross-kingdom synergy).
- *Bacteria only*: BacDive experimental ± predicted bacterial phenotypes.
- *Fungi only*: FUNG-GROWTH experimental + genomic predicted fungi.

**Layers**
- *Experimental bacteria*: BacDive Step3 production/utilization (positive tests only).
- *Predicted bacteria*: `phenotype_confidence.csv` from the genome follow-up pipeline.
- *Experimental fungi*: FUNG-GROWTH carbon-source utilization (`fungi_data/experimental/`).
- *Predicted fungi*: antiSMASH production + CAZyme utilization predictions.

**Synergy** = bidirectional complementarity: metabolites species A produces that species B utilizes.
Organisms with only production *or* only utilization evidence still participate.

**Degradation cross-feeding**: BGC-based production and carbon-source utilization use
different vocabularies; polymer degraders are treated as releasing monomers partners
can consume (`fungi_data/metabolite_bridges.csv`).

**Fungi data not ready?** The tab still runs with bacteria only. See `fungi_investigation/README.md`.
            """
        )
