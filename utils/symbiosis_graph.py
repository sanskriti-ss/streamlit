"""
Graph construction and synergy analysis for the symbiosis network app.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from utils.symbiosis_data import prettify_metabolite, species_display_name

KINGDOM_COLORS = {
    "bacteria": "#4C78A8",
    "fungi": "#54A24B",
}

LAYER_SHAPES = {
    "bacteria_experimental": "circle",
    "bacteria_predicted": "circle-open",
    "fungi_experimental": "diamond",
    "fungi_predicted": "diamond-open",
}

ACTIVITY_COLORS = {
    "production": "#7B9E89",
    "utilization": "#D4858C",
}


@dataclass
class GraphEdge:
    source: str
    target: str
    weight: int
    metabolites: List[str]
    edge_type: str  # prod_util, util_link, prod_link, synergy


@dataclass
class GraphNode:
    node_id: str
    label: str
    node_type: str  # species, metabolite
    kingdom: str = ""
    source_layer: str = ""
    activity: str = ""


def synergy_key(species: str, kingdom: str) -> str:
    """Canonical key for prod/util synergy across layers with different entity_ids."""
    return f"{str(species).strip().lower()}|{str(kingdom).strip().lower()}"


def synergy_key_from_row(row: pd.Series) -> str:
    return synergy_key(str(row.get("species", "")), str(row.get("kingdom", "")))


def _binarize_phenotypes(df: pd.DataFrame, *, index_col: str = "entity_key") -> pd.DataFrame:
    """Pivot to index_col × metabolite for one activity."""
    if df.empty:
        return pd.DataFrame()
    sub = df.copy()
    sub["value"] = 1
    wide = sub.pivot_table(
        index=index_col,
        columns="metabolite",
        values="value",
        aggfunc="max",
        fill_value=0,
    )
    return (wide > 0).astype(np.int8)


def _align_activity_matrices(
    prod: pd.DataFrame, util: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict prod/util to shared entity keys and metabolite columns."""
    if prod.empty or util.empty:
        return prod, util
    shared_mets = sorted(set(prod.columns) & set(util.columns))
    shared_keys = sorted(set(prod.index) & set(util.index))
    if not shared_mets or not shared_keys:
        return prod.iloc[0:0], util.iloc[0:0]
    return prod.loc[shared_keys, shared_mets], util.loc[shared_keys, shared_mets]


def build_activity_matrices(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, dict]]:
    """Return production matrix, utilization matrix, and entity metadata."""
    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), {}

    work = df.copy()
    work["_synergy_key"] = work.apply(synergy_key_from_row, axis=1)

    meta: Dict[str, dict] = {}
    for key, grp in work.groupby("_synergy_key"):
        row0 = grp.iloc[0]
        meta[str(key)] = {
            "species": species_display_name(grp["species"]),
            "genus": str(row0.get("genus", "")),
            "kingdom": str(row0["kingdom"]),
            "source_layer": str(row0["source_layer"]),
        }

    prod = _binarize_phenotypes(work[work["activity"] == "production"], index_col="_synergy_key")
    util = _binarize_phenotypes(work[work["activity"] == "utilization"], index_col="_synergy_key")
    prod, util = _align_activity_matrices(prod, util)
    return prod, util, meta


def compute_synergy_pairs(
    prod: pd.DataFrame,
    util: pd.DataFrame,
    meta: Dict[str, dict],
    *,
    min_edge_weight: int = 1,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Directed prod→util synergy between species."""
    if prod.empty or util.empty:
        return pd.DataFrame(), np.array([])

    prod, util = _align_activity_matrices(prod, util)
    if prod.empty or util.empty:
        return pd.DataFrame(), np.array([])

    keys = list(prod.index)
    p = prod.to_numpy(dtype=np.int8)
    u = util.to_numpy(dtype=np.int8)
    mets = list(prod.columns)

    directed = p @ u.T
    synergy = directed + directed.T
    np.fill_diagonal(synergy, 0)

    pair_rows: List[dict] = []
    n = len(keys)
    for i in range(n):
        for j in range(i + 1, n):
            score = int(synergy[i, j])
            if score < min_edge_weight:
                continue
            ab = int(directed[i, j])
            ba = int(directed[j, i])
            if ab <= 0 and ba <= 0:
                continue
            a_mets = [prettify_metabolite(mets[k]) for k in np.where((p[i] > 0) & (u[j] > 0))[0]]
            b_mets = [prettify_metabolite(mets[k]) for k in np.where((p[j] > 0) & (u[i] > 0))[0]]
            pair_rows.append(
                {
                    "species_a_key": keys[i],
                    "species_b_key": keys[j],
                    "species_a": meta.get(keys[i], {}).get("species", keys[i]),
                    "species_b": meta.get(keys[j], {}).get("species", keys[j]),
                    "kingdom_a": meta.get(keys[i], {}).get("kingdom", ""),
                    "kingdom_b": meta.get(keys[j], {}).get("kingdom", ""),
                    "synergy_score": score,
                    "a_produces_b_utilizes_n": ab,
                    "b_produces_a_utilizes_n": ba,
                    "a_produces_b_utilizes_mets": ", ".join(a_mets[:20]),
                    "b_produces_a_utilizes_mets": ", ".join(b_mets[:20]),
                }
            )

    pairs_df = pd.DataFrame(pair_rows)
    if not pairs_df.empty:
        pairs_df = pairs_df.sort_values("synergy_score", ascending=False)
    return pairs_df, directed


def synergy_directed_edges(
    pairs_df: pd.DataFrame, *, max_edges: int = 80
) -> List[GraphEdge]:
    edges: List[GraphEdge] = []
    if pairs_df.empty:
        return edges
    for _, r in pairs_df.head(max_edges).iterrows():
        ab = int(r.get("a_produces_b_utilizes_n", 0))
        ba = int(r.get("b_produces_a_utilizes_n", 0))
        if ab > 0:
            edges.append(
                GraphEdge(
                    source=str(r["species_a_key"]),
                    target=str(r["species_b_key"]),
                    weight=ab,
                    metabolites=str(r.get("a_produces_b_utilizes_mets", "")).split(", "),
                    edge_type="synergy",
                )
            )
        if ba > 0:
            edges.append(
                GraphEdge(
                    source=str(r["species_b_key"]),
                    target=str(r["species_a_key"]),
                    weight=ba,
                    metabolites=str(r.get("b_produces_a_utilizes_mets", "")).split(", "),
                    edge_type="synergy",
                )
            )
    return edges


def _bfs_neighbors(
    adj: Dict[str, Set[str]], start: str, max_hops: int
) -> Set[str]:
    if start not in adj:
        return {start}
    visited = {start}
    queue: deque = deque([(start, 0)])
    while queue:
        node, depth = queue.popleft()
        if depth >= max_hops:
            continue
        for nb in adj.get(node, set()):
            if nb not in visited:
                visited.add(nb)
                queue.append((nb, depth + 1))
    return visited


def build_metabolite_focus_graph(
    df: pd.DataFrame,
    metabolite: str,
    *,
    max_hops: int = 2,
) -> Tuple[List[GraphNode], List[GraphEdge]]:
    """Center metabolite; link to species that produce or utilize it."""
    sub = df[df["metabolite"].str.lower() == metabolite.lower()]
    if sub.empty:
        return [], []

    nodes: Dict[str, GraphNode] = {
        metabolite: GraphNode(metabolite, prettify_metabolite(metabolite), "metabolite")
    }
    edges: List[GraphEdge] = []

    sub = sub.copy()
    sub["_sk"] = sub.apply(synergy_key_from_row, axis=1)
    for sk_val, grp in sub.groupby("_sk"):
        sk_val = str(sk_val)
        activities = set(grp["activity"].astype(str))
        row0 = grp.iloc[0]
        if sk_val not in nodes:
            nodes[sk_val] = GraphNode(
                sk_val,
                species_display_name(grp["species"]),
                "species",
                kingdom=str(row0["kingdom"]),
                source_layer=str(row0["source_layer"]),
            )
        if "production" in activities:
            edges.append(GraphEdge(sk_val, metabolite, 1, [metabolite], "prod_link"))
        if "utilization" in activities:
            edges.append(GraphEdge(metabolite, sk_val, 1, [metabolite], "util_link"))

    adj: Dict[str, Set[str]] = defaultdict(set)
    for e in edges:
        adj[e.source].add(e.target)
        adj[e.target].add(e.source)

    keep = _bfs_neighbors(adj, metabolite, max_hops)
    nodes = {k: v for k, v in nodes.items() if k in keep}
    edges = [e for e in edges if e.source in keep and e.target in keep]
    return list(nodes.values()), edges


def build_species_focus_graph(
    df: pd.DataFrame,
    focal_key: str,
    pairs_df: pd.DataFrame,
    *,
    max_hops: int = 2,
    max_partner_edges: int = 40,
) -> Tuple[List[GraphNode], List[GraphEdge]]:
    """Center one organism (by synergy key); metabolite links + synergy partners.

    ``focal_key`` is a ``species|kingdom`` synergy key, so all layers/strains of the
    organism are merged into a single focal node (matching the synergy matrices).
    """
    sk = df.apply(synergy_key_from_row, axis=1)
    sub = df[sk == focal_key]
    if sub.empty:
        return [], []

    row0 = sub.iloc[0]
    nodes: Dict[str, GraphNode] = {
        focal_key: GraphNode(
            focal_key,
            species_display_name(sub["species"]),
            "species",
            kingdom=str(row0["kingdom"]),
            source_layer=str(row0["source_layer"]),
        )
    }
    edges: List[GraphEdge] = []

    for _, row in sub.iterrows():
        met = str(row["metabolite"])
        activity = str(row["activity"])
        if met not in nodes:
            nodes[met] = GraphNode(met, prettify_metabolite(met), "metabolite", activity=activity)
        if activity == "production":
            edges.append(GraphEdge(focal_key, met, 1, [met], "prod_link"))
        elif activity == "utilization":
            edges.append(GraphEdge(met, focal_key, 1, [met], "util_link"))

    partner_edges = synergy_directed_edges(pairs_df, max_edges=max_partner_edges)
    for pe in partner_edges:
        if focal_key not in (pe.source, pe.target):
            continue
        for nk in (pe.source, pe.target):
            if nk == focal_key or nk in nodes:
                continue
            match = df[sk == nk]
            if not match.empty:
                nodes[nk] = GraphNode(
                    nk,
                    species_display_name(match["species"]),
                    "species",
                    kingdom=str(match.iloc[0]["kingdom"]),
                    source_layer=str(match.iloc[0]["source_layer"]),
                )
        edges.append(
            GraphEdge(pe.source, pe.target, pe.weight, pe.metabolites, pe.edge_type)
        )

    adj: Dict[str, Set[str]] = defaultdict(set)
    for e in edges:
        adj[e.source].add(e.target)
        adj[e.target].add(e.source)

    keep = _bfs_neighbors(adj, focal_key, max_hops)
    nodes = {k: v for k, v in nodes.items() if k in keep}
    edges = [e for e in edges if e.source in keep and e.target in keep]
    return list(nodes.values()), edges


def _layout_circle(nodes: List[GraphNode], center_id: Optional[str] = None) -> Dict[str, Tuple[float, float]]:
    if not nodes:
        return {}
    n = len(nodes)
    positions: Dict[str, Tuple[float, float]] = {}
    if center_id and any(nd.node_id == center_id for nd in nodes):
        positions[center_id] = (0.0, 0.0)
        others = [nd for nd in nodes if nd.node_id != center_id]
        radius = 1.2 + 0.05 * min(len(others), 30)
        for i, nd in enumerate(others):
            angle = 2 * np.pi * i / max(len(others), 1)
            positions[nd.node_id] = (radius * np.cos(angle), radius * np.sin(angle))
        return positions

    radius = 1.0 + 0.05 * min(n, 30)
    for i, nd in enumerate(nodes):
        angle = 2 * np.pi * i / max(n, 1)
        positions[nd.node_id] = (radius * np.cos(angle), radius * np.sin(angle))
    return positions


def plot_network(
    nodes: List[GraphNode],
    edges: List[GraphEdge],
    *,
    center_id: Optional[str] = None,
    title: str = "Symbiosis network",
) -> go.Figure:
    if not nodes:
        fig = go.Figure()
        fig.update_layout(title=title, annotations=[dict(text="No nodes to display", showarrow=False)])
        return fig

    pos = _layout_circle(nodes, center_id=center_id)
    node_idx = {nd.node_id: nd for nd in nodes}

    edge_traces = []
    for e in edges:
        if e.source not in pos or e.target not in pos:
            continue
        x0, y0 = pos[e.source]
        x1, y1 = pos[e.target]
        color = "#888888" if e.edge_type in ("prod_link", "util_link") else "#E45756"
        edge_traces.append(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode="lines",
                line=dict(width=1 + min(e.weight, 5), color=color),
                hoverinfo="text",
                text=f"{e.source} → {e.target} ({e.weight}): {', '.join(e.metabolites[:5])}",
                showlegend=False,
            )
        )

    xs, ys, texts, colors, sizes, symbols = [], [], [], [], [], []
    for nd in nodes:
        x, y = pos[nd.node_id]
        xs.append(x)
        ys.append(y)
        if nd.node_type == "metabolite":
            colors.append(ACTIVITY_COLORS.get(nd.activity, "#9B8BA8"))
            symbols.append("square")
            sizes.append(14)
            texts.append(f"Metabolite: {nd.label}")
        else:
            colors.append(KINGDOM_COLORS.get(nd.kingdom, "#666666"))
            base_layer = nd.source_layer.split("+", 1)[0]
            symbols.append(LAYER_SHAPES.get(base_layer, "circle"))
            sizes.append(18 if nd.node_id == center_id else 12)
            texts.append(
                f"{nd.label}<br>kingdom={nd.kingdom}<br>layer={nd.source_layer}"
            )

    node_trace = go.Scatter(
        x=xs,
        y=ys,
        mode="markers+text",
        text=[nd.label[:28] for nd in nodes],
        textposition="top center",
        textfont=dict(size=8),
        marker=dict(color=colors, size=sizes, symbol=symbols, line=dict(width=1, color="#333")),
        hovertext=texts,
        hoverinfo="text",
        showlegend=False,
    )

    fig = go.Figure(data=edge_traces + [node_trace])
    fig.update_layout(
        title=title,
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=20, r=20, t=50, b=20),
        height=620,
    )
    return fig
