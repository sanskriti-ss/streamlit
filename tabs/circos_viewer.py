"""
Upload Circos-style chromosome and link `.txt` files and render an interactive circle graph.
"""
import streamlit as st
import matplotlib.pyplot as plt

from utils.circos_graph import (
    aggregate_undirected_edges,
    build_circle_graph_figure,
    parse_chromosome_text,
    parse_links_text,
)
from utils.tooltip_title import display_title_with_tooltip

_DEFAULT_NODE_COLORS = [
    "pred",
    "orange",
    "yellow",
    "green",
    "blue",
    "purple",
    "grey",
    "pred",
    "porange",
    "pyellow",
    "pgreen",
    "pblue",
    "ppurple",
]


def _merge_node_colors_in_order(chromosome_colors: dict, ordered_node_list: list) -> dict:
    """Preserve ring order; assign default Circos colour names for nodes missing from the chromosome file."""
    out = {}
    for i, name in enumerate(ordered_node_list):
        if name in chromosome_colors:
            out[name] = chromosome_colors[name]
        else:
            out[name] = _DEFAULT_NODE_COLORS[i % len(_DEFAULT_NODE_COLORS)]
    return out


def display():
    display_title_with_tooltip(
        title_text="Circle graph (from Circos files)",
        sample_image_filename="circos_example.png",
        description_text=(
            "Upload the chromosome (node ring) file and the links file produced by the **Circos** tab "
            "or your external pipeline.\n\n"
            "- **Chromosome file** lines look like: `chr - Genus Genus 0 100 colourname`\n"
            "- **Links file** lines look like: `GenusA 100 0 GenusB 100 0 color=orange,thickness=15p`\n\n"
            "Genera appear as **coloured arc segments** on the ring (no dot markers). "
            "Each label sits a **fixed radial offset** past its arc; optional faint spokes match that ray. "
            "Chords use **semi-transparent** strokes, **sqrt-scaled** width from `thickness`, and draw order "
            "so thinner links stay visible on top."
        ),
    )
    st.markdown("---")

    col_a, col_b = st.columns(2)
    with col_a:
        chr_file = st.file_uploader(
            "Chromosome file (`*_chromosome_file.txt`)",
            type=["txt"],
            help="Defines each genus on the ring and its Circos colour name.",
        )
    with col_b:
        links_file = st.file_uploader(
            "Links file (`*_alloutputs.txt` or filtered `*_Triples.txt`, etc.)",
            type=["txt"],
            help="Defines edges: genus pairs with color= and thickness=.",
        )

    c1, c2, c3 = st.columns(3)
    with c1:
        agg = st.checkbox(
            "Merge duplicate pairs (undirected, keep strongest link)",
            value=True,
        )
    with c2:
        show_spokes = st.checkbox(
            "Show faint radial guides (center → label)",
            value=True,
            help="Light gray lines along the same direction as the label rotation.",
        )
    with c3:
        st.caption("Thickness uses a sqrt spread; tune chords below.")

    fig_w = st.slider("Figure width (px)", 480, 1000, 800, 20)
    fig_h = st.slider("Figure height (px)", 480, 1000, 800, 20)

    s1, s2, s3 = st.columns(3)
    with s1:
        thick_scale = st.slider(
            "Chord width multiplier",
            0.35,
            1.6,
            1.0,
            0.05,
            help="Scales mapped line width (after sqrt spread from file thickness).",
        )
    with s2:
        edge_opacity = st.slider(
            "Chord opacity",
            0.12,
            0.75,
            0.42,
            0.02,
            help="Lower = more see-through overlap; higher = bolder lines.",
        )
    with s3:
        label_font = st.slider("Genus label size", 7, 14, 10, 1)

    r4, r5 = st.columns(2)
    with r4:
        label_gap = st.slider(
            "Space from arc outer edge to label (× ring radius)",
            0.18,
            0.58,
            0.32,
            0.02,
            help="Same offset for every genus: label centre sits this far outside the coloured arc.",
        )
    with r5:
        arc_fill = st.slider(
            "Arc angular fill",
            0.62,
            0.92,
            0.78,
            0.02,
            help="How much of each angular slot the arc uses; lower leaves wider gaps (no overlap).",
        )

    if not links_file:
        st.info("Upload a **links** file to draw connections. The chromosome file sets node colours; if omitted, colours cycle automatically.")
        return

    links_bytes = links_file.getvalue()
    try:
        links_text = links_bytes.decode("utf-8")
    except UnicodeDecodeError:
        links_text = links_bytes.decode("latin-1")

    edges = parse_links_text(links_text)
    if not edges:
        st.error("No valid link lines found. Expected format: `GenusA 100 0 GenusB 100 0 color=…,thickness=…p`")
        return

    chromosome_colors = {}
    if chr_file:
        cbytes = chr_file.getvalue()
        try:
            ctext = cbytes.decode("utf-8")
        except UnicodeDecodeError:
            ctext = cbytes.decode("latin-1")
        chromosome_colors = parse_chromosome_text(ctext)
        if not chromosome_colors:
            st.warning("Chromosome file had no `chr - …` lines; using automatic colours for all nodes.")

    nodes_from_edges = set()
    for e in edges:
        nodes_from_edges.add(e["source"])
        nodes_from_edges.add(e["target"])

    ordered_nodes = list(chromosome_colors.keys())
    seen = set(ordered_nodes)
    for g in sorted(nodes_from_edges - seen):
        ordered_nodes.append(g)

    chromosome_colors = _merge_node_colors_in_order(chromosome_colors, ordered_nodes)

    st.caption(
        f"{len(ordered_nodes)} nodes, {len(edges)} link lines"
        + (" (aggregated)" if agg else "")
    )

    fig = build_circle_graph_figure(
        chromosome_colors,
        edges,
        aggregate_edges=agg,
        thickness_to_linewidth_scale=thick_scale,
        edge_opacity=edge_opacity,
        font_size=label_font,
        width=fig_w,
        height=fig_h,
        title=f"{links_file.name}",
        show_radial_guides=show_spokes,
        label_radial_gap=label_gap,
        arc_angular_fill=arc_fill,
    )
    st.pyplot(fig, clear_figure=True)
    plt.close(fig)

    if agg:
        n_pairs = len(aggregate_undirected_edges(edges))
        st.caption(f"After merging undirected pairs: {n_pairs} unique connections.")
