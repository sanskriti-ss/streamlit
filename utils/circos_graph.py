"""
Parse Circos-style chromosome / link text files and build an interactive circle (chord-style) graph.
"""
from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Wedge

# Circos-style names used in chromosome files and link color= attributes → hex
CIRCOS_COLOR_TO_HEX: Dict[str, str] = {
    "pred": "#E41A1C",
    "orange": "#FF7F00",
    "yellow": "#FFFF33",
    "green": "#4DAF4A",
    "blue": "#377EB8",
    "purple": "#984EA3",
    "grey": "#999999",
    "gray": "#999999",
    "porange": "#FDB462",
    "pyellow": "#FFF7BC",
    "pgreen": "#B3E2CD",
    "pblue": "#B3CDE3",
    "ppurple": "#DECBE4",
    "vvlorange": "#FFF5EB",
    "vlorange": "#FFE0C2",
    "lorange": "#FFC080",
    "dorange": "#FF8C00",
    "vdorange": "#CC7000",
    "vvdorange": "#994400",
    "vvlred": "#FFE8E8",
    "vlred": "#FFC0C0",
    "lred": "#FF8080",
    "red": "#E41A1C",
    "dred": "#B20000",
    "vdred": "#800000",
    "vvdred": "#500000",
    "black": "#222222",
    "white": "#FFFFFF",
}


def resolve_color(name: Optional[str], fallback: str = "#888888") -> str:
    if not name:
        return fallback
    s = name.strip()
    if s.startswith("#") and len(s) in (4, 7, 9):
        return s
    return CIRCOS_COLOR_TO_HEX.get(s.lower(), fallback)


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    """rgba() string (e.g. for any remaining string-based APIs)."""
    h = hex_color.lstrip("#")
    if len(h) == 6:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"
    return hex_color


def _hex_to_rgba_tuple(hex_color: str, alpha: float) -> Tuple[float, float, float, float]:
    h = hex_color.lstrip("#")
    if len(h) == 6:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        return (r / 255.0, g / 255.0, b / 255.0, alpha)
    return (0.5, 0.5, 0.5, alpha)


def _radial_label_rotation_deg(lx: float, ly: float) -> float:
    """
    Degrees CCW from +x so the label baseline lies on the ray from the origin
    through (lx, ly). Matplotlib rotates around the anchor when rotation_mode='anchor'.
    """
    return float(math.degrees(math.atan2(ly, lx)))


_LINK_LINE_RE = re.compile(
    r"^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(.+)$"
)


def parse_chromosome_text(text: str) -> Dict[str, str]:
    """
    Parse lines like: chr - Pseudomonas Pseudomonas 0 100 pred
    Returns genus_name -> circos color name (for resolve_color).
    """
    out: Dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 7 or parts[0] != "chr" or parts[1] != "-":
            continue
        genus = parts[2]
        color_name = parts[-1]
        out[genus] = color_name
    return out


def parse_links_text(text: str) -> List[Dict[str, Any]]:
    """
    Parse lines like:
    Acinetobacter 100 0 Aeromonas 100 0 color=orange,thickness=15p
    """
    edges: List[Dict[str, Any]] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        m = _LINK_LINE_RE.match(line)
        if not m:
            continue
        g1, g2, tail = m.group(1), m.group(4), m.group(7)
        color_name = "black"
        thickness = 4
        if "color=" in tail:
            cm = re.search(r"color=([^,\s]+)", tail)
            if cm:
                color_name = cm.group(1).strip()
        if "thickness=" in tail:
            tm = re.search(r"thickness=(\d+)p", tail)
            if tm:
                thickness = int(tm.group(1))
        edges.append(
            {
                "source": g1,
                "target": g2,
                "color_name": color_name,
                "thickness": thickness,
            }
        )
    return edges


def _node_positions(
    node_list: List[str], radius: float = 1.0, start_angle: float = np.pi / 2
) -> Dict[str, Tuple[float, float, float]]:
    """Map node -> (x, y, theta). Even spacing starting from top."""
    n = len(node_list)
    if n == 0:
        return {}
    positions: Dict[str, Tuple[float, float, float]] = {}
    for i, name in enumerate(node_list):
        theta = start_angle + 2 * np.pi * i / n
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        positions[name] = (float(x), float(y), float(theta))
    return positions


def _quadratic_bezier(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, n: int = 48
) -> Tuple[np.ndarray, np.ndarray]:
    t = np.linspace(0.0, 1.0, n)
    one_m = 1.0 - t
    xy = (
        np.outer(one_m**2, p0)
        + np.outer(2 * one_m * t, p1)
        + np.outer(t**2, p2)
    )
    return xy[:, 0], xy[:, 1]


def _edge_curve_points(
    theta_a: float,
    theta_b: float,
    radius: float,
    bulge: float = 0.35,
) -> Tuple[np.ndarray, np.ndarray]:
    p0 = np.array([radius * np.cos(theta_a), radius * np.sin(theta_a)])
    p2 = np.array([radius * np.cos(theta_b), radius * np.sin(theta_b)])
    mid = (p0 + p2) / 2.0
    p1 = mid * (1.0 - bulge)
    return _quadratic_bezier(p0, p1, p2)


def _canonical_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def aggregate_undirected_edges(
    edges: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Merge duplicate undirected pairs; keep max thickness and last color (prefer stronger)."""
    best: Dict[Tuple[str, str], Dict[str, Any]] = {}
    for e in edges:
        ca = _canonical_pair(e["source"], e["target"])
        t = e["thickness"]
        if ca not in best or t > best[ca]["thickness"]:
            best[ca] = {
                "source": ca[0],
                "target": ca[1],
                "color_name": e["color_name"],
                "thickness": t,
            }
    return list(best.values())


def _map_thickness_to_linewidth(
    t: float,
    t_min: float,
    t_max: float,
    *,
    min_linewidth: float,
    max_linewidth: float,
    scale: float,
) -> float:
    """Spread thickness values using sqrt so small differences stay visible."""
    if t_max <= t_min:
        lw = (min_linewidth + max_linewidth) / 2.0
    else:
        st = math.sqrt(max(t, 0.0))
        smin = math.sqrt(t_min)
        smax = math.sqrt(t_max)
        u = (st - smin) / (smax - smin)
        lw = min_linewidth + u * (max_linewidth - min_linewidth)
    lw = lw * scale
    return max(min_linewidth * 0.5, min(max_linewidth * scale, lw))


def build_circle_graph_figure(
    chromosome_colors: Dict[str, str],
    edges: List[Dict[str, Any]],
    *,
    node_radius: float = 1.0,
    edge_bulge: float = 0.38,
    thickness_to_linewidth_scale: float = 1.0,
    max_linewidth: float = 8.0,
    min_linewidth: float = 0.4,
    edge_opacity: float = 0.42,
    font_size: int = 10,
    height: int = 720,
    width: int = 720,
    aggregate_edges: bool = True,
    title: Optional[str] = None,
    show_radial_guides: bool = True,
    radial_guide_alpha: float = 0.18,
    # Ring: coloured arcs (annular wedges), non-overlapping in angle
    arc_band_rel: float = 0.15,
    arc_angular_fill: float = 0.78,
    # Labels: constant radial gap from outer arc edge to label centre (same for every genus)
    label_radial_gap: float = 0.32,
) -> Figure:
    """
    Chromosome file gives arc colours; edges list gives links with color + thickness.
    Each genus is a **curved arc** (annular wedge) on the ring—no dot markers, no black rims.
    Arcs are spaced with angular gaps so they do not overlap. Chords attach at the outer ring.

    Each label sits at ``ring_outer + label_radial_gap`` on the same ray as its arc centre,
    so every genus name is the same radial distance from its arc’s outer edge.
    """
    if aggregate_edges:
        edges = aggregate_undirected_edges(edges)

    nodes_from_edges = set()
    for e in edges:
        nodes_from_edges.add(e["source"])
        nodes_from_edges.add(e["target"])

    ordered: List[str] = []
    seen = set()
    for g in chromosome_colors.keys():
        if g not in seen:
            ordered.append(g)
            seen.add(g)
    for g in sorted(nodes_from_edges - seen):
        ordered.append(g)
        seen.add(g)

    dpi = 100
    fig_w_in = width / dpi
    fig_h_in = height / dpi
    fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in), dpi=dpi)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if not ordered:
        ax.text(0.5, 0.5, title or "No nodes to display", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        fig.tight_layout()
        return fig

    pos = _node_positions(ordered, radius=node_radius)

    drawn: List[Dict[str, Any]] = []
    for e in edges:
        src, tgt = e["source"], e["target"]
        if src not in pos or tgt not in pos:
            continue
        drawn.append(e)

    thickness_vals = [float(e["thickness"]) for e in drawn]
    t_min = min(thickness_vals) if thickness_vals else 0.0
    t_max = max(thickness_vals) if thickness_vals else 1.0

    groups: Dict[Tuple[str, float], List[Tuple[str, str]]] = {}
    for e in drawn:
        src, tgt = e["source"], e["target"]
        lw = _map_thickness_to_linewidth(
            float(e["thickness"]),
            t_min,
            t_max,
            min_linewidth=min_linewidth,
            max_linewidth=max_linewidth,
            scale=thickness_to_linewidth_scale,
        )
        hex_c = resolve_color(e.get("color_name"))
        key = (hex_c, round(lw, 3))
        groups.setdefault(key, []).append((src, tgt))

    sorted_keys = sorted(groups.keys(), key=lambda k: k[1], reverse=True)

    n = len(ordered)
    ring_outer_r = float(node_radius)
    ring_inner_r = ring_outer_r * max(0.02, min(0.95, 1.0 - arc_band_rel))
    slot = 2.0 * math.pi / n
    arc_span = slot * max(0.1, min(1.0, arc_angular_fill))
    half = arc_span / 2.0

    # Faint spokes: center → label point (same ray as arc centre)
    label_center_r = ring_outer_r + float(label_radial_gap)
    if show_radial_guides:
        for g in ordered:
            x0, y0, th = pos[g]
            lx = label_center_r * math.cos(th)
            ly = label_center_r * math.sin(th)
            ax.plot(
                [0.0, lx],
                [0.0, ly],
                color=(0.55, 0.55, 0.55, radial_guide_alpha),
                linewidth=0.6,
                zorder=0,
                clip_on=False,
            )

    for hex_c, lw in sorted_keys:
        pairs = groups[(hex_c, lw)]
        rgba = _hex_to_rgba_tuple(hex_c, edge_opacity)
        for src, tgt in pairs:
            _, _, ta = pos[src]
            _, _, tb = pos[tgt]
            xc, yc = _edge_curve_points(ta, tb, ring_outer_r, bulge=edge_bulge)
            ax.plot(
                xc,
                yc,
                color=rgba,
                linewidth=lw,
                solid_capstyle="round",
                zorder=1,
            )

    # Coloured arc segments (annular wedges), gaps between neighbours → no overlap
    for g in ordered:
        _, _, th = pos[g]
        t1 = th - half
        t2 = th + half
        d1 = float(np.degrees(t1))
        d2 = float(np.degrees(t2))
        col = resolve_color(chromosome_colors.get(g))
        wedge = Wedge(
            (0.0, 0.0),
            ring_outer_r,
            d1,
            d2,
            width=ring_outer_r - ring_inner_r,
            facecolor=col,
            edgecolor="none",
            linewidth=0,
            zorder=2,
        )
        ax.add_patch(wedge)

    lx = [label_center_r * math.cos(pos[g][2]) for g in ordered]
    ly = [label_center_r * math.sin(pos[g][2]) for g in ordered]
    for name, x, y in zip(ordered, lx, ly):
        rot = _radial_label_rotation_deg(x, y)
        ax.text(
            x,
            y,
            name,
            rotation=rot,
            rotation_mode="anchor",
            ha="center",
            va="center",
            fontsize=font_size,
            color="#1a1a1a",
            zorder=3,
            clip_on=False,
        )

    lim = label_center_r * 1.38
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    if title:
        fig.suptitle(title, fontsize=12, y=0.98)

    fig.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.02)
    return fig
