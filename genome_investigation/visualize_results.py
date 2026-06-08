"""
Paper-ready figures and interpretation report for genome investigation results.

Usage:
  python -m genome_investigation.visualize_results
  python -m genome_investigation.visualize_results --selected-species genome_investigation/selected_species.yaml
"""

from __future__ import annotations

import argparse
import textwrap
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from genome_investigation.io_utils import load_selected_species, normalize_bacid

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(__file__).resolve().parent / "results"
PAPER_DIR = RESULTS_DIR / "paper"
SELECTED_YAML = Path(__file__).resolve().parent / "selected_species.yaml"
SPECIES_DATA = REPO_ROOT / "species_data"

# Publication style
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.facecolor": "white",
    }
)

CATEGORY_COLORS = {
    "Metabolic generalists": "#4C78A8",
    "Metabolic specialists": "#F58518",
    "Resistance outliers": "#E45756",
}


def species_two_line_label(name: str) -> str:
    """Binomial name on two lines: genus on line 1, epithet(s) on line 2."""
    parts = str(name).strip().split()
    if len(parts) >= 2:
        return f"{parts[0]}\n{' '.join(parts[1:])}"
    return str(name)


def _two_line_ticklabels(ax: plt.Axes, axis: str, names: Sequence[str], *, fontsize: float = 8) -> None:
    labels = [species_two_line_label(n) for n in names]
    if axis == "y":
        ax.set_yticklabels(labels, fontsize=fontsize)
        plt.setp(ax.get_yticklabels(), multialignment="right")
    else:
        ax.set_xticklabels(labels, fontsize=fontsize)


def _bacdive_breadth_by_species(species_names: List[str]) -> pd.DataFrame:
    """Max breadth per activity for each species from local Step3 zips."""
    import zipfile as zf_mod

    files = {
        "utilization": SPECIES_DATA / "step3_met_util_exploded.csv.zip",
        "production": SPECIES_DATA / "step3_met_prod_exploded.csv.zip",
        "resistance": SPECIES_DATA / "step3_met_res_exploded.csv.zip",
    }
    meta = {"BacID", "species", "genus", "order", "type_strain", "is_strain", "species_with_id"}
    rows = []
    for activity, path in files.items():
        if not path.exists():
            continue
        with zf_mod.ZipFile(path, "r") as zf:
            csvs = [n for n in zf.namelist() if n.endswith(".csv")]
            df = pd.read_csv(zf.open(csvs[0]), low_memory=False)
        mets = [c for c in df.columns if c not in meta]
        sub = df[df["species"].isin(species_names)]
        if sub.empty:
            continue
        x = sub[mets].apply(pd.to_numeric, errors="coerce")
        if activity == "resistance":
            breadth = (x == 1).sum(axis=1)
        else:
            breadth = (x.replace(-1, 0).fillna(0) > 0).sum(axis=1)
        part = sub[["species"]].copy()
        part["activity"] = activity
        part["breadth"] = breadth.values
        rows.append(part)
    if not rows:
        return pd.DataFrame()
    long = pd.concat(rows, ignore_index=True)
    return long.groupby(["species", "activity"], as_index=False)["breadth"].max()


def load_category_map(yaml_path: Path) -> Dict[str, str]:
    cfg = load_selected_species(yaml_path)
    mapping: Dict[str, str] = {}
    for block in cfg.get("raw", {}).get("categories") or []:
        if isinstance(block, dict):
            name = block.get("name", "Other")
            for sp in block.get("species") or []:
                mapping[str(sp)] = name
    return mapping


def build_integrated_table(
    genome: pd.DataFrame,
    ranked: pd.DataFrame,
    category_map: Dict[str, str],
) -> pd.DataFrame:
    genome = genome.copy()
    genome["BacID_norm"] = genome["BacID"].map(normalize_bacid)
    ranked = ranked.copy()
    if not ranked.empty and "best_BacID" in ranked.columns:
        ranked["BacID_norm"] = ranked["best_BacID"].map(normalize_bacid)
        rank_best = ranked.sort_values("priority_score", ascending=False).drop_duplicates("species")
        rank_by_sp = rank_best.set_index("species")
    else:
        rank_by_sp = pd.DataFrame()

    breadth = _bacdive_breadth_by_species(genome["species"].astype(str).tolist())
    rows = []
    for _, g in genome.iterrows():
        sp = str(g["species"])
        row = {
            "Species": sp,
            "Research category": category_map.get(sp, "—"),
            "BacID": g.get("BacID"),
            "Strain": g.get("strain") or "—",
            "Genome accession": g.get("genome_accession") or "—",
            "Assembly level": g.get("assembly_level") or "—",
            "GC (%)": g.get("gc_percent"),
            "Match confidence": g.get("match_confidence"),
            "Data source": g.get("source_database") or "—",
        }
        if sp in rank_by_sp.index:
            r = rank_by_sp.loc[sp]
            row["Priority score"] = r.get("priority_score")
            row["Outlier category (auto)"] = r.get("category")
            row["BacDive z (auto rank)"] = r.get("bacdive_z")
        for act in ("utilization", "production", "resistance"):
            if not breadth.empty:
                sub = breadth[(breadth["species"] == sp) & (breadth["activity"] == act)]
                row[f"BacDive {act[:4]} breadth"] = int(sub["breadth"].iloc[0]) if len(sub) else "—"
            else:
                row[f"BacDive {act[:4]} breadth"] = "—"
        conf = float(g.get("match_confidence") or 0)
        if conf >= 0.95:
            row["Interpretation"] = "High-confidence strain-level genome link; suitable for download / antiSMASH."
        elif conf >= 0.7:
            row["Interpretation"] = "Species-level NCBI match; verify strain before mechanistic claims."
        else:
            row["Interpretation"] = "Weak or missing genome link; phenotype may reflect under-testing."
        rows.append(row)
    return pd.DataFrame(rows)


def plot_match_confidence(genome: pd.DataFrame, category_map: Dict[str, str], out: Path) -> None:
    """Standalone figure 1: strain-level (green) vs NCBI-only (red) grid."""
    d = genome.copy()
    if "Research category" not in d.columns and category_map:
        d["Research category"] = d["species"].map(category_map).fillna("Other")
    fig, ax = plt.subplots(figsize=(7, max(4, 0.45 * len(d))))
    plot_match_quality_grid(
        d,
        ax,
        sp_col="species" if "species" in d.columns else "Species",
        acc_col="genome_accession" if "genome_accession" in d.columns else "Genome accession",
        conf_col="match_confidence" if "match_confidence" in d.columns else "Match confidence",
        src_col="source_database" if "source_database" in d.columns else "Data source",
        cat_col="Research category" if "Research category" in d.columns else "category",
    )
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_phenotype_vs_confidence(genome: pd.DataFrame, category_map: Dict[str, str], out: Path) -> None:
    breadth = _bacdive_breadth_by_species(genome["species"].astype(str).tolist())
    if breadth.empty:
        return
    util = breadth[breadth["activity"] == "utilization"][["species", "breadth"]].rename(columns={"breadth": "util"})
    prod = breadth[breadth["activity"] == "production"][["species", "breadth"]].rename(columns={"breadth": "prod"})
    res = breadth[breadth["activity"] == "resistance"][["species", "breadth"]].rename(columns={"breadth": "res"})
    m = genome.merge(util, on="species", how="left").merge(prod, on="species", how="left").merge(res, on="species", how="left")
    m["category"] = m["species"].map(category_map)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for ax, col, title in zip(axes, ["util", "prod", "res"], ["Utilization", "Production", "Resistance"]):
        for cat, sub in m.groupby("category"):
            ax.scatter(
                sub["match_confidence"],
                sub[col],
                s=70,
                alpha=0.85,
                label=cat,
                c=CATEGORY_COLORS.get(cat, "#888"),
                edgecolors="white",
                linewidths=0.5,
            )
        ax.set_xlabel("Genome match confidence")
        ax.set_title(title)
        ax.set_xlim(0, 1.05)
    axes[0].set_ylabel("BacDive metabolite breadth")
    axes[1].legend(loc="best", fontsize=7, frameon=True)
    fig.suptitle("Phenotypic breadth vs. genome evidence strength", y=1.02)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def _col(integrated: pd.DataFrame, *names: str) -> str:
    """Resolve first matching column name."""
    lower = {c.lower(): c for c in integrated.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    raise KeyError(f"None of {names} in {list(integrated.columns)}")


def _is_strain_level_match(row: pd.Series, conf_col: str, src_col: str) -> bool:
    conf = float(pd.to_numeric(row.get(conf_col), errors="coerce") or 0)
    src = str(row.get(src_col, "")).lower()
    return conf >= 0.95 or "bacdive" in src


def plot_match_quality_grid(
    d: pd.DataFrame,
    ax: plt.Axes,
    *,
    sp_col: str,
    acc_col: str,
    conf_col: str,
    src_col: str,
    cat_col: str,
) -> None:
    """
    Grid: one row per species; green = strain-level (BacDive), red = NCBI species-level only.
    """
    d = d.copy()
    d["_strain"] = d.apply(lambda r: _is_strain_level_match(r, conf_col, src_col), axis=1)
    d = d.sort_values([cat_col, sp_col], ascending=[True, True]).reset_index(drop=True)

    n = len(d)
    # Matrix: columns = [strain-level, NCBI-only]; value 1 = active cell, 0 = inactive
    mat = np.zeros((n, 2))
    labels = np.empty((n, 2), dtype=object)
    for i, (_, row) in enumerate(d.iterrows()):
        acc = str(row[acc_col]) if pd.notna(row[acc_col]) else "—"
        if row["_strain"]:
            mat[i, 0] = 1
            mat[i, 1] = 0
            labels[i, 0] = acc[:18]
            labels[i, 1] = ""
        else:
            mat[i, 0] = 0
            mat[i, 1] = 2
            labels[i, 0] = ""
            labels[i, 1] = acc[:18]

    # 0=empty, 1=strain (green), 2=ncbi (red)
    cmap = ListedColormap(["#f5f5f5", "#2ca02c", "#d62728"])
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Strain-level\n(BacDive)", "NCBI species-level\n(fallback)"], fontsize=9)
    ax.set_yticks(np.arange(n))
    _two_line_ticklabels(ax, "y", d[sp_col].astype(str), fontsize=7)
    ax.set_title("A  Genome match quality")

    for i in range(n):
        for j in range(2):
            if mat[i, j] >= 1 and labels[i, j]:
                ax.text(
                    j,
                    i,
                    labels[i, j],
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="white",
                    fontweight="bold",
                )

    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.tick_params(which="minor", size=0)

    ax.legend(
        handles=[
            Patch(facecolor="#2ca02c", edgecolor="white", label="Strain-level (BacDive accession)"),
            Patch(facecolor="#d62728", edgecolor="white", label="NCBI species-level only"),
            Patch(facecolor="#f5f5f5", edgecolor="#ccc", label="Not applicable"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=3,
        fontsize=7,
        frameon=True,
    )


def plot_breadth_grouped_bars(
    d: pd.DataFrame,
    ax: plt.Axes,
    *,
    sp_col: str,
    cat_col: str,
) -> None:
    """Horizontal grouped bars: one row per species, util / prod / resistance breadth."""
    n = len(d)
    y = np.arange(n)
    h = 0.24

    ax.barh(y + h, d["_util"], height=h, label="Utilization", color="#4C78A8", edgecolor="white", linewidth=0.4)
    ax.barh(y, d["_prod"], height=h, label="Production", color="#F58518", edgecolor="white", linewidth=0.4)
    ax.barh(y - h, d["_res"], height=h, label="Resistance", color="#E45756", edgecolor="white", linewidth=0.4)

    ax.set_yticks(y)
    _two_line_ticklabels(ax, "y", d[sp_col].astype(str), fontsize=7)
    for tick_label, (_, row) in zip(ax.get_yticklabels(), d.iterrows()):
        tick_label.set_color(CATEGORY_COLORS.get(str(row[cat_col]), "#333333"))

    # Category group dividers
    prev_cat = None
    for yi, (_, row) in enumerate(d.iterrows()):
        cat = str(row[cat_col])
        if prev_cat is not None and cat != prev_cat:
            ax.axhline(yi - 0.5, color="#cccccc", linewidth=1.0, linestyle="--")
        prev_cat = cat

    xmax = max(float(d[["_util", "_prod", "_res"]].max().max()), 1.0) * 1.12
    for yi, (_, row) in enumerate(d.iterrows()):
        for val, yoff in ((row["_util"], h), (row["_prod"], 0), (row["_res"], -h)):
            if val > 0:
                ax.text(
                    float(val) + xmax * 0.02,
                    yi + yoff,
                    str(int(val)),
                    va="center",
                    ha="left",
                    fontsize=6.5,
                    color="#333333",
                )

    ax.set_xlim(0, xmax)
    ax.invert_yaxis()
    ax.set_xlabel("# BacDive metabolites (breadth)")
    ax.set_title("B  BacDive phenotypic breadth")
    ax.legend(loc="lower right", fontsize=7, frameon=True, title="Activity")


def _scatter_label_clusters(d: pd.DataFrame, *, y_gap: float = 4.0) -> List[pd.DataFrame]:
    """Cluster points that share x (rounded) and nearly the same y for label stacking."""
    clusters: List[pd.DataFrame] = []
    for _, xgrp in d.sort_values(["_conf", "_total_breadth"]).groupby(d["_conf"].round(2)):
        pending: List[pd.DataFrame] = []
        for _, row in xgrp.iterrows():
            row_df = row.to_frame().T
            if not pending:
                pending.append(row_df)
                continue
            last_y = float(pending[-1]["_total_breadth"].iloc[-1])
            if abs(float(row["_total_breadth"]) - last_y) <= y_gap:
                pending[-1] = pd.concat([pending[-1], row_df], ignore_index=True)
            else:
                clusters.append(pending.pop())
                pending.append(row_df)
        for chunk in pending:
            clusters.append(chunk)
    return clusters


def plot_scatter_phenotype_vs_confidence(
    d: pd.DataFrame,
    ax: plt.Axes,
    *,
    sp_col: str,
    cat_col: str,
) -> None:
    """Scatter with non-overlapping two-line species labels."""
    for cat, sub in d.groupby(cat_col):
        ax.scatter(
            sub["_conf"],
            sub["_total_breadth"],
            s=40 + 8 * sub["_res"],
            c=CATEGORY_COLORS.get(str(cat), "#888"),
            label=str(cat),
            edgecolors="white",
            linewidths=0.6,
            alpha=0.9,
            zorder=2,
        )

    y_span = max(float(d["_total_breadth"].max()), 1.0)
    y_stack = max(3.0, y_span * 0.035)

    for grp in _scatter_label_clusters(d.copy()):
        grp = grp.sort_values("_total_breadth").reset_index(drop=True)
        n = len(grp)
        for i, (_, row) in enumerate(grp.iterrows()):
            y_label = float(row["_total_breadth"]) + (i - (n - 1) / 2) * y_stack
            x_label = float(row["_conf"]) + 0.028
            ax.annotate(
                species_two_line_label(row[sp_col]),
                (row["_conf"], row["_total_breadth"]),
                xytext=(x_label, y_label),
                textcoords="data",
                fontsize=6,
                ha="left",
                va="center",
                color=CATEGORY_COLORS.get(str(row[cat_col]), "#333333"),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#dddddd", linewidth=0.4, alpha=0.9),
                arrowprops=dict(arrowstyle="-", color="#bbbbbb", lw=0.5, shrinkA=3, shrinkB=3),
                zorder=3,
            )

    ax.axvline(0.7, color="#999", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Match confidence")
    ax.set_ylabel("Total BacDive breadth (util+prod+res)")
    ax.set_title("D  Phenotype load vs. genome evidence (size ∝ resistance breadth)")
    ax.legend(fontsize=6, loc="upper left", frameon=True)


def plot_composite_dashboard(integrated: pd.DataFrame, out: Path) -> None:
    """
    Multi-panel figure synthesizing Table 1 (integrated summary).
    Panels A (match grid), B (breadth bars), C (heatmap), D (phenotype vs confidence).
    """
    d = integrated.copy()
    sp_col = _col(d, "Species", "species")
    cat_col = _col(d, "Research category", "category")
    conf_col = _col(d, "Match confidence", "match_confidence")
    acc_col = _col(d, "Genome accession", "genome_accession")
    src_col = _col(d, "Data source", "source_database")
    util_col = _col(d, "BacDive util breadth", "util")
    prod_col = _col(d, "BacDive prod breadth", "prod")
    res_col = _col(d, "BacDive resi breadth", "BacDive res breadth", "res")

    d["_cat"] = d[cat_col].map(CATEGORY_COLORS).fillna("#888888")
    d["_conf"] = pd.to_numeric(d[conf_col], errors="coerce").fillna(0)
    d["_util"] = pd.to_numeric(d[util_col], errors="coerce").fillna(0)
    d["_prod"] = pd.to_numeric(d[prod_col], errors="coerce").fillna(0)
    d["_res"] = pd.to_numeric(d[res_col], errors="coerce").fillna(0)
    d["_total_breadth"] = d["_util"] + d["_prod"] + d["_res"]
    d = d.sort_values([cat_col, sp_col], ascending=[True, True]).reset_index(drop=True)

    n_species = len(d)
    fig_h = max(8.5, 0.62 * n_species + 3.5)
    fig = plt.figure(figsize=(12, fig_h))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.2, 1.0],
        width_ratios=[1.0, 1.05],
        hspace=0.38,
        wspace=0.38,
    )

    ax0 = fig.add_subplot(gs[0, 0])
    plot_match_quality_grid(d, ax0, sp_col=sp_col, acc_col=acc_col, conf_col=conf_col, src_col=src_col, cat_col=cat_col)

    ax1 = fig.add_subplot(gs[0, 1])
    plot_breadth_grouped_bars(d, ax1, sp_col=sp_col, cat_col=cat_col)

    ax2 = fig.add_subplot(gs[1, 0])
    heat = d[[sp_col, util_col, prod_col, res_col]].set_index(sp_col)
    heat.columns = ["Utilization", "Production", "Resistance"]
    heat = heat.apply(pd.to_numeric, errors="coerce").fillna(0)
    heat.index = [species_two_line_label(str(s)) for s in heat.index]
    sns.heatmap(heat, annot=True, fmt=".0f", cmap="YlGnBu", ax=ax2, cbar_kws={"shrink": 0.75})
    ax2.set_title("C  Breadth heatmap")
    ax2.set_ylabel("")
    plt.setp(ax2.get_yticklabels(), fontsize=7, multialignment="right")

    ax3 = fig.add_subplot(gs[1, 1])
    plot_scatter_phenotype_vs_confidence(d, ax3, sp_col=sp_col, cat_col=cat_col)

    fig.suptitle("Integrated genome & BacDive phenotype summary (focal species)", fontsize=13, y=0.98)
    fig.subplots_adjust(left=0.14, right=0.96)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_breadth_heatmap(integrated: pd.DataFrame, out: Path) -> None:
    cols = [c for c in integrated.columns if c.startswith("BacDive")]
    if not cols:
        return
    d = integrated.set_index("Species")[cols].apply(pd.to_numeric, errors="coerce")
    d.columns = [c.replace("BacDive ", "").replace(" breadth", "") for c in d.columns]
    d.index = [species_two_line_label(str(i)) for i in d.index]
    fig, ax = plt.subplots(figsize=(5.5, max(3, 0.5 * len(d))))
    sns.heatmap(d, annot=True, fmt=".0f", cmap="YlGnBu", ax=ax, cbar_kws={"label": "# metabolites"})
    ax.set_title("BacDive phenotypic breadth (selected species)")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def write_interpretation_md(
    genome: pd.DataFrame,
    integrated: pd.DataFrame,
    failures: pd.DataFrame,
    out: Path,
) -> None:
    n = len(genome)
    n_high = int((genome["match_confidence"].astype(float) >= 0.95).sum())
    n_mid = int(((genome["match_confidence"].astype(float) >= 0.7) & (genome["match_confidence"].astype(float) < 0.95)).sum())
    n_bacdive = int(genome["source_database"].astype(str).str.contains("bacdive", case=False, na=False).sum())

    lines = [
        "# Genome investigation — results interpretation",
        "",
        "This report summarizes the **optional genome evidence layer** for nine focal species ",
        "(metabolic generalists, specialists, and resistance outliers). ",
        "Phenotypes come from local BacDive Step3 matrices; genome accessions were resolved via ",
        "**BacDive API v2** (strain-level) with **NCBI Datasets** fallback (species-level).",
        "",
        "## Key findings",
        "",
        f"- **{n_high}/{n}** strains have *high-confidence* genome links (score ≥ 0.95; BacDive strain record with INSDC accession).",
        f"- **{n_mid}/{n}** strains rely on species-level NCBI assemblies (score ≈ 0.72); strain designation was not verified in assembly metadata.",
        f"- **{n_bacdive}/{n}** records include direct BacDive genome metadata.",
        "",
        "## How to read match confidence",
        "",
        "| Score | Meaning | Recommended use |",
        "|-------|---------|-----------------|",
        "| **1.0** | BacDive strain lists INSDC accession (e.g. GCA_…) | Download genome; run antiSMASH; cite accession |",
        "| **0.7–0.9** | NCBI species match; strain uncertain | Compare phenotype to reference assembly; note limitation in text |",
        "| **< 0.7** | Weak / missing link | Do not infer genomic mechanisms; phenotype may reflect sparse testing |",
        "",
        "## Per-species summary",
        "",
    ]
    for _, r in integrated.iterrows():
        lines.append(f"### {r['Species']} ({r['Research category']})")
        lines.append(
            f"- **Genome:** `{r['Genome accession']}` ({r['Assembly level']}) — confidence **{r['Match confidence']}** via {r['Data source']}."
        )
        lines.append(f"- **BacDive breadth:** util {r.get('BacDive util breadth', '—')}, prod {r.get('BacDive prod breadth', '—')}, res {r.get('BacDive res breadth', '—')}.")
        lines.append(f"- **Interpretation:** {r['Interpretation']}")
        lines.append("")

    if not failures.empty:
        lines.extend(
            [
                "## Strains logged as low-confidence",
                "",
                "These rows are listed in `logs/genome_lookup_failures.csv` (species-level match without verified strain):",
                "",
            ]
        )
        for _, f in failures.iterrows():
            lines.append(f"- *{f.get('species')}* (BacID {f.get('BacID')}): {f.get('match_notes')}")
        lines.append("")

    lines.extend(
        [
            "## Figures (paper-ready)",
            "",
            "| File | Description |",
            "|------|-------------|",
            "| `paper/fig1_genome_match_confidence.png` | Match confidence by species, colored by research category |",
            "| `paper/fig2_phenotype_vs_genome_confidence.png` | BacDive breadth vs. genome confidence (util/prod/res) |",
            "| `paper/fig3_bacdive_breadth_heatmap.png` | Heatmap of phenotypic breadth |",
            "| `paper/fig4_composite_integrated_dashboard.png` | **Multi-panel summary** (A: match grid; B–C, D: phenotype) |",
            "| `paper/table1_integrated_summary.csv` | Main supplementary table |",
            "| `paper/table1_integrated_summary.md` | Markdown version of Table 1 |",
            "",
            "## Suggested text for Methods (snippet)",
            "",
            "> Genome accessions were retrieved programmatically from BacDive API v2 by BacDive ID, ",
            "> supplemented by NCBI Datasets taxon queries when strain-level links were absent. ",
            "> Match confidence scores reflect strain-level agreement (1.0) or species-level fallback (≈0.72). ",
            "> antiSMASH was not run by default; biosynthetic gene cluster analysis is reserved for ",
            "> a manually selected subset with high-confidence assemblies.",
            "",
        ]
    )
    out.write_text("\n".join(lines), encoding="utf-8")


def dataframe_to_markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        cells = []
        for c in cols:
            v = row[c]
            if pd.isna(v):
                cells.append("—")
            else:
                s = str(v).replace("|", "\\|")
                cells.append(s)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def generate_paper_outputs(
    *,
    genome_path: Path,
    ranked_path: Path,
    failures_path: Path,
    selected_yaml: Path,
    out_dir: Path,
) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    genome = pd.read_csv(genome_path)
    ranked = pd.read_csv(ranked_path) if ranked_path.exists() else pd.DataFrame()
    failures = pd.read_csv(failures_path) if failures_path.exists() else pd.DataFrame()
    category_map = load_category_map(selected_yaml)

    integrated = build_integrated_table(genome, ranked, category_map)

    paths = {}
    t_csv = out_dir / "table1_integrated_summary.csv"
    integrated.to_csv(t_csv, index=False)
    paths["table_csv"] = t_csv

    t_md = out_dir / "table1_integrated_summary.md"
    t_md.write_text("# Table 1. Integrated BacDive phenotype and genome metadata\n\n" + dataframe_to_markdown_table(integrated) + "\n", encoding="utf-8")
    paths["table_md"] = t_md

    f1 = out_dir / "fig1_genome_match_confidence.png"
    plot_match_confidence(genome, category_map, f1)
    paths["fig1"] = f1

    f2 = out_dir / "fig2_phenotype_vs_genome_confidence.png"
    plot_phenotype_vs_confidence(genome, category_map, f2)
    paths["fig2"] = f2

    f3 = out_dir / "fig3_bacdive_breadth_heatmap.png"
    plot_breadth_heatmap(integrated, f3)
    paths["fig3"] = f3

    f4 = out_dir / "fig4_composite_integrated_dashboard.png"
    plot_composite_dashboard(integrated, f4)
    paths["fig4_composite"] = f4

    interp = out_dir / "RESULTS_INTERPRETATION.md"
    write_interpretation_md(genome, integrated, failures, interp)
    paths["interpretation"] = interp

    return paths


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Generate paper-ready genome investigation figures")
    parser.add_argument("--genome-enriched", default=str(RESULTS_DIR / "Step2_5_genome_enriched.csv"))
    parser.add_argument("--ranked", default=str(RESULTS_DIR / "ranked_species_candidates.csv"))
    parser.add_argument("--failures", default=str(REPO_ROOT / "logs" / "genome_lookup_failures.csv"))
    parser.add_argument("--selected-species", default=str(SELECTED_YAML))
    parser.add_argument("--output-dir", default=str(PAPER_DIR))
    parser.add_argument(
        "--composite-only",
        action="store_true",
        help="Regenerate only fig4 from existing table1_integrated_summary.csv",
    )
    args = parser.parse_args(argv)

    if args.composite_only:
        table_path = Path(args.output_dir) / "table1_integrated_summary.csv"
        if not table_path.exists():
            print(f"[error] missing {table_path}; run full visualize_results first")
            return 1
        out = Path(args.output_dir) / "fig4_composite_integrated_dashboard.png"
        plot_composite_dashboard(pd.read_csv(table_path), out)
        print(f"[done] wrote {out}")
        return 0

    paths = generate_paper_outputs(
        genome_path=Path(args.genome_enriched),
        ranked_path=Path(args.ranked),
        failures_path=Path(args.failures),
        selected_yaml=Path(args.selected_species),
        out_dir=Path(args.output_dir),
    )
    print(f"[done] paper-ready outputs in {args.output_dir}:")
    for k, p in paths.items():
        print(f"  {k}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
