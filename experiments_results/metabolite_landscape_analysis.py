"""
Metabolite landscape + outlier analysis (local, reproducible).

Generates plots + a short report for:
3.1 Landscape of Bacterial Metabolite Profiles
3.2 Outlier Species (generalists, specialists, resistance outliers)

Data inputs (already in this repo):
- ./species_data/step3_met_{util,prod,res,sen}_exploded.csv.zip (or .csv)

Optional (for ecological context + genome size correlations):
- ./experiments_results/metadata_species.csv
  Must include a join key column named `species_with_id` (recommended) or (`BacID`,`species`).
  Suggested columns: isolation_source, habitat_complexity, genome_size_mb

Run:
    python experiments_results/metabolite_landscape_analysis.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import textwrap
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

try:
    from bokeh.io import output_file, save
    from bokeh.layouts import column, row
    from bokeh.models import (
        BasicTicker,
        ColorBar,
        ColumnDataSource,
        Div,
        HoverTool,
        LinearColorMapper,
    )
    from bokeh.palettes import Category20, RdBu11
    from bokeh.plotting import figure

    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False

REPO_ROOT = Path(__file__).resolve().parents[1]
SPECIES_DATA = REPO_ROOT / "species_data"

ACTIVITY_FILES = {
    "utilization": SPECIES_DATA / "step3_met_util_exploded.csv.zip",
    "production": SPECIES_DATA / "step3_met_prod_exploded.csv.zip",
    "resistance": SPECIES_DATA / "step3_met_res_exploded.csv.zip",
    "sensitivity": SPECIES_DATA / "step3_met_sen_exploded.csv.zip",
}

METADATA_COLS = ["BacID", "species", "genus", "order", "type_strain", "is_strain", "species_with_id"]


def _now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def _ensure_outdir(base: Path) -> Path:
    out = base / f"outputs_{_now_stamp()}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _read_csv_or_zip(path: Path, *, low_memory: bool = False) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(str(path))
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path, "r") as zf:
            csvs = [n for n in zf.namelist() if n.lower().endswith(".csv")]
            if not csvs:
                raise ValueError(f"No CSV found inside zip: {path}")
            with zf.open(csvs[0]) as f:
                return pd.read_csv(f, low_memory=low_memory)
    return pd.read_csv(path, low_memory=low_memory)


def _get_metabolite_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c not in METADATA_COLS]


def _binarize_activity(df: pd.DataFrame, activity: str) -> pd.DataFrame:
    """
    Return a metabolite matrix as 0/1 for the given activity.

    - utilization/production: treat -1 as 0; any positive -> 1
    - resistance/sensitivity: keep only 1 as positive, treat -1/0/NaN as 0
    """
    mets = _get_metabolite_columns(df)
    x = df[mets].copy()
    x = x.apply(pd.to_numeric, errors="coerce")

    if activity in ("resistance", "sensitivity"):
        out = (x == 1).astype(np.int8)
    else:
        out = x.replace(-1, 0)
        out = (out.fillna(0) > 0).astype(np.int8)
    out.columns = mets
    return out


def _species_key(df: pd.DataFrame) -> pd.Series:
    if "species_with_id" in df.columns:
        return df["species_with_id"].astype(str)
    if "BacID" in df.columns and "species" in df.columns:
        return df["BacID"].astype(str) + " | " + df["species"].astype(str)
    if "species" in df.columns:
        return df["species"].astype(str)
    return pd.Series([f"row_{i}" for i in range(len(df))])


@dataclass
class ActivityData:
    name: str
    df: pd.DataFrame
    mets_bin: pd.DataFrame
    key: pd.Series

    @property
    def breadth(self) -> pd.Series:
        return self.mets_bin.sum(axis=1)


def load_activity_data() -> Dict[str, ActivityData]:
    out: Dict[str, ActivityData] = {}
    for activity, path in ACTIVITY_FILES.items():
        df = _read_csv_or_zip(path, low_memory=False)
        key = _species_key(df)
        mets_bin = _binarize_activity(df, activity)
        out[activity] = ActivityData(activity, df, mets_bin, key)
    return out


def _lookup_phylum_gbif(genus: str, timeout_s: int = 10) -> Tuple[str, str]:
    """
    Query GBIF species match endpoint for a genus -> phylum mapping.
    Returns (phylum, source_note).
    """
    q = urlparse.urlencode({"name": genus, "rank": "GENUS"})
    url = f"https://api.gbif.org/v1/species/match?{q}"
    req = urlrequest.Request(
        url,
        headers={"User-Agent": "streamlit-metabolite-landscape/1.0 (+local-analysis-script)"},
    )
    with urlrequest.urlopen(req, timeout=timeout_s) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    phylum = payload.get("phylum")
    if isinstance(phylum, str) and phylum.strip():
        return phylum.strip(), "gbif"
    return "Unknown", "gbif_no_phylum"


def build_or_load_phylum_cache(
    activities: Dict[str, ActivityData],
    cache_csv: Path,
    *,
    force_refresh: bool = False,
    max_new_queries: Optional[int] = None,
) -> pd.DataFrame:
    """
    Create/reuse a genus->phylum cache CSV.
    - Reuses existing rows if present.
    - Queries missing genera from GBIF and appends them.
    """
    genera = set()
    for ad in activities.values():
        if "genus" in ad.df.columns:
            vals = ad.df["genus"].dropna().astype(str).str.strip()
            genera.update(v for v in vals if v and v.lower() != "nan")
    genera = sorted(genera)

    if cache_csv.exists() and not force_refresh:
        cache = pd.read_csv(cache_csv)
    else:
        cache = pd.DataFrame(columns=["genus", "phylum", "source", "updated_at"])

    existing = {}
    if not cache.empty and "genus" in cache.columns and "phylum" in cache.columns:
        for _, r in cache.iterrows():
            g = str(r["genus"]).strip()
            p = str(r["phylum"]).strip() if pd.notna(r["phylum"]) else "Unknown"
            if g:
                existing[g] = p or "Unknown"

    missing = [g for g in genera if g not in existing or force_refresh]
    if max_new_queries is not None:
        missing = missing[: max(0, int(max_new_queries))]

    new_rows = []
    for i, g in enumerate(missing, start=1):
        phylum = "Unknown"
        source = "lookup_failed"
        # light retry loop for transient HTTP/network issues
        for attempt in range(3):
            try:
                phylum, source = _lookup_phylum_gbif(g)
                break
            except (urlerror.URLError, TimeoutError, ValueError):
                time.sleep(0.6 * (attempt + 1))
        new_rows.append(
            {
                "genus": g,
                "phylum": phylum,
                "source": source,
                "updated_at": dt.datetime.now().isoformat(timespec="seconds"),
            }
        )
        # polite pacing to public API
        if i % 10 == 0:
            time.sleep(0.5)

    if new_rows:
        cache = pd.concat([cache, pd.DataFrame(new_rows)], ignore_index=True)
        cache = cache.drop_duplicates(subset=["genus"], keep="last").sort_values("genus")
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(cache_csv, index=False)
    elif not cache_csv.exists():
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        cache.to_csv(cache_csv, index=False)

    return cache


def attach_phylum_to_activities(activities: Dict[str, ActivityData], phylum_cache: pd.DataFrame) -> None:
    if phylum_cache.empty:
        return
    g2p = {
        str(r["genus"]).strip(): str(r["phylum"]).strip()
        for _, r in phylum_cache.iterrows()
        if pd.notna(r.get("genus"))
    }
    for ad in activities.values():
        if "genus" not in ad.df.columns:
            continue
        ad.df["phylum"] = (
            ad.df["genus"].fillna("Unknown").astype(str).str.strip().map(g2p).fillna("Unknown")
        )


def try_load_metadata(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    md = pd.read_csv(path)
    if "species_with_id" not in md.columns:
        if "BacID" in md.columns and "species" in md.columns:
            md["species_with_id"] = md["BacID"].astype(str) + " | " + md["species"].astype(str)
    if "species_with_id" not in md.columns:
        return None
    return md


def _save_fig(fig: plt.Figure, outdir: Path, name: str) -> Path:
    p = outdir / name
    fig.savefig(p, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return p


def _bokeh_taxonomy_palette(labels: Sequence[str]) -> Dict[str, str]:
    """Map taxonomy labels to hex colours (Bokeh Category20)."""
    uniq = list(dict.fromkeys(labels))
    palette = Category20[max(3, min(20, len(uniq)))]
    return {u: palette[i % len(palette)] for i, u in enumerate(uniq)}


def _cluster_metabolite_order(mat: np.ndarray, met_names: List[str]) -> List[str]:
    """Hierarchical cluster metabolite columns (Jaccard + average linkage)."""
    if mat.shape[1] < 2:
        return met_names
    col_dist = pdist(mat.T, metric="jaccard")
    col_link = linkage(col_dist, method="average")
    order = leaves_list(col_link)
    return [met_names[i] for i in order]


def _prepare_big_picture_matrix(
    ad: ActivityData,
    *,
    taxonomy_col: str = "order",
    top_mets: int = 240,
    max_rows: int = 1200,
) -> Tuple[np.ndarray, List[str], List[str], List[str], List[str], List[int]]:
    """
    Shared row/col prep for static + interactive heatmaps.
    Returns (matrix, row_species, row_tax, col_mets, col_mets_clustered, row_breadth).
    """
    tax = (
        ad.df[taxonomy_col].fillna("Unknown").astype(str)
        if taxonomy_col in ad.df.columns
        else pd.Series(["Unknown"] * len(ad.df))
    )
    species = ad.key.astype(str).tolist()
    prevalence = ad.mets_bin.mean(axis=0).sort_values(ascending=False)
    keep_cols = prevalence.head(top_mets).index.tolist()

    d = ad.mets_bin[keep_cols].copy()
    d["__tax__"] = tax.values
    d["__breadth__"] = ad.breadth.values
    d["__species__"] = species
    d = d.sort_values(["__tax__", "__breadth__"], ascending=[True, False])

    if len(d) > max_rows:
        parts = []
        per = max(10, int(max_rows / max(1, d["__tax__"].nunique())))
        for _, grp in d.groupby("__tax__", sort=False):
            parts.append(grp.head(per))
        d = pd.concat(parts, axis=0)

    mat = d[keep_cols].to_numpy(dtype=np.float32)
    clustered_cols = _cluster_metabolite_order(mat, keep_cols)
    col_idx = [keep_cols.index(c) for c in clustered_cols]
    mat = mat[:, col_idx]

    return (
        mat,
        d["__species__"].tolist(),
        d["__tax__"].tolist(),
        keep_cols,
        clustered_cols,
        d["__breadth__"].astype(int).tolist(),
    )


def bokeh_correlation_heatmap(
    ad: ActivityData,
    outdir: Path,
    *,
    top_mets: int = 120,
) -> Optional[Path]:
    """Interactive phi-correlation heatmap with hover for metabolite pairs."""
    if not BOKEH_AVAILABLE:
        print("[warn] bokeh not installed; skipping interactive correlation heatmap")
        return None

    prevalence = ad.mets_bin.mean(axis=0).sort_values(ascending=False)
    keep = prevalence.head(top_mets).index.tolist()
    x = ad.mets_bin[keep].to_numpy(dtype=np.float32)
    phi = np.corrcoef(x, rowvar=False)
    order = _cluster_metabolite_order(x, keep)
    idx = [keep.index(m) for m in order]
    phi = phi[np.ix_(idx, idx)]
    names = order
    n = len(names)

    xs, ys, phis, ma, mb = [], [], [], [], []
    for i in range(n):
        for j in range(n):
            xs.append(i)
            ys.append(j)
            phis.append(float(phi[i, j]))
            ma.append(names[i])
            mb.append(names[j])

    source = ColumnDataSource(
        data={
            "x": xs,
            "y": ys,
            "phi": phis,
            "met_a": ma,
            "met_b": mb,
            "label": [f"{a} ↔ {b}" for a, b in zip(ma, mb)],
        }
    )

    mapper = LinearColorMapper(palette=list(reversed(RdBu11)), low=-1.0, high=1.0)
    p = figure(
        title=f"Correlation structure (phi): top {top_mets} {ad.name} metabolites (clustered)",
        width=920,
        height=880,
        tools="pan,wheel_zoom,box_zoom,reset,save,tap",
        x_range=(-0.5, n - 0.5),
        y_range=(-0.5, n - 0.5),
        toolbar_location="above",
    )
    p.grid.grid_line_color = None
    p.axis.visible = False
    p.rect(
        x="x",
        y="y",
        width=1,
        height=1,
        source=source,
        line_color=None,
        fill_color={"field": "phi", "transform": mapper},
    )
    p.add_tools(
        HoverTool(
            tooltips=[
                ("pair", "@label"),
                ("metabolite (row)", "@met_a"),
                ("metabolite (col)", "@met_b"),
                ("φ", "@phi{0.000}"),
            ],
            mode="mouse",
        )
    )
    color_bar = ColorBar(
        color_mapper=mapper,
        ticker=BasicTicker(desired_num_ticks=9),
        label_standoff=8,
        border_line_color=None,
        location=(0, 0),
        title="φ",
    )
    p.add_layout(color_bar, "right")

    # Metabolite axis labels on hover bands (top + left strips)
    div = Div(
        text=(
            f"<b>How to explore:</b> hover any cell for metabolite names and φ. "
            f"Columns/rows follow the same hierarchical cluster order as the static PNG "
            f"({n} metabolites). Scroll/zoom to inspect branches."
        ),
        width=900,
    )

    out = outdir / f"3_1_correlation_heatmap_{ad.name}_interactive.html"
    output_file(out, title=f"Correlation heatmap — {ad.name}")
    save(column(div, p))
    return out


def bokeh_big_picture_heatmap(
    ad: ActivityData,
    outdir: Path,
    *,
    taxonomy_col: str = "order",
    top_mets: int = 240,
    max_rows: int = 1200,
) -> Optional[Path]:
    """Interactive utilization landscape heatmap with taxonomy legend + rich hover."""
    if not BOKEH_AVAILABLE:
        print("[warn] bokeh not installed; skipping interactive big-picture heatmap")
        return None

    mat, row_species, row_tax, _keep, col_mets, row_breadth = _prepare_big_picture_matrix(
        ad,
        taxonomy_col=taxonomy_col,
        top_mets=top_mets,
        max_rows=max_rows,
    )
    n_rows, n_cols = mat.shape
    color_map = _bokeh_taxonomy_palette(row_tax)

    xs, ys, vals, sp, tx, met, br = [], [], [], [], [], [], []
    for i in range(n_rows):
        for j in range(n_cols):
            xs.append(j)
            ys.append(i)
            vals.append(float(mat[i, j]))
            sp.append(row_species[i])
            tx.append(row_tax[i])
            met.append(col_mets[j])
            br.append(row_breadth[i])

    source = ColumnDataSource(
        data={
            "x": xs,
            "y": ys,
            "value": vals,
            "species": sp,
            "taxonomy": tx,
            "metabolite": met,
            "breadth": br,
        }
    )

    mapper = LinearColorMapper(palette=["#ffffff", "#1a1a1a"], low=0.0, high=1.0)
    p = figure(
        title=f"Big picture: {ad.name} metabolite landscape (rows by {taxonomy_col}, cols clustered)",
        width=1100,
        height=720,
        tools="pan,wheel_zoom,box_zoom,reset,save,tap",
        x_range=(-0.5, n_cols - 0.5),
        y_range=(-0.5, n_rows - 0.5),
        toolbar_location="above",
    )
    p.grid.grid_line_color = None
    p.axis.visible = False
    p.rect(
        x="x",
        y="y",
        width=1,
        height=1,
        source=source,
        line_color=None,
        fill_color={"field": "value", "transform": mapper},
    )
    p.add_tools(
        HoverTool(
            tooltips=[
                ("species", "@species"),
                (taxonomy_col, "@taxonomy"),
                ("metabolite", "@metabolite"),
                ("utilized", "@value"),
                ("species breadth", "@breadth"),
            ],
            mode="mouse",
        )
    )
    p.add_layout(
        ColorBar(color_mapper=mapper, ticker=BasicTicker(), location=(0, 0), title="signal"),
        "right",
    )

    # Taxonomy colour strip + legend
    strip_x, strip_y, strip_c = [], [], []
    for i, t in enumerate(row_tax):
        strip_x.append(-1.2)
        strip_y.append(i)
        strip_c.append(color_map.get(t, "#888888"))
    strip = figure(
        width=40,
        height=720,
        x_range=(-1.5, -0.5),
        y_range=p.y_range,
        toolbar_location=None,
        tools="",
    )
    strip.axis.visible = False
    strip.grid.grid_line_color = None
    strip.rect(x=strip_x, y=strip_y, width=0.8, height=1, color=strip_c, line_color=None)

    legend_items = "".join(
        f'<span style="display:inline-block;width:14px;height:14px;background:{c};margin-right:6px;vertical-align:middle;"></span>'
        f'<span style="margin-right:14px;">{t}</span>'
        for t, c in sorted(color_map.items(), key=lambda kv: kv[0])
    )
    legend = Div(
        text=f"<h4>{taxonomy_col.title()} colours</h4><div style='line-height:1.8'>{legend_items}</div>",
        width=1050,
    )
    help_div = Div(
        text=(
            "<b>How to explore:</b> hover cells for species, taxonomy, metabolite, and signal. "
            "Left strip shows row taxonomy colour. Column order = hierarchical cluster of metabolites "
            "(Jaccard + average linkage)."
        ),
        width=1050,
    )

    out = outdir / f"3_1_big_picture_heatmap_{ad.name}_interactive.html"
    output_file(out, title=f"Big picture heatmap — {ad.name}")
    save(column(help_div, legend, row(strip, p)))
    return out


def bokeh_taxonomic_trends(
    activities: Dict[str, ActivityData],
    outdir: Path,
    *,
    taxonomy_col: str = "order",
    top_k_taxa: int = 20,
) -> Optional[Path]:
    """Interactive point plot: mean breadth by taxonomy × activity."""
    if not BOKEH_AVAILABLE:
        print("[warn] bokeh not installed; skipping interactive taxonomic trends")
        return None

    parts = []
    for a, ad in activities.items():
        if taxonomy_col not in ad.df.columns:
            continue
        parts.append(
            pd.DataFrame(
                {
                    "activity": a,
                    taxonomy_col: ad.df[taxonomy_col].fillna("Unknown").astype(str),
                    "breadth": ad.breadth,
                }
            )
        )
    d = pd.concat(parts, ignore_index=True)
    top_taxa = (
        d.groupby(taxonomy_col)["breadth"].mean().sort_values(ascending=False).head(top_k_taxa).index.tolist()
    )
    d = d[d[taxonomy_col].isin(top_taxa)]
    summary = d.groupby([taxonomy_col, "activity"], as_index=False)["breadth"].mean()
    summary["tax_rank"] = summary[taxonomy_col].map({t: i for i, t in enumerate(top_taxa[::-1])})

    activities_list = sorted(d["activity"].unique())
    dodge = 0.22
    offsets = {a: (i - (len(activities_list) - 1) / 2) * dodge for i, a in enumerate(activities_list)}

    p = figure(
        title=f"Taxonomic patterns: mean breadth by {taxonomy_col} (top {top_k_taxa})",
        width=980,
        height=520,
        tools="pan,wheel_zoom,box_zoom,reset,save,tap",
        y_range=(-0.5, len(top_taxa) - 0.5),
        toolbar_location="above",
    )
    p.yaxis.ticker = list(range(len(top_taxa)))
    p.yaxis.major_label_overrides = {i: t for i, t in enumerate(top_taxa[::-1])}

    for act in activities_list:
        sub = summary[summary["activity"] == act]
        src = ColumnDataSource(
            sub.assign(x=sub["breadth"], y=sub["tax_rank"] + offsets[act], activity=act)
        )
        p.scatter(
            x="x",
            y="y",
            size=10,
            alpha=0.85,
            source=src,
            legend_label=act,
        )
        p.add_tools(
            HoverTool(
                renderers=[p.renderers[-1]],
                tooltips=[
                    (taxonomy_col, f"@{{{taxonomy_col}}}"),
                    ("activity", "@activity"),
                    ("mean breadth", "@breadth{0.0}"),
                ],
            )
        )
    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    p.xaxis.axis_label = "Mean # metabolites (binarized)"

    out = outdir / f"3_1_taxonomic_trends_{taxonomy_col}_interactive.html"
    output_file(out, title=f"Taxonomic trends — {taxonomy_col}")
    save(p)
    return out


def _align_prod_util(
    activities: Dict[str, ActivityData],
    *,
    min_breadth: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], int]:
    """
    Align production/utilization matrices on shared species keys and metabolite columns.

    Drops species with production or utilization breadth below ``min_breadth``
    (default 2 = exclude species with ≤1 metabolite in either activity).
    Returns (prod_bin, util_bin, keys, n_species_before_filter).
    """
    prod = activities["production"]
    util = activities["utilization"]
    keys_prod = prod.key.astype(str)
    keys_util = util.key.astype(str)

    common_keys = sorted(set(keys_prod) & set(keys_util))
    if not common_keys:
        raise ValueError("No overlapping species keys between production and utilization tables")

    prod_idx = keys_prod[keys_prod.isin(common_keys)].index
    util_idx = keys_util[keys_util.isin(common_keys)].index
    prod_sub = prod.df.loc[prod_idx].copy()
    util_sub = util.df.loc[util_idx].copy()
    prod_sub["_key"] = keys_prod.loc[prod_idx].values
    util_sub["_key"] = keys_util.loc[util_idx].values

    prod_sub = prod_sub.drop_duplicates(subset=["_key"], keep="first").set_index("_key").loc[common_keys]
    util_sub = util_sub.drop_duplicates(subset=["_key"], keep="first").set_index("_key").loc[common_keys]

    prod_mets = set(_get_metabolite_columns(prod_sub))
    util_mets = set(_get_metabolite_columns(util_sub))
    shared_mets = sorted(prod_mets & util_mets)
    if not shared_mets:
        raise ValueError("No shared metabolite columns between production and utilization")

    prod_bin = _binarize_activity(prod_sub, "production")
    util_bin = _binarize_activity(util_sub, "utilization")
    prod_bin.index = prod_sub.index
    util_bin.index = util_sub.index
    prod_bin = prod_bin[shared_mets]
    util_bin = util_bin[shared_mets]

    n_before = len(common_keys)
    if min_breadth > 1:
        prod_b = prod_bin.sum(axis=1)
        util_b = util_bin.sum(axis=1)
        keep = (prod_b >= min_breadth) & (util_b >= min_breadth)
        prod_bin = prod_bin.loc[keep]
        util_bin = util_bin.loc[keep]
        common_keys = prod_bin.index.tolist()

    if not common_keys:
        raise ValueError(
            f"No species left after breadth filter (min_breadth={min_breadth}); "
            "lower --synergy-min-breadth or disable filtering with 1"
        )
    return prod_bin, util_bin, common_keys, n_before


def export_synergy_groups(
    activities: Dict[str, ActivityData],
    outdir: Path,
    *,
    top_n: int = 20,
    max_triplet_seeds: int = 200,
    min_breadth: int = 2,
) -> Dict[str, Path]:
    """
    Export complementary species groupings ranked by bidirectional prod↔util synergy.

    Pair score = |A prod ∩ B util| + |B prod ∩ A util|.
    Triplet score = sum of all 6 directed edges among the three species.

    Species with production or utilization breadth < ``min_breadth`` are excluded
    (default 2 removes species with ≤1 metabolite in either activity).
    """
    prod_bin, util_bin, keys, n_before = _align_prod_util(activities, min_breadth=min_breadth)
    print(
        f"[info] synergy species pool: {len(keys):,} / {n_before:,} "
        f"(min prod & util breadth ≥ {min_breadth})"
    )
    key_to_species: Dict[str, str] = {}
    for ad_name in ("production", "utilization"):
        ad = activities[ad_name]
        for k, sp in zip(
            ad.key.astype(str),
            ad.df.get("species", pd.Series([""] * len(ad.df))).astype(str),
        ):
            key_to_species.setdefault(k, sp)

    p = prod_bin.to_numpy(dtype=np.int8)
    u = util_bin.to_numpy(dtype=np.int8)
    met_cols = list(prod_bin.columns)
    # directed[i,j] = metabolites produced by i utilized by j
    directed = p @ u.T
    synergy = directed + directed.T
    np.fill_diagonal(synergy, 0)

    n = len(keys)
    key_to_idx = {k: i for i, k in enumerate(keys)}

    # Score all pairs via matrix multiply; only expand metabolite lists for top candidates.
    tri_i, tri_j = np.triu_indices(n, k=1)
    scores_flat = synergy[tri_i, tri_j]
    positive = scores_flat > 0
    tri_i = tri_i[positive]
    tri_j = tri_j[positive]
    scores_flat = scores_flat[positive]

    # Rank all positive pairs; materialize details for top_n + triplet seeds only.
    rank_cap = min(len(scores_flat), max(top_n, max_triplet_seeds) * 4)
    if rank_cap < len(scores_flat):
        top_local = np.argpartition(scores_flat, -rank_cap)[-rank_cap:]
        order = top_local[np.argsort(scores_flat[top_local])[::-1]]
    else:
        order = np.argsort(scores_flat)[::-1]

    pair_rows: List[dict] = []
    for loc in order:
        i, j = int(tri_i[loc]), int(tri_j[loc])
        score = int(scores_flat[loc])
        ka, kb = keys[i], keys[j]
        n_ab = int(directed[i, j])
        n_ba = int(directed[j, i])
        mask_ab = (p[i] == 1) & (u[j] == 1)
        mask_ba = (p[j] == 1) & (u[i] == 1)
        mets_ab = [met_cols[k] for k in np.flatnonzero(mask_ab)]
        mets_ba = [met_cols[k] for k in np.flatnonzero(mask_ba)]
        pair_rows.append(
            {
                "species_a_key": ka,
                "species_b_key": kb,
                "species_a": key_to_species.get(ka, ka),
                "species_b": key_to_species.get(kb, kb),
                "synergy_score": score,
                "a_produces_b_utilizes_n": n_ab,
                "b_produces_a_utilizes_n": n_ba,
                "a_produces_b_utilizes_mets": "; ".join(mets_ab[:80])
                + ("; …" if len(mets_ab) > 80 else ""),
                "b_produces_a_utilizes_mets": "; ".join(mets_ba[:80])
                + ("; …" if len(mets_ba) > 80 else ""),
            }
        )

    pairs_df = pd.DataFrame(pair_rows).sort_values("synergy_score", ascending=False).head(top_n)
    pairs_path = outdir / "3_3_synergy_top_pairs.csv"
    pairs_df.to_csv(pairs_path, index=False)

    # Triplets: extend top pair seeds with best third partner (vectorized directed sums)
    pair_rows_sorted = sorted(pair_rows, key=lambda r: r["synergy_score"], reverse=True)
    seed_pairs = pair_rows_sorted[:max_triplet_seeds]
    triplet_rows: List[dict] = []
    seen_triplets: set = set()

    for pr in seed_pairs:
        i = key_to_idx[pr["species_a_key"]]
        j = key_to_idx[pr["species_b_key"]]
        # total synergy of adding k to pair (i,j): all 6 directed edges
        scores = directed[i, :] + directed[:, i] + directed[j, :] + directed[:, j]
        scores[i] = -1
        scores[j] = -1
        k_best = int(np.argmax(scores))
        if scores[k_best] <= 0:
            continue
        trip = tuple(sorted([keys[i], keys[j], keys[k_best]]))
        if trip in seen_triplets:
            continue
        seen_triplets.add(trip)
        ki, kj, kk = key_to_idx[trip[0]], key_to_idx[trip[1]], key_to_idx[trip[2]]
        n_ab = int(directed[ki, kj])
        n_ba = int(directed[kj, ki])
        n_ac = int(directed[ki, kk])
        n_ca = int(directed[kk, ki])
        n_bc = int(directed[kj, kk])
        n_cb = int(directed[kk, kj])
        score = n_ab + n_ba + n_ac + n_ca + n_bc + n_cb

        def _mets(a_idx: int, b_idx: int) -> str:
            m = [met_cols[t] for t in np.flatnonzero((p[a_idx] == 1) & (u[b_idx] == 1))]
            return "; ".join(m[:25])

        triplet_rows.append(
            {
                "species_a_key": trip[0],
                "species_b_key": trip[1],
                "species_c_key": trip[2],
                "species_a": key_to_species.get(trip[0], trip[0]),
                "species_b": key_to_species.get(trip[1], trip[1]),
                "species_c": key_to_species.get(trip[2], trip[2]),
                "synergy_score": score,
                "directed_edge_counts": f"ab={n_ab},ba={n_ba},ac={n_ac},ca={n_ca},bc={n_bc},cb={n_cb}",
                "example_mets_a_to_b": _mets(ki, kj),
                "example_mets_a_to_c": _mets(ki, kk),
                "example_mets_b_to_c": _mets(kj, kk),
            }
        )

    triplets_df = (
        pd.DataFrame(triplet_rows)
        .sort_values("synergy_score", ascending=False)
        .drop_duplicates(subset=["species_a_key", "species_b_key", "species_c_key"])
        .head(top_n)
    )
    triplets_path = outdir / "3_3_synergy_top_triplets.csv"
    triplets_df.to_csv(triplets_path, index=False)

    quartet_rows: List[dict] = []
    for _, tr in triplets_df.head(min(50, len(triplets_df))).iterrows():
        trip_keys = [tr["species_a_key"], tr["species_b_key"], tr["species_c_key"]]
        idxs = [key_to_idx[k] for k in trip_keys]
        scores = np.zeros(n, dtype=np.int64)
        for a in idxs:
            scores += directed[a, :] + directed[:, a]
        for a in idxs:
            scores[a] = -1
        kd = int(np.argmax(scores))
        if scores[kd] <= 0:
            continue
        q_keys = tuple(sorted(trip_keys + [keys[kd]]))
        quartet_rows.append(
            {
                "species_keys": " | ".join(q_keys),
                "species_names": " | ".join(key_to_species.get(k, k) for k in q_keys),
                "synergy_score": int(scores[kd]),
            }
        )
    quartets_df = (
        pd.DataFrame(quartet_rows)
        .sort_values("synergy_score", ascending=False)
        .drop_duplicates(subset=["species_keys"])
        .head(top_n)
    )
    quartets_path = outdir / "3_3_synergy_top_quartets.csv"
    quartets_df.to_csv(quartets_path, index=False)

    return {"pairs": pairs_path, "triplets": triplets_path, "quartets": quartets_path}


def _truncate_label(text: str, max_len: int = 42) -> str:
    s = str(text).strip()
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + "…"


def _parse_directed_edge_counts(
    counts_str: str,
    species: Tuple[str, str, str],
) -> List[Tuple[str, str, int]]:
    """Parse triplet edge string 'ab=4,ba=0,...' into (source, target, weight) tuples."""
    mapping = {
        "ab": (species[0], species[1]),
        "ba": (species[1], species[0]),
        "ac": (species[0], species[2]),
        "ca": (species[2], species[0]),
        "bc": (species[1], species[2]),
        "cb": (species[2], species[1]),
    }
    edges: List[Tuple[str, str, int]] = []
    for key, (src, dst) in mapping.items():
        m = re.search(rf"\b{key}=(-?\d+)", str(counts_str))
        if m:
            w = int(m.group(1))
            if w > 0:
                edges.append((src, dst, w))
    return edges


def _synergy_edges_from_pairs(pairs_df: pd.DataFrame) -> List[dict]:
    rows: List[dict] = []
    for _, r in pairs_df.iterrows():
        a, b = str(r["species_a"]), str(r["species_b"])
        n_ab = int(r.get("a_produces_b_utilizes_n", 0) or 0)
        n_ba = int(r.get("b_produces_a_utilizes_n", 0) or 0)
        m_ab = str(r.get("a_produces_b_utilizes_mets", "") or "")
        m_ba = str(r.get("b_produces_a_utilizes_mets", "") or "")
        if n_ab > 0:
            rows.append(
                {
                    "source": a,
                    "target": b,
                    "weight": n_ab,
                    "metabolites": m_ab,
                    "group": "pair",
                }
            )
        if n_ba > 0:
            rows.append(
                {
                    "source": b,
                    "target": a,
                    "weight": n_ba,
                    "metabolites": m_ba,
                    "group": "pair",
                }
            )
    return rows


def _synergy_edges_from_triplets(triplets_df: pd.DataFrame) -> List[dict]:
    rows: List[dict] = []
    for _, r in triplets_df.iterrows():
        sp = (str(r["species_a"]), str(r["species_b"]), str(r["species_c"]))
        for src, dst, w in _parse_directed_edge_counts(r.get("directed_edge_counts", ""), sp):
            rows.append({"source": src, "target": dst, "weight": w, "metabolites": "", "group": "triplet"})
    return rows


def plot_synergy_pairs_bidirectional(pairs_df: pd.DataFrame, outdir: Path, *, top_n: int = 15) -> Path:
    """Stacked horizontal bars: metabolites A→B vs B→A for top synergy pairs."""
    d = pairs_df.sort_values("synergy_score", ascending=True).tail(top_n).copy()
    d["pair_label"] = [
        _truncate_label(f"{a} ↔ {b}", 50) for a, b in zip(d["species_a"], d["species_b"])
    ]

    fig, ax = plt.subplots(figsize=(11, max(5, 0.38 * len(d) + 1.5)))
    y = np.arange(len(d))
    a2b = d["a_produces_b_utilizes_n"].astype(int)
    b2a = d["b_produces_a_utilizes_n"].astype(int)
    ax.barh(y, a2b, height=0.72, label="A produces → B utilizes", color="#4C78A8", alpha=0.92)
    ax.barh(y, b2a, height=0.72, left=a2b, label="B produces → A utilizes", color="#F58518", alpha=0.92)
    for yi, score in enumerate(d["synergy_score"]):
        ax.text(
            float(a2b.values[yi] + b2a.values[yi]) + 0.15,
            yi,
            f"Σ {int(score)}",
            va="center",
            fontsize=8,
            color="#333",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(d["pair_label"], fontsize=8)
    ax.set_xlabel("# shared metabolites (directed)")
    ax.set_title(f"Top {len(d)} complementary species pairs (prod ↔ util synergy)")
    ax.legend(loc="lower right", frameon=True)
    sns.despine(ax=ax, left=False)
    fig.tight_layout()
    return _save_fig(fig, outdir, "3_3_synergy_pairs_bidirectional.png")


def plot_synergy_rank_summary(
    pairs_df: pd.DataFrame,
    triplets_df: pd.DataFrame,
    quartets_df: pd.DataFrame,
    outdir: Path,
    *,
    top_n: int = 12,
) -> Path:
    """Three-panel lollipop chart of synergy scores for pairs, triplets, quartets."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2), gridspec_kw={"width_ratios": [1.1, 1.15, 1.0]})

    def _panel(ax, df, label_col: str, title: str, color: str) -> None:
        d = df.sort_values("synergy_score", ascending=True).tail(top_n)
        labels = [_truncate_label(x, 36) for x in d[label_col]]
        y = np.arange(len(d))
        ax.hlines(y, 0, d["synergy_score"], color=color, alpha=0.35, linewidth=2)
        ax.scatter(d["synergy_score"], y, color=color, s=48, zorder=3, edgecolors="white", linewidths=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Synergy score")
        ax.set_title(title, fontsize=11)
        sns.despine(ax=ax, left=False)

    pairs_for_panel = pairs_df.assign(
        pair_label=lambda x: x["species_a"].astype(str) + " ↔ " + x["species_b"].astype(str)
    )
    _panel(axes[0], pairs_for_panel, "pair_label", f"Pairs (top {min(top_n, len(pairs_df))})", "#4C78A8")

    triplets_for_panel = triplets_df.assign(
        triplet_label=lambda x: x["species_a"].astype(str)
        + " + "
        + x["species_b"].astype(str)
        + " + "
        + x["species_c"].astype(str)
    )
    _panel(
        axes[1],
        triplets_for_panel,
        "triplet_label",
        f"Triplets (top {min(top_n, len(triplets_df))})",
        "#54A24B",
    )

    quartets_for_panel = quartets_df.assign(
        quartet_label=lambda x: x["species_names"].astype(str).str.replace(" | ", "\n", regex=False)
    )
    _panel(
        axes[2],
        quartets_for_panel,
        "quartet_label",
        f"Quartets (top {min(top_n, len(quartets_df))})",
        "#B279A2",
    )

    fig.suptitle("Complementary metabolite synergy (production ↔ utilization)", y=1.02, fontsize=13)
    fig.tight_layout()
    return _save_fig(fig, outdir, "3_3_synergy_rank_summary.png")


def bokeh_synergy_network(
    pairs_df: pd.DataFrame,
    triplets_df: pd.DataFrame,
    outdir: Path,
    *,
    max_edges: int = 80,
) -> Optional[Path]:
    """Interactive directed network of top pair + triplet prod→util edges."""
    if not BOKEH_AVAILABLE:
        print("[warn] bokeh not installed; skipping synergy network")
        return None

    edge_rows = _synergy_edges_from_pairs(pairs_df) + _synergy_edges_from_triplets(triplets_df)
    if not edge_rows:
        return None

    edges_df = pd.DataFrame(edge_rows)
    edges_df = edges_df.sort_values("weight", ascending=False).head(max_edges)

    nodes = sorted(set(edges_df["source"]) | set(edges_df["target"]))
    n = len(nodes)
    node_idx = {name: i for i, name in enumerate(nodes)}
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    radius = 1.0 + 0.08 * min(n, 30)
    xs = (np.cos(angles) * radius).tolist()
    ys = (np.sin(angles) * radius).tolist()
    degree = {nd: 0 for nd in nodes}
    for _, e in edges_df.iterrows():
        degree[e["source"]] += int(e["weight"])
        degree[e["target"]] += int(e["weight"])

    curve_xs, curve_ys, ew, elabel, emets = [], [], [], [], []
    for _, e in edges_df.iterrows():
        i, j = node_idx[e["source"]], node_idx[e["target"]]
        x0, y0, x1, y1 = xs[i], ys[i], xs[j], ys[j]
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        bend = 0.22 * radius
        cx = mx + bend * (y1 - y0)
        cy = my - bend * (x1 - x0)
        t_vals = np.linspace(0, 1, 24)
        curve_x = (1 - t_vals) ** 2 * x0 + 2 * (1 - t_vals) * t_vals * cx + t_vals**2 * x1
        curve_y = (1 - t_vals) ** 2 * y0 + 2 * (1 - t_vals) * t_vals * cy + t_vals**2 * y1
        curve_xs.append(curve_x.tolist())
        curve_ys.append(curve_y.tolist())
        w = int(e["weight"])
        ew.append(w)
        emets.append(str(e.get("metabolites") or ""))
        elabel.append(f"{e['source']} → {e['target']} ({w} met.)")

    edge_source = ColumnDataSource(
        data={
            "xs": curve_xs,
            "ys": curve_ys,
            "weight": ew,
            "label": elabel,
            "metabolites": emets,
        }
    )

    node_degrees = [degree[nd] for nd in nodes]
    node_sizes = [12 + 5 * min(d, 24) for d in node_degrees]
    node_source = ColumnDataSource(
        data={
            "x": xs,
            "y": ys,
            "name": nodes,
            "degree": node_degrees,
            "size": node_sizes,
            "label": [_truncate_label(nd, 38) for nd in nodes],
        }
    )

    p = figure(
        title="Complementary metabolite synergy network (directed prod → util)",
        width=980,
        height=820,
        x_range=(-1.55, 1.55),
        y_range=(-1.55, 1.55),
        tools="pan,wheel_zoom,box_zoom,reset,save,tap",
        toolbar_location="above",
    )
    p.axis.visible = False
    p.grid.grid_line_color = None

    edge_renderer = p.multi_line(
        xs="xs",
        ys="ys",
        source=edge_source,
        line_color="#7A7A7A",
        line_alpha=0.6,
        line_width="weight",
    )
    node_renderer = p.scatter(
        x="x",
        y="y",
        size="size",
        source=node_source,
        fill_color="#2C7FB8",
        line_color="white",
        line_width=1.5,
        fill_alpha=0.92,
    )
    p.add_tools(
        HoverTool(
            renderers=[edge_renderer],
            tooltips=[
                ("edge", "@label"),
                ("metabolites", "@metabolites"),
            ],
        ),
        HoverTool(
            renderers=[node_renderer],
            tooltips=[("species", "@name"), ("weighted degree", "@degree")],
        ),
    )

    div = Div(
        text=(
            "<b>How to explore:</b> nodes = species; edge width ∝ # metabolites produced by source "
            "and utilized by target. Hover edges for direction and metabolites (pairs). "
            "Includes top pair + triplet synergies."
        ),
        width=940,
    )

    out = outdir / "3_3_synergy_network_interactive.html"
    output_file(out, title="Synergy network")
    save(column(div, p))
    return out


def visualize_synergy_outputs(
    synergy_paths: Dict[str, Path],
    outdir: Path,
    *,
    skip_bokeh: bool = False,
) -> List[str]:
    """Build static + interactive figures from synergy CSV exports."""
    pairs_df = pd.read_csv(synergy_paths["pairs"])
    triplets_df = pd.read_csv(synergy_paths["triplets"])
    quartets_df = pd.read_csv(synergy_paths["quartets"])

    written: List[str] = []
    written.append(plot_synergy_pairs_bidirectional(pairs_df, outdir).name)
    written.append(
        plot_synergy_rank_summary(pairs_df, triplets_df, quartets_df, outdir).name
    )
    if not skip_bokeh:
        net = bokeh_synergy_network(pairs_df, triplets_df, outdir)
        if net is not None:
            written.append(net.name)
    return written


def plot_overview_breadth(activities: Dict[str, ActivityData], outdir: Path) -> Path:
    rows = []
    for a, ad in activities.items():
        rows.append(pd.DataFrame({"activity": a, "breadth": ad.breadth}))
    df = pd.concat(rows, ignore_index=True)
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    sns.violinplot(data=df, x="activity", y="breadth", inner="quartile", ax=ax, cut=0)
    ax.set_title("Overview: distribution of metabolite breadth across activities (species-level)")
    ax.set_xlabel("")
    ax.set_ylabel("# metabolites with positive signal (binarized)")
    return _save_fig(fig, outdir, "3_1_overview_breadth_violin.png")


def plot_taxonomic_trends(
    activities: Dict[str, ActivityData],
    outdir: Path,
    *,
    taxonomy_col: str = "order",
    top_k_taxa: int = 20,
) -> Path:
    """
    Show mean breadth by taxonomy (order) for each activity.
    Note: phylum is not present in the provided species_data tables; `order` is the highest rank available.
    """
    parts = []
    for a, ad in activities.items():
        if taxonomy_col not in ad.df.columns:
            continue
        tmp = pd.DataFrame(
            {
                "activity": a,
                taxonomy_col: ad.df[taxonomy_col].fillna("Unknown").astype(str),
                "breadth": ad.breadth,
            }
        )
        parts.append(tmp)
    d = pd.concat(parts, ignore_index=True)
    top_taxa = (
        d.groupby(taxonomy_col)["breadth"].mean().sort_values(ascending=False).head(top_k_taxa).index.tolist()
    )
    d = d[d[taxonomy_col].isin(top_taxa)]

    fig, ax = plt.subplots(figsize=(12.5, 6.0))
    sns.pointplot(
        data=d,
        y=taxonomy_col,
        x="breadth",
        hue="activity",
        dodge=0.5,
        ax=ax,
        errorbar=("ci", 95),
    )
    ax.set_title(f"Taxonomic patterns: mean breadth by {taxonomy_col} (top {top_k_taxa})")
    ax.set_xlabel("Mean # metabolites (binarized)")
    ax.set_ylabel("")
    ax.legend(title="activity", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return _save_fig(fig, outdir, f"3_1_taxonomic_trends_{taxonomy_col}_pointplot.png")


def _top_correlated_pairs(phi: np.ndarray, names: List[str], n: int = 25) -> pd.DataFrame:
    tri = np.triu_indices_from(phi, k=1)
    vals = phi[tri]
    order = np.argsort(vals)[::-1]
    rows = []
    for idx in order[:n]:
        i = tri[0][idx]
        j = tri[1][idx]
        rows.append({"metabolite_a": names[i], "metabolite_b": names[j], "phi": float(vals[idx])})
    return pd.DataFrame(rows)


def correlation_structure(
    activities: Dict[str, ActivityData],
    outdir: Path,
    *,
    top_mets: int = 220,
) -> Tuple[Path, Path]:
    """
    Compute metabolite co-occurrence (phi correlation ≈ corr of 0/1) within each activity,
    on the most prevalent metabolites.
    Saves:
    - heatmap of correlations (utilization)
    - CSV of top correlated pairs per activity
    """
    all_top = []
    heat_path: Optional[Path] = None
    for a, ad in activities.items():
        prevalence = ad.mets_bin.mean(axis=0).sort_values(ascending=False)
        keep = prevalence.head(top_mets).index.tolist()
        x = ad.mets_bin[keep].to_numpy(dtype=np.float32)
        phi = np.corrcoef(x, rowvar=False)
        pairs = _top_correlated_pairs(phi, keep, n=30)
        pairs.insert(0, "activity", a)
        all_top.append(pairs)

        if a == "utilization":
            fig, ax = plt.subplots(figsize=(10.5, 9.5))
            sns.heatmap(phi, cmap="vlag", center=0, vmin=-1, vmax=1, ax=ax, cbar_kws={"shrink": 0.6})
            ax.set_title(f"Correlation structure (phi): top {top_mets} utilization metabolites")
            ax.set_xticks([])
            ax.set_yticks([])
            heat_path = _save_fig(fig, outdir, "3_1_correlation_heatmap_util_top_mets.png")

    top_pairs = pd.concat(all_top, ignore_index=True)
    csv_path = outdir / "3_1_top_correlated_pairs_phi.csv"
    top_pairs.to_csv(csv_path, index=False)
    return heat_path or (outdir / "3_1_correlation_heatmap_util_top_mets.png"), csv_path


def big_picture_heatmap_clustered_by_taxonomy(
    ad: ActivityData,
    outdir: Path,
    *,
    taxonomy_col: str = "order",
    top_mets: int = 240,
    max_rows: int = 1200,
) -> Path:
    """
    “Big picture” figure: cluster metabolites while grouping rows by taxonomy.
    """
    mat, _species, tax_labels, _keep, clustered_cols, _breadth = _prepare_big_picture_matrix(
        ad,
        taxonomy_col=taxonomy_col,
        top_mets=top_mets,
        max_rows=max_rows,
    )

    col_dist = pdist(mat.T, metric="jaccard")
    col_link = linkage(col_dist, method="average")

    uniq = pd.Series(tax_labels).unique().tolist()
    palette = sns.color_palette("tab20", n_colors=min(20, len(uniq)))
    color_map = {u: palette[i % len(palette)] for i, u in enumerate(uniq)}
    row_colors = pd.Series(tax_labels).map(color_map).to_numpy()

    g = sns.clustermap(
        mat,
        col_linkage=col_link,
        row_cluster=False,
        row_colors=row_colors,
        cmap="Greys",
        xticklabels=False,
        yticklabels=False,
        figsize=(12.5, 10.5),
        cbar_pos=(0.02, 0.85, 0.02, 0.12),
    )
    g.fig.suptitle(
        f"Big picture: {ad.name} metabolite landscape (top {top_mets} mets), rows grouped by {taxonomy_col}",
        y=1.02,
        fontsize=12,
    )
    out = outdir / f"3_1_big_picture_heatmap_{ad.name}_clustered.png"
    g.fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(g.fig)
    return out


def _zscore(s: pd.Series) -> pd.Series:
    s = s.astype(float)
    mu = float(s.mean())
    sd = float(s.std(ddof=0)) or 1.0
    return (s - mu) / sd


def find_outliers(ad: ActivityData, *, top_n: int = 3, mode: str = "high") -> pd.DataFrame:
    z = _zscore(ad.breadth)
    df = pd.DataFrame(
        {
            "species_key": ad.key,
            "species": ad.df.get("species", pd.Series(["Unknown"] * len(ad.df))).astype(str),
            "genus": ad.df.get("genus", pd.Series(["Unknown"] * len(ad.df))).astype(str),
            "order": ad.df.get("order", pd.Series(["Unknown"] * len(ad.df))).astype(str),
            "type_strain": ad.df.get("type_strain", pd.Series(["Unknown"] * len(ad.df))).astype(str),
            "breadth": ad.breadth.astype(int),
            "z": z.astype(float),
        }
    )
    valid = df["species_key"].notna() & (df["species"].astype(str).str.lower() != "nan")
    df = df[valid]
    df = df.sort_values("z", ascending=(mode != "high")).head(top_n)
    return df


def write_outlier_report(
    activities: Dict[str, ActivityData],
    outdir: Path,
    metadata: Optional[pd.DataFrame],
    *,
    top_n: int = 3,
) -> Path:
    lines: List[str] = []
    lines.append("# Outlier species summary\n")

    def df_to_md_table(df: pd.DataFrame, max_cols: int = 18) -> str:
        """
        Minimal markdown table renderer to avoid the optional 'tabulate' dependency.
        Truncates very wide tables for readability.
        """
        if df is None or df.empty:
            return "_(no rows)_\n"
        d = df.copy()
        if len(d.columns) > max_cols:
            d = d.iloc[:, :max_cols]
        cols = [str(c) for c in d.columns]

        def esc(x: object) -> str:
            s = "" if pd.isna(x) else str(x)
            return s.replace("|", "\\|").replace("\n", " ")

        header = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        rows = ["| " + " | ".join(esc(v) for v in row) + " |" for row in d.to_numpy()]
        return "\n".join([header, sep] + rows) + "\n"

    def enrich(df: pd.DataFrame) -> pd.DataFrame:
        if metadata is None or df.empty:
            return df
        md = metadata.drop_duplicates(subset=["species_with_id"]).copy()
        return df.merge(md, how="left", left_on="species_key", right_on="species_with_id")

    util_hi = enrich(find_outliers(activities["utilization"], top_n=top_n, mode="high"))
    prod_hi = enrich(find_outliers(activities["production"], top_n=top_n, mode="high"))
    lines.append("## 3.2.1 Metabolic generalists\n")
    lines.append(f"### Utilization breadth outliers (top {top_n})\n")
    lines.append(df_to_md_table(util_hi) + "\n")
    lines.append(f"### Production breadth outliers (top {top_n})\n")
    lines.append(df_to_md_table(prod_hi) + "\n")

    util_lo = enrich(find_outliers(activities["utilization"], top_n=top_n, mode="low"))
    prod_lo = enrich(find_outliers(activities["production"], top_n=top_n, mode="low"))
    lines.append("## 3.2.2 Metabolic specialists\n")
    lines.append(f"### Utilization breadth (bottom {top_n})\n")
    lines.append(df_to_md_table(util_lo) + "\n")
    lines.append(f"### Production breadth (bottom {top_n})\n")
    lines.append(df_to_md_table(prod_lo) + "\n")

    res_hi = enrich(find_outliers(activities["resistance"], top_n=top_n, mode="high"))
    lines.append("## 3.2.3 Resistance outliers\n")
    lines.append(f"### Resistance breadth outliers (top {top_n})\n")
    lines.append(df_to_md_table(res_hi) + "\n")

    if metadata is None:
        lines.append("## Ecological context (optional metadata)\n")
        lines.append(
            textwrap.dedent(
                """
                The requested ecological context (isolation source, habitat complexity, genome size correlations)
                is **not present** in the default `species_data/*exploded*.csv(.zip)` tables in this repo.

                If you add `experiments_results/metadata_species.csv` with a `species_with_id` column (or BacID+species),
                this script will automatically join it into the outlier tables so you can write the paragraph-level
                exemplar descriptions.
                """
            ).strip()
            + "\n"
        )
    else:
        lines.append("## Ecological context (from `metadata_species.csv`)\n")
        for col in ["isolation_source", "habitat_complexity", "genome_size_mb"]:
            if col in metadata.columns:
                lines.append(f"- Found metadata column: `{col}`\n")
        lines.append(
            "\nAdd your own narrative paragraphs using the enriched tables above (now including metadata columns).\n"
        )

    out = outdir / "3_2_outlier_species_report.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default=str(REPO_ROOT / "experiments_results"), help="Base output directory")
    parser.add_argument("--top-mets-corr", type=int, default=220)
    parser.add_argument("--top-mets-corr-interactive", type=int, default=120)
    parser.add_argument("--top-mets-big", type=int, default=240)
    parser.add_argument("--max-rows-big", type=int, default=1200)
    parser.add_argument("--synergy-top-n", type=int, default=20)
    parser.add_argument(
        "--synergy-min-breadth",
        type=int,
        default=2,
        help="Exclude species with fewer than this many production OR utilization metabolites (default 2 = drop ≤1)",
    )
    parser.add_argument("--skip-bokeh", action="store_true")
    parser.add_argument("--outlier-top-n", type=int, default=20, help="Outliers per category in 3_2 report")
    parser.add_argument(
        "--outlier-report-only",
        action="store_true",
        help="Only regenerate 3_2_outlier_species_report.md (no new outputs_* folder)",
    )
    parser.add_argument(
        "--output-run-dir",
        type=str,
        default=None,
        help="Existing outputs_YYYYMMDD_HHMMSS folder for --outlier-report-only (default: latest)",
    )
    parser.add_argument("--force-refresh-phylum-cache", action="store_true")
    parser.add_argument("--max-new-phylum-queries", type=int, default=None)
    args = parser.parse_args()

    base = Path(args.outdir)
    base.mkdir(parents=True, exist_ok=True)

    print(f"[info] loading activity data from {SPECIES_DATA} …")
    activities = load_activity_data()
    md = try_load_metadata(REPO_ROOT / "experiments_results" / "metadata_species.csv")

    if args.outlier_report_only:
        if args.output_run_dir:
            outdir = Path(args.output_run_dir)
        else:
            runs = sorted(base.glob("outputs_*"), key=lambda p: p.stat().st_mtime, reverse=True)
            if not runs:
                raise SystemExit(f"No outputs_* folders under {base}")
            outdir = runs[0]
        if not outdir.is_dir():
            raise SystemExit(f"Output run directory not found: {outdir}")
        print(f"[info] writing outlier report (top {args.outlier_top_n} per category) → {outdir} …")
        report_path = write_outlier_report(activities, outdir, md, top_n=args.outlier_top_n)
        print(f"[done] wrote {report_path}")
        return 0

    outdir = _ensure_outdir(base)
    phylum_cache_path = REPO_ROOT / "experiments_results" / "genus_phylum_cache.csv"
    print(f"[info] building/reusing phylum cache at {phylum_cache_path} …")
    phylum_cache = build_or_load_phylum_cache(
        activities,
        phylum_cache_path,
        force_refresh=args.force_refresh_phylum_cache,
        max_new_queries=args.max_new_phylum_queries,
    )
    attach_phylum_to_activities(activities, phylum_cache)

    resolved_phylum = int((phylum_cache["phylum"] != "Unknown").sum()) if "phylum" in phylum_cache.columns else 0
    taxonomy_col = "phylum" if resolved_phylum > 0 else "order"

    print("[info] plotting overview breadth …")
    plot_overview_breadth(activities, outdir)

    print(f"[info] plotting taxonomic trends ({taxonomy_col}-level) …")
    plot_taxonomic_trends(activities, outdir, taxonomy_col=taxonomy_col, top_k_taxa=20)

    print("[info] computing correlation structure …")
    correlation_structure(activities, outdir, top_mets=args.top_mets_corr)

    print(f"[info] generating big-picture heatmap (utilization, grouped by {taxonomy_col}) …")
    big_picture_heatmap_clustered_by_taxonomy(
        activities["utilization"],
        outdir,
        taxonomy_col=taxonomy_col,
        top_mets=args.top_mets_big,
        max_rows=args.max_rows_big,
    )

    bokeh_paths: List[str] = []
    if not args.skip_bokeh:
        print("[info] writing interactive Bokeh plots …")
        for fn, label in [
            (
                lambda: bokeh_correlation_heatmap(
                    activities["utilization"],
                    outdir,
                    top_mets=args.top_mets_corr_interactive,
                ),
                "correlation heatmap",
            ),
            (
                lambda: bokeh_big_picture_heatmap(
                    activities["utilization"],
                    outdir,
                    taxonomy_col=taxonomy_col,
                    top_mets=args.top_mets_big,
                    max_rows=args.max_rows_big,
                ),
                "big-picture heatmap",
            ),
            (
                lambda: bokeh_taxonomic_trends(
                    activities,
                    outdir,
                    taxonomy_col=taxonomy_col,
                    top_k_taxa=20,
                ),
                "taxonomic trends",
            ),
        ]:
            try:
                p = fn()
                if p is not None:
                    bokeh_paths.append(p.name)
            except Exception as exc:
                print(f"[warn] Bokeh {label} failed: {exc}")

    print("[info] computing complementary prod↔util synergy groups …")
    synergy_paths = export_synergy_groups(
        activities,
        outdir,
        top_n=args.synergy_top_n,
        min_breadth=args.synergy_min_breadth,
    )

    print("[info] plotting synergy visualizations …")
    synergy_viz = visualize_synergy_outputs(synergy_paths, outdir, skip_bokeh=args.skip_bokeh)

    print("[info] writing outlier report …")
    write_outlier_report(activities, outdir, md, top_n=args.outlier_top_n)

    (outdir / "README.txt").write_text(
        textwrap.dedent(
            f"""
            Outputs written to:
              {outdir}

            Taxonomy grouping used:
              {taxonomy_col}
            Cached genus->phylum map:
              {phylum_cache_path}
            Phylum mapped genera:
              {resolved_phylum} / {len(phylum_cache) if isinstance(phylum_cache, pd.DataFrame) else 0}

            Key files:
              - 3_1_overview_breadth_violin.png
              - 3_1_taxonomic_trends_{taxonomy_col}_pointplot.png
              - 3_1_correlation_heatmap_util_top_mets.png
              - 3_1_top_correlated_pairs_phi.csv
              - 3_1_big_picture_heatmap_utilization_clustered.png
              - 3_2_outlier_species_report.md

            Interactive (Bokeh HTML):
              {chr(10).join("  - " + p for p in bokeh_paths) if bokeh_paths else "  (none — install bokeh or omit --skip-bokeh)"}

            Complementary synergy groups (prod ↔ util):
              - {synergy_paths["pairs"].name}
              - {synergy_paths["triplets"].name}
              - {synergy_paths["quartets"].name}

            Synergy visualizations:
              {chr(10).join("  - " + v for v in synergy_viz) if synergy_viz else "  (none)"}
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    print(f"[done] wrote outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

