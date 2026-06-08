"""
Figures for genomic follow-up pipeline (BGC, AMR, BacDive phenotypes, gene hits).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from genome_investigation.visualize_results import CATEGORY_COLORS, species_two_line_label, _bacdive_breadth_by_species

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
        "figure.facecolor": "white",
    }
)


def _load_optional(path: Optional[Path]) -> pd.DataFrame:
    if path and Path(path).exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def fig_pipeline_status(integrated: pd.DataFrame, out: Path) -> None:
    """Horizontal status bars: download / antiSMASH / AMRFinder / gene hits."""
    d = integrated.copy()
    sp_col = "species"
    labels = [species_two_line_label(s) for s in d[sp_col]]

    def _status(col: str, ok_vals=("success", "parsed", "completed")) -> np.ndarray:
        if col not in d.columns:
            return np.zeros(len(d))
        return np.array([1.0 if str(v).lower() in ok_vals else 0.0 for v in d[col]])

    has_dl = np.array([1.0 if str(d.iloc[i].get("genome_accession", "")) else 0.0 for i in range(len(d))])
    asm = _status("antismash_status")
    amr = _status("amrfinder_status")
    hits = np.array([min(1.0, float(d.iloc[i].get("genomic_hit_n", 0) or 0) / 5.0) for i in range(len(d))])

    mat = np.column_stack([has_dl, asm, amr, hits])
    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(d))))
    im = ax.imshow(mat, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["Genome\naccession", "antiSMASH", "AMRFinder", "Gene hits"], fontsize=9)
    ax.set_yticks(np.arange(len(d)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_title("Genomic follow-up pipeline status by species")
    plt.colorbar(im, ax=ax, shrink=0.5, label="Step complete / hits (scaled)")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_bgc_and_amr(integrated: pd.DataFrame, out: Path) -> None:
    d = integrated.copy()
    if "species" not in d.columns:
        return
    if "bgc_count_total" not in d.columns:
        d["bgc_count_total"] = 0
    if "amr_gene_count" not in d.columns:
        d["amr_gene_count"] = 0

    d["bgc_count_total"] = pd.to_numeric(d["bgc_count_total"], errors="coerce").fillna(0)
    d["amr_gene_count"] = pd.to_numeric(d["amr_gene_count"], errors="coerce").fillna(0)
    d = d.sort_values("bgc_count_total", ascending=True)

    labels = [species_two_line_label(s) for s in d["species"]]
    y = np.arange(len(d))
    h = 0.35

    fig, ax = plt.subplots(figsize=(9, max(4, 0.5 * len(d))))
    ax.barh(y - h / 2, d["bgc_count_total"], height=h, color="#6A4C93", label="BGC regions (antiSMASH)")
    ax.barh(y + h / 2, d["amr_gene_count"], height=h, color="#E45756", label="AMR genes (AMRFinder)")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Count")
    ax.set_title("Biosynthetic clusters vs. acquired resistance genes")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_phenotype_breadth_vs_genomic(focal: pd.DataFrame, integrated: pd.DataFrame, out: Path) -> None:
    breadth = _bacdive_breadth_by_species(focal["species"].astype(str).tolist())
    if breadth.empty:
        return

    util = breadth[breadth["activity"] == "utilization"][["species", "breadth"]].rename(columns={"breadth": "util"})
    prod = breadth[breadth["activity"] == "production"][["species", "breadth"]].rename(columns={"breadth": "prod"})
    res = breadth[breadth["activity"] == "resistance"][["species", "breadth"]].rename(columns={"breadth": "res"})
    m = focal[["species", "match_confidence"]].merge(util, on="species", how="left").merge(prod, on="species", how="left").merge(res, on="species", how="left")
    extra_cols = ["species"] + [c for c in ("bgc_count_total", "amr_gene_count", "genomic_hit_n") if c in integrated.columns]
    if len(extra_cols) > 1:
        m = m.merge(integrated[extra_cols].drop_duplicates("species"), on="species", how="left")
    else:
        m["bgc_count_total"] = 0
        m["amr_gene_count"] = 0
        m["genomic_hit_n"] = 0
    for col in ("bgc_count_total", "amr_gene_count", "genomic_hit_n"):
        if col not in m.columns:
            m[col] = 0
        m[col] = pd.to_numeric(m[col], errors="coerce").fillna(0)

    m["total_breadth"] = m[["util", "prod", "res"]].fillna(0).sum(axis=1)
    m["genomic_score"] = m["bgc_count_total"] + m["amr_gene_count"] * 2 + m["genomic_hit_n"] * 0.5

    fig, ax = plt.subplots(figsize=(7, 5))
    sizes = 40 + 6 * pd.to_numeric(m["res"], errors="coerce").fillna(0)
    sc = ax.scatter(
        m["total_breadth"],
        m["genomic_score"],
        s=sizes,
        c=pd.to_numeric(m["match_confidence"], errors="coerce").fillna(0),
        cmap="viridis",
        edgecolors="white",
        linewidths=0.6,
    )
    for _, row in m.iterrows():
        ax.annotate(
            species_two_line_label(row["species"]),
            (row["total_breadth"], row["genomic_score"]),
            fontsize=6,
            xytext=(4, 4),
            textcoords="offset points",
        )
    plt.colorbar(sc, ax=ax, label="Match confidence")
    ax.set_xlabel("BacDive Step3 total breadth (util + prod + res)")
    ax.set_ylabel("Genomic evidence score (BGC + 2×AMR + 0.5×gene hits)")
    ax.set_title("Phenotype breadth vs. genomic follow-up evidence")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_metabolite_heatmap(pheno: pd.DataFrame, out: Path, top_n: int = 25) -> None:
    if pheno.empty:
        return
    pos = pheno[
        ((pheno["activity"] == "utilization") & (pheno["result"].isin(["+", "positive", "yes"])))
        | ((pheno["activity"] == "production") & (pheno["result"].isin(["yes", "+", "positive"])))
        | ((pheno["activity"] == "resistance") & (pheno["result"] == "resistant"))
    ].copy()
    if pos.empty:
        return

    pos["label"] = pos["activity"].str[:4] + ": " + pos["metabolite"].str.slice(0, 40)
    counts = pos.groupby(["species", "label"]).size().reset_index(name="n")
    wide = counts.pivot(index="species", columns="label", values="n").fillna(0)
    if wide.shape[1] > top_n:
        top_cols = wide.sum(axis=0).sort_values(ascending=False).head(top_n).index
        wide = wide[top_cols]

    wide.index = [species_two_line_label(str(i)) for i in wide.index]
    fig_h = max(4, 0.45 * len(wide))
    fig, ax = plt.subplots(figsize=(min(14, 0.35 * wide.shape[1] + 4), fig_h))
    sns.heatmap(wide, cmap="YlOrRd", ax=ax, cbar_kws={"label": "BacDive API record"})
    ax.set_title("BacDive metabolite phenotypes (API, focal strains)")
    ax.set_ylabel("")
    plt.xticks(rotation=60, ha="right", fontsize=7)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_gene_hits_by_activity(hits: pd.DataFrame, out: Path) -> None:
    if hits.empty:
        return
    d = hits.groupby(["species", "activity"]).size().reset_index(name="n_hits")
    d["species"] = d["species"].map(species_two_line_label)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.4 * d["species"].nunique())))
    sns.barplot(data=d, y="species", x="n_hits", hue="activity", ax=ax, palette="Set2")
    ax.set_xlabel("Targeted gene / AMR hits")
    ax.set_title("Genomic search hits by phenotype category")
    ax.legend(title="Activity", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_composite_dashboard(
    focal: pd.DataFrame,
    integrated: pd.DataFrame,
    pheno: pd.DataFrame,
    hits: pd.DataFrame,
    out: Path,
) -> None:
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    # A: BGC vs AMR (inline)
    ax0 = fig.add_subplot(gs[0, 0])
    d = integrated.copy()
    for col in ("bgc_count_total", "amr_gene_count", "genomic_hit_n"):
        if col not in d.columns:
            d[col] = 0
        d[col] = pd.to_numeric(d[col], errors="coerce").fillna(0)
    y = np.arange(len(d))
    ax0.barh(y - 0.2, d["bgc_count_total"], height=0.35, color="#6A4C93", label="BGC")
    ax0.barh(y + 0.2, d["amr_gene_count"], height=0.35, color="#E45756", label="AMR")
    ax0.set_yticks(y)
    ax0.set_yticklabels([species_two_line_label(s) for s in d["species"]], fontsize=6)
    ax0.set_title("A  BGC & AMR counts")
    ax0.legend(fontsize=7)

    # B: phenotype counts from API
    ax1 = fig.add_subplot(gs[0, 1])
    if not pheno.empty and "bacdive_util_n" in integrated.columns and len(integrated) > 0:
        cols = ["bacdive_util_n", "bacdive_prod_n", "bacdive_resistance_n"]
        sub = integrated[["species"] + [c for c in cols if c in integrated.columns]].copy()
        for c in cols:
            sub[c] = pd.to_numeric(sub[c], errors="coerce").fillna(0)
        x = np.arange(len(sub))
        w = 0.25
        ax1.bar(x - w, sub.get("bacdive_util_n", 0), width=w, label="Util", color="#4C78A8")
        ax1.bar(x, sub.get("bacdive_prod_n", 0), width=w, label="Prod", color="#F58518")
        ax1.bar(x + w, sub.get("bacdive_resistance_n", 0), width=w, label="Res", color="#E45756")
        ax1.set_xticks(x)
        ax1.set_xticklabels(
            [species_two_line_label(s) for s in sub["species"]],
            rotation=45,
            ha="right",
            fontsize=6,
        )
        ax1.set_title("B  BacDive API metabolite records")
        ax1.legend(fontsize=7)

    # C: gene hits
    ax2 = fig.add_subplot(gs[1, 0])
    if not hits.empty:
        hc = hits.groupby("species").size().sort_values(ascending=True)
        ax2.barh(
            [species_two_line_label(s) for s in hc.index],
            hc.values,
            color="#54A24B",
        )
        ax2.set_title("C  Targeted genomic hits")
        ax2.set_xlabel("# hits")

    # D: confidence vs genomic score
    ax3 = fig.add_subplot(gs[1, 1])
    m = integrated.copy()
    for col in ("bgc_count_total", "amr_gene_count", "genomic_hit_n"):
        if col not in m.columns:
            m[col] = 0
        m[col] = pd.to_numeric(m[col], errors="coerce").fillna(0)
    m["genomic_score"] = m["bgc_count_total"] + m["amr_gene_count"] * 2 + m["genomic_hit_n"] * 0.5
    ax3.scatter(
        pd.to_numeric(m["match_confidence"], errors="coerce"),
        m["genomic_score"],
        s=60,
        c="#4C78A8",
        edgecolors="white",
    )
    for _, row in m.iterrows():
        ax3.annotate(
            species_two_line_label(row["species"]),
            (row["match_confidence"], row["genomic_score"]),
            fontsize=5,
            xytext=(3, 3),
            textcoords="offset points",
        )
    ax3.set_xlabel("Match confidence")
    ax3.set_ylabel("Genomic evidence score")
    ax3.set_title("D  Evidence strength")

    fig.suptitle("Genomic follow-up dashboard (focal species)", fontsize=13, y=0.98)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def fig_confidence_heatmap(conf: pd.DataFrame, activity: str, out: Path) -> None:
    """Species x metabolite heatmap of phenotype confidence for one activity.

    BacDive-observed cells are pinned to 1.0/0.0 and marked with a dot.
    """
    if conf.empty:
        return
    d = conf[conf["activity"] == activity]
    if d.empty:
        return

    wide = d.pivot_table(
        index="metabolite", columns="species", values="confidence", aggfunc="max"
    ).fillna(0.0)
    if wide.empty:
        return
    # order metabolites by total evidence so the strongest rows are on top
    wide = wide.loc[wide.sum(axis=1).sort_values(ascending=False).index]

    obs = d[d["observed"]]
    obs_wide = (
        obs.pivot_table(index="metabolite", columns="species", values="confidence", aggfunc="max")
        if not obs.empty
        else pd.DataFrame()
    )

    n_rows, n_cols = wide.shape
    fig_h = max(4.0, 0.35 * n_rows + 1.5)
    fig_w = max(6.0, 1.5 * n_cols + 3.0)
    fig, ax = plt.subplots(figsize=(min(18, fig_w), min(24, fig_h)))
    sns.heatmap(
        wide,
        cmap="RdYlGn",
        vmin=0.0,
        vmax=1.0,
        ax=ax,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Confidence (BacDive=observed, else genomic estimate)"},
    )
    ax.set_xticklabels(
        [species_two_line_label(t.get_text()) for t in ax.get_xticklabels()],
        fontsize=7,
        rotation=0,
    )
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

    # mark BacDive-observed cells
    if not obs_wide.empty:
        col_pos = {c: i for i, c in enumerate(wide.columns)}
        row_pos = {r: i for i, r in enumerate(wide.index)}
        for met in obs_wide.index:
            for sp in obs_wide.columns:
                val = obs_wide.loc[met, sp]
                if pd.notna(val) and met in row_pos and sp in col_pos:
                    ax.text(
                        col_pos[sp] + 0.5,
                        row_pos[met] + 0.5,
                        "\u25cf",
                        ha="center",
                        va="center",
                        fontsize=5,
                        color="black",
                    )

    ax.set_title(
        f"Predicted {activity} confidence per metabolite  (\u25cf = BacDive-tested)",
        fontsize=11,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Metabolite")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def generate_followup_figures(
    focal: pd.DataFrame,
    run_dir: Path,
    fig_dir: Path,
    *,
    phenotype_path: Optional[Path] = None,
    hits_path: Optional[Path] = None,
) -> Dict[str, str]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    integrated = _load_optional(run_dir / "integrated_genomic_followup.csv")
    pheno = _load_optional(phenotype_path or run_dir / "bacdive_phenotype_metabolites.csv")
    hits = _load_optional(hits_path or run_dir / "targeted_gene_hits.csv")

    paths: Dict[str, str] = {}
    p1 = fig_dir / "fig1_pipeline_status.png"
    fig_pipeline_status(integrated if not integrated.empty else focal, p1)
    paths["fig1"] = str(p1)

    p2 = fig_dir / "fig2_bgc_vs_amr.png"
    fig_bgc_and_amr(integrated if not integrated.empty else focal, p2)
    paths["fig2"] = str(p2)

    p3 = fig_dir / "fig3_phenotype_vs_genomic.png"
    fig_phenotype_breadth_vs_genomic(focal, integrated if not integrated.empty else focal, p3)
    paths["fig3"] = str(p3)

    if not pheno.empty:
        p4 = fig_dir / "fig4_bacdive_metabolite_heatmap.png"
        fig_metabolite_heatmap(pheno, p4)
        paths["fig4"] = str(p4)

    if not hits.empty:
        p5 = fig_dir / "fig5_gene_hits_by_activity.png"
        fig_gene_hits_by_activity(hits, p5)
        paths["fig5"] = str(p5)

    p6 = fig_dir / "fig6_composite_genomic_dashboard.png"
    fig_composite_dashboard(focal, integrated if not integrated.empty else focal, pheno, hits, p6)
    paths["fig6"] = str(p6)

    conf = _load_optional(run_dir / "phenotype_confidence.csv")
    if not conf.empty:
        for key, activity, fname in (
            ("fig7", "resistance", "fig7_confidence_resistance.png"),
            ("fig8", "production", "fig8_confidence_production.png"),
            ("fig9", "utilization", "fig9_confidence_utilization.png"),
        ):
            fp = fig_dir / fname
            fig_confidence_heatmap(conf, activity, fp)
            if fp.exists():
                paths[key] = str(fp)

    import shutil

    stable_dir = fig_dir.parent.parent
    for name in paths.values():
        shutil.copy2(name, fig_dir.parent / Path(name).name)
        if stable_dir.name == "genomic_followup":
            shutil.copy2(name, stable_dir / Path(name).name)

    return paths
