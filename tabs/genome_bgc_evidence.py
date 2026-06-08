"""
Streamlit: Genome and BGC Evidence (optional genome investigation layer).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOME_CSV = REPO_ROOT / "genome_investigation" / "results" / "Step2_5_genome_enriched.csv"
DEFAULT_ANTISMASH_CSV = REPO_ROOT / "genome_investigation" / "results" / "antismash_summary.csv"
DEFAULT_AMR_CSV = REPO_ROOT / "genome_investigation" / "results" / "amrfinder_summary.csv"
DEFAULT_RANKED_CSV = REPO_ROOT / "genome_investigation" / "results" / "ranked_species_candidates.csv"
FOLLOWUP_DIR = REPO_ROOT / "genome_investigation" / "results" / "genomic_followup"
PAPER_DIR = REPO_ROOT / "genome_investigation" / "results" / "paper"
SPECIES_DATA = REPO_ROOT / "species_data"


@st.cache_data(show_spinner=False)
def _load_csv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p)


@st.cache_data(show_spinner=False)
def _bacdive_phenotype_summary(species_name: str) -> dict:
    """Breadth summary from local Step3 tables (no API calls in Streamlit)."""
    import zipfile

    summary = {}
    files = {
        "utilization": SPECIES_DATA / "step3_met_util_exploded.csv.zip",
        "production": SPECIES_DATA / "step3_met_prod_exploded.csv.zip",
        "resistance": SPECIES_DATA / "step3_met_res_exploded.csv.zip",
        "sensitivity": SPECIES_DATA / "step3_met_sen_exploded.csv.zip",
    }
    meta = {"BacID", "species", "genus", "order", "type_strain", "is_strain", "species_with_id"}
    for activity, fpath in files.items():
        if not fpath.exists():
            continue
        with zipfile.ZipFile(fpath, "r") as zf:
            csvs = [n for n in zf.namelist() if n.endswith(".csv")]
            df = pd.read_csv(zf.open(csvs[0]), low_memory=False)
        sub = df[df["species"].astype(str) == species_name]
        if sub.empty:
            summary[activity] = {"breadth": 0, "n_strains": 0}
            continue
        mets = [c for c in df.columns if c not in meta]
        x = sub[mets].apply(pd.to_numeric, errors="coerce")
        if activity in ("resistance", "sensitivity"):
            breadth = int((x == 1).sum(axis=1).max())
        else:
            breadth = int((x.replace(-1, 0).fillna(0) > 0).sum(axis=1).max())
        summary[activity] = {"breadth": breadth, "n_strains": len(sub)}
    return summary


def display(_data_frames=None) -> None:
    st.title("Genome and BGC Evidence")
    st.caption(
        "Optional genome metadata and biosynthetic gene cluster (BGC) evidence. "
        "Run `genome_investigation` CLI tools to generate the CSVs below."
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        genome_path = st.text_input("Step2_5 genome enriched CSV", value=str(DEFAULT_GENOME_CSV))
    with col2:
        asm_path = st.text_input("antiSMASH summary CSV", value=str(DEFAULT_ANTISMASH_CSV))
    with col3:
        ranked_path = st.text_input("Ranked candidates CSV", value=str(DEFAULT_RANKED_CSV))

    genome_df = _load_csv(genome_path)
    asm_df = _load_csv(asm_path)
    ranked_df = _load_csv(ranked_path)

    if genome_df.empty:
        st.warning(
            "No genome enrichment file found. Generate one with:\n\n"
            "`python -m genome_investigation.genome_enrichment --input genome_investigation/results/selected_test_input.csv "
            "--output genome_investigation/results/Step2_5_genome_enriched.csv`"
        )
        return

    # Overview visualizations (all species in enrichment file)
    st.subheader("Overview")
    if PAPER_DIR.exists() and (PAPER_DIR / "RESULTS_INTERPRETATION.md").exists():
        with st.expander("Results interpretation (markdown)", expanded=False):
            st.markdown((PAPER_DIR / "RESULTS_INTERPRETATION.md").read_text(encoding="utf-8"))
        composite = PAPER_DIR / "fig4_composite_integrated_dashboard.png"
        if composite.exists():
            st.image(str(composite), use_container_width=True, caption="Integrated dashboard (Table 1)")
        with st.expander("Individual figures", expanded=False):
            c1, c2, c3 = st.columns(3)
            for col, name in zip(
                [c1, c2, c3],
                [
                    "fig1_genome_match_confidence.png",
                    "fig2_phenotype_vs_genome_confidence.png",
                    "fig3_bacdive_breadth_heatmap.png",
                ],
            ):
                p = PAPER_DIR / name
                if p.exists():
                    col.image(str(p), use_container_width=True)
    else:
        st.info(
            "Generate paper-ready figures: "
            "`python -m genome_investigation.visualize_results`"
        )
        chart = genome_df.sort_values("match_confidence", ascending=False)
        st.bar_chart(chart.set_index("species")["match_confidence"])

    table_path = PAPER_DIR / "table1_integrated_summary.csv"
    if table_path.exists():
        st.dataframe(pd.read_csv(table_path), use_container_width=True, hide_index=True)

    st.subheader("Genomic follow-up pipeline")
    st.caption(
        "Run: `python -m genome_investigation.genome_followup_pipeline` "
        "(download, antiSMASH, AMRFinder, BacDive metabolites, gene search, figures)."
    )
    followup_fig = FOLLOWUP_DIR / "fig6_composite_genomic_dashboard.png"
    if not followup_fig.exists():
        runs = sorted(FOLLOWUP_DIR.glob("run_*/figures/fig6_composite_genomic_dashboard.png"))
        followup_fig = runs[-1] if runs else followup_fig
    if followup_fig.exists():
        st.image(str(followup_fig), use_container_width=True, caption="Genomic follow-up composite")
    integrated_path = FOLLOWUP_DIR / "integrated_genomic_followup.csv"
    if integrated_path.exists():
        st.dataframe(pd.read_csv(integrated_path), use_container_width=True, hide_index=True)
    amr_df = _load_csv(str(DEFAULT_AMR_CSV))
    if not amr_df.empty:
        st.write("**AMRFinder summary**")
        st.dataframe(amr_df, use_container_width=True, hide_index=True)

    st.divider()
    species_list = sorted(genome_df["species"].dropna().astype(str).unique())
    default_ix = 0
    if ranked_df is not None and not ranked_df.empty:
        top_sp = str(ranked_df.iloc[0]["species"])
        if top_sp in species_list:
            default_ix = species_list.index(top_sp)

    species = st.selectbox("Species", species_list, index=default_ix)

    g_rows = genome_df[genome_df["species"].astype(str) == species]
    if g_rows.empty:
        st.error("No genome row for selected species")
        return
    g = g_rows.sort_values("match_confidence", ascending=False).iloc[0]

    st.subheader("1. BacDive phenotype summary")
    pheno = _bacdive_phenotype_summary(species)
    if pheno:
        st.dataframe(pd.DataFrame(pheno).T, use_container_width=True)
    else:
        st.info("No local Step3 phenotype rows for this species.")

    st.subheader("2. Genome metadata")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.metric("Accession", g.get("genome_accession", "—"))
        st.metric("Assembly level", g.get("assembly_level", "—"))
    with c2:
        st.metric("Genome size (bp)", g.get("genome_size_bp", "—"))
        st.metric("GC %", g.get("gc_percent", "—"))
    with c3:
        st.metric("Gene count", g.get("gene_count", "—"))
        st.metric("Source", g.get("source_database", "—"))

    st.subheader("3. Genome match confidence")
    conf = float(g.get("match_confidence") or 0)
    st.progress(min(1.0, max(0.0, conf)))
    st.write(f"**Confidence:** {conf:.3f}")
    st.write(f"**Notes:** {g.get('match_notes', '')}")
    if conf < 0.7:
        st.warning("Low-confidence genome match — interpret genomic evidence cautiously.")

    asm_row = None
    if not asm_df.empty:
        asm_matches = asm_df[asm_df["species"].astype(str) == species]
        if not asm_matches.empty:
            asm_row = asm_matches.iloc[0]

    amr_row = None
    if not amr_df.empty:
        amr_matches = amr_df[amr_df["species"].astype(str) == species]
        if not amr_matches.empty:
            amr_row = amr_matches.iloc[0]

    st.subheader("4. AMRFinder (resistance genes)")
    if amr_row is None:
        st.info("No AMRFinder summary. Run `genome_followup_pipeline` after downloading genomes.")
    else:
        st.metric("AMR genes", int(amr_row.get("amr_gene_count") or 0))
        st.write(f"**Classes:** {amr_row.get('amr_classes', '—')}")
        st.write(f"**Genes:** {amr_row.get('amr_genes', '—')}")

    st.subheader("5. antiSMASH BGC summary")
    if asm_row is None:
        st.info("No antiSMASH summary for this species. Run antiSMASH on selected genomes only (see README).")
    else:
        st.metric("Total BGCs", int(asm_row.get("bgc_count_total") or 0))
        st.write(f"**Status:** {asm_row.get('antismash_status', '')}")
        notes = str(asm_row.get("antismash_notes", ""))
        if "fasta" in notes.lower():
            st.warning("antiSMASH was run on FASTA — annotation quality may be lower than GenBank/GBFF input.")
        st.write(f"**Notes:** {notes}")

        st.subheader("6. BGC type counts")
        type_cols = [
            ("NRPS", "nrps_count"),
            ("PKS", "pks_count"),
            ("Terpene", "terpene_count"),
            ("RiPP / ribosomal", "ribosomal_peptide_count"),
            ("Saccharide", "saccharide_count"),
            ("Siderophore-related", "siderophore_related_count"),
            ("Other", "other_bgc_count"),
        ]
        counts = {label: int(asm_row.get(col) or 0) for label, col in type_cols}
        st.bar_chart(pd.Series(counts))

        st.subheader("7. Known-cluster hits")
        st.write(asm_row.get("knownclusterblast_hits", "—"))
        st.write(asm_row.get("most_similar_known_clusters", "—"))

    st.subheader("8. Prioritization rationale")
    if not ranked_df.empty:
        r = ranked_df[ranked_df["species"].astype(str) == species]
        if not r.empty:
            for _, row in r.iterrows():
                st.markdown(f"**{row.get('category', '')}** — priority {row.get('priority_score', '')}")
                st.write(row.get("rationale", ""))
                if row.get("safety_flag"):
                    st.error(str(row["safety_flag"]))
        else:
            st.info("Species not in ranked candidates table.")
    else:
        st.info("Generate ranked_species_candidates.csv with species_prioritization.py")
