import streamlit as st
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import sys
from pathlib import Path

# Add parent directory to path to import utils
sys.path.append(str(Path(__file__).parent.parent))
from utils.tooltip_title import display_title_with_tooltip

def summary_bargraphs(data_frames):
    display_title_with_tooltip(
        "Metabolite Counts Bar Graph",
        sample_image_filename="bythenumbers_metcount.png"
    )
    
    # Widget commands (outside cached functions)
    test_type = st.radio("Choose Test Type:", ["Negatively Tested", "Positively Tested"], index=1)
    strain_option = st.radio("Include Strains?", ["Isolate", "Strain"], index=0)
    sort_category = st.selectbox(
        "Select category to sort top 15 by:",
        ["Production", "Utilization", "Resistance", "Sensitivity"],
        index=0
    )
    
    # Override option: let user choose specific genera instead of auto top 15.
    override = st.checkbox("Override Top 15 Selection (choose specific genera)")
    
    # Compute the combined DataFrame to extract available genera.
    combined = _compute_combined_df(data_frames, test_type, strain_option)
    available_genera = sorted(list(combined.index))
    
    override_genera = None
    if override:
        override_genera = st.multiselect(
            "Select up to fifteen genera to plot:",
            options=available_genera,
            default=[]
        )
        if len(override_genera) > 15:
            st.warning("Please select no more than fifteen genera. Only the first 15 will be used.")
            override_genera = override_genera[:15]
    
    # Call the cached computation/plotting function with the widget selections as parameters.
    fig = _compute_summary_bargraphs(data_frames, test_type, strain_option, sort_category, override_genera)
    st.pyplot(fig)

@st.cache_data
def _compute_combined_df(data_frames, test_type, strain_option):
    # Process widget input values
    test_type_short = "negatively" if test_type == "Negatively Tested" else "positively"
    strain_short = "isolate" if strain_option == "Isolate" else "strain"
    
    # Mapping for the four categories
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }
    
    category_counts = {}
    all_genera = set()
    
    for cat, cat_short in category_mapping.items():
        file_key = f"step4_{test_type_short}_tested_by_genera_{cat_short}_{strain_short}.csv"
        if file_key in data_frames:
            df = data_frames[file_key].copy()
            # Exclude the non-metabolite columns.
            metabolite_cols = df.columns.difference(['genus', 'species_count'])
            # Count each metabolite as 1 if its value is > 0.
            df['unique_metabolites'] = (df[metabolite_cols] > 0).sum(axis=1)
            counts_series = df.set_index('genus')['unique_metabolites']
            category_counts[cat] = counts_series
            all_genera.update(counts_series.index.tolist())
        else:
            raise FileNotFoundError(f"File not found: {file_key}")
    
    combined = pd.DataFrame(index=list(all_genera))
    for cat in category_mapping.keys():
        combined[cat] = category_counts.get(cat, pd.Series()).reindex(combined.index, fill_value=0)
    return combined

@st.cache_data
def _compute_summary_bargraphs(data_frames, test_type, strain_option, sort_category, override_genera):
    # Retrieve the combined DataFrame from the cached helper.
    combined = _compute_combined_df(data_frames, test_type, strain_option)
    
    # Use override genera if provided; otherwise, select the top 15 by sort_category.
    if override_genera is not None and len(override_genera) > 0:
        top_selection = combined.loc[combined.index.intersection(override_genera)]
        top_selection = top_selection.sort_values(sort_category, ascending=False)
    else:
        top_selection = combined.sort_values(sort_category, ascending=False).head(15)
    
    labels = top_selection.index.tolist()
    x = np.arange(len(labels))
    n_categories = len(combined.columns)  # Should be 4
    total_width = 0.8  # Total width for all bars at one x value
    bar_width = total_width / n_categories
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Colors for the four categories.
    colors = {
        "Production": "#4c72b0",
        "Utilization": "#55a868",
        "Resistance": "#c44e52",
        "Sensitivity": "#8172b2"
    }
    
    # Plot each category as a side-by-side bar.
    for i, cat in enumerate(["Production", "Utilization", "Resistance", "Sensitivity"]):
        offset = -total_width/2 + i*bar_width + bar_width/2
        ax.bar(x + offset, top_selection[cat], bar_width, 
               color=colors[cat], alpha=0.8, label=cat)
    
    ax.set_xlabel("Genera")
    ax.set_ylabel("Unique Metabolite Count")
    ax.set_title(f"Genera by Unique Metabolite Counts (Sorted by {sort_category})\n({test_type}, {strain_option})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.legend()
    
    return fig

def display(data_frames):
    display_title_with_tooltip(
        "By the Numbers",
        description_text="Summary statistics:"
    )
    
    # Choose a file that represents 'Isolate' and 'Strain'.
    file_isolate = next((fname for fname in data_frames if "isolate" in fname), None)
    file_strain = next((fname for fname in data_frames if "strain" in fname and "isolate" not in fname), None)
    
    if not file_isolate:
        st.error("No file with 'isolate' found in the data.")
        return
    
    df_isolate = data_frames[file_isolate]
    metabolite_columns = df_isolate.columns.difference(['genus', 'species_count'])
    num_metabolites = len(metabolite_columns)
    num_genera = df_isolate["genus"].nunique()
    total_species_isolate = df_isolate["species_count"].sum()
    
    if file_strain:
        df_strain = data_frames[file_strain]
        total_species_strain = df_strain["species_count"].sum()
    else:
        total_species_strain = "Not Available"
    
    st.write(f"**Number of Metabolites:** {num_metabolites}")
    st.write(f"**Number of Genera:** {num_genera}")
    st.write(f"**Total Number of (Isolate) Species:** {total_species_isolate}")
    st.write(f"**Total Number (Strain) Species:** {total_species_strain}")
    
    # Draw the grouped bar chart below the summaries.
    summary_bargraphs(data_frames)

    
    ###########################################
    # New Section: Homogeneous Metabolite Summary
    ###########################################
    st.markdown("---")
    display_title_with_tooltip(
        "Homogeneous Metabolite Summary by Genus",
        sample_image_filename="bythenumbers_metsummary.png",
        description_text=(
            "For a selected metabolite category, this section lists genera where, for at least one metabolite, "
            "all species (with a minimum of 5 species) tested uniformly positive or uniformly negative. "
            "It also shows which metabolite(s) met that criteria. If a genus appears with both positive and negative "
            "results (i.e. mixed), it is omitted."
        )
    )
    
    # New selection options for metabolite category and strain.
    homo_category = st.selectbox(
        "Select Metabolite Category:",
        ["Production", "Utilization", "Resistance", "Sensitivity"],
        index=0
    )
    homo_strain = st.selectbox(
        "Select Strain Option:",
        ["Isolate", "Strain"],
        index=0
    )
    
    if st.button("Generate Homogeneous Summary"):
        summary_df = _compute_homogeneous_summary(data_frames, homo_category, homo_strain)
        if summary_df.empty:
            st.write("No genera found meeting the criteria.")
        else:
            st.dataframe(summary_df)


@st.cache_data
def _compute_homogeneous_summary(data_frames, homo_category, homo_strain):
    """
    For the selected metabolite category and strain option, this function examines both the 
    positively_tested and negatively_tested CSV files. For each genus (with species_count >= 5),
    it checks each metabolite column (excluding 'genus' and 'species_count') and records the metabolite(s)
    for which the test is homogeneous (i.e. in the positive file, the value equals species_count; in the
    negative file, the value equals species_count). If a genus qualifies in both positive and negative, it is omitted.
    
    Returns a DataFrame indexed by genus with the following columns:
        - Result: "Positive" or "Negative"
        - Metabolites: a comma-separated list of metabolite names for which the homogeneous result is observed.
        - Species Count: the species_count for that genus.
    """
    cat_mapping = {"Production": "prod", "Utilization": "util", "Resistance": "res", "Sensitivity": "sen"}
    if homo_category not in cat_mapping:
        raise ValueError("Invalid metabolite category selection.")
    cat_short = cat_mapping[homo_category]
    strain_short = "isolate" if homo_strain == "Isolate" else "strain"
    
    # Construct file keys.
    file_key_pos = f"step4_positively_tested_by_genera_{cat_short}_{strain_short}.csv"
    file_key_neg = f"step4_negatively_tested_by_genera_{cat_short}_{strain_short}.csv"
    
    pos_df = data_frames.get(file_key_pos, None)
    neg_df = data_frames.get(file_key_neg, None)
    
    pos_summary = {}
    neg_summary = {}
    
    # Process positive file.
    if pos_df is not None:
        # Assume columns: 'genus', 'species_count', plus one column per metabolite.
        for idx, row in pos_df.iterrows():
            genus = row["genus"]
            species_count = row["species_count"]
            if species_count < 5:
                continue
            # For each metabolite column, check if value equals species_count.
            homogeneous_metabolites = []
            for col in pos_df.columns.difference(["genus", "species_count"]):
                if row[col] == species_count:
                    homogeneous_metabolites.append(col)
            if homogeneous_metabolites:
                pos_summary[genus] = {
                    "Result": "Positive",
                    "Metabolites": ", ".join(homogeneous_metabolites),
                    "Species Count": species_count
                }
    
    # Process negative file.
    if neg_df is not None:
        for idx, row in neg_df.iterrows():
            genus = row["genus"]
            species_count = row["species_count"]
            if species_count < 5:
                continue
            homogeneous_metabolites = []
            for col in neg_df.columns.difference(["genus", "species_count"]):
                if row[col] == species_count:
                    homogeneous_metabolites.append(col)
            if homogeneous_metabolites:
                neg_summary[genus] = {
                    "Result": "Negative",
                    "Metabolites": ", ".join(homogeneous_metabolites),
                    "Species Count": species_count
                }
    
    # Remove genera that appear in both summaries (mixed results).
    mixed = set(pos_summary.keys()) & set(neg_summary.keys())
    for genus in mixed:
        del pos_summary[genus]
        del neg_summary[genus]
    
    # Combine the remaining dictionaries.
    final_summary = {**pos_summary, **neg_summary}
    summary_df = pd.DataFrame.from_dict(final_summary, orient='index')
    if not summary_df.empty:
        summary_df = summary_df.sort_values(by="Species Count", ascending=False)
    return summary_df
