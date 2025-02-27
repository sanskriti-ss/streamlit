import streamlit as st
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO

def summary_bargraphs(data_frames):
    st.header("Metabolite Counts Bar Graph")
    
    # User choices for test type and strain option
    test_type = st.radio("Choose Test Type:", ["Negatively Tested", "Positively Tested"], index=1)
    strain_option = st.radio("Include Strains?", ["No Strains", "Yes Strains"], index=0)
    
    test_type_short = "negatively" if test_type == "Negatively Tested" else "positively"
    strain_short = "nostrain" if strain_option == "No Strains" else "yesstrain"
    
    # Mapping for the four categories
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }
    
    # For each category, load the corresponding file and compute unique metabolite counts per genus.
    category_counts = {}
    all_genera = set()
    
    for cat, cat_short in category_mapping.items():
        file_key = f"step4_{test_type_short}_tested_by_genera_{cat_short}_{strain_short}.csv"
        if file_key in data_frames:
            df = data_frames[file_key]
            # Exclude the non-metabolite columns
            metabolite_cols = df.columns.difference(['genus', 'species_count'])
            # Instead of summing, count each metabolite as 1 if its value is > 0
            df['unique_metabolites'] = (df[metabolite_cols] > 0).sum(axis=1)
            counts_series = df.set_index('genus')['unique_metabolites']
            category_counts[cat] = counts_series
            all_genera.update(counts_series.index.tolist())
        else:
            st.error(f"File not found: {file_key}")
            return
    
    # Combine data from all categories into a single DataFrame.
    combined = pd.DataFrame(index=list(all_genera))
    for cat in category_mapping.keys():
        # Some genera might be missing from a category; fill those with 0.
        combined[cat] = category_counts.get(cat, pd.Series()).reindex(combined.index, fill_value=0)
    
    # Compute the total unique metabolites across all categories for ranking
    # (i.e. if a metabolite is present in multiple categories, it is counted in each)
    combined['total_all'] = combined.sum(axis=1)
    top20 = combined.sort_values('total_all', ascending=False).head(20)
    
    # Plotting the overlayed bar chart.
    labels = top20.index.tolist()
    x = np.arange(len(labels))
    width = 0.8  # same x-position for each genus
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Colors for the four categories with translucency
    colors = {
        "Production": "#4c72b0",
        "Utilization": "#55a868",
        "Resistance": "#c44e52",
        "Sensitivity": "#8172b2"
    }
    
    # Plot each category as an overlapping bar at the same x positions.
    for cat in category_mapping.keys():
        ax.bar(x, top20[cat], width, 
               color=colors[cat], alpha=0.6, label=cat)
    
    ax.set_xlabel("Genera")
    ax.set_ylabel("Unique Metabolite Count")
    ax.set_title(f"Top 20 Genera by Unique Metabolite Counts ({test_type}, {strain_option})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.legend()
    
    st.pyplot(fig)

def display(data_frames):
    st.header("By the Numbers")
    st.write("Summary statistics:")
    
    # Choose a file that represents 'No Strain'
    file_nostrain = next((fname for fname in data_frames if "nostrain" in fname), None)
    # Choose a file that represents 'Yes Strain'
    file_yesstrain = next((fname for fname in data_frames if "yesstrain" in fname), None)
    
    if not file_nostrain:
        st.error("No file with 'nostrain' found in the data.")
        return
    
    # Use the 'No Strain' file for common summaries.
    df_nostrain = data_frames[file_nostrain]
    
    # Count metabolite columns (all columns except 'genus' and 'species_count')
    metabolite_columns = df_nostrain.columns.difference(['genus', 'species_count'])
    num_metabolites = len(metabolite_columns)
    
    # Number of unique genera
    num_genera = df_nostrain["genus"].nunique()
    
    # Total species count for 'No Strain'
    total_species_nostrain = df_nostrain["species_count"].sum()
    
    # If available, get total species count for 'Yes Strain'
    if file_yesstrain:
        df_yesstrain = data_frames[file_yesstrain]
        total_species_yesstrain = df_yesstrain["species_count"].sum()
    else:
        total_species_yesstrain = "Not Available"
    
    st.write(f"**Number of Metabolites:** {num_metabolites}")
    st.write(f"**Number of Genera:** {num_genera}")
    st.write(f"**Total Number of (No Strain) Species:** {total_species_nostrain}")
    st.write(f"**Total Number (Yes Strain) Species:** {total_species_yesstrain}")
    
    # Draw the overlayed bar chart below the summaries.
    summary_bargraphs(data_frames)
