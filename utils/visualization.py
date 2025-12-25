import streamlit as st
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import textwrap
from io import BytesIO
import pandas as pd

def create_heatmap_by_numbers(df, num_genera, num_metabolites, fmt="d", test_type=None, strain_option=None, category=None):
    """Creates a heatmap by numbers, filtering out all-zero columns."""
    metabolite_columns = df.columns.difference(['genus', 'species_count']).tolist()
    df[metabolite_columns] = df[metabolite_columns].apply(pd.to_numeric, errors='coerce')

    # Select top genera first
    top_genera = df[metabolite_columns].sum(axis=1).nlargest(num_genera).index
    df_top_genera = df.iloc[top_genera]
    
    # Filter metabolites: only keep those with at least one non-zero value in the selected genera
    non_zero_metabolites = df_top_genera[metabolite_columns].columns[(df_top_genera[metabolite_columns] > 0).any()]
    
    # If we have more non-zero metabolites than requested, select top N by sum
    if len(non_zero_metabolites) > num_metabolites:
        top_metabolites = df_top_genera[non_zero_metabolites].sum().nlargest(num_metabolites).index
    else:
        # Use all non-zero metabolites if we have fewer than requested
        top_metabolites = non_zero_metabolites

    heatmap_data = df_top_genera[top_metabolites].set_index(df['genus'].iloc[top_genera])

    plt.figure(figsize=(12, 8))
    
    # Determine format automatically if not specified
    if fmt == "auto":
        fmt = "d" if np.all(heatmap_data.dropna().astype(int) == heatmap_data.dropna()) else ".2f"
    
    # Build descriptive title
    if test_type and strain_option and category:
        title = f"{test_type} - {strain_option} - {category} - Showing {len(top_metabolites)} metabolites"
    else:
        title = f"Heatmap - Showing {len(top_metabolites)} metabolites"
    
    sns.heatmap(heatmap_data, annot=True, fmt=fmt, cmap="viridis", linewidths=.5)
    plt.title(title)
    
    buf = BytesIO()
    plt.savefig(buf, format='jpg')
    buf.seek(0)
    
    st.pyplot(plt)
    
    return buf, title

def create_heatmap_by_proportions(df, num_genera, num_metabolites, test_type=None, strain_option=None, category=None):
    """Creates a heatmap by proportions."""
    df['species_count'] = pd.to_numeric(df['species_count'], errors='coerce')
    
    # Avoid division by zero by replacing NaN and inf values
    df.iloc[:, 1:] = df.iloc[:, 1:].div(df['species_count'].replace(0, np.nan), axis=0)
    df = df.fillna(0)  # Replace NaNs with zero after division
    
    return create_heatmap_by_numbers(df, num_genera, num_metabolites, fmt=".2f", test_type=test_type, strain_option=strain_option, category=category)
