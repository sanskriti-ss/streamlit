import streamlit as st
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import textwrap
from io import BytesIO
import pandas as pd

def create_heatmap_by_numbers(df, num_genera, num_metabolites, fmt="d"):
    """Creates a heatmap by numbers."""
    metabolite_columns = df.columns.difference(['genus', 'species_count']).tolist()
    df[metabolite_columns] = df[metabolite_columns].apply(pd.to_numeric, errors='coerce')

    top_metabolites = df[metabolite_columns].sum().nlargest(num_metabolites).index
    top_genera = df[metabolite_columns].sum(axis=1).nlargest(num_genera).index

    heatmap_data = df.iloc[top_genera][top_metabolites].set_index(df['genus'].iloc[top_genera])

    plt.figure(figsize=(12, 8))
    
    # Determine format automatically if not specified
    if fmt == "auto":
        fmt = "d" if np.all(heatmap_data.dropna().astype(int) == heatmap_data.dropna()) else ".2f"
    
    sns.heatmap(heatmap_data, annot=True, fmt=fmt, cmap="viridis", linewidths=.5)
    plt.title("Heatmap")
    
    buf = BytesIO()
    plt.savefig(buf, format='jpg')
    buf.seek(0)
    
    st.pyplot(plt)
    
    return buf, "Heatmap"

def create_heatmap_by_proportions(df, num_genera, num_metabolites):
    """Creates a heatmap by proportions."""
    df['species_count'] = pd.to_numeric(df['species_count'], errors='coerce')
    
    # Avoid division by zero by replacing NaN and inf values
    df.iloc[:, 1:] = df.iloc[:, 1:].div(df['species_count'].replace(0, np.nan), axis=0)
    df = df.fillna(0)  # Replace NaNs with zero after division
    
    return create_heatmap_by_numbers(df, num_genera, num_metabolites, fmt=".2f")
