import streamlit as st
st.set_page_config(page_title="Visualizing BacDive")

import os
import pandas as pd
from tabs import (
    general_overview,
    circos,
    circos_viewer,
    trends,
    cards2,
    by_the_numbers,
    comparison,
    species_analysis,
    genome_bgc_evidence,
    symbiosis_network,
)
from utils.data_loader import load_data

# Define available tabs
available_tabs = [
    "General Overview",
    "Circos",
    "Circle graph",
    "Trends",
    "Cards",
    "By the Numbers",
    "Comparison",
    "Species Analysis",
    "Genome and BGC Evidence",
    "Symbiosis Network",
]

# Load CSV files from the data_files folder (do this before sidebar to avoid reloading)
data_folder = "data_files"
data_frames = load_data(data_folder)

# Sidebar Navigation - must be at the top before any content
st.sidebar.title("Navigation")

# Initialize tab selection in session state if not present
if 'selected_tab' not in st.session_state:
    # Check query parameters for persisted tab
    query_tab = st.query_params.get("tab", None)
    if query_tab and query_tab in available_tabs:
        st.session_state.selected_tab = query_tab
    else:
        st.session_state.selected_tab = "General Overview"

# Create radio button with on_change callback
def on_tab_change():
    st.session_state.selected_tab = st.session_state.tab_radio
    # Update query params when tab changes
    st.query_params["tab"] = st.session_state.tab_radio

tab = st.sidebar.radio(
    "Go to", 
    available_tabs,
    key="tab_radio",
    index=available_tabs.index(st.session_state.selected_tab),
    on_change=on_tab_change
)

if data_frames:
    st.sidebar.success("CSV files loaded successfully.")

# Route to the correct tab based on session state
if st.session_state.selected_tab == "General Overview":
    general_overview.display(data_frames)
elif st.session_state.selected_tab == "Circos":
    circos.display()
elif st.session_state.selected_tab == "Circle graph":
    circos_viewer.display()
elif st.session_state.selected_tab == "Trends":
    trends.display(data_frames)
elif st.session_state.selected_tab == "Cards":
    cards2.display(data_frames)
elif st.session_state.selected_tab == "By the Numbers":
    by_the_numbers.display(data_frames)
elif st.session_state.selected_tab == "Comparison":
    comparison.display(data_frames)
elif st.session_state.selected_tab == "Species Analysis":
    species_analysis.display(data_frames)
elif st.session_state.selected_tab == "Genome and BGC Evidence":
    genome_bgc_evidence.display(data_frames)
elif st.session_state.selected_tab == "Symbiosis Network":
    symbiosis_network.display(data_frames)
