import streamlit as st
st.set_page_config(page_title="Visualizing BacDive")

import os
import pandas as pd
from tabs import general_overview, circos, trends, cards2, by_the_numbers, comparison, species_analysis
from utils.data_loader import load_data

# Initialize session state for tab persistence
if 'active_tab' not in st.session_state:
    st.session_state.active_tab = "General Overview"

# Load CSV files from the data_files folder
data_folder = "data_files"
data_frames = load_data(data_folder)

if data_frames:
    st.sidebar.success("CSV files loaded successfully.")

# Sidebar Navigation - use session state for default value
st.sidebar.title("Navigation")
tab = st.sidebar.radio(
    "Go to", 
    ["General Overview", "Circos", "Trends", "Cards", "By the Numbers", "Comparison", "Species Analysis"],
    index=["General Overview", "Circos", "Trends", "Cards", "By the Numbers", "Comparison", "Species Analysis"].index(st.session_state.active_tab),
    key="tab_selector"
)

# Update session state when tab changes
if tab != st.session_state.active_tab:
    st.session_state.active_tab = tab

# Route to the correct tab
if tab == "General Overview":
    general_overview.display(data_frames)
elif tab == "Circos":
    circos.display()
elif tab == "Trends":
    trends.display(data_frames)
elif tab == "Cards":
    cards2.display(data_frames)
elif tab == "By the Numbers":
    by_the_numbers.display(data_frames)
elif tab == "Comparison":
    comparison.display(data_frames)
elif tab == "Species Analysis":
    species_analysis.display(data_frames)
