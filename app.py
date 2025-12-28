import streamlit as st
st.set_page_config(page_title="Visualizing BacDive")

import os
import pandas as pd
from tabs import general_overview, circos, trends, cards2, by_the_numbers, comparison, species_analysis
from utils.data_loader import load_data

# Define available tabs
available_tabs = ["General Overview", "Circos", "Trends", "Cards", "By the Numbers", "Comparison", "Species Analysis"]

# Load CSV files from the data_files folder
data_folder = "data_files"
data_frames = load_data(data_folder)

if data_frames:
    st.sidebar.success("CSV files loaded successfully.")

# Sidebar Navigation
st.sidebar.title("Navigation")

# Get saved tab from query params for initial load only
query_tab = st.query_params.get("tab", None)
if query_tab in available_tabs:
    default_idx = available_tabs.index(query_tab)
else:
    default_idx = 0

# Radio button with default index set once at load
tab = st.sidebar.radio(
    "Go to", 
    available_tabs,
    index=default_idx
)

# Save the selected tab to query params for persistence
st.query_params["tab"] = tab

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
