import streamlit as st
st.set_page_config(page_title="Visualizing BacDive")

import os
import pandas as pd
from tabs import general_overview, circos, trends, cards, by_the_numbers, comparison
from utils.data_loader import load_data



# Load CSV files from the data_files folder
data_folder = "data_files"
data_frames = load_data(data_folder)

if data_frames:
    st.sidebar.success("CSV files loaded successfully.")

# Sidebar Navigation
st.sidebar.title("Navigation")
tab = st.sidebar.radio("Go to", ["General Overview", "Circos", "Trends", "Cards", "By the Numbers", "Comparison"])

# Route to the correct tab
if tab == "General Overview":
    general_overview.display(data_frames)
elif tab == "Circos":
    circos.display()
elif tab == "Trends":
    trends.display(data_frames)
elif tab == "Cards":
    cards.display(data_frames)
elif tab == "By the Numbers":
    by_the_numbers.display(data_frames)
elif tab == "Comparison":
    comparison.display(data_frames)
