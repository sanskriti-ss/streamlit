import streamlit as st
import os
import pandas as pd
from tabs import general_overview, circos, trends, cards
from utils.data_loader import load_data

# App title
st.title("Visualizing BacDive")
st.write("Work in progress.")

# Load CSV files from the data_files folder
data_folder = "data_files"
data_frames = load_data(data_folder)

if data_frames:
    st.sidebar.success("CSV files loaded successfully.")

# Sidebar Navigation
st.sidebar.title("Navigation")
tab = st.sidebar.radio("Go to", ["General Overview", "Circos", "Trends", "Cards"])

# Route to the correct tab
if tab == "General Overview":
    general_overview.display(data_frames)
elif tab == "Circos":
    circos.display()
elif tab == "Trends":
    trends.display(data_frames)
elif tab == "Cards":
    cards.display()
   
