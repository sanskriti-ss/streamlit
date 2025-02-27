import streamlit as st
import textwrap
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO
from utils.visualization import create_heatmap_by_numbers, create_heatmap_by_proportions

def display(data_frames):
    st.title("General Overview")
    st.write("This section will show a summary of the antibiotic data.")

    # User choices
    test_type = st.radio("Choose Test Type:", ["Negatively Tested", "Positively Tested"])
    strain_option = st.radio("Include Strains?", ["No Strains", "Yes Strains"])
    category = st.selectbox("Select Category:", ["Production", "Utilization", "Resistance", "Sensitivity"])

    # Mapping selections to filename
    test_type_short = "negatively" if test_type == "Negatively Tested" else "positively"
    strain_short = "nostrain" if strain_option == "No Strains" else "yesstrain"
    category_mapping = {"Production": "prod", "Utilization": "util", "Resistance": "res", "Sensitivity": "sen"}
    category_short = category_mapping[category]

    file_key = f"step4_{test_type_short}_tested_by_genera_{category_short}_{strain_short}.csv"

    # Check if file exists
    if file_key in data_frames:
        selected_file = file_key
        st.write(f"Selected File: **{selected_file}**")
    else:
        st.error("File not found! Please check your selections.")
        return

    # Sliders
    num_genera = st.slider("Number of Genera:", 5, 20, 10)
    num_metabolites = st.slider("Number of Metabolites:", 5, 40, 20)

    # Generate heatmap by numbers
    st.subheader("By Numbers")
    if st.button("Generate Heatmap by Numbers"):
        buf, title = create_heatmap_by_numbers(data_frames[selected_file], num_genera, num_metabolites)
        st.download_button("Download Heatmap", buf, file_name=f"{title.replace(' ', '_')}.jpg", mime="image/jpeg")

    st.markdown("---")

    # Generate heatmap by proportions
    st.subheader("By Proportions")
    if st.button("Generate Heatmap by Proportions"):
        buf, title = create_heatmap_by_proportions(data_frames[selected_file], num_genera, num_metabolites)
        st.download_button("Download Heatmap", buf, file_name=f"{title.replace(' ', '_')}.jpg", mime="image/jpeg")
