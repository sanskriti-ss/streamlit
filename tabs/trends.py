import streamlit as st
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import plotly.graph_objects as go
from io import BytesIO

def display(data_frames):
    st.title("Parallel Coordinates Plot for Genera Trends")
    st.write("Compare up to 5 genera across Production, Utilization, Resistance, and Sensitivity.")

    # Step 1: Store selections in session state to prevent auto-re-running
    if "test_type" not in st.session_state:
        st.session_state.test_type = "Positively Tested"
    if "strain_option" not in st.session_state:
        st.session_state.strain_option = "No Strains"
    if "selected_genera" not in st.session_state:
        st.session_state.selected_genera = []

    # Step 2: User selections (prevent re-run on change)
    st.session_state.test_type = st.radio(
        "Choose Test Type:", ["Negatively Tested", "Positively Tested"],
        key="test_type_radio"
    )
    st.session_state.strain_option = st.radio(
        "Include Strains?", ["No Strains", "Yes Strains"],
        key="strain_option_radio"
    )

    # Step 3: Map selections to filenames
    test_type_short = "negatively" if st.session_state.test_type == "Negatively Tested" else "positively"
    strain_short = "nostrain" if st.session_state.strain_option == "No Strains" else "yesstrain"
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }

    # Step 4: Load the four relevant files
    relevant_files = {
        category: f"step4_{test_type_short}_tested_by_genera_{short}_{strain_short}.csv"
        for category, short in category_mapping.items()
    }

    # Step 5: Check if all files exist in data
    missing_files = [f for f in relevant_files.values() if f not in data_frames]
    if missing_files:
        st.error(f"Missing files: {missing_files}")
        return

    # Step 6: Extract the genus names (assuming they exist in all files)
    sample_df = data_frames[relevant_files["Production"]]
    genus_list = sample_df["genus"].tolist()

    # Step 7: Multi-selection for genera (up to 5) without triggering re-runs
    st.session_state.selected_genera = st.multiselect(
        "Search and Select Up to 5 Genera:",
        genus_list,
        default=st.session_state.selected_genera,
        key="selected_genera_multiselect"
    )

    # Step 8: Button to trigger execution of the parallel coordinates plot
    if st.button("Generate Parallel Coordinates Diagram"):
        categories = ["Production", "Utilization", "Resistance", "Sensitivity"]
        genus_values = {}

        # Extract data from all four files for selected genera
        for genus in st.session_state.selected_genera:
            values = []
            for category in categories:
                file_name = relevant_files[category]
                df = data_frames[file_name]
                # Find the selected genus row
                row = df[df["genus"] == genus]
                if not row.empty:
                    # Sum all chemical presence counts (assuming columns 2 onward)
                    values.append(row.iloc[0, 1:].sum())
                else:
                    values.append(0)
            genus_values[genus] = values

        # Step 9: Generate the Parallel Coordinates Diagram
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["b", "r", "g", "c", "m", "y"]
        random.shuffle(colors)

        for i, (genus, values) in enumerate(genus_values.items()):
            color = colors[i % len(colors)]
            ax.plot(categories, values, marker='o', linestyle='-', color=color, label=genus)
            ax.fill_between(categories, values, color=color, alpha=0.2)

        ax.set_yscale("log")
        ax.set_xlabel("Category")
        ax.set_ylabel("Number of Species (Log Scale)")
        ax.set_title("Parallel Coordinates Plot for Selected Genera")
        ax.grid(True, which="both", linestyle="--", linewidth=0.5)
        ax.legend(title="Genus", loc="upper right")

        st.pyplot(fig)

    # --- Sankey Diagram Section ---
    st.markdown("---")
    st.title("Sankey Diagram: Genus → Strains → Metabolite Testing")
    st.write("Visualize how a genus splits into strains and how they are tested for a specific metabolite.")

    # Step 1: Load metabolite list from file
    data_folder = "data_files"
    metabolite_file = os.path.join(data_folder, "step3_overall_unique_mets.txt")

    if os.path.exists(metabolite_file):
        with open(metabolite_file, "r") as f:
            metabolite_list = f.read().splitlines()
    else:
        st.error("Metabolite file not found!")
        metabolite_list = []

    # Step 2: User selects a genus for Sankey (single selection)
    selected_genus_sankey = st.selectbox("Select a Genus for Sankey Diagram:", genus_list, key="sankey_genus")

    # Step 3: User selects a metabolite
    selected_metabolite = st.selectbox("Select a Metabolite:", metabolite_list, key="sankey_metabolite")

    # Step 4: Button to trigger Sankey plot generation
    if st.button("Generate Sankey Diagram"):
        # Step 5: Extract data for Yes Strains / No Strains
        strain_counts = {"Yes Strain": 0, "No Strain": 0}
        test_counts = {
            "Yes Strain → Positive": 0, "Yes Strain → Negative": 0,
            "No Strain → Positive": 0, "No Strain → Negative": 0
        }

        for strain_status in ["yesstrain", "nostrain"]:
            for test_status in ["positively", "negatively"]:
                file_name = f"step4_{test_status}_tested_by_genera_prod_{strain_status}.csv"
                if file_name in data_frames:
                    df = data_frames[file_name]
                    genus_row = df[df["genus"] == selected_genus_sankey]
                    if not genus_row.empty:
                        strain_counts["Yes Strain" if strain_status == "yesstrain" else "No Strain"] += genus_row.iloc[:, 1:].sum().sum()
                        if selected_metabolite in df.columns:
                            test_counts[f"{'Yes Strain' if strain_status == 'yesstrain' else 'No Strain'} → {'Positive' if test_status == 'positively' else 'Negative'}"] += genus_row[selected_metabolite].sum()

        # Step 6: Build Sankey Data
        labels = ["Genus", "Yes Strain", "No Strain", "Positive Test", "Negative Test"]
        sources = [0, 0, 1, 1, 2, 2]
        targets = [1, 2, 3, 4, 3, 4]

        def extract_numeric(value):
            if isinstance(value, pd.Series):
                return value.values[0] if len(value) == 1 else value.sum()
            return value

        values = [
            extract_numeric(strain_counts["Yes Strain"]),
            extract_numeric(strain_counts["No Strain"]),
            extract_numeric(test_counts["Yes Strain → Positive"]),
            extract_numeric(test_counts["Yes Strain → Negative"]),
            extract_numeric(test_counts["No Strain → Positive"]),
            extract_numeric(test_counts["No Strain → Negative"])
        ]

        link_colors = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#937860", "#da8bc3"]

        fig_sankey = go.Figure(go.Sankey(
            node=dict(
                pad=15, thickness=20, line=dict(color="black", width=0.5),
                label=labels, color="#636363"
            ),
            link=dict(
                source=sources, target=targets, value=values,
                color=link_colors,
                label=[f"{v} species" if v > 0 else "" for v in values]
            )
        ))

        fig_sankey.update_layout(
            title_text=f"Sankey Diagram for {selected_genus_sankey} → {selected_metabolite}",
            font_size=12
        )

        st.plotly_chart(fig_sankey)

        # Step 7: Provide download for the Sankey diagram as PNG
        buffer = BytesIO()
        fig_sankey.write_image(buffer, format="png")
        buffer.seek(0)

        st.download_button(
            label="Download Sankey Diagram as PNG",
            data=buffer,
            file_name=f"Sankey_{selected_genus_sankey}_{selected_metabolite}.png",
            mime="image/png"
        )
