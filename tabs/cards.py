import streamlit as st
import pandas as pd

def display(data_frames):
    st.header("Cards")
    st.write("Click on a card below to view details for each genus. (Details will be added soon.)")
    
    # Use one CSV file (e.g., one with "nostrain") to get the list of genera.
    file_key = next((fname for fname in data_frames if "nostrain" in fname), None)
    if not file_key:
        st.error("No data available for cards.")
        return
    
    df = data_frames[file_key]
    # Get a sorted list of unique genera.
    genera = sorted(df["genus"].unique())
    
    # Arrange the cards in a grid with 3 cards per row.
    num_cols = 3
    rows = [genera[i:i+num_cols] for i in range(0, len(genera), num_cols)]
    
    # Default file naming parameters: positively tested, no strain.
    # For species count, we'll use the Production file.
    prod_key = "step4_positively_tested_by_genera_prod_nostrain.csv"
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }
    
    for row in rows:
        cols = st.columns(num_cols)
        for idx, genus in enumerate(row):
            with cols[idx]:
                with st.expander(genus, expanded=False):
                    # Get the number of species from the Production file.
                    if prod_key in data_frames:
                        df_prod = data_frames[prod_key]
                        row_prod = df_prod[df_prod["genus"] == genus]
                        if not row_prod.empty:
                            species_count = row_prod["species_count"].iloc[0]
                        else:
                            species_count = "N/A"
                    else:
                        species_count = "N/A"
                    st.write(f"**Number of species in {genus}:** {species_count}")
                    
                    # For each category, count unique metabolites (without double counting!!!).
                    for cat, short in category_mapping.items():
                        file_key_cat = f"step4_positively_tested_by_genera_{short}_nostrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            if not row_cat.empty:
                                # Exclude the non-metabolite columns.
                                metabolite_cols = row_cat.columns.difference(["genus", "species_count"])
                                # Count each metabolite as 1 if its value > 0.
                                unique_count = (row_cat[metabolite_cols] > 0).sum(axis=1).iloc[0]
                            else:
                                unique_count = 0
                        else:
                            unique_count = "N/A"
                        st.write(f"**Unique {cat} metabolites:** {unique_count}")
