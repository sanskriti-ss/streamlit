import streamlit as st
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import plotly.graph_objects as go
from io import BytesIO
import sys
from pathlib import Path

# Add parent directory to path to import utils
sys.path.append(str(Path(__file__).parent.parent))
from utils.tooltip_title import display_title_with_tooltip

def display(data_frames):
    display_title_with_tooltip(
        "Parallel Coordinates Plot for Genera Trends",
        sample_image_filename="trends_5genera_comparison.png",
        description_text="Compare up to 5 genera across Production, Utilization, Resistance, and Sensitivity."
    )

    ### Step 1: store selections so it stops re-running on every change
    if "test_type" not in st.session_state:
        st.session_state.test_type = "Positively Tested"
    if "strain_option" not in st.session_state:
        st.session_state.strain_option = "Isolate"
    if "selected_genera" not in st.session_state:
        st.session_state.selected_genera = []

    ### Step 2: allowing use to select
    st.session_state.test_type = st.radio(
        "Choose Test Type:", ["Negatively Tested", "Positively Tested"],
        key="test_type_radio"
    )
    st.session_state.strain_option = st.radio(
        "Include Strains?", ["Isolate", "Strain"],
        key="strain_option_radio"
    )

    ### Step 3: mapping selections to file names
    test_type_short = "negatively" if st.session_state.test_type == "Negatively Tested" else "positively"
    strain_short = "isolate" if st.session_state.strain_option == "Isolate" else "strain"
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }
    relevant_files = {
        category: f"step4_{test_type_short}_tested_by_genera_{short}_{strain_short}.csv"
        for category, short in category_mapping.items()
    }

    # Step 5: Extract the genus names (assuming they exist in all files)
    sample_df = data_frames[relevant_files["Production"]]
    genus_list = sample_df["genus"].tolist()

    # Step 6: Multi-selection for genera (up to 5)
    st.session_state.selected_genera = st.multiselect(
        "Search and Select Up to 5 Genera:",
        genus_list,
        default=st.session_state.selected_genera,
        key="selected_genera_multiselect"
    )

    ### Step 7: parallel plot
    if st.button("Generate Parallel Coordinates Diagram"):
        genus_values = _compute_parallel_coordinates(data_frames, relevant_files, st.session_state.selected_genera)
        categories = ["Production", "Utilization", "Resistance", "Sensitivity"]

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
    display_title_with_tooltip(
        "Comprehensive Sankey Diagram: Genus → Strains → Categories → Test Results",
        sample_image_filename="trends_sankey.png",
        description_text="Visualize the complete flow: how a genus splits by strain status, then by metabolite categories (Production, Utilization, Resistance, Sensitivity), and finally by test results (Positive/Negative)."
    )

    # Checkbox for filtering (default unchecked)
    filter_enabled = st.checkbox(
        "Only show non-zero genus-metabolite combinations",
        value=False,
        help="Check this to filter the dropdowns to only show genus/metabolite pairs that have actual test data."
    )

    # Step 1: load in metabolite list
    data_folder = "data_files"
    metabolite_file = os.path.join(data_folder, "step3_overall_unique_mets.txt")

    if os.path.exists(metabolite_file):
        with open(metabolite_file, "r") as f:
            metabolite_list = f.read().splitlines()
    else:
        st.error("Metabolite file not found!")
        metabolite_list = []

    # Step 2: If filtering is enabled, build index and filter options
    if filter_enabled:
        # Import the functions from species_analysis
        import sys
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
        from tabs.species_analysis import load_all_activity_dfs, build_genus_metabolite_index, build_metabolite_genus_index, get_filtered_options
        
        with st.spinner("Loading activity files and building index..."):
            activity_dfs = load_all_activity_dfs()
            
            # Check if any files loaded
            loaded_count = sum(1 for v in activity_dfs.values() if v is not None)
            if loaded_count == 0:
                st.error("No activity files found in species_data folder.")
                return
            
            # Create a hash for cache invalidation
            dfs_hash = str(hash(tuple(k for k, v in activity_dfs.items() if v is not None)))
            
            # Build indices
            genus_index = build_genus_metabolite_index(dfs_hash, activity_dfs)
            metabolite_index = build_metabolite_genus_index(dfs_hash, activity_dfs)
        
        if not genus_index:
            st.warning("No genus-metabolite data found.")
            return
        
        # Use session state to track selections
        if 'sankey_genus_filtered' not in st.session_state:
            st.session_state.sankey_genus_filtered = "-- choose --"
        if 'sankey_metabolite_filtered' not in st.session_state:
            st.session_state.sankey_metabolite_filtered = "-- choose --"
        
        # Get filtered options
        filtered_genera, filtered_metabolites = get_filtered_options(
            genus_index,
            metabolite_index,
            st.session_state.sankey_genus_filtered,
            st.session_state.sankey_metabolite_filtered
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            genus_options = ["-- choose --"] + filtered_genera
            # Get current value from session state
            current_genus_idx = 0
            if st.session_state.sankey_genus_filtered in genus_options:
                current_genus_idx = genus_options.index(st.session_state.sankey_genus_filtered)
            
            selected_genus_sankey = st.selectbox(
                "Select a Genus for Sankey Diagram:",
                genus_options,
                index=current_genus_idx,
                key="sankey_genus_filtered_select"
            )
            
            # Update session state if changed
            if selected_genus_sankey != st.session_state.sankey_genus_filtered:
                st.session_state.sankey_genus_filtered = selected_genus_sankey
                # Reset metabolite if genus changed
                if selected_genus_sankey != "-- choose --":
                    valid_metabolites = genus_index.get(selected_genus_sankey, set())
                    if st.session_state.sankey_metabolite_filtered not in valid_metabolites:
                        st.session_state.sankey_metabolite_filtered = "-- choose --"
        
        with col2:
            if selected_genus_sankey and selected_genus_sankey != "-- choose --":
                metabolite_options = ["-- choose --"] + sorted(genus_index.get(selected_genus_sankey, []))
            else:
                metabolite_options = ["-- choose --"] + filtered_metabolites
            
            # Get current value from session state
            current_met_idx = 0
            if st.session_state.sankey_metabolite_filtered in metabolite_options:
                current_met_idx = metabolite_options.index(st.session_state.sankey_metabolite_filtered)
            
            selected_metabolite = st.selectbox(
                "Select a Metabolite:",
                metabolite_options,
                index=current_met_idx,
                key="sankey_metabolite_filtered_select"
            )
            
            # Update session state if changed
            if selected_metabolite != st.session_state.sankey_metabolite_filtered:
                st.session_state.sankey_metabolite_filtered = selected_metabolite
        
        # Show info
        if selected_genus_sankey != "-- choose --":
            metabolite_count = len(genus_index.get(selected_genus_sankey, []))
            st.info(f"{selected_genus_sankey} has test data for {metabolite_count} metabolites across all activity types.")
        
        # Show message if not both selected, but don't return - continue to button
        if selected_genus_sankey == "-- choose --" or selected_metabolite == "-- choose --":
            if selected_genus_sankey != "-- choose --" or selected_metabolite != "-- choose --":
                st.info("Please select both a genus and a metabolite to generate the diagram.")
    else:
        # Unfiltered mode - original behavior
        # Step 2: select one genus
        selected_genus_sankey = st.selectbox("Select a Genus for Sankey Diagram:", genus_list, key="sankey_genus")

        # Step 3: select one metabolite
        selected_metabolite = st.selectbox("Select a Metabolite:", metabolite_list, key="sankey_metabolite")

    # Step 4: button for sankey
    if st.button("Generate Sankey Diagram"):
        labels, sources, targets, values, link_colors, debug_info = _compute_sankey_data(
            data_frames, relevant_files, selected_genus_sankey, selected_metabolite
        )
        # does data exist?
        total_value = sum(values)
        st.write(f"**Total flow value: {total_value}**")
        
        if total_value == 0:
            st.warning(f"No data found for {selected_genus_sankey} and {selected_metabolite}. Please try a different combination.")
            fig_sankey = None
        else:
            # Filter out zero values to clean up the diagram
            filtered_sources = []
            filtered_targets = []
            filtered_values = []
            filtered_colors = []
            
            for i, value in enumerate(values):
                if value > 0:
                    filtered_sources.append(sources[i])
                    filtered_targets.append(targets[i])
                    filtered_values.append(value)
                    filtered_colors.append(link_colors[i])
            
            fig_sankey = go.Figure(go.Sankey(
                node=dict(
                    pad=15, 
                    thickness=20, 
                    line=dict(color="#000000", width=0.5),
                    label=labels, 
                    color="#4A4A4A"
                ),
                link=dict(
                    source=filtered_sources, target=filtered_targets, value=filtered_values,
                    color=filtered_colors,
                    label=[f"{v} specimens" if v > 0 else "" for v in filtered_values]
                ),
                textfont=dict(
                    family="Arial, sans-serif",
                    size=18,
                    color="#FFFFFF"
                )
            ))
            fig_sankey.update_layout(
                title_text=f"Comprehensive Sankey: {selected_genus_sankey} → {selected_metabolite}",
                title_font=dict(
                    family="Arial, sans-serif",
                    size=20,
                    color="#000000"
                ),
                font=dict(
                    family="Arial, sans-serif",
                    size=18,
                    color="#000000"
                ),
                height=600,
                paper_bgcolor="#FAFAFA",
                plot_bgcolor="#FAFAFA"
            )

            st.plotly_chart(fig_sankey, use_container_width=True)

        # Download functionality - only show if diagram was created
        if 'fig_sankey' in locals() and fig_sankey is not None:
            try:
                buffer = BytesIO()
                fig_sankey.write_image(buffer, format="png")
                buffer.seek(0)

                st.download_button(
                    label="Download Sankey Diagram as PNG",
                    data=buffer,
                    file_name=f"Comprehensive_Sankey_{selected_genus_sankey}_{selected_metabolite}.png",
                    mime="image/png"
                )
            except Exception as e:
                st.warning(f"Download feature temporarily unavailable: {str(e)}")
                st.info("You can still right-click on the diagram and save it manually.")

@st.cache_data
def _compute_parallel_coordinates(data_frames, relevant_files, selected_genera):
    """
    Compute the values for the parallel coordinates plot.
    """
    categories = ["Production", "Utilization", "Resistance", "Sensitivity"]
    genus_values = {}

    for genus in selected_genera:
        values = []
        for category in categories:
            file_name = relevant_files[category]
            df = data_frames[file_name]
            row = df[df["genus"] == genus]
            if not row.empty:
                # Sum all chemical presence counts (assuming columns 2 onward)
                values.append(row.iloc[0, 1:].sum())
            else:
                values.append(0)
        genus_values[genus] = values

    return genus_values

@st.cache_data
def _compute_sankey_data(data_frames, relevant_files, selected_genus_sankey, selected_metabolite):
    """
    Compute the values for the Sankey diagram with complete flow:
    Genus → Strain Status → Category Type → Test Result
    """
    
    # Initialize data structure for the comprehensive flow
    category_mapping = {
        "Production": "prod",
        "Utilization": "util",
        "Resistance": "res",
        "Sensitivity": "sen"
    }
    
    # Track counts for each flow path
    strain_counts = {"Strain": 0, "Isolate": 0}
    category_counts = {}
    test_counts = {}
    debug_info = []
    
    # Initialize category and test count dictionaries
    for strain in ["Strain", "Isolate"]:
        for category in category_mapping.keys():
            category_key = f"{strain} → {category}"
            category_counts[category_key] = 0
            
            for test_result in ["Positive", "Negative"]:
                test_key = f"{strain} → {category} → {test_result}"
                test_counts[test_key] = 0

    # Process each file to collect data
    for strain_status in ["strain", "isolate"]:
        strain_label = "Strain" if strain_status == "strain" else "Isolate"
        
        for category_name, category_short in category_mapping.items():
            for test_status in ["positively", "negatively"]:
                test_label = "Positive" if test_status == "positively" else "Negative"
                
                file_name = f"step4_{test_status}_tested_by_genera_{category_short}_{strain_status}.csv"
                debug_info.append(f"Looking for file: {file_name}")
                
                if file_name in data_frames:
                    df = data_frames[file_name]
                    debug_info.append(f"  - File found with {len(df)} rows and {len(df.columns)} columns")
                    
                    # Find genus row
                    genus_rows = df[df["genus"] == selected_genus_sankey]
                    debug_info.append(f"  - Found {len(genus_rows)} rows for genus '{selected_genus_sankey}'")
                    
                    if not genus_rows.empty:
                        # Check if metabolite exists (case-insensitive search)
                        metabolite_cols = [col for col in df.columns if col.lower() == selected_metabolite.lower()]
                        
                        if metabolite_cols:
                            metabolite_col = metabolite_cols[0]  # Use the first match
                            metabolite_value = genus_rows[metabolite_col].iloc[0]
                            debug_info.append(f"  - Found metabolite '{metabolite_col}' with value: {metabolite_value}")
                            
                            if pd.notna(metabolite_value) and metabolite_value > 0:
                                value = int(metabolite_value)
                                debug_info.append(f"  - Using value: {value}")
                                
                                # Add to strain totals
                                strain_counts[strain_label] += value
                                
                                # Add to category totals
                                category_key = f"{strain_label} → {category_name}"
                                category_counts[category_key] += value
                                
                                # Add to test result totals
                                test_key = f"{strain_label} → {category_name} → {test_label}"
                                test_counts[test_key] += value
                                
                                debug_info.append(f"  - Added {value} to: {test_key}")
                        else:
                            # Show available metabolites for debugging
                            available_mets = [col for col in df.columns if col not in ['genus', 'species_count']][:10]
                            debug_info.append(f"  - Metabolite '{selected_metabolite}' not found. Available: {available_mets}...")
                    else:
                        available_genera = df['genus'].unique()[:10]
                        debug_info.append(f"  - Genus not found. Available genera: {available_genera.tolist()}...")
                else:
                    debug_info.append(f"  - File not found in data_frames")

    debug_info.append(f"\nFinal strain counts: {strain_counts}")
    debug_info.append(f"Final category counts: {category_counts}")
    debug_info.append(f"Final test counts: {test_counts}")

    # Build comprehensive Sankey diagram with 4 levels
    labels = [
        # Level 0: Genus
        selected_genus_sankey,  # 0
        
        # Level 1: Strain Status
        "Strain",          # 1
        "Isolate",           # 2
        
        # Level 2: Categories (Strain)
        "YS Production",       # 3
        "YS Utilization",      # 4  
        "YS Resistance",       # 5
        "YS Sensitivity",      # 6
        
        # Level 2: Categories (Isoalte)  
        "NS Production",       # 7
        "NS Utilization",      # 8
        "NS Resistance",       # 9
        "NS Sensitivity",      # 10
        
        # Level 3: Test Results
        "Positive Test",       # 11
        "Negative Test"        # 12
    ]
    
    sources = []
    targets = []
    values = []
    
    # Level 0 → Level 1: Genus to Strain Status
    sources.extend([0, 0])
    targets.extend([1, 2])
    values.extend([
        strain_counts["Strain"],
        strain_counts["Isolate"]
    ])
    
    # Level 1 → Level 2: Strain Status to Categories
    # Yes Strain to categories
    sources.extend([1, 1, 1, 1])
    targets.extend([3, 4, 5, 6])
    values.extend([
        category_counts["Strain → Production"],
        category_counts["Strain → Utilization"],
        category_counts["Strain → Resistance"],
        category_counts["Strain → Sensitivity"]
    ])

    # Isolate to categories
    sources.extend([2, 2, 2, 2])
    targets.extend([7, 8, 9, 10])
    values.extend([
        category_counts["Isolate → Production"],
        category_counts["Isolate → Utilization"],
        category_counts["Isolate → Resistance"],
        category_counts["Isolate → Sensitivity"]
    ])
    
    # Level 2 → Level 3: Categories to Test Results
    category_indices = [3, 4, 5, 6, 7, 8, 9, 10]
    category_names = ["Production", "Utilization", "Resistance", "Sensitivity"] * 2
    strain_names = ["Strain"] * 4 + ["Isolate"] * 4

    for i, (category_idx, category_name, strain_name) in enumerate(zip(category_indices, category_names, strain_names)):
        # To Positive Test
        sources.append(category_idx)
        targets.append(11)
        values.append(test_counts[f"{strain_name} → {category_name} → Positive"])
        
        # To Negative Test
        sources.append(category_idx)
        targets.append(12)
        values.append(test_counts[f"{strain_name} → {category_name} → Negative"])
    
    # Color scheme
    link_colors = (
        ["#4c72b0", "#55a868"] +  # Genus to strain (blue, green)
        ["#87ceeb"] * 4 +         # Strain to categories (light blue)
        ["#90ee90"] * 4 +         # Isolate to categories (light green)
        ["#2ca02c"] * 8 +         # Categories to positive (bright green)
        ["#d62728"] * 8           # Categories to negative (red)
    )
    
    return labels, sources, targets, values, link_colors, debug_info
