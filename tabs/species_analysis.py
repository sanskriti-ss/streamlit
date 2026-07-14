##########################################################
### 1) Imports 
##########################################################

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import zipfile
import io
import numpy as np
import os
from typing import Dict, List, Optional, Set
from utils.species_utils import (
    load_species_data_from_path, 
    calculate_species_stats, 
    process_multiple_files,
    create_comparison_table,
    get_top_species_across_files
)

##########################################################
### Helper Functions
##########################################################

##########################################################
### Cached Functions for Filtered Selectboxes
##########################################################

@st.cache_data(show_spinner="Loading activity files...")
def load_all_activity_dfs() -> Dict[str, pd.DataFrame]:
    """
    Load all 4 activity dataframes from species_data folder.
    Returns dict: {'production': df, 'utilization': df, 'resistance': df, 'sensitivity': df}
    """
    species_folder = "species_data"
    activity_file_map = {
        'production': os.path.join(species_folder, "step3_met_prod_exploded.csv.zip"),
        'utilization': os.path.join(species_folder, "step3_met_util_exploded.csv.zip"),
        'resistance': os.path.join(species_folder, "step3_met_res_exploded.csv.zip"),
        'sensitivity': os.path.join(species_folder, "step3_met_sen_exploded.csv.zip")
    }
    
    dfs = {}
    for activity, path in activity_file_map.items():
        if os.path.exists(path):
            try:
                dfs[activity] = load_species_data_from_path(path)
            except Exception as e:
                st.warning(f"Could not load {activity} file: {e}")
                dfs[activity] = None
        else:
            dfs[activity] = None
    return dfs


@st.cache_data(show_spinner="Building genus-metabolite index...")
def build_genus_metabolite_index(_dfs_hash: str, dfs: Dict[str, pd.DataFrame]) -> Dict[str, Set[str]]:
    """
    Build genus -> set(metabolites) where there is any non-zero value
    across any activity dataframe. Vectorized for speed.
    
    _dfs_hash is used for cache invalidation (pass a hash of file paths or timestamps).
    """
    genus_index = {}
    
    for activity, df in dfs.items():
        if df is None or df.empty:
            continue
        
        # Detect metadata cols
        metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain', 'is_strain', 'species_with_id']
        metabolite_cols = [c for c in df.columns if c not in metadata_cols]
        
        if not metabolite_cols or 'genus' not in df.columns:
            continue
        
        # Vectorized: tested = non-NaN and not-equal-to-0 (counts 1 and -1)
        metabolite_df = df[metabolite_cols]
        tested = metabolite_df.notna() & (metabolite_df != 0)
        tested_with_genus = tested.copy()
        tested_with_genus['genus'] = df['genus'].fillna('Unknown')
        
        # Group by genus and check if any strain has nonzero for each metabolite
        grp = tested_with_genus.groupby('genus').any()
        
        for genus in grp.index:
            # Get metabolite columns where at least one strain has nonzero
            cols_true = set(grp.columns[grp.loc[genus].values])
            genus_index.setdefault(genus, set()).update(cols_true)
    
    return genus_index


@st.cache_data(show_spinner="Building metabolite-genus index...")
def build_metabolite_genus_index(_dfs_hash: str, dfs: Dict[str, pd.DataFrame]) -> Dict[str, Set[str]]:
    """
    Build metabolite -> set(genera) where there is any non-zero value.
    This is the reverse index for filtering genera after metabolite selection.
    """
    metabolite_index = {}
    
    for activity, df in dfs.items():
        if df is None or df.empty:
            continue
        
        metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain', 'is_strain', 'species_with_id']
        metabolite_cols = [c for c in df.columns if c not in metadata_cols]
        
        if not metabolite_cols or 'genus' not in df.columns:
            continue
        
        # For each metabolite, find genera with nonzero values
        for met in metabolite_cols:
            if met not in df.columns:
                continue
            # Get rows where metabolite is nonzero
            mask = df[met].notna() & (df[met] != 0)
            genera_with_data = set(df.loc[mask, 'genus'].dropna().unique())
            metabolite_index.setdefault(met, set()).update(genera_with_data)
    
    return metabolite_index


def get_filtered_options(
    genus_index: Dict[str, Set[str]],
    metabolite_index: Dict[str, Set[str]],
    selected_genus: Optional[str] = None,
    selected_metabolite: Optional[str] = None
) -> tuple:
    """
    Get filtered genus and metabolite lists based on current selections.
    
    Returns: (genus_list, metabolite_list)
    """
    all_genera = sorted(genus_index.keys())
    all_metabolites = sorted(set().union(*genus_index.values())) if genus_index else []
    
    # Filter metabolites based on selected genus
    if selected_genus and selected_genus != "-- choose --":
        filtered_metabolites = sorted(genus_index.get(selected_genus, []))
    else:
        filtered_metabolites = all_metabolites
    
    # Filter genera based on selected metabolite
    if selected_metabolite and selected_metabolite != "-- choose --":
        filtered_genera = sorted(metabolite_index.get(selected_metabolite, []))
    else:
        filtered_genera = all_genera
    
    return filtered_genera, filtered_metabolites


def get_all_options_from_df(df: pd.DataFrame) -> tuple:
    """
    Get all genus and metabolite options from a single dataframe (unfiltered mode).
    """
    if df is None or df.empty:
        return [], []
    
    metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain', 'is_strain', 'species_with_id']
    metabolite_cols = [c for c in df.columns if c not in metadata_cols]
    
    genera = sorted(df['genus'].dropna().unique()) if 'genus' in df.columns else []
    metabolites = sorted(metabolite_cols)
    
    return genera, metabolites


def load_species_data(file_path):
    return load_species_data_from_path(file_path)

def calculate_metabolite_utilization(df, file_path=""):
    # Detect activity type from file path
    activity_type = "resistance"  # default
    file_path_lower = file_path.lower()
    
    if 'res' in file_path_lower or 'resistance' in file_path_lower:
        activity_type = "resistance"
    elif 'prod' in file_path_lower or 'production' in file_path_lower:
        activity_type = "production"
    elif 'sen' in file_path_lower or 'sensitivity' in file_path_lower:
        activity_type = "sensitivity"
    elif 'util' in file_path_lower or 'utilization' in file_path_lower:
        activity_type = "utilization"
    
    return calculate_species_stats(df, activity_type)

def create_top_species_chart(stats_df, top_n=20, analysis_title=""):
    """
    Create a bar chart showing top species by metabolite utilization.
    
    Args:
        stats_df (pd.DataFrame): Species utilization statistics
        top_n (int): Number of top species to show
        analysis_title (str): Title to determine chart labels
        
    Returns:
        plotly.graph_objects.Figure: Bar chart
    """
    if stats_df is None or stats_df.empty:
        return None
    
    # Determine which activity column to use based on analysis title
    activity_col = 'metabolites_utilized'  # default
    rate_col = 'utilization_rate'  # default
    
    title_lower = analysis_title.lower()
    
    if ('res' in title_lower or 'resistance' in title_lower) and 'metabolites_resistant' in stats_df.columns:
        activity_col = 'metabolites_resistant'
        rate_col = 'resistance_rate'
        chart_title = f"Top {top_n} Species by Metabolite Resistance"
        x_title = "Number of Metabolites Resistant To"
        rate_title = "Resistance Rate (%)"
    elif ('prod' in title_lower or 'production' in title_lower) and 'metabolites_produced' in stats_df.columns:
        activity_col = 'metabolites_produced'
        rate_col = 'production_rate'
        chart_title = f"Top {top_n} Species by Metabolite Production"
        x_title = "Number of Metabolites Produced"
        rate_title = "Production Rate (%)"
    elif ('sen' in title_lower or 'sensitivity' in title_lower) and 'metabolites_sensitive' in stats_df.columns:
        activity_col = 'metabolites_sensitive'
        rate_col = 'sensitivity_rate'
        chart_title = f"Top {top_n} Species by Metabolite Sensitivity"
        x_title = "Number of Metabolites Sensitive To"
        rate_title = "Sensitivity Rate (%)"
    elif ('util' in title_lower or 'utilization' in title_lower) and 'metabolites_utilized' in stats_df.columns:
        activity_col = 'metabolites_utilized'
        rate_col = 'utilization_rate'
        chart_title = f"Top {top_n} Species by Metabolite Utilization"
        x_title = "Number of Metabolites Utilized"
        rate_title = "Utilization Rate (%)"
    else:
        # Fallback: find the column with the highest total activity
        potential_cols = {
            'metabolites_resistant': ('resistance_rate', 'Resistance', 'Resistant To'),
            'metabolites_produced': ('production_rate', 'Production', 'Produced'),
            'metabolites_sensitive': ('sensitivity_rate', 'Sensitivity', 'Sensitive To'),
            'metabolites_utilized': ('utilization_rate', 'Utilization', 'Utilized')
        }
        
        max_total = 0
        activity_type = 'Utilization'
        x_action = 'Utilized'
        
        for col, (rate, act_type, x_act) in potential_cols.items():
            if col in stats_df.columns:
                total = stats_df[col].sum()
                if total > max_total:
                    max_total = total
                    activity_col = col
                    rate_col = rate
                    activity_type = act_type
                    x_action = x_act
        
        chart_title = f"Top {top_n} Species by Metabolite {activity_type}"
        x_title = f"Number of Metabolites {x_action}"
        rate_title = f"{activity_type} Rate (%)"
    
    # Sort by the appropriate activity column and get top N
    top_species = stats_df.nlargest(top_n, activity_col)
    
    # Create bar chart
    fig = go.Figure()
    
    # Use species_with_id if available, otherwise fall back to species
    if 'species_with_id' in top_species.columns:
        y_labels = [f"{row['species_with_id']}<br>({row['genus']})" for _, row in top_species.iterrows()]
    else:
        y_labels = [f"{row['species']}<br>({row['genus']})" for _, row in top_species.iterrows()]
    
    fig.add_trace(go.Bar(
        x=top_species[activity_col],
        y=y_labels,
        orientation='h',
        text=top_species[activity_col],
        textposition='auto',
        marker=dict(
            color=top_species[rate_col],  # Use the correct rate column
            colorscale='Viridis',
            colorbar=dict(title=rate_title)
        ),
        hovertemplate="<b>%{y}</b><br>" +
                     f"{x_title}: %{{x}}<br>" +
                     "Order: %{customdata[0]}<br>" +
                     "Strain Type: %{customdata[1]}<br>" +
                     f"{rate_title}: %{{customdata[2]:.1f}}%<br>" +
                     "<extra></extra>",
        customdata=[[row['order'], 
                    row['type_strain'] if 'type_strain' in row else 'Unknown',
                    row[rate_col]] 
                   for _, row in top_species.iterrows()]
    ))
    
    fig.update_layout(
        title=chart_title,
        xaxis_title=x_title,
        yaxis_title="Species (Genus)",
        height=max(400, top_n * 25),
        font=dict(size=12),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    return fig

def create_utilization_heatmap(df, top_species_list, top_metabolites=50):
    """
    Create a heatmap showing metabolite utilization for top species.
    
    Args:
        df (pd.DataFrame): Raw species data
        top_species_list (list): List of top species names
        top_metabolites (int): Number of most variable metabolites to show
        
    Returns:
        plotly.graph_objects.Figure: Heatmap
    """
    if df is None or df.empty:
        return None
    
    # Filter to top species - use species column for filtering but species_with_id for display
    df_top = df[df['species'].isin(top_species_list)].copy()
    
    # Add species_with_id column if it doesn't exist
    if 'species_with_id' not in df_top.columns and 'BacID' in df_top.columns:
        df_top['species_with_id'] = df_top['species'] + ' (ID: ' + df_top['BacID'].astype(str) + ')'
    
    # Identify metabolite columns
    metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain', 'species_with_id']
    if 'is_strain' in df.columns:
        metadata_cols.append('is_strain')
    metabolite_cols = [col for col in df.columns if col not in metadata_cols]
    
    # Convert -1 to 0
    for col in metabolite_cols:
        df_top[col] = df_top[col].replace(-1, 0).clip(0, 1)
    
    # Find most variable metabolites
    metabolite_variance = df_top[metabolite_cols].var().sort_values(ascending=False)
    top_metabolite_cols = metabolite_variance.head(top_metabolites).index.tolist()
    
    # Create heatmap data - use species_with_id if available for y-axis labels
    if 'species_with_id' in df_top.columns:
        heatmap_data = df_top.set_index('species_with_id')[top_metabolite_cols]
    else:
        heatmap_data = df_top.set_index('species')[top_metabolite_cols]
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='RdYlBu_r',
        zmid=0.5,
        colorbar=dict(title="Utilization<br>(0=No, 1=Yes)")
    ))
    
    fig.update_layout(
        title=f"Metabolite Utilization Heatmap - Top {len(top_species_list)} Species",
        xaxis_title="Metabolites",
        yaxis_title="Species",
        height=max(400, len(top_species_list) * 30),
        xaxis=dict(tickangle=45),
        font=dict(size=10)
    )
    
    return fig

def display_top_species_table(stats_df, top_n=5, analysis_title=""):
    """
    Display a formatted table of top species.
    
    Args:
        stats_df (pd.DataFrame): Species utilization statistics
        top_n (int): Number of top species to show
        analysis_title (str): Title to determine the type of analysis
    """
    if stats_df is None or stats_df.empty:
        st.error("No data available for analysis")
        return []
    
    # Determine which activity column to use based on analysis title
    activity_col = 'metabolites_utilized'  # default
    rate_col = 'utilization_rate'  # default
    metabolite_type = "Utilization"
    action_verb = "Utilized"
    
    title_lower = analysis_title.lower()
    
    if ('res' in title_lower or 'resistance' in title_lower) and 'metabolites_resistant' in stats_df.columns:
        activity_col = 'metabolites_resistant'
        rate_col = 'resistance_rate'
        metabolite_type = "Resistance"
        action_verb = "Resistant to"
    elif ('prod' in title_lower or 'production' in title_lower) and 'metabolites_produced' in stats_df.columns:
        activity_col = 'metabolites_produced'
        rate_col = 'production_rate'
        metabolite_type = "Production"
        action_verb = "Produced"
    elif ('sen' in title_lower or 'sensitivity' in title_lower) and 'metabolites_sensitive' in stats_df.columns:
        activity_col = 'metabolites_sensitive'
        rate_col = 'sensitivity_rate'
        metabolite_type = "Sensitivity"
        action_verb = "Sensitive to"
    elif ('util' in title_lower or 'utilization' in title_lower) and 'metabolites_utilized' in stats_df.columns:
        activity_col = 'metabolites_utilized'
        rate_col = 'utilization_rate'
        metabolite_type = "Utilization"
        action_verb = "Utilized"
    else:
        # Fallback: find the column with the highest total activity
        potential_cols = {
            'metabolites_resistant': ('resistance_rate', 'Resistance', 'Resistant to'),
            'metabolites_produced': ('production_rate', 'Production', 'Produced'),
            'metabolites_sensitive': ('sensitivity_rate', 'Sensitivity', 'Sensitive to'),
            'metabolites_utilized': ('utilization_rate', 'Utilization', 'Utilized')
        }
        
        max_total = 0
        for col, (rate, met_type, verb) in potential_cols.items():
            if col in stats_df.columns:
                total = stats_df[col].sum()
                if total > max_total:
                    max_total = total
                    activity_col = col
                    rate_col = rate
                    metabolite_type = met_type
                    action_verb = verb
    
    # Get top species by the appropriate activity column
    top_species = stats_df.nlargest(top_n, activity_col)
    
    # Format the table - include strain/isolate information if available
    display_columns = ['species_with_id', 'genus', 'order', activity_col, rate_col]
    if 'type_strain' in top_species.columns:
        display_columns.insert(-2, 'type_strain')  # Insert before activity column
    
    display_df = top_species[display_columns].copy()
    display_df[rate_col] = display_df[rate_col].round(2)
    
    # Set column names
    column_names = ['Species (BacID)', 'Genus', 'Order', f'Metabolites {action_verb}', f'{metabolite_type} Rate (%)']
    if 'type_strain' in display_columns:
        column_names.insert(-2, 'Strain Type')
    
    display_df.columns = column_names
    display_df.index = range(1, len(display_df) + 1)
    
    st.subheader(f"Top {top_n} Species by Metabolite {metabolite_type}")
    st.dataframe(display_df, use_container_width=True)
    
    return top_species['species'].tolist()

def display(data_frames=None):
    """
    Main display function for the Species Analysis tab.
    
    Args:
        data_frames (dict): Dictionary of loaded data frames (for compatibility)
    """
    st.title("Species Metabolite Resistance Analysis")
    st.caption(
        "Bacterial BacDive deep-dive. For bacteria ↔ fungi metabolic partners, "
        "use the Symbiosis Network tab."
    )
    if st.button(
        "Open Symbiosis Network (bacteria + fungi)",
        key="species_analysis_to_symbiosis",
        help="Switch to the cross-kingdom synergy explorer",
    ):
        st.session_state.selected_tab = "Symbiosis Network"
        st.query_params["tab"] = "Symbiosis Network"
        st.rerun()
    st.markdown("---")
    
    # Analysis mode selection
    st.subheader("Analysis Mode")
    analysis_mode = st.radio(
        "Choose analysis type:",
        ["Single File Analysis", "Multi-File Comparison (Prod/Res/Sen/Util)"],
        help="Single file analyzes one dataset, Multi-file compares across production, resistance, sensitivity, and utilization data"
    )
    
    if analysis_mode == "Single File Analysis":
        display_single_file_analysis()
    else:
        display_multi_file_analysis()

def display_single_file_analysis():
    """Display single file analysis interface."""
    st.subheader("Single File Data Upload")
    
    # Check for available files
    import os
    species_folder = "species_data"
    available_files = []
    
    if os.path.exists(species_folder):
        expected_files = [
            "step3_met_prod_exploded.csv.zip",
            "step3_met_util_exploded.csv.zip", 
            "step3_met_res_exploded.csv.zip",
            "step3_met_sen_exploded.csv.zip"
        ]
        
        for filename in expected_files:
            file_path = os.path.join(species_folder, filename)
            if os.path.exists(file_path):
                available_files.append((file_path, filename))
    
    if available_files:
        use_default = st.checkbox("Use file from species_data folder", value=True)
        
        if use_default:
            selected_file = st.selectbox(
                "Select a data file:",
                available_files,
                format_func=lambda x: x[1],
                index=0
            )
            file_path = selected_file[0]
            file_name = selected_file[1]
            st.info(f"Using: {file_name}")
        else:
            uploaded_file = st.file_uploader(
                "Upload a compressed CSV file (.zip)", 
                type=['zip'],
                help="Upload a ZIP file containing CSV data with species and metabolite information"
            )
            
            if uploaded_file is not None:
                # Save uploaded file temporarily
                file_path = f"temp_{uploaded_file.name}"
                file_name = uploaded_file.name
                with open(file_path, "wb") as f:
                    f.write(uploaded_file.getbuffer())
            else:
                st.warning("Please upload a file or select a default file")
                return
    else:
        st.warning("No default files found. Please upload a file.")
        uploaded_file = st.file_uploader(
            "Upload a compressed CSV file (.zip)", 
            type=['zip'],
            help="Upload a ZIP file containing CSV data with species and metabolite information"
        )
        
        if uploaded_file is not None:
            # Save uploaded file temporarily
            file_path = f"temp_{uploaded_file.name}"
            file_name = uploaded_file.name
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
        else:
            st.warning("Please upload a file")
            return
    
    # Load and process data
    with st.spinner("Loading and processing data..."):
        df = load_species_data(file_path)
        
        if df is None:
            return
        
        # Calculate utilization statistics
        stats_df = calculate_metabolite_utilization(df, file_path)
        
        stats_df = stats_df[
            (~stats_df['species'].str.contains("Unknown", case=False, na=False)) &
            (~stats_df['genus'].str.contains("Unknown", case=False, na=False)) &
            (~stats_df['order'].str.contains("Unknown", case=False, na=False))
        ]
        
        if stats_df is None:
            st.error("Failed to process utilization statistics")
            return
    
    # Determine analysis title from file name
    analysis_title = "Single File"
    if 'res' in file_name.lower():
        analysis_title = "Resistance"
    elif 'prod' in file_name.lower():
        analysis_title = "Production"
    elif 'sen' in file_name.lower():
        analysis_title = "Sensitivity"
    elif 'util' in file_name.lower():
        analysis_title = "Utilization"
    
    display_analysis_results(df, stats_df, analysis_title)
    
def display_multi_file_analysis():
    """Display multi-file analysis interface."""
    st.subheader("Multi-File Data Upload")
    st.info("Upload production, resistance, sensitivity, and utilization data files for comparison")
    
    # File upload options
    use_default = st.checkbox("Use files from species_data folder", value=True)
    
    if use_default:
        # Look for default files
        import os
        species_folder = "species_data"
        
        # Define expected files
        expected_files = {
            "step3_met_prod_exploded.csv.zip": "Production",
            "step3_met_util_exploded.csv.zip": "Utilization", 
            "step3_met_res_exploded.csv.zip": "Resistance",
            "step3_met_sen_exploded.csv.zip": "Sensitivity"
        }
        
        available_files = []
        file_types = []
        
        if os.path.exists(species_folder):
            for filename, file_type in expected_files.items():
                file_path = os.path.join(species_folder, filename)
                if os.path.exists(file_path):
                    available_files.append(file_path)
                    file_types.append(file_type)
        
        if available_files:
            st.success(f"Found {len(available_files)} data files:")
            for file, file_type in zip(available_files, file_types):
                st.write(f"- **{file_type}**: {os.path.basename(file)}")
            
            # Option to select specific files or use all
            analysis_option = st.radio(
                "Analysis options:",
                ["Analyze all files", "Select specific files"],
                index=0
            )
            
            if analysis_option == "Select specific files":
                selected_indices = st.multiselect(
                    "Select files to analyze:",
                    range(len(available_files)),
                    format_func=lambda x: f"{file_types[x]}: {os.path.basename(available_files[x])}",
                    default=list(range(len(available_files)))
                )
                
                if selected_indices:
                    selected_files = [available_files[i] for i in selected_indices]
                    selected_types = [file_types[i] for i in selected_indices]
                    process_multiple_file_analysis(selected_files, selected_types)
                else:
                    st.warning("Please select at least one file to analyze")
            else:
                # Analyze all available files
                process_multiple_file_analysis(available_files, file_types)
        else:
            st.warning("No files found in species_data folder")
    
    else:
        # Manual file upload
        col1, col2, col3 = st.columns(3)
        
        uploaded_files = []
        file_types = []
        
        with col1:
            st.write("**Production Data**")
            prod_file = st.file_uploader("Upload production file", type=['zip', 'csv'], key="prod")
            if prod_file:
                uploaded_files.append(prod_file)
                file_types.append('Production')
        
        with col2:
            st.write("**Resistance Data**")  
            res_file = st.file_uploader("Upload resistance file", type=['zip', 'csv'], key="res")
            if res_file:
                uploaded_files.append(res_file)
                file_types.append('Resistance')
        
        with col3:
            st.write("**Sensitivity Data**")
            sen_file = st.file_uploader("Upload sensitivity file", type=['zip', 'csv'], key="sen")
            if sen_file:
                uploaded_files.append(sen_file)
                file_types.append('Sensitivity')
        
        if uploaded_files:
            # Save uploaded files temporarily
            file_paths = []
            for i, file in enumerate(uploaded_files):
                temp_path = f"temp_{file_types[i]}_{file.name}"
                with open(temp_path, "wb") as f:
                    f.write(file.getbuffer())
                file_paths.append(temp_path)
            
            process_multiple_file_analysis(file_paths, file_types)

def process_multiple_file_analysis(file_paths, file_types):
    """Process and display results for multiple files."""
    with st.spinner("Processing multiple files..."):
        results = process_multiple_files(file_paths, file_types)
        
        if results is None:
            st.error("Failed to process files")
            return
        
        display_multi_file_results(results)

def display_multi_file_results(results):
    """Display results from multi-file analysis."""
    st.markdown("---")
    st.subheader("Multi-File Analysis Results")
    
    combined_stats = results['combined_stats']
    individual_data = results['individual_data'] 
    top_species_overall = results['top_species_overall']
    
    # Summary statistics
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Files", len(individual_data))
    with col2:
        st.metric("Total Species", len(combined_stats['species'].unique()))
    with col3: 
        st.metric("File Types", len(combined_stats['file_type'].unique()))
    
    # Top species across all files
    st.subheader("Top Species Across All File Types")
    if not top_species_overall.empty:
        display_df = top_species_overall[['species', 'genus', 'order', 'metabolites_utilized', 'overall_utilization_rate', 'file_type']].copy()
        display_df.columns = ['Species', 'Genus', 'Order', 'Total Metabolites Utilized', 'Overall Utilization Rate (%)', 'File Types']
        display_df['Overall Utilization Rate (%)'] = display_df['Overall Utilization Rate (%)'].round(2)
        display_df.index = range(1, len(display_df) + 1)
        st.dataframe(display_df, use_container_width=True)
    
    # Comparison table
    st.subheader("Top 5 Species by File Type")
    comparison_table = create_comparison_table(results, top_n=5)
    if not comparison_table.empty:
        st.dataframe(comparison_table, use_container_width=True)
    
    # Individual file analysis
    st.subheader("Individual File Analysis")
    selected_file_type = st.selectbox("Select file type to analyze in detail:", list(individual_data.keys()))
    
    if selected_file_type:
        df = individual_data[selected_file_type]
        stats_df = combined_stats[combined_stats['file_type'] == selected_file_type].copy()
        display_analysis_results(df, stats_df, selected_file_type)


def display_analysis_results(df, stats_df, analysis_title):
    """Display analysis results for a dataset."""
    
    # Debug: Print available columns
    st.write("**Debug Info:**")
    st.write(f"Analysis Title: {analysis_title}")
    st.write(f"Stats DataFrame columns: {list(stats_df.columns)}")
    st.write(f"Stats DataFrame shape: {stats_df.shape}")
    if not stats_df.empty:
        st.write("First few rows of stats_df:")
        st.dataframe(stats_df.head())
        
        # Additional debug: Check the distribution of values
        st.write("**Value Distribution Debug:**")
        for col in ['metabolites_resistant', 'metabolites_produced', 'metabolites_sensitive', 'metabolites_utilized']:
            if col in stats_df.columns:
                st.write(f"{col}: min={stats_df[col].min()}, max={stats_df[col].max()}, mean={stats_df[col].mean():.2f}")
        
        # Check metabolites_tested column
        if 'total_metabolites' in stats_df.columns:
            st.write(f"**total_metabolites issue:** All values are {stats_df['total_metabolites'].iloc[0]} - this should vary per species!")
            st.write(f"Unique total_metabolites values: {stats_df['total_metabolites'].unique()}")
        
        # Debug raw data values
        st.write("**Raw Data Sample:**")
        if df is not None and not df.empty:
            # Show metabolite columns sample
            metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain', 'species_with_id']
            metabolite_cols = [col for col in df.columns if col not in metadata_cols][:10]  # First 10 metabolite columns
            sample_df = df[['species'] + metabolite_cols].head()
            st.dataframe(sample_df)
            
            # Show unique values in metabolite columns
            st.write("**Unique values in metabolite columns:**")
            for col in metabolite_cols[:5]:  # Check first 5 metabolite columns
                unique_vals = df[col].unique()
                st.write(f"{col}: {unique_vals}")
                
            # Calculate correct metabolites_tested for a few species manually
            st.write("**Manual calculation of metabolites_tested:**")
            for i in range(min(3, len(df))):
                species_row = df.iloc[i]
                species_name = species_row['species']
                # Count non-zero values in metabolite columns
                non_zero_count = sum(1 for col in metabolite_cols if species_row[col] != 0 and species_row[col] != -1)
                st.write(f"{species_name}: {non_zero_count} metabolites tested (non-zero)")
    
    # Determine which activity column to use based on analysis title and available data
    activity_col = 'metabolites_utilized'  # default
    rate_col = 'utilization_rate'  # default
    activity_name = 'Utilized'  # default
    
    # Check analysis title first to determine the expected activity type
    title_lower = analysis_title.lower()
    
    if ('res' in title_lower or 'resistance' in title_lower) and 'metabolites_resistant' in stats_df.columns:
        activity_col = 'metabolites_resistant'
        rate_col = 'resistance_rate'
        activity_name = 'Resistant'
    elif ('prod' in title_lower or 'production' in title_lower) and 'metabolites_produced' in stats_df.columns:
        activity_col = 'metabolites_produced'
        rate_col = 'production_rate'
        activity_name = 'Produced'
    elif ('sen' in title_lower or 'sensitivity' in title_lower) and 'metabolites_sensitive' in stats_df.columns:
        activity_col = 'metabolites_sensitive'
        rate_col = 'sensitivity_rate'
        activity_name = 'Sensitive'
    elif ('util' in title_lower or 'utilization' in title_lower) and 'metabolites_utilized' in stats_df.columns:
        activity_col = 'metabolites_utilized'
        rate_col = 'utilization_rate'
        activity_name = 'Utilized'
    else:
        # Fallback: find the column with the highest total activity
        potential_cols = {
            'metabolites_resistant': ('resistance_rate', 'Resistant'),
            'metabolites_produced': ('production_rate', 'Produced'),
            'metabolites_sensitive': ('sensitivity_rate', 'Sensitive'),
            'metabolites_utilized': ('utilization_rate', 'Utilized')
        }
        
        max_total = 0
        for col, (rate, name) in potential_cols.items():
            if col in stats_df.columns:
                total = stats_df[col].sum()
                if total > max_total:
                    max_total = total
                    activity_col = col
                    rate_col = rate
                    activity_name = name
    
    # Debug: Show what columns were selected
    st.write(f"**Selected columns:** Activity: {activity_col}, Rate: {rate_col}, Name: {activity_name}")
    
    # Check if we have any meaningful data
    if activity_col in stats_df.columns:
        max_activity = stats_df[activity_col].max()
        total_activity = stats_df[activity_col].sum()
        st.write(f"**Activity Data:** Max: {max_activity}, Total: {total_activity}")
        
        if max_activity == 0:
            st.error(f"⚠️ All species have 0 {activity_name.lower()} metabolites. This suggests a data processing issue.")
            st.info("This could be because:")
            st.info("1. The resistance data uses different value encoding (e.g., -1 for resistant instead of 1)")
            st.info("2. The data needs different preprocessing for resistance analysis")
            st.info("3. The file contains different data than expected")
            st.info("4. The 'total_metabolites' calculation is wrong - it should count metabolites tested per species, not total columns")
            return
    
    # Only proceed if we have meaningful data
    if activity_col in stats_df.columns and stats_df[activity_col].max() > 0:
        # Display summary statistics
        st.subheader(f"Data Summary - {analysis_title}")
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Species", len(stats_df))
        with col2:
            # Show range of metabolites tested instead of fixed total
            if 'total_metabolites' in stats_df.columns:
                min_tested = stats_df['total_metabolites'].min()
                max_tested = stats_df['total_metabolites'].max()
                avg_tested = stats_df['total_metabolites'].mean()
                st.metric("Metabolites Tested", f"{min_tested}-{max_tested} (avg: {avg_tested:.0f})")
            else:
                st.metric("Metabolites Tested", "Unknown")
        with col3:
            st.metric(f"Max {activity_name}", stats_df[activity_col].max() if not stats_df.empty and activity_col in stats_df.columns else 0)
        with col4:
            rate_mean = stats_df[rate_col].mean() if not stats_df.empty and rate_col in stats_df.columns else 0
            st.metric(f"Avg {activity_name} Rate", f"{rate_mean:.1f}%")
        
        # Interactive controls
        st.subheader("Analysis Controls")
        col1, col2 = st.columns(2)
        
        with col1:
            top_n_table = st.slider("Number of top species in table", 5, 20, 5, key=f"table_{analysis_title}")
            top_n_chart = st.slider("Number of species in chart", 10, 50, 20, key=f"chart_{analysis_title}")
        
        with col2:
            show_heatmap = st.checkbox("Show utilization heatmap", value=True, key=f"heatmap_{analysis_title}")
            top_metabolites = st.slider("Metabolites in heatmap", 20, 100, 50, key=f"metabolites_{analysis_title}") if show_heatmap else 50
        
        # Display results
        st.markdown("---")
        
        # Top species table
        top_species_list = display_top_species_table(stats_df, top_n_table, analysis_title)
        
        # Bar chart
        st.subheader("Top Species Visualization")
        chart_fig = create_top_species_chart(stats_df, top_n_chart, analysis_title)
        if chart_fig:
            st.plotly_chart(chart_fig, use_container_width=True)
        else:
            st.warning("No chart data available")
        
        # Heatmap
        if show_heatmap and top_species_list:
            st.subheader("Metabolite Activity Heatmap")
            heatmap_fig = create_utilization_heatmap(df, top_species_list[:10], top_metabolites)
            if heatmap_fig:
                st.plotly_chart(heatmap_fig, use_container_width=True)
            else:
                st.warning("No heatmap data available")
        
        # Download options
        st.markdown("---")
        st.subheader("Download Results")
        
        if st.button(f"Generate Download Data - {analysis_title}", key=f"download_{analysis_title}"):
            # Use the determined activity column for download
            if activity_col in stats_df.columns:
                download_df = stats_df.nlargest(50, activity_col)
                csv = download_df.to_csv(index=False)
                
                st.download_button(
                    label="Download Top 50 Species (CSV)",
                    data=csv,
                    file_name=f"top_species_metabolite_{activity_name.lower()}_{analysis_title.lower().replace(' ', '_')}.csv",
                    mime="text/csv",
                    key=f"download_btn_{analysis_title}"
                )
            else:
                st.error(f"Column {activity_col} not found in data")
        
        # Sankey Diagram Section
        # display_sankey_section(df, analysis_title)
    else:
        st.warning("No meaningful activity data found to display.")
        st.info("**Possible fixes needed in utils/species_utils.py:**")
        st.code("""
# In calculate_species_stats function, change this:
stats_df['total_metabolites'] = len(metabolite_cols)  # WRONG - same for all species

# To this:
for idx, row in df.iterrows():
    # Count non-zero values for each species
    metabolites_tested = sum(1 for col in metabolite_cols if row[col] != 0 and row[col] != -1)
    stats_df.loc[stats_df['species'] == row['species'], 'total_metabolites'] = metabolites_tested
        """)


##########################################################
### Sankey Diagram Section (with filtered selectboxes)
##########################################################

def display_sankey_section(df: pd.DataFrame = None, analysis_title: str = ""):
    """
    Display the Sankey diagram section with optional filtered selectboxes.
    """
    st.markdown("---")
    st.subheader("Comprehensive Sankey Diagram: Genus -> Strains -> Categories -> Test Results")
    st.write("Visualize the complete flow: how a genus splits by strain status, then by metabolite categories (Production, Utilization, Resistance, Sensitivity), and finally by test results (Positive/Negative).")
    
    # Checkbox for filtering (default unchecked)
    filter_enabled = st.checkbox(
        "Only show non-zero genus-metabolite combinations",
        value=False,
        help="Check this to filter the dropdowns to only show genus/metabolite pairs that have actual test data. This requires loading all activity files and may take a moment the first time."
    )
    
    if filter_enabled:
        # Load all activity files and build indices
        with st.spinner("Loading activity files and building index..."):
            activity_dfs = load_all_activity_dfs()
            
            # Check if any files loaded
            loaded_count = sum(1 for v in activity_dfs.values() if v is not None)
            if loaded_count == 0:
                st.error("No activity files found in species_data folder. Please ensure the files exist.")
                return
            
            # Create a hash for cache invalidation (based on which files are loaded)
            dfs_hash = str(hash(tuple(k for k, v in activity_dfs.items() if v is not None)))
            
            # Build indices
            genus_index = build_genus_metabolite_index(dfs_hash, activity_dfs)
            metabolite_index = build_metabolite_genus_index(dfs_hash, activity_dfs)
        
        if not genus_index:
            st.warning("No genus-metabolite data found in the loaded files.")
            return
        
        # Use session state to track selections and enable bidirectional filtering
        if 'sankey_genus' not in st.session_state:
            st.session_state.sankey_genus = "-- choose --"
        if 'sankey_metabolite' not in st.session_state:
            st.session_state.sankey_metabolite = "-- choose --"
        
        # Get filtered options based on current selections
        filtered_genera, filtered_metabolites = get_filtered_options(
            genus_index,
            metabolite_index,
            st.session_state.sankey_genus,
            st.session_state.sankey_metabolite
        )
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Genus selectbox
            genus_options = ["-- choose --"] + filtered_genera
            selected_genus = st.selectbox(
                "Select a Genus for Sankey Diagram:",
                genus_options,
                key="sankey_genus_select"
            )
            if selected_genus != st.session_state.sankey_genus:
                st.session_state.sankey_genus = selected_genus
                # Reset metabolite if it's no longer valid
                if selected_genus != "-- choose --":
                    valid_metabolites = genus_index.get(selected_genus, set())
                    if st.session_state.sankey_metabolite not in valid_metabolites and st.session_state.sankey_metabolite != "-- choose --":
                        st.session_state.sankey_metabolite = "-- choose --"
                st.rerun()
        
        with col2:
            # Metabolite selectbox (filtered based on genus if selected)
            if selected_genus and selected_genus != "-- choose --":
                metabolite_options = ["-- choose --"] + sorted(genus_index.get(selected_genus, []))
            else:
                metabolite_options = ["-- choose --"] + filtered_metabolites
            
            selected_metabolite = st.selectbox(
                "Select a Metabolite:",
                metabolite_options,
                key="sankey_metabolite_select"
            )
            if selected_metabolite != st.session_state.sankey_metabolite:
                st.session_state.sankey_metabolite = selected_metabolite
                st.rerun()
        
        # Show info about current selection
        if selected_genus != "-- choose --":
            metabolite_count = len(genus_index.get(selected_genus, []))
            st.info(f"**{selected_genus}** has test data for **{metabolite_count}** metabolites across all activity types.")
        
    else:
        # Unfiltered mode - use single df or load a default one
        if df is not None:
            genera, metabolites = get_all_options_from_df(df)
        else:
            # Try to load a default file
            default_path = "species_data/step3_met_res_exploded.csv.zip"
            if os.path.exists(default_path):
                df = load_species_data_from_path(default_path)
                genera, metabolites = get_all_options_from_df(df)
            else:
                st.warning("No data available. Please load a file first.")
                return
        
        col1, col2 = st.columns(2)
        
        with col1:
            selected_genus = st.selectbox(
                "Select a Genus for Sankey Diagram:",
                ["-- choose --"] + genera,
                key="sankey_genus_unfiltered"
            )
        
        with col2:
            selected_metabolite = st.selectbox(
                "Select a Metabolite:",
                ["-- choose --"] + metabolites,
                key="sankey_metabolite_unfiltered"
            )
        
        # Load all activity dfs for Sankey generation
        activity_dfs = load_all_activity_dfs()
    
    # Generate Sankey if both selections are made
    if selected_genus != "-- choose --" and selected_metabolite != "-- choose --":
        st.markdown("---")
        
        # Configuration
        top_k = st.slider("Number of top strains to show:", min_value=5, max_value=50, value=20, key="sankey_top_k")
        
        # Create Sankey
        with st.spinner("Generating Sankey diagram..."):
            fig = create_genus_sankey(
                activity_dfs if filter_enabled else load_all_activity_dfs(),
                selected_genus,
                selected_metabolite,
                top_k=top_k
            )
        
        if fig:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No data available for the selected genus/metabolite combination.")
    elif selected_genus != "-- choose --" or selected_metabolite != "-- choose --":
        st.info("Please select both a genus and a metabolite to generate the Sankey diagram.")


def create_genus_sankey(
    activity_dfs: Dict[str, pd.DataFrame],
    genus: str,
    metabolite: str,
    top_k: int = 20,
    positive_values_map: Dict[str, set] = None
) -> Optional[go.Figure]:
    """
    Build Sankey nodes/links for genus/metabolite across activity_dfs.
    
    positive_values_map: mapping activity -> set(values considered positive)
    e.g. {'resistance': {1}, ...}
    """
    if positive_values_map is None:
        positive_values_map = {
            'resistance': {1},      # 1 = resistant
            'sensitivity': {1},     # 1 = sensitive
            'production': {1},      # 1 = produces
            'utilization': {1}      # 1 = utilizes
        }
    
    # Collect rows for this genus across activities
    species_rows = {}  # species_with_id -> dict: category -> (tested, positive)
    categories = ['production', 'utilization', 'resistance', 'sensitivity']
    
    for cat in categories:
        df = activity_dfs.get(cat)
        if df is None or df.empty:
            continue
        if 'genus' not in df.columns or metabolite not in df.columns:
            continue
        
        df_gen = df[df['genus'] == genus].copy()
        if df_gen.empty:
            continue
        
        # Identify strain label
        if 'species_with_id' in df_gen.columns:
            strain_label = df_gen['species_with_id'].fillna(df_gen.get('species', 'unknown'))
        else:
            strain_label = df_gen['species'].fillna('unknown') if 'species' in df_gen.columns else pd.Series(['unknown'] * len(df_gen))
        
        values = df_gen[metabolite]
        tested = values.notna() & (values != 0)  # tested if nonzero
        positives = values.isin(positive_values_map.get(cat, {1}))
        
        for lab, t, p in zip(strain_label, tested, positives):
            entry = species_rows.setdefault(lab, {c: (0, 0) for c in categories})
            # Accumulate if multiple rows for same strain
            prev_t, prev_p = entry[cat]
            entry[cat] = (prev_t + int(t), prev_p + int(p))
    
    if not species_rows:
        return None
    
    # Select top_k strains by total tested across categories
    strain_counts = {s: sum(v[c][0] for c in categories) for s, v in species_rows.items()}
    top_strains = sorted(strain_counts.keys(), key=lambda s: strain_counts[s], reverse=True)[:top_k]
    
    # Filter out strains with no tests
    top_strains = [s for s in top_strains if strain_counts[s] > 0]
    
    if not top_strains:
        return None
    
    # Build nodes
    nodes = []
    node_index = {}
    
    # Root genus node
    node_index['genus'] = 0
    nodes.append({'label': genus})
    
    # Strain nodes
    for s in top_strains:
        node_index[f"strain::{s}"] = len(nodes)
        # Truncate long strain names for display
        display_name = s if len(s) <= 40 else s[:37] + "..."
        nodes.append({'label': display_name})
    
    # Category nodes
    for c in categories:
        node_index[f"cat::{c}"] = len(nodes)
        nodes.append({'label': c.capitalize()})
    
    # Result nodes
    node_index['res::Positive'] = len(nodes)
    nodes.append({'label': 'Positive'})
    node_index['res::Negative'] = len(nodes)
    nodes.append({'label': 'Negative'})
    
    # Build links
    source, target, value = [], [], []
    
    # Genus -> strain (value = total tests for that strain)
    for s in top_strains:
        tests = strain_counts[s]
        if tests <= 0:
            continue
        source.append(node_index['genus'])
        target.append(node_index[f"strain::{s}"])
        value.append(tests)
    
    # Strain -> category (value = number of tests in that category)
    for s in top_strains:
        for c in categories:
            t, p = species_rows[s][c]
            if t > 0:
                source.append(node_index[f"strain::{s}"])
                target.append(node_index[f"cat::{c}"])
                value.append(t)
    
    # Category -> result (sum positives/negatives across strains)
    for c in categories:
        pos_sum = sum(species_rows[s][c][1] for s in top_strains)
        tested_sum = sum(species_rows[s][c][0] for s in top_strains)
        neg_sum = tested_sum - pos_sum
        
        if pos_sum > 0:
            source.append(node_index[f"cat::{c}"])
            target.append(node_index['res::Positive'])
            value.append(pos_sum)
        
        if neg_sum > 0:
            source.append(node_index[f"cat::{c}"])
            target.append(node_index['res::Negative'])
            value.append(neg_sum)
    
    # Build Sankey figure
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            label=[n['label'] for n in nodes],
            pad=15,
            thickness=15,
            line=dict(color="#000000", width=0.5),
            color=[
                '#1f77b4',  # genus - blue
                *['#2ca02c' for _ in top_strains],  # strains - green
                '#ff7f0e', '#9467bd', '#d62728', '#8c564b',  # categories - various
                '#17becf', '#bcbd22'  # results - cyan, olive
            ][:len(nodes)]
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color='rgba(200,200,200,0.4)'
        ),
        textfont=dict(
            family="Arial, sans-serif",
            size=18,
            color="#FFFFFF"
        )
    )])
    
    fig.update_layout(
        title_text=f"Sankey: {genus} -> top {len(top_strains)} strains -> categories -> results<br><sub>Metabolite: {metabolite}</sub>",
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
    
    return fig
