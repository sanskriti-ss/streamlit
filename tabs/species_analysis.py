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
