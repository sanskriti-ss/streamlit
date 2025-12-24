"""
Species Analysis Utilities
Helper functions for analyzing species metabolite utilization data.
"""

import pandas as pd
import zipfile
import os
from typing import Dict, List, Tuple, Optional

def process_multiple_files(file_paths: List[str], file_types: List[str] = None) -> Dict[str, pd.DataFrame]:
    """
    Process multiple species data files and combine statistics.
    
    Args:
        file_paths (List[str]): List of file paths to process
        file_types (List[str]): Optional list of file type labels (e.g., ['prod', 'res', 'sen'])
        
    Returns:
        Dict[str, pd.DataFrame]: Dictionary with combined results
    """
    combined_stats = []
    all_species_data = {}
    
    for i, file_path in enumerate(file_paths):
        file_type = file_types[i] if file_types and i < len(file_types) else f"file_{i+1}"
        
        # Load data
        df = load_species_data_from_path(file_path)
        if df is None:
            continue
            
        # Calculate stats with activity type
        stats = calculate_species_stats(df, file_type)
        if stats is None:
            continue
            
        # Add file type information
        stats['file_type'] = file_type
        combined_stats.append(stats)
        all_species_data[file_type] = df
    
    if not combined_stats:
        return None
    
    # Combine all statistics
    combined_df = pd.concat(combined_stats, ignore_index=True)
    
    return {
        'combined_stats': combined_df,
        'individual_data': all_species_data,
        'top_species_overall': get_top_species_across_files(combined_df)
    }

def load_species_data_from_path(file_path: str) -> Optional[pd.DataFrame]:
    """
    Load species data from various file formats.
    
    Args:
        file_path (str): Path to the data file
        
    Returns:
        pd.DataFrame: Loaded data or None if failed
    """
    try:
        if file_path.endswith('.zip'):
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                csv_files = [f for f in zip_ref.namelist() if f.endswith('.csv')]
                if not csv_files:
                    return None
                
                with zip_ref.open(csv_files[0]) as csv_file:
                    df = pd.read_csv(csv_file)
        elif file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        else:
            return None
            
        return df
    except Exception as e:
        print(f"Error loading {file_path}: {str(e)}")
        return None

def calculate_species_stats(df: pd.DataFrame, activity_type: str = "utilization") -> Optional[pd.DataFrame]:
    """
    Calculate utilization statistics for species in a dataframe.
    
    Args:
        df (pd.DataFrame): Species data
        activity_type (str): Type of activity (utilization, resistance, production, sensitivity)
        
    Returns:
        pd.DataFrame: Statistics dataframe
    """
    if df is None or df.empty:
        return None
    
    # Debug: Print available columns for troubleshooting
    print(f"DEBUG - File columns: {list(df.columns)}")
    print(f"DEBUG - First few rows of key columns:")
    key_cols = ['BacID', 'species', 'genus', 'order', 'type_strain']
    existing_key_cols = [col for col in key_cols if col in df.columns]
    if existing_key_cols:
        print(df[existing_key_cols].head(3))
    else:
        print("No expected key columns found!")
        print("First 5 columns of the dataframe:")
        print(df.iloc[:3, :5])
    
    # Identify columns - check for both 'type_strain' and 'is_strain'
    # Also check for common variations in column names
    metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain']
    if 'is_strain' in df.columns:
        metadata_cols.append('is_strain')
    
    # Check for actual columns that exist in the dataframe
    actual_metadata_cols = [col for col in metadata_cols if col in df.columns]
    metabolite_cols = [col for col in df.columns if col not in actual_metadata_cols]
    
    if not metabolite_cols:
        return None
    

    # Process data - remember don't convert -1 to 0 for resistance/sensitivity data
    df_processed = df.copy()
    
    # Only clip values for production/utilization data
    if not ('res' in activity_type.lower() or 'resistance' in activity_type.lower() or 
            'sen' in activity_type.lower() or 'sensitivity' in activity_type.lower()):
        for col in metabolite_cols:
            df_processed[col] = df_processed[col].replace(-1, 0).clip(0, 1)
    # For resistance/sensitivity data, keep -1 values as they represent negative test results

    # Calculate statistics
    stats_list = []

    # Determine appropriate column names based on activity type
    if 'res' in activity_type.lower() or 'resistance' in activity_type.lower():
        activity_col = 'metabolites_resistant'
        rate_col = 'resistance_rate'
        activity_name = 'resistant to'
    elif 'prod' in activity_type.lower() or 'production' in activity_type.lower():
        activity_col = 'metabolites_produced'
        rate_col = 'production_rate'
        activity_name = 'produced'
    elif 'sen' in activity_type.lower() or 'sensitivity' in activity_type.lower():
        activity_col = 'metabolites_sensitive'
        rate_col = 'sensitivity_rate'
        activity_name = 'sensitive to'
    else:  # utilization or default
        activity_col = 'metabolites_utilized'
        rate_col = 'utilization_rate'
        activity_name = 'utilized'

    for _, row in df_processed.iterrows():
        # Count positives and tested based on activity type
        if 'res' in activity_type.lower() or 'resistance' in activity_type.lower():
            # 1 = resistant; tested = {1, -1}
            metabolites_with_activity = sum(1 for col in metabolite_cols if pd.notna(row[col]) and row[col] == 1)
            metabolites_tested = sum(1 for col in metabolite_cols if pd.notna(row[col]) and row[col] in (1, -1))
        elif 'sen' in activity_type.lower() or 'sensitivity' in activity_type.lower():
            # 1 = sensitive; tested = {1, -1}
            metabolites_with_activity = sum(1 for col in metabolite_cols if pd.notna(row[col]) and row[col] == 1)
            metabolites_tested = sum(1 for col in metabolite_cols if pd.notna(row[col]) and row[col] in (1, -1))
        else:
            # production/utilization: 1 = positive; tested = any non-NaN (counts 0 and 1)
            metabolites_with_activity = sum(1 for col in metabolite_cols if pd.notna(row[col]) and row[col] == 1)
            metabolites_tested = sum(1 for col in metabolite_cols if pd.notna(row[col]))

        # Rate uses metabolites_tested (prevents 100% when zeros exist)
        activity_rate = (metabolites_with_activity / metabolites_tested * 100) if metabolites_tested > 0 else 0

        # Strain/type_strain mapping (unchanged, but normalized)
        strain_info = 'unknown'
        if 'type_strain' in df_processed.columns and pd.notna(row['type_strain']):
            ts = str(row['type_strain']).strip().lower()
            if ts == 'yes':
                strain_info = 'Type Strain'
            elif ts == 'no':
                strain_info = 'Strain'
            else:
                strain_info = str(row['type_strain'])
        elif 'is_strain' in df_processed.columns and pd.notna(row.get('is_strain')):
            strain_info = 'Strain' if row['is_strain'] == 1 else 'Isolate'

        stats_list.append({
            'BacID': row.get('BacID'),
            'species': row.get('species'),
            'species_with_id': f"{row.get('species')} (ID: {row.get('BacID')})" if 'BacID' in df_processed.columns else row.get('species'),
            'genus': row.get('genus'),
            'order': row.get('order'),
            'type_strain': strain_info,
            activity_col: metabolites_with_activity,
            'metabolites_tested': metabolites_tested,
            'total_metabolites': len(metabolite_cols),  # keep for reference; don’t use for rate
            rate_col: activity_rate,
            'activity_type': activity_type.capitalize()
        })
    
    return pd.DataFrame(stats_list)

def get_top_species_across_files(combined_stats: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    """
    Get top species by total utilization across all file types.
    
    Args:
        combined_stats (pd.DataFrame): Combined statistics from multiple files
        top_n (int): Number of top species to return
        
    Returns:
        pd.DataFrame: Top species with aggregated statistics
    """
    # Group by species and aggregate
    species_totals = combined_stats.groupby(['species', 'genus', 'order']).agg({
        'metabolites_utilized': 'sum',
        'metabolites_tested': 'sum',
        'total_metabolites': 'sum',
        'file_type': lambda x: ', '.join(x.unique())
    }).reset_index()
    
    # Recalculate utilization rate
    species_totals['overall_utilization_rate'] = (
        species_totals['metabolites_utilized'] / species_totals['total_metabolites'] * 100
    )
    
    # Sort and return top N
    return species_totals.nlargest(top_n, 'metabolites_utilized')

def create_comparison_table(results: Dict[str, pd.DataFrame], top_n: int = 5) -> pd.DataFrame:
    """
    Create a comparison table showing top species across different file types.
    
    Args:
        results (Dict): Results from process_multiple_files
        top_n (int): Number of top species per file type
        
    Returns:
        pd.DataFrame: Comparison table
    """
    if not results or 'combined_stats' not in results:
        return pd.DataFrame()
    
    combined_stats = results['combined_stats']
    comparison_data = []
    
    # Get top species for each file type
    for file_type in combined_stats['file_type'].unique():
        file_stats = combined_stats[combined_stats['file_type'] == file_type]
        top_species = file_stats.nlargest(top_n, 'metabolites_utilized')
        
        for i, (_, row) in enumerate(top_species.iterrows(), 1):
            comparison_data.append({
                'File Type': file_type.upper(),
                'Rank': i,
                'Species': row['species'],
                'Genus': row['genus'],
                'Order': row['order'],
                'Metabolites Utilized': row['metabolites_utilized'],
                'Utilization Rate (%)': round(row['utilization_rate'], 2)
            })
    
    return pd.DataFrame(comparison_data)

def get_unique_metabolites_by_species(df: pd.DataFrame, species_name: str) -> List[str]:
    """
    Get list of metabolites utilized by a specific species.
    
    Args:
        df (pd.DataFrame): Species data
        species_name (str): Name of the species
        
    Returns:
        List[str]: List of utilized metabolites
    """
    species_data = df[df['species'] == species_name]
    if species_data.empty:
        return []
    
    metadata_cols = ['BacID', 'species', 'genus', 'order', 'type_strain']
    metabolite_cols = [col for col in df.columns if col not in metadata_cols]
    
    utilized_metabolites = []
    for _, row in species_data.iterrows():
        for metabolite in metabolite_cols:
            if pd.notna(row[metabolite]) and row[metabolite] == 1:
                utilized_metabolites.append(metabolite)
    
    return list(set(utilized_metabolites))
