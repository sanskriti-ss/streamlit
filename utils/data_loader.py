import os
import pandas as pd

def load_data(data_folder):
    """Loads CSV files from the data folder and returns a dictionary of dataframes."""
    if not os.path.exists(data_folder):
        return {}

    csv_files = [f for f in os.listdir(data_folder) if f.endswith('.csv')]
    return {file: pd.read_csv(os.path.join(data_folder, file)) for file in csv_files}
