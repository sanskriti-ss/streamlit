# make sure you have streamlit. pip install streamlit
import streamlit as st
import seaborn as sns
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import textwrap

def main():
    st.title("Visualizing BacDive")
    st.write("Work in progress.")

    # Load CSV files from the data_files subfolder
    data_folder = 'data_files'
    csv_files = [f for f in os.listdir(data_folder) if f.endswith('.csv')]
    data_frames = {file: pd.read_csv(os.path.join(data_folder, file)) for file in csv_files}

    if data_frames:
        st.write("CSV files have been successfully loaded.")

    # Set up the sidebar navigation
    st.sidebar.title("Navigation")
    tab = st.sidebar.radio("Go to", ["General Overview", "Circos", "Trends"])

    # Display content based on selected tab
    if tab == "General Overview":
        st.title("General Overview")
        st.write("This section will show a summary of the antibiotic data.")

        # Dropdown for file selection
        selected_file = st.selectbox("Select File:", list(data_frames.keys()))

        # Sliders for number of genera and metabolites
        num_genera = st.slider("Number of Genera:", min_value=5, max_value=20, value=10)
        num_metabolites = st.slider("Number of Metabolites:", min_value=5, max_value=40, value=20)

        # Function to create a heatmap by numbers
        def create_heatmap_by_numbers(selected_file, num_genera, num_metabolites):
            df = data_frames[selected_file]

            # Exclude 'species_count' and define metabolite columns
            metabolite_columns = df.columns.difference(['genus', 'species_count']).tolist()

            # Convert metabolite columns to numeric, coerce errors to NaN
            df[metabolite_columns] = df[metabolite_columns].apply(pd.to_numeric, errors='coerce')

            # Select top metabolites and genera based on user input
            top_metabolites = df[metabolite_columns].sum().nlargest(num_metabolites).index
            top_genera = df[metabolite_columns].sum(axis=1).nlargest(num_genera).index

            # Wrap the x-axis labels
            wrapped_labels = [textwrap.fill(label, width=20) for label in top_metabolites]

            heatmap_data = df.iloc[top_genera][top_metabolites]
            heatmap_data = heatmap_data.set_index(df['genus'].iloc[top_genera])

            # Apply log transformation for heatmap color scale
            heatmap_data_log = heatmap_data.applymap(lambda x: np.log10(x + 1) if x > 0 else 0)

            # Create the heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(heatmap_data_log, annot=heatmap_data.astype(int), fmt='d', cmap="viridis", cbar=True, linewidths=.5,
                        annot_kws={"size": 8, "color": "white"})

            # Title with the selected number of genera and metabolites
            title = f'Heatmap of Top {num_metabolites} Metabolites and Top {num_genera} Genera for {selected_file}'
            plt.title(title)
            plt.xlabel('Top Metabolites')
            plt.ylabel('Top Genera')
            plt.xticks(ticks=np.arange(len(wrapped_labels)) + 0.5, labels=wrapped_labels, rotation=90, ha='center', fontsize=8)
            plt.yticks(rotation=0)

            st.pyplot(plt)

        # Function to create a heatmap by proportions
        def create_heatmap_by_proportions(selected_file, num_genera, num_metabolites):
            df = data_frames[selected_file]

            # Exclude 'species_count' and define metabolite columns
            metabolite_columns = df.columns.difference(['genus', 'species_count']).tolist()

            # Convert metabolite columns and species_count to numeric, coerce errors to NaN
            df[metabolite_columns] = df[metabolite_columns].apply(pd.to_numeric, errors='coerce')
            df['species_count'] = pd.to_numeric(df['species_count'], errors='coerce')

            # Calculate proportions
            df[metabolite_columns] = df[metabolite_columns].div(df['species_count'], axis=0)

            # Select top metabolites and genera based on user input
            top_metabolites = df[metabolite_columns].sum().nlargest(num_metabolites).index
            top_genera = df[metabolite_columns].sum(axis=1).nlargest(num_genera).index

            # Wrap the x-axis labels
            wrapped_labels = [textwrap.fill(label, width=20) for label in top_metabolites]

            heatmap_data = df.iloc[top_genera][top_metabolites]
            heatmap_data = heatmap_data.set_index(df['genus'].iloc[top_genera])

            # Create the heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(heatmap_data, annot=True, fmt='.2f', cmap="viridis", cbar=True, linewidths=.5,
                        annot_kws={"size": 8, "color": "white"})

            # Title with the selected number of genera and metabolites
            title = f'Heatmap of Top {num_metabolites} Metabolites and Top {num_genera} Genera Proportions for {selected_file}'
            plt.title(title)
            plt.xlabel('Top Metabolites')
            plt.ylabel('Top Genera')
            plt.xticks(ticks=np.arange(len(wrapped_labels)) + 0.5, labels=wrapped_labels, rotation=90, ha='center', fontsize=8)
            plt.yticks(rotation=0)

            st.pyplot(plt)

        # Button to generate the heatmap by numbers
        st.subheader("By Numbers")
        if st.button("Generate Heatmap by Numbers"):
            create_heatmap_by_numbers(selected_file, num_genera, num_metabolites)

        # Divider
        st.markdown("---")

        # Button to generate the heatmap by proportions
        st.subheader("By Proportions")
        if st.button("Generate Heatmap by Proportions"):
            create_heatmap_by_proportions(selected_file, num_genera, num_metabolites)

    elif tab == "Circos":
        st.title("Circos Visualization")
        st.write("This section will include a Circos plot for antibiotic relationships.")

    elif tab == "Trends":
        st.title("Trends Over Time")
        st.write("This section will display trends in antibiotic resistance/utilization over time.")

if __name__ == "__main__":
    main()
