# make sure you have streamlit. pip install streamlit
import streamlit as st
import pandas as pd
import numpy as np
import os
import re
import json
import csv

def main():
    st.title("Visualizing BacDive")
    st.write("Work in progress.")

    # Set up the sidebar navigation
    st.sidebar.title("Navigation")
    tab = st.sidebar.radio("Go to", ["General Overview", "Circos", "Trends"])

    # Display content based on selected tab
    if tab == "General Overview":
        st.title("General Overview")
        st.write("This section will show a summary of the antibiotic data.")

    elif tab == "Circos":
        st.title("Circos Visualization")
        st.write("This section will include a Circos plot for antibiotic relationships.")

    elif tab == "Trends":
        st.title("Trends Over Time")
        st.write("This section will display trends in antibiotic resistance/utilization over time.")

main()



