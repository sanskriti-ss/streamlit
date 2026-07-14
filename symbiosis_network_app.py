"""
Standalone launcher for the symbiosis network tab.

Prefer the main app:
    streamlit run app.py
then open **Symbiosis Network** in the sidebar.

This entrypoint remains for quick standalone demos:
    streamlit run symbiosis_network_app.py
"""

from __future__ import annotations

import streamlit as st

from tabs.symbiosis_network import display

st.set_page_config(page_title="Symbiosis Network", layout="wide")
display()
