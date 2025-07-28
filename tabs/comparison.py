import streamlit as st
import pandas as pd
import numpy as np
import re
from collections import Counter, defaultdict
from itertools import combinations
from collections import Counter
import io
from typing import Tuple, Dict, List
import random

##### loading in our data files
##### need to replace if we get new files.
df_pos_res_nostrain = pd.read_csv("data_files/step4_positively_tested_by_genera_res_nostrain.csv")
df_pos_res_yesstrain = pd.read_csv("data_files/step4_positively_tested_by_genera_res_yesstrain.csv")
df_pos_sen_nostrain = pd.read_csv("data_files/step4_positively_tested_by_genera_sen_nostrain.csv")
df_pos_sen_yesstrain = pd.read_csv("data_files/step4_positively_tested_by_genera_sen_yesstrain.csv")

@st.cache_data

def load_category_data(strain_suffix: str) -> dict[str, pd.DataFrame]:
    """Load the four genus×metabolite tables  for the given  suffix."""
    base = "data_files/step4_positively_tested_by_genera_"
    return {
        "Pos Resistance": pd.read_csv(f"{base}res_{strain_suffix}.csv"),
        "Pos Sensitive": pd.read_csv(f"{base}sen_{strain_suffix}.csv"),
        "Pos Utilization": pd.read_csv(f"{base}util_{strain_suffix}.csv"),
        "Pos Production": pd.read_csv(f"{base}prod_{strain_suffix}.csv"),
    }

def display(data_frames):
    st.title("Comparison of Genera")
    st.write("Select up to 10 genera and see which metabolites fall into each category.")

    # 1) Choose strain
    strain_option = st.radio("Strain Option", ["Isolates","Strains"])
    suffix = "nostrain" if strain_option=="Isolates" else "yesstrain"
    dfs = load_category_data(suffix)

    # 2) Genus selector
    all_genera = sorted(dfs["Pos Resistance"]["genus"].unique())
    selected = st.multiselect("Select genera (max 10)", all_genera, [])
    if not selected or len(selected)>10:
        st.info("Pick 1–10 genera."); return

    # helper to get metabolites (exclude species_count!)
    def mets_for(genus: str, df: pd.DataFrame) -> set[str]:
        row = df.loc[df["genus"]==genus]
        if row.empty: return set()
        row = row.iloc[0]
        return {col for col in df.columns 
                if col not in ("genus","species_count") and row[col]>0}

    # 3) Build per-genus × category table
    table_data = {cat: [] for cat in dfs}
    index_labels = []
    for genus in selected:
        # grab species_count from any df (they should all agree)
        count = int(dfs["Pos Resistance"].loc[
            dfs["Pos Resistance"]["genus"]==genus, "species_count"
        ].iloc[0])
        index_labels.append(f"{genus} ({count})")

        for cat, df in dfs.items():
            mets = mets_for(genus, df)
            table_data[cat].append(", ".join(sorted(mets)) or "—")

    # 5) display the df
    out = pd.DataFrame(table_data, index=index_labels)
    out.index.name = "Genus (species count)"
    st.dataframe(out, use_container_width=True)
    

    # 4) Shared & Synergy summary (also skipping species_count)
    # shared = intersection across all selected
    shared_summary: dict[str, str] = {}
    for cat, df in dfs.items():
        sets = [mets_for(g, df) for g in selected]
        shared = set.intersection(*sets) if sets else set()
        header = cat.replace("Pos ", "Shared ")
        shared_summary[header] = ", ".join(sorted(shared)) or "—"

    # prod/util synergy
        # --- Build a metabolite → {producers, utilizers} map ---
    from collections import defaultdict

    prod_df = dfs["Pos Production"]
    util_df = dfs["Pos Utilization"]

    synergy_map = defaultdict(lambda: {"prod": set(), "util": set()})

    # collect producers
    for genus in selected:
        for met in mets_for(genus, prod_df):
            synergy_map[met]["prod"].add(genus)

    # collect utilizers
    for genus in selected:
        for met in mets_for(genus, util_df):
            synergy_map[met]["util"].add(genus)

    # now format only those with both non‐empty prod & util
    synergy_entries = []
    for met, groups in synergy_map.items():
        prods = groups["prod"]
        utils = groups["util"]
        if prods and utils:
            prod_list = ", ".join(sorted(prods))
            util_list = ", ".join(sorted(utils))
            synergy_entries.append(
                f"{met} (prod: {prod_list}; util: {util_list})"
            )

    shared_summary["Prod/Util Synergy"] = "; ".join(synergy_entries) or "—"


    summary_df = pd.DataFrame([shared_summary], index=["Summary"])
    st.subheader("Shared & Synergy Summary")

    st.table(summary_df)
    csv = summary_df.to_csv(index=True).encode("utf-8")
    st.download_button(
        "Download Shared & Synergy Summary",
        csv,
        file_name="shared_synergy_summary.csv",
        mime="text/csv",
    )


