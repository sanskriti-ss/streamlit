<<<<<<< Updated upstream
import streamlit as st
import pandas as pd
from pathlib import Path
from PIL import Image

# point to your icons folder
ICON_DIR = Path(__file__).parent / "icons"

ICON_PATHS = {
    "Production": ICON_DIR / "production.png",
    "Utilization": ICON_DIR / "utilization.png",
    "Resistance": ICON_DIR / "resistance.png",
    "Sensitivity": ICON_DIR / "sensitivity.png",
}

# Helper: Render a rectangle with given count and rank.
def render_shape(count, rank, color, width=80, height=40):
    # Determine border style based on rank.
    if isinstance(rank, int):
        if rank <= 10:
            border = "2px solid black"
        elif rank <= 100:
            border = "2px solid grey"
        else:
            border = "none"
    else:
        border = "none"
    
    # Use a wrapper div for spacing.
    wrapper_style = "display:inline-block; margin:6px;"
    
    html = (
        f"<div style='{wrapper_style} background-color: {color}; width:{width}px; height:{height}px; "
        f"border: {border}; text-align: center; line-height: {height}px; border-radius: 4px; "
        f"color: white; font-size: 14px;'>"
        f"{count} #{rank}"
        f"</div>"
    )
    return html

# Cached helper to compute ranking based on species_count from the Production file.
@st.cache_data
def get_genus_ranking(data_frames, strain_option="Isolate"):
    prod_key = f"step4_positively_tested_by_genera_prod_{'nostrain' if strain_option=='Isolate' else 'yesstrain'}.csv"
    if prod_key not in data_frames:
        return {}
    df_prod = data_frames[prod_key]
    ranking_df = df_prod[['genus', 'species_count']].dropna().sort_values(by='species_count', ascending=False)
    # Use standard competition ranking: method='min'
    ranking_df['rank'] = ranking_df['species_count'].rank(method='min', ascending=False).astype(int)
    return dict(zip(ranking_df['genus'], ranking_df['rank']))

# Cached helper to compute unique metabolite counts ranking for a given category.
@st.cache_data
def get_unique_mets_ranking(data_frames, cat, strain_option="Isolate"):
    cat_mapping = {"Production": "prod", "Utilization": "util", "Resistance": "res", "Sensitivity": "sen"}
    if cat not in cat_mapping:
        return {}
    file_key = f"step4_positively_tested_by_genera_{cat_mapping[cat]}_{'nostrain' if strain_option=='Isolate' else 'yesstrain'}.csv"
    if file_key not in data_frames:
        return {}
    df_cat = data_frames[file_key].copy()
    metabolite_cols = df_cat.columns.difference(["genus", "species_count"])
    df_cat['unique_mets'] = (df_cat[metabolite_cols] > 0).sum(axis=1)
    ranking_df = df_cat[['genus', 'unique_mets']].dropna().sort_values(by='unique_mets', ascending=False)
    # Use standard competition ranking.
    ranking_df['rank'] = ranking_df['unique_mets'].rank(method='min', ascending=False).astype(int)
    return dict(zip(ranking_df['genus'], ranking_df['rank']))

def display(data_frames):
    st.header("Cards")
    st.write("Click on a card below to view details for each genus. (Details will be added soon.)")
    
    # Use one CSV file (with "nostrain") to get the list of genera.
    file_key = next((fname for fname in data_frames if "nostrain" in fname), None)
    if not file_key:
        st.error("No data available for cards.")
        return
    
    df = data_frames[file_key]
    genera = sorted(df["genus"].unique())
    
    # Compute rankings for each strain option separately.
    species_ranking_no = get_genus_ranking(data_frames, strain_option="Isolate")
    species_ranking_yes = get_genus_ranking(data_frames, strain_option="Strain")
    
    unique_rankings_no = {}
    unique_rankings_yes = {}
    for cat in ["Production", "Utilization", "Resistance", "Sensitivity"]:
        unique_rankings_no[cat] = get_unique_mets_ranking(data_frames, cat, strain_option="Isolate")
        unique_rankings_yes[cat] = get_unique_mets_ranking(data_frames, cat, strain_option="Strain")
    
    # Updated pastel color mapping.
    pastel_colors = {
        "Production": "#A8E6CF",    # pastel green
        "Utilization": "#FF8B94",    # pastel red
        "Resistance": "#89CFF0",     # pastel blue (more vibrant)
        "Sensitivity": "#CDB4DB"     # pastel purple
    }
    
    # Arrange cards in a grid with 2 cards per row.
    num_cols = 2
    rows = [genera[i:i+num_cols] for i in range(0, len(genera), num_cols)]
    
    # File naming parameters.
    prod_key_no = "step4_positively_tested_by_genera_prod_nostrain.csv"
    prod_key_yes = "step4_positively_tested_by_genera_prod_yesstrain.csv"
    cat_mapping = {"Production": "prod", "Utilization": "util", "Resistance": "res", "Sensitivity": "sen"}
    
    for row in rows:
        cols = st.columns(num_cols)
        for idx, genus in enumerate(row):
            with cols[idx]:
                with st.expander(genus, expanded=False):
                    st.write("**No Strain**")
                    # Species count for No Strains.
                    if prod_key_no in data_frames:
                        df_prod_no = data_frames[prod_key_no]
                        row_prod_no = df_prod_no[df_prod_no["genus"] == genus]
                        species_count_no = row_prod_no["species_count"].iloc[0] if not row_prod_no.empty else "N/A"
                    else:
                        species_count_no = "N/A"
                    rank_no = species_ranking_no.get(genus, "N/A")
                    st.write(f"**# of species:** {species_count_no} (#{rank_no})")
                    
                    # Render unique metabolite counts in a 2x2 grid for No Strains.
                    cols_grid1 = st.columns(2)
                    cols_grid2 = st.columns(2)
                    
                    # Production.
                    with cols_grid1[0]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Production']}_nostrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_no.get("Production", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Production"]), unsafe_allow_html=True)
                    
                    # Utilization.
                    with cols_grid1[1]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Utilization']}_nostrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_no.get("Utilization", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Utilization"]), unsafe_allow_html=True)
                    
                    # Resistance.
                    with cols_grid2[0]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Resistance']}_nostrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_no.get("Resistance", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Resistance"]), unsafe_allow_html=True)
                    
                    # Sensitivity.
                    with cols_grid2[1]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Sensitivity']}_nostrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_no.get("Sensitivity", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Sensitivity"]), unsafe_allow_html=True)
                    
                    st.markdown("<hr style='margin-top:5px; margin-bottom:5px;'>", unsafe_allow_html=True)
                    
                    st.write("**Yes Strain**")
                    # Species count for Yes Strains.
                    if prod_key_yes in data_frames:
                        df_prod_yes = data_frames[prod_key_yes]
                        row_prod_yes = df_prod_yes[df_prod_yes["genus"] == genus]
                        species_count_yes = row_prod_yes["species_count"].iloc[0] if not row_prod_yes.empty else "N/A"
                    else:
                        species_count_yes = "N/A"
                    rank_yes = species_ranking_yes.get(genus, "N/A")
                    st.write(f"**# of species:** {species_count_yes} (#{rank_yes})")
                    
                    # Grid for Yes Strain unique metabolite counts.
                    cols_grid1 = st.columns(2)
                    cols_grid2 = st.columns(2)
                    
                    # Production.
                    with cols_grid1[0]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Production']}_yesstrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_yes.get("Production", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Production"]), unsafe_allow_html=True)
                    
                    # Utilization.
                    with cols_grid1[1]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Utilization']}_yesstrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_yes.get("Utilization", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Utilization"]), unsafe_allow_html=True)
                    
                    # Resistance.
                    with cols_grid2[0]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Resistance']}_yesstrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_yes.get("Resistance", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Resistance"]), unsafe_allow_html=True)
                    
                    # Sensitivity.
                    with cols_grid2[1]:
                        file_key_cat = f"step4_positively_tested_by_genera_{cat_mapping['Sensitivity']}_yesstrain.csv"
                        if file_key_cat in data_frames:
                            df_cat = data_frames[file_key_cat]
                            row_cat = df_cat[df_cat["genus"] == genus]
                            unique_count = (row_cat[row_cat.columns.difference(["genus", "species_count"])] > 0).sum(axis=1).iloc[0] if not row_cat.empty else 0
                        else:
                            unique_count = "N/A"
                        rank_val = unique_rankings_yes.get("Sensitivity", {}).get(genus, "")
                        st.markdown(render_shape(unique_count, rank_val if rank_val else "", pastel_colors["Sensitivity"]), unsafe_allow_html=True)
=======
import streamlit as st
import pandas as pd
import base64
import os

# Helper: load icon as base64 string (cached)
@st.cache_data
def get_icon_base64(type_key: str, level: str) -> str:
    path = os.path.join('icons', f"{type_key}_{level}.png")
    if not os.path.exists(path):
        return ''
    with open(path, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode('utf-8')

# Helper: Render a "card" with inline background-image icon
def render_shape(count, rank, type_key, width=80, height=80):
    # Determine level based on rank
    if isinstance(rank, int):
        if rank <= 10:
            level = 'high'
        elif rank <= 100:
            level = 'mid'
        else:
            level = 'low'
    else:
        level = 'low'

    b64 = get_icon_base64(type_key, level)
    if b64:
        bg_style = f"background-image: url('data:image/png;base64,{b64}');"
    else:
        bg_style = 'background-color: #444;'

    wrapper_style = (
        "display:inline-block; margin:6px; "
        f"width:{width}px; height:{height}px; "
        f"{bg_style} background-size: contain; background-repeat: no-repeat; background-position: center; "
        "text-align:center; "
        f"line-height:{height}px; color:white; font-size:14px;"
    )

    html = f"<div style=\"{wrapper_style}\">{count} #{rank}</div>"
    return html

# Cached helper to compute genus ranking
@st.cache_data
def get_genus_ranking(data_frames, strain_option='No Strains'):
    prod_key = f"step4_positively_tested_by_genera_prod_{'nostrain' if strain_option=='No Strains' else 'yesstrain'}.csv"
    if prod_key not in data_frames:
        return {}
    df = data_frames[prod_key]
    ranking = (
        df[['genus','species_count']]
          .dropna()
          .sort_values('species_count', ascending=False)
    )
    ranking['rank'] = ranking['species_count'].rank(method='min', ascending=False).astype(int)
    return dict(zip(ranking['genus'], ranking['rank']))

# Cached helper to compute unique metabolite rankings
def get_unique_mets_ranking(data_frames, cat, strain_option='No Strains'):
    cat_map = {'Production':'prod','Utilization':'util','Resistance':'res','Sensitivity':'sen'}
    if cat not in cat_map:
        return {}
    key = cat_map[cat]
    file_key = f"step4_positively_tested_by_genera_{key}_{'nostrain' if strain_option=='No Strains' else 'yesstrain'}.csv"
    if file_key not in data_frames:
        return {}
    df = data_frames[file_key].copy()
    mets = df.columns.difference(['genus','species_count'])
    df['unique_mets'] = (df[mets] > 0).sum(axis=1)
    ranking = (
        df[['genus','unique_mets']]
          .dropna()
          .sort_values('unique_mets', ascending=False)
    )
    ranking['rank'] = ranking['unique_mets'].rank(method='min', ascending=False).astype(int)
    return dict(zip(ranking['genus'], ranking['rank']))

# Main display function

def display(data_frames):
    st.header('Cards')
    st.write('Click on a card below to view details for each genus.')

    base_key = next((f for f in data_frames if 'nostrain' in f), None)
    if not base_key:
        st.error('No data available for cards.')
        return
    df_base = data_frames[base_key]
    genera = sorted(df_base['genus'].unique())

    # Compute rankings
    species_no = get_genus_ranking(data_frames, 'No Strains')
    species_yes = get_genus_ranking(data_frames, 'Yes Strains')
    unique_no = {cat: get_unique_mets_ranking(data_frames, cat, 'No Strains') for cat in ['Production','Utilization','Resistance','Sensitivity']}
    unique_yes = {cat: get_unique_mets_ranking(data_frames, cat, 'Yes Strains') for cat in ['Production','Utilization','Resistance','Sensitivity']}

    cat_map = {'Production':'prod','Utilization':'util','Resistance':'res','Sensitivity':'sen'}

    # Layout
    cols_per_row = 2
    for i in range(0, len(genera), cols_per_row):
        cols = st.columns(cols_per_row)
        for j, genus in enumerate(genera[i:i+cols_per_row]):
            with cols[j]:
                with st.expander(genus):
                    # No Strains block
                    st.write('**No Strain**')
                    df_pn = data_frames.get('step4_positively_tested_by_genera_prod_nostrain.csv', pd.DataFrame())
                    r = df_pn[df_pn['genus']==genus]
                    sp_n = r['species_count'].iloc[0] if not r.empty else 'N/A'
                    rk_n = species_no.get(genus, 'N/A')
                    st.write(f"**# of species:** {sp_n} (#{rk_n})")

                    g1, g2 = st.columns(2), st.columns(2)
                    for idx, cat in enumerate(['Production','Utilization','Resistance','Sensitivity']):
                        key = cat_map[cat]
                        grp = g1[idx] if idx<2 else g2[idx-2]
                        with grp:
                            df_c = data_frames.get(f'step4_positively_tested_by_genera_{key}_nostrain.csv', pd.DataFrame())
                            rc = df_c[df_c['genus']==genus]
                            uc = (rc[rc.columns.difference(['genus','species_count'])]>0).sum(axis=1).iloc[0] if not rc.empty else 0
                            rv = unique_no[cat].get(genus, '')
                            st.markdown(render_shape(uc, rv, key), unsafe_allow_html=True)

                    st.markdown("<hr style='margin:5px 0;'>", unsafe_allow_html=True)

                    # Yes Strains block
                    st.write('**Yes Strain**')
                    df_py = data_frames.get('step4_positively_tested_by_genera_prod_yesstrain.csv', pd.DataFrame())
                    r2 = df_py[df_py['genus']==genus]
                    sp_y = r2['species_count'].iloc[0] if not r2.empty else 'N/A'
                    rk_y = species_yes.get(genus, 'N/A')
                    st.write(f"**# of species:** {sp_y} (#{rk_y})")

                    g1, g2 = st.columns(2), st.columns(2)
                    for idx, cat in enumerate(['Production','Utilization','Resistance','Sensitivity']):
                        key = cat_map[cat]
                        grp = g1[idx] if idx<2 else g2[idx-2]
                        with grp:
                            df_c = data_frames.get(f'step4_positively_tested_by_genera_{key}_yesstrain.csv', pd.DataFrame())
                            rc = df_c[df_c['genus']==genus]
                            uc = (rc[rc.columns.difference(['genus','species_count'])]>0).sum(axis=1).iloc[0] if not rc.empty else 0
                            rv = unique_yes[cat].get(genus, '')
                            st.markdown(render_shape(uc, rv, key), unsafe_allow_html=True)
>>>>>>> Stashed changes
