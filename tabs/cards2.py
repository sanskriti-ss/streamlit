
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
def get_genus_ranking(data_frames, strain_option='Isolate'):
    prod_key = f"step4_positively_tested_by_genera_prod_{'isolate' if strain_option=='Isolate' else 'strain'}.csv"
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
def get_unique_mets_ranking(data_frames, cat, strain_option='Isolate'):
    cat_map = {'Production':'prod','Utilization':'util','Resistance':'res','Sensitivity':'sen'}
    if cat not in cat_map:
        return {}
    key = cat_map[cat]
    file_key = f"step4_positively_tested_by_genera_{key}_{'isolate' if strain_option=='Isolate' else 'strain'}.csv"
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

    base_key = next((f for f in data_frames if 'isolate' in f), None)
    if not base_key:
        st.error('No data available for cards.')
        return
    df_base = data_frames[base_key]
    genera = sorted(df_base['genus'].unique())

    # Compute rankings
    species_no = get_genus_ranking(data_frames, 'Isolate')
    species_yes = get_genus_ranking(data_frames, 'Strain')
    unique_no = {cat: get_unique_mets_ranking(data_frames, cat, 'Isolate') for cat in ['Production','Utilization','Resistance','Sensitivity']}
    unique_yes = {cat: get_unique_mets_ranking(data_frames, cat, 'Strain') for cat in ['Production','Utilization','Resistance','Sensitivity']}

    cat_map = {'Production':'prod','Utilization':'util','Resistance':'res','Sensitivity':'sen'}

    # Layout
    cols_per_row = 2
    for i in range(0, len(genera), cols_per_row):
        cols = st.columns(cols_per_row)
        for j, genus in enumerate(genera[i:i+cols_per_row]):
            with cols[j]:
                with st.expander(genus):
                    # Isolate block
                    st.write('**Isolates**')
                    df_pn = data_frames.get('step4_positively_tested_by_genera_prod_isolate.csv', pd.DataFrame())
                    r = df_pn[df_pn['genus']==genus]
                    sp_n = r['species_count'].iloc[0] if not r.empty else 'N/A'
                    rk_n = species_no.get(genus, 'N/A')
                    st.write(f"**# of species:** {sp_n} (#{rk_n})")

                    g1, g2 = st.columns(2), st.columns(2)
                    for idx, cat in enumerate(['Production','Utilization','Resistance','Sensitivity']):
                        key = cat_map[cat]
                        grp = g1[idx] if idx<2 else g2[idx-2]
                        with grp:
                            df_c = data_frames.get(f'step4_positively_tested_by_genera_{key}_isolate.csv', pd.DataFrame())
                            rc = df_c[df_c['genus']==genus]
                            uc = (rc[rc.columns.difference(['genus','species_count'])]>0).sum(axis=1).iloc[0] if not rc.empty else 0
                            rv = unique_no[cat].get(genus, '')
                            st.markdown(render_shape(uc, rv, key), unsafe_allow_html=True)

                    st.markdown("<hr style='margin:5px 0;'>", unsafe_allow_html=True)

                    # Yes Strains block
                    st.write('**Strains**')
                    df_py = data_frames.get('step4_positively_tested_by_genera_prod_strain.csv', pd.DataFrame())
                    r2 = df_py[df_py['genus']==genus]
                    sp_y = r2['species_count'].iloc[0] if not r2.empty else 'N/A'
                    rk_y = species_yes.get(genus, 'N/A')
                    st.write(f"**# of species:** {sp_y} (#{rk_y})")

                    g1, g2 = st.columns(2), st.columns(2)
                    for idx, cat in enumerate(['Production','Utilization','Resistance','Sensitivity']):
                        key = cat_map[cat]
                        grp = g1[idx] if idx<2 else g2[idx-2]
                        with grp:
                            df_c = data_frames.get(f'step4_positively_tested_by_genera_{key}_strain.csv', pd.DataFrame())
                            rc = df_c[df_c['genus']==genus]
                            uc = (rc[rc.columns.difference(['genus','species_count'])]>0).sum(axis=1).iloc[0] if not rc.empty else 0
                            rv = unique_yes[cat].get(genus, '')
                            st.markdown(render_shape(uc, rv, key), unsafe_allow_html=True)
