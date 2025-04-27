import streamlit as st
import pandas as pd
import numpy as np
import re
from collections import Counter, defaultdict
from itertools import combinations
from collections import Counter
import io
from typing import Tuple, Dict, List

##### loading in our data files
##### need to replace if we get new files.
df_pos_res_nostrain = pd.read_csv("data_files/step4_positively_tested_by_genera_res_nostrain.csv")
df_pos_res_yesstrain = pd.read_csv("data_files/step4_positively_tested_by_genera_res_yesstrain.csv")
df_pos_sen_nostrain = pd.read_csv("data_files/step4_positively_tested_by_genera_sen_nostrain.csv")
df_pos_sen_yesstrain = pd.read_csv("data_files/step4_positively_tested_by_genera_sen_yesstrain.csv")

@st.cache_data

##### generates files we can download from dictionary, for modularity
def write_final_file(file_name: str, processed_lines: list) -> bytes:
    buffer = io.StringIO()
    buffer.writelines(processed_lines)
    # Convert the in-memory text file to bytes (using UTF-8 encoding)
    return buffer.getvalue().encode("utf-8")

##### generates combinations from dictionary, used once, for modularity
def generate_combinations_from_dict(data_dict: dict, line_format: str = "black") -> list:
    lines = []
    for key, items in data_dict.items():
        for a, b in combinations(items, 2):
            line = f"{a} 100 0 {b} 100 0 color={line_format}"
            lines.append(line)
    return lines

##### generates combinations from df, used once, for modularity
def generate_combinations_from_df(df, line_format="black") -> list:
    lines = []
    for column in df.columns:
        # Drop NaN rows to get the actual elements in this column
        genera = df[column].dropna()
        # Build combinations
        for a, b in combinations(genera, 2):
            line = f"{a} 100 0 {b} 100 0 color={line_format}"
            lines.append(line)
    return lines

##### removes spaces from column names to match edges list
def crop_to_edges(values, primary_data_file, edges_list):
    st.write("**Debug**: Entering crop_to_edges")
    st.write(f"DataFrame type: {type(primary_data_file)}")

    if values['gen_or_met'] == "Metabolites":
        # Remove spaces from column names so they match 'edges_list'
        primary_data_file.columns = primary_data_file.columns.str.replace(' ', '', regex=False)

        # first column is 'genus'; skip it
        columns_to_keep = [col for col in primary_data_file.columns[1:] if col in edges_list]

        st.write("DataFrame columns:", primary_data_file.columns.tolist())
        st.write("Edges list:", edges_list)
        st.write("Columns to keep:", columns_to_keep)

        # crop columns
        cropped_df = primary_data_file[columns_to_keep]
        return cropped_df

    else:  # 'Genera'
        st.write("**Debug**: 'Genera' path in crop_to_edges")
        # The first column in your DataFrame is presumably the genus name
        # Keep only rows whose genus is in edges_list
        # (This assumes the genus column is the very first column index=0)
        cropped_df = primary_data_file[primary_data_file.iloc[:, 0].isin(edges_list)]
        st.write("Cropped DataFrame (Genera):")
        st.write(cropped_df.head())
        return cropped_df

##### changes the colours of the lines based on thickness, so thicker lines get bolder colours.
def change_colors(processed_lines):
        # Gather distinct thicknesses
    thicknesses = sorted({
        int(line.split('thickness=')[1].split('p')[0])
        for line in processed_lines
        if 'thickness=' in line
    })

    # Shades of orange
    orange_colours = [
        'vvlorange', 'vlorange', 'lorange', 'orange',
        'dorange', 'vdorange', 'vvdorange'
    ]
    # If we exceed the length of orange_colours, we stay at the last color
    darkest_orange = orange_colours[-1]

    thickness_to_colour = {}
    for i, t in enumerate(thicknesses):
        if i < len(orange_colours):
            thickness_to_colour[t] = orange_colours[i]
        else:
            thickness_to_colour[t] = darkest_orange

    coloured_lines = []
    for line in processed_lines:
        if 'thickness=' in line:
            t = int(line.split('thickness=')[1].split('p')[0])
            color = thickness_to_colour.get(t, 'black')
            new_line = line.replace('color=black', f'color={color}')
            coloured_lines.append(new_line)
        else:
            coloured_lines.append(line)

    return coloured_lines

##### does same as above, but for secondary file.
def change_colors_secondary(processed_lines):
    """
    For SECONDARY data: assign shades of GREEN based on thickness.
    Once we run out of green list, stay on the darkest green.
    """
    # Gather distinct thicknesses
    thicknesses = sorted({
        int(line.split('thickness=')[1].split('p')[0])
        for line in processed_lines
        if 'thickness=' in line
    })

    # Shades of green
    green_colours = [
        'vvlred', 'vlred', 'lred', 'red',
        'dred', 'vdred', 'vvdred'
    ]
    # If we exceed the length of green_colours, we stay at the last color
    darkest_green = green_colours[-1]

    thickness_to_colour = {}
    for i, t in enumerate(thicknesses):
        if i < len(green_colours):
            thickness_to_colour[t] = green_colours[i]
        else:
            thickness_to_colour[t] = darkest_green

    coloured_lines = []
    for line in processed_lines:
        if 'thickness=' in line:
            t = int(line.split('thickness=')[1].split('p')[0])
            color = thickness_to_colour.get(t, 'black')
            new_line = line.replace('color=black', f'color={color}')
            coloured_lines.append(new_line)
        else:
            coloured_lines.append(line)

    return coloured_lines

##### changes thickness of lines depending on the frequency, deleting rest.
def process_lines(combination_count, lines, filter_strength_int):
    """Process the lines, adjusting thickness and removing low-frequency pairs."""
    processed_lines = []
    processed_pairs = set()

    for pair, line in lines:
        if pair in processed_pairs:
            # Already processed this pair, skip
            continue

        # Count how many times this pair appears
        count = combination_count[pair]
        if count <= filter_strength_int:
            # Skip pairs that don't meet the frequency threshold
            continue

        processed_pairs.add(pair)

        # Remove any existing thickness specification
        line = re.sub(r',thickness=\d+p', '', line)

        # Append the new thickness (count * 4)
        thickness = f",thickness={count * 4}p"
        modified_line = line.strip() + thickness + "\n"
        processed_lines.append(modified_line)

    return processed_lines


##### writing the first 'blank' lines with 100 0 100 0 format for links
def write_initial_file(pairs, line_format):
    lines = []
    for pair in pairs:
        line = f"{pair[0]} 100 0 {pair[1]} 100 0 color={line_format}"
        lines.append(line)
    return lines

##### counts the number of instances of each line.
def count_combinations(lines: List[str]) -> Tuple[Dict[tuple, int], List[tuple]]:
    combination_count = defaultdict(int)
    line_tuples = []

    for line in lines:
        # Extract only the relevant pair, ignoring extra numeric parts
        parts = line.split()
        if len(parts) >= 4:
            genus1 = parts[0]
            genus2 = parts[3]
            pair = (genus1, genus2)
            combination_count[pair] += 1
            line_tuples.append((pair, line))

    return combination_count, line_tuples

##### we need to process the input files! we need to put it into 0 100 format, do thickness, colour, etc.
def process_data_frame(file_name: str, line_format: str, df: pd.DataFrame) -> bytes:
    # Step 1: get all combination lines
    pairs = generate_combinations_from_df(df)
    # Step 2: do the initial file write logic (in memory or disk)
    write_initial_file(file_name, pairs, line_format)  
    # Step 3: read lines, count combos, filter, etc.
    combination_count, lines = count_combinations(file_name)
    processed_lines = process_lines(combination_count, lines)
    coloured_lines = change_colors(processed_lines)
    # Step 4: use 'write_final_file' to return bytes
    final_data = write_final_file(file_name, coloured_lines)
    return final_data


##### this function starts us on our journey to generate all the links.
def generate_links(primary_data_file: pd.DataFrame,
                   txt_file_name: str,
                   values: dict,
                   line_format: str) -> dict:


    st.write("**Debug**: Entering 'generate_links_in_memory' function")
    st.write(f"txt_file_name = {txt_file_name}")
    st.write("Preview of primary_data_file (first 5 rows):")
    st.write(primary_data_file.head())

    # 1. Map filter strength strings to a numeric cutoff
    filter_strength_map = {
        'Doubles': 1,
        'Triples': 2,
        'Quadruples': 3,
        'Quintuples': 4
    }
    filter_strength_int = filter_strength_map.get(values['filter_strength'], 0)

    # Dictionary of {filename: bytes} to return
    result_files = {}

    ###############################################################################
    # (A) Unfiltered Test Content
    unfiltered_name = f"{txt_file_name}_top_{values['number_to_include']}_unfilteredtest.txt"
    unfiltered_buffer = io.StringIO()

    if values['gen_or_met'] == 'Genera':
        st.write("**Debug**: 'Genera' path for unfiltered test")
        # Build a dictionary: metabolite -> list of genera
        resistances = {}
        for metabolite in primary_data_file.columns[1:]:
            resistant_genera = primary_data_file[primary_data_file[metabolite] != 0]['genus'].tolist()
            if len(resistant_genera) > 0:
                resistances[metabolite] = resistant_genera

        # Show them in the unfiltered “test” file
        for k, v in resistances.items():
            unfiltered_buffer.write(f"{k}: {v}\n")

    else:  # Metabolites
        st.write("**Debug**: 'Metabolites' path for unfiltered test")
        # Remove spaces
        primary_data_file.columns = primary_data_file.columns.str.replace(' ', '')

        resistances = {}
        for genus in primary_data_file.index:
            resistant_cols = primary_data_file.columns[primary_data_file.loc[genus] != 0].tolist()
            if resistant_cols:
                resistances[genus] = resistant_cols

        for k, v in resistances.items():
            unfiltered_buffer.write(f"{k}: {v}\n")

    # Store unfiltered test
    result_files[unfiltered_name] = unfiltered_buffer.getvalue().encode("utf-8")

    ###############################################################################
    # (B) Build the "all outputs" lines (color=black, no thickness yet)
    alloutputs_name = f"{txt_file_name}_top_{values['number_to_include']}_alloutputs.txt"
    alloutputs_buffer = io.StringIO()

    # We'll gather these lines in a list too for frequency counting
    alloutputs_list = []

    if values['gen_or_met'] == 'Genera':
        # Rebuild 'resistances_df' to identify combos
        resistances_df = {}
        for metabolite in primary_data_file.columns[1:]:
            gen_list = primary_data_file[primary_data_file[metabolite] != 0]['genus'].tolist()
            if len(gen_list) > 1:
                resistances_df[metabolite] = gen_list

        # For each metabolite, print combos
        for metabolite, gen_list in resistances_df.items():
            pairs = combinations(gen_list, 2)
            for pair in pairs:
                line_str = f"{pair[0]} 100 0 {pair[1]} 100 0 color=black\n"
                # Write to buffer
                alloutputs_buffer.write(line_str)
                alloutputs_list.append(line_str.strip())

    else:  # Metabolites
        resistances_dict = {}
        for genus in primary_data_file.index:
            meta_list = primary_data_file.columns[primary_data_file.loc[genus] != 0].tolist()
            if len(meta_list) > 1:
                resistances_dict[genus] = meta_list

        for genus, meta_list in resistances_dict.items():
            pairs = combinations(meta_list, 2)
            for pair in pairs:
                line_str = f"{pair[0]} 100 0 {pair[1]} 100 0 color=black\n"
                alloutputs_buffer.write(line_str)
                alloutputs_list.append(line_str.strip())

    # Store "all outputs" as-is
    alloutputs_data = alloutputs_buffer.getvalue().encode("utf-8")
    result_files[alloutputs_name] = alloutputs_data

    ###############################################################################
    # (C) Build Final File with Frequency-based thickness & coloring
    # 1) Count combos
    combination_count, line_tuples = count_combinations(alloutputs_list)

    # 2) Process lines: set thickness = count * 4, remove lines w/ freq <= filter_strength_int
    processed_lines = process_lines(combination_count, line_tuples, filter_strength_int)

    # 3) Color lines based on thickness
    colored_lines = change_colors(processed_lines)

    # 4) Save final lines
    final_name = f"{txt_file_name}_top_{values['number_to_include']}{values['filter_strength']}.txt"
    final_buffer = io.StringIO()
    for line in colored_lines:
        final_buffer.write(line)

    result_files[final_name] = final_buffer.getvalue().encode("utf-8")

    st.write("**Debug**: 'generate_links_in_memory' completed.")
    return result_files

##### for our 'secondary' files, if selected.
def generate_secondary_files(
    values: dict,
    txt_file_name: str,
    edges_list: list,
    line_format: str,
    alloutputs_content: bytes
) -> dict:
    result_files = {}

    secondary_df_map = {
        'Resistance + Nostrains': df_pos_res_nostrain,
        'Resistance + Onlystrains': df_pos_res_yesstrain,
        'Sensitivity + Nostrains': df_pos_sen_nostrain,
        'Sensitivity + Onlystrains': df_pos_sen_yesstrain
    }

    # We'll iterate over each secondary dataset the user selected
    for secondary_label in values.get('res_or_sens_checkboxes', []):
        st.write(f"**Debug**: Secondary dataset = {secondary_label}")

        if secondary_label not in secondary_df_map:
            st.write(f"**Error**: Unrecognized secondary label: {secondary_label}")
            continue

        # Load the DataFrame
        secondary_df = secondary_df_map[secondary_label]
        # Possibly crop to the same edges as the primary
        overlay_df = crop_to_edges(values, secondary_df, edges_list)

        # Build lines for overlay in memory, with color=black initially
        overlay_lines = []
        if values['gen_or_met'] == 'Genera':
            # Similar approach as in generate_links
            combos_map = {}
            for metabolite in overlay_df.columns[1:]:
                # e.g. a row has genus in col[0], 1/0 in others
                gen_list = overlay_df[overlay_df[metabolite] != 0].iloc[:, 0].tolist()
                if len(gen_list) > 1:
                    combos_map[metabolite] = gen_list

            for metabolite, g_list in combos_map.items():
                pairs = combinations(g_list, 2)
                for pair in pairs:
                    line_str = f"{pair[0]} 100 0 {pair[1]} 100 0 color=black"
                    overlay_lines.append(line_str)

        else:  # Metabolites
            # Remove spaces, set index
            overlay_df.columns = overlay_df.columns.str.replace(' ', '')
            # Typically col[0] is 'genus'
            if 'genus' in overlay_df.columns:
                overlay_df.set_index('genus', inplace=True, drop=False)

            combos_map = {}
            for genus in overlay_df.index:
                meta_list = overlay_df.columns[overlay_df.loc[genus] != 0].tolist()
                if len(meta_list) > 1:
                    combos_map[genus] = meta_list

            for genus, m_list in combos_map.items():
                pairs = combinations(m_list, 2)
                for pair in pairs:
                    line_str = f"{pair[0]} 100 0 {pair[1]} 100 0 color=black"
                    overlay_lines.append(line_str)

        # Now for each filter strength the user selected, build a final file
        filter_str_map = {
            'Doubles': 1,
            'Triples': 2,
            'Quadruples': 3,
            'Quintuples': 4
        }

        for strength_label in values.get('filter_strength_checkboxes', []):
            st.write(f"**Debug**: Secondary filter strength = {strength_label}")
            filter_strength_int = filter_str_map.get(strength_label, 0)

            # 1) Count combos
            combination_count, line_tuples = count_combinations(overlay_lines)
            # 2) process lines (remove freq <= filter_strength_int, set thickness)
            processed_lines = process_lines(combination_count, line_tuples, filter_strength_int)
            # 3) color them green
            colored_lines = change_colors_secondary(processed_lines)

            # Build a final filename, e.g. "top_genera_res_1strain_Resistance + Onlystrains_Triples.txt"
            # or something similar
            safe_secondary_label = secondary_label.replace(' ', '_')
            out_filename = f"{txt_file_name}_{safe_secondary_label}_{strength_label}.txt"
            final_buf = io.StringIO()
            for line in colored_lines:
                final_buf.write(line + "\n" if not line.endswith("\n") else line)

            result_files[out_filename] = final_buf.getvalue().encode("utf-8")

    return result_files

##### basically the start of our journey, to generate the chromosome file.
def generate_edges_file(primary_data_file, txt_file_name, values):
    edges_list = []
    # We'll create a StringIO for the chromosome (edges) data
    edges_buffer = io.StringIO()

    # Build the output filename using 'number_to_include'
    output_file1 = f"{txt_file_name}_top_{values['number_to_include']}_chromosome_file.txt"

    # We'll store all resulting files in this dictionary
    result_files = {}

    # Rotate edge colors if needed
    colours = [
        'pred', 'orange', 'yellow', 'green', 'blue', 'purple', 'grey',
        'pred', 'porange', 'pyellow', 'pgreen', 'pblue', 'ppurple'
    ]

    # Decide the line_format for primary data
    if values['res_or_sens'] in ['Resistance + Nostrains', 'Resistance + Onlystrains']:
        line_format = 'black'
    else:
        line_format = 'black'

    # 1) If Genera is chosen
    if values['gen_or_met'] == 'Genera':
        st.write("Genera is chosen, code should work")

        # Collect the Genera from the first column
        for index, row in primary_data_file.iterrows():
            genera = row[0]
            edges_list.append(genera)

        st.write("The genera list is:")
        st.write(edges_list)

        # Write each genus as a 'chr' line in memory
        for i, genus in enumerate(edges_list):
            colour = colours[i % len(colours)]
            edges_buffer.write(f"chr - {genus} {genus} 0 100 {colour}\n")

    # 2) If Metabolites is chosen
    elif values['gen_or_met'] == 'Metabolites':
        st.write("Metabolites is chosen, code may work")

        # Remove spaces from column names, because Circos can choke on those
        primary_data_file.columns = primary_data_file.columns.str.replace(' ', '')

        # The first column is 'genus'? The rest are metabolite columns
        edges_list = primary_data_file.columns.tolist()[1:]  # skip index 0
        st.write("Printing primary_data_file:")
        st.write(primary_data_file)

        # Write each metabolite as a 'chr' line in memory
        for i, metabolite in enumerate(edges_list):
            colour = colours[i % len(colours)]
            edges_buffer.write(f"chr - {metabolite} {metabolite} 0 100 {colour}\n")

        # Assuming the first column is 'genus', set that as index
        primary_data_file.set_index('genus', inplace=True)

    # Convert the chromosome text to bytes and store in our results
    edges_file_bytes = edges_buffer.getvalue().encode("utf-8")
    result_files[output_file1] = edges_file_bytes

    # Now generate link file(s) in memory
    links_files = generate_links(primary_data_file, txt_file_name, values, line_format)
    result_files.update(links_files)

    # Generate secondary files, if any, in memory
    #   The 'alloutputs' file is presumably in links_files.
    #   We'll fetch it by name:
    alloutputs_filename = f"{txt_file_name}_top_{values['number_to_include']}_alloutputs.txt"
    alloutputs_content = links_files.get(alloutputs_filename, b"")

    secondary_files = generate_secondary_files(values, txt_file_name, edges_list, line_format, alloutputs_content)
    result_files.update(secondary_files)

    return result_files

##### crops the dataframe properly, because we just want the total resistances
def crop_df(values, df):
    # we just want to sum up the total resistances in the genera, and if user has selected metabolites, for that.
    if values['gen_or_met'] == 'Genera':
        df['Total Resistances'] = df.iloc[:, 1:].astype(int).sum(axis=1)
        top_genera = df.sort_values(by='Total Resistances', ascending=False).head(values['number_to_include'])
        return top_genera
    else:
        metabolite_resistances = df.iloc[:, 1:].astype(int).sum(axis=0)
        top_metabolites = metabolite_resistances.sort_values(ascending=False).head(values['number_to_include'])
        columns_to_keep = [df.columns[0]] + top_metabolites.index.tolist()
        filtered_df = df[columns_to_keep]
        return filtered_df
    
##### leads us to the generate_edges_file, which does the heavy lifting 
def generate_files(values):
    # based on dropdown selection, first creating the map, selecting it
    res_or_sens_map = {
    'Resistance + Nostrains': ('res_1strain', df_pos_res_nostrain),
    'Resistance + Onlystrains': ('res_0strain',  df_pos_res_yesstrain),
    'Sensitivity + Nostrains':  ('sen_1strain',  df_pos_sen_nostrain),
    'Sensitivity + Onlystrains': ('sen_0strain', df_pos_sen_yesstrain)
    }

    # determining file name for edge_file
    # first creating the map
    gen_or_met_map = {
        'Genera': 'top_genera',
        'Metabolites': 'top_metabolites'
    }

    # get labels and choose df from the maps
    short_label, df = res_or_sens_map[values['res_or_sens']]
    prefix = gen_or_met_map[values['gen_or_met']]
    # getting the txt file name 
    txt_file_name = f"{prefix}_{short_label}"

    ########## cropping the df
    primary_cropped_df = crop_df(values, df)

    # 3) Generate all in-memory files (edges + links + secondaries) --> the edges file function leads to the rest
    result_files = generate_edges_file(primary_cropped_df, txt_file_name, values)

    # 4) Return the dictionary so display() can do st.download_button
    return result_files

#### main function
def display():
    st.header("Circos Settings")
    values = {}
    values['res_or_sens'] = st.selectbox("Primary Dataset:", 
        ['Resistance + Nostrains', 'Resistance + Onlystrains', 'Sensitivity + Nostrains', 'Sensitivity + Onlystrains'], index=0)
    values['gen_or_met'] = st.selectbox("Sort by Metabolite or Genus:", 
        ['Metabolites', 'Genera'], index=1)
    # values['choose_colours'] = st.selectbox("Selection for Dynamic Line Colour:", 
    #     ['Red', 'Orange', 'Green', 'Blue', 'Purple'], index=4)
    values['filter_strength'] = st.selectbox("Selection for Pairs:", 
        ['Doubles', 'Triples', 'Quadruples', 'Quintuples'], index=1)
    values['number_to_include'] = st.slider("Number of Edges", min_value=1, max_value=50, step=1, value=50)
    values['res_or_sens_checkboxes'] = st.multiselect("Secondary generated files (Res/Sens, Nostrains/Onlystrains):", 
        ['Resistance + Nostrains', 'Resistance + Onlystrains', 'Sensitivity + Nostrains', 'Sensitivity + Onlystrains'])
    values['filter_strength_checkboxes'] = st.multiselect("Secondary generated files (Selection Options):", 
        ['Doubles', 'Triples', 'Quadruples', 'Quintuples'])
    # New widget for secondary dataset naming:
    values['secondary_dataset'] = st.selectbox("Secondary Dataset:", ['Nostrains', 'Onlystrains'], index=0)
    
    if st.button("Generate Files"):
        result_files = generate_files(values)
        st.success("Files generated!")
        st.subheader("Download Generated Files")

        for file_name, file_data in result_files.items():
            st.download_button(
                label=f"Download {file_name}",
                data=file_data,
                file_name=file_name,
                mime="text/plain"
            )