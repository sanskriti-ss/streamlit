# Species Analysis Features

## Known Issues and Notes

total_metabolites is still len(metabolite_cols) (same for all species). That's ok for reference but not used for rate -- rates use metabolites_tested (good). However get_top_species_across_files uses total_metabolites in overall_utilization_rate; you may want to use aggregated metabolites_tested instead.
Performance: the function iterates with iterrows over large dataframes (34k x ~1.5k cols) -- expect slow execution; vectorized implementation would be much faster.

---
# Adding Features
## Filtered Genus-Metabolite Selection for Sankey Diagrams
Users can optionally filter the genus and metabolite dropdowns to only show combinations that have actual test data.

### How it works:
1. **Default (unchecked):** All genera and metabolites are shown in the dropdowns, even if they have no intersection. Selecting a combination with no data will show a warning.

2. **Filtered mode (checked):** 
   - When the user checks "Only show non-zero genus-metabolite combinations", the app loads all 4 activity files (production, utilization, resistance, sensitivity).
   - A precomputed index maps each genus to the set of metabolites for which any strain has nonzero test data.
   - The dropdowns are filtered so only valid combinations are shown.
   - If a genus is selected first, only metabolites with data for that genus appear.
   - If a metabolite is selected first, only genera with data for that metabolite appear.

#### Backend:

- The activity files are loaded once and cached using `@st.cache_data`.
- The genus-metabolite index is built using vectorized pandas operations (groupby + any), making it fast even for large datasets.
- Subsequent interactions use the cached data, so filtering is instant after the initial load.
