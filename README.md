### Intro
BacDive is the largest aggregate database for bacterial information. From here, we’ve built a visualization tool that allows the user to see resistance, sensitivity, production, and utilization trends between strains, species, genera, and more. The Bacdive Streamlit can be used to generate a great number of graphs and tables. The user by default gets information on metabolite resistances, production, utilization, and sensation. The user can further select whether they want to view strains of bacteria species, or exclude them. Furthermore, instead of only showing positive results, negative results can also be chosen. It can detail the shared antibiotic resistances, production, utilization, and sensitivity across different genera; show the relative ranking of total unique metabolites of a genus in aforementioned categories; give a general overview of the BacDive database stats, and much more. 

For more information about the data available from BacDive, visit the [BacDive Overview](BacDive_Overview.md).

# How to use the streamlit :)
The webapp is split into six tabs: General Overview, Circos, Trends, Cards, By the Numbers, and Comparison.
They each produce different types of visualizations, tables, and/or graphs.

## General Overview
The general overview allows you to do selections for negatively/positively tested, strain/isolate, production/utilization/resistance/sensitivity

You can then select genera by highest absolute values or proportions: that is, genera that have the highest number of species that are ‘resistant’ to the metabolite, or genera that have the highest proportion of species that are. However, we advise the user to take this section with a grain of salt, as it is just meant to provide a quick overview.

## Circos
Users are able to choose a ‘primary’ category of interest — whether the edges be metabolites or genera. They can also choose a ‘secondary’ categor(ies) of interest in order to overlap multiple graphs on a single one. The output will be in .txt files, that you will need to run on circos yourself. You must install it for this to work! Instructions available in [Circos Instructions](Circos_Instructions.md).

## Trends
Two possible visualizations on this page:
1) A parallel diagram. Another way of looking at the number of species in each category (prod/res/sen/util), across genera (that you can choose). 

2) A sankey diagram: select a genus, and you'll see which metabolites it has been tested for, in which category, and whether it was tested positively or negatively.
    
## Cards
Default load 'playing cards' for each of the genera, where you can see how many metabolites the isolate/strain species in the genus have been positively tested for production/utilization/resistance/sensitivity.

## By the Numbers
There are three main features on this page. First is the summary statistics section at the top, which details how many metabolites, genera, isolates, and strains there are in your dataset. Following that, there is a “Metabolite Counts Bar Graph” section, where you can plot the top 15 genera based on neg/pos, isolate/strain, and prod/util/res/sen. This gives an overview of which genera have the most data in that category. See the figure below for an example. The third main feature is the “Homogeneous Metabolite Summary by Genus” section. For a selected metabolite category, this section lists genera where, for at least one metabolite, all species (with a minimum of 5 species) tested uniformly positive or uniformly negative. It also shows which metabolite(s) met that criteria. If a genus appears with both positive and negative results (i.e. mixed), it is omitted. 

## Comparison
This is the section where interesting synergies between genera can be found. 
The user can select up to 10 genera to compare. Upon selecting to include/exclude strains, and the top genera, a table is generated with each genus and its species count as the row name, with columns for which metabolites each genera tested positive for in resistance, sensitivity, production, and utilization. 

More important is the generated shared and synergy summary, which shows the prod/util synergy between the genera, i.e., do any of the genera have species that produce metabolites that the others utilize? Do they have shared utilization, which could indicate that certain metabolites are promoting their growth? Is there any metabolite they are all resistant to?
 
(Note that it is very slow to run, currently. Working on it!)


## Genome investigation (optional)

An optional pipeline under `genome_investigation/` enriches strains with genome accessions (BacDive → NCBI), supports **selected-only** genome download, optional antiSMASH, and ranked follow-up candidates. It does **not** modify Step1–Step3 CSV generation or existing tabs unless you run the new tools.

- Streamlit tab: **Genome and BGC Evidence**
- Details: [genome_investigation/README.md](genome_investigation/README.md)

### Optional antiSMASH biosynthetic gene cluster analysis

antiSMASH is **off by default** (`enable_antismash: false` in `genome_investigation/config/genome_config.yaml`). It runs only for species listed in `genome_investigation/selected_species.yaml` after genomes are downloaded to `data/genomes/`.

1. Install antiSMASH via conda (recommended for Apple Silicon):
   ```bash
   conda create -n antismash_env -c bioconda -c conda-forge antismash
   conda activate antismash_env
   pip install nrpys  # builds ARM-native Rust extension
   download-antismash-databases
   ```
   The pipeline auto-discovers the binary at `~/miniconda3/envs/antismash_env/bin/antismash`.
2. Download selected genomes: `python -m genome_investigation.genome_download --download-genomes ...`
3. Run or import results:

```bash
# Run locally (requires antiSMASH on PATH)
python -m genome_investigation.antismash_runner --run-antismash

# Or parse results produced elsewhere
python -m genome_investigation.antismash_runner --parse-antismash-only \
  --antismash-output-dir results/antismash
```

Outputs: `results/antismash/{species_slug}/{accession}/` and `genome_investigation/results/antismash_summary.csv`. Prefer GenBank/GBFF input; FASTA is allowed but logged as lower-quality annotation input. Extra antiSMASH flags are **not** passed unless you supply `--antismash-extra-args`.

### Automated genomic follow-up (one command)

Runs enrichment, NCBI download, antiSMASH (confidence ≥ 1.0), AMRFinder (*B. cereus* / *M. luteus*), BacDive metabolite export, gene search, and figures:

```bash
python -m genome_investigation.genome_followup_pipeline
```

Requires optional CLIs: NCBI `datasets`, `antismash`, `amrfinder`. Results: `genome_investigation/results/genomic_followup/` (see [genome_investigation/README.md](genome_investigation/README.md)).

## Symbiosis network (bacteria + fungi)

A separate Streamlit app explores cross-kingdom production ↔ utilization synergy networks. It works with partial data (bacteria-only is fine).

```bash
streamlit run symbiosis_network_app.py
```

- Configure data paths in [`symbiosis_data_paths.yaml`](symbiosis_data_paths.yaml)
- Fungi pipelines: [`fungi_investigation/README.md`](fungi_investigation/README.md)

Layers (multi-select in the app):

| Layer | Source |
|-------|--------|
| Experimental bacteria | BacDive Step3 (`species_data/step3_met_*`) |
| Predicted bacteria | `genome_investigation/results/genomic_followup/phenotype_confidence.csv` |
| Experimental fungi | FUNG-GROWTH → `fungi_data/experimental/fungi_phenotypes_long.csv` |
| Predicted fungi | antiSMASH + CAZyme → `fungi_data/predicted/` |

# If you want to make local changes
make sure requirements are installed :)))
to run, do streamlit run app.py

# Notes to devs:
Notes:
Species Analysis is a little slow to load. Still figuring out why?
