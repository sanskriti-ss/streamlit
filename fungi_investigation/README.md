# Fungi phenotype data pipelines

Scripts to build experimental and predicted fungi layers for the symbiosis network app.

## Experimental — FUNG-GROWTH (manual)

Carbon-source utilization for ~398 fungal species: [FUNG-GROWTH](https://www.fung-growth.org/)

FUNG-GROWTH has no public API — search and export records manually on the website, then place a CSV in `fungi_data/raw/fung_growth/`.

1. Export or download a CSV into `fungi_data/raw/fung_growth/`
2. Run ingest:

```bash
python -m fungi_investigation.fung_growth_ingest \
  --input fungi_data/raw/fung_growth/export.csv
```

A sample file is provided for testing:

```bash
python -m fungi_investigation.fung_growth_ingest \
  --input fungi_data/raw/fung_growth/sample_export.csv
```

Output: `fungi_data/experimental/fungi_phenotypes_long.csv`

## Predicted — production + utilization (automated)

Download reference genomes, run antiSMASH, and scan protein FASTA for utilization:

```bash
/Users/sanskriti/miniconda3/bin/python -m fungi_investigation.fetch_fungi_predicted
```

Species list: `fungi_investigation/selected_fungi.yaml`  
Genomes: `data/genomes_fungi/`  
antiSMASH output: `results/antismash_fungi/`

## Predicted — production (antiSMASH only)

**Option A — local antiSMASH** (requires `antismash_env` conda env on PATH):

```bash
PATH="$HOME/miniconda3/envs/antismash_env/bin:$PATH" \
  python -m fungi_investigation.fetch_fungi_predicted --skip-download --skip-utilization
```

**Option B — antiSMASH-DB bulk** (421+ fungal genomes):

Download JSON subset from [antiSMASH-DB v5](https://dl.secondarymetabolites.org/database/5.0/), then:

```bash
python -m fungi_investigation.antismash_db_ingest \
  --json-root /path/to/fungal_antismash_jsons
```

Output: `fungi_data/predicted/fungi_phenotype_confidence.csv`

## Predicted — utilization (CAZyme keywords)

Scans protein FASTA under `data/genomes/` using the same keyword map as the bacteria genome pipeline:

```bash
python -m fungi_investigation.fungi_utilization_predict \
  --fasta-root data/genomes
```

Output: `fungi_data/predicted/fungi_utilization_confidence.csv`

For MycoCosm annotations, download protein FASTA via the [JGI Data Portal API](https://sites.google.com/lbl.gov/data-portal-help/home/tips_tutorials/api-tutorial) into a directory and pass `--fasta-root`.

## Metabolite name harmonization

Cross-kingdom edges require aligned metabolite names. Edit:

`fungi_data/metabolite_aliases.csv`

## Production ↔ utilization ontology bridge

antiSMASH *production* labels (BGC classes: siderophores, terpenes, polyketides) and
carbon-source *utilization* labels (glucose, cellulose, …) share no vocabulary, so raw
synergy is always zero. The loader bridges them via **extracellular degradation
cross-feeding**: a species that utilizes/degrades a polymer releases monomers into the
shared environment (inferred production). The polymer→monomer map lives in:

`fungi_data/metabolite_bridges.csv`

Columns: `polymer,monomer,release_confidence` (e.g. `cellulose,glucose,0.7`). Derived
production rows carry a `+degradation` suffix on their source layer and confidence equal
to the utilization confidence × `release_confidence`. Add rows to extend the bridge.

## Symbiosis app

```bash
streamlit run symbiosis_network_app.py
```

Paths are configured in `symbiosis_data_paths.yaml` at the repo root.
