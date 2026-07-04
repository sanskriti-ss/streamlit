# Genome investigation (optional layer)

Adds genome metadata enrichment, **selected-only** genome download, optional antiSMASH, and species prioritization—without scanning all ~30k BacDive strains by default.

## Phase 1 — Genome metadata enrichment

```bash
python -m genome_investigation.genome_enrichment \
  --input species_data/step3_met_util_exploded.csv.zip \
  --output genome_investigation/results/Step2_5_genome_enriched.csv \
  --limit 100 \
  --allow-species-only
```

| Flag | Description |
|------|-------------|
| `--input` | Step1/Step2/Step3 CSV or zip |
| `--output` | `Step2_5_genome_enriched.csv` |
| `--limit` | Max rows |
| `--force-refresh` | Ignore API cache |
| `--min-confidence` | Threshold for failure log |
| `--allow-species-only` | Accept weaker NCBI species-level matches |

Lookup order: **BacDive API v2** (`/v2/fetch/{BacID}`) → **NCBI Datasets** by species. Responses cached under `genome_investigation/cache/api/`. Failures: `logs/genome_lookup_failures.csv`.

## Phase 1b — Selected genome download

Edit `selected_species.yaml`, then:

```bash
python -m genome_investigation.genome_download \
  --selected-species genome_investigation/selected_species.yaml \
  --genome-enriched genome_investigation/results/Step2_5_genome_enriched.csv \
  --download-genomes \
  --min-confidence 0.7 \
  --max-genomes 20
```

Files go to `data/genomes/{species_slug}/{accession}/` (requires NCBI `datasets` CLI).

## Phase 1c — Optional antiSMASH

See main [README](../README.md#optional-antismash-biosynthetic-gene-cluster-analysis). Disabled by default (`enable_antismash: false` in `config/genome_config.yaml`).

## Phase 1d — Species prioritization

```bash
python -m genome_investigation.species_prioritization \
  --genome-enriched genome_investigation/results/Step2_5_genome_enriched.csv \
  --antismash-summary genome_investigation/results/antismash_summary.csv \
  --output genome_investigation/results/ranked_species_candidates.csv
```

Tune weights in `config/prioritization_weights.yaml`.

## Paper-ready visualizations

After enrichment, generate figures and Table 1 for manuscripts:

```bash
python -m genome_investigation.visualize_results
```

Outputs in `genome_investigation/results/paper/` (300 dpi PNGs, CSV/Markdown table, `RESULTS_INTERPRETATION.md`).  
Also viewable in the Streamlit tab **Genome and BGC Evidence**.

## Automated genomic follow-up (download + antiSMASH + AMRFinder + figures)

One command for the nine focal species in `results/selected_test_input.csv`:

```bash
python -m genome_investigation.genome_followup_pipeline
```

**Prerequisites (optional per step):**

| Tool | Step |
|------|------|
| NCBI `datasets` CLI | Genome download |
| `antismash` | BGC detection (confidence ≥ 1.0 strains) |
| `amrfinder` (AMRFinderPlus) | Resistance genes (*B. cereus*, *M. luteus*) |

**Outputs** (timestamped under `results/genomic_followup/run_YYYYMMDD_HHMMSS/`):

| File | Content |
|------|---------|
| `focal_genome_enriched.csv` | Per-strain accessions |
| `download_manifest.csv` | Download status |
| `antismash_summary.csv` | BGC counts / types |
| `amrfinder_summary.csv` | AMR gene calls |
| `bacdive_phenotype_metabolites.csv` | Util / prod / resistance lists from BacDive API |
| `targeted_gene_hits.csv` | Protein/GFF keyword hits + AMRFinder genes |
| `phenotype_confidence.csv` | Per-(species, activity, metabolite) confidence (BacDive-observed pinned to 1.0/0.0, untested metabolites from genomic evidence) |
| `integrated_genomic_followup.csv` | Merged summary |
| `figures/` | Pipeline status, BGC vs AMR, heatmaps, composite dashboard, and `fig7/8/9_confidence_{resistance,production,utilization}.png` |

**Phenotype confidence** (`phenotype_confidence.py`) estimates how likely each species is resistant to / produces / utilizes metabolites **beyond those tested in BacDive**. Tested phenotypes are pinned (resistant/positive → 1.0, sensitive/negative → 0.0, marked with ● in figures); untested metabolites are scored from genomic evidence: AMRFinder gene `Class`/`Subclass` × `%identity`×`%coverage` (resistance), antiSMASH BGC product types (production, capped at 0.65 since BGC presence is putative), and CAZyme/transporter annotation counts in the proteome (utilization, saturating to 0.9). These are heuristic genomic-evidence scores, not calibrated probabilities.

Stable copies also written to `results/genomic_followup/` and `results/antismash_summary.csv` / `amrfinder_summary.csv`.

Useful flags:

```bash
# Re-run figures only (after tools finished)
python -m genome_investigation.genome_followup_pipeline \
  --skip-enrichment --skip-download --skip-antismash --skip-amrfinder --skip-phenotype --skip-gene-search

# Phenotype + enrichment only (no external tools)
python -m genome_investigation.genome_followup_pipeline --skip-download --skip-antismash --skip-amrfinder --skip-gene-search
```

## Tests

```bash
pytest genome_investigation/tests -q
```
