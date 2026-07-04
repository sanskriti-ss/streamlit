# FUNG-GROWTH raw exports

Place bulk exports from [FUNG-GROWTH](https://www.fung-growth.org/) here.

## Obtaining data

1. Browse https://www.fung-growth.org/ and export growth profiles, or
2. Contact the FUNG-GROWTH authors for a bulk CSV dump.

## Expected CSV format

```csv
species,genus,substrate,growth_positive,genome_accession
Aspergillus niger,Aspergillus,D-glucose,1,GCA_000002395.2
Aspergillus niger,Aspergillus,cellulose,1,GCA_000002395.2
```

`growth_positive` may also be `growth_score` (numeric; >0 = growth).

## Ingest

```bash
python -m fungi_investigation.fung_growth_ingest \
  --input fungi_data/raw/fung_growth/export.csv
```

Output: `fungi_data/experimental/fungi_phenotypes_long.csv`

A small **sample** file `sample_export.csv` is included for pipeline testing.
