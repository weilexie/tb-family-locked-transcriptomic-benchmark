# Family-locked TB blood-transcriptome benchmark

Analysis code for the Intelligent Healthcare manuscript: cohort-family-locked validation of public blood-transcriptomic classifiers for active tuberculosis.

Public expression matrices come from `curatedTBData` 2.4.0 (`assay_curated`) via ExperimentHub. Source GEO accessions are listed in Supplementary Table S1.

## Frozen split and seed

- Primary random seed: `20260501` (`01_config/config.yml`, `01_config/external_validation_lockfile.json`)
- Family lock date: 2026-05-01
- Family-grouped 5-fold CV and model fitting use `20260501`
- Combination-wise locked-external evaluation uses `20260501 + i`
- Permutation negative controls in script 13 use `20260502`

Do not overwrite the lockfile (`overwrite_external_validation_lockfile: false`).

## Run

From the repository root, with R 4.5.x:

```r
Rscript run_all.R --mode prepare_data
Rscript run_all.R --mode optimize_final
Rscript run_all.R --mode enhance_final
Rscript run_all.R --mode stress_test
Rscript run_all.R --mode manuscript_final
```

`prepare_data` needs internet for ExperimentHub. Derived matrices are written to `03_data_curated/` and `07_results/`, which are gitignored.

Package versions are recorded in `renv.lock`. The pipeline currently runs with `runtime.use_renv: false`; restore from the lockfile if you need the same library set.

## Layout

Folders are numbered in run order. `03_data_curated/` is created by `prepare_data` and is not stored in this repository.

- `01_config/`: phenotype rules, family lock, seed
- `02_data_meta/`: harmonized sample and family tables
- `04_R_functions/`: helpers `00`–`30` (figures `21`–`30`)
- `05_pipeline/`: scripts in order
  - `00`–`07` prepare data and score reference signatures
  - `08_model_optimization.R` calls `08a` development screen, `08b` locked external validation, `08c` summary
  - `09`–`13` verification, ensembles, operating points, stress tests
  - `14_refresh_manuscript_assets.R` calls `14a_manuscript_package_updates.R`, then `15`
  - `16` sample-flow check, `17` statistics, `18` published formulas, `19` 10,000-permutation test
- `06_supplementary_tables/`: S1–S8. S3 includes `published_formula_same_lock`. S8 has 16 selective-prediction states (primary elastic-net Sens>=0.90 was infeasible in development).

Nested cross-validation is not in this deposit. The second model role is penalized-linear (primary elastic-net, expanded ridge).
