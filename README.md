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

- `01_config/`: phenotype rules, family lock, seed
- `02_data_meta/`: harmonized sample and family tables
- `04_reusable_R_functions/`: helpers
- `05_pipeline_scripts/`: numbered scripts 00-15
- `revision_extras/`: revision-only mapping, permutation, and original-signature checks
- `supplementary_tables/`: S1-S8 workbooks. S8 has 16 selective-prediction states (one primary elastic-net Sens>=0.90 rule infeasible in development).
