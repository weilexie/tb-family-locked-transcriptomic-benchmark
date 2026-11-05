# Final manuscript-facing package sync.
# This mode consumes existing final results and does not train or select models.

if (!exists("project_path", mode = "function")) {
  Sys.setenv(TBML_PROJECT_ROOT = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  source(file.path("04_R_functions", "00_load_all.R"))
}

required <- c(
  "data_meta/sample_harmonization_master.tsv",
  "config/external_validation_lockfile.json",
  "results/ml_final/external_pooled_summary.tsv",
  "results/ml_final/external_predictions.tsv.gz",
  "results/ml_final/paired_family_bootstrap_model_differences.tsv",
  "results/ml_final/calibrated_operating_points_external.tsv",
  "results/ml_final/selective_prediction_error_asymmetry.tsv",
  "results/ml_final/selective_prediction_family_level_coverage.tsv"
)
missing <- required[!file.exists(project_path(required))]
if (length(missing) > 0) {
  stop(sprintf("Cannot refresh manuscript-facing package; missing required outputs: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

Sys.setenv(TBML_PACKAGE_UPDATE_STAGE = "pre_figure")
source(project_path("scripts/14a_manuscript_package_updates.R"), local = TRUE)
source(project_path("scripts/15_integrate_validation_checks.R"), local = TRUE)

sf_build_submission_figures()
message("Submission figure rebuild completed under 09_submission/.")

sync_dir <- function(src, dst) {
  src_dir <- project_path(src)
  dst_dir <- project_path(dst)
  if (!dir.exists(src_dir)) stop(sprintf("Missing source directory: %s", src), call. = FALSE)
  if (dir.exists(dst_dir)) unlink(dst_dir, recursive = TRUE, force = TRUE)
  dir.create(dirname(dst_dir), recursive = TRUE, showWarnings = FALSE)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(src_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(files) > 0) {
    file.copy(files, dst_dir, recursive = TRUE, overwrite = TRUE)
  }
}

sync_dir("figures_submission", "manuscript/figures")
sync_dir("plot_data_submission", "manuscript/figure_data")
sync_dir("captions_submission", "manuscript/captions")
file.copy(
  project_path("captions_submission/main_figure_captions.md"),
  project_path("manuscript/figure_legends.md"),
  overwrite = TRUE
)

Sys.setenv(TBML_PACKAGE_UPDATE_STAGE = "post_sync")
source(project_path("scripts/14a_manuscript_package_updates.R"), local = TRUE)
Sys.unsetenv("TBML_PACKAGE_UPDATE_STAGE")

message("Final manuscript-facing figures, plot data, and captions refreshed under manuscript/.")
