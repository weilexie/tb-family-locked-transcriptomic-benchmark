# Sync/verify current final ML optimization outputs.
# Current manuscript-facing analysis uses results/ml_final.
dst_dir <- project_path("results/ml_final")
if (!dir.exists(dst_dir)) stop("Missing final optimization outputs: results/ml_final", call. = FALSE)
required <- c(
  "algorithm_combination_grid.tsv",
  "development_combination_performance.tsv",
  "selected_combinations.tsv",
  "external_pooled_summary.tsv",
  "external_predictions.tsv.gz",
  "external_lockfile_integrity.tsv"
)
missing <- required[!file.exists(file.path(dst_dir, required))]
if (length(missing) > 0) stop(sprintf("Missing final optimization files: %s", paste(missing, collapse = ", ")), call. = FALSE)
message("Final ML outputs verified under results/ml_final/.")
