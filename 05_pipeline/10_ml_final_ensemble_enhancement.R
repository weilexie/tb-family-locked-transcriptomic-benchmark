# Development-selected ensemble enhancement for final ML results.
# This is a secondary/sensitivity layer: ensemble selection uses development OOF predictions only.

out_dir <- "results/ml_final"
dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)

grid <- read_tsv(file.path(out_dir, "algorithm_combination_grid.tsv"))
dev_perf <- read_tsv(file.path(out_dir, "development_combination_performance.tsv"))
oof <- read_tsv(file.path(out_dir, "development_oof_predictions.tsv.gz"))
lock_hash_before <- hash_file_or_object("config/external_validation_lockfile.json")

candidate_rows <- list()
for (task in automl_final_tasks()) {
  candidate_rows[[length(candidate_rows) + 1]] <- mlfinal_ensemble_candidate_sets(dev_perf, task, top_n = 12)
}
candidates <- do.call(rbind, candidate_rows)
write_tsv_atomic(candidates, file.path(out_dir, "ensemble_candidate_grid.tsv"))

message(sprintf("Final ensemble candidate grid: %s candidates", nrow(candidates)))

dev_rows <- lapply(seq_len(nrow(candidates)), function(i) {
  if (i %% 4 == 0) message(sprintf("Development ensemble candidate %s/%s", i, nrow(candidates)))
  mlfinal_ensemble_oof_summary(oof, candidates[i, , drop = FALSE])
})
ensemble_dev <- do.call(rbind, dev_rows)
ensemble_dev <- ensemble_dev[order(ensemble_dev$task, ensemble_dev$development_AUROC, ensemble_dev$development_AUPRC, decreasing = TRUE), , drop = FALSE]
write_tsv_atomic(ensemble_dev, file.path(out_dir, "ensemble_development_performance.tsv"))

selected <- do.call(rbind, lapply(split(ensemble_dev, ensemble_dev$task), function(d) {
  d[order(d$development_AUROC, d$development_AUPRC, d$development_median_family_AUROC, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
}))
write_tsv_atomic(selected, file.path(out_dir, "selected_ensemble_development_only.tsv"))

external_prediction_rows <- list()
for (task in automl_final_tasks()) {
  task_candidates <- candidates[candidates$task == task, , drop = FALSE]
  selected_task <- selected[selected$task == task, , drop = FALSE]
  base_ids <- unique(unlist(strsplit(task_candidates$combination_ids[task_candidates$ensemble_id %in% selected_task$ensemble_id], ";", fixed = TRUE)))
  base_ids <- base_ids[nzchar(base_ids)]
  message(sprintf("External base prediction for %s selected ensemble: %s base models", task, length(base_ids)))
  external_base <- mlfinal_external_base_predictions(task, base_ids, grid, lock_hash_before, seed = 20260501)
  write_tsv_atomic(external_base, file.path(out_dir, sprintf("ensemble_external_base_predictions_%s.tsv.gz", automl_final_task_short(task))))
  cand <- candidates[candidates$ensemble_id == selected_task$ensemble_id[1], , drop = FALSE]
  external_prediction_rows[[length(external_prediction_rows) + 1]] <- mlfinal_ensemble_external_predict(oof[oof$task == task, , drop = FALSE], external_base, cand)
}
ensemble_external_pred <- do.call(rbind, external_prediction_rows)
write_tsv_atomic(ensemble_external_pred, file.path(out_dir, "ensemble_external_predictions.tsv.gz"))
ensemble_external_summary <- mlfinal_ensemble_external_summary(ensemble_external_pred)
write_tsv_atomic(ensemble_external_summary, file.path(out_dir, "ensemble_external_pooled_summary.tsv"))

current <- read_tsv(file.path(out_dir, "external_pooled_summary.tsv"))
current_rows <- data.frame(
  task = current$task,
  model_layer = paste("current", current$selected_role, sep = "_"),
  model = current$selector_classifier,
  AUROC = current$AUROC,
  AUPRC = current$AUPRC,
  sensitivity = current$sensitivity,
  specificity = current$specificity,
  Brier = current$Brier,
  calibration_slope = current$calibration_slope,
  stringsAsFactors = FALSE
)
ensemble_rows <- data.frame(
  task = ensemble_external_summary$task,
  model_layer = "development_selected_ensemble",
  model = ensemble_external_summary$ensemble_name,
  AUROC = ensemble_external_summary$AUROC,
  AUPRC = ensemble_external_summary$AUPRC,
  sensitivity = ensemble_external_summary$sensitivity,
  specificity = ensemble_external_summary$specificity,
  Brier = ensemble_external_summary$Brier,
  calibration_slope = ensemble_external_summary$calibration_slope,
  stringsAsFactors = FALSE
)
comparison <- rbind(current_rows, ensemble_rows)
write_tsv_atomic(comparison, file.path(out_dir, "current_vs_ensemble_external_summary.tsv"))

lock_hash_after <- hash_file_or_object("config/external_validation_lockfile.json")
if (!identical(lock_hash_before, lock_hash_after)) stop("External validation lockfile changed during ensemble enhancement.", call. = FALSE)
message("Final ensemble enhancement completed.")
