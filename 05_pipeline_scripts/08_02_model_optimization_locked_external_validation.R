# ML final selected-model external validation stage
# Sourced by scripts/08_ml_final_optimization.R.

external_prediction_rows <- list()
for (i in seq_len(nrow(selected))) {
  combo <- selected[i, , drop = FALSE]
  td <- task_data_store[[combo$task]]
  combo_grid <- grid[grid$combination_id == combo$combination_id, , drop = FALSE][1, , drop = FALSE]
  combo_grid$selected_role <- combo$selected_role
  ext_pred <- try(automl_final_external_evaluate_combination(td, combo_grid, td$all_sig, lock_hash_before, seed = (config$project$seed %||% 20260501) + i), silent = TRUE)
  if (inherits(ext_pred, "try-error")) {
    warning(sprintf("External evaluation failed for %s: %s", combo$combination_id, as.character(ext_pred)))
  } else {
    external_prediction_rows[[length(external_prediction_rows) + 1]] <- ext_pred
  }
}
if (length(external_prediction_rows) == 0) stop("No selected final model could be externally evaluated.", call. = FALSE)
external_predictions <- do.call(rbind, external_prediction_rows)
write_tsv_atomic(external_predictions, file.path(out_dir, "external_predictions.tsv.gz"))

ext_groups <- split(external_predictions, paste(external_predictions$task, external_predictions$selected_role, external_predictions$combination_id, sep = "__"))
ext_summary <- do.call(rbind, lapply(ext_groups, function(d) {
  m <- compute_binary_metrics(d$true_label, d$score)
  cal <- calibration_summary(d$true_label, d$score)
  data.frame(
    task = d$task[1], selected_role = d$selected_role[1], normalization = d$normalization[1], selector = d$selector[1], classifier = d$classifier[1],
    selector_classifier = d$selector_classifier[1], combination_id = d$combination_id[1], n = m$n, n_pos = m$n_pos,
    AUROC = m$AUROC, AUPRC = m$AUPRC, sensitivity = m$sensitivity, specificity = m$specificity,
    balanced_accuracy = m$balanced_accuracy, Brier = m$Brier,
    calibration_intercept = cal$calibration_intercept, calibration_slope = cal$calibration_slope,
    stringsAsFactors = FALSE
  )
}))
write_tsv_atomic(ext_summary, file.path(out_dir, "external_pooled_summary.tsv"))

ext_family <- do.call(rbind, lapply(ext_groups, function(d) {
  fam <- summarize_per_study_metrics(d, "cohort_family")
  fam$task <- d$task[1]
  fam$selected_role <- d$selected_role[1]
  fam$normalization <- d$normalization[1]
  fam$selector <- d$selector[1]
  fam$classifier <- d$classifier[1]
  fam$selector_classifier <- d$selector_classifier[1]
  fam$combination_id <- d$combination_id[1]
  fam
}))
write_tsv_atomic(ext_family, file.path(out_dir, "external_per_family_performance.tsv"))

best_published <- data.frame()
bench_path <- project_path("results/benchmark/benchmark_rank_table.tsv")
if (file.exists(bench_path)) {
  bench <- read_tsv("results/benchmark/benchmark_rank_table.tsv")
  bench_task <- bench[bench$task %in% tasks, , drop = FALSE]
  best_published <- do.call(rbind, lapply(split(bench_task, bench_task$task), function(d) {
    d <- d[order(d$AUROC, d$AUPRC, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
    data.frame(task = d$task, layer = "best_reference_signature_benchmark", model = d$signature_id, external_AUROC = NA_real_, external_AUPRC = NA_real_, development_or_benchmark_AUROC = d$AUROC, development_or_benchmark_AUPRC = d$AUPRC, stringsAsFactors = FALSE)
  }))
}

final_rows <- data.frame(
  task = ext_summary$task,
  layer = paste0("final_", ext_summary$selected_role),
  model = ext_summary$selector_classifier,
  external_AUROC = ext_summary$AUROC,
  external_AUPRC = ext_summary$AUPRC,
  development_or_benchmark_AUROC = dev_perf$AUROC[match(ext_summary$combination_id, dev_perf$combination_id)],
  development_or_benchmark_AUPRC = dev_perf$AUPRC[match(ext_summary$combination_id, dev_perf$combination_id)],
  stringsAsFactors = FALSE
)
comparison <- rbind(best_published, final_rows)
write_tsv_atomic(comparison, file.path(out_dir, "final_vs_published_signature_summary.tsv"))

lock_hash_after <- hash_file_or_object("config/external_validation_lockfile.json")
write_tsv_atomic(data.frame(
  lockfile_hash_before = lock_hash_before,
  lockfile_hash_after = lock_hash_after,
  unchanged = identical(lock_hash_before, lock_hash_after),
  stringsAsFactors = FALSE
), file.path(out_dir, "external_lockfile_integrity.tsv"))
if (!identical(lock_hash_before, lock_hash_after)) stop("External validation lockfile changed during final optimization.", call. = FALSE)

