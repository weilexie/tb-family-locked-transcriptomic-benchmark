# ML final summary-report stage
# Sourced by scripts/08_model_optimization.R.

completed_n <- sum(grid_status$status == "completed")
skipped_n <- sum(grid_status$status == "skipped")
primary_ATB_vs_OD_ext <- ext_summary[ext_summary$task == "primary_ATB_vs_OD" & ext_summary$selected_role == "best_performance", , drop = FALSE]
expanded_ATB_vs_nonATB_ext <- ext_summary[ext_summary$task == "expanded_ATB_vs_nonATB" & ext_summary$selected_role == "best_performance", , drop = FALSE]
claim <- if (nrow(primary_ATB_vs_OD_ext) > 0 && is.finite(primary_ATB_vs_OD_ext$AUROC[1]) && primary_ATB_vs_OD_ext$AUROC[1] >= 0.75) {
  "primary ATB-vs-OD analysis reached the pre-specified AUROC upgrade threshold for an optimized ATB-vs-OD multicohort ML model."
} else if (nrow(expanded_ATB_vs_nonATB_ext) > 0 && is.finite(expanded_ATB_vs_nonATB_ext$AUROC[1]) && expanded_ATB_vs_nonATB_ext$AUROC[1] >= 0.80) {
  "expanded ATB-vs-nonATB analysis reached the pre-specified AUROC threshold for an expanded-sample screening/general discrimination model; primary ATB-vs-OD analysis remains the stricter clinical-discrimination line."
} else {
  "The controlled large-scale ML screen supports a benchmark/generalizability framing rather than a strong superiority claim."
}

write_tsv_atomic(data.frame(
  item = c("planned_combinations", "completed_combinations", "skipped_combinations", "tasks", "normalizations", "claim", "external_lock_unchanged", "development_only_selection", "separate_primary_and_expanded_claims"),
  value = c(nrow(grid), completed_n, skipped_n, paste(tasks, collapse = ";"), paste(normalizations, collapse = ";"), claim, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
), "report/ml_algorithm_combinations_final_summary.tsv")

message(sprintf("ML final optimization completed: %s completed, %s skipped", completed_n, skipped_n))

