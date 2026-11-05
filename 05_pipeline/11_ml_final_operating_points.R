# Final operating-point, uncertainty, and robustness layer.
# Thresholds are selected on development OOF predictions and applied once to locked external predictions.

out_dir <- "results/ml_final"
man_dir <- "manuscript"
dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(file.path(man_dir, "figures")), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(file.path(man_dir, "tables")), recursive = TRUE, showWarnings = FALSE)

selected <- read_tsv(file.path(out_dir, "selected_combinations.tsv"))
oof <- read_tsv(file.path(out_dir, "development_oof_predictions.tsv.gz"))
ext <- read_tsv(file.path(out_dir, "external_predictions.tsv.gz"))

threshold_rows <- list()
external_op_rows <- list()
for (i in seq_len(nrow(selected))) {
  combo <- selected[i, , drop = FALSE]
  dev <- oof[oof$combination_id == combo$combination_id, , drop = FALSE]
  ext_d <- ext[ext$combination_id == combo$combination_id & ext$selected_role == combo$selected_role, , drop = FALSE]
  if (nrow(dev) == 0 || nrow(ext_d) == 0) next
  th <- mlfinal_select_thresholds(dev$true_label, dev$score)
  th$task <- combo$task
  th$selected_role <- combo$selected_role
  th$model <- combo$selector_classifier
  th$combination_id <- combo$combination_id
  threshold_rows[[length(threshold_rows) + 1]] <- th
  for (j in seq_len(nrow(th))) {
    m <- compute_binary_metrics(ext_d$true_label, ext_d$score, threshold = th$threshold[j])
    m$task <- combo$task
    m$selected_role <- combo$selected_role
    m$model <- combo$selector_classifier
    m$combination_id <- combo$combination_id
    m$threshold_name <- th$threshold_name[j]
    m$locked_threshold <- th$threshold[j]
    external_op_rows[[length(external_op_rows) + 1]] <- m
  }
}
threshold_summary <- do.call(rbind, threshold_rows)
external_ops <- do.call(rbind, external_op_rows)
write_tsv_atomic(threshold_summary, file.path(out_dir, "development_locked_thresholds.tsv"))
write_tsv_atomic(external_ops, file.path(out_dir, "external_operating_point_summary.tsv"))
mlfinal_write_display_tsv(external_ops, file.path(man_dir, "tables", "Supplementary_Table_external_operating_points.tsv"))

ext_ci_source <- ext
ext_ci_source$model_id <- paste(ext_ci_source$task, ext_ci_source$selected_role, ext_ci_source$selector_classifier, sep = "__")
ci <- mlfinal_bootstrap_ci_rows(ext_ci_source, "model_id", n_boot = 500, seed = 20260502)
parts <- do.call(rbind, strsplit(ci$model_id, "__", fixed = TRUE))
ci$task <- parts[, 1]
ci$selected_role <- parts[, 2]
ci$model <- parts[, 3]
ci <- ci[, c("task", "selected_role", "model", "AUROC", "AUROC_ci_low", "AUROC_ci_high", "AUPRC", "AUPRC_ci_low", "AUPRC_ci_high"), drop = FALSE]
write_tsv_atomic(ci, file.path(out_dir, "external_auc_auprc_bootstrap_ci.tsv"))
mlfinal_write_display_tsv(ci, file.path(man_dir, "tables", "Supplementary_Table_external_auc_auprc_ci.tsv"))

if (file.exists(project_path(file.path(out_dir, "current_vs_ensemble_external_summary.tsv")))) {
  ens <- read_tsv(file.path(out_dir, "current_vs_ensemble_external_summary.tsv"))
  mlfinal_write_display_tsv(ens, file.path(man_dir, "tables", "Supplementary_Table_current_vs_ensemble_external.tsv"))
}

plot_pdf_png(file.path(man_dir, "figures", "Figure7_operating_points"), function() {
  d <- external_ops[external_ops$selected_role %in% c("best_performance", "best_sparse") & external_ops$threshold_name %in% c("fixed_0.50", "development_youden", "development_sensitivity_0.90", "development_specificity_0.70"), , drop = FALSE]
  d$key <- paste(d$task, d$selected_role, d$threshold_name, sep = "\n")
  graphics::par(mfrow = c(1, 2), mar = c(10, 4, 4, 1))
  vals <- d$sensitivity
  names(vals) <- d$key
  safe_barplot(vals, "Locked external sensitivity by development-locked threshold", ylab = "Sensitivity", col = "steelblue")
  graphics::abline(h = 0.90, lty = 2, col = "firebrick")
  vals <- d$specificity
  names(vals) <- d$key
  safe_barplot(vals, "Locked external specificity by development-locked threshold", ylab = "Specificity", col = "darkseagreen")
  graphics::abline(h = 0.70, lty = 2, col = "firebrick")
}, width = 14, height = 6)

plot_pdf_png(file.path(man_dir, "figures", "Figure8_ensemble_robustness"), function() {
  if (!file.exists(project_path(file.path(out_dir, "current_vs_ensemble_external_summary.tsv")))) {
    graphics::plot.new(); graphics::text(0.5, 0.5, "No ensemble results available"); return(invisible())
  }
  d <- read_tsv(file.path(out_dir, "current_vs_ensemble_external_summary.tsv"))
  d$key <- paste(d$task, d$model_layer, sep = "\n")
  graphics::par(mfrow = c(1, 3), mar = c(9, 4, 4, 1))
  vals <- d$AUROC; names(vals) <- d$key
  safe_barplot(vals, "Current vs development-selected ensemble: AUROC", ylab = "AUROC", col = "steelblue")
  graphics::abline(h = 0.75, lty = 2, col = "gray40")
  vals <- d$AUPRC; names(vals) <- d$key
  safe_barplot(vals, "Current vs development-selected ensemble: AUPRC", ylab = "AUPRC", col = "tan")
  vals <- d$Brier; names(vals) <- d$key
  safe_barplot(vals, "Current vs development-selected ensemble: Brier", ylab = "Brier", col = "tomato")
}, width = 14, height = 5)

write_tsv_atomic(data.frame(
  item = c("ensemble_note", "operating_point_note"),
  value = c(
    "Development-selected ensembles were screened using development OOF predictions only and then evaluated once in locked external validation.",
    "Thresholds were selected on development OOF predictions and applied to locked external predictions."
  ),
  stringsAsFactors = FALSE
), "report/final_enhancement_summary.tsv")

message("Final operating-point and robustness layer completed.")
