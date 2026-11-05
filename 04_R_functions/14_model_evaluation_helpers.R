# ML final model selection and external evaluation
# Split from modeling_algorithm_combinations_final.R for maintainability.

automl_final_task_short <- function(task) {
  if (identical(task, "primary_ATB_vs_OD")) return("primary_ATB_vs_OD")
  if (identical(task, "expanded_ATB_vs_nonATB")) return("expanded_ATB_vs_nonATB")
  gsub("[^A-Za-z0-9]+", "_", task)
}

automl_final_pooled_summary <- function(pred) {
  if (nrow(pred) == 0) return(data.frame())
  groups <- split(pred, pred$combination_id)
  out <- lapply(groups, function(d) {
    m <- compute_binary_metrics(d$true_label, d$score)
    fam <- summarize_per_study_metrics(d, "cohort_family")
    data.frame(
      task = d$task[1], normalization = d$normalization[1], selector = d$selector[1], classifier = d$classifier[1],
      selector_classifier = d$selector_classifier[1], combination_id = d$combination_id[1],
      n_predictions = nrow(d), n_families = length(unique(d$cohort_family)),
      AUROC = m$AUROC, AUPRC = m$AUPRC, balanced_accuracy = m$balanced_accuracy, Brier = m$Brier,
      median_family_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
      n_evaluable_families = sum(!is.na(fam$AUROC)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

automl_final_select_best <- function(perf, sparse_classifiers = c("ridge_logistic", "lasso_logistic", "elastic_net_025", "elastic_net_050", "elastic_net_075")) {
  perf <- perf[is.finite(perf$AUROC), , drop = FALSE]
  if (nrow(perf) == 0) return(list(best = data.frame(), sparse = data.frame()))
  perf <- perf[order(perf$AUROC, perf$AUPRC, perf$median_family_AUROC, decreasing = TRUE), , drop = FALSE]
  best <- perf[1, , drop = FALSE]
  sparse <- perf[perf$classifier %in% sparse_classifiers, , drop = FALSE]
  if (nrow(sparse) > 0) sparse <- sparse[order(sparse$AUROC, sparse$AUPRC, sparse$median_family_AUROC, decreasing = TRUE), , drop = FALSE][1, , drop = FALSE]
  list(best = best, sparse = sparse)
}

automl_final_external_evaluate_combination <- function(task_data, combo_row, all_sample_sig, lock_hash, seed = 20260501) {
  meta <- task_data$meta
  X_raw <- task_data$X_raw
  dev_meta <- task_data$dev_meta
  ext_meta <- task_data$ext_meta
  if (nrow(ext_meta) == 0) stop("no_external_samples", call. = FALSE)
  sig_names <- colnames(all_sample_sig)
  selected <- automl_final_select_features(combo_row$selector, X_raw[, dev_meta$standard_sample_id, drop = FALSE], dev_meta$task_label, sig_names)
  if (!identical(selected$status, "completed")) stop(selected$reason, call. = FALSE)
  sig_dev <- all_sample_sig[dev_meta$standard_sample_id, , drop = FALSE]
  sig_ext <- all_sample_sig[ext_meta$standard_sample_id, , drop = FALSE]
  prepared <- automl_final_prepare_xy(
    X_raw[, dev_meta$standard_sample_id, drop = FALSE], X_raw[, ext_meta$standard_sample_id, drop = FALSE],
    dev_meta$standard_sample_id, ext_meta$standard_sample_id, selected, combo_row$normalization, sig_dev, sig_ext, combo_row$classifier
  )
  score <- automl_final_fit_predict(combo_row$classifier, prepared$train, dev_meta$task_label, prepared$test, seed = seed)
  data.frame(
    task = combo_row$task[1], selected_role = combo_row$selected_role %||% "best_performance",
    standard_sample_id = ext_meta$standard_sample_id, study = ext_meta$study, cohort_family = ext_meta$cohort_family,
    true_label = ext_meta$task_label, score = score,
    normalization = combo_row$normalization, selector = combo_row$selector, classifier = combo_row$classifier,
    selector_classifier = combo_row$selector_classifier, combination_id = combo_row$combination_id,
    n_gene_features = length(selected$genes), n_signature_features = length(selected$signatures),
    n_model_features = ncol(prepared$train), lockfile_hash = lock_hash,
    stringsAsFactors = FALSE
  )
}

