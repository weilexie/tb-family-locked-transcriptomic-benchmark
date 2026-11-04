# Minimal ledger-driven submission updates.
# Uses existing locked results only; no model training, tuning, or selection.

if (!exists("project_path", mode = "function")) {
  Sys.setenv(TBML_PROJECT_ROOT = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  source(file.path("04_reusable_R_functions", "00_load_all.R"))
}

stage <- Sys.getenv("TBML_PACKAGE_UPDATE_STAGE", unset = "standalone")

dir.create(project_path("manuscript/tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("plot_inputs_submission"), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path("quality_control"), recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(as.numeric(x)), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

ci_text <- function(lo, hi) paste0(fmt(lo), "-", fmt(hi))

display_model <- function(x) {
  y <- as.character(x)
  y[y %in% c("reference_signature_scores+svm_rbf", "published_signature_scores+svm_rbf")] <- "Signature-score SVM"
  y[y %in% c("reference_signature_scores+elastic_net_075", "published_signature_scores+elastic_net_075")] <- "Sparse primary model"
  y[y == "stable_top100_plus_signature_scores+random_forest"] <- "Expanded random forest"
  y[y == "stable_top100_plus_signature_scores+ridge_logistic"] <- "Sparse expanded model"
  y[y %in% c("mean_top10", "mean_top10_reference_signatures")] <- "Top-10 signature mean ensemble"
  y[y == "mean_top5"] <- "Top-5 signature mean ensemble"
  y
}

if (stage == "post_sync") {
  message("Manuscript-facing files refreshed.")
} else {
  required <- c(
    "manuscript/tables/Table1_included_cohorts.tsv",
    "manuscript/tables/Table2_reference_signature_locked_external_head_to_head.tsv",
    "manuscript/tables/Table3_final_selected_models.tsv",
    "manuscript/tables/Supplementary_Table_locked_external_performance.tsv",
    "manuscript/tables/Supplementary_Table_family_cluster_bootstrap_ci.tsv",
    "manuscript/tables/Supplementary_Table_paired_family_bootstrap_model_differences.tsv",
    "results/ml_final/external_error_by_sample.tsv",
    "results/ml_final/negative_control_summary.tsv",
    "results/ml_final/negative_control_family_stratified_permutation.tsv"
  )
  missing <- required[!file.exists(project_path(required))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required revision inputs: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  table1 <- read_tsv("manuscript/tables/Table1_included_cohorts.tsv")
  perf <- read_tsv("manuscript/tables/Supplementary_Table_locked_external_performance.tsv")
  errors <- read_tsv("results/ml_final/external_error_by_sample.tsv")
  fam_ci <- read_tsv("manuscript/tables/Supplementary_Table_family_cluster_bootstrap_ci.tsv")

  count_eval <- function(task, role) {
    d <- errors[errors$task == task & errors$selected_role == role, , drop = FALSE]
    if (nrow(d) == 0) return(data.frame(cohort_family = character(), n = integer()))
    out <- stats::aggregate(standard_sample_id ~ cohort_family, d, function(z) length(unique(z)))
    names(out)[names(out) == "standard_sample_id"] <- "n"
    out
  }

  lookup <- function(keys, tab) {
    out <- rep(0L, length(keys))
    m <- match(keys, tab$cohort_family)
    hit <- !is.na(m)
    out[hit] <- as.integer(tab$n[m[hit]])
    out
  }

  primary_counts <- count_eval("primary_ATB_vs_OD", "best_performance")
  expanded_counts <- count_eval("expanded_ATB_vs_nonATB", "best_sparse")
  is_external <- table1$split_primary_ATB_vs_OD == "external_validation"
  table1$primary_model_evaluable_n <- ifelse(is_external, lookup(table1$cohort_family, primary_counts), NA_integer_)
  table1$primary_final_evaluated_n <- table1$primary_model_evaluable_n
  table1$expanded_model_evaluable_n <- ifelse(is_external, lookup(table1$cohort_family, expanded_counts), NA_integer_)
  table1$expanded_final_evaluated_n <- table1$expanded_model_evaluable_n
  table1$sample_count_note <- ifelse(
    is_external,
    "Task-ready counts indicate phenotype eligibility; final evaluated counts indicate locked model-evaluable samples with prediction outputs.",
    "Development or nonlocked family: final evaluated counts are not applicable to locked external performance."
  )
  write_tsv_atomic(table1, "manuscript/tables/Table1_included_cohorts.tsv")

  flow <- data.frame(
    cohort_family = table1$cohort_family,
    split = table1$split_primary_ATB_vs_OD,
    n_samples_harmonized = table1$n_samples_harmonized,
    primary_task_ready_n = table1$primary_ATB_vs_OD_ready_n,
    primary_model_evaluable_n = table1$primary_model_evaluable_n,
    primary_final_evaluated_n = table1$primary_final_evaluated_n,
    expanded_task_ready_n = table1$expanded_ATB_vs_nonATB_ready_n,
    expanded_model_evaluable_n = table1$expanded_model_evaluable_n,
    expanded_final_evaluated_n = table1$expanded_final_evaluated_n,
    stringsAsFactors = FALSE
  )
  flow$primary_exclusion_n <- ifelse(flow$split == "external_validation", flow$primary_task_ready_n - flow$primary_final_evaluated_n, NA_integer_)
  flow$expanded_exclusion_n <- ifelse(flow$split == "external_validation", flow$expanded_task_ready_n - flow$expanded_final_evaluated_n, NA_integer_)
  flow$primary_exclusion_reason <- ifelse(is.na(flow$primary_exclusion_n), "not_applicable",
    ifelse(flow$primary_exclusion_n > 0, "phenotype-eligible but absent from locked prediction outputs after expression/signature feature filtering", "no_attrition"))
  flow$expanded_exclusion_reason <- ifelse(is.na(flow$expanded_exclusion_n), "not_applicable",
    ifelse(flow$expanded_exclusion_n > 0, "phenotype-eligible but absent from locked prediction outputs after expression/signature feature filtering", "no_attrition"))
  write_tsv_atomic(flow, "manuscript/tables/Supplementary_Table_sample_flow_by_family.tsv")

  primary_n <- sum(flow$primary_final_evaluated_n[flow$split == "external_validation"], na.rm = TRUE)
  expanded_n <- sum(flow$expanded_final_evaluated_n[flow$split == "external_validation"], na.rm = TRUE)
  locked_primary_n <- perf$n[perf$task == "primary_ATB_vs_OD" & perf$selected_role == "best_performance"][1]
  locked_expanded_n <- perf$n[perf$task == "expanded_ATB_vs_nonATB" & perf$selected_role == "best_sparse"][1]
  if (!isTRUE(all.equal(primary_n, locked_primary_n))) stop("Primary final evaluated n does not match locked performance n.", call. = FALSE)
  if (!isTRUE(all.equal(expanded_n, locked_expanded_n))) stop("Expanded final evaluated n does not match locked performance n.", call. = FALSE)
  if (any(flow$primary_final_evaluated_n > flow$primary_task_ready_n, na.rm = TRUE) ||
      any(flow$expanded_final_evaluated_n > flow$expanded_task_ready_n, na.rm = TRUE)) {
    stop("Final evaluated n exceeds task-ready n.", call. = FALSE)
  }

  flow_plot <- data.frame(
    stage = c("All harmonized samples", "Primary task-ready", "Primary external task-ready", "Primary final evaluated", "Expanded task-ready", "Expanded external task-ready", "Expanded final evaluated"),
    value = c(sum(table1$n_samples_harmonized, na.rm = TRUE), sum(table1$primary_ATB_vs_OD_ready_n, na.rm = TRUE), sum(flow$primary_task_ready_n[flow$split == "external_validation"], na.rm = TRUE), primary_n, sum(table1$expanded_ATB_vs_nonATB_ready_n, na.rm = TRUE), sum(flow$expanded_task_ready_n[flow$split == "external_validation"], na.rm = TRUE), expanded_n),
    unit = "samples",
    x = c(1, 2, 3, 4, 2, 3, 4),
    note = c("harmonized public cohort set", "phenotype-eligible primary samples", "phenotype-eligible locked-external primary samples", "locked-external primary samples with final prediction outputs", "phenotype-eligible expanded samples", "phenotype-eligible locked-external expanded samples", "locked-external expanded samples with final prediction outputs"),
    stringsAsFactors = FALSE
  )
  flow_plot$label <- paste0(flow_plot$stage, "\n", flow_plot$value, " ", flow_plot$unit)
  write_tsv_atomic(flow_plot, "plot_inputs_submission/Figure2_task_ready_sample_flow.tsv")

  neg <- read_tsv("results/ml_final/negative_control_summary.tsv")
  neg$negative_control <- gsub("random_published_signature_benchmark_spread|random_available_signature_id_score", "random_reference_signature_benchmark_spread", neg$negative_control)
  if ("interpretation" %in% names(neg)) {
    neg$interpretation <- gsub("published", "reference", neg$interpretation, fixed = TRUE)
    neg$interpretation <- gsub("available signature-id", "reference-signature", neg$interpretation, fixed = TRUE)
  }
  neg <- unique(neg)
  write_tsv_atomic(neg, "results/ml_final/negative_control_summary.tsv")
  write_tsv_atomic(neg, "manuscript/tables/Supplementary_Table_negative_control_summary.tsv")

  table2 <- read_tsv("manuscript/tables/Table2_reference_signature_locked_external_head_to_head.tsv")
  paired <- read_tsv("manuscript/tables/Supplementary_Table_paired_family_bootstrap_model_differences.tsv")
  for (nm in c("AUROC_ci_low", "AUROC_ci_high", "AUPRC_ci_low", "AUPRC_ci_high")) {
    if (!nm %in% names(table2)) table2[[nm]] <- NA_real_
  }
  for (i in seq_len(nrow(table2))) {
    if (grepl("^final_", table2$model_layer[i])) {
      role <- sub("^final_", "", table2$model_layer[i])
      hit <- fam_ci[fam_ci$task == table2$analysis_id[i] & fam_ci$selected_role == role, , drop = FALSE]
      if (nrow(hit) > 0) {
        table2$AUROC_ci_low[i] <- hit$AUROC_ci_low[1]
        table2$AUROC_ci_high[i] <- hit$AUROC_ci_high[1]
        table2$AUPRC_ci_low[i] <- hit$AUPRC_ci_low[1]
        table2$AUPRC_ci_high[i] <- hit$AUPRC_ci_high[1]
      }
    }
  }
  table2$n_shared_samples <- table2$n_shared_families <- NA_integer_
  table2$final_AUROC_on_shared_samples <- table2$comparator_AUROC_on_shared_samples <- NA_real_
  table2$final_AUPRC_on_shared_samples <- table2$comparator_AUPRC_on_shared_samples <- NA_real_
  table2$paired_delta_definition <- NA_character_
  for (i in seq_len(nrow(table2))) {
    if (table2$model_layer[i] == "final_best_performance") {
      table2$n_shared_samples[i] <- table2$n_samples[i]
      table2$n_shared_families[i] <- table2$n_families[i]
      table2$final_AUROC_on_shared_samples[i] <- table2$AUROC[i]
      table2$final_AUPRC_on_shared_samples[i] <- table2$AUPRC[i]
      table2$paired_delta_definition[i] <- "reference row for final selected model"
    } else {
      hit <- paired[paired$task == table2$analysis_id[i] & paired$model_b == table2$signature_or_model[i] & paired$model_b_layer == table2$model_layer[i], , drop = FALSE]
      if (nrow(hit) > 0) {
        table2$n_shared_samples[i] <- hit$n_shared_samples[1]
        table2$n_shared_families[i] <- hit$n_shared_families[1]
        table2$comparator_AUROC_on_shared_samples[i] <- table2$AUROC[i]
        table2$comparator_AUPRC_on_shared_samples[i] <- table2$AUPRC[i]
        table2$final_AUROC_on_shared_samples[i] <- table2$AUROC[i] + hit$delta_AUROC[1]
        table2$final_AUPRC_on_shared_samples[i] <- table2$AUPRC[i] + hit$delta_AUPRC[1]
        table2$paired_delta_definition[i] <- "final selected model minus comparator on shared evaluable locked-external samples; negative values favor comparator"
      }
    }
  }
  write_tsv_atomic(table2, "manuscript/tables/Table2_reference_signature_locked_external_head_to_head.tsv")

  table3 <- read_tsv("manuscript/tables/Table3_final_selected_models.tsv")
  table3$external_AUROC_ci_low_family_cluster <- table3$external_AUROC_ci_high_family_cluster <- NA_real_
  table3$external_AUPRC_ci_low_family_cluster <- table3$external_AUPRC_ci_high_family_cluster <- NA_real_
  table3$external_AUROC_family_cluster_95CI <- table3$external_AUPRC_family_cluster_95CI <- NA_character_
  for (i in seq_len(nrow(table3))) {
    hit <- fam_ci[fam_ci$task == table3$analysis_id[i] & fam_ci$selected_role == table3$selected_role[i], , drop = FALSE]
    if (nrow(hit) > 0) {
      table3$external_AUROC_ci_low_family_cluster[i] <- hit$AUROC_ci_low[1]
      table3$external_AUROC_ci_high_family_cluster[i] <- hit$AUROC_ci_high[1]
      table3$external_AUPRC_ci_low_family_cluster[i] <- hit$AUPRC_ci_low[1]
      table3$external_AUPRC_ci_high_family_cluster[i] <- hit$AUPRC_ci_high[1]
      table3$external_AUROC_family_cluster_95CI[i] <- ci_text(hit$AUROC_ci_low[1], hit$AUROC_ci_high[1])
      table3$external_AUPRC_family_cluster_95CI[i] <- ci_text(hit$AUPRC_ci_low[1], hit$AUPRC_ci_high[1])
    }
  }
  write_tsv_atomic(table3, "manuscript/tables/Table3_final_selected_models.tsv")

  metrics <- perf
  metrics$task_label <- ifelse(metrics$task == "primary_ATB_vs_OD", "Primary clinical-differential task", "Expanded screening/generalization task")
  metrics$display_name <- display_model(metrics$selector_classifier)
  metrics$prevalence <- metrics$n_pos / metrics$n
  metrics$threshold <- if ("threshold" %in% names(metrics)) metrics$threshold else 0.5
  metrics$PPV <- if ("PPV" %in% names(metrics)) metrics$PPV else NA_real_
  metrics$NPV <- if ("NPV" %in% names(metrics)) metrics$NPV else NA_real_
  metric_cols <- c("task", "task_label", "display_name", "n", "n_pos", "prevalence", "threshold", "AUROC", "AUPRC", "sensitivity", "specificity", "PPV", "NPV", "balanced_accuracy", "Brier", "selected_role", "selector_classifier", "calibration_slope", "calibration_intercept")
  metrics <- metrics[, metric_cols, drop = FALSE]
  if (file.exists(project_path("results/ml_final/current_vs_ensemble_external_summary.tsv"))) {
    ens <- read_tsv("results/ml_final/current_vs_ensemble_external_summary.tsv")
    ens <- ens[ens$model_layer == "development_selected_ensemble", , drop = FALSE]
    if (nrow(ens) > 0) {
      ens$task_label <- ifelse(ens$task == "primary_ATB_vs_OD", "Primary clinical-differential task", "Expanded screening/generalization task")
      ens$display_name <- display_model(ens$model)
      ens$n <- ifelse(ens$task == "primary_ATB_vs_OD", primary_n, expanded_n)
      ens$n_pos <- perf$n_pos[match(ens$task, perf$task)]
      ens$prevalence <- ens$n_pos / ens$n
      ens$threshold <- 0.5
      ens$selected_role <- ens$model_layer
      ens$selector_classifier <- ens$model
      ens$PPV <- ens$NPV <- ens$calibration_intercept <- NA_real_
      ens$balanced_accuracy <- rowMeans(ens[, c("sensitivity", "specificity")], na.rm = TRUE)
      metrics <- rbind(metrics, ens[, metric_cols, drop = FALSE])
    }
  }
  write_tsv_atomic(metrics, "plot_inputs_submission/Figure5_calibration_summary.tsv")

  audit <- data.frame(
    issue_id = c("P0-01", "P0-02", "P0-03", "P0-04", "P1-07", "P1-08", "P2-01"),
    status = "implemented",
    validation_result = c(
      sprintf("primary=%s; expanded=%s", primary_n, expanded_n),
      ifelse(any(duplicated(neg)), "failed", "passed"),
      "shared-sample paired columns added",
      sprintf("primary AUROC CI %s; expanded AUROC CI %s",
        ci_text(fam_ci$AUROC_ci_low[fam_ci$task == "primary_ATB_vs_OD" & fam_ci$selected_role == "best_performance"][1], fam_ci$AUROC_ci_high[fam_ci$task == "primary_ATB_vs_OD" & fam_ci$selected_role == "best_performance"][1]),
        ci_text(fam_ci$AUROC_ci_low[fam_ci$task == "expanded_ATB_vs_nonATB" & fam_ci$selected_role == "best_sparse"][1], fam_ci$AUROC_ci_high[fam_ci$task == "expanded_ATB_vs_nonATB" & fam_ci$selected_role == "best_sparse"][1])),
      "Figure5 input cleaned",
      "Figure6 generated from selected-rule inputs",
      "manifest generated during manuscript_final post-sync"
    ),
    stringsAsFactors = FALSE
  )
  if (file.exists(project_path("manuscript/tables/Table_index.tsv"))) {
    index <- read_tsv("manuscript/tables/Table_index.tsv")
    additions <- data.frame(
      table_file = "Supplementary_Table_sample_flow_by_family.tsv",
      placement = "supplementary",
      title = "Sample flow by cohort family",
      purpose = "Reconciles task-ready and final evaluated counts.",
      source_table = "Table1_included_cohorts.tsv; external_error_by_sample.tsv",
      notes = "Included in Supplementary Table S1.",
      exists_in_package = TRUE,
      stringsAsFactors = FALSE
    )
    for (nm in names(index)) if (!nm %in% names(additions)) additions[[nm]] <- NA
    for (nm in names(additions)) if (!nm %in% names(index)) index[[nm]] <- NA
    index <- index[!index$table_file %in% additions$table_file, , drop = FALSE]
    index <- rbind(index[, names(additions), drop = FALSE], additions)
    index$exists_in_package <- file.exists(project_path(file.path("manuscript/tables", index$table_file)))
    write_tsv_atomic(index, "manuscript/tables/Table_index.tsv")
  }

  message(sprintf("Submission revision table inputs completed (%s).", stage))
}
