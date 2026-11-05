# Generalization stress-test analyses for the conservative revision.
# This stage does not change the external lockfile or select thresholds on external data.

out_dir <- "results/ml_final"
man_table_dir <- "manuscript/tables"
report_dir <- "report"
dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(man_table_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(report_dir), recursive = TRUE, showWarnings = FALSE)

config <- read_config("config/config.yml")
lock_hash_before <- hash_file_or_object("config/external_validation_lockfile.json")
lock <- read_external_lockfile("config/external_validation_lockfile.json")
development_families <- unlist(lock$development_families)
external_families <- unlist(lock$external_validation_families)
tasks <- c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB")

selected <- read_tsv(file.path(out_dir, "selected_combinations.tsv"))
dev_perf <- read_tsv(file.path(out_dir, "development_combination_performance.tsv"))
ext_summary <- read_tsv(file.path(out_dir, "external_pooled_summary.tsv"))
ext_pred <- read_tsv(file.path(out_dir, "external_predictions.tsv.gz"))
ext_family <- read_tsv(file.path(out_dir, "external_per_family_performance.tsv"))

# 1. Published signatures on the exact same locked external families.
signature_external <- data.frame()
signature_external_family <- data.frame()
head_to_head <- data.frame()
score_path <- "results/benchmark/signature_scores_by_sample.tsv.gz"
bench_path <- "results/benchmark/benchmark_rank_table.tsv"
if (file.exists(project_path(score_path)) && file.exists(project_path(bench_path))) {
  signature_scores <- read_tsv(score_path)
  bench <- read_tsv(bench_path)
  registry <- if (file.exists(project_path("results/benchmark/signature_registry_resolved.tsv"))) {
    read_tsv("results/benchmark/signature_registry_resolved.tsv")
  } else {
    data.frame(signature_id = unique(signature_scores$signature_id), stringsAsFactors = FALSE)
  }
  overlap <- if (file.exists(project_path("results/benchmark/signature_training_overlap_by_study.tsv"))) {
    read_tsv("results/benchmark/signature_training_overlap_by_study.tsv")
  } else {
    data.frame(signature_id = character(), cohort_family = character(), stringsAsFactors = FALSE)
  }
  overlap$overlap_key <- paste(overlap$signature_id, overlap$cohort_family, sep = "__")

  rows <- list()
  fam_rows <- list()
  for (task in tasks) {
    task_scores <- signature_scores[signature_scores$task == task & signature_scores$cohort_family %in% external_families, , drop = FALSE]
    for (sig in sort(unique(task_scores$signature_id))) {
      d_all <- task_scores[task_scores$signature_id == sig, , drop = FALSE]
      sig_overlap_keys <- paste(sig, unique(d_all$cohort_family), sep = "__")
      overlap_families <- unique(d_all$cohort_family[sig_overlap_keys %in% overlap$overlap_key])
      scopes <- list(
        all_locked_external = d_all,
        external_no_signature_training_family_overlap = d_all[!d_all$cohort_family %in% overlap_families, , drop = FALSE]
      )
      for (scope_name in names(scopes)) {
        d <- scopes[[scope_name]]
        m <- mlaudit_auc_row(d$true_label, d$score)
        rows[[length(rows) + 1]] <- data.frame(
          task = task,
          signature_id = sig,
          evaluation_scope = scope_name,
          n_external_families = length(unique(d$cohort_family)),
          n_external_training_overlap_families = length(overlap_families),
          n = m$n,
          n_pos = m$n_pos,
          prevalence = m$prevalence,
          external_AUROC = m$AUROC,
          external_AUPRC = m$AUPRC,
          development_or_benchmark_AUROC = bench$AUROC[match(paste(task, sig), paste(bench$task, bench$signature_id))],
          development_or_benchmark_AUPRC = bench$AUPRC[match(paste(task, sig), paste(bench$task, bench$signature_id))],
          signature_level = registry$signature_level[match(sig, registry$signature_id)] %||% NA_character_,
          intended_task = registry$intended_task[match(sig, registry$signature_id)] %||% NA_character_,
          stringsAsFactors = FALSE
        )
        for (fam in sort(unique(d$cohort_family))) {
          df <- d[d$cohort_family == fam, , drop = FALSE]
          mf <- mlaudit_auc_row(df$true_label, df$score)
          fam_rows[[length(fam_rows) + 1]] <- data.frame(
            task = task,
            signature_id = sig,
            evaluation_scope = scope_name,
            cohort_family = fam,
            n = mf$n,
            n_pos = mf$n_pos,
            AUROC = mf$AUROC,
            AUPRC = mf$AUPRC,
            signature_training_family_overlap = paste(sig, fam, sep = "__") %in% overlap$overlap_key,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  signature_external <- do.call(rbind, rows)
  signature_external <- signature_external[order(signature_external$task, signature_external$evaluation_scope, -signature_external$external_AUROC, -signature_external$external_AUPRC, na.last = TRUE), , drop = FALSE]
  signature_external_family <- do.call(rbind, fam_rows)
  write_tsv_atomic(signature_external, file.path(out_dir, "published_signature_external_head_to_head.tsv"))
  write_tsv_atomic(signature_external_family, file.path(out_dir, "published_signature_external_per_family.tsv"))
  unlink(project_path(file.path(man_table_dir, "Supplementary_Table_published_signature_external_head_to_head.tsv")), force = TRUE)
  mlfinal_write_display_tsv(signature_external, file.path(man_table_dir, "Supplementary_Table_reference_signature_external_head_to_head.tsv"))

  comparison_path <- file.path(out_dir, "final_vs_published_signature_summary.tsv")
  if (file.exists(project_path(comparison_path))) {
    comparison <- read_tsv(comparison_path)
    idx <- which(comparison$layer %in% c("best_reference_signature_benchmark", "best_published_signature_benchmark"))
    for (i in idx) {
      key <- signature_external$task == comparison$task[i] &
        signature_external$signature_id == comparison$model[i] &
        signature_external$evaluation_scope == "all_locked_external"
      if (any(key)) {
        comparison$external_AUROC[i] <- signature_external$external_AUROC[which(key)[1]]
        comparison$external_AUPRC[i] <- signature_external$external_AUPRC[which(key)[1]]
      }
    }
    write_tsv_atomic(comparison, comparison_path)
  }

  final_h2h <- data.frame(
    task = ext_summary$task,
    model_layer = paste("final", ext_summary$selected_role, sep = "_"),
    model = ext_summary$selector_classifier,
    evaluation_scope = "all_locked_external",
    n = ext_summary$n,
    n_pos = ext_summary$n_pos,
    external_AUROC = ext_summary$AUROC,
    external_AUPRC = ext_summary$AUPRC,
    development_or_benchmark_AUROC = dev_perf$AUROC[match(ext_summary$combination_id, dev_perf$combination_id)],
    development_or_benchmark_AUPRC = dev_perf$AUPRC[match(ext_summary$combination_id, dev_perf$combination_id)],
    stringsAsFactors = FALSE
  )
  prespecified_signatures <- c("Sweeney_OD_3", "Kaforou_OD_44", "Kaforou_OD_53", "PennNich_RISK_6", "Zak_RISK_16")
  published_top <- do.call(rbind, lapply(split(signature_external[signature_external$evaluation_scope == "all_locked_external", ], signature_external$task[signature_external$evaluation_scope == "all_locked_external"]), function(d) {
    d_ranked <- d[is.finite(d$development_or_benchmark_AUROC), , drop = FALSE]
    d_ranked <- d_ranked[order(d_ranked$development_or_benchmark_AUROC, d_ranked$development_or_benchmark_AUPRC, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    selected_ids <- unique(c(head(d_ranked$signature_id, min(10, nrow(d_ranked))), intersect(prespecified_signatures, d$signature_id)))
    out <- d[d$signature_id %in% selected_ids, , drop = FALSE]
    out$selection_basis <- ifelse(out$signature_id %in% head(d_ranked$signature_id, min(10, nrow(d_ranked))), "top_development_benchmark", "prespecified_signature")
    out
  }))
  published_h2h <- data.frame(
    task = published_top$task,
    model_layer = "published_signature",
    model = published_top$signature_id,
    evaluation_scope = published_top$evaluation_scope,
    n = published_top$n,
    n_pos = published_top$n_pos,
    external_AUROC = published_top$external_AUROC,
    external_AUPRC = published_top$external_AUPRC,
    development_or_benchmark_AUROC = published_top$development_or_benchmark_AUROC,
    development_or_benchmark_AUPRC = published_top$development_or_benchmark_AUPRC,
    selection_basis = published_top$selection_basis,
    stringsAsFactors = FALSE
  )
  final_h2h$selection_basis <- "development_selected_final_model"
  head_to_head <- rbind(final_h2h, published_h2h)
  task_best_published <- stats::aggregate(external_AUROC ~ task, published_h2h, max, na.rm = TRUE)
  names(task_best_published)[2] <- "best_published_external_AUROC"
  head_to_head <- merge(head_to_head, task_best_published, by = "task", all.x = TRUE)
  head_to_head$delta_AUROC_vs_best_published <- head_to_head$external_AUROC - head_to_head$best_published_external_AUROC
  head_to_head <- head_to_head[order(head_to_head$task, -head_to_head$external_AUROC, -head_to_head$external_AUPRC), , drop = FALSE]
  write_tsv_atomic(head_to_head, file.path(out_dir, "head_to_head_final_vs_published_external.tsv"))
  mlfinal_write_display_tsv(head_to_head, file.path(man_table_dir, "Supplementary_Table_locked_external_head_to_head.tsv"))
  compact_comparison <- data.frame(
    task = head_to_head$task,
    layer = head_to_head$model_layer,
    model = head_to_head$model,
    external_AUROC = head_to_head$external_AUROC,
    external_AUPRC = head_to_head$external_AUPRC,
    development_or_benchmark_AUROC = head_to_head$development_or_benchmark_AUROC,
    development_or_benchmark_AUPRC = head_to_head$development_or_benchmark_AUPRC,
    selection_basis = head_to_head$selection_basis,
    stringsAsFactors = FALSE
  )
  write_tsv_atomic(compact_comparison, file.path(out_dir, "final_vs_published_signature_summary.tsv"))
}

# 2. primary ATB-vs-OD analysis mixed-family-only sensitivity based on development OOF predictions.
split_balance <- read_tsv("results/tables/development_validation_split_balance.tsv")
split_balance$is_mixed_ATB_OD_family <- split_balance$ATB > 0 & split_balance$OD > 0
split_balance$is_single_class_family <- xor(split_balance$ATB > 0, split_balance$OD > 0)
split_balance$positive_fraction <- split_balance$ATB / pmax(split_balance$ATB + split_balance$OD, 1)
split_balance$risk_note <- ifelse(split_balance$split == "development" & split_balance$is_single_class_family,
  "development_single_class_family_supports_separation_but_not_family_AUROC",
  ifelse(split_balance$is_mixed_ATB_OD_family, "mixed_ATB_OD_family_evaluable", "empty_or_not_task_ready")
)
write_tsv_atomic(split_balance, file.path(out_dir, "primary_ATB_vs_OD_family_label_structure_check.tsv"))

oof <- read_tsv(file.path(out_dir, "development_oof_predictions.tsv.gz"))
primary_ATB_vs_OD_mixed_dev <- split_balance$cohort_family[split_balance$split == "development" & split_balance$is_mixed_ATB_OD_family]
primary_ATB_vs_OD_mixed_ext <- split_balance$cohort_family[split_balance$split == "external_validation" & split_balance$is_mixed_ATB_OD_family]
write_tsv_atomic(data.frame(
  split = c("development", "external_validation"),
  n_mixed_families = c(length(primary_ATB_vs_OD_mixed_dev), length(primary_ATB_vs_OD_mixed_ext)),
  mixed_families = c(paste(primary_ATB_vs_OD_mixed_dev, collapse = ";"), paste(primary_ATB_vs_OD_mixed_ext, collapse = ";")),
  stringsAsFactors = FALSE
), file.path(out_dir, "primary_ATB_vs_OD_mixed_family_inventory.tsv"))

primary_ATB_vs_OD_oof_mixed <- oof[oof$task == "primary_ATB_vs_OD" & oof$cohort_family %in% primary_ATB_vs_OD_mixed_dev, , drop = FALSE]
mixed_groups <- split(primary_ATB_vs_OD_oof_mixed, primary_ATB_vs_OD_oof_mixed$combination_id)
mixed_perf <- do.call(rbind, lapply(mixed_groups, function(d) {
  m <- compute_binary_metrics(d$true_label, d$score)
  fam <- summarize_per_study_metrics(d, "cohort_family")
  data.frame(
    combination_id = d$combination_id[1],
    task = d$task[1],
    selector_classifier = d$selector_classifier[1],
    mixed_development_n = m$n,
    mixed_development_n_pos = m$n_pos,
    mixed_development_AUROC = m$AUROC,
    mixed_development_AUPRC = m$AUPRC,
    mixed_development_Brier = m$Brier,
    mixed_development_median_family_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
    mixed_development_n_evaluable_families = sum(is.finite(fam$AUROC)),
    stringsAsFactors = FALSE
  )
}))
all_primary_ATB_vs_OD <- dev_perf[dev_perf$task == "primary_ATB_vs_OD", , drop = FALSE]
all_primary_ATB_vs_OD <- all_primary_ATB_vs_OD[order(all_primary_ATB_vs_OD$AUROC, all_primary_ATB_vs_OD$AUPRC, all_primary_ATB_vs_OD$median_family_AUROC, decreasing = TRUE), , drop = FALSE]
all_primary_ATB_vs_OD$all_development_rank <- seq_len(nrow(all_primary_ATB_vs_OD))
mixed_perf <- merge(mixed_perf, all_primary_ATB_vs_OD[, c("combination_id", "AUROC", "AUPRC", "median_family_AUROC", "all_development_rank")], by = "combination_id", all.x = TRUE)
names(mixed_perf)[names(mixed_perf) == "AUROC"] <- "all_development_AUROC"
names(mixed_perf)[names(mixed_perf) == "AUPRC"] <- "all_development_AUPRC"
names(mixed_perf)[names(mixed_perf) == "median_family_AUROC"] <- "all_development_median_family_AUROC"
mixed_perf <- mixed_perf[order(mixed_perf$mixed_development_AUROC, mixed_perf$mixed_development_AUPRC, mixed_perf$mixed_development_median_family_AUROC, decreasing = TRUE), , drop = FALSE]
mixed_perf$mixed_development_rank <- seq_len(nrow(mixed_perf))
mixed_perf$currently_selected <- mixed_perf$combination_id %in% selected$combination_id[selected$task == "primary_ATB_vs_OD"]
write_tsv_atomic(mixed_perf, file.path(out_dir, "primary_ATB_vs_OD_mixed_family_development_rank_shift.tsv"))
write_tsv_atomic(mixed_perf[mixed_perf$currently_selected | mixed_perf$mixed_development_rank <= 30, , drop = FALSE], file.path(man_table_dir, "Supplementary_Table_primary_ATB_vs_OD_mixed_family_rank_shift.tsv"))

selected_primary_ATB_vs_OD <- selected[selected$task == "primary_ATB_vs_OD", , drop = FALSE]
lmfo_rows <- list()
for (i in seq_len(nrow(selected_primary_ATB_vs_OD))) {
  combo_id <- selected_primary_ATB_vs_OD$combination_id[i]
  d <- primary_ATB_vs_OD_oof_mixed[primary_ATB_vs_OD_oof_mixed$combination_id == combo_id, , drop = FALSE]
  for (fam in primary_ATB_vs_OD_mixed_dev) {
    held <- d[d$cohort_family == fam, , drop = FALSE]
    pooled_except <- d[d$cohort_family != fam, , drop = FALSE]
    mh <- compute_binary_metrics(held$true_label, held$score)
    mp <- compute_binary_metrics(pooled_except$true_label, pooled_except$score)
    lmfo_rows[[length(lmfo_rows) + 1]] <- data.frame(
      selected_role = selected_primary_ATB_vs_OD$selected_role[i],
      model = selected_primary_ATB_vs_OD$selector_classifier[i],
      combination_id = combo_id,
      heldout_mixed_family = fam,
      heldout_family_n = mh$n,
      heldout_family_AUROC = mh$AUROC,
      heldout_family_AUPRC = mh$AUPRC,
      pooled_except_family_n = mp$n,
      pooled_except_family_AUROC = mp$AUROC,
      pooled_except_family_AUPRC = mp$AUPRC,
      stringsAsFactors = FALSE
    )
  }
}
primary_ATB_vs_OD_lmfo <- do.call(rbind, lmfo_rows)
write_tsv_atomic(primary_ATB_vs_OD_lmfo, file.path(out_dir, "primary_ATB_vs_OD_leave_one_mixed_family_out_sensitivity.tsv"))
write_tsv_atomic(primary_ATB_vs_OD_lmfo, file.path(man_table_dir, "Supplementary_Table_primary_ATB_vs_OD_leave_one_mixed_family_out.tsv"))

# 3. Negative controls and confounding checks.
selected_groups <- split(ext_pred, paste(ext_pred$task, ext_pred$selected_role, ext_pred$combination_id, sep = "__"))
n_perm <- config$negative_controls$n_permutations %||% 100
perm_rows <- do.call(rbind, lapply(selected_groups, function(d) {
  p <- mlaudit_score_label_permutation(d, n_perm = n_perm, seed = 20260502)
  data.frame(
    task = d$task[1],
    selected_role = d$selected_role[1],
    model = d$selector_classifier[1],
    combination_id = d$combination_id[1],
    p,
    stringsAsFactors = FALSE
  )
}))
write_tsv_atomic(perm_rows, file.path(out_dir, "negative_control_score_label_permutation.tsv"))
write_tsv_atomic(perm_rows, file.path(man_table_dir, "Supplementary_Table_negative_control_score_label_permutation.tsv"))

platform_rows <- list()
for (task in tasks) {
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  meta$split <- ifelse(meta$cohort_family %in% external_families, "external_validation",
    ifelse(meta$cohort_family %in% development_families, "development", "not_in_lock"))
  if (!"platform_type" %in% names(meta)) meta$platform_type <- "unknown"
  meta$platform_type[is.na(meta$platform_type) | !nzchar(as.character(meta$platform_type))] <- "unknown"
  dev <- meta[meta$split == "development", , drop = FALSE]
  platform_prev <- tapply(as_binary(dev$task_label), dev$platform_type, mean)
  default_prev <- mean(as_binary(dev$task_label))
  for (sp in c("development", "external_validation")) {
    d <- meta[meta$split == sp, , drop = FALSE]
    if (nrow(d) == 0) next
    score <- as.numeric(platform_prev[as.character(d$platform_type)])
    score[!is.finite(score) | is.na(score)] <- default_prev
    m <- mlaudit_auc_row(d$task_label, score)
    platform_rows[[length(platform_rows) + 1]] <- data.frame(
      task = task,
      split = sp,
      n = m$n,
      n_pos = m$n_pos,
      n_platform_levels = length(unique(d$platform_type)),
      platform_label_cramers_v = mlaudit_cramers_v(d$platform_type, d$task_label),
      development_platform_prevalence_score_AUROC = m$AUROC,
      development_platform_prevalence_score_AUPRC = m$AUPRC,
      platform_levels = paste(sort(unique(d$platform_type)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
platform_check <- do.call(rbind, platform_rows)
write_tsv_atomic(platform_check, file.path(out_dir, "platform_label_confounding_check.tsv"))
write_tsv_atomic(platform_check, file.path(man_table_dir, "Supplementary_Table_platform_label_confounding_check.tsv"))

negative_summary <- data.frame(
  negative_control = "selected_external_score_label_permutation",
  n_iterations = n_perm,
  AUROC_median = stats::median(perm_rows$permuted_AUROC_median, na.rm = TRUE),
  AUROC_q025 = stats::median(perm_rows$permuted_AUROC_q025, na.rm = TRUE),
  AUROC_q975 = stats::median(perm_rows$permuted_AUROC_q975, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_tsv_atomic(negative_summary, file.path(out_dir, "negative_control_summary.tsv"))
write_tsv_atomic(negative_summary, file.path(man_table_dir, "Supplementary_Table_negative_control_summary.tsv"))

# 4. Corrected primary ATB-vs-OD analysis error composition from the final selected external model.
primary_ATB_vs_OD_best <- ext_pred[ext_pred$task == "primary_ATB_vs_OD" & ext_pred$selected_role == "best_performance", , drop = FALSE]
if (nrow(primary_ATB_vs_OD_best) > 0) {
  y <- as_binary(primary_ATB_vs_OD_best$true_label)
  pred_lab <- ifelse(primary_ATB_vs_OD_best$score >= 0.5, 1L, 0L)
  primary_ATB_vs_OD_best$prediction_status <- ifelse(pred_lab == 1L & y == 1L, "TP",
    ifelse(pred_lab == 0L & y == 0L, "TN", ifelse(pred_lab == 1L & y == 0L, "FP", "FN"))
  )
  err_comp <- as.data.frame(table(primary_ATB_vs_OD_best$true_label, primary_ATB_vs_OD_best$prediction_status), stringsAsFactors = FALSE)
  names(err_comp) <- c("true_label", "prediction_status", "n")
  err_comp <- err_comp[err_comp$n > 0, , drop = FALSE]
  write_tsv_atomic(err_comp, file.path(out_dir, "external_error_composition_primary_ATB_vs_OD_best_performance.tsv"))
  err_by_family <- as.data.frame(table(primary_ATB_vs_OD_best$cohort_family, primary_ATB_vs_OD_best$true_label, primary_ATB_vs_OD_best$prediction_status), stringsAsFactors = FALSE)
  names(err_by_family) <- c("cohort_family", "true_label", "prediction_status", "n")
  err_by_family <- err_by_family[err_by_family$n > 0, , drop = FALSE]
  write_tsv_atomic(err_by_family, file.path(out_dir, "external_error_composition_by_family_primary_ATB_vs_OD_best_performance.tsv"))
}

lock_hash_after <- hash_file_or_object("config/external_validation_lockfile.json")
if (!identical(lock_hash_before, lock_hash_after)) stop("External validation lockfile changed during generalization stress-test stage.", call. = FALSE)

primary_ATB_vs_OD_head <- if (nrow(head_to_head) > 0) head_to_head[head_to_head$task == "primary_ATB_vs_OD", , drop = FALSE] else data.frame()
primary_ATB_vs_OD_final <- primary_ATB_vs_OD_head[primary_ATB_vs_OD_head$model_layer == "final_best_performance", , drop = FALSE]
primary_ATB_vs_OD_pub <- primary_ATB_vs_OD_head[primary_ATB_vs_OD_head$model_layer == "published_signature", , drop = FALSE]
primary_ATB_vs_OD_mixed_selected <- mixed_perf[mixed_perf$currently_selected, , drop = FALSE]
family_ci_for_report <- if (file.exists(project_path("results/ml_final/family_cluster_bootstrap_ci.tsv"))) {
  read_tsv("results/ml_final/family_cluster_bootstrap_ci.tsv")
} else {
  data.frame()
}
primary_ATB_vs_OD_ci <- family_ci_for_report[family_ci_for_report$task == "primary_ATB_vs_OD" & family_ci_for_report$selected_role == "best_performance", , drop = FALSE]
write_tsv_atomic(data.frame(
  item = c("primary_final_external_AUROC", "primary_final_external_AUPRC", "best_published_signature_external_AUROC", "mixed_development_families", "family_cluster_AUROC_ci_low", "family_cluster_AUROC_ci_high", "manuscript_implication"),
  value = c(
    if (nrow(primary_ATB_vs_OD_final) > 0) primary_ATB_vs_OD_final$external_AUROC[1] else NA_real_,
    if (nrow(primary_ATB_vs_OD_final) > 0) primary_ATB_vs_OD_final$external_AUPRC[1] else NA_real_,
    if (nrow(primary_ATB_vs_OD_pub) > 0) max(primary_ATB_vs_OD_pub$external_AUROC, na.rm = TRUE) else NA_real_,
    paste(primary_ATB_vs_OD_mixed_dev, collapse = ";"),
    if (nrow(primary_ATB_vs_OD_ci) > 0) primary_ATB_vs_OD_ci$AUROC_ci_low[1] else NA_real_,
    if (nrow(primary_ATB_vs_OD_ci) > 0) primary_ATB_vs_OD_ci$AUROC_ci_high[1] else NA_real_,
    "Benchmark plus optimization plus generalization ceiling; not a deployable diagnostic signature."
  ),
  stringsAsFactors = FALSE
), "report/generalization_stress_test_response.tsv")

message("Generalization stress-test analyses completed.")
