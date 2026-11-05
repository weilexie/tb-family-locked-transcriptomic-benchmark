# Additional stress-test analyses for the conservative submission revision.
# This stage does not train new classifiers on external data and does not use
# locked external labels for model or threshold selection.

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
oof <- read_tsv(file.path(out_dir, "development_oof_predictions.tsv.gz"))
thresholds <- if (file.exists(project_path(file.path(out_dir, "development_locked_thresholds.tsv")))) {
  read_tsv(file.path(out_dir, "development_locked_thresholds.tsv"))
} else {
  data.frame()
}

gst_ece <- function(y_true, prob, n_bins = 10) {
  y <- as_binary(y_true)
  p <- pmin(pmax(as.numeric(prob), 1e-6), 1 - 1e-6)
  keep <- !is.na(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  if (length(y) == 0) return(NA_real_)
  br <- cut(p, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  bins <- split(seq_along(p), br)
  sum(vapply(bins, function(idx) {
    if (length(idx) == 0) return(0)
    length(idx) / length(p) * abs(mean(y[idx]) - mean(p[idx]))
  }, numeric(1)))
}

gst_calibration_metrics <- function(y_true, prob, threshold = 0.5) {
  m <- compute_binary_metrics(y_true, prob, threshold = threshold)
  cal <- calibration_summary(y_true, prob)
  data.frame(
    n = m$n,
    n_pos = m$n_pos,
    threshold = threshold,
    AUROC = m$AUROC,
    AUPRC = m$AUPRC,
    sensitivity = m$sensitivity,
    specificity = m$specificity,
    Brier = m$Brier,
    ECE = gst_ece(y_true, prob),
    calibration_intercept = cal$calibration_intercept,
    calibration_slope = cal$calibration_slope,
    stringsAsFactors = FALSE
  )
}

gst_fit_isotonic <- function(y_true, score) {
  y <- as_binary(y_true)
  x <- as.numeric(score)
  keep <- !is.na(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  if (length(y) < 10 || length(unique(y)) < 2 || length(unique(x)) < 2) return(list(ok = FALSE))
  ord <- order(x)
  fit <- try(stats::isoreg(x[ord], y[ord]), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(ok = FALSE))
  xx <- fit$x
  yy <- pmin(pmax(fit$yf, 1e-6), 1 - 1e-6)
  u <- !duplicated(xx, fromLast = TRUE)
  list(ok = TRUE, x = xx[u], y = yy[u])
}

gst_predict_isotonic <- function(fit, score) {
  if (!isTRUE(fit$ok)) return(pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6))
  as.numeric(stats::approx(fit$x, fit$y, xout = as.numeric(score), rule = 2, ties = "ordered")$y)
}

gst_fit_temperature <- function(y_true, score) {
  y <- as_binary(y_true)
  z <- stats::qlogis(pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6))
  keep <- !is.na(y) & is.finite(z)
  y <- y[keep]
  z <- z[keep]
  if (length(y) < 10 || length(unique(y)) < 2 || length(unique(z)) < 2) return(list(ok = FALSE, temperature = NA_real_))
  nll <- function(log_t) {
    t <- exp(log_t)
    p <- stats::plogis(z / t)
    -sum(y * log(pmin(pmax(p, 1e-6), 1 - 1e-6)) + (1 - y) * log(pmin(pmax(1 - p, 1e-6), 1 - 1e-6)))
  }
  opt <- try(stats::optim(0, nll, method = "Brent", lower = log(0.05), upper = log(20)), silent = TRUE)
  if (inherits(opt, "try-error")) return(list(ok = FALSE, temperature = NA_real_))
  list(ok = TRUE, temperature = exp(opt$par))
}

gst_predict_temperature <- function(fit, score) {
  if (!isTRUE(fit$ok) || !is.finite(fit$temperature)) return(pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6))
  z <- stats::qlogis(pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6))
  stats::plogis(z / fit$temperature)
}

gst_fit_beta <- function(y_true, score) {
  y <- as_binary(y_true)
  p <- pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6)
  x1 <- log(p)
  x2 <- log(1 - p)
  keep <- !is.na(y) & is.finite(x1) & is.finite(x2)
  y <- y[keep]
  x1 <- x1[keep]
  x2 <- x2[keep]
  if (length(y) < 10 || length(unique(y)) < 2) return(list(ok = FALSE, fit = NULL))
  fit <- try(stats::glm(y ~ x1 + x2, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(ok = FALSE, fit = NULL))
  list(ok = TRUE, fit = fit)
}

gst_predict_beta <- function(fit, score) {
  if (!isTRUE(fit$ok) || is.null(fit$fit)) return(pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6))
  p <- pmin(pmax(as.numeric(score), 1e-6), 1 - 1e-6)
  x1 <- log(p)
  x2 <- log(1 - p)
  as.numeric(stats::predict(fit$fit, newdata = data.frame(x1 = x1, x2 = x2), type = "response"))
}

gst_group_bootstrap <- function(pred, threshold = 0.5, group_col = "cohort_family", n_boot = 1000, seed = 20260502) {
  groups <- sort(unique(pred[[group_col]]))
  obs <- compute_binary_metrics(pred$true_label, pred$score, threshold = threshold)
  if (length(groups) < 2) {
    return(data.frame(
      n = obs$n, n_pos = obs$n_pos, n_clusters = length(groups),
      AUROC = obs$AUROC, AUROC_ci_low = NA_real_, AUROC_ci_high = NA_real_,
      AUPRC = obs$AUPRC, AUPRC_ci_low = NA_real_, AUPRC_ci_high = NA_real_,
      sensitivity = obs$sensitivity, sensitivity_ci_low = NA_real_, sensitivity_ci_high = NA_real_,
      specificity = obs$specificity, specificity_ci_low = NA_real_, specificity_ci_high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  set.seed(seed)
  vals <- replicate(n_boot, {
    sampled <- sample(groups, length(groups), replace = TRUE)
    boot <- do.call(rbind, lapply(seq_along(sampled), function(i) {
      d <- pred[pred[[group_col]] == sampled[i], , drop = FALSE]
      d$cluster_draw <- i
      d
    }))
    m <- compute_binary_metrics(boot$true_label, boot$score, threshold = threshold)
    c(AUROC = m$AUROC, AUPRC = m$AUPRC, sensitivity = m$sensitivity, specificity = m$specificity)
  })
  data.frame(
    n = obs$n,
    n_pos = obs$n_pos,
    n_clusters = length(groups),
    AUROC = obs$AUROC,
    AUROC_ci_low = as.numeric(stats::quantile(vals["AUROC", ], 0.025, na.rm = TRUE)),
    AUROC_ci_high = as.numeric(stats::quantile(vals["AUROC", ], 0.975, na.rm = TRUE)),
    AUPRC = obs$AUPRC,
    AUPRC_ci_low = as.numeric(stats::quantile(vals["AUPRC", ], 0.025, na.rm = TRUE)),
    AUPRC_ci_high = as.numeric(stats::quantile(vals["AUPRC", ], 0.975, na.rm = TRUE)),
    sensitivity = obs$sensitivity,
    sensitivity_ci_low = as.numeric(stats::quantile(vals["sensitivity", ], 0.025, na.rm = TRUE)),
    sensitivity_ci_high = as.numeric(stats::quantile(vals["sensitivity", ], 0.975, na.rm = TRUE)),
    specificity = obs$specificity,
    specificity_ci_low = as.numeric(stats::quantile(vals["specificity", ], 0.025, na.rm = TRUE)),
    specificity_ci_high = as.numeric(stats::quantile(vals["specificity", ], 0.975, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

gst_equal_family_summary <- function(pred, threshold = 0.5) {
  fam <- summarize_per_study_metrics(pred, "cohort_family")
  data.frame(
    n_families = length(unique(fam$cohort_family)),
    n_evaluable_auc_families = sum(is.finite(fam$AUROC)),
    pooled_AUROC = compute_roc_auc(pred$true_label, pred$score),
    pooled_AUPRC = compute_pr_auc(pred$true_label, pred$score),
    mean_family_AUROC = mean(fam$AUROC, na.rm = TRUE),
    median_family_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
    family_AUROC_iqr_low = as.numeric(stats::quantile(fam$AUROC, 0.25, na.rm = TRUE)),
    family_AUROC_iqr_high = as.numeric(stats::quantile(fam$AUROC, 0.75, na.rm = TRUE)),
    family_weighted_sensitivity = mean(fam$sensitivity, na.rm = TRUE),
    family_weighted_specificity = mean(fam$specificity, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

gst_within_family_score_permutation <- function(pred, n_perm = 1000, seed = 20260502) {
  families <- sort(unique(pred$cohort_family))
  set.seed(seed)
  vals <- replicate(n_perm, {
    d <- pred
    for (fam in families) {
      idx <- which(d$cohort_family == fam)
      if (length(idx) > 1) d$score[idx] <- sample(d$score[idx], length(idx), replace = FALSE)
    }
    compute_roc_auc(d$true_label, d$score)
  })
  data.frame(
    negative_control = "within_family_score_permutation",
    n_iterations = n_perm,
    AUROC_median = stats::median(vals, na.rm = TRUE),
    AUROC_q025 = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
    AUROC_q975 = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
    interpretation = "permutes scores within external cohort families to preserve family score distributions",
    stringsAsFactors = FALSE
  )
}

gst_threshold_for <- function(pred_row) {
  if (nrow(thresholds) == 0) return(0.5)
  hit <- thresholds[thresholds$combination_id == pred_row$combination_id[1] &
    thresholds$selected_role == pred_row$selected_role[1] &
    thresholds$threshold_name == "fixed_0.50", , drop = FALSE]
  if (nrow(hit) == 0) return(0.5)
  as.numeric(hit$threshold[1])
}

gst_join_metadata <- function(pred) {
  rows <- list()
  for (task in unique(pred$task)) {
    meta_path <- file.path("data_curated/task_matrices", task, "metadata.tsv")
    if (!file.exists(project_path(meta_path))) next
    meta <- read_tsv(meta_path)
    keep_cols <- intersect(c(
      "standard_sample_id", "platform", "platform_type", "country", "geographical_region",
      "age", "age_group", "sex", "hiv_status", "diabetes_status", "tissue",
      "control_type_detail", "original_diagnosis", "original_control_type",
      "harmonized_label_4class", "mapping_confidence"
    ), names(meta))
    d <- pred[pred$task == task, , drop = FALSE]
    d <- merge(d, meta[, keep_cols, drop = FALSE], by = "standard_sample_id", all.x = TRUE)
    rows[[length(rows) + 1]] <- d
  }
  out <- do.call(rbind, rows)
  for (nm in c("platform", "platform_type", "country", "geographical_region", "age_group", "sex", "hiv_status", "diabetes_status", "tissue", "control_type_detail")) {
    if (!nm %in% names(out)) out[[nm]] <- NA_character_
    out[[nm]][is.na(out[[nm]]) | !nzchar(as.character(out[[nm]]))] <- "unknown"
  }
  out$OD_subtype <- ifelse(out$true_label == "OD", out$control_type_detail, NA_character_)
  out
}

gst_selective_metrics <- function(y_true, score, lower, upper) {
  y <- as_binary(y_true)
  s <- as.numeric(score)
  keep <- !is.na(y) & is.finite(s)
  y <- y[keep]
  s <- s[keep]
  answered <- s <= lower | s >= upper
  pred <- ifelse(s >= upper, 1L, ifelse(s <= lower, 0L, NA_integer_))
  tp <- sum(answered & pred == 1L & y == 1L)
  tn <- sum(answered & pred == 0L & y == 0L)
  fp <- sum(answered & pred == 1L & y == 0L)
  fn <- sum(answered & pred == 0L & y == 1L)
  answered_n <- sum(answered)
  uncertain_n <- sum(!answered)
  data.frame(
    n = length(y),
    n_pos = sum(y == 1L),
    lower_threshold = lower,
    upper_threshold = upper,
    answered_n = answered_n,
    uncertain_n = uncertain_n,
    coverage = answered_n / length(y),
    abstention_rate = uncertain_n / length(y),
    error_rate_among_answered = ifelse(answered_n > 0, (fp + fn) / answered_n, NA_real_),
    sensitivity_answered = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
    specificity_answered = ifelse(tn + fp > 0, tn / (tn + fp), NA_real_),
    PPV_answered = ifelse(tp + fp > 0, tp / (tp + fp), NA_real_),
    NPV_answered = ifelse(tn + fn > 0, tn / (tn + fn), NA_real_),
    false_negative_rate_answered = ifelse(tp + fn > 0, fn / (tp + fn), NA_real_),
    false_positive_rate_answered = ifelse(tn + fp > 0, fp / (tn + fp), NA_real_),
    tp_answered = tp, tn_answered = tn, fp_answered = fp, fn_answered = fn,
    stringsAsFactors = FALSE
  )
}

gst_pick_strategy <- function(grid, strategy) {
  d <- grid
  if (strategy == "dev_min_error_at_coverage_0.70") {
    d <- d[d$dev_coverage >= 0.70, , drop = FALSE]
    if (nrow(d) == 0) return(integer())
    return(which(grid$key == d$key[order(d$dev_error_rate_among_answered, -d$dev_coverage, -d$dev_balanced_accuracy_answered)[1]]))
  }
  if (strategy == "dev_specificity_ge_0.70_max_coverage") {
    d <- d[d$dev_specificity_answered >= 0.70, , drop = FALSE]
    if (nrow(d) == 0) return(integer())
    return(which(grid$key == d$key[order(-d$dev_coverage, d$dev_error_rate_among_answered, -d$dev_sensitivity_answered)[1]]))
  }
  if (strategy == "dev_sensitivity_ge_0.90_max_coverage") {
    d <- d[d$dev_sensitivity_answered >= 0.90, , drop = FALSE]
    if (nrow(d) == 0) return(integer())
    return(which(grid$key == d$key[order(-d$dev_coverage, d$dev_error_rate_among_answered, -d$dev_specificity_answered)[1]]))
  }
  if (strategy == "dev_balanced_accuracy_max") {
    return(which.max(grid$dev_balanced_accuracy_answered))
  }
  integer()
}

gst_entropy <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  p <- as.numeric(table(x)) / length(x)
  -sum(p * log(p))
}

gst_family_oracle_row <- function(d, family_col, task, split_name, definition_label) {
  y <- as_binary(d$task_label)
  fam <- as.character(d[[family_col]])
  keep <- !is.na(y) & !is.na(fam) & nzchar(fam)
  d <- d[keep, , drop = FALSE]
  y <- y[keep]
  fam <- fam[keep]
  fam_prev <- ave(y, fam, FUN = mean)
  fam_tab <- tapply(y, fam, function(z) length(unique(z[!is.na(z)])) == 2)
  data.frame(
    task = task,
    split = split_name,
    family_definition = definition_label,
    n = length(y),
    n_families = length(unique(fam)),
    n_mixed_families = sum(fam_tab, na.rm = TRUE),
    n_single_class_families = sum(!fam_tab, na.rm = TRUE),
    family_oracle_AUROC = compute_roc_auc(y, fam_prev),
    family_label_cramers_v = mlaudit_cramers_v(fam, d$task_label),
    median_family_phenotype_entropy = stats::median(tapply(d$task_label, fam, gst_entropy), na.rm = TRUE),
    analysis_scope = "post_hoc_regrouping_no_model_refit",
    stringsAsFactors = FALSE
  )
}

gst_family_phenotype_rows <- function(d, family_col, task, split_name, definition_label) {
  fam <- as.character(d[[family_col]])
  keep <- !is.na(fam) & nzchar(fam)
  d <- d[keep, , drop = FALSE]
  fam <- fam[keep]
  y <- as_binary(d$task_label)
  rows <- lapply(sort(unique(fam)), function(ff) {
    idx <- fam == ff
    pos_n <- sum(y[idx] == 1L, na.rm = TRUE)
    neg_n <- sum(y[idx] == 0L, na.rm = TRUE)
    fam_type <- if (pos_n > 0 && neg_n > 0) {
      if (min(pos_n, neg_n) < 10) "mixed_low_volume" else "mixed"
    } else if (pos_n > 0) {
      "positive_only"
    } else if (neg_n > 0) {
      "negative_only"
    } else {
      "not_task_ready"
    }
    data.frame(
      task = task,
      split = split_name,
      family_definition = definition_label,
      family_id = ff,
      n = sum(idx),
      positive_n = pos_n,
      negative_n = neg_n,
      positive_fraction = ifelse(pos_n + neg_n > 0, pos_n / (pos_n + neg_n), NA_real_),
      phenotype_entropy = gst_entropy(d$task_label[idx]),
      family_type = fam_type,
      low_volume_flag = fam_type == "mixed_low_volume",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gst_map_family_definitions <- function(d, crosswalk) {
  keep <- c(
    "standard_sample_id",
    "family_geo_accession_level",
    "family_parent_study_level_main",
    "family_conservative_exclusion",
    "family_maximally_conservative_merge"
  )
  merge(d, unique(crosswalk[, keep, drop = FALSE]), by = "standard_sample_id", all.x = TRUE)
}

gst_within_family_perm_for_col <- function(pred, family_col, n_perm = 1000, seed = 20260502) {
  d <- pred[!is.na(pred[[family_col]]) & nzchar(as.character(pred[[family_col]])), , drop = FALSE]
  if (nrow(d) == 0) {
    return(data.frame(
      negative_control = "within_family_score_permutation",
      n_iterations = n_perm,
      AUROC_median = NA_real_,
      AUROC_q025 = NA_real_,
      AUROC_q975 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  d$cohort_family <- as.character(d[[family_col]])
  gst_within_family_score_permutation(d, n_perm = n_perm, seed = seed)[
    ,
    c("negative_control", "n_iterations", "AUROC_median", "AUROC_q025", "AUROC_q975"),
    drop = FALSE
  ]
}

gst_lite_split_sensitivity <- function(task, selected_role, n_repeats = 500, seed = 20260530) {
  sel <- selected[selected$task == task & selected$selected_role == selected_role, , drop = FALSE]
  if (nrow(sel) != 1) return(list(results = data.frame(), summary = data.frame()))
  combo_id <- sel$combination_id[1]
  model_oof <- oof[oof$task == task & oof$combination_id == combo_id & oof$cohort_family %in% development_families, , drop = FALSE]
  if (nrow(model_oof) == 0) return(list(results = data.frame(), summary = data.frame()))
  bench <- if (file.exists(project_path("results/benchmark/benchmark_rank_table_development.tsv"))) {
    read_tsv("results/benchmark/benchmark_rank_table_development.tsv")
  } else {
    read_tsv("results/benchmark/benchmark_rank_table.tsv")
  }
  task_bench <- bench[bench$task == task, , drop = FALSE]
  task_bench <- task_bench[order(task_bench$AUROC, task_bench$AUPRC, decreasing = TRUE), , drop = FALSE]
  top_sigs <- head(task_bench$signature_id, 10)
  task_sig <- sig_scores[sig_scores$task == task & sig_scores$signature_id %in% top_sigs & sig_scores$cohort_family %in% development_families, , drop = FALSE]
  families <- sort(unique(model_oof$cohort_family))
  split_size <- min(length(external_families), length(families))
  if (split_size < 2 || length(unique(as_binary(model_oof$true_label))) < 2) return(list(results = data.frame(), summary = data.frame()))
  set.seed(seed)
  rows <- list()
  attempts <- 0L
  while (length(rows) < n_repeats && attempts < n_repeats * 80L) {
    attempts <- attempts + 1L
    held <- sort(sample(families, split_size, replace = FALSE))
    dm <- model_oof[model_oof$cohort_family %in% held, , drop = FALSE]
    if (length(unique(as_binary(dm$true_label))) < 2) next
    model_auc <- compute_roc_auc(dm$true_label, dm$score)
    model_auprc <- compute_pr_auc(dm$true_label, dm$score)
    best_ref_auc <- NA_real_
    best_ref_auprc <- NA_real_
    best_ref <- NA_character_
    best_shared_n <- NA_integer_
    for (sig in top_sigs) {
      ds <- task_sig[task_sig$signature_id == sig & task_sig$standard_sample_id %in% dm$standard_sample_id, , drop = FALSE]
      if (nrow(ds) == 0 || length(unique(as_binary(ds$true_label))) < 2) next
      auc <- compute_roc_auc(ds$true_label, ds$score)
      auprc <- compute_pr_auc(ds$true_label, ds$score)
      if (!is.finite(best_ref_auc) || (is.finite(auc) && auc > best_ref_auc)) {
        best_ref_auc <- auc
        best_ref_auprc <- auprc
        best_ref <- sig
        best_shared_n <- nrow(ds)
      }
    }
    if (!is.finite(best_ref_auc)) next
    rows[[length(rows) + 1]] <- data.frame(
      task = task,
      selected_role = selected_role,
      model = sel$selector_classifier[1],
      combination_id = combo_id,
      split_id = sprintf("lite_%03d", length(rows) + 1L),
      heldout_family_count = length(held),
      heldout_families = paste(held, collapse = ";"),
      n = nrow(dm),
      n_pos = sum(as_binary(dm$true_label) == 1L, na.rm = TRUE),
      model_AUROC = model_auc,
      model_AUPRC = model_auprc,
      best_reference_signature = best_ref,
      best_reference_AUROC = best_ref_auc,
      best_reference_AUPRC = best_ref_auprc,
      n_shared_reference_samples = best_shared_n,
      delta_ML_vs_best_reference_AUROC = model_auc - best_ref_auc,
      analysis_scope = "development_oof_family_block_lite_no_refit",
      stringsAsFactors = FALSE
    )
  }
  results <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  summary <- if (nrow(results) > 0) {
    data.frame(
      task = task,
      selected_role = selected_role,
      model = sel$selector_classifier[1],
      combination_id = combo_id,
      n_splits = nrow(results),
      heldout_family_count = split_size,
      AUROC_median = stats::median(results$model_AUROC, na.rm = TRUE),
      AUROC_q025 = as.numeric(stats::quantile(results$model_AUROC, 0.025, na.rm = TRUE)),
      AUROC_q975 = as.numeric(stats::quantile(results$model_AUROC, 0.975, na.rm = TRUE)),
      delta_ML_vs_best_ref_median = stats::median(results$delta_ML_vs_best_reference_AUROC, na.rm = TRUE),
      delta_ML_vs_best_ref_q025 = as.numeric(stats::quantile(results$delta_ML_vs_best_reference_AUROC, 0.025, na.rm = TRUE)),
      delta_ML_vs_best_ref_q975 = as.numeric(stats::quantile(results$delta_ML_vs_best_reference_AUROC, 0.975, na.rm = TRUE)),
      prop_ML_beats_best_ref = mean(results$delta_ML_vs_best_reference_AUROC > 0, na.rm = TRUE),
      prop_delta_gt_0_02 = mean(results$delta_ML_vs_best_reference_AUROC > 0.02, na.rm = TRUE),
      prop_delta_gt_0_05 = mean(results$delta_ML_vs_best_reference_AUROC > 0.05, na.rm = TRUE),
      analysis_scope = "development_oof_family_block_lite_no_refit",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }
  list(results = results, summary = summary)
}

gst_family_calibration_threshold_audit <- function(d) {
  th <- gst_threshold_for(d)
  rows <- lapply(sort(unique(d$cohort_family)), function(fam) {
    df <- d[d$cohort_family == fam, , drop = FALSE]
    m <- compute_binary_metrics(df$true_label, df$score, threshold = th)
    cal <- calibration_summary(df$true_label, df$score)
    pred_pos <- mean(as.numeric(df$score) >= th, na.rm = TRUE)
    data.frame(
      task = df$task[1],
      selected_role = df$selected_role[1],
      model = df$selector_classifier[1],
      combination_id = df$combination_id[1],
      cohort_family = fam,
      n = m$n,
      n_pos = m$n_pos,
      prevalence = m$n_pos / max(m$n, 1),
      threshold = th,
      AUROC = m$AUROC,
      AUPRC = m$AUPRC,
      sensitivity = m$sensitivity,
      specificity = m$specificity,
      PPV = m$PPV,
      NPV = m$NPV,
      Brier = m$Brier,
      ECE = gst_ece(df$true_label, df$score),
      calibration_intercept = cal$calibration_intercept,
      calibration_slope = cal$calibration_slope,
      calibration_estimable = is.finite(cal$calibration_slope),
      mean_predicted_probability = mean(as.numeric(df$score), na.rm = TRUE),
      predicted_positive_fraction = pred_pos,
      all_positive_like = pred_pos >= 0.95,
      all_negative_like = pred_pos <= 0.05,
      severe_positive_bias = pred_pos >= 0.85,
      severe_negative_bias = pred_pos <= 0.15,
      analysis_scope = "locked_external_family_specific_no_threshold_refit",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gst_update_table_index <- function(additions) {
  index_path <- "manuscript/tables/Table_index.tsv"
  if (!file.exists(project_path(index_path))) return(invisible(FALSE))
  index <- read_tsv(index_path)
  for (nm in names(index)) if (!nm %in% names(additions)) additions[[nm]] <- NA
  for (nm in names(additions)) if (!nm %in% names(index)) index[[nm]] <- NA
  additions <- additions[, names(index), drop = FALSE]
  index <- index[!index$table_file %in% additions$table_file, , drop = FALSE]
  index <- rbind(index, additions)
  index$exists_in_package <- file.exists(project_path(file.path("manuscript/tables", index$table_file)))
  write_tsv_atomic(index, index_path)
  invisible(TRUE)
}

# 1. Sample-level external error table with task metadata.
ext_meta <- gst_join_metadata(ext_pred)
threshold_rows <- do.call(rbind, lapply(split(ext_meta, paste(ext_meta$task, ext_meta$selected_role, ext_meta$combination_id, sep = "__")), function(d) {
  th <- gst_threshold_for(d)
  y <- as_binary(d$true_label)
  pred <- ifelse(d$score >= th, 1L, 0L)
  d$sample_id <- d$standard_sample_id
  d$predicted_probability <- d$score
  d$locked_threshold <- th
  d$predicted_label <- ifelse(pred == 1L, "ATB", ifelse(d$task[1] == "primary_ATB_vs_OD", "OD", "nonATB"))
  d$error_type <- ifelse(pred == 1L & y == 1L, "TP",
    ifelse(pred == 0L & y == 0L, "TN", ifelse(pred == 1L & y == 0L, "FP", "FN"))
  )
  d
}))
err_cols <- intersect(c(
  "sample_id", "standard_sample_id", "study", "cohort_family", "task", "selected_role",
  "true_label", "predicted_probability", "locked_threshold", "predicted_label",
  "error_type", "platform", "platform_type", "country", "geographical_region",
  "hiv_status", "age", "age_group", "sex", "diabetes_status", "tissue", "OD_subtype",
  "control_type_detail", "original_diagnosis", "original_control_type",
  "selector_classifier", "combination_id"
), names(threshold_rows))
external_error_by_sample <- threshold_rows[, err_cols, drop = FALSE]
write_tsv_atomic(external_error_by_sample, file.path(out_dir, "external_error_by_sample.tsv"))
write_tsv_atomic(external_error_by_sample, file.path(man_table_dir, "Supplementary_Table_external_error_by_sample.tsv"))

primary_ATB_vs_OD_best_err <- external_error_by_sample[
  external_error_by_sample$task == "primary_ATB_vs_OD" &
    external_error_by_sample$selected_role == "best_performance", , drop = FALSE
]
write_tsv_atomic(primary_ATB_vs_OD_best_err, file.path(out_dir, "external_error_by_sample_primary_ATB_vs_OD_best_performance.tsv"))
err_comp <- as.data.frame(table(primary_ATB_vs_OD_best_err$true_label, primary_ATB_vs_OD_best_err$error_type), stringsAsFactors = FALSE)
names(err_comp) <- c("true_label", "prediction_status", "n")
err_comp <- err_comp[err_comp$n > 0, , drop = FALSE]
write_tsv_atomic(err_comp, file.path(out_dir, "external_error_composition_primary_ATB_vs_OD_best_performance.tsv"))

subgroup_vars <- intersect(c("cohort_family", "platform_type", "geographical_region", "hiv_status", "age_group", "OD_subtype"), names(primary_ATB_vs_OD_best_err))
sub_rows <- list()
for (v in subgroup_vars) {
  tmp <- primary_ATB_vs_OD_best_err
  tmp[[v]][is.na(tmp[[v]]) | !nzchar(as.character(tmp[[v]]))] <- "not_applicable_or_unknown"
  tab <- as.data.frame(table(tmp[[v]], tmp$error_type), stringsAsFactors = FALSE)
  names(tab) <- c("level", "error_type", "n")
  tab <- tab[tab$n > 0, , drop = FALSE]
  tab$variable <- v
  sub_rows[[length(sub_rows) + 1]] <- tab[, c("variable", "level", "error_type", "n"), drop = FALSE]
}
err_meta <- do.call(rbind, sub_rows)
write_tsv_atomic(err_meta, file.path(out_dir, "external_error_composition_by_metadata_primary_ATB_vs_OD_best_performance.tsv"))

score_dist <- do.call(rbind, lapply(split(primary_ATB_vs_OD_best_err, primary_ATB_vs_OD_best_err$error_type), function(d) {
  data.frame(
    error_type = d$error_type[1],
    n = nrow(d),
    score_median = stats::median(d$predicted_probability, na.rm = TRUE),
    score_iqr_low = as.numeric(stats::quantile(d$predicted_probability, 0.25, na.rm = TRUE)),
    score_iqr_high = as.numeric(stats::quantile(d$predicted_probability, 0.75, na.rm = TRUE)),
    score_min = min(d$predicted_probability, na.rm = TRUE),
    score_max = max(d$predicted_probability, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv_atomic(score_dist, file.path(out_dir, "external_score_distribution_by_error_primary_ATB_vs_OD_best_performance.tsv"))

# 2. Published-signature locked external head-to-head with family metrics and overlap status.
sig_scores <- read_tsv("results/benchmark/signature_scores_by_sample.tsv.gz")
registry <- read_tsv("results/benchmark/signature_registry_resolved.tsv")
overlap <- read_tsv("results/benchmark/signature_training_overlap_by_study.tsv")
bench <- if (file.exists(project_path("results/benchmark/benchmark_rank_table_development.tsv"))) {
  read_tsv("results/benchmark/benchmark_rank_table_development.tsv")
} else {
  read_tsv("results/benchmark/benchmark_rank_table.tsv")
}
overlap$overlap_key <- paste(overlap$signature_id, overlap$cohort_family, sep = "__")
known_training <- unique(c(overlap$signature_id, registry$signature_id[nzchar(as.character(registry$training_studies_mapped))]))
prespecified_signatures <- c("Sweeney_OD_3", "Kaforou_OD_44", "Kaforou_OD_53", "PennNich_RISK_6", "Zak_RISK_16")

head_rows <- list()
for (task in tasks) {
  task_bench <- bench[bench$task == task, , drop = FALSE]
  task_bench <- task_bench[order(task_bench$AUROC, task_bench$AUPRC, decreasing = TRUE), , drop = FALSE]
  sig_ids <- unique(c(head(task_bench$signature_id, 10), prespecified_signatures))
  task_scores <- sig_scores[sig_scores$task == task & sig_scores$cohort_family %in% external_families & sig_scores$signature_id %in% sig_ids, , drop = FALSE]
  for (sig in sort(unique(task_scores$signature_id))) {
    d <- task_scores[task_scores$signature_id == sig, , drop = FALSE]
    m <- gst_calibration_metrics(d$true_label, d$score, threshold = 0.5)
    fam <- summarize_per_study_metrics(d, "cohort_family")
    n_overlap <- length(unique(d$cohort_family[paste(sig, d$cohort_family, sep = "__") %in% overlap$overlap_key]))
    training_status <- ifelse(sig %in% known_training, "known_training_cohorts", "unknown_training_cohorts")
    overlap_status <- ifelse(n_overlap > 0, "known_external_training_overlap",
      ifelse(training_status == "known_training_cohorts", "no_known_external_training_overlap", "unknown_training_cohorts"))
    head_rows[[length(head_rows) + 1]] <- data.frame(
      signature_or_model = sig,
      model_layer = "published_signature",
      task = task,
      n_samples = m$n,
      n_families = length(unique(d$cohort_family)),
      AUROC = m$AUROC,
      AUPRC = m$AUPRC,
      family_median_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
      family_IQR_AUROC = paste0(sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.25, na.rm = TRUE))), "-", sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.75, na.rm = TRUE)))),
      sensitivity_at_locked_threshold = m$sensitivity,
      specificity_at_locked_threshold = m$specificity,
      Brier = m$Brier,
      calibration_slope = m$calibration_slope,
      overlap_status = overlap_status,
      n_external_training_overlap_families = n_overlap,
      selection_basis = ifelse(sig %in% head(task_bench$signature_id, 10), "top10_development_benchmark", "prespecified_signature"),
      stringsAsFactors = FALSE
    )
  }

  top_for_ens <- head(task_bench$signature_id[task_bench$signature_id %in% unique(sig_scores$signature_id)], 10)
  if (length(top_for_ens) >= 2) {
    wide <- reshape(sig_scores[sig_scores$task == task & sig_scores$signature_id %in% top_for_ens, ],
      idvar = c("standard_sample_id", "study", "cohort_family", "true_label"),
      timevar = "signature_id", direction = "wide")
    score_cols <- grep("^score\\.", names(wide), value = TRUE)
    dev_w <- wide[wide$cohort_family %in% development_families, , drop = FALSE]
    ext_w <- wide[wide$cohort_family %in% external_families, , drop = FALSE]
    ext_scores <- list(
      mean_top10_published_signatures = rowMeans(ext_w[, score_cols, drop = FALSE], na.rm = TRUE),
      max_top10_published_signatures = apply(ext_w[, score_cols, drop = FALSE], 1, max, na.rm = TRUE)
    )
    complete_dev <- stats::complete.cases(dev_w[, score_cols, drop = FALSE])
    complete_ext <- stats::complete.cases(ext_w[, score_cols, drop = FALSE])
    if (sum(complete_dev) >= 20 && length(unique(as_binary(dev_w$true_label[complete_dev]))) == 2 && sum(complete_ext) > 0) {
      form <- stats::as.formula(paste("as_binary(true_label) ~", paste(sprintf("`%s`", score_cols), collapse = " + ")))
      fit <- try(stats::glm(form, data = dev_w[complete_dev, c("true_label", score_cols), drop = FALSE], family = stats::binomial()), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        ext_scores$logistic_top10_published_signatures <- as.numeric(stats::predict(fit, newdata = ext_w[, score_cols, drop = FALSE], type = "response"))
      }
    }
    for (ens_name in names(ext_scores)) {
      de <- ext_w
      de$score <- ext_scores[[ens_name]]
      de <- de[is.finite(de$score), , drop = FALSE]
      m <- gst_calibration_metrics(de$true_label, de$score, threshold = 0.5)
      fam <- summarize_per_study_metrics(de, "cohort_family")
      head_rows[[length(head_rows) + 1]] <- data.frame(
        signature_or_model = ens_name,
        model_layer = "published_signature_ensemble",
        task = task,
        n_samples = m$n,
        n_families = length(unique(de$cohort_family)),
        AUROC = m$AUROC,
        AUPRC = m$AUPRC,
        family_median_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
        family_IQR_AUROC = paste0(sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.25, na.rm = TRUE))), "-", sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.75, na.rm = TRUE)))),
        sensitivity_at_locked_threshold = m$sensitivity,
        specificity_at_locked_threshold = m$specificity,
        Brier = m$Brier,
        calibration_slope = m$calibration_slope,
        overlap_status = "development_fit_or_rule_based_signature_ensemble",
        n_external_training_overlap_families = NA_integer_,
        selection_basis = "top10_development_signature_ensemble",
        stringsAsFactors = FALSE
      )
    }
  }
}
for (i in seq_len(nrow(ext_summary))) {
  d <- ext_pred[ext_pred$combination_id == ext_summary$combination_id[i] & ext_pred$selected_role == ext_summary$selected_role[i], , drop = FALSE]
  fam <- summarize_per_study_metrics(d, "cohort_family")
  head_rows[[length(head_rows) + 1]] <- data.frame(
    signature_or_model = ext_summary$selector_classifier[i],
    model_layer = paste("final", ext_summary$selected_role[i], sep = "_"),
    task = ext_summary$task[i],
    n_samples = ext_summary$n[i],
    n_families = length(unique(d$cohort_family)),
    AUROC = ext_summary$AUROC[i],
    AUPRC = ext_summary$AUPRC[i],
    family_median_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
    family_IQR_AUROC = paste0(sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.25, na.rm = TRUE))), "-", sprintf("%.3f", as.numeric(stats::quantile(fam$AUROC, 0.75, na.rm = TRUE)))),
    sensitivity_at_locked_threshold = ext_summary$sensitivity[i],
    specificity_at_locked_threshold = ext_summary$specificity[i],
    Brier = ext_summary$Brier[i],
    calibration_slope = ext_summary$calibration_slope[i],
    overlap_status = "not_applicable_final_model",
    n_external_training_overlap_families = NA_integer_,
    selection_basis = "development_selected_final_model",
    stringsAsFactors = FALSE
  )
}
published_locked_h2h <- do.call(rbind, head_rows)
published_locked_h2h <- published_locked_h2h[order(published_locked_h2h$task, -published_locked_h2h$AUROC, -published_locked_h2h$AUPRC), , drop = FALSE]
write_tsv_atomic(published_locked_h2h, file.path(out_dir, "published_signature_locked_external_head_to_head.tsv"))
unlink(project_path(file.path(man_table_dir, "Supplementary_Table_published_signature_locked_external_head_to_head.tsv")), force = TRUE)
mlfinal_write_display_tsv(published_locked_h2h, file.path(man_table_dir, "Supplementary_Table_reference_signature_locked_external_head_to_head.tsv"))

# 3. Discovery/training-overlap-aware external benchmark.
overlap_rows <- list()
for (task in tasks) {
  task_scores <- sig_scores[sig_scores$task == task & sig_scores$cohort_family %in% external_families, , drop = FALSE]
  for (sig in sort(unique(task_scores$signature_id))) {
    d0 <- task_scores[task_scores$signature_id == sig, , drop = FALSE]
    overlap_fam <- unique(d0$cohort_family[paste(sig, d0$cohort_family, sep = "__") %in% overlap$overlap_key])
    scopes <- list(
      secondary_all_locked_external = d0,
      primary_overlap_free_locked_external = d0[!d0$cohort_family %in% overlap_fam, , drop = FALSE]
    )
    if (sig %in% known_training) {
      scopes$sensitivity_known_training_only_overlap_free <- scopes$primary_overlap_free_locked_external
    }
    for (scope_name in names(scopes)) {
      d <- scopes[[scope_name]]
      m <- mlaudit_auc_row(d$true_label, d$score)
      overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
        task = task,
        signature_id = sig,
        benchmark_scope = scope_name,
        n_samples = m$n,
        n_pos = m$n_pos,
        n_families = length(unique(d$cohort_family)),
        n_excluded_overlap_families = length(overlap_fam),
        training_cohort_status = ifelse(sig %in% known_training, "known_training_cohorts", "unknown_training_cohorts"),
        AUROC = m$AUROC,
        AUPRC = m$AUPRC,
        stringsAsFactors = FALSE
      )
    }
  }
}
overlap_bench <- do.call(rbind, overlap_rows)
write_tsv_atomic(overlap_bench, file.path(out_dir, "overlap_aware_signature_benchmark.tsv"))
write_tsv_atomic(overlap_bench, file.path(man_table_dir, "Supplementary_Table_overlap_aware_signature_benchmark.tsv"))
all_scope <- overlap_bench[overlap_bench$benchmark_scope == "secondary_all_locked_external", c("task", "signature_id", "AUROC", "AUPRC")]
primary_scope <- overlap_bench[overlap_bench$benchmark_scope == "primary_overlap_free_locked_external", c("task", "signature_id", "AUROC", "AUPRC", "n_samples", "n_families", "n_excluded_overlap_families", "training_cohort_status")]
infl <- merge(all_scope, primary_scope, by = c("task", "signature_id"), suffixes = c("_all_available", "_overlap_free"), all = TRUE)
infl$overlap_inflation_AUROC <- infl$AUROC_all_available - infl$AUROC_overlap_free
infl$overlap_inflation_AUPRC <- infl$AUPRC_all_available - infl$AUPRC_overlap_free
write_tsv_atomic(infl, file.path(out_dir, "overlap_inflation_index.tsv"))

# 4. Family-cluster uncertainty, leave-family deltas, and equal-family summaries.
selected_groups <- split(ext_pred, paste(ext_pred$task, ext_pred$selected_role, ext_pred$combination_id, sep = "__"))
fam_ci <- do.call(rbind, lapply(selected_groups, function(d) {
  th <- gst_threshold_for(d)
  ci <- gst_group_bootstrap(d, threshold = th, n_boot = config$metrics$bootstrap_n %||% 1000)
  data.frame(task = d$task[1], selected_role = d$selected_role[1], model = d$selector_classifier[1], combination_id = d$combination_id[1], threshold = th, ci, stringsAsFactors = FALSE)
}))
write_tsv_atomic(fam_ci, file.path(out_dir, "family_cluster_bootstrap_ci.tsv"))
write_tsv_atomic(fam_ci, file.path(man_table_dir, "Supplementary_Table_family_cluster_bootstrap_ci.tsv"))

loo <- do.call(rbind, lapply(selected_groups, function(d) {
  full <- compute_binary_metrics(d$true_label, d$score, threshold = gst_threshold_for(d))
  rows <- do.call(rbind, lapply(sort(unique(d$cohort_family)), function(fam) {
    dd <- d[d$cohort_family != fam, , drop = FALSE]
    m <- compute_binary_metrics(dd$true_label, dd$score, threshold = gst_threshold_for(d))
    data.frame(
      removed_family = fam,
      AUROC_without_family = m$AUROC,
      delta_AUROC = m$AUROC - full$AUROC,
      AUPRC_without_family = m$AUPRC,
      delta_AUPRC = m$AUPRC - full$AUPRC,
      sensitivity_without_family = m$sensitivity,
      specificity_without_family = m$specificity,
      stringsAsFactors = FALSE
    )
  }))
  data.frame(task = d$task[1], selected_role = d$selected_role[1], model = d$selector_classifier[1], combination_id = d$combination_id[1], rows, stringsAsFactors = FALSE)
}))
write_tsv_atomic(loo, file.path(out_dir, "leave_one_external_family_out.tsv"))
write_tsv_atomic(loo, file.path(man_table_dir, "Supplementary_Table_leave_one_external_family_out.tsv"))

eq <- do.call(rbind, lapply(selected_groups, function(d) {
  sm <- gst_equal_family_summary(d, threshold = gst_threshold_for(d))
  data.frame(task = d$task[1], selected_role = d$selected_role[1], model = d$selector_classifier[1], combination_id = d$combination_id[1], sm, stringsAsFactors = FALSE)
}))
write_tsv_atomic(eq, file.path(out_dir, "equal_family_weighted_summary.tsv"))
write_tsv_atomic(eq, file.path(man_table_dir, "Supplementary_Table_equal_family_weighted_summary.tsv"))

# 5. Mixed-family-only primary ATB-vs-OD analysis sensitivity summary.
split_balance <- read_tsv("results/tables/development_validation_split_balance.tsv")
mixed_dev <- split_balance$cohort_family[split_balance$split == "development" & split_balance$ATB > 0 & split_balance$OD > 0]
mixed_ext <- split_balance$cohort_family[split_balance$split == "external_validation" & split_balance$ATB > 0 & split_balance$OD > 0]
mixed_rows <- list()
for (i in seq_len(nrow(selected[selected$task == "primary_ATB_vs_OD", , drop = FALSE]))) {
  combo <- selected[selected$task == "primary_ATB_vs_OD", , drop = FALSE][i, ]
  dev_d <- oof[oof$combination_id == combo$combination_id & oof$cohort_family %in% mixed_dev, , drop = FALSE]
  ext_d <- ext_pred[ext_pred$combination_id == combo$combination_id & ext_pred$selected_role == combo$selected_role & ext_pred$cohort_family %in% mixed_ext, , drop = FALSE]
  dm <- compute_binary_metrics(dev_d$true_label, dev_d$score, threshold = 0.5)
  em <- compute_binary_metrics(ext_d$true_label, ext_d$score, threshold = gst_threshold_for(ext_d))
  ef <- summarize_per_study_metrics(ext_d, "cohort_family")
  mixed_rows[[length(mixed_rows) + 1]] <- data.frame(
    task = combo$task,
    selected_role = combo$selected_role,
    model = combo$selector_classifier,
    combination_id = combo$combination_id,
    n_mixed_development_families = length(unique(dev_d$cohort_family)),
    development_mixed_AUROC = dm$AUROC,
    development_mixed_AUPRC = dm$AUPRC,
    n_mixed_external_families = length(unique(ext_d$cohort_family)),
    external_mixed_AUROC = em$AUROC,
    external_mixed_AUPRC = em$AUPRC,
    external_mixed_sensitivity = em$sensitivity,
    external_mixed_specificity = em$specificity,
    external_mixed_family_median_AUROC = stats::median(ef$AUROC, na.rm = TRUE),
    external_mixed_family_iqr_AUROC = paste0(sprintf("%.3f", as.numeric(stats::quantile(ef$AUROC, 0.25, na.rm = TRUE))), "-", sprintf("%.3f", as.numeric(stats::quantile(ef$AUROC, 0.75, na.rm = TRUE)))),
    stringsAsFactors = FALSE
  )
}
mixed_summary <- do.call(rbind, mixed_rows)
write_tsv_atomic(mixed_summary, file.path(out_dir, "mixed_family_only_primary_ATB_vs_OD_summary.tsv"))
write_tsv_atomic(mixed_summary, file.path(man_table_dir, "Supplementary_Table_mixed_family_only_primary_ATB_vs_OD_summary.tsv"))

# 6. Model-shopping audit with development-vs-external scatter-ready data.
shopping_rows <- list()
for (short in c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB")) {
  p <- file.path(out_dir, sprintf("ensemble_external_base_predictions_%s.tsv.gz", short))
  if (!file.exists(project_path(p))) next
  base <- read_tsv(p)
  for (combo_id in sort(unique(base$combination_id))) {
    d <- base[base$combination_id == combo_id, , drop = FALSE]
    m <- compute_binary_metrics(d$true_label, d$score)
    dp <- dev_perf[dev_perf$combination_id == combo_id, , drop = FALSE]
    shopping_rows[[length(shopping_rows) + 1]] <- data.frame(
      task = d$task[1],
      combination_id = combo_id,
      selector_classifier = d$selector_classifier[1],
      feature_source = d$selector[1],
      classifier = d$classifier[1],
      model_complexity = d$n_model_features[1] %||% NA_real_,
      development_AUROC = dp$AUROC[1] %||% NA_real_,
      development_AUPRC = dp$AUPRC[1] %||% NA_real_,
      external_AUROC = m$AUROC,
      external_AUPRC = m$AUPRC,
      external_Brier = m$Brier,
      external_calibration_slope = calibration_summary(d$true_label, d$score)$calibration_slope,
      external_sensitivity_at_0_50 = m$sensitivity,
      external_specificity_at_0_50 = m$specificity,
      selected_final_model = combo_id %in% selected$combination_id,
      stringsAsFactors = FALSE
    )
  }
}
if (length(shopping_rows) > 0) {
  shopping <- do.call(rbind, shopping_rows)
  shopping <- shopping[order(shopping$task, -shopping$development_AUROC, -shopping$development_AUPRC), , drop = FALSE]
  shopping$model_rank_by_development <- ave(-shopping$development_AUROC, shopping$task, FUN = function(x) rank(x, ties.method = "first"))
  write_tsv_atomic(shopping, file.path(out_dir, "development_vs_external_model_selection_check.tsv"))
  write_tsv_atomic(shopping, file.path(man_table_dir, "Supplementary_Table_development_vs_external_model_selection_check.tsv"))
}

# 7. Expanded development-only calibration repair audit.
cal_rows <- list()
for (i in seq_len(nrow(selected))) {
  combo <- selected[i, , drop = FALSE]
  dev <- oof[oof$combination_id == combo$combination_id, , drop = FALSE]
  ext <- ext_pred[ext_pred$combination_id == combo$combination_id & ext_pred$selected_role == combo$selected_role, , drop = FALSE]
  if (nrow(dev) == 0 || nrow(ext) == 0) next
  method_scores <- list(raw_uncalibrated = ext$score)
  platt <- mlaudit_fit_platt(dev$true_label, dev$score)
  method_scores$development_oof_platt <- mlaudit_predict_platt(platt, ext$score)
  iso <- gst_fit_isotonic(dev$true_label, dev$score)
  method_scores$development_oof_isotonic <- gst_predict_isotonic(iso, ext$score)
  temp <- gst_fit_temperature(dev$true_label, dev$score)
  method_scores$development_oof_temperature <- gst_predict_temperature(temp, ext$score)
  beta <- gst_fit_beta(dev$true_label, dev$score)
  method_scores$development_oof_beta <- gst_predict_beta(beta, ext$score)
  for (method in names(method_scores)) {
    cm <- gst_calibration_metrics(ext$true_label, method_scores[[method]], threshold = 0.5)
    cal_rows[[length(cal_rows) + 1]] <- data.frame(
      task = combo$task,
      selected_role = combo$selected_role,
      model = combo$selector_classifier,
      combination_id = combo$combination_id,
      calibration_method = method,
      cm,
      calibration_fit_note = ifelse(method == "raw_uncalibrated", "none", "fit_on_development_oof_only"),
      stringsAsFactors = FALSE
    )
  }
}
cal_repair <- do.call(rbind, cal_rows)
write_tsv_atomic(cal_repair, file.path(out_dir, "calibration_repair_external_summary.tsv"))
write_tsv_atomic(cal_repair, file.path(man_table_dir, "Supplementary_Table_calibration_repair_external_summary.tsv"))

# 8. Selective prediction / abstention using development OOF thresholds.
sel_grid_rows <- list()
sel_point_rows <- list()
strategy_names <- c(
  "dev_min_error_at_coverage_0.70",
  "dev_specificity_ge_0.70_max_coverage",
  "dev_sensitivity_ge_0.90_max_coverage",
  "dev_balanced_accuracy_max"
)
for (i in seq_len(nrow(selected))) {
  combo <- selected[i, , drop = FALSE]
  dev <- oof[oof$combination_id == combo$combination_id, , drop = FALSE]
  ext <- ext_pred[ext_pred$combination_id == combo$combination_id & ext_pred$selected_role == combo$selected_role, , drop = FALSE]
  if (nrow(dev) == 0 || nrow(ext) == 0) next
  qs_low <- unique(as.numeric(stats::quantile(dev$score, probs = seq(0.05, 0.45, 0.05), na.rm = TRUE)))
  qs_high <- unique(as.numeric(stats::quantile(dev$score, probs = seq(0.55, 0.95, 0.05), na.rm = TRUE)))
  grid <- expand.grid(lower_threshold = qs_low, upper_threshold = qs_high)
  grid <- grid[grid$lower_threshold < grid$upper_threshold, , drop = FALSE]
  if (nrow(grid) == 0) next
  met <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    dm <- gst_selective_metrics(dev$true_label, dev$score, grid$lower_threshold[j], grid$upper_threshold[j])
    em <- gst_selective_metrics(ext$true_label, ext$score, grid$lower_threshold[j], grid$upper_threshold[j])
    data.frame(
      lower_threshold = grid$lower_threshold[j],
      upper_threshold = grid$upper_threshold[j],
      dev_coverage = dm$coverage,
      dev_abstention_rate = dm$abstention_rate,
      dev_error_rate_among_answered = dm$error_rate_among_answered,
      dev_sensitivity_answered = dm$sensitivity_answered,
      dev_specificity_answered = dm$specificity_answered,
      dev_balanced_accuracy_answered = mean(c(dm$sensitivity_answered, dm$specificity_answered), na.rm = TRUE),
      external_coverage = em$coverage,
      external_abstention_rate = em$abstention_rate,
      external_error_rate_among_answered = em$error_rate_among_answered,
      external_sensitivity_answered = em$sensitivity_answered,
      external_specificity_answered = em$specificity_answered,
      external_false_negative_rate_answered = em$false_negative_rate_answered,
      external_false_positive_rate_answered = em$false_positive_rate_answered,
      external_PPV_answered = em$PPV_answered,
      external_NPV_answered = em$NPV_answered,
      n = em$n,
      answered_n = em$answered_n,
      uncertain_n = em$uncertain_n,
      external_tp_answered = em$tp_answered,
      external_tn_answered = em$tn_answered,
      external_fp_answered = em$fp_answered,
      external_fn_answered = em$fn_answered,
      stringsAsFactors = FALSE
    )
  }))
  met$key <- paste(sprintf("%.6f", met$lower_threshold), sprintf("%.6f", met$upper_threshold), sep = "__")
  met$task <- combo$task
  met$selected_role <- combo$selected_role
  met$model <- combo$selector_classifier
  met$combination_id <- combo$combination_id
  met$is_selected_strategy <- FALSE
  met$selected_strategy <- NA_character_
  for (strategy in strategy_names) {
    idx <- gst_pick_strategy(met, strategy)
    if (length(idx) == 1 && is.finite(idx)) {
      row <- met[idx, , drop = FALSE]
      row$is_selected_strategy <- TRUE
      row$selected_strategy <- strategy
      sel_point_rows[[length(sel_point_rows) + 1]] <- row
      met$is_selected_strategy[idx] <- TRUE
      met$selected_strategy[idx] <- paste(na.omit(c(met$selected_strategy[idx], strategy)), collapse = ";")
    }
  }
  sel_grid_rows[[length(sel_grid_rows) + 1]] <- met
}
selective_grid <- do.call(rbind, sel_grid_rows)
selective_points <- do.call(rbind, sel_point_rows)
write_tsv_atomic(selective_grid, file.path(out_dir, "selective_prediction_risk_coverage.tsv"))
write_tsv_atomic(selective_points, file.path(out_dir, "selective_prediction_operating_points.tsv"))
write_tsv_atomic(selective_points, file.path(man_table_dir, "Supplementary_Table_selective_prediction_operating_points.tsv"))

# 9. Negative controls: permutation, random score, random signature-id score, platform, family oracle.
n_perm <- config$negative_controls$n_permutations %||% 100
perm_rows <- do.call(rbind, lapply(selected_groups, function(d) {
  p <- mlaudit_score_label_permutation(d, n_perm = n_perm, seed = 20260502)
  data.frame(task = d$task[1], selected_role = d$selected_role[1], model = d$selector_classifier[1], combination_id = d$combination_id[1], negative_control = "selected_external_score_label_permutation", p, stringsAsFactors = FALSE)
}))
write_tsv_atomic(perm_rows, file.path(out_dir, "negative_control_score_label_permutation.tsv"))

within_family_perm <- do.call(rbind, lapply(seq_along(selected_groups), function(i) {
  d <- selected_groups[[i]]
  p <- gst_within_family_score_permutation(d, n_perm = n_perm, seed = 20260502 + i)
  data.frame(task = d$task[1], selected_role = d$selected_role[1], model = d$selector_classifier[1], combination_id = d$combination_id[1], p, stringsAsFactors = FALSE)
}))
write_tsv_atomic(within_family_perm, file.path(out_dir, "negative_control_family_stratified_permutation.tsv"))
write_tsv_atomic(within_family_perm, file.path(man_table_dir, "Supplementary_Table_negative_control_family_stratified_permutation.tsv"))

set.seed(20260502)
rand_rows <- list()
for (task in tasks) {
  task_ext <- ext_pred[ext_pred$task == task & ext_pred$selected_role == "best_performance", , drop = FALSE]
  if (nrow(task_ext) > 0) {
    vals <- replicate(n_perm, compute_roc_auc(task_ext$true_label, stats::runif(nrow(task_ext))))
    rand_rows[[length(rand_rows) + 1]] <- data.frame(
      task = task,
      negative_control = "random_uniform_score",
      n_iterations = n_perm,
      AUROC_median = stats::median(vals, na.rm = TRUE),
      AUROC_q025 = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
      AUROC_q975 = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  task_sig <- sig_scores[sig_scores$task == task & sig_scores$cohort_family %in% external_families, , drop = FALSE]
  sig_ids <- sort(unique(task_sig$signature_id))
  if (length(sig_ids) > 0) {
    sampled <- sample(sig_ids, n_perm, replace = length(sig_ids) < n_perm)
    vals <- vapply(sampled, function(sig) {
      d <- task_sig[task_sig$signature_id == sig, , drop = FALSE]
      compute_roc_auc(d$true_label, d$score)
    }, numeric(1))
    rand_rows[[length(rand_rows) + 1]] <- data.frame(
      task = task,
      negative_control = "random_available_signature_id_score",
      n_iterations = length(vals),
      AUROC_median = stats::median(vals, na.rm = TRUE),
      AUROC_q025 = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
      AUROC_q975 = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
}
random_control <- do.call(rbind, rand_rows)
write_tsv_atomic(random_control, file.path(out_dir, "negative_control_random_signature.tsv"))
write_tsv_atomic(random_control, file.path(man_table_dir, "Supplementary_Table_negative_control_random_signature.tsv"))

family_oracle_rows <- list()
for (task in tasks) {
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  meta$split <- ifelse(meta$cohort_family %in% external_families, "external_validation",
    ifelse(meta$cohort_family %in% development_families, "development", "not_in_lock"))
  for (sp in c("development", "external_validation")) {
    d <- meta[meta$split == sp, , drop = FALSE]
    if (nrow(d) == 0) next
    fam_prev <- ave(as_binary(d$task_label), d$cohort_family, FUN = mean)
    family_oracle_rows[[length(family_oracle_rows) + 1]] <- data.frame(
      task = task,
      split = sp,
      n = nrow(d),
      n_families = length(unique(d$cohort_family)),
      family_oracle_AUROC = compute_roc_auc(d$task_label, fam_prev),
      family_label_cramers_v = mlaudit_cramers_v(d$cohort_family, d$task_label),
      interpretation = "oracle_uses_family_labels_within_split_to_expose_label_family_confounding_risk_not_a_deployable_model",
      stringsAsFactors = FALSE
    )
  }
}
family_oracle <- do.call(rbind, family_oracle_rows)
write_tsv_atomic(family_oracle, file.path(out_dir, "family_label_confounding_check.tsv"))
write_tsv_atomic(family_oracle, file.path(man_table_dir, "Supplementary_Table_family_label_confounding_check.tsv"))

# 9b. Ledger-driven integrated strengthening tables, excluding optional biology Module 11.
family_rules <- if (file.exists(project_path("results/tables/uncertain_family_assignment.tsv"))) {
  read_tsv("results/tables/uncertain_family_assignment.tsv")
} else {
  data.frame(study = character(), cohort_family = character(), family_assignment_rule = character(), notes = character(), stringsAsFactors = FALSE)
}
master <- read_tsv("data_meta/sample_harmonization_master.tsv")
crosswalk <- merge(
  master[, intersect(c("standard_sample_id", "study", "cohort_family", "platform", "platform_type", "geographical_region"), names(master)), drop = FALSE],
  family_rules[, intersect(c("study", "cohort_family", "family_assignment_rule", "publication_key", "confidence", "notes"), names(family_rules)), drop = FALSE],
  by = c("study", "cohort_family"),
  all.x = TRUE
)
crosswalk$family_assignment_rule[is.na(crosswalk$family_assignment_rule)] <- "study_as_family_default"
crosswalk$publication_key[is.na(crosswalk$publication_key)] <- crosswalk$study[is.na(crosswalk$publication_key)]
crosswalk$confidence[is.na(crosswalk$confidence)] <- "medium"
crosswalk$notes[is.na(crosswalk$notes)] <- ""
crosswalk$family_geo_accession_level <- crosswalk$study
crosswalk$family_parent_study_level_main <- crosswalk$cohort_family
crosswalk$family_conservative_exclusion <- ifelse(crosswalk$family_assignment_rule == "study_as_family_default", crosswalk$cohort_family, NA_character_)
crosswalk$family_maximally_conservative_merge <- crosswalk$cohort_family
crosswalk$ambiguity_flag <- crosswalk$family_assignment_rule != "study_as_family_default"
crosswalk$merge_reason <- ifelse(crosswalk$ambiguity_flag, crosswalk$family_assignment_rule, "study_as_family_default")
crosswalk_out <- unique(crosswalk[, c(
  "standard_sample_id", "study", "cohort_family", "family_geo_accession_level",
  "family_parent_study_level_main", "family_conservative_exclusion",
  "family_maximally_conservative_merge", "publication_key", "merge_reason",
  "ambiguity_flag", "confidence", "notes"
), drop = FALSE])
write_tsv_atomic(crosswalk_out, file.path(out_dir, "family_definition_crosswalk.tsv"))
write_tsv_atomic(crosswalk_out, file.path(man_table_dir, "Supplementary_Table_family_definition_crosswalk.tsv"))

family_defs <- data.frame(
  family_definition = c("geo_accession_level", "parent_study_level_main", "conservative_exclusion", "maximally_conservative_merge"),
  family_col = c("family_geo_accession_level", "family_parent_study_level_main", "family_conservative_exclusion", "family_maximally_conservative_merge"),
  interpretation = c(
    "each source accession is treated as its own family",
    "current main analysis family definition",
    "non-default merged accessions are excluded from grouping-sensitive summaries",
    "current related-accession merges retained as the maximal encoded conservative merge"
  ),
  stringsAsFactors = FALSE
)

oracle_def_rows <- list()
phenotype_rows <- list()
for (task in tasks) {
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  meta <- gst_map_family_definitions(meta, crosswalk_out)
  meta$split <- ifelse(meta$cohort_family %in% external_families, "external_validation",
    ifelse(meta$cohort_family %in% development_families, "development", "not_in_lock"))
  for (sp in c("development", "external_validation")) {
    d <- meta[meta$split == sp, , drop = FALSE]
    if (nrow(d) == 0) next
    for (i in seq_len(nrow(family_defs))) {
      oracle_def_rows[[length(oracle_def_rows) + 1]] <- gst_family_oracle_row(d, family_defs$family_col[i], task, sp, family_defs$family_definition[i])
      phenotype_rows[[length(phenotype_rows) + 1]] <- gst_family_phenotype_rows(d, family_defs$family_col[i], task, sp, family_defs$family_definition[i])
    }
  }
}
family_definition_oracle <- do.call(rbind, oracle_def_rows)
family_phenotype_structure <- do.call(rbind, phenotype_rows)
write_tsv_atomic(family_definition_oracle, file.path(out_dir, "family_definition_oracle_stability.tsv"))
write_tsv_atomic(family_definition_oracle, file.path(man_table_dir, "Supplementary_Table_family_definition_oracle_stability.tsv"))
write_tsv_atomic(family_phenotype_structure, file.path(out_dir, "family_phenotype_structure.tsv"))
write_tsv_atomic(family_phenotype_structure, file.path(man_table_dir, "Supplementary_Table_family_phenotype_structure.tsv"))

family_def_perm_rows <- list()
for (g in seq_along(selected_groups)) {
  d0 <- gst_map_family_definitions(selected_groups[[g]], crosswalk_out)
  for (i in seq_len(nrow(family_defs))) {
    p <- gst_within_family_perm_for_col(d0, family_defs$family_col[i], n_perm = n_perm, seed = 20260530 + 100L * g + i)
    family_def_perm_rows[[length(family_def_perm_rows) + 1]] <- data.frame(
      task = d0$task[1],
      selected_role = d0$selected_role[1],
      model = d0$selector_classifier[1],
      combination_id = d0$combination_id[1],
      family_definition = family_defs$family_definition[i],
      p,
      analysis_scope = "post_hoc_regrouping_no_model_refit",
      stringsAsFactors = FALSE
    )
  }
}
family_definition_perm <- do.call(rbind, family_def_perm_rows)
write_tsv_atomic(family_definition_perm, file.path(out_dir, "family_definition_permutation_stability.tsv"))
write_tsv_atomic(family_definition_perm, file.path(man_table_dir, "Supplementary_Table_family_definition_permutation_stability.tsv"))

reference_h2h_for_stability <- if (file.exists(project_path(file.path(out_dir, "reference_signature_locked_external_head_to_head.tsv")))) {
  read_tsv(file.path(out_dir, "reference_signature_locked_external_head_to_head.tsv"))
} else if (file.exists(project_path(file.path(out_dir, "published_signature_locked_external_head_to_head.tsv")))) {
  read_tsv(file.path(out_dir, "published_signature_locked_external_head_to_head.tsv"))
} else {
  data.frame()
}
stability_rows <- list()
for (i in seq_len(nrow(selected))) {
  task <- selected$task[i]
  role <- selected$selected_role[i]
  combo <- selected$combination_id[i]
  ext_hit <- ext_summary[ext_summary$task == task & ext_summary$selected_role == role & ext_summary$combination_id == combo, , drop = FALSE]
  refs <- reference_h2h_for_stability[reference_h2h_for_stability$task == task & !grepl("^final_", reference_h2h_for_stability$model_layer), , drop = FALSE]
  best_ref_auc <- if (nrow(refs) > 0) max(refs$AUROC, na.rm = TRUE) else NA_real_
  best_ref <- if (nrow(refs) > 0) refs$signature_or_model[which.max(refs$AUROC)] else NA_character_
  for (fd in family_defs$family_definition) {
    dev_o <- family_definition_oracle[family_definition_oracle$task == task & family_definition_oracle$split == "development" & family_definition_oracle$family_definition == fd, , drop = FALSE]
    ext_o <- family_definition_oracle[family_definition_oracle$task == task & family_definition_oracle$split == "external_validation" & family_definition_oracle$family_definition == fd, , drop = FALSE]
    perm <- family_definition_perm[family_definition_perm$task == task & family_definition_perm$selected_role == role & family_definition_perm$combination_id == combo & family_definition_perm$family_definition == fd, , drop = FALSE]
    stability_rows[[length(stability_rows) + 1]] <- data.frame(
      task = task,
      selected_role = role,
      model = selected$selector_classifier[i],
      combination_id = combo,
      family_definition = fd,
      development_n_families = if (nrow(dev_o) > 0) dev_o$n_families[1] else NA_integer_,
      development_n_mixed_families = if (nrow(dev_o) > 0) dev_o$n_mixed_families[1] else NA_integer_,
      development_family_oracle_AUROC = if (nrow(dev_o) > 0) dev_o$family_oracle_AUROC[1] else NA_real_,
      external_n_families = if (nrow(ext_o) > 0) ext_o$n_families[1] else NA_integer_,
      external_n_mixed_families = if (nrow(ext_o) > 0) ext_o$n_mixed_families[1] else NA_integer_,
      external_family_oracle_AUROC = if (nrow(ext_o) > 0) ext_o$family_oracle_AUROC[1] else NA_real_,
      observed_external_AUROC = if (nrow(ext_hit) > 0) ext_hit$AUROC[1] else NA_real_,
      best_reference_signature = best_ref,
      best_reference_AUROC = best_ref_auc,
      delta_ML_vs_best_reference_AUROC = if (nrow(ext_hit) > 0 && is.finite(best_ref_auc)) ext_hit$AUROC[1] - best_ref_auc else NA_real_,
      within_family_permutation_AUROC_median = if (nrow(perm) > 0) perm$AUROC_median[1] else NA_real_,
      within_family_permutation_AUROC_q025 = if (nrow(perm) > 0) perm$AUROC_q025[1] else NA_real_,
      within_family_permutation_AUROC_q975 = if (nrow(perm) > 0) perm$AUROC_q975[1] else NA_real_,
      analysis_scope = "post_hoc_regrouping_no_model_refit",
      stability_note = "supports sensitivity of family-structure summaries; not a full model-refit family-definition rerun",
      stringsAsFactors = FALSE
    )
  }
}
family_definition_stability <- do.call(rbind, stability_rows)
write_tsv_atomic(family_definition_stability, file.path(out_dir, "family_definition_stability_matrix.tsv"))
write_tsv_atomic(family_definition_stability, file.path(man_table_dir, "Supplementary_Table_family_definition_stability_matrix.tsv"))

split_repeats <- config$split_sensitivity_repeats %||% 500
split_lite <- lapply(seq_len(nrow(selected)), function(i) {
  gst_lite_split_sensitivity(selected$task[i], selected$selected_role[i], n_repeats = split_repeats, seed = 20260530 + i)
})
split_lite_results <- do.call(rbind, lapply(split_lite, `[[`, "results"))
split_lite_summary <- do.call(rbind, lapply(split_lite, `[[`, "summary"))
write_tsv_atomic(split_lite_results, file.path(out_dir, "family_block_split_lite_results.tsv"))
write_tsv_atomic(split_lite_summary, file.path(out_dir, "family_block_split_lite_summary.tsv"))
write_tsv_atomic(split_lite_summary, file.path(man_table_dir, "Supplementary_Table_family_block_split_lite_summary.tsv"))

family_calibration <- do.call(rbind, lapply(selected_groups, gst_family_calibration_threshold_audit))
write_tsv_atomic(family_calibration, file.path(out_dir, "family_specific_calibration_thresholds.tsv"))
write_tsv_atomic(family_calibration, file.path(man_table_dir, "Supplementary_Table_family_specific_calibration_thresholds.tsv"))

artifact_status <- if (file.exists(project_path(file.path(out_dir, "model_external_artifact_status.tsv")))) read_tsv(file.path(out_dir, "model_external_artifact_status.tsv")) else data.frame()
dev_ext_check <- if (file.exists(project_path(file.path(out_dir, "development_vs_external_model_selection_check.tsv")))) read_tsv(file.path(out_dir, "development_vs_external_model_selection_check.tsv")) else data.frame()
lock_hash_now <- hash_file_or_object("config/external_validation_lockfile.json")
selected_external_ok <- all(selected$combination_id %in% unique(ext_pred$combination_id))
artifact_bool <- function(x) x %in% TRUE | toupper(as.character(x)) == "TRUE"
artifact_selected <- if (nrow(artifact_status) > 0) artifact_bool(artifact_status$selected_for_external) else logical()
artifact_available <- if (nrow(artifact_status) > 0) artifact_bool(artifact_status$external_prediction_available) else logical()
nonselected_no_artifact_ok <- if (nrow(artifact_status) > 0) {
  idx <- !artifact_selected & !artifact_available
  all(artifact_status$reason_no_external[idx] == "not_selected_for_locked_external_prediction", na.rm = TRUE)
} else {
  FALSE
}
nonselected_artifact_labeled_ok <- if (nrow(artifact_status) > 0) {
  idx <- !artifact_selected & artifact_available
  if (!any(idx)) TRUE else all(artifact_status$reason_no_external[idx] == "external_artifact_available", na.rm = TRUE)
} else {
  FALSE
}
no_leakage_audit <- data.frame(
  audit_item = c(
    "external_lockfile_unchanged_during_stage",
    "external_artifact_status_covers_full_grid",
    "selected_models_have_external_predictions",
    "nonselected_external_artifact_scope_is_labeled",
    "development_external_selection_scope_is_explicit",
    "thresholds_are_development_locked",
    "no_leakage_grouped_split_tests_present"
  ),
  status = c(
    ifelse(identical(lock_hash_before, lock_hash_now), "pass", "fail"),
    ifelse(nrow(artifact_status) == 952, "pass", "fail"),
    ifelse(selected_external_ok, "pass", "fail"),
    ifelse(nonselected_no_artifact_ok && nonselected_artifact_labeled_ok, "pass", "fail"),
    ifelse(nrow(dev_ext_check) > 0 && all(dev_ext_check$selected_final_model %in% c(TRUE, FALSE)), "pass", "fail"),
    ifelse(nrow(thresholds) > 0 && all(grepl("^development_|^fixed_", thresholds$threshold_name)), "pass", "fail"),
    ifelse(all(file.exists(project_path(c("tests/testthat/test_01_no_leakage.R", "tests/testthat/test_02_splits_are_grouped.R", "tests/testthat/test_03_thresholds_locked.R")))), "pass", "fail")
  ),
  evidence_file = c(
    "config/external_validation_lockfile.json",
    "results/ml_final/model_external_artifact_status.tsv",
    "results/ml_final/model_external_artifact_status.tsv",
    "results/ml_final/model_external_artifact_status.tsv",
    "results/ml_final/development_vs_external_model_selection_check.tsv",
    "results/ml_final/development_locked_thresholds.tsv",
    "06_validation_tests/01_testthat/test_01_no_leakage.R;06_validation_tests/01_testthat/test_02_splits_are_grouped.R;06_validation_tests/01_testthat/test_03_thresholds_locked.R"
  ),
  interpretation = c(
    "external lockfile hash is stable within the current analysis run",
    "all 952 development grid combinations are represented in the external-evaluation scope check",
    "selected final models have locked external prediction rows in external_predictions.tsv.gz",
    "nonselected rows are either explicitly not external-tested or labeled as top-development check outputs, not selected final models",
    "development and external metrics are displayed as a check, not as external model selection",
    "threshold names and threshold table document development-locked or fixed thresholds",
    "automated tests cover no leakage, grouped splits, and locked thresholds"
  ),
  stringsAsFactors = FALSE
)
write_tsv_atomic(no_leakage_audit, file.path(out_dir, "data_separation_checks.tsv"))
write_tsv_atomic(no_leakage_audit, file.path(man_table_dir, "Supplementary_Table_data_separation_checks.tsv"))

claim_evidence_map <- data.frame(
  claim_id = c(
    "C01_family_identity_carries_phenotype_information",
    "C02_limited_incremental_gain_over_reference_signatures",
    "C03_primary_mixed_family_support_is_sparse",
    "C04_family_definition_sensitivity_is_bounded",
    "C05_primary_calibration_and_threshold_transport_are_unstable",
    "C06_permutation_is_stress_test_not_decomposition",
    "C07_external_no_leakage_scope",
    "C08_selective_prediction_does_not_solve_deployment",
    "C09_reference_signature_benchmark_is_standardized_not_exact_replication"
  ),
  manuscript_section = c("Results", "Abstract/Results", "Results", "Supplementary Results", "Results/Discussion", "Methods/Results", "Methods", "Results/Discussion", "Methods/Limitations"),
  risk_level = c("low", "medium", "low", "medium", "low", "medium", "low", "medium", "medium"),
  supporting_file = c(
    "Supplementary_Table_family_label_confounding_check.tsv;Supplementary_Table_family_phenotype_structure.tsv",
    "Table2_reference_signature_locked_external_head_to_head.tsv;Supplementary_Table_family_block_split_lite_summary.tsv",
    "Supplementary_Table_mixed_family_only_primary_ATB_vs_OD_summary.tsv;Supplementary_Table_family_phenotype_structure.tsv",
    "Supplementary_Table_family_definition_stability_matrix.tsv",
    "Supplementary_Table_family_specific_calibration_thresholds.tsv;Supplementary_Table_calibrated_operating_points_external.tsv",
    "Supplementary_Table_negative_control_family_stratified_permutation.tsv",
    "Supplementary_Table_data_separation_checks.tsv;Supplementary_Table_model_external_artifact_status.tsv",
    "Supplementary_Table_selective_prediction_family_level_coverage.tsv;Table4_operating_points_selective_prediction.tsv",
    "Supplementary_Table_reference_signature_implementation_and_gene_coverage.tsv;Supplementary_Table_signature_coverage_sensitivity.tsv"
  ),
  revised_claim = c(
    "Family labels carry phenotype information and therefore sample-level validation can overstate generalization.",
    "Optimized ML did not show a stable advantage over the strongest reference-signature comparators.",
    "The primary ATB-versus-OD development contrast is supported by few mixed families.",
    "Post-hoc regrouping of encoded family definitions supports the family-structure interpretation, but it is not a full refit sensitivity.",
    "Primary probability calibration and family-level threshold transfer remain unstable despite development-only repair checks.",
    "Within-family score permutation is a structural stress test and not an additive AUROC decomposition.",
    "The 952-combination grid is a development search; locked external artifacts are generated only for selected/top comparators.",
    "Selective prediction reduces coverage and does not fully rescue false-negative or family-level operating risk.",
    "Reference signatures are standardized gene-score benchmarks because original coefficients are not always available."
  ),
  status = "supported_or_softened",
  stringsAsFactors = FALSE
)
write_tsv_atomic(claim_evidence_map, file.path(out_dir, "statement_support_summary.tsv"))

module_coverage <- data.frame(
  module = c("Module 01", "Module 02", "Module 03 Lite", "Module 04", "Module 05", "Module 06", "Module 07", "Module 08", "Module 09", "Module 10", "Module 12"),
  integrated_status = c("implemented", "implemented_scope_limited", "implemented_lite", "implemented", "implemented", "implemented", "implemented", "implemented", "implemented", "implemented", "implemented"),
  primary_outputs = c(
    "negative_control_summary.tsv; check_submission_consistency.R",
    "family_definition_crosswalk.tsv; family_definition_stability_matrix.tsv",
    "family_block_split_lite_results.tsv; family_block_split_lite_summary.tsv",
    "family_phenotype_structure.tsv; family_label_confounding_check.tsv",
    "external_per_family_performance.tsv; equal_family_weighted_summary.tsv; leave_one_external_family_out.tsv",
    "reference_signature_implementation_and_gene_coverage.tsv; signature_coverage_sensitivity.tsv; overlap_aware_signature_benchmark.tsv",
    "negative_control_summary.tsv; negative_control_score_label_permutation.tsv; negative_control_family_stratified_permutation.tsv; negative_control_random_signature.tsv",
    "family_specific_calibration_thresholds.tsv; calibration_repair_external_summary.tsv; calibrated_operating_points_external.tsv",
    "selective_prediction_risk_coverage.tsv; selective_prediction_family_level_coverage.tsv",
    "data_separation_checks.tsv; model_external_artifact_status.tsv; external_lockfile_integrity.tsv",
    "statement_support_summary.tsv; manuscript_wording_checks"
  ),
  manuscript_outputs = c(
    "Supplementary_Table_negative_control_summary.tsv",
    "Supplementary_Table_family_definition_stability_matrix.tsv",
    "Supplementary_Table_family_block_split_lite_summary.tsv",
    "Supplementary_Table_family_phenotype_structure.tsv",
    "Supplementary_Table_equal_family_weighted_summary.tsv",
    "Supplementary_Table_reference_signature_implementation_and_gene_coverage.tsv",
    "Supplementary_Table_negative_control_family_stratified_permutation.tsv",
    "Supplementary_Table_family_specific_calibration_thresholds.tsv",
    "Supplementary_Table_selective_prediction_family_level_coverage.tsv",
    "Supplementary_Table_data_separation_checks.tsv",
    "source-data only"
  ),
  limitation = c(
    "Package consistency is checked by scripts rather than a separate narrative audit directory.",
    "Family definitions are evaluated by post-hoc regrouping; this is not a full model-refit F0-F3 rerun.",
    "Lite uses development OOF family-block rotations with fixed candidates; it is not a full 952-grid rerun per split.",
    "Phenotype entropy and mixed-family scarcity support structure, not mechanism.",
    "Family-level AUROC is only estimable when both classes occur in a family.",
    "Reference signatures are standardized benchmarks, not exact original-model replications.",
    "Random reference-signature spread is a benchmark spread, not a null random gene-set distribution.",
    "Family-specific calibration slopes may be non-estimable in small or separated families.",
    "Selective prediction is an operating-point sensitivity, not a deployment solution.",
    "Data-separation checks are machine-readable and do not replace independent code review.",
    "Statement-support notes soften wording without adding biological mechanism claims."
  ),
  stringsAsFactors = FALSE
)
write_tsv_atomic(module_coverage, file.path(out_dir, "supporting_analysis_coverage.tsv"))

gst_update_table_index(data.frame(
  table_file = c(
    "Supplementary_Table_family_definition_crosswalk.tsv",
    "Supplementary_Table_family_definition_oracle_stability.tsv",
    "Supplementary_Table_family_definition_permutation_stability.tsv",
    "Supplementary_Table_family_definition_stability_matrix.tsv",
    "Supplementary_Table_family_phenotype_structure.tsv",
    "Supplementary_Table_family_block_split_lite_summary.tsv",
    "Supplementary_Table_family_specific_calibration_thresholds.tsv",
    "Supplementary_Table_data_separation_checks.tsv"
  ),
  placement = c(
    "source_data",
    "supplementary",
    "supplementary",
    "supplementary",
    "supplementary",
    "supplementary",
    "supplementary",
    "source_data"
  ),
  title = c(
    "Family-definition crosswalk",
    "Family-definition oracle stability",
    "Family-definition permutation stability",
    "Family-definition stability matrix",
    "Family phenotype structure",
    "Family-block split sensitivity Lite",
    "Family-specific calibration and threshold transport",
    "Data-separation checks"
  ),
  purpose = c(
    "Maps samples from accession-level to encoded parent-study family definitions.",
    "Quantifies family-label oracle behavior under encoded family definitions.",
    "Checks within-family score-permutation behavior under encoded family definitions.",
    "Summarizes family-definition sensitivity for selected models and reference comparators.",
    "Quantifies family phenotype mixing, entropy, and low-volume flags.",
    "Summarizes development OOF family-block rotations with fixed candidate models.",
    "Reports family-specific calibration, threshold behavior, and collapse flags.",
    "Compiles machine-readable data-separation and reproducibility checks."
  ),
  source_table = c(
    "data_meta/sample_harmonization_master.tsv; results/tables/uncertain_family_assignment.tsv",
    "data_curated/task_matrices/*/metadata.tsv",
    "results/ml_final/external_predictions.tsv.gz",
    "family_definition_oracle_stability.tsv; family_definition_permutation_stability.tsv",
    "data_curated/task_matrices/*/metadata.tsv",
    "development_oof_predictions.tsv.gz; signature_scores_by_sample.tsv.gz",
    "external_predictions.tsv.gz",
    "model_external_artifact_status.tsv; 06_validation_tests/01_testthat"
  ),
  notes = c(
    "Source-data only.",
    "Included in Supplementary Table S2.",
    "Included in Supplementary Table S2.",
    "Included in Supplementary Table S2.",
    "Included in Supplementary Table S2.",
    "Included in Supplementary Table S5.",
    "Included in Supplementary Table S7.",
    "Included in Supplementary Table S4 as the data_separation_checks sheet."
  ),
  exists_in_package = TRUE,
  stringsAsFactors = FALSE
))

platform_check <- if (file.exists(project_path(file.path(out_dir, "platform_label_confounding_check.tsv")))) {
  read_tsv(file.path(out_dir, "platform_label_confounding_check.tsv"))
} else {
  data.frame()
}
negative_summary <- rbind(
  data.frame(
    negative_control = "selected_external_score_label_permutation",
    n_iterations = n_perm,
    AUROC_median = stats::median(perm_rows$permuted_AUROC_median, na.rm = TRUE),
    AUROC_q025 = stats::median(perm_rows$permuted_AUROC_q025, na.rm = TRUE),
    AUROC_q975 = stats::median(perm_rows$permuted_AUROC_q975, na.rm = TRUE),
    interpretation = "score-label permutation centered near chance",
    stringsAsFactors = FALSE
  ),
  data.frame(
    negative_control = "random_uniform_score",
    n_iterations = n_perm,
    AUROC_median = stats::median(random_control$AUROC_median[random_control$negative_control == "random_uniform_score"], na.rm = TRUE),
    AUROC_q025 = stats::median(random_control$AUROC_q025[random_control$negative_control == "random_uniform_score"], na.rm = TRUE),
    AUROC_q975 = stats::median(random_control$AUROC_q975[random_control$negative_control == "random_uniform_score"], na.rm = TRUE),
    interpretation = "random score centered near chance",
    stringsAsFactors = FALSE
  ),
  data.frame(
    negative_control = "random_available_signature_id_score",
    n_iterations = n_perm,
    AUROC_median = stats::median(random_control$AUROC_median[random_control$negative_control == "random_available_signature_id_score"], na.rm = TRUE),
    AUROC_q025 = stats::median(random_control$AUROC_q025[random_control$negative_control == "random_available_signature_id_score"], na.rm = TRUE),
    AUROC_q975 = stats::median(random_control$AUROC_q975[random_control$negative_control == "random_available_signature_id_score"], na.rm = TRUE),
    interpretation = "random reference-signature IDs show benchmark spread, not a null gene-signature distribution",
    stringsAsFactors = FALSE
  )
)
write_tsv_atomic(negative_summary, file.path(out_dir, "negative_control_summary.tsv"))
write_tsv_atomic(negative_summary, file.path(man_table_dir, "Supplementary_Table_negative_control_summary.tsv"))

# 10. Submission-lean main tables.
reference_locked_h2h_path <- file.path(out_dir, "reference_signature_locked_external_head_to_head.tsv")
table2_source <- if (file.exists(project_path(reference_locked_h2h_path))) {
  read_tsv(reference_locked_h2h_path)
} else {
  published_locked_h2h
}
table2 <- table2_source[table2_source$task == "primary_ATB_vs_OD", , drop = FALSE]
top_overlap_free_for_main <- c("Hoang_OD_3", "Darboe_RISK_11", "Tabone_OD_11", "Chen_HIV_4")
table2_keep <- grepl("^final_", table2$model_layer) |
  table2$model_layer == "published_signature_ensemble" |
  grepl("canonical_literature", table2$comparator_reason %||% "", fixed = TRUE) |
  grepl("development_top10", table2$selection_basis %||% "", fixed = TRUE) |
  table2$signature_or_model %in% top_overlap_free_for_main
table2 <- table2[table2_keep, , drop = FALSE]
table2 <- table2[order(table2$model_layer != "final_best_performance", -table2$AUROC), , drop = FALSE]
table2_out <- table2[, c("signature_or_model", "model_layer", "task", "evaluation_scope", "comparator_reason", "selection_basis", "overlap_status", "n_samples", "n_families", "AUROC", "AUPRC", "family_median_AUROC", "family_IQR_AUROC"), drop = FALSE]
table2_out$AUROC_ci_low <- table2_out$AUROC_ci_high <- NA_real_
table2_out$AUPRC_ci_low <- table2_out$AUPRC_ci_high <- NA_real_
for (i in seq_len(nrow(table2_out))) {
  if (grepl("^final_", table2_out$model_layer[i])) {
    role <- sub("^final_", "", table2_out$model_layer[i])
    hit <- fam_ci[fam_ci$task == table2_out$task[i] & fam_ci$selected_role == role, , drop = FALSE]
    if (nrow(hit) > 0) {
      table2_out$AUROC_ci_low[i] <- hit$AUROC_ci_low[1]
      table2_out$AUROC_ci_high[i] <- hit$AUROC_ci_high[1]
      table2_out$AUPRC_ci_low[i] <- hit$AUPRC_ci_low[1]
      table2_out$AUPRC_ci_high[i] <- hit$AUPRC_ci_high[1]
    }
  }
}
names(table2_out)[names(table2_out) == "task"] <- "analysis_id"
unlink(project_path("manuscript/tables/Table2_published_signature_locked_external_head_to_head.tsv"), force = TRUE)
mlfinal_write_display_tsv(table2_out,
  "manuscript/tables/Table2_reference_signature_locked_external_head_to_head.tsv")

table3 <- merge(ext_summary, selected[, c("task", "selected_role", "selector_classifier", "selector", "classifier", "combination_id")], by = c("task", "selected_role", "selector_classifier", "combination_id"), all.x = TRUE)
table3 <- merge(table3, dev_perf[, c("combination_id", "AUROC", "AUPRC")], by = "combination_id", all.x = TRUE, suffixes = c("_external", "_development"))
feature_col <- if ("selector.x" %in% names(table3)) "selector.x" else "selector"
classifier_col <- if ("classifier.x" %in% names(table3)) "classifier.x" else "classifier"
external_auc_col <- if ("AUROC_external" %in% names(table3)) "AUROC_external" else "AUROC.x"
external_auprc_col <- if ("AUPRC_external" %in% names(table3)) "AUPRC_external" else "AUPRC.x"
development_auc_col <- if ("AUROC_development" %in% names(table3)) "AUROC_development" else "AUROC.y"
table3_out <- data.frame(
  analysis_id = table3$task,
  model = table3$selector_classifier,
  feature_source = table3[[feature_col]],
  classifier = table3[[classifier_col]],
  development_AUROC = table3[[development_auc_col]],
  external_AUROC = table3[[external_auc_col]],
  external_AUPRC = table3[[external_auprc_col]],
  Brier = table3$Brier,
  calibration_slope = table3$calibration_slope,
  n_families_external = vapply(table3$combination_id, function(id) length(unique(ext_pred$cohort_family[ext_pred$combination_id == id])), integer(1)),
  stringsAsFactors = FALSE
)
write_tsv_atomic(table3_out, "manuscript/tables/Table3_final_selected_models.tsv")

table4_rows <- list()
for (i in seq_len(nrow(selected))) {
  combo <- selected[i, , drop = FALSE]
  for (strategy in strategy_names) {
    point <- selective_points[
      selective_points$combination_id == combo$combination_id &
        selective_points$selected_role == combo$selected_role &
        selective_points$selected_strategy == strategy, , drop = FALSE
    ]
    if (nrow(point) == 1L) {
      table4_rows[[length(table4_rows) + 1L]] <- data.frame(
        analysis_id = point$task,
        selected_role = point$selected_role,
        model = point$model,
        threshold_strategy = strategy,
        rule_status = "selected_in_development",
        lower_threshold = point$lower_threshold,
        upper_threshold = point$upper_threshold,
        sensitivity = point$external_sensitivity_answered,
        specificity = point$external_specificity_answered,
        PPV = point$external_PPV_answered %||% NA_real_,
        NPV = point$external_NPV_answered %||% NA_real_,
        coverage = point$external_coverage,
        abstention_rate = point$external_abstention_rate,
        error_rate_among_answered = point$external_error_rate_among_answered,
        false_negative_rate_among_answered = point$external_false_negative_rate_answered,
        false_positive_rate_among_answered = point$external_false_positive_rate_answered,
        max_development_sensitivity = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      grid <- selective_grid[
        selective_grid$combination_id == combo$combination_id &
          selective_grid$selected_role == combo$selected_role, , drop = FALSE
      ]
      max_dev_sens <- if (nrow(grid) > 0L) max(grid$dev_sensitivity_answered, na.rm = TRUE) else NA_real_
      status <- if (strategy == "dev_sensitivity_ge_0.90_max_coverage" &&
                    is.finite(max_dev_sens) && max_dev_sens < 0.90) {
        "infeasible_in_development"
      } else {
        "not_evaluable"
      }
      table4_rows[[length(table4_rows) + 1L]] <- data.frame(
        analysis_id = combo$task,
        selected_role = combo$selected_role,
        model = combo$selector_classifier,
        threshold_strategy = strategy,
        rule_status = status,
        lower_threshold = NA_real_,
        upper_threshold = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_,
        PPV = NA_real_,
        NPV = NA_real_,
        coverage = NA_real_,
        abstention_rate = NA_real_,
        error_rate_among_answered = NA_real_,
        false_negative_rate_among_answered = NA_real_,
        false_positive_rate_among_answered = NA_real_,
        max_development_sensitivity = max_dev_sens,
        stringsAsFactors = FALSE
      )
    }
  }
}
table4 <- do.call(rbind, table4_rows)
stopifnot(nrow(table4) == 16L,
          sum(table4$rule_status == "infeasible_in_development") == 1L,
          !any(table4$rule_status == "not_evaluable"))
write_tsv_atomic(table4, "manuscript/tables/Table4_operating_points_selective_prediction.tsv")
write_tsv_atomic(table4, file.path(out_dir, "selective_prediction_rule_feasibility_and_external_results.tsv"))
write_tsv_atomic(table4, file.path(man_table_dir, "Supplementary_Table_selective_prediction_rule_feasibility_and_external_results.tsv"))

lock_hash_after <- hash_file_or_object("config/external_validation_lockfile.json")
if (!identical(lock_hash_before, lock_hash_after)) stop("External validation lockfile changed during generalization stress-test stage.", call. = FALSE)

write_tsv_atomic(data.frame(
  item = c("external_error_by_sample", "same_lock_head_to_head", "overlap_aware_benchmark", "family_uncertainty", "mixed_family_sensitivity", "model_shopping_check", "calibration_repair", "selective_prediction", "negative_controls", "interpretation"),
  value = c(
    "joins locked external predictions to metadata and reports TP/TN/FP/FN",
    "final model, top signatures, and simple signature ensembles on the same locked external families",
    "separates all-available and overlap-free signature evaluation",
    "family-cluster bootstrap, leave-family jackknife, and equal-family summaries",
    "exposes limited mixed-family support for primary ATB-vs-OD analysis",
    "scatter-plot-ready development versus external performance review",
    "raw, Platt, isotonic, temperature, and beta calibration fitted only on development OOF predictions",
    "development-locked abstention/selective-use analysis",
    "score-label permutation, random score, random signature IDs, platform checks, and family-label oracle checks",
    "strengthens generalization-ceiling framing without tuning locked external AUROC"
  ),
  stringsAsFactors = FALSE
), "report/generalization_stress_test_enhancements.tsv")

message("Generalization stress-test enhancement analyses completed.")
