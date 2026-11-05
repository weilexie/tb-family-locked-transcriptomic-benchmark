mlfinal_ensemble_candidate_sets <- function(dev_perf, task, top_n = 12) {
  d <- dev_perf[dev_perf$task == task & is.finite(dev_perf$AUROC), , drop = FALSE]
  d <- d[order(d$AUROC, d$AUPRC, d$median_family_AUROC, decreasing = TRUE), , drop = FALSE]
  d <- head(d, min(top_n, nrow(d)))
  rows <- list()
  add <- function(name, ids, method) {
    ids <- unique(ids)
    ids <- ids[nzchar(ids)]
    if (length(ids) > 0) rows[[length(rows) + 1]] <<- data.frame(task = task, ensemble_id = paste(task, name, sep = "__"), ensemble_name = name, method = method, combination_ids = paste(ids, collapse = ";"), n_base_models = length(ids), stringsAsFactors = FALSE)
  }
  add("single_best_development", d$combination_id[1], "single")
  for (k in c(3, 5, 10)) {
    if (nrow(d) >= k) {
      add(sprintf("mean_top%s", k), d$combination_id[seq_len(k)], "mean")
      add(sprintf("stack_ridge_top%s", k), d$combination_id[seq_len(k)], "stack_ridge")
    }
  }
  do.call(rbind, rows)
}

mlfinal_wide_oof <- function(oof, combination_ids) {
  d <- oof[oof$combination_id %in% combination_ids, , drop = FALSE]
  key <- unique(d[, c("standard_sample_id", "study", "cohort_family", "true_label", "task"), drop = FALSE])
  wide <- stats::reshape(d[, c("standard_sample_id", "combination_id", "score")], idvar = "standard_sample_id", timevar = "combination_id", direction = "wide")
  names(wide) <- sub("^score\\.", "", names(wide))
  out <- merge(key, wide, by = "standard_sample_id", all.x = TRUE)
  out
}

mlfinal_prepare_score_matrix <- function(wide, ids) {
  x <- as.matrix(wide[, ids, drop = FALSE])
  storage.mode(x) <- "numeric"
  for (j in seq_len(ncol(x))) {
    med <- stats::median(x[, j], na.rm = TRUE)
    if (!is.finite(med)) med <- 0.5
    x[!is.finite(x[, j]) | is.na(x[, j]), j] <- med
  }
  x
}

mlfinal_fit_stack_ridge <- function(x_train, y_train, x_test) {
  y <- as_binary(y_train)
  if (!requireNamespace("glmnet", quietly = TRUE) || length(unique(y)) < 2) {
    return(rowMeans(x_test))
  }
  fit <- glmnet::glmnet(x_train, y, family = "binomial", alpha = 0, weights = class_weights(y))
  pred_train <- stats::predict(fit, newx = x_train, type = "response")
  aucs <- apply(pred_train, 2, function(p) compute_roc_auc(y, as.numeric(p)))
  idx <- which.max(ifelse(is.na(aucs), -Inf, aucs))
  as.numeric(stats::predict(fit, newx = x_test, s = fit$lambda[idx], type = "response"))
}

mlfinal_meta_cv_predict <- function(wide, ids, method, seed = 20260501) {
  meta <- unique(wide[, c("standard_sample_id", "cohort_family"), drop = FALSE])
  folds <- automl_final_group_folds(meta, "cohort_family", k = 5, seed = seed)
  out <- list()
  for (fold in folds) {
    train_idx <- wide$standard_sample_id %in% fold$train_ids
    test_idx <- wide$standard_sample_id %in% fold$test_ids
    x_train <- mlfinal_prepare_score_matrix(wide[train_idx, , drop = FALSE], ids)
    x_test <- mlfinal_prepare_score_matrix(wide[test_idx, , drop = FALSE], ids)
    score <- if (method == "single") {
      as.numeric(x_test[, 1])
    } else if (method == "mean") {
      rowMeans(x_test)
    } else if (method == "stack_ridge") {
      mlfinal_fit_stack_ridge(x_train, wide$true_label[train_idx], x_test)
    } else {
      stop(sprintf("unknown_ensemble_method_%s", method), call. = FALSE)
    }
    out[[length(out) + 1]] <- data.frame(
      task = wide$task[test_idx],
      standard_sample_id = wide$standard_sample_id[test_idx],
      study = wide$study[test_idx],
      cohort_family = wide$cohort_family[test_idx],
      true_label = wide$true_label[test_idx],
      score = score,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

mlfinal_ensemble_oof_summary <- function(oof, candidate_row) {
  ids <- strsplit(candidate_row$combination_ids[1], ";", fixed = TRUE)[[1]]
  wide <- mlfinal_wide_oof(oof, ids)
  pred <- mlfinal_meta_cv_predict(wide, ids, candidate_row$method[1])
  m <- compute_binary_metrics(pred$true_label, pred$score)
  fam <- summarize_per_study_metrics(pred, "cohort_family")
  data.frame(
    task = candidate_row$task[1],
    ensemble_id = candidate_row$ensemble_id[1],
    ensemble_name = candidate_row$ensemble_name[1],
    method = candidate_row$method[1],
    n_base_models = candidate_row$n_base_models[1],
    development_AUROC = m$AUROC,
    development_AUPRC = m$AUPRC,
    development_Brier = m$Brier,
    development_median_family_AUROC = stats::median(fam$AUROC, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

mlfinal_external_base_predictions <- function(task, combination_ids, grid, lock_hash_before, seed = 20260501) {
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  xdf <- read_tsv(file.path("data_curated/task_matrices", task, "X_raw_gene_by_sample.tsv.gz"))
  genes <- xdf$gene
  X_raw <- as.matrix(xdf[, setdiff(names(xdf), "gene"), drop = FALSE])
  storage.mode(X_raw) <- "numeric"
  rownames(X_raw) <- genes
  lock <- read_external_lockfile("config/external_validation_lockfile.json")
  dev_meta <- meta[meta$cohort_family %in% unlist(lock$development_families), , drop = FALSE]
  ext_meta <- meta[meta$cohort_family %in% unlist(lock$external_validation_families), , drop = FALSE]
  all_sig <- automl_final_load_signature_matrix(task, meta$standard_sample_id, max_sigs = 30)
  td <- list(meta = meta, X_raw = X_raw, dev_meta = dev_meta, ext_meta = ext_meta, all_sig = all_sig)
  rows <- list()
  for (i in seq_along(combination_ids)) {
    combo <- grid[grid$combination_id == combination_ids[i], , drop = FALSE][1, , drop = FALSE]
    combo$selected_role <- "ensemble_base"
    pred <- automl_final_external_evaluate_combination(td, combo, all_sig, lock_hash_before, seed = seed + i)
    rows[[length(rows) + 1]] <- pred
  }
  do.call(rbind, rows)
}

mlfinal_ensemble_external_predict <- function(oof, external_base, candidate_row) {
  ids <- strsplit(candidate_row$combination_ids[1], ";", fixed = TRUE)[[1]]
  dev_wide <- mlfinal_wide_oof(oof, ids)
  ext_key <- unique(external_base[, c("standard_sample_id", "study", "cohort_family", "true_label", "task"), drop = FALSE])
  ext_wide <- stats::reshape(external_base[, c("standard_sample_id", "combination_id", "score")], idvar = "standard_sample_id", timevar = "combination_id", direction = "wide")
  names(ext_wide) <- sub("^score\\.", "", names(ext_wide))
  ext_wide <- merge(ext_key, ext_wide, by = "standard_sample_id", all.x = TRUE)
  x_dev <- mlfinal_prepare_score_matrix(dev_wide, ids)
  x_ext <- mlfinal_prepare_score_matrix(ext_wide, ids)
  score <- if (candidate_row$method[1] == "single") {
    as.numeric(x_ext[, 1])
  } else if (candidate_row$method[1] == "mean") {
    rowMeans(x_ext)
  } else if (candidate_row$method[1] == "stack_ridge") {
    mlfinal_fit_stack_ridge(x_dev, dev_wide$true_label, x_ext)
  } else {
    stop(sprintf("unknown_ensemble_method_%s", candidate_row$method[1]), call. = FALSE)
  }
  data.frame(
    task = ext_wide$task,
    ensemble_id = candidate_row$ensemble_id[1],
    ensemble_name = candidate_row$ensemble_name[1],
    method = candidate_row$method[1],
    n_base_models = candidate_row$n_base_models[1],
    standard_sample_id = ext_wide$standard_sample_id,
    study = ext_wide$study,
    cohort_family = ext_wide$cohort_family,
    true_label = ext_wide$true_label,
    score = score,
    stringsAsFactors = FALSE
  )
}

mlfinal_ensemble_external_summary <- function(pred) {
  groups <- split(pred, pred$ensemble_id)
  out <- lapply(groups, function(d) {
    m <- compute_binary_metrics(d$true_label, d$score)
    cal <- calibration_summary(d$true_label, d$score)
    data.frame(
      task = d$task[1],
      ensemble_id = d$ensemble_id[1],
      ensemble_name = d$ensemble_name[1],
      method = d$method[1],
      n_base_models = d$n_base_models[1],
      n = m$n,
      n_pos = m$n_pos,
      AUROC = m$AUROC,
      AUPRC = m$AUPRC,
      sensitivity = m$sensitivity,
      specificity = m$specificity,
      balanced_accuracy = m$balanced_accuracy,
      Brier = m$Brier,
      calibration_intercept = cal$calibration_intercept,
      calibration_slope = cal$calibration_slope,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
