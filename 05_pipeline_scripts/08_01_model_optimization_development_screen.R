# ML final development-screen stage
# Sourced by scripts/08_ml_final_optimization.R.

config <- read_config("config/config.yml")
out_dir <- "results/ml_final"
dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)

tasks <- automl_final_tasks()
normalizations <- automl_final_normalizations()
selectors <- automl_final_selectors()
classifiers <- automl_final_classifiers()
grid <- automl_final_build_grid(tasks, normalizations, selectors, classifiers)
lock <- read_external_lockfile("config/external_validation_lockfile.json")
assert_no_overlap_between_dev_and_external(lock)
lock_hash_before <- hash_file_or_object("config/external_validation_lockfile.json")
package_status <- automl_final_package_status()
write_tsv_atomic(package_status, file.path(out_dir, "package_status.tsv"))

message(sprintf("ML final optimization grid: %s combinations", nrow(grid)))

resume_existing <- all(file.exists(project_path(file.path(out_dir, c(
  "algorithm_combination_grid.tsv",
  "development_combination_performance.tsv",
  "development_oof_predictions.tsv.gz",
  "selected_combinations.tsv"
))))) && !identical(Sys.getenv("TBML_FORCE_OPTIMIZE_FINAL", unset = "0"), "1")

if (resume_existing) {
  message("Existing final development-screen outputs found; resuming from selected models to external validation/reporting.")
  grid_status <- read_tsv(file.path(out_dir, "algorithm_combination_grid.tsv"))
  dev_perf <- read_tsv(file.path(out_dir, "development_combination_performance.tsv"))
  selected <- read_tsv(file.path(out_dir, "selected_combinations.tsv"))
  task_data_store <- list()
  for (task in tasks) {
    meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
    xdf <- read_tsv(file.path("data_curated/task_matrices", task, "X_raw_gene_by_sample.tsv.gz"))
    genes <- xdf$gene
    X_raw <- as.matrix(xdf[, setdiff(names(xdf), "gene"), drop = FALSE])
    storage.mode(X_raw) <- "numeric"
    rownames(X_raw) <- genes
    dev_meta <- meta[meta$cohort_family %in% unlist(lock$development_families), , drop = FALSE]
    ext_meta <- meta[meta$cohort_family %in% unlist(lock$external_validation_families), , drop = FALSE]
    all_sig <- automl_final_load_signature_matrix(task, meta$standard_sample_id, max_sigs = 30)
    task_data_store[[task]] <- list(meta = meta, X_raw = X_raw, dev_meta = dev_meta, ext_meta = ext_meta, all_sig = all_sig)
  }
} else {
all_prediction_rows <- list()
grid_rows <- list()
task_data_store <- list()

for (task in tasks) {
  message(sprintf("Starting algorithm-combination screen for %s", task))
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  xdf <- read_tsv(file.path("data_curated/task_matrices", task, "X_raw_gene_by_sample.tsv.gz"))
  genes <- xdf$gene
  X_raw <- as.matrix(xdf[, setdiff(names(xdf), "gene"), drop = FALSE])
  storage.mode(X_raw) <- "numeric"
  rownames(X_raw) <- genes
  dev_meta <- meta[meta$cohort_family %in% unlist(lock$development_families), , drop = FALSE]
  ext_meta <- meta[meta$cohort_family %in% unlist(lock$external_validation_families), , drop = FALSE]
  all_sig <- automl_final_load_signature_matrix(task, meta$standard_sample_id, max_sigs = 30)
  task_data_store[[task]] <- list(meta = meta, X_raw = X_raw, dev_meta = dev_meta, ext_meta = ext_meta, all_sig = all_sig)
  folds <- automl_final_group_folds(dev_meta, "cohort_family", k = 5, seed = config$project$seed %||% 20260501)
  task_grid <- grid[grid$task == task, , drop = FALSE]
  selector_cache <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(task_grid))) {
    combo <- task_grid[i, , drop = FALSE]
    if (i %% 25 == 0) message(sprintf("%s: %s/%s combinations processed", task, i, nrow(task_grid)))
    combo_status <- "completed"
    combo_reason <- ""
    combo_preds <- list()
    n_model_features <- integer()
    n_gene_features <- integer()
    n_signature_features <- integer()

    if (!automl_final_classifier_available(combo$classifier)) {
      combo_status <- "skipped"
      combo_reason <- sprintf("missing_package_for_%s", combo$classifier)
    }

    if (identical(combo_status, "completed")) {
      for (fold in folds) {
        train_meta <- dev_meta[dev_meta$standard_sample_id %in% fold$train_ids, , drop = FALSE]
        test_meta <- dev_meta[dev_meta$standard_sample_id %in% fold$test_ids, , drop = FALSE]
        cache_key <- paste(fold$fold_id, combo$selector, sep = "__")
        if (!exists(cache_key, envir = selector_cache, inherits = FALSE)) {
          selected <- try(automl_final_select_features(
            combo$selector,
            X_raw[, train_meta$standard_sample_id, drop = FALSE],
            train_meta$task_label,
            colnames(all_sig)
          ), silent = TRUE)
          if (inherits(selected, "try-error")) {
            selected <- list(genes = character(), signatures = character(), status = "skipped", reason = as.character(selected))
          }
          assign(cache_key, selected, envir = selector_cache)
        } else {
          selected <- get(cache_key, envir = selector_cache, inherits = FALSE)
        }
        if (!identical(selected$status, "completed")) {
          combo_status <- "skipped"
          combo_reason <- selected$reason %||% "selector_failed"
          break
        }
        if (length(selected$genes) + length(selected$signatures) == 0) {
          combo_status <- "skipped"
          combo_reason <- "selector_returned_no_features"
          break
        }
        sig_train <- all_sig[train_meta$standard_sample_id, , drop = FALSE]
        sig_test <- all_sig[test_meta$standard_sample_id, , drop = FALSE]
        one_fold <- try({
          prepared <- automl_final_prepare_xy(
            X_raw[, train_meta$standard_sample_id, drop = FALSE],
            X_raw[, test_meta$standard_sample_id, drop = FALSE],
            train_meta$standard_sample_id,
            test_meta$standard_sample_id,
            selected,
            combo$normalization,
            sig_train,
            sig_test,
            combo$classifier
          )
          score <- suppressWarnings(automl_final_fit_predict(
            combo$classifier,
            prepared$train,
            train_meta$task_label,
            prepared$test,
            seed = (config$project$seed %||% 20260501) + i
          ))
          data.frame(
            task = task,
            standard_sample_id = test_meta$standard_sample_id,
            study = test_meta$study,
            cohort_family = test_meta$cohort_family,
            true_label = test_meta$task_label,
            score = score,
            normalization = combo$normalization,
            selector = combo$selector,
            classifier = combo$classifier,
            selector_classifier = combo$selector_classifier,
            combination_id = combo$combination_id,
            fold_id = fold$fold_id,
            stringsAsFactors = FALSE
          )
        }, silent = TRUE)
        if (inherits(one_fold, "try-error")) {
          combo_status <- "skipped"
          combo_reason <- paste0("fit_error: ", substr(gsub("[\r\n\t]+", " ", as.character(one_fold)), 1, 180))
          break
        }
        combo_preds[[length(combo_preds) + 1]] <- one_fold
        n_model_features <- c(n_model_features, length(unique(c(colnames(prepared$train)))))
        n_gene_features <- c(n_gene_features, length(selected$genes))
        n_signature_features <- c(n_signature_features, length(selected$signatures))
      }
    }

    if (identical(combo_status, "completed")) {
      pred <- do.call(rbind, combo_preds)
      if (nrow(pred) != nrow(dev_meta) || length(unique(pred$standard_sample_id)) != nrow(dev_meta)) {
        combo_status <- "skipped"
        combo_reason <- "incomplete_oof_predictions"
      } else {
        all_prediction_rows[[length(all_prediction_rows) + 1]] <- pred
      }
    }

    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      task = task,
      normalization = combo$normalization,
      selector = combo$selector,
      classifier = combo$classifier,
      selector_classifier = combo$selector_classifier,
      combination_id = combo$combination_id,
      status = combo_status,
      skip_reason = combo_reason,
      n_folds = length(folds),
      n_oof_predictions = if (identical(combo_status, "completed")) nrow(dev_meta) else 0,
      median_gene_features = if (length(n_gene_features)) stats::median(n_gene_features) else NA_real_,
      median_signature_features = if (length(n_signature_features)) stats::median(n_signature_features) else NA_real_,
      median_model_features = if (length(n_model_features)) stats::median(n_model_features) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

grid_status <- do.call(rbind, grid_rows)
write_tsv_atomic(grid_status, file.path(out_dir, "algorithm_combination_grid.tsv"))
if (length(all_prediction_rows) == 0) stop("No final algorithm-combination predictions were generated.", call. = FALSE)
all_predictions <- do.call(rbind, all_prediction_rows)
write_tsv_atomic(all_predictions, file.path(out_dir, "development_oof_predictions.tsv.gz"))

dev_perf <- automl_final_pooled_summary(all_predictions)
dev_perf <- merge(dev_perf, grid_status[, c("combination_id", "status", "skip_reason", "median_gene_features", "median_signature_features", "median_model_features")], by = "combination_id", all.x = TRUE)
dev_perf <- dev_perf[order(dev_perf$task, dev_perf$AUROC, dev_perf$AUPRC, decreasing = TRUE), , drop = FALSE]
write_tsv_atomic(dev_perf, file.path(out_dir, "development_combination_performance.tsv"))

selected_rows <- list()
for (task in tasks) {
  perf_task <- dev_perf[dev_perf$task == task, , drop = FALSE]
  sel <- automl_final_select_best(perf_task)
  if (nrow(sel$best) == 0) stop(sprintf("No selected model for %s", task), call. = FALSE)
  best <- sel$best
  best$selected_role <- "best_performance"
  selected_rows[[length(selected_rows) + 1]] <- best
  sparse <- sel$sparse
  if (nrow(sparse) > 0) {
    sparse$selected_role <- "best_sparse"
    selected_rows[[length(selected_rows) + 1]] <- sparse
  }
  short <- automl_final_task_short(task)
  payload <- list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    task = task,
    selection_provenance = "development_family_grouped_oof_only",
    lockfile_hash = lock_hash_before,
    best_performance = as.list(best[1, setdiff(names(best), "selected_role"), drop = FALSE]),
    best_sparse = if (nrow(sparse) > 0) as.list(sparse[1, setdiff(names(sparse), "selected_role"), drop = FALSE]) else NULL
  )
  write_json_atomic(payload, file.path(out_dir, sprintf("selected_combination_%s.json", short)))
}
selected <- do.call(rbind, selected_rows)
write_tsv_atomic(selected, file.path(out_dir, "selected_combinations.tsv"))
}

