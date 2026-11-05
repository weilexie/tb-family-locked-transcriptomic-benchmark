# ML final task grid, package status, and grouped folds
# Split from modeling_algorithm_combinations_final.R for maintainability.

automl_final_tasks <- function() {
  c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB")
}

automl_final_normalizations <- function() {
  c("deployment_locked_training_scaler", "samplewise_rank_percentile")
}

automl_final_selectors <- function() {
  c(
    "stable_top20", "stable_top50", "stable_top100", "stable_top200", "stable_top500",
    "lasso_selector", "ridge_rank_selector", "elastic_net_selector",
    "stepglm_forward", "stepglm_backward", "stepglm_both",
    "gbm_importance_selector", "xgboost_importance_selector", "rf_importance_selector",
    "published_signature_scores", "stable_top100_plus_signature_scores", "stable_top200_plus_signature_scores"
  )
}

automl_final_classifiers <- function() {
  c(
    "ridge_logistic", "lasso_logistic", "elastic_net_025", "elastic_net_050", "elastic_net_075",
    "stepglm_logistic", "lda", "naive_bayes", "pls_glm", "svm_linear", "svm_rbf",
    "gbm", "xgboost", "random_forest"
  )
}

automl_final_package_status <- function() {
  pkgs <- c("glmnet", "gbm", "xgboost", "ranger", "randomForest", "e1071", "MASS", "pls", "pROC")
  data.frame(package = pkgs, available = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE), stringsAsFactors = FALSE)
}

automl_final_build_grid <- function(tasks = automl_final_tasks(), normalizations = automl_final_normalizations(), selectors = automl_final_selectors(), classifiers = automl_final_classifiers()) {
  grid <- expand.grid(task = tasks, normalization = normalizations, selector = selectors, classifier = classifiers, stringsAsFactors = FALSE)
  grid$combination_id <- sprintf("%s__%s__%s+%s", grid$task, grid$normalization, grid$selector, grid$classifier)
  grid$selector_classifier <- sprintf("%s+%s", grid$selector, grid$classifier)
  grid$status <- "pending"
  grid$skip_reason <- ""
  grid
}

automl_final_group_folds <- function(meta, group_col = "cohort_family", k = 5, seed = 20260501) {
  groups <- sort(unique(meta[[group_col]]))
  sizes <- as.integer(table(meta[[group_col]])[groups])
  ord <- order(sizes, decreasing = TRUE)
  groups <- groups[ord]
  fold_load <- rep(0, min(k, length(groups)))
  fold_groups <- vector("list", length(fold_load))
  for (g in groups) {
    idx <- which.min(fold_load)
    fold_groups[[idx]] <- c(fold_groups[[idx]], g)
    fold_load[idx] <- fold_load[idx] + sum(meta[[group_col]] == g)
  }
  lapply(seq_along(fold_groups), function(i) {
    test_groups <- fold_groups[[i]]
    list(
      fold_id = sprintf("groupfold_%02d", i),
      heldout_groups = test_groups,
      train_ids = meta$standard_sample_id[!meta[[group_col]] %in% test_groups],
      test_ids = meta$standard_sample_id[meta[[group_col]] %in% test_groups]
    )
  })
}

