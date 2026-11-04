# ML final feature and signature selectors
# Split from modeling_algorithm_combinations_final.R for maintainability.

automl_final_univariate_rank <- function(X_train_raw, y_train, max_rank = 700) {
  y <- as_binary(y_train)
  present <- rowMeans(!is.na(X_train_raw))
  keep <- present >= 0.60
  X <- X_train_raw[keep, , drop = FALSE]
  if (nrow(X) == 0 || length(unique(y)) < 2) return(character())
  pos_mean <- rowMeans(X[, y == 1, drop = FALSE], na.rm = TRUE)
  neg_mean <- rowMeans(X[, y == 0, drop = FALSE], na.rm = TRUE)
  pos_var <- apply(X[, y == 1, drop = FALSE], 1, stats::var, na.rm = TRUE)
  neg_var <- apply(X[, y == 0, drop = FALSE], 1, stats::var, na.rm = TRUE)
  se <- sqrt(pos_var / max(sum(y == 1), 1) + neg_var / max(sum(y == 0), 1))
  se[is.na(se) | se == 0] <- stats::median(se[is.finite(se) & se > 0], na.rm = TRUE)
  se[is.na(se) | se == 0] <- 1
  score <- abs(pos_mean - neg_mean) / se
  score[!is.finite(score)] <- 0
  ranked <- names(sort(score, decreasing = TRUE))
  head(ranked, min(max_rank, length(ranked)))
}

automl_final_safe_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- NA_real_
  x
}

automl_final_impute_scale_samples <- function(train, test) {
  train <- automl_final_safe_matrix(train)
  test <- automl_final_safe_matrix(test)
  med <- apply(train, 2, stats::median, na.rm = TRUE)
  med[is.na(med)] <- 0
  for (j in seq_len(ncol(train))) {
    train[is.na(train[, j]), j] <- med[j]
    test[is.na(test[, j]), j] <- med[j]
  }
  center <- colMeans(train, na.rm = TRUE)
  scale <- apply(train, 2, stats::sd, na.rm = TRUE)
  scale[is.na(scale) | scale == 0] <- 1
  train <- sweep(sweep(train, 2, center, "-"), 2, scale, "/")
  test <- sweep(sweep(test, 2, center, "-"), 2, scale, "/")
  list(train = train, test = test, center = center, scale = scale)
}

automl_final_reduce_zero_variance <- function(train, test, max_features = Inf) {
  if (ncol(train) == 0) return(list(train = train, test = test, features = character()))
  vars <- apply(train, 2, stats::var, na.rm = TRUE)
  keep <- is.finite(vars) & vars > 1e-10
  if (!any(keep)) return(list(train = train[, FALSE, drop = FALSE], test = test[, FALSE, drop = FALSE], features = character()))
  train <- train[, keep, drop = FALSE]
  test <- test[, colnames(train), drop = FALSE]
  if (is.finite(max_features) && ncol(train) > max_features) {
    train <- train[, seq_len(max_features), drop = FALSE]
    test <- test[, colnames(train), drop = FALSE]
  }
  list(train = train, test = test, features = colnames(train))
}

automl_final_load_signature_matrix <- function(task, sample_ids, max_sigs = 30) {
  path <- project_path("results/benchmark/signature_scores_by_sample.tsv.gz")
  if (!file.exists(path)) return(matrix(nrow = length(sample_ids), ncol = 0, dimnames = list(sample_ids, character())))
  sig <- read_tsv("results/benchmark/signature_scores_by_sample.tsv.gz")
  sig <- sig[sig$task == task & sig$standard_sample_id %in% sample_ids, , drop = FALSE]
  if (nrow(sig) == 0) return(matrix(nrow = length(sample_ids), ncol = 0, dimnames = list(sample_ids, character())))
  preferred <- sig[sig$scoring_method == "signed_mean_rank", , drop = FALSE]
  if (nrow(preferred) > 0) sig <- preferred
  summary_rel <- "results/benchmark/benchmark_rank_table_development.tsv"
  if (!file.exists(project_path(summary_rel))) {
    summary_rel <- "results/benchmark/benchmark_rank_table.tsv"
  }
  if (file.exists(project_path(summary_rel))) {
    rank <- read_tsv(summary_rel)
    rank <- rank[rank$task == task, , drop = FALSE]
    rank <- rank[order(rank$AUROC, rank$AUPRC, decreasing = TRUE), , drop = FALSE]
    top <- head(unique(rank$signature_id), max_sigs)
  } else {
    top <- head(sort(unique(sig$signature_id)), max_sigs)
  }
  sig <- sig[sig$signature_id %in% top, , drop = FALSE]
  if (nrow(sig) == 0) return(matrix(nrow = length(sample_ids), ncol = 0, dimnames = list(sample_ids, character())))
  wide <- stats::reshape(sig[, c("standard_sample_id", "signature_id", "score")], idvar = "standard_sample_id", timevar = "signature_id", direction = "wide")
  rn <- wide$standard_sample_id
  wide$standard_sample_id <- NULL
  mat <- as.matrix(wide)
  storage.mode(mat) <- "numeric"
  colnames(mat) <- sub("^score\\.", "sig__", colnames(mat))
  out <- matrix(NA_real_, nrow = length(sample_ids), ncol = ncol(mat), dimnames = list(sample_ids, colnames(mat)))
  out[match(rn, sample_ids), ] <- mat
  out
}

automl_final_glmnet_select <- function(X_train_raw, y_train, alpha = 1, top_n = 80) {
  if (!requireNamespace("glmnet", quietly = TRUE)) return(character())
  ranked <- automl_final_univariate_rank(X_train_raw, y_train, max_rank = 300)
  if (length(ranked) < 2) return(ranked)
  genes <- head(ranked, min(300, length(ranked)))
  x <- t(X_train_raw[genes, , drop = FALSE])
  x[is.na(x)] <- 0
  y <- as_binary(y_train)
  if (length(unique(y)) < 2) return(head(ranked, top_n))
  fit <- try(glmnet::glmnet(x, y, family = "binomial", alpha = alpha, weights = class_weights(y)), silent = TRUE)
  if (inherits(fit, "try-error")) return(head(ranked, top_n))
  pred <- try(stats::predict(fit, newx = x, type = "response"), silent = TRUE)
  if (inherits(pred, "try-error")) return(head(ranked, top_n))
  aucs <- apply(pred, 2, function(p) compute_roc_auc(y, as.numeric(p)))
  idx <- which.max(ifelse(is.na(aucs), -Inf, aucs))
  co <- as.matrix(stats::coef(fit, s = fit$lambda[idx]))[-1, 1]
  ranked_coef <- names(sort(abs(co), decreasing = TRUE))
  ranked_coef <- ranked_coef[is.finite(abs(co[ranked_coef])) & abs(co[ranked_coef]) > 1e-8]
  out <- unique(c(ranked_coef, ranked))
  head(out, min(top_n, length(out)))
}

automl_final_step_select <- function(X_train_raw, y_train, direction = "both", top_n = 40) {
  ranked <- automl_final_univariate_rank(X_train_raw, y_train, max_rank = 80)
  if (length(ranked) < 3 || !requireNamespace("MASS", quietly = TRUE)) return(head(ranked, top_n))
  genes <- head(ranked, min(25, length(ranked)))
  x <- as.data.frame(t(X_train_raw[genes, , drop = FALSE]), check.names = FALSE)
  names(x) <- make.names(names(x), unique = TRUE)
  x[is.na(x)] <- 0
  x$y <- as_binary(y_train)
  if (length(unique(x$y)) < 2) return(head(ranked, top_n))
  fit0 <- try(suppressWarnings(stats::glm(y ~ 1, data = x, family = stats::binomial())), silent = TRUE)
  fit1 <- try(suppressWarnings(stats::glm(y ~ ., data = x, family = stats::binomial())), silent = TRUE)
  if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) return(head(ranked, top_n))
  st <- try(suppressWarnings(stats::step(fit0, scope = list(lower = fit0, upper = fit1), direction = direction, trace = 0, steps = 12)), silent = TRUE)
  if (inherits(st, "try-error")) return(head(ranked, top_n))
  vars <- names(stats::coef(st))[-1]
  original <- genes[match(vars, make.names(genes, unique = TRUE))]
  original <- original[!is.na(original)]
  out <- unique(c(original, ranked))
  head(out, min(top_n, length(out)))
}

automl_final_importance_select <- function(X_train_raw, y_train, method = "gbm", top_n = 80) {
  ranked <- automl_final_univariate_rank(X_train_raw, y_train, max_rank = 250)
  if (length(ranked) < 3) return(ranked)
  genes <- head(ranked, min(180, length(ranked)))
  x <- as.data.frame(t(X_train_raw[genes, , drop = FALSE]), check.names = FALSE)
  names(x) <- make.names(names(x), unique = TRUE)
  x[is.na(x)] <- 0
  y <- as_binary(y_train)
  if (length(unique(y)) < 2) return(head(ranked, top_n))
  imp <- NULL
  if (method == "gbm" && requireNamespace("gbm", quietly = TRUE)) {
    df <- x; df$y <- y
    fit <- try(suppressWarnings(gbm::gbm(y ~ ., data = df, distribution = "bernoulli", n.trees = 40, interaction.depth = 1, shrinkage = 0.05, bag.fraction = 0.8, train.fraction = 1, verbose = FALSE)), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      rel <- try(summary(fit, plotit = FALSE), silent = TRUE)
      if (!inherits(rel, "try-error")) imp <- setNames(rel$rel.inf, rel$var)
    }
  } else if (method == "xgboost" && requireNamespace("xgboost", quietly = TRUE)) {
    fit <- try(xgboost::xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2, eta = 0.08, subsample = 0.8, colsample_bytree = 0.8, nthread = 1), data = xgboost::xgb.DMatrix(as.matrix(x), label = y), nrounds = 35, verbose = 0), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      rel <- try(xgboost::xgb.importance(model = fit), silent = TRUE)
      if (!inherits(rel, "try-error") && nrow(rel) > 0) imp <- setNames(rel$Gain, rel$Feature)
    }
  } else if (method == "rf" && requireNamespace("ranger", quietly = TRUE)) {
    df <- x; df$y <- factor(ifelse(y == 1, "ATB", "nonATB"), levels = c("nonATB", "ATB"))
    fit <- try(ranger::ranger(y ~ ., data = df, probability = TRUE, num.trees = 80, mtry = max(1, floor(sqrt(ncol(x)))), importance = "impurity", seed = 20260501, num.threads = 1), silent = TRUE)
    if (!inherits(fit, "try-error")) imp <- fit$variable.importance
  }
  if (is.null(imp) || length(imp) == 0) return(head(ranked, top_n))
  imp <- sort(imp, decreasing = TRUE)
  selected <- genes[match(names(imp), make.names(genes, unique = TRUE))]
  selected <- selected[!is.na(selected)]
  head(unique(c(selected, ranked)), min(top_n, length(ranked)))
}

automl_final_select_features <- function(selector, X_train_raw, y_train, sig_feature_names = character()) {
  ranked <- automl_final_univariate_rank(X_train_raw, y_train, max_rank = 700)
  empty <- list(genes = character(), signatures = character(), status = "completed", reason = "")
  if (selector == "stable_top20") return(modifyList(empty, list(genes = head(ranked, 20))))
  if (selector == "stable_top50") return(modifyList(empty, list(genes = head(ranked, 50))))
  if (selector == "stable_top100") return(modifyList(empty, list(genes = head(ranked, 100))))
  if (selector == "stable_top200") return(modifyList(empty, list(genes = head(ranked, 200))))
  if (selector == "stable_top500") return(modifyList(empty, list(genes = head(ranked, 500))))
  if (selector == "lasso_selector") return(modifyList(empty, list(genes = automl_final_glmnet_select(X_train_raw, y_train, alpha = 1, top_n = 80))))
  if (selector == "ridge_rank_selector") return(modifyList(empty, list(genes = automl_final_glmnet_select(X_train_raw, y_train, alpha = 0, top_n = 100))))
  if (selector == "elastic_net_selector") return(modifyList(empty, list(genes = automl_final_glmnet_select(X_train_raw, y_train, alpha = 0.5, top_n = 80))))
  if (selector == "stepglm_forward") return(modifyList(empty, list(genes = automl_final_step_select(X_train_raw, y_train, direction = "forward", top_n = 40))))
  if (selector == "stepglm_backward") return(modifyList(empty, list(genes = automl_final_step_select(X_train_raw, y_train, direction = "backward", top_n = 40))))
  if (selector == "stepglm_both") return(modifyList(empty, list(genes = automl_final_step_select(X_train_raw, y_train, direction = "both", top_n = 40))))
  if (selector == "gbm_importance_selector") return(modifyList(empty, list(genes = automl_final_importance_select(X_train_raw, y_train, method = "gbm", top_n = 80))))
  if (selector == "xgboost_importance_selector") return(modifyList(empty, list(genes = automl_final_importance_select(X_train_raw, y_train, method = "xgboost", top_n = 80))))
  if (selector == "rf_importance_selector") {
    if (!requireNamespace("ranger", quietly = TRUE) && !requireNamespace("randomForest", quietly = TRUE)) return(modifyList(empty, list(status = "skipped", reason = "missing_random_forest_package")))
    return(modifyList(empty, list(genes = automl_final_importance_select(X_train_raw, y_train, method = "rf", top_n = 80))))
  }
  if (selector == "published_signature_scores") return(modifyList(empty, list(signatures = sig_feature_names)))
  if (selector == "stable_top100_plus_signature_scores") return(modifyList(empty, list(genes = head(ranked, 100), signatures = sig_feature_names)))
  if (selector == "stable_top200_plus_signature_scores") return(modifyList(empty, list(genes = head(ranked, 200), signatures = sig_feature_names)))
  modifyList(empty, list(status = "skipped", reason = "unknown_selector"))
}

