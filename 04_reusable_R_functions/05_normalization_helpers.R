fit_normalizer <- function(X_train, mode, gene_universe = NULL) {
  gene_universe <- gene_universe %||% rownames(X_train)
  X_train <- X_train[gene_universe, , drop = FALSE]
  if (mode == "deployment_locked_training_scaler") {
    center <- apply(X_train, 1, median, na.rm = TRUE)
    scale <- apply(X_train, 1, mad, na.rm = TRUE)
    scale[is.na(scale) | scale == 0] <- 1
    list(
      mode = mode,
      gene_universe = gene_universe,
      center = center,
      scale = scale,
      fit_sample_ids = colnames(X_train)
    )
  } else if (mode == "samplewise_rank_percentile") {
    list(
      mode = mode,
      gene_universe = gene_universe,
      fit_sample_ids = colnames(X_train)
    )
  } else if (mode == "cohortwise_zscore_benchmark_only") {
    center <- rowMeans(X_train, na.rm = TRUE)
    scale <- apply(X_train, 1, sd, na.rm = TRUE)
    scale[is.na(scale) | scale == 0] <- 1
    list(
      mode = mode,
      gene_universe = gene_universe,
      center = center,
      scale = scale,
      fit_sample_ids = colnames(X_train),
      benchmark_only = TRUE
    )
  } else {
    stop(sprintf("Unknown normalization mode: %s", mode), call. = FALSE)
  }
}

transform_normalizer <- function(fit, X) {
  genes <- fit$gene_universe
  missing <- setdiff(genes, rownames(X))
  if (length(missing) > 0) stop(sprintf("Missing genes for normalization: %s", paste(missing, collapse = ", ")), call. = FALSE)
  X <- X[genes, , drop = FALSE]
  if (fit$mode %in% c("deployment_locked_training_scaler", "cohortwise_zscore_benchmark_only")) {
    out <- sweep(sweep(X, 1, fit$center, "-"), 1, fit$scale, "/")
    out[is.na(out)] <- 0
    out
  } else if (fit$mode == "samplewise_rank_percentile") {
    out <- normalize_samplewise_rank(X, genes)
    out[is.na(out)] <- 0
    out
  } else {
    stop(sprintf("Unknown normalizer mode: %s", fit$mode), call. = FALSE)
  }
}

normalize_samplewise_rank <- function(X, gene_universe = rownames(X)) {
  X <- X[gene_universe, , drop = FALSE]
  apply(X, 2, function(v) rank(v, ties.method = "average", na.last = "keep") / sum(!is.na(v)))
}

assert_normalizer_no_test_leakage <- function(fit, train_ids, test_ids) {
  leaked <- intersect(fit$fit_sample_ids %||% character(), test_ids)
  if (length(leaked) > 0) {
    stop(sprintf("Normalizer fit contains held-out/test samples: %s", paste(leaked, collapse = ", ")), call. = FALSE)
  }
  if (!all((fit$fit_sample_ids %||% character()) %in% train_ids)) {
    stop("Normalizer fit contains samples not present in declared training IDs.", call. = FALSE)
  }
  TRUE
}

detect_and_log2_transform <- function(x, q99_gt = 100, max_gt = 1000, offset = 1) {
  q99 <- as.numeric(stats::quantile(x, 0.99, na.rm = TRUE))
  mx <- max(x, na.rm = TRUE)
  has_negative <- any(x < 0, na.rm = TRUE)
  if (!has_negative && (q99 > q99_gt || mx > max_gt)) {
    list(expression = log2(x + offset), decision = "log2_transformed", q99 = q99, max = mx, has_negative = has_negative)
  } else if (has_negative || mx < max_gt) {
    list(expression = x, decision = "already_log", q99 = q99, max = mx, has_negative = has_negative)
  } else {
    list(expression = x, decision = "uncertain_no_transform", q99 = q99, max = mx, has_negative = has_negative)
  }
}
