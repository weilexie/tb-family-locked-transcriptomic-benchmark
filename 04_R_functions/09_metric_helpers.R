as_binary <- function(y_true, positive = "ATB") {
  if (is.numeric(y_true)) return(as.integer(y_true > 0))
  as.integer(as.character(y_true) == positive)
}

class_weights <- function(y_true, positive = "ATB") {
  y <- as_binary(y_true, positive = positive)
  if (length(y) == 0 || length(unique(y[!is.na(y)])) < 2) return(rep(1, length(y)))
  n_pos <- sum(y == 1, na.rm = TRUE)
  n_neg <- sum(y == 0, na.rm = TRUE)
  w <- ifelse(y == 1, length(y) / (2 * n_pos), length(y) / (2 * n_neg))
  w[!is.finite(w) | is.na(w)] <- 1
  w
}

compute_roc_auc <- function(y_true, score) {
  y <- as_binary(y_true)
  keep <- !is.na(y) & !is.na(score) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (length(y) == 0) return(NA_real_)
  if (length(unique(y)) < 2) return(NA_real_)
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(tryCatch(as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE, direction = "<"))), error = function(e) NA_real_))
  }
  r <- rank(score)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

compute_pr_auc <- function(y_true, score) {
  y <- as_binary(y_true)
  keep <- !is.na(y) & !is.na(score) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (length(y) == 0) return(NA_real_)
  if (length(unique(y)) < 2) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  recall <- tp / sum(y == 1)
  precision <- tp / pmax(tp + fp, 1)
  recall <- c(0, recall)
  precision <- c(1, precision)
  sum(diff(recall) * (head(precision, -1) + tail(precision, -1)) / 2)
}

compute_binary_metrics <- function(y_true, score, threshold = NULL, positive = "ATB") {
  y <- as_binary(y_true, positive = positive)
  keep <- !is.na(y) & !is.na(score) & is.finite(score)
  y <- y[keep]
  score <- score[keep]
  if (is.null(threshold)) threshold <- 0.5
  pred <- as.integer(score >= threshold)
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA_real_)
  ppv <- ifelse(tp + fp > 0, tp / (tp + fp), NA_real_)
  npv <- ifelse(tn + fn > 0, tn / (tn + fn), NA_real_)
  data.frame(
    n = length(y),
    n_pos = sum(y == 1),
    prevalence = mean(y == 1),
    threshold = threshold,
    AUROC = compute_roc_auc(y, score),
    AUPRC = compute_pr_auc(y, score),
    sensitivity = sensitivity,
    specificity = specificity,
    PPV = ppv,
    NPV = npv,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    Brier = mean((score - y)^2),
    stringsAsFactors = FALSE
  )
}

summarize_per_study_metrics <- function(per_sample_predictions, group_col) {
  groups <- sort(unique(per_sample_predictions[[group_col]]))
  do.call(rbind, lapply(groups, function(g) {
    dat <- per_sample_predictions[per_sample_predictions[[group_col]] == g, , drop = FALSE]
    out <- compute_binary_metrics(dat$true_label, dat$score)
    out[[group_col]] <- g
    out
  }))
}
