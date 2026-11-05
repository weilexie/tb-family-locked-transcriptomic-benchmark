mlaudit_clamp01 <- function(x, eps = 1e-6) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

mlaudit_auc_row <- function(y_true, score) {
  y <- as_binary(y_true)
  keep <- !is.na(y) & is.finite(score)
  y <- y[keep]
  score <- as.numeric(score[keep])
  data.frame(
    n = length(y),
    n_pos = sum(y == 1),
    prevalence = ifelse(length(y) > 0, mean(y == 1), NA_real_),
    AUROC = compute_roc_auc(y, score),
    AUPRC = compute_pr_auc(y, score),
    stringsAsFactors = FALSE
  )
}

mlaudit_fit_platt <- function(y_true, score) {
  y <- as_binary(y_true)
  x <- stats::qlogis(mlaudit_clamp01(score))
  keep <- !is.na(y) & is.finite(x)
  y <- y[keep]
  x <- x[keep]
  if (length(y) < 10 || length(unique(y)) < 2 || length(unique(x)) < 2) {
    return(list(ok = FALSE, intercept = NA_real_, slope = NA_real_, fit = NULL))
  }
  fit <- try(stats::glm(y ~ x, family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(ok = FALSE, intercept = NA_real_, slope = NA_real_, fit = NULL))
  }
  co <- stats::coef(fit)
  list(ok = TRUE, intercept = unname(co[1]), slope = unname(co[2]), fit = fit)
}

mlaudit_predict_platt <- function(platt, score) {
  if (!isTRUE(platt$ok) || is.null(platt$fit)) return(mlaudit_clamp01(score))
  x <- stats::qlogis(mlaudit_clamp01(score))
  pred <- as.numeric(stats::predict(platt$fit, newdata = data.frame(x = x), type = "response"))
  pred[!is.finite(pred)] <- NA_real_
  pred
}

mlaudit_cramers_v <- function(x, y) {
  tab <- table(x, y, useNA = "no")
  if (any(dim(tab) < 2) || sum(tab) == 0) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

mlaudit_score_label_permutation <- function(pred, n_perm = 100, seed = 20260502) {
  observed <- compute_roc_auc(pred$true_label, pred$score)
  set.seed(seed)
  vals <- replicate(n_perm, compute_roc_auc(sample(pred$true_label), pred$score))
  data.frame(
    observed_AUROC = observed,
    n_permutations = n_perm,
    permuted_AUROC_median = stats::median(vals, na.rm = TRUE),
    permuted_AUROC_q025 = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE)),
    permuted_AUROC_q975 = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE)),
    empirical_p_ge_observed = (1 + sum(vals >= observed, na.rm = TRUE)) / (n_perm + 1),
    stringsAsFactors = FALSE
  )
}
