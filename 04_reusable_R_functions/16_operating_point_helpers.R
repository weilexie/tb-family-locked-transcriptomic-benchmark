mlfinal_select_thresholds <- function(y_true, score) {
  keep <- is.finite(score) & !is.na(score)
  y_true <- y_true[keep]
  score <- score[keep]
  ths <- sort(unique(c(0, 0.5, 1, score)), decreasing = TRUE)
  perf <- do.call(rbind, lapply(ths, function(th) {
    m <- compute_binary_metrics(y_true, score, threshold = th)
    m$threshold <- th
    m
  }))
  pick <- function(name, idx, fallback_col) {
    if (length(idx) == 0) idx <- which.max(perf[[fallback_col]])
    perf[idx[which.max(perf$balanced_accuracy[idx])], , drop = FALSE]
  }
  fixed <- perf[which.min(abs(perf$threshold - 0.5)), , drop = FALSE]
  fixed$threshold_name <- "fixed_0.50"
  youden <- perf[which.max(perf$sensitivity + perf$specificity - 1), , drop = FALSE]
  youden$threshold_name <- "development_youden"
  sens90 <- pick("development_sensitivity_0.90", which(perf$sensitivity >= 0.90), "sensitivity")
  sens90$threshold_name <- "development_sensitivity_0.90"
  spec70 <- pick("development_specificity_0.70", which(perf$specificity >= 0.70), "specificity")
  spec70$threshold_name <- "development_specificity_0.70"
  spec80 <- pick("development_specificity_0.80", which(perf$specificity >= 0.80), "specificity")
  spec80$threshold_name <- "development_specificity_0.80"
  out <- rbind(fixed, youden, sens90, spec70, spec80)
  out[, c("threshold_name", "threshold", "n", "n_pos", "AUROC", "AUPRC", "sensitivity", "specificity", "PPV", "NPV", "balanced_accuracy", "Brier"), drop = FALSE]
}

mlfinal_bootstrap_ci_rows <- function(pred, id_cols, n_boot = 500, seed = 20260502) {
  set.seed(seed)
  groups <- split(pred, pred[[id_cols[1]]])
  rows <- list()
  for (nm in names(groups)) {
    d <- groups[[nm]]
    y <- d$true_label
    score <- d$score
    n <- length(score)
    auc_vals <- replicate(n_boot, {
      idx <- sample(seq_len(n), replace = TRUE)
      compute_roc_auc(y[idx], score[idx])
    })
    auprc_vals <- replicate(n_boot, {
      idx <- sample(seq_len(n), replace = TRUE)
      compute_pr_auc(y[idx], score[idx])
    })
    rows[[length(rows) + 1]] <- data.frame(
      model_id = nm,
      AUROC = compute_roc_auc(y, score),
      AUROC_ci_low = as.numeric(stats::quantile(auc_vals, 0.025, na.rm = TRUE)),
      AUROC_ci_high = as.numeric(stats::quantile(auc_vals, 0.975, na.rm = TRUE)),
      AUPRC = compute_pr_auc(y, score),
      AUPRC_ci_low = as.numeric(stats::quantile(auprc_vals, 0.025, na.rm = TRUE)),
      AUPRC_ci_high = as.numeric(stats::quantile(auprc_vals, 0.975, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
