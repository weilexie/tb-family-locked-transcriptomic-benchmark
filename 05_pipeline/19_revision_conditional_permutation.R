options(stringsAsFactors = FALSE)
root <- normalizePath(Sys.getenv("TBML_PROJECT_ROOT", unset = getwd()), winslash = "/")
Sys.setenv(TBML_PROJECT_ROOT = root)
source(file.path(root, "04_R_functions", "00_load_all.R"))

out_dir <- file.path(root, "07_results", "04_revision", "statistics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lock_path <- file.path(root, "01_config", "external_validation_lockfile.json")
lock_hash <- unname(tools::md5sum(lock_path))
set.seed(20261103)
n_perm <- 10000L

pred <- read_tsv(file.path(root, "07_results", "03_ml_final", "external_predictions.tsv.gz"))
d <- pred[pred$task == "primary_ATB_vs_OD" & pred$selected_role == "best_performance", , drop = FALSE]
stopifnot(nrow(d) == 827L, length(unique(d$cohort_family)) == 6L,
          !anyDuplicated(d$standard_sample_id))
y <- as_binary(d$true_label)
score <- as.numeric(d$score)
family <- as.character(d$cohort_family)
observed <- compute_roc_auc(y, score)
groups <- split(seq_along(score), family)

null_auc <- replicate(n_perm, {
  permuted <- score
  for (idx in groups) permuted[idx] <- sample(permuted[idx], replace = FALSE)
  compute_roc_auc(y, permuted)
})
exact <- read_tsv(file.path(out_dir, "conditional_permutation_expectation.tsv"))
exact <- exact[exact$task == "primary_ATB_vs_OD" & exact$selected_role == "best_performance", , drop = FALSE]
stopifnot(nrow(exact) == 1L)

summary <- data.frame(
  task = "primary_ATB_vs_OD",
  selected_role = "best_performance",
  model = unique(d$selector_classifier),
  n_samples = nrow(d),
  n_families = length(groups),
  n_permutations = n_perm,
  seed = 20261103L,
  observed_AUROC = observed,
  conditional_expected_AUROC = exact$conditional_permutation_expected_AUROC,
  observed_minus_conditional_expected = observed - exact$conditional_permutation_expected_AUROC,
  simulated_null_mean_AUROC = mean(null_auc),
  simulated_null_median_AUROC = median(null_auc),
  null_q025 = unname(stats::quantile(null_auc, 0.025)),
  null_q975 = unname(stats::quantile(null_auc, 0.975)),
  empirical_one_sided_p = (1 + sum(null_auc >= observed)) / (n_perm + 1),
  interpretation = "Conditional score-label randomization within the six fixed external families; tests within-family alignment, not AUROC versus 0.5 or causal attribution.",
  stringsAsFactors = FALSE
)
write_tsv_atomic(summary, file.path(out_dir, "conditional_permutation_test.tsv"))
write_tsv_atomic(data.frame(iteration = seq_len(n_perm), permuted_AUROC = null_auc),
                 gzfile(file.path(out_dir, "conditional_permutation_null_distribution.tsv.gz")))
stopifnot(identical(lock_hash, unname(tools::md5sum(lock_path))))
print(summary)
