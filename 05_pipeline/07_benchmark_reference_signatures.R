config <- read_config("config/config.yml")
registry <- load_signature_registry()
if (nrow(registry) == 0) {
  stop("No reference signature registry available. Install TBSignatureProfiler or fill config/reference_signature_formula_registry.yml.", call. = FALSE)
}
write_tsv_atomic(registry, "results/benchmark/signature_registry_resolved.tsv")

family_map <- read_tsv("data_meta/cohort_family_map.tsv")
overlap <- flag_signature_training_overlap(registry, family_map)
write_tsv_atomic(overlap, "results/benchmark/signature_training_overlap_by_study.tsv")

task_names <- c("primary_ATB_vs_OD", "secondary_ATB_vs_LTBI", "expanded_ATB_vs_nonATB")
coverage_rows <- list()
score_rows <- list()
study_perf_rows <- list()
family_perf_rows <- list()
orientation_rows <- list()

for (task in task_names) {
  message(sprintf("Benchmarking reference signatures for %s", task))
  meta <- read_tsv(file.path("data_curated/task_matrices", task, "metadata.tsv"))
  xdf <- read_tsv(file.path("data_curated/task_matrices", task, "X_raw_gene_by_sample.tsv.gz"))
  genes <- xdf$gene
  X <- as.matrix(xdf[, setdiff(names(xdf), "gene"), drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- genes
  X_rank <- normalize_samplewise_rank(X)
  colnames(X_rank) <- colnames(X)

  lock <- if (file.exists(project_path("config/external_validation_lockfile.json"))) read_external_lockfile("config/external_validation_lockfile.json") else NULL
  dev_families <- if (!is.null(lock)) unlist(lock$development_families) else unique(meta$cohort_family)
  dev_samples <- meta$standard_sample_id[meta$cohort_family %in% dev_families]
  dev_samples <- intersect(dev_samples, colnames(X_rank))
  y_dev <- meta$task_label[match(dev_samples, meta$standard_sample_id)]

  for (sig in registry$signature_id) {
    sig_genes <- get_signature_genes(sig)
    covered <- intersect(sig_genes$genes, rownames(X_rank))
    coverage <- if (length(sig_genes$genes) == 0) 0 else length(covered) / length(sig_genes$genes)
    reg_row <- registry[registry$signature_id == sig, , drop = FALSE]
    signature_level <- ifelse(coverage < (config$reference_signature_benchmark$minimum_gene_coverage_primary %||% 0.60), "Level 3", reg_row$signature_level)
    primary_eligible <- signature_level %in% c("Level 0", "Level 1") && coverage >= (config$reference_signature_benchmark$minimum_gene_coverage_primary %||% 0.60)
    coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
      task = task,
      signature_id = sig,
      n_signature_genes = length(sig_genes$genes),
      n_covered = length(covered),
      coverage = coverage,
      signature_level = signature_level,
      primary_gene_coverage_eligible = primary_eligible,
      stringsAsFactors = FALSE
    )
    if (length(covered) == 0) next
    if (sig_genes$has_direction) {
      raw_score <- score_signature_signed_mean(X_rank, sig_genes$up_genes, sig_genes$down_genes)
      scoring_method <- "signed_mean_rank"
    } else {
      raw_score <- colMeans(X_rank[covered, , drop = FALSE], na.rm = TRUE)
      scoring_method <- "unsigned_mean_rank_exploratory"
    }
    names(raw_score) <- colnames(X_rank)
    orient <- orient_signature_score(raw_score[dev_samples], y_dev, original_direction = NULL)
    oriented_score <- raw_score * orient$orientation
    orientation_rows[[length(orientation_rows) + 1]] <- data.frame(
      task = task,
      signature_id = sig,
      orientation = orient$orientation,
      orientation_source = "development_families_only",
      scoring_method = scoring_method,
      stringsAsFactors = FALSE
    )
    pred <- data.frame(
      task = task,
      signature_id = sig,
      standard_sample_id = names(oriented_score),
      score = as.numeric(oriented_score),
      scoring_method = scoring_method,
      stringsAsFactors = FALSE
    )
    pred <- merge(pred, meta[, c("standard_sample_id", "study", "cohort_family", "task_label")], by = "standard_sample_id", all.x = TRUE)
    names(pred)[names(pred) == "task_label"] <- "true_label"
    score_rows[[length(score_rows) + 1]] <- pred
    sp <- summarize_per_study_metrics(pred, "study")
    sp$task <- task
    sp$signature_id <- sig
    sp$signature_level <- signature_level
    sp$coverage <- coverage
    study_perf_rows[[length(study_perf_rows) + 1]] <- sp
    fp <- summarize_per_study_metrics(pred, "cohort_family")
    fp$task <- task
    fp$signature_id <- sig
    fp$signature_level <- signature_level
    fp$coverage <- coverage
    family_perf_rows[[length(family_perf_rows) + 1]] <- fp
  }
}

coverage_table <- do.call(rbind, coverage_rows)
scores <- do.call(rbind, score_rows)
study_perf <- do.call(rbind, study_perf_rows)
family_perf <- do.call(rbind, family_perf_rows)
orientation_table <- do.call(rbind, orientation_rows)

write_tsv_atomic(coverage_table, "results/benchmark/signature_gene_coverage.tsv")
write_tsv_atomic(scores, "results/benchmark/signature_scores_by_sample.tsv.gz")
write_tsv_atomic(study_perf, "results/benchmark/signature_per_study_performance.tsv")
write_tsv_atomic(family_perf, "results/benchmark/signature_per_family_performance.tsv")
write_tsv_atomic(orientation_table, "results/benchmark/signature_score_orientation.tsv")

overlap_keys <- unique(overlap[, c("signature_id", "cohort_family")])
family_perf$overlap_key <- paste(family_perf$signature_id, family_perf$cohort_family, sep = "__")
overlap_key <- paste(overlap_keys$signature_id, overlap_keys$cohort_family, sep = "__")
family_perf$include_in_primary_summary <- family_perf$signature_level %in% c("Level 0", "Level 1") &
  family_perf$coverage >= (config$reference_signature_benchmark$minimum_gene_coverage_primary %||% 0.60) &
  !family_perf$overlap_key %in% overlap_key

primary <- family_perf[family_perf$include_in_primary_summary, , drop = FALSE]
summary_primary <- aggregate(cbind(AUROC, AUPRC, balanced_accuracy) ~ task + signature_id, primary, median, na.rm = TRUE)
summary_primary$n_families <- aggregate(AUROC ~ task + signature_id, primary, function(x) sum(!is.na(x)))$AUROC
summary_primary <- summary_primary[order(summary_primary$task, -summary_primary$AUROC, -summary_primary$AUPRC), ]
write_tsv_atomic(summary_primary, "results/benchmark/signature_summary_primary_no_discovery_overlap.tsv")

if (!is.null(lock)) {
  development_families <- unlist(lock$development_families)
  development_primary <- primary[primary$cohort_family %in% development_families, , drop = FALSE]
} else {
  development_primary <- primary
}
summary_development <- aggregate(cbind(AUROC, AUPRC, balanced_accuracy) ~ task + signature_id, development_primary, median, na.rm = TRUE)
summary_development$n_families <- aggregate(AUROC ~ task + signature_id, development_primary, function(x) sum(!is.na(x)))$AUROC
summary_development <- summary_development[order(summary_development$task, -summary_development$AUROC, -summary_development$AUPRC), ]
summary_development$selection_scope <- "development_families_only"
write_tsv_atomic(summary_development, "results/benchmark/signature_summary_development_no_discovery_overlap.tsv")
write_tsv_atomic(summary_development, "results/benchmark/benchmark_rank_table_development.tsv")

summary_all <- aggregate(cbind(AUROC, AUPRC, balanced_accuracy) ~ task + signature_id, family_perf, median, na.rm = TRUE)
summary_all$n_families <- aggregate(AUROC ~ task + signature_id, family_perf, function(x) sum(!is.na(x)))$AUROC
summary_all <- summary_all[order(summary_all$task, -summary_all$AUROC, -summary_all$AUPRC), ]
write_tsv_atomic(summary_all, "results/benchmark/signature_summary_all_studies.tsv")
write_tsv_atomic(summary_primary, "results/benchmark/benchmark_rank_table.tsv")

top_primary_ATB_vs_OD <- head(summary_primary[summary_primary$task == "primary_ATB_vs_OD", ], 10)
write_tsv_atomic(data.frame(
  item = c("scoring_note", "primary_summary_rule", "top_primary_ATB_vs_OD_rows"),
  value = c(
    "Reference signatures scored using TBSignatureProfiler gene sets when available, with manual rank-based signed/unsigned fallback.",
    "Primary summaries exclude discovery/training family overlap and signatures below the configured gene coverage threshold.",
    nrow(top_primary_ATB_vs_OD)
  ),
  stringsAsFactors = FALSE
), "report/benchmark_interpretation.tsv")
