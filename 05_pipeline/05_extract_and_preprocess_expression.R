config <- read_config("config/config.yml")
ds <- read_tsv("results/tables/DataSummary_full.tsv")
studies <- unique(ds$GEOAccession)

qc_rows <- list()
failures <- list()
for (study in studies) {
  out_path <- project_path("data_curated", "by_study", study, "expression_gene_symbol.tsv.gz")
  if (file.exists(out_path)) {
    log_path <- project_path("data_curated", "by_study", study, "preprocessing_log.tsv")
    if (file.exists(log_path)) {
      qc_rows[[length(qc_rows) + 1]] <- read_tsv(file.path("data_curated", "by_study", study, "preprocessing_log.tsv"))
      next
    }
  }
  message(sprintf("Preprocessing expression for %s", study))
  study_row <- ds[match(study, ds$GEOAccession), , drop = FALSE]
  res <- tryCatch(
    preprocess_one_study_expression(study, config, study_row),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    failures[[length(failures) + 1]] <- data.frame(study = study, error_message = conditionMessage(res), stringsAsFactors = FALSE)
  } else {
    qc_rows[[length(qc_rows) + 1]] <- res
  }
}

qc <- if (length(qc_rows)) do.call(rbind, qc_rows) else data.frame()
write_tsv_atomic(qc, "results/tables/expression_qc_by_study.tsv")
coverage <- if (nrow(qc) > 0) {
  data.frame(
    study = qc$study,
    n_genes = qc$n_genes_after_missing_filter,
    n_samples = qc$n_samples_final,
    assay_name = qc$assay_name,
    log2_decision = qc$log2_decision,
    stringsAsFactors = FALSE
  )
} else {
  data.frame(study = character(), n_genes = integer(), n_samples = integer(), assay_name = character(), log2_decision = character())
}
write_tsv_atomic(coverage, "results/tables/gene_symbol_coverage_by_study.tsv")

if (length(failures)) {
  write_tsv_atomic(do.call(rbind, failures), "results/tables/study_download_failures.tsv")
} else {
  write_tsv_atomic(data.frame(study = character(), error_message = character()), "results/tables/study_download_failures.tsv")
}
