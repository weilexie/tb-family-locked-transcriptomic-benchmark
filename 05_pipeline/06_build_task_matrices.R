master <- read_tsv("data_meta/sample_harmonization_master.tsv")
tasks <- list(
  primary_ATB_vs_OD = list(flag = "inclusion_flag_primary_ATB_vs_OD", label = "harmonized_label_binary_primary_ATB_vs_OD"),
  secondary_ATB_vs_LTBI = list(flag = "inclusion_flag_secondary_ATB_vs_LTBI", label = "harmonized_label_binary_secondary_ATB_vs_LTBI"),
  expanded_ATB_vs_nonATB = list(flag = "inclusion_flag_expanded_ATB_vs_nonATB", label = "harmonized_label_binary_expanded_ATB_vs_nonATB")
)

for (task in names(tasks)) {
  message(sprintf("Building task matrix for %s", task))
  spec <- tasks[[task]]
  meta <- master[master[[spec$flag]] == TRUE & !is.na(master[[spec$label]]) & nzchar(master[[spec$label]]), , drop = FALSE]
  meta$task_label <- meta[[spec$label]]
  studies <- sort(unique(meta$study))
  available <- studies[file.exists(project_path("data_curated", "by_study", studies, "expression_gene_symbol.tsv.gz"))]
  meta <- meta[meta$study %in% available, , drop = FALSE]
  X <- build_task_expression_matrix(
    available,
    meta,
    gene_presence_min_fraction = 1 - (read_config("config/config.yml")$expression$missing_gene_max_fraction_training %||% 0.20)
  )
  meta <- meta[match(colnames(X), meta$standard_sample_id), , drop = FALSE]
  out_dir <- file.path("data_curated", "task_matrices", task)
  dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)
  write_tsv_atomic(meta, file.path(out_dir, "metadata.tsv"))
  gzip_tsv_matrix(X, file.path(out_dir, "X_raw_gene_by_sample.tsv.gz"))
  gene_presence <- sapply(rownames(X), function(g) {
    by_study <- tapply(!is.na(X[g, ]), meta$study, any)
    sum(by_study, na.rm = TRUE)
  })
  write_tsv_atomic(data.frame(
    gene = rownames(X),
    n_studies_measured = as.integer(gene_presence),
    n_studies_total = length(unique(meta$study)),
    fraction_studies_measured = as.numeric(gene_presence) / length(unique(meta$study)),
    stringsAsFactors = FALSE
  ), file.path(out_dir, "gene_universe.tsv"))
}
