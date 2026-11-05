update_and_clean_gene_symbols <- function(expr, feature_metadata = NULL) {
  keep <- !is.na(rownames(expr)) & nzchar(rownames(expr)) & !grepl("[;/|]", rownames(expr))
  expr[keep, , drop = FALSE]
}

collapse_duplicate_genes <- function(expr, rule = "highest_mean_expression") {
  genes <- rownames(expr)
  if (!anyDuplicated(genes)) {
    return(list(expression = expr, resolution = data.frame(gene = genes, action = "kept_unique")))
  }
  rows <- lapply(unique(genes), function(g) {
    idx <- which(genes == g)
    chosen <- idx[which.max(rowMeans(expr[idx, , drop = FALSE], na.rm = TRUE))]
    list(gene = g, chosen_row = chosen, dropped_rows = paste(setdiff(idx, chosen), collapse = ","))
  })
  chosen <- vapply(rows, `[[`, integer(1), "chosen_row")
  out <- expr[chosen, , drop = FALSE]
  rownames(out) <- genes[chosen]
  list(expression = out, resolution = do.call(rbind, lapply(rows, as.data.frame)))
}

preprocess_one_study_expression <- function(study, config, study_row = NULL) {
  out_dir <- project_path("data_curated", "by_study", study)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  platform_type <- if (!is.null(study_row) && grepl("RNA|seq", study_row$GeneralType %||% "", ignore.case = TRUE)) "rnaseq" else "microarray"
  assay_priority <- config$primary_data_source$assay_priority[[platform_type]] %||% c("assay_curated")
  loaded <- load_curatedtb_study(
    study = study,
    assay_priority = assay_priority,
    curated_only = config$primary_data_source$use_curated_only_first %||% TRUE
  )
  expr0 <- loaded$expression
  original_dim <- dim(expr0)
  colnames(expr0) <- paste(study, colnames(expr0), sep = "__")

  expr1 <- update_and_clean_gene_symbols(expr0, loaded$feature_metadata)
  symbol_log <- data.frame(
    original_features = original_dim[[1]],
    retained_after_symbol_clean = nrow(expr1),
    dropped = original_dim[[1]] - nrow(expr1),
    stringsAsFactors = FALSE
  )
  collapse <- collapse_duplicate_genes(expr1, config$expression$duplicate_gene_rule %||% "highest_mean_expression")
  expr2 <- collapse$expression
  log2_result <- detect_and_log2_transform(
    expr2,
    q99_gt = config$expression$log2_detection$q99_gt,
    max_gt = config$expression$log2_detection$max_gt,
    offset = config$expression$log2_detection$offset
  )
  expr3 <- log2_result$expression
  na_fraction <- rowMeans(is.na(expr3))
  keep <- na_fraction <= (config$expression$missing_gene_max_fraction_training %||% 0.20)
  expr4 <- expr3[keep, , drop = FALSE]
  if (anyNA(expr4)) {
    for (i in which(rowSums(is.na(expr4)) > 0)) {
      expr4[i, is.na(expr4[i, ])] <- stats::median(expr4[i, ], na.rm = TRUE)
    }
  }

  sample_annotation <- loaded$sample_metadata
  sample_annotation$original_sample_id <- sub(paste0("^", study, "__"), "", sample_annotation$standard_sample_id)
  sample_annotation$standard_sample_id <- paste(study, sample_annotation$original_sample_id, sep = "__")

  feature_annotation <- data.frame(gene = rownames(expr4), stringsAsFactors = FALSE)
  qc <- data.frame(
    study = study,
    assay_name = loaded$assay_name,
    object_class = loaded$object_class,
    n_features_original = original_dim[[1]],
    n_samples_original = original_dim[[2]],
    n_genes_after_clean = nrow(expr1),
    n_genes_after_collapse = nrow(expr2),
    n_genes_after_missing_filter = nrow(expr4),
    n_samples_final = ncol(expr4),
    log2_decision = log2_result$decision,
    q99 = log2_result$q99,
    max = log2_result$max,
    has_negative = log2_result$has_negative,
    notes = loaded$notes,
    stringsAsFactors = FALSE
  )

  gzip_tsv_matrix(expr4, file.path("data_curated", "by_study", study, "expression_gene_symbol.tsv.gz"))
  write_tsv_atomic(sample_annotation, file.path("data_curated", "by_study", study, "sample_annotation.tsv"))
  write_tsv_atomic(feature_annotation, file.path("data_curated", "by_study", study, "feature_annotation.tsv"))
  write_tsv_atomic(qc, file.path("data_curated", "by_study", study, "preprocessing_log.tsv"))
  write_tsv_atomic(symbol_log, file.path("data_curated", "by_study", study, "gene_symbol_update_log.tsv"))
  write_tsv_atomic(collapse$resolution, file.path("data_curated", "by_study", study, "duplicate_gene_resolution.tsv"))
  qc
}

read_study_expression <- function(study) {
  path <- project_path("data_curated", "by_study", study, "expression_gene_symbol.tsv.gz")
  if (!file.exists(path)) stop(sprintf("Missing preprocessed expression for %s", study), call. = FALSE)
  df <- read_tsv(file.path("data_curated", "by_study", study, "expression_gene_symbol.tsv.gz"))
  genes <- df$gene
  mat <- as.matrix(df[, setdiff(names(df), "gene"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat
}

build_task_expression_matrix <- function(studies, sample_meta, gene_universe = NULL, gene_presence_min_fraction = 1.0) {
  matrices <- list()
  for (study in studies) {
    mat <- read_study_expression(study)
    keep_samples <- intersect(sample_meta$standard_sample_id[sample_meta$study == study], colnames(mat))
    if (length(keep_samples) == 0) next
    matrices[[study]] <- mat[, keep_samples, drop = FALSE]
  }
  if (length(matrices) == 0) stop("No study matrices available for task.", call. = FALSE)
  if (is.null(gene_universe)) {
    gene_counts <- table(unlist(lapply(matrices, rownames), use.names = FALSE))
    min_studies <- ceiling(length(matrices) * gene_presence_min_fraction)
    gene_universe <- names(gene_counts)[gene_counts >= min_studies]
    gene_universe <- sort(gene_universe)
  }
  combined <- do.call(cbind, lapply(matrices, function(m) {
    out <- matrix(NA_real_, nrow = length(gene_universe), ncol = ncol(m), dimnames = list(gene_universe, colnames(m)))
    present <- intersect(gene_universe, rownames(m))
    out[present, ] <- m[present, , drop = FALSE]
    out
  }))
  combined
}
