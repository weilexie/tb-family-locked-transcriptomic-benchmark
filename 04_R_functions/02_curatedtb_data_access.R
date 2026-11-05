load_datasummary <- function() {
  if (!requireNamespace("curatedTBData", quietly = TRUE)) {
    stop("curatedTBData is not installed; install Bioconductor curatedTBData before prepare_data.", call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  utils::data(list = "DataSummary", package = "curatedTBData", envir = env)
  obj <- get("DataSummary", envir = env)
  as.data.frame(obj)
}

load_signature_info_training <- function() {
  if (!requireNamespace("curatedTBData", quietly = TRUE)) {
    stop("curatedTBData is not installed; install Bioconductor curatedTBData before prepare_data.", call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  utils::data(list = "SignatureInfoTraining", package = "curatedTBData", envir = env)
  obj <- get("SignatureInfoTraining", envir = env)
  as.data.frame(obj)
}

curatedtb_dry_run <- function(studies, curated_only = FALSE) {
  if (!requireNamespace("curatedTBData", quietly = TRUE)) {
    stop("curatedTBData is not installed; cannot run dry-run.", call. = FALSE)
  }
  rows <- lapply(studies, function(study) {
    resources <- tryCatch(
      curatedTBData::curatedTBData(study, dry.run = TRUE, curated.only = curated_only),
      error = function(e) structure(character(), error_message = conditionMessage(e))
    )
    if (length(resources) == 0) {
      data.frame(study = study, curated_only = curated_only, resource = NA_character_, status = "failed", error_message = attr(resources, "error_message") %||% "no resources", stringsAsFactors = FALSE)
    } else {
      data.frame(study = study, curated_only = curated_only, resource = as.character(resources), status = "ok", error_message = "", stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}

load_curatedtb_study <- function(study, assay_priority, curated_only = TRUE) {
  if (!requireNamespace("curatedTBData", quietly = TRUE)) {
    stop("curatedTBData is not installed.", call. = FALSE)
  }
  obj_list <- curatedTBData::curatedTBData(study, dry.run = FALSE, curated.only = curated_only)
  if (is.list(obj_list) && study %in% names(obj_list)) {
    study_object <- obj_list[[study]]
  } else if (is.list(obj_list) && length(obj_list) >= 1) {
    study_object <- obj_list[[1]]
  } else {
    study_object <- obj_list
  }
  selected <- select_best_assay(study_object, assay_priority)
  expr <- extract_expression_matrix(study_object, selected$assay_name)
  sample_metadata <- extract_sample_metadata(study_object, study)
  feature_metadata <- extract_feature_metadata(study_object, selected$assay_name)
  list(
    study = study,
    assay_name = selected$assay_name,
    expression = expr,
    sample_metadata = sample_metadata,
    feature_metadata = feature_metadata,
    object_class = paste(class(study_object), collapse = ";"),
    notes = selected$notes
  )
}

select_best_assay <- function(mae_or_list, assay_priority) {
  assay_names <- character()
  if (inherits(mae_or_list, "MultiAssayExperiment")) {
    assay_names <- names(MultiAssayExperiment::experiments(mae_or_list))
  } else if (is.list(mae_or_list)) {
    assay_names <- names(mae_or_list)
  }
  if (length(assay_names) == 0) {
    stop("Could not determine assay names from curatedTBData object.", call. = FALSE)
  }
  priority <- unlist(assay_priority, use.names = FALSE)
  selected <- priority[priority %in% assay_names][1]
  if (is.na(selected) || length(selected) == 0) selected <- assay_names[[1]]
  list(
    assay_name = selected,
    notes = sprintf("Selected assay '%s' from available assays: %s", selected, paste(assay_names, collapse = ", "))
  )
}

extract_expression_matrix <- function(study_object, assay_name) {
  if (inherits(study_object, "MultiAssayExperiment")) {
    expr <- MultiAssayExperiment::experiments(study_object)[[assay_name]]
  } else if (is.list(study_object)) {
    expr <- study_object[[assay_name]]
  } else if (is.matrix(study_object)) {
    expr <- study_object
  } else {
    stop(sprintf("Unsupported study object class: %s", paste(class(study_object), collapse = ";")), call. = FALSE)
  }
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr
}

extract_sample_metadata <- function(study_object, study) {
  if (inherits(study_object, "MultiAssayExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(study_object), stringsAsFactors = FALSE)
  } else {
    meta <- data.frame(row.names = colnames(extract_expression_matrix(study_object, names(study_object)[1])))
  }
  meta$original_sample_id <- rownames(meta)
  meta$standard_sample_id <- paste(study, meta$original_sample_id, sep = "__")
  meta$study <- study
  rownames(meta) <- NULL
  meta
}

extract_feature_metadata <- function(study_object, assay_name) {
  expr <- extract_expression_matrix(study_object, assay_name)
  data.frame(
    feature_id = rownames(expr),
    gene_symbol = rownames(expr),
    stringsAsFactors = FALSE
  )
}
