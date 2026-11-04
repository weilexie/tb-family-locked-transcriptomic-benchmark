load_signature_registry <- function(path = "config/reference_signature_formula_registry.yml") {
  if (requireNamespace("TBSignatureProfiler", quietly = TRUE)) {
    env <- new.env(parent = emptyenv())
    utils::data(list = "TBsignatures", package = "TBSignatureProfiler", envir = env)
    utils::data(list = "TBsignaturesSplit", package = "TBSignatureProfiler", envir = env)
    utils::data(list = "sigAnnotData", package = "TBSignatureProfiler", envir = env)
    signatures <- get("TBsignatures", envir = env)
    split <- get("TBsignaturesSplit", envir = env)
    annot <- get("sigAnnotData", envir = env)
    rows <- lapply(names(signatures), function(sig) {
      split_entry <- split[[sig]]
      has_direction <- is.list(split_entry) && length(split_entry) >= 2
      up <- down <- character()
      if (has_direction) {
        up_name <- grep("_up$", names(split_entry), value = TRUE)
        dn_name <- grep("_dn$|_down$", names(split_entry), value = TRUE)
        up <- unique(unlist(split_entry[up_name], use.names = FALSE))
        down <- unique(unlist(split_entry[dn_name], use.names = FALSE))
      }
      data.frame(
        signature_id = sig,
        source_package = "TBSignatureProfiler",
        source_object = ifelse(has_direction, "TBsignaturesSplit", "TBsignatures"),
        signature_level = ifelse(has_direction, "Level 1", "Level 2"),
        has_formula = FALSE,
        has_coefficients = FALSE,
        has_up_down_direction = has_direction,
        n_up_genes = length(up),
        n_down_genes = length(down),
        n_total_genes = length(unique(signatures[[sig]])),
        intended_task = ifelse(grepl("_OD|OD_", sig), "ATB_vs_OD", "TB_host_response"),
        original_score_direction = "unknown_or_development_oriented",
        training_studies_raw = "",
        training_studies_mapped = "",
        registry_status = "resolved_from_TBSignatureProfiler",
        notes = "",
        stringsAsFactors = FALSE
      )
    })
    registry <- do.call(rbind, rows)
    attr(registry, "TBsignatures") <- signatures
    attr(registry, "TBsignaturesSplit") <- split
    return(registry)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml is required.", call. = FALSE)
  registry <- yaml::read_yaml(project_path(path))
  data.frame(signature_id = character(), source_package = character(), source_object = character())
}

flag_signature_training_overlap <- function(signature_registry, study_family_map) {
  if (is.null(signature_registry)) {
    return(data.frame(
      signature_id = c("synthetic_overlap_signature", "synthetic_independent_signature"),
      study = c("SYN_STUDY_1", "SYN_STUDY_4"),
      cohort_family = c("FAM_1", "FAM_4"),
      is_training_study_overlap = c(TRUE, FALSE),
      is_training_family_overlap = c(TRUE, FALSE),
      include_in_primary_summary = c(FALSE, TRUE),
      reason = c("synthetic discovery overlap", "synthetic no overlap"),
      stringsAsFactors = FALSE
    ))
  }
  sig_train <- load_signature_info_training()
  names(sig_train) <- sub("^TBSignature$", "signature_id", names(sig_train))
  names(sig_train) <- sub("^Study$", "study", names(sig_train))
  merged <- merge(sig_train, study_family_map[, c("study", "cohort_family")], by = "study", all.x = TRUE)
  if (nrow(merged) == 0) {
    return(data.frame(signature_id = character(), study = character(), cohort_family = character(), is_training_study_overlap = logical(), is_training_family_overlap = logical(), include_in_primary_summary = logical(), reason = character()))
  }
  merged$is_training_study_overlap <- TRUE
  merged$is_training_family_overlap <- !is.na(merged$cohort_family)
  merged$include_in_primary_summary <- FALSE
  merged$reason <- "signature discovery/training overlap from SignatureInfoTraining"
  merged
}
