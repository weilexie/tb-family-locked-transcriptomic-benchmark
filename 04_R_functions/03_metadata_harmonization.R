first_present_value <- function(row, candidates) {
  cols <- intersect(candidates, names(row))
  if (length(cols) == 0) return(NA_character_)
  vals <- as.character(unlist(row[cols], use.names = FALSE))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) == 0) NA_character_ else vals[[1]]
}

classify_tb_label <- function(row) {
  text_cols <- grep("status|diagnosis|disease|group|class|control|pneumonia|tb", names(row), ignore.case = TRUE, value = TRUE)
  text <- tolower(paste(as.character(unlist(row[text_cols], use.names = FALSE)), collapse = " "))
  tbstatus <- tolower(first_present_value(row, c("TBStatus", "tb_status", "TB_status")) %||% "")
  pneumonia <- tolower(first_present_value(row, c("PneumoniaStatus", "pneumonia_status")) %||% "")
  if (grepl("ptb|active tb|active tuberculosis|pulmonary tuberculosis|tb disease|\\batb\\b", tbstatus) ||
      grepl("ptb|active tb|active tuberculosis|pulmonary tuberculosis|tb disease|\\batb\\b", text)) {
    return(list(label = "ATB", confidence = "high", rule = "label_regex_atb"))
  }
  if (grepl("ltbi|latent", tbstatus) || grepl("ltbi|latent", text)) {
    return(list(label = "LTBI", confidence = "high", rule = "label_regex_ltbi"))
  }
  if (grepl("positive", pneumonia) || grepl("pneumonia|other disease|other respiratory|\\bod\\b|copd|asthma|sarcoidosis", text)) {
    return(list(label = "OD", confidence = "high", rule = "label_regex_od"))
  }
  if (grepl("control|healthy|uninfected|negative|no tb|normal", tbstatus) ||
      grepl("healthy control|uninfected control|\\bhc\\b|control", text)) {
    return(list(label = "HC", confidence = "high", rule = "label_regex_hc"))
  }
  list(label = "OTHER_EXCLUDE", confidence = "low", rule = "no_high_confidence_label")
}

classify_treatment_status <- function(row) {
  time_cols <- grep("time|treat|week|month|visit|follow", names(row), ignore.case = TRUE, value = TRUE)
  if (length(time_cols) == 0) {
    return(list(status = "pretreatment_baseline", evidence = "no_timepoint_column_curated_metadata_assumed_baseline_needs_review", review = TRUE))
  }
  text <- tolower(paste(as.character(unlist(row[time_cols], use.names = FALSE)), collapse = " "))
  text <- gsub("\\bna\\b|nan|unknown", "", text)
  if (!nzchar(trimws(text))) {
    return(list(status = "pretreatment_baseline", evidence = "empty_timepoint_value_curated_metadata_assumed_baseline_needs_review", review = TRUE))
  }
  if (grepl("baseline|pre.?treatment|diagnosis|enrollment|day 0|week 0|month 0|\\b0\\b", text)) {
    return(list(status = "pretreatment_baseline", evidence = "explicit_baseline_or_pretreatment_timepoint", review = FALSE))
  }
  if (grepl("post|after|follow|week|month|treatment|end", text)) {
    return(list(status = "on_treatment_or_post_treatment", evidence = "nonbaseline_timepoint_detected", review = FALSE))
  }
  list(status = "unknown", evidence = "timepoint_present_but_unmapped", review = TRUE)
}

standardize_column_data <- function(column_data, study, study_row) {
  df <- as.data.frame(column_data, stringsAsFactors = FALSE)
  sample_ids <- rownames(df)
  if (is.null(sample_ids) || any(!nzchar(sample_ids))) sample_ids <- paste0(study, "_sample_", seq_len(nrow(df)))
  df$original_sample_id <- sample_ids
  out_rows <- lapply(seq_len(nrow(df)), function(i) {
    row <- df[i, , drop = FALSE]
    label <- classify_tb_label(row)
    tx <- classify_treatment_status(row)
    harmonized <- label$label
    primary_ATB_vs_OD <- if (harmonized %in% c("ATB", "OD")) harmonized else NA_character_
    secondary_ATB_vs_LTBI <- if (harmonized %in% c("ATB", "LTBI")) harmonized else NA_character_
    expanded_ATB_vs_nonATB <- if (harmonized %in% c("ATB", "OD", "LTBI", "HC")) ifelse(harmonized == "ATB", "ATB", "nonATB") else NA_character_
    include_baseline <- identical(tx$status, "pretreatment_baseline")
    include_label <- label$confidence == "high"
    data.frame(
      original_sample_id = row$original_sample_id,
      standard_sample_id = paste(study, row$original_sample_id, sep = "__"),
      subject_id_if_available = first_present_value(row, c("PatientID", "SubjectID", "subject_id", "patient_id")),
      study = study,
      cohort_family = study,
      platform = study_row$Platform %||% NA_character_,
      platform_type = ifelse(grepl("RNA|seq", study_row$GeneralType %||% "", ignore.case = TRUE), "RNA-seq", "microarray"),
      tissue = first_present_value(row, c("Tissue", "tissue")) %||% study_row$Tissue %||% NA_character_,
      country = study_row[["Country/Region"]] %||% NA_character_,
      geographical_region = study_row[["Country/Region"]] %||% NA_character_,
      age = first_present_value(row, c("Age", "age")),
      age_group = "unknown",
      sex = first_present_value(row, c("Gender", "Sex", "gender", "sex")),
      hiv_status = first_present_value(row, c("HIVStatus", "HIV_status", "hiv_status")),
      diabetes_status = first_present_value(row, c("DiabetesStatus", "diabetes_status")),
      bcg_status_if_available = first_present_value(row, c("BcgVaccinated", "BCG", "bcg")),
      treatment_status = tx$status,
      timepoint = first_present_value(row, grep("time|week|month|visit", names(row), ignore.case = TRUE, value = TRUE)),
      baseline_evidence = tx$evidence,
      original_diagnosis = paste(as.character(unlist(row[grep("status|diagnosis|disease|group|control|pneumonia|tb", names(row), ignore.case = TRUE, value = TRUE)], use.names = FALSE)), collapse = " | "),
      original_TBStatus = first_present_value(row, c("TBStatus", "tb_status")),
      original_control_type = first_present_value(row, c("Control", "control_type", "PneumoniaStatus")),
      harmonized_label_4class = harmonized,
      harmonized_label_binary_primary_ATB_vs_OD = primary_ATB_vs_OD,
      harmonized_label_binary_secondary_ATB_vs_LTBI = secondary_ATB_vs_LTBI,
      harmonized_label_binary_expanded_ATB_vs_nonATB = expanded_ATB_vs_nonATB,
      control_type_detail = ifelse(harmonized == "OD", "OD_or_ORD_from_metadata", ifelse(harmonized == "HC", "HC", "")),
      inclusion_flag_primary_ATB_vs_OD = include_label && include_baseline && harmonized %in% c("ATB", "OD"),
      inclusion_flag_secondary_ATB_vs_LTBI = include_label && include_baseline && harmonized %in% c("ATB", "LTBI"),
      inclusion_flag_expanded_ATB_vs_nonATB = include_label && include_baseline && harmonized %in% c("ATB", "OD", "LTBI", "HC"),
      exclusion_reason = ifelse(include_label && include_baseline, "", paste(c(if (!include_label) "low_confidence_label" else NULL, if (!include_baseline) "not_baseline" else NULL), collapse = ";")),
      mapping_rule_id = label$rule,
      mapping_confidence = label$confidence,
      manual_override_flag = FALSE,
      source_metadata_columns_used = paste(names(df), collapse = ";"),
      notes = ifelse(tx$review, "baseline assumption needs review", ""),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out_rows)
}
