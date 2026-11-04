if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
  stop("ExperimentHub is required for sample metadata extraction.", call. = FALSE)
}
config <- read_config("config/config.yml")
ds <- if (file.exists(project_path("results/tables/DataSummary_full.tsv"))) read_tsv("results/tables/DataSummary_full.tsv") else load_datasummary()
eh <- ExperimentHub::ExperimentHub()
q <- AnnotationHub::query(eh, "curatedTBData")
md <- as.data.frame(S4Vectors::mcols(q))

all_meta <- list()
failures <- list()
for (study in unique(ds$GEOAccession)) {
  title <- paste0(study, "_column_data")
  hit <- rownames(md)[md$title == title]
  if (length(hit) == 0) {
    failures[[length(failures) + 1]] <- data.frame(study = study, step = "column_data_lookup", error_message = "resource_not_found")
    next
  }
  column_data <- tryCatch(q[[hit[[1]]]], error = function(e) e)
  if (inherits(column_data, "error")) {
    failures[[length(failures) + 1]] <- data.frame(study = study, step = "column_data_download", error_message = conditionMessage(column_data))
    next
  }
  study_row <- ds[match(study, ds$GEOAccession), , drop = FALSE]
  all_meta[[length(all_meta) + 1]] <- standardize_column_data(column_data, study, study_row)
}

master <- if (length(all_meta)) do.call(rbind, all_meta) else data.frame()
write_tsv_atomic(master, "data_meta/sample_harmonization_master.tsv")

needs_review <- master[
  master$mapping_confidence != "high" | grepl("needs_review|assumed", master$baseline_evidence, ignore.case = TRUE),
  ,
  drop = FALSE
]
if (nrow(needs_review) > 0) {
  review <- data.frame(
    study = needs_review$study,
    sample_id = needs_review$standard_sample_id,
    original_columns = needs_review$source_metadata_columns_used,
    suspected_label = needs_review$harmonized_label_4class,
    suspected_treatment_status = needs_review$treatment_status,
    reason_for_uncertainty = paste(needs_review$mapping_confidence, needs_review$baseline_evidence, needs_review$notes, sep = " | "),
    current_action = ifelse(needs_review$inclusion_flag_primary_ATB_vs_OD | needs_review$inclusion_flag_secondary_ATB_vs_LTBI | needs_review$inclusion_flag_expanded_ATB_vs_nonATB, "provisionally_include", "exclude"),
    how_to_override = "Edit config/manual_sample_exclusions.yml or add a future manual phenotype override table.",
    stringsAsFactors = FALSE
  )
} else {
  review <- data.frame(study = character(), sample_id = character(), original_columns = character(), suspected_label = character(), suspected_treatment_status = character(), reason_for_uncertainty = character(), current_action = character(), how_to_override = character())
}
write_tsv_atomic(review, "results/tables/phenotype_needs_manual_review.tsv")

study_counts <- aggregate(
  cbind(inclusion_flag_primary_ATB_vs_OD, inclusion_flag_secondary_ATB_vs_LTBI, inclusion_flag_expanded_ATB_vs_nonATB) ~ study,
  data = transform(master, inclusion_flag_primary_ATB_vs_OD = as.integer(inclusion_flag_primary_ATB_vs_OD), inclusion_flag_secondary_ATB_vs_LTBI = as.integer(inclusion_flag_secondary_ATB_vs_LTBI), inclusion_flag_expanded_ATB_vs_nonATB = as.integer(inclusion_flag_expanded_ATB_vs_nonATB)),
  sum
)
label_counts <- as.data.frame.matrix(table(master$study, master$harmonized_label_4class))
label_counts$study <- rownames(label_counts)
study_counts <- merge(study_counts, label_counts, by = "study", all = TRUE)
write_tsv_atomic(study_counts, "results/tables/study_inclusion_exclusion_preliminary.tsv")

if (length(failures)) {
  write_tsv_atomic(do.call(rbind, failures), "results/tables/sample_metadata_download_failures.tsv")
} else {
  write_tsv_atomic(data.frame(study = character(), step = character(), error_message = character()), "results/tables/sample_metadata_download_failures.tsv")
}

