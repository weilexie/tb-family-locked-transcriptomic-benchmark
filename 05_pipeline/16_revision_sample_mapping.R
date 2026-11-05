# Reviewer R2.1: distinguish metadata records from selected-assay samples.
# This analysis preserves the submitted matrices, model predictions, and lock.
root <- normalizePath(Sys.getenv("TBML_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = TRUE)
out <- file.path(root, "07_results", "04_revision", "sample_mapping")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tab <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tab <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
master <- read_tab(file.path(root, "02_data_meta", "sample_harmonization_master.tsv"))
table1 <- read_tab(file.path(root, "08_manuscript", "01_tables", "Table1_included_cohorts.tsv"))
lock <- jsonlite::read_json(file.path(root, "01_config", "external_validation_lockfile.json"))
lock_path <- file.path(root, "01_config", "external_validation_lockfile.json")
lock_before <- unname(tools::md5sum(lock_path))
tasks <- c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB")
assay_ids <- character()
assay_by_study <- list()
for (study in unique(master$study)) {
  path <- file.path(root, "03_data_curated", "01_by_study", study, "expression_gene_symbol.tsv.gz")
  if (!file.exists(path)) next
  con <- gzfile(path, "rt")
  header <- readLines(con, n = 1, warn = FALSE)
  close(con)
  ids <- strsplit(header, "\t", fixed = TRUE)[[1]][-1]
  ids <- sub('^"(.*)"$', '\\1', ids)
  assay_by_study[[study]] <- ids
  assay_ids <- c(assay_ids, ids)
}
stopifnot(!anyDuplicated(master$standard_sample_id), !anyDuplicated(assay_ids))

# SampleCode identifies the GSE89403 biological specimen at a particular visit.
# Upstream curatedTBData retains one representative GSM for each DescriptionID.
raw89403 <- read_tab(file.path(root, "03_data_curated", "01_by_study", "GSE89403", "sample_annotation.tsv"))
g <- master[master$study == "GSE89403", , drop = FALSE]
idx <- match(g$standard_sample_id, raw89403$standard_sample_id)
stopifnot(!anyNA(idx))
g$SampleCode <- raw89403$SampleCode[idx]
g$PatientID <- raw89403$PatientID[idx]
g$MeasurementTime <- raw89403$MeasurementTime[idx]
g$represented_in_selected_assay <- g$standard_sample_id %in% assay_by_study[["GSE89403"]]
specimens <- split(seq_len(nrow(g)), g$SampleCode)
specimen_rows <- lapply(specimens, function(ix) {
  z <- g[ix, , drop = FALSE]
  representatives <- z$standard_sample_id[z$represented_in_selected_assay]
  if (length(representatives) != 1) stop("Expected one selected-assay representative per GSE89403 SampleCode.")
  if (length(unique(z$harmonized_label_4class)) != 1 || length(unique(z$treatment_status)) != 1 || length(unique(z$MeasurementTime)) != 1 || length(unique(z$PatientID)) != 1) stop("Discordant technical-replicate annotation.")
  data.frame(SampleCode = z$SampleCode[1], PatientID = z$PatientID[1], MeasurementTime = z$MeasurementTime[1],
    harmonized_label_4class = z$harmonized_label_4class[1], treatment_status = z$treatment_status[1],
    metadata_records_n = nrow(z), selected_assay_representative = representatives,
    primary_task_ready = z$inclusion_flag_primary_ATB_vs_OD[1],
    expanded_task_ready = z$inclusion_flag_expanded_ATB_vs_nonATB[1], stringsAsFactors = FALSE)
})
specimen_table <- do.call(rbind, specimen_rows)
g$selected_assay_representative <- specimen_table$selected_assay_representative[match(g$SampleCode, specimen_table$SampleCode)]
g$metadata_records_per_specimen <- specimen_table$metadata_records_n[match(g$SampleCode, specimen_table$SampleCode)]
g$record_to_assay_reason <- ifelse(g$represented_in_selected_assay, "selected_assay_representative", "technical_replicate_already_collapsed_upstream")
cross_cols <- c("original_sample_id", "standard_sample_id", "SampleCode", "PatientID", "MeasurementTime", "harmonized_label_4class", "treatment_status", "baseline_evidence", "inclusion_flag_primary_ATB_vs_OD", "inclusion_flag_expanded_ATB_vs_nonATB", "represented_in_selected_assay", "selected_assay_representative", "metadata_records_per_specimen", "record_to_assay_reason")
write_tab(g[, cross_cols], "GSE89403_record_to_specimen_crosswalk.tsv")
write_tab(specimen_table, "GSE89403_unique_specimens.tsv")
stage_groups <- split(seq_len(nrow(g)), interaction(g$harmonized_label_4class, g$treatment_status, drop = TRUE))
stage_counts <- do.call(rbind, lapply(stage_groups, function(ix) {
  z <- g[ix, , drop = FALSE]
  data.frame(harmonized_label_4class = z$harmonized_label_4class[1], treatment_status = z$treatment_status[1],
    metadata_records_n = nrow(z), unique_SampleCode_n = length(unique(z$SampleCode)),
    selected_assay_samples_n = sum(z$represented_in_selected_assay),
    primary_eligible_records_n = sum(z$inclusion_flag_primary_ATB_vs_OD),
    expanded_eligible_records_n = sum(z$inclusion_flag_expanded_ATB_vs_nonATB))
}))
write_tab(stage_counts, "GSE89403_flow_by_label_and_treatment.tsv")

flows <- list()
for (task in tasks) {
  flag <- paste0("inclusion_flag_", task)
  taskmeta <- read_tab(file.path(root, "03_data_curated", "02_task_matrices", task, "metadata.tsv"))
  for (family in unique(master$cohort_family)) {
    z <- master[master$cohort_family == family, , drop = FALSE]
    eligible <- z[!is.na(z[[flag]]) & z[[flag]], , drop = FALSE]
    matched <- eligible[eligible$standard_sample_id %in% assay_ids, , drop = FALSE]
    in_matrix <- matched[matched$standard_sample_id %in% taskmeta$standard_sample_id, , drop = FALSE]
    collapsed <- sum(eligible$study == "GSE89403" & !eligible$standard_sample_id %in% assay_ids)
    missing_assay <- nrow(eligible) - nrow(matched) - collapsed
    split <- if (family %in% unlist(lock$external_validation_families)) "locked_external" else if (family %in% unlist(lock$development_families)) "development" else "not_task_eligible"
    flows[[length(flows) + 1]] <- data.frame(task = task, cohort_family = family, analysis_split = split,
      metadata_records_all_visits_n = nrow(z), baseline_metadata_records_n = sum(z$treatment_status == "pretreatment_baseline"),
      baseline_phenotype_eligible_records_n = nrow(eligible), technical_replicate_records_collapsed_upstream_n = collapsed,
      other_eligible_records_without_selected_assay_n = missing_assay, selected_assay_eligible_samples_n = nrow(matched),
      task_matrix_samples_n = nrow(in_matrix), ATB_task_matrix_n = sum(in_matrix$harmonized_label_4class == "ATB"),
      OD_task_matrix_n = sum(in_matrix$harmonized_label_4class == "OD"), LTBI_task_matrix_n = sum(in_matrix$harmonized_label_4class == "LTBI"),
      HC_task_matrix_n = sum(in_matrix$harmonized_label_4class == "HC"),
      count_unit_note = if (family == "GSE89403") "Metadata GSM records include technical replicates; selected-assay baseline samples are unique SampleCode specimens." else "Metadata records and selected-assay samples; subject independence is not inferred from sample IDs alone.", stringsAsFactors = FALSE)
  }
}
flow <- do.call(rbind, flows)
write_tab(flow, "family_task_sample_flow_reconciled.tsv")
totals <- aggregate(flow[, c("metadata_records_all_visits_n", "baseline_metadata_records_n", "baseline_phenotype_eligible_records_n", "technical_replicate_records_collapsed_upstream_n", "other_eligible_records_without_selected_assay_n", "selected_assay_eligible_samples_n", "task_matrix_samples_n", "ATB_task_matrix_n", "OD_task_matrix_n", "LTBI_task_matrix_n", "HC_task_matrix_n")], flow[, c("task", "analysis_split")], sum)
write_tab(totals, "task_split_sample_flow_totals.tsv")

for (task in tasks) {
  f <- flow[flow$task == task, ]
  prefix <- if (task == tasks[1]) "primary" else "expanded"
  for (name in c("metadata_records_all_visits_n", "baseline_metadata_records_n", "baseline_phenotype_eligible_records_n", "technical_replicate_records_collapsed_upstream_n", "other_eligible_records_without_selected_assay_n", "task_matrix_samples_n")) {
    table1[[paste(prefix, name, sep = "_")]] <- f[[name]][match(table1$cohort_family, f$cohort_family)]
  }
}
names(table1)[names(table1) == "n_samples_harmonized"] <- "metadata_records_all_visits_n"
for (label in c("ATB", "OD", "LTBI", "HC")) names(table1)[names(table1) == label] <- paste0(label, "_metadata_records_all_visits_n")
table1$sample_count_note[table1$cohort_family == "GSE89403"] <- "914 GSM metadata records include longitudinal visits and technical replicates. Baseline phenotype-eligible records (216 primary, 266 expanded) correspond to 108 and 129 unique selected-assay specimens; upstream technical-replicate collapse explains the difference. Existing model predictions already use these unique specimens."
write_tab(table1, "Table1_included_cohorts_reconciled.tsv")
unmatched <- master[master$standard_sample_id %in% master$standard_sample_id[master$inclusion_flag_primary_ATB_vs_OD | master$inclusion_flag_expanded_ATB_vs_nonATB] & !master$standard_sample_id %in% assay_ids & master$study != "GSE89403", , drop = FALSE]
write_tab(unmatched, "other_eligible_records_without_selected_assay.tsv")

# Source provenance supports the harmonization description without claiming that
# all studies used a common raw-data normalization or that local ComBat was run.
pkg <- find.package("curatedTBData")
description <- read.dcf(file.path(pkg, "DESCRIPTION"))
upstream <- readLines(file.path(pkg, "scripts", "make-data_GSE89403.R"), warn = FALSE)
source_hits <- grep("technical replicates|distinct\\(ID_table|SYMBOL_NEW|probesetsToGenes|normalized_log2|assay_curated", upstream)
provenance <- data.frame(source = "curatedTBData/scripts/make-data_GSE89403.R", source_line = source_hits, source_text = upstream[source_hits], package_version = as.character(utils::packageVersion("curatedTBData")), package_commit = unname(description[1, "git_last_commit"]), stringsAsFactors = FALSE)
write_tab(provenance, "GSE89403_upstream_processing_evidence.tsv")
study_provenance <- list()
for (study in unique(master$study)) {
  source_path <- file.path(pkg, "scripts", paste0("make-data_", study, ".R"))
  if (!file.exists(source_path)) {
    study_provenance[[length(study_provenance) + 1]] <- data.frame(study = study, source_file = basename(source_path), source_line = NA_integer_, evidence_type = "source_not_packaged", source_text = "Installed package does not contain this study's curation script.")
    next
  }
  lines <- readLines(source_path, warn = FALSE)
  for (pattern in c(mapping = "map_gene_symbol|match_gene_symbol|SYMBOL_NEW|probesetsToGenes", normalization = "rma\\(|normalizeBetweenArrays\\(|data_normalized|normalized_log2|readRDS.*normal|read\\.delim.*normal")) {
    hits <- grep(pattern, lines, ignore.case = TRUE)
    if (length(hits)) study_provenance[[length(study_provenance) + 1]] <- data.frame(study = study, source_file = basename(source_path), source_line = hits, evidence_type = if (grepl("map_gene", pattern)) "mapping" else "normalization_input_or_call", source_text = lines[hits])
  }
}
write_tab(do.call(rbind, study_provenance), "study_upstream_curation_source_evidence.tsv")
local_processing <- data.frame(stage = c("input", "upstream_mapping", "GSE89403_source_scale", "technical_replicates", "local_gene_cleanup", "local_duplicate_genes", "local_missing_values", "task_gene_universe", "model_scaling", "signature_scaling", "ComBat_SVA", "raw_read_processing"),
  implementation = c("curatedTBData 2.4.0 assay_curated, not assay_reprocess; normalized curated public matrices", "Study-specific upstream curation maps probes/features to genes; GSE89403 uses supplied SYMBOL_NEW and probesetsToGenes(..., median). No common annotation release was independently reconstructed by this project.", "GEO supplied normalized_log2 expression; local log2 detection retained already_log (q99=9.71804556, maximum=15.82038167)", "GSE89403 upstream distinct(DescriptionID) retains one representative GSM per expression specimen; all GSM metadata records remain available", "Remove blank or ambiguous symbols containing semicolon, slash, or vertical bar; no local HGNChelper gene-symbol update call", "Retain duplicate gene row with highest whole-study mean expression", "Whole-study gene missingness filter <=20%; fill remaining missing values with gene median within that study", "Genes available in >=80% of studies in the whole task, including locked external studies; this availability filter is not training-fold-only", "Training-only median/MAD frozen scaling, or selected-gene samplewise percentile ranks, depending on the model", "Training-only signature median imputation then mean/SD scaling", "No local ComBat/SVA/removeBatchEffect call was applied; package vignette mentions SVA but this project did not call it", "No project-level raw CEL/FASTQ reprocessing, TMM/DESeq or common quantile-normalization step was applied"), stringsAsFactors = FALSE)
write_tab(local_processing, "harmonization_implementation.tsv")

if ("--verify-source" %in% commandArgs(trailingOnly = TRUE)) {
  source_mae <- curatedTBData::curatedTBData("GSE89403", dry.run = FALSE, curated.only = TRUE)[[1]]
  source_map <- as.data.frame(MultiAssayExperiment::sampleMap(source_mae))
  selected_map <- source_map[as.character(source_map$assay) == "assay_curated", , drop = FALSE]
  source_meta <- as.data.frame(SummarizedExperiment::colData(source_mae))
  stopifnot(nrow(source_meta) == nrow(g), nrow(selected_map) == 453L,
    setequal(paste0("GSE89403__", selected_map$colname), assay_by_study[["GSE89403"]]),
    setequal(rownames(source_meta), g$original_sample_id))
  compare <- match(g$original_sample_id, rownames(source_meta))
  canonical_id <- function(x) { x <- as.character(x); x[is.na(x) | x %in% c("NA", "NaN", "")] <- NA_character_; x }
  stopifnot(identical(as.character(source_meta$SampleCode[compare]), as.character(g$SampleCode)),
    identical(canonical_id(source_meta$PatientID[compare]), canonical_id(g$PatientID)))
  write_tab(source_map, "GSE89403_source_MAE_sampleMap.tsv")
  write_tab(data.frame(check = c("source_metadata_records", "source_curated_assay_samples", "source_all_assay_map_rows", "source_local_selected_sample_IDs_equal", "source_local_SampleCode_and_PatientID_equal"),
    value = c(nrow(source_meta), nrow(selected_map), nrow(source_map), TRUE, TRUE)), "GSE89403_source_MAE_validation.tsv")
}

primary <- specimen_table[specimen_table$primary_task_ready, ]
expanded <- specimen_table[specimen_table$expanded_task_ready, ]
stopifnot(nrow(g) == 914L, nrow(specimen_table) == 453L, nrow(primary) == 108L, nrow(expanded) == 129L,
  !anyDuplicated(primary$PatientID), !anyDuplicated(expanded$PatientID),
  sum(g$inclusion_flag_primary_ATB_vs_OD) == 216L, sum(g$inclusion_flag_expanded_ATB_vs_nonATB) == 266L,
  sum(primary$harmonized_label_4class == "ATB") == 91L, sum(primary$harmonized_label_4class == "OD") == 17L,
  identical(lock_before, unname(tools::md5sum(lock_path))))
message("SAMPLE_FLOW_RECONCILED: GSE89403 primary216 records->108 specimens; expanded266->129; lock unchanged")
print(totals[totals$analysis_split == "locked_external", ])
