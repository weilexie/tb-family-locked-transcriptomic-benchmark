config <- read_config("config/config.yml")
if (!requireNamespace("curatedTBData", quietly = TRUE)) {
  stop("curatedTBData is required for 01_pull_curatedTBData_metadata.R. Install it with BiocManager::install('curatedTBData').", call. = FALSE)
}

ds <- load_datasummary()
sig <- load_signature_info_training()
write_tsv_atomic(ds, "results/tables/DataSummary_full.tsv")
write_tsv_atomic(sig, "results/tables/SignatureInfoTraining.tsv")

raw_cols <- data.frame(
  table = c(rep("DataSummary", ncol(ds)), rep("SignatureInfoTraining", ncol(sig))),
  column = c(names(ds), names(sig)),
  stringsAsFactors = FALSE
)
write_tsv_atomic(raw_cols, "results/tables/study_metadata_raw_columns.tsv")

studies <- unique(ds$GEOAccession)
dry <- curatedtb_dry_run(studies, curated_only = config$primary_data_source$use_curated_only_first %||% TRUE)
write_tsv_atomic(dry, "results/tables/curatedTBData_resource_dryrun.tsv")

write_tsv_atomic(data.frame(
  item = c("run_date", "R_version", "Bioconductor_version", "curatedTBData_version", "TBSignatureProfiler_version", "DataSummary_rows", "DataSummary_unique_studies", "SignatureInfoTraining_available", "SignatureInfoTraining_rows", "Dry_run_resources_queried"),
  value = c(
    as.character(Sys.time()),
    as.character(getRversion()),
    if (requireNamespace("BiocManager", quietly = TRUE)) as.character(BiocManager::version()) else "not_available",
    as.character(utils::packageVersion("curatedTBData")),
    if (requireNamespace("TBSignatureProfiler", quietly = TRUE)) as.character(utils::packageVersion("TBSignatureProfiler")) else "not_installed",
    nrow(ds),
    length(unique(ds$GEOAccession)),
    nrow(sig) > 0,
    nrow(sig),
    nrow(dry)
  ),
  stringsAsFactors = FALSE
), "report/references_checked.tsv")
