ds <- read_tsv("results/tables/DataSummary_full.tsv")
master <- read_tsv("data_meta/sample_harmonization_master.tsv")

family_for_study <- function(study) {
  if (study %in% c("GSE19435", "GSE19439", "GSE19442", "GSE19443", "GSE19444")) {
    return(c(family = "FAM_GSE19435_19444_Berry_related", rule = "known_same_publication_series", confidence = "high"))
  }
  if (study %in% c("GSE42825", "GSE42826", "GSE42827", "GSE42830", "GSE42831", "GSE42832")) {
    return(c(family = "FAM_GSE42825_42832_related", rule = "known_same_publication_series", confidence = "high"))
  }
  if (study %in% c("GSE107991", "GSE107992", "GSE107993", "GSE107994")) {
    return(c(family = "FAM_GSE107991_107994_related", rule = "adjacent_accession_series", confidence = "medium"))
  }
  c(family = study, rule = "study_as_family_default", confidence = "medium")
}

rows <- lapply(unique(ds$GEOAccession), function(study) {
  fam <- family_for_study(study)
  d <- ds[match(study, ds$GEOAccession), , drop = FALSE]
  data.frame(
    study = study,
    cohort_family = unname(fam[["family"]]),
    parent_accession = study,
    publication_key = unname(fam[["family"]]),
    platform = d$Platform,
    region = d[["Country/Region"]],
    family_assignment_rule = unname(fam[["rule"]]),
    manual_override_flag = FALSE,
    confidence = unname(fam[["confidence"]]),
    notes = ifelse(unname(fam[["rule"]]) == "study_as_family_default", "No same-family rule detected automatically.", "Grouped to reduce hidden cohort leakage."),
    stringsAsFactors = FALSE
  )
})
family_map <- do.call(rbind, rows)
write_tsv_atomic(family_map, "data_meta/cohort_family_map.tsv")
write_tsv_atomic(family_map[family_map$confidence != "high", , drop = FALSE], "results/tables/uncertain_family_assignment.tsv")

master$cohort_family <- family_map$cohort_family[match(master$study, family_map$study)]
write_tsv_atomic(master, "data_meta/sample_harmonization_master.tsv")

