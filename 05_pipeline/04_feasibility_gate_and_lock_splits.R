config <- read_config("config/config.yml")
if (identical(Sys.getenv("TBML_REFRESH_EXTERNAL_LOCKFILE"), "1")) {
  config$project$overwrite_external_validation_lockfile <- TRUE
}
master <- read_tsv("data_meta/sample_harmonization_master.tsv")

task_defs <- config$labels$task_definitions
task_rows <- list()
family_rows <- list()
for (task in names(task_defs)) {
  flag <- switch(task,
    primary_ATB_vs_OD = "inclusion_flag_primary_ATB_vs_OD",
    secondary_ATB_vs_LTBI = "inclusion_flag_secondary_ATB_vs_LTBI",
    expanded_ATB_vs_nonATB = "inclusion_flag_expanded_ATB_vs_nonATB"
  )
  label_col <- switch(task,
    primary_ATB_vs_OD = "harmonized_label_binary_primary_ATB_vs_OD",
    secondary_ATB_vs_LTBI = "harmonized_label_binary_secondary_ATB_vs_LTBI",
    expanded_ATB_vs_nonATB = "harmonized_label_binary_expanded_ATB_vs_nonATB"
  )
  tm <- master[master[[flag]] == TRUE & !is.na(master[[label_col]]) & nzchar(master[[label_col]]), , drop = FALSE]
  tm$task_label <- tm[[label_col]]
  tab <- as.data.frame.matrix(table(tm$cohort_family, tm$task_label))
  if (nrow(tab) == 0) tab <- data.frame()
  tab$cohort_family <- rownames(tab)
  pos_name <- task_defs[[task]]$positive[[1]]
  neg_names <- if (task == "expanded_ATB_vs_nonATB") "nonATB" else unlist(task_defs[[task]]$negative)
  pos_count <- if (pos_name %in% names(tab)) tab[[pos_name]] else rep(0, nrow(tab))
  neg_count <- rowSums(tab[, intersect(neg_names, names(tab)), drop = FALSE])
  eligible_external <- sum(pos_count >= config$external_validation$min_per_class_external & neg_count >= config$external_validation$min_per_class_external)
  independent_families <- length(unique(tm$cohort_family))
  total_samples <- nrow(tm)
  decision <- "red"
  if (task == config$feasibility_gate$primary_task &&
      independent_families >= config$feasibility_gate$green_min_independent_families &&
      eligible_external >= config$feasibility_gate$green_min_external_families &&
      total_samples >= config$feasibility_gate$green_min_total_samples) {
    decision <- "green"
  } else if (independent_families >= config$feasibility_gate$yellow_min_independent_families &&
             eligible_external >= config$feasibility_gate$yellow_min_external_families &&
             total_samples >= config$feasibility_gate$yellow_min_total_samples) {
    decision <- "yellow"
  }
  task_rows[[task]] <- data.frame(
    task = task,
    n_samples = total_samples,
    n_independent_families = independent_families,
    n_external_eligible_families = eligible_external,
    decision = decision,
    stringsAsFactors = FALSE
  )
  if (nrow(tab) > 0) {
    count_cols <- setdiff(names(tab), "cohort_family")
    fam_long <- do.call(rbind, lapply(count_cols, function(label) {
      data.frame(
        task = task,
        cohort_family = tab$cohort_family,
        label = label,
        n = tab[[label]],
        stringsAsFactors = FALSE
      )
    }))
    family_rows[[task]] <- fam_long
  }
}
feas_by_task <- do.call(rbind, task_rows)
write_tsv_atomic(feas_by_task, "results/tables/feasibility_by_task.tsv")
write_tsv_atomic(do.call(rbind, family_rows), "results/tables/feasibility_by_study_family.tsv")

primary_ATB_vs_OD <- master[master$inclusion_flag_primary_ATB_vs_OD == TRUE & !is.na(master$harmonized_label_binary_primary_ATB_vs_OD), , drop = FALSE]
primary_ATB_vs_OD$task_label <- primary_ATB_vs_OD$harmonized_label_binary_primary_ATB_vs_OD
decision <- feas_by_task[feas_by_task$task == config$feasibility_gate$primary_task, , drop = FALSE]
write_tsv_atomic(decision, "results/tables/feasibility_decision.tsv")

if (nrow(primary_ATB_vs_OD) > 0 && decision$decision %in% c("green", "yellow")) {
  split <- make_external_validation_split(primary_ATB_vs_OD, "cohort_family", config)
  lock <- write_external_lockfile(split, config, hash_file_or_object(primary_ATB_vs_OD))
  split <- list(
    task = lock$task,
    split_unit = lock$split_unit,
    development_families = unlist(lock$development_families),
    external_validation_families = unlist(lock$external_validation_families),
    external_selection_algorithm = lock$external_selection_algorithm
  )
  assert_no_overlap_between_dev_and_external(split)
  plan <- data.frame(
    cohort_family = c(split$development_families, split$external_validation_families),
    split = c(rep("development", length(split$development_families)), rep("external_validation", length(split$external_validation_families))),
    stringsAsFactors = FALSE
  )
  write_tsv_atomic(plan, "results/tables/external_validation_plan.tsv")
  bal <- as.data.frame.matrix(table(primary_ATB_vs_OD$cohort_family, primary_ATB_vs_OD$task_label))
  bal$cohort_family <- rownames(bal)
  bal <- merge(plan, bal, by = "cohort_family", all.x = TRUE)
  write_tsv_atomic(bal, "results/tables/development_validation_split_balance.tsv")
} else {
  write_tsv_atomic(data.frame(cohort_family = character(), split = character()), "results/tables/external_validation_plan.tsv")
  write_tsv_atomic(data.frame(), "results/tables/development_validation_split_balance.tsv")
}

write_tsv_atomic(data.frame(
  item = "primary_ATB_vs_OD_decision",
  value = decision$decision,
  stringsAsFactors = FALSE
), "report/feasibility_summary.tsv")
