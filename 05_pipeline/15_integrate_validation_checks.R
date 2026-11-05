# Supporting-analysis tables for the Tuberculosis submission.
# This stage does not fit models. It documents development-only selection scope,
# reference-signature implementation, model external artifacts, and sample flow.

if (!exists("project_path", mode = "function")) {
  Sys.setenv(TBML_PROJECT_ROOT = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  source(file.path("04_R_functions", "00_load_all.R"))
}

out_dir <- "results/ml_final"
bench_dir <- "results/benchmark"
man_table_dir <- "manuscript/tables"
dir.create(project_path(out_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(bench_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(project_path(man_table_dir), recursive = TRUE, showWarnings = FALSE)

config <- read_config("config/config.yml")
lock <- read_external_lockfile("config/external_validation_lockfile.json")
development_families <- unlist(lock$development_families)
external_families <- unlist(lock$external_validation_families)
coverage_threshold <- config$reference_signature_benchmark$minimum_gene_coverage_primary %||% 0.60

write_both <- function(x, result_path, manuscript_name, display = FALSE) {
  write_tsv_atomic(x, result_path)
  if (isTRUE(display)) {
    mlfinal_write_display_tsv(x, file.path(man_table_dir, manuscript_name))
  } else {
    write_tsv_atomic(x, file.path(man_table_dir, manuscript_name))
  }
}

# 1. Development-only reference-signature ranking table.
if (file.exists(project_path(file.path(bench_dir, "signature_per_family_performance.tsv")))) {
  family_perf <- read_tsv(file.path(bench_dir, "signature_per_family_performance.tsv"))
  overlap <- read_tsv(file.path(bench_dir, "signature_training_overlap_by_study.tsv"))
  overlap_keys <- unique(overlap[, c("signature_id", "cohort_family")])
  family_perf$overlap_key <- paste(family_perf$signature_id, family_perf$cohort_family, sep = "__")
  overlap_key <- paste(overlap_keys$signature_id, overlap_keys$cohort_family, sep = "__")
  include <- family_perf$signature_level %in% c("Level 0", "Level 1") &
    family_perf$coverage >= coverage_threshold &
    !family_perf$overlap_key %in% overlap_key &
    family_perf$cohort_family %in% development_families
  dev_primary <- family_perf[include, , drop = FALSE]
  if (nrow(dev_primary) > 0) {
    rank_dev <- aggregate(cbind(AUROC, AUPRC, balanced_accuracy) ~ task + signature_id, dev_primary, median, na.rm = TRUE)
    rank_dev$n_families <- aggregate(AUROC ~ task + signature_id, dev_primary, function(x) sum(!is.na(x)))$AUROC
    rank_dev$selection_scope <- "development_families_only"
    rank_dev <- rank_dev[order(rank_dev$task, -rank_dev$AUROC, -rank_dev$AUPRC), , drop = FALSE]
    write_tsv_atomic(rank_dev, file.path(bench_dir, "benchmark_rank_table_development.tsv"))
    write_tsv_atomic(rank_dev, file.path(bench_dir, "signature_summary_development_no_discovery_overlap.tsv"))
  }
}

# 2. Reference signature implementation and coverage table.
registry <- read_tsv(file.path(bench_dir, "signature_registry_resolved.tsv"))
coverage <- read_tsv(file.path(bench_dir, "signature_gene_coverage.tsv"))
orientation <- read_tsv(file.path(bench_dir, "signature_score_orientation.tsv"))
impl <- merge(coverage, registry, by = "signature_id", all.x = TRUE, sort = FALSE)
impl <- merge(impl, orientation, by = c("task", "signature_id"), all.x = TRUE, sort = FALSE)
direction_available <- impl$has_up_down_direction %in% TRUE | tolower(as.character(impl$has_up_down_direction)) == "true"
coefficients_available <- impl$has_coefficients %in% TRUE | tolower(as.character(impl$has_coefficients)) == "true"
training_mapped <- if ("training_studies_mapped" %in% names(impl)) as.character(impl$training_studies_mapped) else rep("", nrow(impl))
impl_out <- data.frame(
  task = impl$task,
  signature_id = impl$signature_id,
  source_package = impl$source_package,
  source_object = impl$source_object,
  registry_status = impl$registry_status,
  intended_comparison = impl$intended_task,
  original_gene_n = impl$n_signature_genes,
  genes_mapped_n = impl$n_covered,
  gene_coverage = impl$coverage,
  direction_available = direction_available,
  coefficients_available = coefficients_available,
  coefficient_source = ifelse(coefficients_available, "source_registry", "unavailable_in_current_registry"),
  scoring_method = impl$scoring_method,
  orientation = impl$orientation,
  orientation_source = impl$orientation_source,
  coverage_threshold = coverage_threshold,
  excluded_due_to_low_coverage = impl$coverage < coverage_threshold,
  signature_level = impl$signature_level.x %||% impl$signature_level.y,
  probe_to_gene_rule = "curatedTBData gene symbols; duplicate genes collapsed within study by highest mean expression",
  missing_gene_rule = "score computed from covered genes; low-coverage signatures flagged rather than silently imputed",
  overlap_with_training_families = ifelse(nzchar(training_mapped), "known_or_mapped", "unknown_or_none_recorded"),
  overlap_handling = "discovery/training-overlap-aware benchmark plus overlap-free locked external table",
  notes = "Current implementation uses TBSignatureProfiler gene sets when available; signed gene sets are scored as mean rank(up)-mean rank(down), unsigned gene sets as exploratory mean rank.",
  stringsAsFactors = FALSE
)
write_both(
  impl_out,
  file.path(out_dir, "reference_signature_implementation_and_gene_coverage.tsv"),
  "Supplementary_Table_reference_signature_implementation_and_gene_coverage.tsv",
  display = TRUE
)

# 3. Signature coverage sensitivity.
sig_scores <- read_tsv(file.path(bench_dir, "signature_scores_by_sample.tsv.gz"))
ext_summary <- read_tsv(file.path(out_dir, "external_pooled_summary.tsv"))
final_primary <- ext_summary[ext_summary$task == "primary_ATB_vs_OD" & ext_summary$selected_role == "best_performance", , drop = FALSE][1, ]
final_pred <- read_tsv(file.path(out_dir, "external_predictions.tsv.gz"))
final_pred <- final_pred[final_pred$task == "primary_ATB_vs_OD" & final_pred$selected_role == "best_performance", , drop = FALSE]
impl_primary <- impl_out[impl_out$task == "primary_ATB_vs_OD", , drop = FALSE]
filter_sets <- list(
  all_scored_signatures = unique(impl_primary$signature_id),
  coverage_ge_0_50 = unique(impl_primary$signature_id[impl_primary$gene_coverage >= 0.50]),
  coverage_ge_0_70 = unique(impl_primary$signature_id[impl_primary$gene_coverage >= 0.70]),
  direction_available = unique(impl_primary$signature_id[impl_primary$direction_available %in% TRUE]),
  coefficient_available = unique(impl_primary$signature_id[impl_primary$coefficients_available %in% TRUE])
)
sensitivity_rows <- lapply(names(filter_sets), function(filter_name) {
  ids <- filter_sets[[filter_name]]
  task_scores <- sig_scores[sig_scores$task == "primary_ATB_vs_OD" &
    sig_scores$cohort_family %in% external_families &
    sig_scores$signature_id %in% ids, , drop = FALSE]
  if (nrow(task_scores) == 0 || length(ids) == 0) {
    return(data.frame(
      coverage_filter = filter_name,
      signature_n = length(ids),
      best_signature = NA_character_,
      best_auc = NA_real_,
      selected_ml_auc = final_primary$AUROC,
      paired_delta_ml_vs_best_signature = NA_real_,
      n_shared_samples = NA_integer_,
      n_shared_families = NA_integer_,
      interpretation = "no signatures available under this filter",
      stringsAsFactors = FALSE
    ))
  }
  perf <- do.call(rbind, lapply(split(task_scores, task_scores$signature_id), function(d) {
    m <- compute_binary_metrics(d$true_label, d$score)
    data.frame(signature_id = d$signature_id[1], AUROC = m$AUROC, AUPRC = m$AUPRC, n = m$n, stringsAsFactors = FALSE)
  }))
  perf <- perf[order(-perf$AUROC, -perf$AUPRC), , drop = FALSE]
  best <- perf[1, , drop = FALSE]
  comp <- task_scores[task_scores$signature_id == best$signature_id[1], , drop = FALSE]
  shared <- merge(
    final_pred[, c("standard_sample_id", "cohort_family", "true_label", "score"), drop = FALSE],
    comp[, c("standard_sample_id", "score"), drop = FALSE],
    by = "standard_sample_id",
    suffixes = c("_ml", "_signature")
  )
  delta <- if (nrow(shared) > 0) compute_roc_auc(shared$true_label, shared$score_ml) - compute_roc_auc(shared$true_label, shared$score_signature) else NA_real_
  data.frame(
    coverage_filter = filter_name,
    signature_n = length(unique(ids)),
    best_signature = best$signature_id[1],
    best_auc = best$AUROC[1],
    selected_ml_auc = final_primary$AUROC,
    paired_delta_ml_vs_best_signature = delta,
    n_shared_samples = nrow(shared),
    n_shared_families = length(unique(shared$cohort_family)),
    interpretation = ifelse(is.finite(delta) && delta > 0, "selected ML higher than best signature under this filter", "no stable ML advantage indicated by this filter"),
    stringsAsFactors = FALSE
  )
})
coverage_sensitivity <- do.call(rbind, sensitivity_rows)
write_both(
  coverage_sensitivity,
  file.path(out_dir, "signature_coverage_sensitivity.tsv"),
  "Supplementary_Table_signature_coverage_sensitivity.tsv",
  display = TRUE
)

# 4. Development grid versus locked-external artifact status.
grid <- read_tsv(file.path(out_dir, "algorithm_combination_grid.tsv"))
dev_perf <- read_tsv(file.path(out_dir, "development_combination_performance.tsv"))
selected <- read_tsv(file.path(out_dir, "selected_combinations.tsv"))
external_check <- if (file.exists(project_path(file.path(out_dir, "all_model_development_external_evaluation.tsv")))) {
  read_tsv(file.path(out_dir, "all_model_development_external_evaluation.tsv"))
} else {
  data.frame(combination_id = character())
}
ext_all <- read_tsv(file.path(out_dir, "external_predictions.tsv.gz"))
ext_summary_for_artifacts <- read_tsv(file.path(out_dir, "external_pooled_summary.tsv"))
artifact <- merge(grid, dev_perf[, c("combination_id", "AUROC", "AUPRC", "median_family_AUROC")], by = "combination_id", all.x = TRUE, suffixes = c("", "_development"))
artifact <- merge(
  artifact,
  external_check[, intersect(c("combination_id", "external_AUROC", "external_AUPRC", "external_Brier", "external_calibration_slope", "rank_by_development", "rank_by_external", "evaluation_scope_note", "audit_scope"), names(external_check)), drop = FALSE],
  by = "combination_id",
  all.x = TRUE
)
artifact$selected_for_external <- artifact$combination_id %in% selected$combination_id
artifact$external_prediction_available <- artifact$combination_id %in% unique(c(external_check$combination_id, ext_all$combination_id))
artifact <- merge(
  artifact,
  ext_summary_for_artifacts[, intersect(c("combination_id", "AUROC", "AUPRC", "Brier", "calibration_slope"), names(ext_summary_for_artifacts)), drop = FALSE],
  by = "combination_id",
  all.x = TRUE,
  suffixes = c("", "_selected_summary")
)
for (pair in list(
  c("external_AUROC", "AUROC_selected_summary"),
  c("external_AUPRC", "AUPRC_selected_summary"),
  c("external_Brier", "Brier_selected_summary"),
  c("external_calibration_slope", "calibration_slope_selected_summary")
)) {
  if (all(pair %in% names(artifact))) {
    miss <- is.na(artifact[[pair[1]]])
    artifact[[pair[1]]][miss] <- artifact[[pair[2]]][miss]
  }
}
task_sample_n <- aggregate(standard_sample_id ~ task, ext_all, function(x) length(unique(x)))
names(task_sample_n)[2] <- "external_sample_n_selected_models"
artifact <- merge(artifact, task_sample_n, by = "task", all.x = TRUE)
artifact$external_family_n <- ifelse(artifact$external_prediction_available, length(external_families), NA_integer_)
artifact$external_sample_n <- ifelse(artifact$external_prediction_available, artifact$external_sample_n_selected_models, NA_integer_)
artifact$reason_no_external <- ifelse(
  artifact$external_prediction_available,
  "external_artifact_available",
  ifelse(artifact$selected_for_external, "selected_but_external_artifact_missing", "not_selected_for_locked_external_prediction")
)
artifact_out <- artifact[, intersect(c(
  "task", "combination_id", "normalization", "selector", "classifier", "selector_classifier",
  "AUROC", "AUPRC", "median_family_AUROC", "selected_for_external", "external_prediction_available",
  "reason_no_external", "external_AUROC", "external_AUPRC", "external_Brier", "external_calibration_slope",
  "rank_by_development", "rank_by_external", "external_family_n", "external_sample_n", "evaluation_scope_note", "audit_scope"
), names(artifact)), drop = FALSE]
names(artifact_out) <- sub("^AUROC$", "development_oof_auc", names(artifact_out))
names(artifact_out) <- sub("^AUPRC$", "development_oof_auprc", names(artifact_out))
write_both(
  artifact_out,
  file.path(out_dir, "model_external_artifact_status.tsv"),
  "Supplementary_Table_model_external_artifact_status.tsv",
  display = TRUE
)

# 5. Long sample-flow table with final-evaluable reasons.
flow <- read_tsv(file.path(man_table_dir, "Supplementary_Table_sample_flow_by_family.tsv"))
make_flow_long <- function(task_prefix, case_label, control_label, final_role) {
  task_name <- if (identical(task_prefix, "primary")) "primary_ATB_vs_OD" else "expanded_ATB_vs_nonATB"
  ready_col <- paste0(task_prefix, "_task_ready_n")
  final_col <- paste0(task_prefix, "_final_evaluated_n")
  excl_col <- paste0(task_prefix, "_exclusion_n")
  reason_col <- paste0(task_prefix, "_exclusion_reason")
  out <- data.frame(
    cohort_family = flow$cohort_family,
    task = task_name,
    analysis_split = ifelse(is.na(flow$split), "inventory_only", ifelse(flow$split == "external_validation", "locked_external", flow$split)),
    task_ready_n = flow[[ready_col]],
    final_evaluable_n = flow[[final_col]],
    excluded_from_final_n = flow[[excl_col]],
    exclusion_reason = flow[[reason_col]],
    reason_not_task_eligible = ifelse(flow[[ready_col]] > 0, "", "no task-ready samples after phenotype harmonization for this task"),
    final_set_note = ifelse(flow$split == "external_validation", "locked external final set uses selected-model/reference-signature evaluable samples", "not applicable to locked external final evaluation"),
    stringsAsFactors = FALSE
  )
  err <- read_tsv(file.path(out_dir, "external_error_by_sample.tsv"))
  err <- err[err$task == task_name & err$selected_role == final_role, , drop = FALSE]
  case_n <- aggregate(standard_sample_id ~ cohort_family, err[err$true_label == "ATB", , drop = FALSE], function(x) length(unique(x)))
  names(case_n)[2] <- "case_n_final"
  control_n <- aggregate(standard_sample_id ~ cohort_family, err[err$true_label != "ATB", , drop = FALSE], function(x) length(unique(x)))
  names(control_n)[2] <- "control_n_final"
  out <- merge(out, case_n, by = "cohort_family", all.x = TRUE)
  out <- merge(out, control_n, by = "cohort_family", all.x = TRUE)
  out$case_label <- case_label
  out$control_label <- control_label
  out$final_count_role <- ifelse(out$analysis_split == "locked_external", final_role, NA_character_)
  out
}
flow_long <- rbind(
  make_flow_long("primary", "ATB", "OD", "best_performance"),
  make_flow_long("expanded", "ATB", "nonATB", "best_sparse")
)
flow_long$task_eligible <- flow_long$task_ready_n > 0
write_both(
  flow_long,
  file.path(out_dir, "sample_flow_final_evaluable_reason.tsv"),
  "Supplementary_Table_sample_flow_final_evaluable_reason.tsv",
  display = TRUE
)

# 6. Claim audit ledger. This table records the high-risk claim classes handled by this reinforcement pass.
claim_audit <- data.frame(
  file = c(
    "09_submission/TB_ML_main_manuscript.md",
    "09_submission/TB_ML_main_manuscript.md",
    "09_submission/TB_ML_main_manuscript.md",
    "09_submission/TB_ML_supplementary_information.md",
    "04_R_functions/22_figure_02_locked_external_benchmark.R",
    "09_submission/build_tuberculosis_submission_package.py"
  ),
  risk_type = c(
    "title_and_abstract_overclaim",
    "reference_signature_implementation_mismatch",
    "952_external_scope_ambiguity",
    "within_family_permutation_additive_decomposition",
    "ambiguous_family_oracle_anchor",
    "cover_letter_and_docx_package_overclaim"
  ),
  revised_sentence = c(
    "Use limited incremental gain / family-locked benchmark wording rather than structural-bound overclaiming.",
    "Describe TBSignatureProfiler gene-set scoring with signed or unsigned rank means and development-only orientation; do not claim original coefficients/formulas unless available.",
    "Bind 952 to development-side screening; locked external evaluation is for selected/top models with saved artifacts.",
    "Interpret within-family permutation as a negative-control stress test consistent with residual family-level score positioning, not an additive AUROC decomposition.",
    "Label the 0.937 anchor as a development family-label oracle and state it is not an external classifier benchmark.",
    "Use the same conservative limited-incremental-gain framing in generated DOCX title, abstract, highlights, and cover letter."
  ),
  status = "implemented_or_source_patched",
  stringsAsFactors = FALSE
)
write_tsv_atomic(claim_audit, file.path(out_dir, "manuscript_wording_checks.tsv"))

# 7. Compact integrated strengthening summary. This replaces a stand-alone
# supplementary figure for the same evidence layer.
fmt_num <- function(x, digits = 3) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return("NA")
  formatC(x[1], format = "f", digits = digits)
}
fmt_range <- function(x, digits = 3) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return("NA")
  r <- range(x, na.rm = TRUE)
  if (abs(r[1] - r[2]) < 10^-digits) return(formatC(r[1], format = "f", digits = digits))
  paste0(formatC(r[1], format = "f", digits = digits), "-", formatC(r[2], format = "f", digits = digits))
}
fmt_pct <- function(x, digits = 1) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return("NA")
  paste0(formatC(100 * x[1], format = "f", digits = digits), "%")
}

family_def_summary <- read_tsv("manuscript/tables/Supplementary_Table_family_definition_stability_matrix.tsv")
split_lite_summary <- read_tsv("manuscript/tables/Supplementary_Table_family_block_split_lite_summary.tsv")
family_calibration_summary <- read_tsv("manuscript/tables/Supplementary_Table_family_specific_calibration_thresholds.tsv")
no_leakage_summary <- read_tsv("manuscript/tables/Supplementary_Table_data_separation_checks.tsv")
negative_control_summary <- read_tsv("manuscript/tables/Supplementary_Table_negative_control_summary.tsv")
within_family_perm_summary <- read_tsv("results/ml_final/negative_control_family_stratified_permutation.tsv")
claim_map_summary <- read_tsv(file.path(out_dir, "statement_support_summary.tsv"))
module_coverage_summary <- read_tsv(file.path(out_dir, "supporting_analysis_coverage.tsv"))

fd_primary <- family_def_summary[
  family_def_summary$task == "primary_ATB_vs_OD" &
    family_def_summary$selected_role == "best_performance",
  ,
  drop = FALSE
]
split_primary <- split_lite_summary[
  split_lite_summary$task == "primary_ATB_vs_OD" &
    split_lite_summary$selected_role == "best_performance",
  ,
  drop = FALSE
]
split_expanded <- split_lite_summary[
  split_lite_summary$task == "expanded_ATB_vs_nonATB" &
    split_lite_summary$selected_role == "best_performance",
  ,
  drop = FALSE
]
cal_primary <- family_calibration_summary[
  family_calibration_summary$task == "primary_ATB_vs_OD" &
    family_calibration_summary$selected_role == "best_performance",
  ,
  drop = FALSE
]
collapse_primary <- cal_primary[cal_primary$all_positive_like %in% TRUE | cal_primary$severe_positive_bias %in% TRUE, , drop = FALSE]
collapse_text <- if (nrow(collapse_primary) > 0) {
  paste(
    paste0(collapse_primary$cohort_family, " predicted-positive ", vapply(collapse_primary$predicted_positive_fraction, fmt_pct, character(1))),
    collapse = "; "
  )
} else {
  "No all-positive-like external family flagged"
}

wf_primary <- within_family_perm_summary[
  within_family_perm_summary$task == "primary_ATB_vs_OD" &
    within_family_perm_summary$selected_role == "best_performance",
  ,
  drop = FALSE
]
global_iterations <- unique(negative_control_summary$n_iterations)
global_iteration_text <- if (length(global_iterations) == 1) as.character(global_iterations) else paste(global_iterations, collapse = "/")

integrated_strengthening_summary <- data.frame(
  strengthening_layer = c(
    "Negative-control resampling",
    "Family-definition sensitivity",
    "Family-block split Lite",
    "Family-specific threshold collapse",
    "No-leakage and artifact-scope audit",
    "Claim and module integration"
  ),
  module_link = c(
    "Module 01; Module 07",
    "Module 02; Module 04",
    "Module 03 Lite",
    "Module 08",
    "Module 10",
    "Module 12"
  ),
  key_result = c(
    paste0(
      "Global negative controls use ", global_iteration_text, " iterations; primary within-family score permutation median AUROC ",
      fmt_num(wf_primary$AUROC_median), " (", fmt_num(wf_primary$AUROC_q025), "-", fmt_num(wf_primary$AUROC_q975), ")."
    ),
    paste0(
      "Across four encoded family definitions, primary mixed development families range ",
      fmt_range(fd_primary$development_n_mixed_families, digits = 0), "; development family-label oracle AUROC ",
      fmt_range(fd_primary$development_family_oracle_AUROC), "; within-family permutation median AUROC ",
      fmt_range(fd_primary$within_family_permutation_AUROC_median), "."
    ),
    paste0(
      "Primary best model: ", split_primary$n_splits[1], " family-block rotations, median ML-minus-best-reference AUROC ",
      fmt_num(split_primary$delta_ML_vs_best_ref_median), ", ML beats best reference in ",
      fmt_pct(split_primary$prop_ML_beats_best_ref), " of rotations; expanded best model: median delta ",
      fmt_num(split_expanded$delta_ML_vs_best_ref_median), ", beats in ", fmt_pct(split_expanded$prop_ML_beats_best_ref), "."
    ),
    paste0(collapse_text, "."),
    paste0(
      sum(tolower(no_leakage_summary$status) == "pass"), "/", nrow(no_leakage_summary),
      " no-leakage checks pass, covering lockfile stability, artifact scope, selected predictions, threshold locks, and grouped-split tests."
    ),
    paste0(
      nrow(claim_map_summary), " manuscript claims mapped to supporting evidence; ",
      paste(module_coverage_summary$module, collapse = ", "), " integrated, with Module 11 intentionally excluded."
    )
  ),
  interpretation_for_submission = c(
    "Permutation and random-score checks support the view that observed performance is not a trivial label/random-score artifact, while preserving family-level caution.",
    "Family-boundary choices do not remove the family-structure signal; the result should be described as sensitivity support, not a full refit analysis.",
    "The primary ML advantage is not stable under family-block rotations; expanded-task behavior is less fragile but still context dependent.",
    "Primary external operating behavior is family-specific, with some families pushed toward near all-positive calls at the locked threshold.",
    "The external validation story remains development-locked and auditable; external artifacts are scoped rather than used for model selection.",
    "The final manuscript claims are explicitly softened and traceable to tables rather than unsupported visual rhetoric."
  ),
  scope_limit = c(
    "Within-family permutation is a structural stress test, not an additive AUROC decomposition.",
    "Family definitions are evaluated by post-hoc regrouping rather than full model refitting.",
    "Lite rotations reuse fixed candidates and development out-of-fold predictions.",
    "Family-specific calibration can be unstable or non-estimable in small/separated families.",
    "Machine-readable checks do not replace independent code review.",
    "The integration table does not add biological mechanism claims."
  ),
  source_tables = c(
    "Supplementary_Table_negative_control_summary.tsv; Supplementary_Table_negative_control_family_stratified_permutation.tsv",
    "Supplementary_Table_family_definition_stability_matrix.tsv; Supplementary_Table_family_phenotype_structure.tsv",
    "Supplementary_Table_family_block_split_lite_summary.tsv",
    "Supplementary_Table_family_specific_calibration_thresholds.tsv; Supplementary_Table_calibrated_operating_points_external.tsv",
    "Supplementary_Table_data_separation_checks.tsv; Supplementary_Table_model_external_artifact_status.tsv",
    "results/ml_final/statement_support_summary.tsv; results/ml_final/supporting_analysis_coverage.tsv"
  ),
  stringsAsFactors = FALSE
)
write_tsv_atomic(integrated_strengthening_summary, file.path(out_dir, "supporting_analysis_summary.tsv"))

# 8. Update manuscript table index with the new reinforcement tables.
index_path <- file.path(man_table_dir, "Table_index.tsv")
if (file.exists(project_path(index_path))) {
  index <- read_tsv(index_path)
  additions <- data.frame(
    table_file = c(
      "Supplementary_Table_reference_signature_implementation_and_gene_coverage.tsv",
      "Supplementary_Table_signature_coverage_sensitivity.tsv",
      "Supplementary_Table_model_external_artifact_status.tsv",
      "Supplementary_Table_sample_flow_final_evaluable_reason.tsv"
    ),
    placement = c("supplementary", "supplementary", "source_data", "supplementary"),
    title = c(
      "Reference-signature implementation and gene coverage",
      "Reference-signature coverage sensitivity",
      "External-evaluation scope",
      "Sample flow and final-evaluable reasons"
    ),
    purpose = c(
      "Documents gene coverage, direction/coefficient availability, scoring method, and overlap handling for each reference signature.",
      "Checks whether primary-task conclusions change under simple signature-coverage filters.",
      "Separates the 952-combination development grid from the subset with saved locked-external prediction outputs.",
      "Explains task-ready versus final-evaluable sample counts by family and task, including GSE89403 attrition."
    ),
    source_table = c(
      "results/benchmark/signature_gene_coverage.tsv; results/benchmark/signature_score_orientation.tsv",
      "results/benchmark/signature_scores_by_sample.tsv.gz; results/ml_final/external_predictions.tsv.gz",
      "results/ml_final/algorithm_combination_grid.tsv; results/ml_final/all_model_development_external_evaluation.tsv",
      "manuscript/tables/Supplementary_Table_sample_flow_by_family.tsv; results/ml_final/external_error_by_sample.tsv"
    ),
    notes = c(
      "Included in Supplementary Table S3.",
      "Included in Supplementary Table S3.",
      "Retained under source_data/table_sources and summarized in Supplementary Table S4.",
      "Included in Supplementary Table S1."
    ),
    exists_in_package = TRUE,
    stringsAsFactors = FALSE
  )
  index <- index[!index$table_file %in% additions$table_file, , drop = FALSE]
  index <- rbind(index, additions)
  write_tsv_atomic(index, index_path)
}

message("Computational reinforcement tables completed.")
