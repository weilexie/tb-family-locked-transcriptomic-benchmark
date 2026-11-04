# Reviewer R1.3: published fixed formulas, separated from submitted gene-set scores.
# Run from the revision-copy root. No training, threshold fitting, or external selection.
options(stringsAsFactors = FALSE)
root <- normalizePath(Sys.getenv("TBML_PROJECT_ROOT", unset = getwd()), winslash = "/")
if (!grepl("TB ML_revision_temp_20261103_2308$", root)) {
  stop("This revision script is restricted to the verified revision copy.")
}
Sys.setenv(TBML_PROJECT_ROOT = root)
source(file.path(root, "04_reusable_R_functions", "00_load_all.R"))
out_dir <- file.path(root, "07_results", "04_revision", "reference_signatures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, name) {
  data.table::fwrite(x, file.path(out_dir, name), sep = "\t", na = "NA")
}
lock_path <- file.path(root, "01_config", "external_validation_lockfile.json")
lock_hash <- unname(tools::md5sum(lock_path))
lock <- jsonlite::fromJSON(lock_path)
tasks <- c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB")

# Fixed formula and direction taken from the publications, never selected on this data.
# Gliddon3 was developed for TB vs LTBI and remains a contextual comparator.
specs <- list(
  Sweeney_OD_3 = list(
    genes = c("GBP5", "DUSP3", "KLF2"),
    weights = c(0.5, 0.5, -1),
    formula = "(GBP5 + DUSP3)/2 - KLF2",
    role = "published_fixed_formula_primary_comparator",
    citation = "Sweeney et al. 2016; DOI:10.1016/S2213-2600(16)00048-5",
    source = paste0("Formula independently tabulated in Turner et al. 2020 Table 2, ",
                    "PMC7113842; log2 means verified in Warsinske et al. 2019 S1 Text, PMC6478271.")),
  Gliddon_OD_4 = list(
    genes = c("GBP6", "PRDM1", "TMCC1", "ARG1"),
    weights = c(1, 1, -1, -1),
    formula = "GBP6 + PRDM1 - TMCC1 - ARG1",
    role = "published_fixed_formula_primary_comparator",
    citation = "Gliddon et al. 2021; DOI:10.3389/fimmu.2021.637164",
    source = "Original Table 1 directions and DRS method, PMC7982854; explicit formula in Turner et al. 2020 Table 2."),
  Gliddon_OD_3 = list(
    genes = c("FCGR1A", "C1QB", "ZNF296"),
    weights = c(1, 1, -1),
    formula = "FCGR1A + C1QB - ZNF296",
    role = "published_fixed_formula_context_TB_vs_LTBI",
    citation = "Gliddon et al. 2021; DOI:10.3389/fimmu.2021.637164",
    source = "Original Table 1 directions and DRS method, PMC7982854; explicit formula in Turner et al. 2020 Table 2.")
)
needed_genes <- unique(c(unlist(lapply(specs, `[[`, "genes")),
                         "UB2L6", "UBE2L6", "FCGR1C", "FCGR1B"))
meta <- lapply(tasks, function(task) {
  x <- as.data.frame(data.table::fread(file.path(root, "03_data_curated", "02_task_matrices", task, "metadata.tsv")))
  stopifnot(!anyDuplicated(x$standard_sample_id))
  x$task <- task
  x$split <- ifelse(x$cohort_family %in% lock$external_validation_families,
                    "locked_external", "development")
  stopifnot(all(x$cohort_family %in% c(lock$development_families, lock$external_validation_families)))
  x
})
names(meta) <- tasks
all_meta <- do.call(rbind, meta)
study_rows <- list()
availability <- list()

# Read only selected genes from each study before task-universe intersection.
# Missing entire genes remain missing; no surrogate gene or cross-study imputation is used.
read_needed <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1, warn = FALSE)
  keep <- character()
  repeat {
    block <- readLines(con, n = 500, warn = FALSE)
    if (!length(block)) break
    first <- sub("\t.*$", "", block)
    first <- gsub('^"|"$', "", first)
    keep <- c(keep, block[first %in% needed_genes])
  }
  d <- as.data.frame(data.table::fread(text = paste(c(header, keep), collapse = "\n")))
  g <- d$gene
  x <- as.matrix(d[, setdiff(names(d), "gene"), drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- g
  stopifnot(!anyDuplicated(g))
  x
}
for (study in sort(unique(all_meta$study))) {
  message("Published fixed-formula inputs: ", study)
  study_dir <- file.path(root, "03_data_curated", "01_by_study", study)
  x <- read_needed(file.path(study_dir, "expression_gene_symbol.tsv.gz"))
  qc <- as.data.frame(data.table::fread(file.path(study_dir, "preprocessing_log.tsv")))
  availability[[length(availability) + 1]] <- data.frame(
    study = study, gene = needed_genes, available_before_task_intersection = needed_genes %in% rownames(x),
    assay_name = qc$assay_name[1], log2_decision = qc$log2_decision[1],
    source_expression = paste0("03_data_curated/01_by_study/", study, "/expression_gene_symbol.tsv.gz"))
  for (task in tasks) {
    m <- meta[[task]]
    m <- m[m$study == study, , drop = FALSE]
    if (!nrow(m)) next
    stopifnot(all(m$standard_sample_id %in% colnames(x)))
    for (sig in names(specs)) {
      sp <- specs[[sig]]
      missing <- setdiff(sp$genes, rownames(x))
      scores <- rep(NA_real_, nrow(m))
      if (!length(missing)) {
        v <- x[sp$genes, m$standard_sample_id, drop = FALSE]
        scores <- as.numeric(crossprod(sp$weights, v))
        scores[colSums(!is.finite(v)) > 0] <- NA_real_
      }
      study_rows[[length(study_rows) + 1]] <- data.frame(
        task = task, signature_id = sig, standard_sample_id = m$standard_sample_id,
        study = study, cohort_family = m$cohort_family, split = m$split,
        true_label = m$task_label, score = scores,
        implementation = "published_fixed_formula_on_harmonized_log_expression",
        formula = sp$formula, comparator_role = sp$role,
        n_required_genes = length(sp$genes), n_available_genes = length(sp$genes) - length(missing),
        missing_genes = paste(missing, collapse = ";"),
        input_scale = "study-level curated expression after the existing log2-detection preprocessing",
        status = ifelse(is.finite(scores), "complete_formula_evaluable", "missing_required_gene"))
    }
  }
}
scores <- do.call(rbind, study_rows)
write_out(scores, "published_fixed_formula_scores.tsv.gz")
write_out(do.call(rbind, availability), "gene_availability_before_task_intersection.tsv")
write_out(do.call(rbind, lapply(names(specs), function(sig) {
  s <- specs[[sig]]
  data.frame(signature_id = sig, formula = s$formula, required_genes = paste(s$genes, collapse = ";"),
             fixed_weights = paste(s$weights, collapse = ";"), role = s$role,
             citation = s$citation, verification_source = s$source,
             reproduction_boundary = "Published arithmetic and fixed gene direction reproduced; upstream platform preprocessing and original diagnostic assay are not reproduced.")
})), "published_formula_sources.tsv")

auc <- function(y, s) {
  y <- as_binary(y)
  k <- is.finite(s) & !is.na(y)
  y <- y[k]; s <- s[k]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(s, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
perf <- function(d) {
  d <- d[is.finite(d$score), , drop = FALSE]
  data.frame(n = nrow(d), n_ATB = sum(d$true_label == "ATB"), n_nonATB = sum(d$true_label != "ATB"),
             n_families = length(unique(d$cohort_family)), AUROC = auc(d$true_label, d$score),
             AUPRC = compute_pr_auc(d$true_label, d$score))
}
summary_rows <- lapply(split(scores, interaction(scores$task, scores$signature_id, scores$split, drop = TRUE)), function(d) {
  cbind(d[1, c("task", "signature_id", "split", "implementation", "comparator_role")], perf(d))
})
write_out(do.call(rbind, summary_rows), "published_fixed_formula_performance.tsv")
family_rows <- lapply(split(scores, interaction(scores$task, scores$signature_id, scores$cohort_family, drop = TRUE)), function(d) {
  cbind(d[1, c("task", "signature_id", "cohort_family", "split", "comparator_role")], perf(d))
})
write_out(do.call(rbind, family_rows), "published_fixed_formula_per_family.tsv")

old <- as.data.frame(data.table::fread(file.path(root, "07_results", "02_reference_signature_benchmark", "signature_scores_by_sample.tsv.gz")))
old <- old[old$signature_id %in% names(specs) & old$task %in% tasks, ]
paired_impl <- merge(scores, old[, c("task", "signature_id", "standard_sample_id", "score", "scoring_method")],
                     by = c("task", "signature_id", "standard_sample_id"), suffixes = c("_formula", "_submitted"))
impl_rows <- lapply(split(paired_impl, interaction(paired_impl$task, paired_impl$signature_id, paired_impl$split, drop = TRUE)), function(d) {
  d <- d[is.finite(d$score_formula) & is.finite(d$score_submitted), ]
  if (!nrow(d)) return(NULL)
  data.frame(task = d$task[1], signature_id = d$signature_id[1], split = d$split[1],
             n_shared = nrow(d), n_families = length(unique(d$cohort_family)),
             formula_AUROC = auc(d$true_label, d$score_formula),
             submitted_gene_set_AUROC = auc(d$true_label, d$score_submitted),
             formula_AUPRC = compute_pr_auc(d$true_label, d$score_formula),
             submitted_gene_set_AUPRC = compute_pr_auc(d$true_label, d$score_submitted),
             submitted_scoring_method = d$scoring_method[1],
             interpretation = "Descriptive comparison of implementations on identical evaluable samples")
})
write_out(do.call(rbind, impl_rows), "formula_vs_submitted_gene_set_same_samples.tsv")

external <- as.data.frame(data.table::fread(file.path(root, "07_results", "03_ml_final", "external_predictions.tsv.gz")))
set.seed(20261103)
paired_rows <- list()
for (task in tasks) for (sig in names(specs)) {
  sd <- scores[scores$task == task & scores$signature_id == sig & scores$split == "locked_external" & is.finite(scores$score), ]
  for (role in c("best_performance", "best_sparse")) {
    md <- external[external$task == task & external$selected_role == role, ]
    d <- merge(md[, c("standard_sample_id", "score", "true_label", "cohort_family")],
               sd[, c("standard_sample_id", "score")], by = "standard_sample_id", suffixes = c("_model", "_reference"))
    if (!nrow(d)) next
    stopifnot(!anyDuplicated(d$standard_sample_id))
    families <- sort(unique(d$cohort_family))
    boot <- replicate(1000, {
      f <- sample(families, length(families), replace = TRUE)
      idx <- unlist(lapply(f, function(x) which(d$cohort_family == x)), use.names = FALSE)
      c(auc(d$true_label[idx], d$score_model[idx]), auc(d$true_label[idx], d$score_reference[idx]))
    })
    deltas <- boot[1, ] - boot[2, ]
    qr <- quantile(boot[2, ], c(.025, .975), na.rm = TRUE, names = FALSE)
    qd <- quantile(deltas, c(.025, .975), na.rm = TRUE, names = FALSE)
    paired_rows[[length(paired_rows) + 1]] <- data.frame(
      task = task, selected_role = role, signature_id = sig, n_shared = nrow(d),
      n_ATB = sum(d$true_label == "ATB"), n_families = length(families),
      selected_model_AUROC = auc(d$true_label, d$score_model),
      published_formula_AUROC = auc(d$true_label, d$score_reference),
      formula_family_ci_low = qr[1], formula_family_ci_high = qr[2],
      delta_model_minus_formula = auc(d$true_label, d$score_model) - auc(d$true_label, d$score_reference),
      paired_family_ci_low = qd[1], paired_family_ci_high = qd[2],
      selected_model_AUPRC = compute_pr_auc(d$true_label, d$score_model),
      published_formula_AUPRC = compute_pr_auc(d$true_label, d$score_reference),
      bootstrap_repetitions = 1000, seed = 20261103,
      inference_scope = "Exploratory unadjusted family-cluster intervals; no stable rank, equivalence, or superiority inference")
  }
}
write_out(do.call(rbind, paired_rows), "same_lock_paired_model_vs_published_formula.tsv")

# Hoang3 has no fixed published LDA weights in the inspected article/appendix.
# The symbol correction is source-based: probe 2070170 is UBE2L6 in the same
# appendix (signature 2 and the DET table), although signature 3 prints UB2L6.
hoang <- data.frame(
  signature_id = "Hoang_OD_3", published_probe = c("1510364", "5570398", "2070170"),
  published_signature3_symbol = c("GBP5", "FCGR1C", "UB2L6"),
  source_supported_symbol = c("GBP5", "FCGR1C", "UBE2L6"),
  source = c("Hoang2021 appendix PDF p5", "Hoang2021 appendix PDF p5", "Hoang2021 appendix PDF p5 signature2 and p10 DET table"),
  algorithm_status = "LDA retrained over repeated source-cohort splits; fixed final coefficients absent from inspected publication and appendix",
  benchmark_status = "Not an exact frozen classifier; submitted one-gene mean-rank score is exploratory and must not headline original-Hoang superiority")
write_out(hoang, "Hoang3_original_source_verification.tsv")
stopifnot(identical(unname(tools::md5sum(lock_path)), lock_hash))
message("REVISION_ORIGINAL_FORMULAS_OK; external lock unchanged: ", lock_hash)
