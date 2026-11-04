args <- commandArgs(trailingOnly = TRUE)
mode <- "manuscript_final"
if ("--mode" %in% args) {
  idx <- match("--mode", args)
  if (length(args) < idx + 1) stop("--mode was provided without a value", call. = FALSE)
  mode <- args[[idx + 1]]
}

valid_modes <- c("prepare_data", "optimize_final", "enhance_final", "stress_test", "manuscript_final")
if (!mode %in% valid_modes) {
  stop(sprintf("Invalid mode '%s'. Expected one of: %s", mode, paste(valid_modes, collapse = ", ")), call. = FALSE)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
Sys.setenv(TBML_PROJECT_ROOT = project_root, TBML_RUN_MODE = mode)

tbml_directory_aliases <- c(
  config = "01_config",
  data_meta = "02_data_meta",
  data_curated = "03_data_curated",
  R = "04_reusable_R_functions",
  scripts = "05_pipeline_scripts",
  tests = "06_validation_tests",
  results = "07_results",
  manuscript = "08_manuscript",
  report = "10_report",
  logs = "11_logs"
)
tbml_resolve <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  first <- sub("/.*$", "", normalized)
  rest <- sub("^[^/]+/?", "", normalized)
  mapped_first <- tbml_directory_aliases[[first]]
  mapped <- if (!is.null(mapped_first)) {
    if (nzchar(rest)) file.path(mapped_first, rest) else mapped_first
  } else {
    normalized
  }
  file.path(project_root, mapped)
}

source(tbml_resolve("R/00_load_all.R"))
config <- read_config("config/config.yml")
ensure_project_dirs(config)

log_path <- project_path("logs", sprintf("run_all_%s.log", mode))
dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

message(sprintf("[%s] run_all start: mode=%s", Sys.time(), mode))

run_script <- function(path) {
  full <- project_path(path)
  if (!file.exists(full)) stop(sprintf("Required script not found: %s", path), call. = FALSE)
  message(sprintf("[%s] START %s", Sys.time(), path))
  source(full, local = new.env(parent = globalenv()))
  message(sprintf("[%s] END %s", Sys.time(), path))
}

if (mode == "prepare_data") {
  for (s in c(
    "scripts/00_setup_packages.R",
    "scripts/01_pull_curatedTBData_metadata.R",
    "scripts/02_extract_sample_metadata_and_harmonize.R",
    "scripts/03_build_cohort_families.R",
    "scripts/04_feasibility_gate_and_lock_splits.R",
    "scripts/05_extract_and_preprocess_expression.R",
    "scripts/06_build_task_matrices.R",
    "scripts/07_benchmark_reference_signatures.R"
  )) run_script(s)
} else if (mode == "optimize_final") {
  for (s in c(
    "scripts/00_setup_packages.R",
    "scripts/08_ml_final_optimization.R",
    "scripts/09_verify_final_results.R"
  )) run_script(s)
} else if (mode == "enhance_final") {
  for (s in c(
    "scripts/09_verify_final_results.R",
    "scripts/10_ml_final_ensemble_enhancement.R",
    "scripts/11_ml_final_operating_points.R",
    "scripts/12_generalization_stress_test.R",
    "scripts/13_generalization_stress_test_enhancements.R"
  )) run_script(s)
} else if (mode == "stress_test") {
  run_script("scripts/12_generalization_stress_test.R")
  run_script("scripts/13_generalization_stress_test_enhancements.R")
} else if (mode == "manuscript_final") {
  run_script("scripts/14_refresh_manuscript_assets.R")
}

message(sprintf("[%s] run_all completed: mode=%s", Sys.time(), mode))
