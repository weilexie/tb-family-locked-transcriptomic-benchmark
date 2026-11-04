project_root <- Sys.getenv("TBML_PROJECT_ROOT", unset = getwd())
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

tbml_directory_aliases <- c(
  config = "01_config",
  data_meta = "02_data_meta",
  data_curated = "03_data_curated",
  R = "04_reusable_R_functions",
  scripts = "05_pipeline_scripts",
  tests = "06_validation_tests",
  results = "07_results",
  manuscript = "08_manuscript",
  submission = "09_submission",
  figures_submission = "09_submission/figures",
  plot_data_submission = "09_submission/plot_data",
  plot_inputs_submission = "09_submission/plot_inputs",
  captions_submission = "09_submission/captions",
  quality_control = "09_submission/quality_control",
  report = "10_report",
  logs = "11_logs"
)

tbml_path_aliases <- c(
  "data_curated/by_study" = "03_data_curated/01_by_study",
  "data_curated/task_matrices" = "03_data_curated/02_task_matrices",
  "results/tables" = "07_results/01_data_preparation",
  "results/benchmark" = "07_results/02_reference_signature_benchmark",
  "results/ml_final" = "07_results/03_ml_final",
  "manuscript/tables" = "08_manuscript/01_tables",
  "manuscript/figure_data" = "08_manuscript/02_figure_data",
  "manuscript/figures" = "08_manuscript/03_figures",
  "manuscript/captions" = "08_manuscript/04_captions",
  "tests/testthat" = "06_validation_tests/01_testthat"
)

tbml_map_path <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  nested_aliases <- tbml_path_aliases[order(nchar(names(tbml_path_aliases)), decreasing = TRUE)]
  for (alias in names(nested_aliases)) {
    if (normalized == alias || startsWith(normalized, paste0(alias, "/"))) {
      rest <- sub(paste0("^", alias, "/?"), "", normalized)
      return(if (nzchar(rest)) file.path(nested_aliases[[alias]], rest) else nested_aliases[[alias]])
    }
  }
  first <- sub("/.*$", "", normalized)
  rest <- sub("^[^/]+/?", "", normalized)
  mapped_first <- tbml_directory_aliases[[first]]
  if (!is.null(mapped_first)) {
    if (nzchar(rest)) file.path(mapped_first, rest) else mapped_first
  } else {
    normalized
  }
}

tbml_resolve <- function(path) {
  file.path(project_root, tbml_map_path(path))
}

extra_libs <- Sys.getenv("TBML_R_LIBS", unset = "")
extra_libs <- if (nzchar(extra_libs)) strsplit(extra_libs, .Platform$path.sep, fixed = TRUE)[[1]] else character()
extra_libs <- normalizePath(extra_libs[nzchar(extra_libs) & dir.exists(extra_libs)], winslash = "/", mustWork = FALSE)
if (length(extra_libs) > 0) .libPaths(unique(c(extra_libs, .libPaths())))

source_files <- list.files(tbml_resolve("R"), pattern = "\\.R$", full.names = TRUE)
source_files <- source_files[basename(source_files) != "00_load_all.R"]
for (f in sort(source_files)) {
  source(f, local = .GlobalEnv)
}
