`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

project_directory_aliases <- function() {
  c(
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
}

project_path_aliases <- function() {
  c(
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
}

map_project_relative_path <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  nested_aliases <- project_path_aliases()
  nested_aliases <- nested_aliases[order(nchar(names(nested_aliases)), decreasing = TRUE)]
  for (alias in names(nested_aliases)) {
    if (normalized == alias || startsWith(normalized, paste0(alias, "/"))) {
      rest <- sub(paste0("^", alias, "/?"), "", normalized)
      return(if (nzchar(rest)) file.path(nested_aliases[[alias]], rest) else nested_aliases[[alias]])
    }
  }
  aliases <- project_directory_aliases()
  first <- sub("/.*$", "", normalized)
  rest <- sub("^[^/]+/?", "", normalized)
  mapped_first <- if (first %in% names(aliases)) aliases[[first]] else NULL
  if (!is.null(mapped_first)) {
    if (nzchar(rest)) file.path(mapped_first, rest) else mapped_first
  } else {
    normalized
  }
}

project_path_one <- function(path) {
  if (!is.character(path) || length(path) != 1) stop("project_path_one expects one path string.", call. = FALSE)
  if (grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]{2}|^/", path)) return(path)
  file.path(Sys.getenv("TBML_PROJECT_ROOT", unset = getwd()), map_project_relative_path(path))
}

project_path <- function(...) {
  args <- list(...)
  if (length(args) == 0) return(Sys.getenv("TBML_PROJECT_ROOT", unset = getwd()))
  if (length(args) == 1 && length(args[[1]]) > 1) {
    return(vapply(args[[1]], project_path_one, character(1), USE.NAMES = FALSE))
  }
  path <- do.call(file.path, args)
  vapply(path, project_path_one, character(1), USE.NAMES = FALSE)
}

read_config <- function(path = "config/config.yml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read config/config.yml", call. = FALSE)
  }
  yaml::read_yaml(project_path(path))
}

ensure_project_dirs <- function(config = NULL) {
  dirs <- c(
    "config", "R", "scripts", "tests/testthat", "data_meta",
    "data_curated", "logs", "results", "results/tables",
    "results/benchmark", "results/ml_final", "manuscript/tables",
    "manuscript/figure_data", "manuscript/figures", "manuscript/captions",
    "submission", "report"
  )
  for (d in dirs) dir.create(project_path(d), recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

write_tsv_atomic <- function(x, path) {
  path <- project_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = if (grepl("\\.gz$", path)) ".tsv.gz" else ".tsv")
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(as.data.frame(x), tmp, sep = "\t", na = "NA")
  } else {
    con <- if (grepl("\\.gz$", tmp)) gzfile(tmp, "wt") else file(tmp, "wt")
    on.exit(close(con), add = TRUE)
    write.table(as.data.frame(x), con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

mlfinal_sanitize_terms <- function(x) {
  x <- as.data.frame(x)
  names(x) <- gsub("published", "reference", names(x), ignore.case = TRUE)
  for (nm in names(x)) {
    if (is.factor(x[[nm]])) x[[nm]] <- as.character(x[[nm]])
    if (is.character(x[[nm]])) {
      x[[nm]] <- gsub("published_signature_scores", "reference_signature_scores", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("published_signature_ensemble", "reference_signature_ensemble", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("published_signature", "reference_signature", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("published_signatures", "reference_signatures", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("Published signature", "Reference signature", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("published signature", "reference signature", x[[nm]], fixed = TRUE)
      x[[nm]] <- gsub("published", "reference", x[[nm]], ignore.case = TRUE)
    }
  }
  x
}

mlfinal_write_display_tsv <- function(x, path) {
  write_tsv_atomic(mlfinal_sanitize_terms(x), path)
}

read_tsv <- function(path) {
  path <- project_path(path)
  if (requireNamespace("data.table", quietly = TRUE)) {
    as.data.frame(data.table::fread(file = path, sep = "\t", data.table = FALSE))
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

write_json_atomic <- function(x, path, auto_unbox = TRUE, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write JSON", call. = FALSE)
  }
  path <- project_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".json")
  jsonlite::write_json(x, tmp, auto_unbox = auto_unbox, pretty = pretty, null = "null")
  if (file.exists(path)) unlink(path)
  file.rename(tmp, path)
  invisible(path)
}

read_json <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to read JSON", call. = FALSE)
  }
  jsonlite::read_json(project_path(path), simplifyVector = TRUE)
}

hash_file_or_object <- function(x) {
  if (is.character(x) && length(x) == 1 && file.exists(project_path(x))) {
    return(unname(tools::md5sum(project_path(x))))
  }
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    return(unname(tools::md5sum(x)))
  }
  tmp <- tempfile(fileext = ".rds")
  saveRDS(x, tmp)
  on.exit(unlink(tmp), add = TRUE)
  unname(tools::md5sum(tmp))
}

write_lines <- function(lines, path) {
  path <- project_path(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

gzip_tsv_matrix <- function(mat, path) {
  df <- data.frame(gene = rownames(mat), mat, check.names = FALSE)
  write_tsv_atomic(df, path)
}
