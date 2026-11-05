make_external_validation_split <- function(task_meta, group_col, config) {
  groups <- sort(unique(task_meta[[group_col]]))
  if (length(groups) < 2) stop("Need at least two groups to create an external split.", call. = FALSE)
  tab <- table(task_meta[[group_col]], task_meta$task_label)
  eligible <- rownames(tab)[rowSums(tab > 0) >= 2]
  if (length(eligible) == 0) eligible <- groups
  target_fraction <- config$external_validation$target_fraction_families %||% 0.25
  n_ext_target <- max(1, ceiling(length(groups) * target_fraction))
  n_ext_min <- config$external_validation$min_external_families_green %||% 2
  max_ext_preserving_dev <- max(1, length(eligible) - 3)
  n_ext <- min(n_ext_target, max(n_ext_min, max_ext_preserving_dev))
  n_ext <- min(n_ext, length(eligible))
  external <- tail(sort(eligible), n_ext)
  development <- setdiff(groups, external)
  if (length(development) < 2) {
    external <- tail(sort(eligible), 1)
    development <- setdiff(groups, external)
  }
  list(
    task = "primary_ATB_vs_OD",
    split_unit = group_col,
    development_families = development,
    external_validation_families = external,
    external_selection_algorithm = "deterministic_sorted_tail_v1"
  )
}

write_external_lockfile <- function(split, config, input_hash, path = "config/external_validation_lockfile.json") {
  full <- project_path(path)
  overwrite <- isTRUE(config$project$overwrite_external_validation_lockfile)
  if (file.exists(full) && !overwrite) {
    return(read_json(path))
  }
  payload <- list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = config$project$seed,
    task = split$task,
    split_unit = split$split_unit,
    development_families = as.character(split$development_families),
    external_validation_families = as.character(split$external_validation_families),
    external_selection_algorithm = split$external_selection_algorithm,
    input_metadata_sha256 = input_hash,
    config_sha256 = hash_file_or_object("config/config.yml"),
    notes = "Real-data external validation lockfile generated before feature screening; do not overwrite unless config explicitly allows it."
  )
  write_json_atomic(payload, path)
  payload
}

read_external_lockfile <- function(path = "config/external_validation_lockfile.json") {
  read_json(path)
}

assert_no_overlap_between_dev_and_external <- function(split) {
  dev <- split$development_families
  ext <- split$external_validation_families
  overlap <- intersect(dev, ext)
  if (length(overlap) > 0) {
    stop(sprintf("Development/external overlap detected: %s", paste(overlap, collapse = ", ")), call. = FALSE)
  }
  TRUE
}

make_leave_group_out_folds <- function(meta, group_col) {
  groups <- sort(unique(meta[[group_col]]))
  lapply(groups, function(g) {
    list(
      heldout_group = g,
      train_ids = meta$standard_sample_id[meta[[group_col]] != g],
      test_ids = meta$standard_sample_id[meta[[group_col]] == g]
    )
  })
}
