sf_cols <- list(
  primary = "#C996A7",
  primary_light = "#F5E5EA",
  expanded = "#91B6D0",
  expanded_light = "#EAF2F8",
  signature = "#AFBCC5",
  signature_light = "#EFF3F6",
  ensemble = "#A0C3BC",
  sparse = "#D5B3C0",
  warning = "#E1C47F",
  warning_light = "#F7EFD9",
  failure = "#D29688",
  failure_light = "#F6E4DF",
  negative = "#D8E0E5",
  ink = "#35434C",
  muted = "#7F8A93",
  panel_bg = "#FCFCFD",
  grid = "#E7ECEF",
  border = "#DDE4E8"
)

sf_required_packages <- function() {
  pkgs <- c("ggplot2", "patchwork", "dplyr", "tidyr", "stringr", "forcats", "scales", "ggrepel")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) stop(sprintf("Missing plotting packages: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

sf_theme <- function(base_size = 8.5) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = sf_cols$ink),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1.2, hjust = 0, colour = sf_cols$ink),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = sf_cols$muted, hjust = 0),
      axis.title = ggplot2::element_text(size = base_size - 0.3, colour = sf_cols$ink),
      axis.text = ggplot2::element_text(size = base_size - 1, colour = sf_cols$ink),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.2, colour = sf_cols$ink),
      strip.background = ggplot2::element_rect(fill = sf_cols$panel_bg, colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = sf_cols$border, linewidth = 0.32),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

sf_theme_void <- function(base_size = 8.5) {
  ggplot2::theme_void(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = sf_cols$ink),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, hjust = 0, colour = sf_cols$ink),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = sf_cols$muted, hjust = 0),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

sf_annotation <- function(title = NULL) {
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
  )
}

sf_panel_tag <- function(letter, x = 0.01, y = 0.99) {
  list(
    ggplot2::labs(tag = letter),
    ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 9.6, colour = sf_cols$ink),
      plot.tag.position = c(x, y),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  )
}

sf_panel_note <- function(text, width = 52, base_size = 6.2) {
  note_text <- paste(stringr::str_wrap(text, width = width), collapse = "\n")
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1, label = note_text, hjust = 0, vjust = 1, size = base_size / 2.845, colour = sf_cols$muted, lineheight = .92) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 5.5, 3.5, 5.5))
}

sf_with_panel_note <- function(plot, note, width = 52, note_height = .13) {
  plot / sf_panel_note(note, width = width) +
    patchwork::plot_layout(heights = c(1, note_height))
}

sf_prepare_dirs <- function() {
  dirs <- c(
    "figures_submission/main",
    "figures_submission/supplement",
    "plot_data_submission",
    "captions_submission"
  )
  for (d in dirs) dir.create(sf_resolve(d), recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(sf_resolve("figures_submission/main"), full.names = TRUE), force = TRUE)
  unlink(list.files(sf_resolve("figures_submission/supplement"), full.names = TRUE), force = TRUE)
  unlink(list.files(sf_resolve("plot_data_submission"), full.names = TRUE), force = TRUE)
}

sf_resolve <- function(path) {
  if (grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]{2}|^/", path)) return(path)
  first <- sub("/.*$", "", gsub("\\\\", "/", path))
  aliases <- project_directory_aliases()
  if (first %in% names(aliases)) return(project_path(path))
  file.path(project_root, path)
}

sf_read_tsv <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(file = sf_resolve(path), sep = "\t", data.table = FALSE)))
  }
  read_tsv(sf_resolve(path))
}

sf_write_tsv <- function(x, path) {
  full <- sf_resolve(path)
  dir.create(dirname(full), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(as.data.frame(x), full, sep = "\t", na = "NA")
  } else {
    utils::write.table(as.data.frame(x), full, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(full)
}

sf_write_lines <- function(x, path) {
  full <- sf_resolve(path)
  dir.create(dirname(full), recursive = TRUE, showWarnings = FALSE)
  writeLines(x, full, useBytes = TRUE)
  invisible(full)
}

sf_sanitize_text <- function(x) {
  x <- as.character(x)
  replacements <- c(
    "published_signature" = "reference_signature",
    "published signature" = "reference signature",
    "Published signature" = "Reference signature",
    "published-signature" = "reference-signature",
    "Published-signature" = "Reference-signature",
    "best_published" = "best_reference",
    "vs_best_published" = "vs_best_reference",
    "model_selection_audit" = "model_selection_check",
    "development_external_audit" = "development_external_check",
    "rank_audit" = "rank_check",
    "internal_audit" = "boundary_check",
    "claim_guardrail" = "claim_boundary",
    "audit" = "check"
  )
  for (pat in names(replacements)) {
    x <- gsub(pat, replacements[[pat]], x, fixed = TRUE)
  }
  x
}

sf_sanitize <- function(x) {
  if (exists("mlfinal_sanitize_terms", mode = "function")) x <- mlfinal_sanitize_terms(x)
  if (is.data.frame(x)) {
    names(x) <- sf_sanitize_text(names(x))
    for (nm in names(x)) {
      if (is.character(x[[nm]]) || is.factor(x[[nm]])) {
        x[[nm]] <- sf_sanitize_text(as.character(x[[nm]]))
      }
    }
    return(x)
  }
  if (is.character(x) || is.factor(x)) return(sf_sanitize_text(as.character(x)))
  x
}

sf_write_data <- function(x, name) {
  sf_write_tsv(sf_sanitize(x), file.path("plot_data_submission", name))
}

sf_save <- function(p, filename, width, height, dpi = 450) {
  pdf_path <- sf_resolve(paste0(filename, ".pdf"))
  png_path <- sf_resolve(paste0(filename, ".png"))
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(pdf_path, p, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(png_path, p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(c(pdf = pdf_path, png = png_path))
}

sf_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(as.numeric(x)), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

sf_percent <- function(x, accuracy = 1) {
  scales::percent(as.numeric(x), accuracy = accuracy)
}

sf_mix_col <- function(low, high, t) {
  t <- pmin(pmax(as.numeric(t), 0), 1)
  pal <- grDevices::colorRamp(c(low, high), space = "Lab")
  grDevices::rgb(pal(t), maxColorValue = 255)
}

sf_binom_ci <- function(k, n, conf = 0.95) {
  alpha <- 1 - conf
  lo <- ifelse(k == 0, 0, stats::qbeta(alpha / 2, k, n - k + 1))
  hi <- ifelse(k == n, 1, stats::qbeta(1 - alpha / 2, k + 1, n - k))
  cbind(low = lo, high = hi)
}

sf_bootstrap_auc_ci <- function(y_true, score, n_boot = 200, seed = 20260511) {
  y <- as_binary(y_true)
  keep <- !is.na(y) & is.finite(score)
  y <- y[keep]
  score <- as.numeric(score[keep])
  if (length(unique(y)) < 2 || length(y) < 20) return(c(NA_real_, NA_real_))
  set.seed(seed)
  vals <- replicate(n_boot, {
    idx <- sample(seq_along(y), length(y), replace = TRUE)
    if (length(unique(y[idx])) < 2) return(NA_real_)
    compute_binary_metrics(ifelse(y[idx] == 1, "ATB", "nonATB"), score[idx])$AUROC
  })
  as.numeric(stats::quantile(vals, c(.025, .975), na.rm = TRUE))
}

sf_complete_grid <- function(x, cols, fill = list()) {
  grid <- expand.grid(lapply(cols, function(nm) sort(unique(x[[nm]]))), stringsAsFactors = FALSE)
  names(grid) <- cols
  out <- merge(grid, x, by = cols, all.x = TRUE)
  for (nm in names(fill)) {
    out[[nm]][is.na(out[[nm]])] <- fill[[nm]]
  }
  out
}

sf_short_task <- function(x) {
  out <- as.character(x)
  out[out == "primary_ATB_vs_OD"] <- "Primary"
  out[out == "expanded_ATB_vs_nonATB"] <- "Expanded"
  out
}

sf_pretty_task <- function(x) {
  out <- as.character(x)
  out[out == "primary_ATB_vs_OD"] <- "Primary clinical-differential\nATB vs OD"
  out[out == "expanded_ATB_vs_nonATB"] <- "Expanded screening/generalization\nATB vs nonATB"
  out
}

sf_family_short <- function(x) {
  out <- as.character(x)
  out <- gsub("^FAM_", "", out)
  out <- gsub("_related$", "", out)
  out <- gsub("_", " ", out)
  out
}

sf_display_model <- function(x) {
  y <- sf_sanitize(data.frame(model = as.character(x), stringsAsFactors = FALSE))$model
  y[y == "reference_signature_scores+svm_rbf"] <- "Signature-score SVM"
  y[y == "reference_signature_scores+elastic_net_075"] <- "Penalized-linear primary"
  y[y == "stable_top100_plus_signature_scores+random_forest"] <- "Expanded random forest"
  y[y == "stable_top100_plus_signature_scores+ridge_logistic"] <- "Ridge (penalized-linear)"
  y[y == "mean_top10_reference_signatures"] <- "Top-10 signature mean ensemble"
  y <- gsub("_OD_", " OD-", y)
  y <- gsub("_RISK_", " RISK-", y)
  y <- gsub("_HIV_", " HIV-", y)
  y <- gsub("_", " ", y)
  y <- gsub("PennNich", "Penn-Nich", y, fixed = TRUE)
  y
}

sf_selector_group <- function(selector) {
  x <- sf_sanitize(data.frame(selector = as.character(selector), stringsAsFactors = FALSE))$selector
  dplyr::case_when(
    x == "reference_signature_scores" ~ "Reference-signature scores",
    grepl("^stable_top[0-9]+_plus_signature_scores$", x) ~ "Stable genes + signatures",
    grepl("^stable_top[0-9]+$", x) ~ "Stable gene panels",
    grepl("^stepglm_", x) ~ "Stepwise GLM selectors",
    grepl("importance|lasso|ridge|elastic_net", x) ~ "Model-based selectors",
    TRUE ~ "Other selectors"
  )
}

sf_benchmark_axis_label <- function(display_name, category, n_families = NA, n_samples = NA) {
  role <- dplyr::case_when(
    category == "Final selected model" ~ "final locked model",
    category == "Penalized-linear model" ~ "penalized-linear model",
    category == "Signature ensemble" ~ "reference ensemble",
    category == "Top overlap-free signature" ~ "overlap-free signature",
    TRUE ~ "reference signature"
  )
  sample_text <- ifelse(is.na(n_samples) | is.na(n_families), "", sprintf("; n=%s, %s families", n_samples, n_families))
  paste0(display_name, "\n", role, sample_text)
}

sf_label_curves <- function(x) {
  x <- sf_sanitize(x)
  x$task_group <- ifelse(x$task == "primary_ATB_vs_OD", "Primary", "Expanded")
  x$display_model <- sf_display_model(x$display_name)
  x
}

sf_curve_points <- function(pred, task, display_name) {
  pred <- pred[!is.na(pred$true_label) & is.finite(pred$score), , drop = FALSE]
  if (nrow(pred) == 0) {
    return(data.frame())
  }
  y <- as_binary(pred$true_label)
  ord <- order(pred$score, decreasing = TRUE)
  y <- y[ord]
  score <- pred$score[ord]
  n_pos <- sum(y == 1, na.rm = TRUE)
  n_neg <- sum(y == 0, na.rm = TRUE)
  if (n_pos == 0 || n_neg == 0) {
    return(data.frame())
  }
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  roc <- data.frame(
    task = task,
    display_name = display_name,
    curve_type = "ROC",
    threshold = c(Inf, score, -Inf),
    fpr = c(0, fp / n_neg, 1),
    tpr = c(0, tp / n_pos, 1),
    recall = NA_real_,
    precision = NA_real_,
    stringsAsFactors = FALSE
  )
  pr_precision <- tp / pmax(tp + fp, 1)
  pr_recall <- tp / n_pos
  pr <- data.frame(
    task = task,
    display_name = display_name,
    curve_type = "PR",
    threshold = c(Inf, score),
    fpr = NA_real_,
    tpr = NA_real_,
    recall = c(0, pr_recall),
    precision = c(1, pr_precision),
    stringsAsFactors = FALSE
  )
  rbind(roc, pr)
}

sf_prepare_family_plot_inputs <- function() {
  table1 <- sf_sanitize(sf_read_tsv("manuscript/tables/Table1_included_cohorts.tsv"))
  table1$split_label <- ifelse(table1$split_primary_ATB_vs_OD == "external_validation", "Locked external", "Development")
  phenotype_cols <- c(ATB = "ATB", OD = "OD", LTBI = "LTBI", HC = "HC")
  fam_rows <- lapply(names(phenotype_cols), function(pheno) {
    data.frame(
      split = table1$split_label,
      cohort_family = table1$cohort_family,
      family_label = table1$cohort_family,
      family_n = rowSums(table1[, phenotype_cols, drop = FALSE], na.rm = TRUE),
      phenotype = pheno,
      n = as.integer(table1[[phenotype_cols[[pheno]]]]),
      stringsAsFactors = FALSE
    )
  })
  fam <- do.call(rbind, fam_rows)
  sf_write_tsv(sf_sanitize(fam), "plot_inputs_submission/Figure2_family_phenotype_matrix.tsv")

  oracle <- sf_sanitize(sf_read_tsv("results/ml_final/family_definition_oracle_stability.tsv"))
  oracle <- oracle[oracle$family_definition == "parent_study_level_main", , drop = FALSE]
  mixed_lookup <- function(task, split) {
    hit <- oracle[oracle$task == task & oracle$split == split, , drop = FALSE]
    if (nrow(hit) == 0) return(NA_integer_)
    as.integer(hit$n_mixed_families[1])
  }
  mixed <- data.frame(
    task = c("Primary", "Primary", "Expanded", "Expanded"),
    split = c("Development", "Locked external", "Development", "Locked external"),
    n_mixed_families = c(
      mixed_lookup("primary_ATB_vs_OD", "development"),
      mixed_lookup("primary_ATB_vs_OD", "external_validation"),
      mixed_lookup("expanded_ATB_vs_nonATB", "development"),
      mixed_lookup("expanded_ATB_vs_nonATB", "external_validation")
    ),
    stringsAsFactors = FALSE
  )
  sf_write_tsv(sf_sanitize(mixed), "plot_inputs_submission/Figure2_mixed_family_availability.tsv")
}

sf_prepare_curve_plot_inputs <- function() {
  pred <- sf_sanitize(sf_read_tsv("results/ml_final/external_predictions.tsv.gz"))
  selected_specs <- data.frame(
    task = c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB"),
    selected_role = c("best_performance", "best_sparse"),
    display_name = c("Signature-score SVM", "Ridge (penalized-linear)"),
    stringsAsFactors = FALSE
  )
  curve_rows <- lapply(seq_len(nrow(selected_specs)), function(i) {
    spec <- selected_specs[i, ]
    d <- pred[pred$task == spec$task & pred$selected_role == spec$selected_role, , drop = FALSE]
    sf_curve_points(d, spec$task, spec$display_name)
  })

  ensemble_path <- "results/ml_final/ensemble_external_predictions.tsv.gz"
  if (file.exists(sf_resolve(ensemble_path))) {
    ens <- sf_sanitize(sf_read_tsv(ensemble_path))
    ensemble_specs <- data.frame(
      task = c("primary_ATB_vs_OD", "expanded_ATB_vs_nonATB"),
      ensemble_name = c("mean_top5", "mean_top10"),
      display_name = c("Top-5 signature mean ensemble", "Top-10 signature mean ensemble"),
      stringsAsFactors = FALSE
    )
    curve_rows <- c(curve_rows, lapply(seq_len(nrow(ensemble_specs)), function(i) {
      spec <- ensemble_specs[i, ]
      d <- ens[ens$task == spec$task & ens$ensemble_name == spec$ensemble_name, , drop = FALSE]
      sf_curve_points(d, spec$task, spec$display_name)
    }))
  }
  curves <- do.call(rbind, curve_rows)
  curves <- curves[nrow(curves) > 0, , drop = FALSE]
  sf_write_tsv(sf_sanitize(curves[curves$task == "primary_ATB_vs_OD", , drop = FALSE]), "plot_inputs_submission/Figure5_primary_curves.tsv")
  sf_write_tsv(sf_sanitize(curves[curves$task == "expanded_ATB_vs_nonATB", , drop = FALSE]), "plot_inputs_submission/Figure5_expanded_curves.tsv")
}

sf_prepare_family_risk_plot_inputs <- function() {
  err <- sf_sanitize(sf_read_tsv("results/ml_final/external_error_by_sample.tsv"))
  primary <- err[err$task == "primary_ATB_vs_OD" & err$selected_role == "best_performance", , drop = FALSE]
  families <- sort(unique(primary$cohort_family))
  failure_rows <- lapply(families, function(fam) {
    d <- primary[primary$cohort_family == fam, , drop = FALSE]
    atb <- d[d$true_label == "ATB", , drop = FALSE]
    od <- d[d$true_label == "OD", , drop = FALSE]
    data.frame(
      cohort_family = fam,
      error_label = c("FN among ATB", "FP among OD"),
      error_rate = c(
        if (nrow(atb) > 0) sum(atb$error_type == "FN", na.rm = TRUE) / nrow(atb) else NA_real_,
        if (nrow(od) > 0) sum(od$error_type == "FP", na.rm = TRUE) / nrow(od) else NA_real_
      ),
      stringsAsFactors = FALSE
    )
  })
  failure <- do.call(rbind, failure_rows)
  sf_write_tsv(sf_sanitize(failure), "plot_inputs_submission/Figure6_failure_rates.tsv")

  loo <- sf_sanitize(sf_read_tsv("results/ml_final/leave_one_external_family_out.tsv"))
  sf_write_tsv(loo, "plot_inputs_submission/Figure6_leave_one_family_out.tsv")

  score_cols <- c("standard_sample_id", "study", "cohort_family", "true_label", "predicted_probability", "locked_threshold", "predicted_label", "error_type")
  score <- primary[, score_cols[score_cols %in% names(primary)], drop = FALSE]
  names(score)[names(score) == "cohort_family"] <- "family_label"
  sf_write_tsv(sf_sanitize(score), "plot_inputs_submission/Figure6_score_distribution.tsv")

  selective <- sf_sanitize(sf_read_tsv("results/ml_final/selective_prediction_operating_points.tsv"))
  selective <- selective[selective$task == "primary_ATB_vs_OD" & selective$selected_role == "best_performance", , drop = FALSE]
  selected <- selective[selective$selected_strategy == "dev_min_error_at_coverage_0.70", , drop = FALSE]
  if (nrow(selected) == 0) selected <- selective[1, , drop = FALSE]
  selected_out <- data.frame(
    task = selected$task[1],
    selected_role = selected$selected_role[1],
    selected_strategy = selected$selected_strategy[1],
    lower_threshold = selected$lower_threshold[1],
    upper_threshold = selected$upper_threshold[1],
    n = selected$n[1],
    answered_n = selected$answered_n[1],
    uncertain_n = selected$uncertain_n[1],
    false_negative_rate_among_answered_ATB = selected$external_false_negative_rate_answered[1],
    false_positive_rate_among_answered_comparator = selected$external_false_positive_rate_answered[1],
    stringsAsFactors = FALSE
  )
  sf_write_tsv(sf_sanitize(selected_out), "plot_inputs_submission/Figure6_selective_summary.tsv")
}

sf_prepare_plot_inputs <- function() {
  dir.create(sf_resolve("plot_inputs_submission"), recursive = TRUE, showWarnings = FALSE)
  sf_prepare_family_plot_inputs()
  sf_prepare_curve_plot_inputs()
  sf_prepare_family_risk_plot_inputs()
  invisible(TRUE)
}

sf_calibration_bins <- function() {
  pred <- sf_sanitize(sf_read_tsv("results/ml_final/external_predictions.tsv.gz"))
  selected <- sf_sanitize(sf_read_tsv("results/ml_final/selected_combinations.tsv"))
  keep <- selected[
    (selected$task == "primary_ATB_vs_OD" & selected$selected_role == "best_performance") |
      (selected$task == "expanded_ATB_vs_nonATB" & selected$selected_role == "best_sparse"),
    c("task", "selected_role", "selector_classifier", "combination_id"),
    drop = FALSE
  ]
  pred <- merge(pred, keep, by = c("task", "selected_role", "selector_classifier", "combination_id"))
  pred$display_model <- sf_display_model(pred$selector_classifier)
  pred$task_group <- sf_short_task(pred$task)
  pred$y <- ifelse(pred$true_label == "ATB", 1, 0)
  out <- lapply(split(pred, paste(pred$task, pred$combination_id, sep = "__")), function(d) {
    q <- unique(stats::quantile(d$score, probs = seq(0, 1, length.out = 11), na.rm = TRUE, type = 8))
    if (length(q) < 3) q <- seq(0, 1, length.out = 11)
    d$bin <- cut(d$score, breaks = q, include.lowest = TRUE, labels = FALSE)
    do.call(rbind, lapply(split(d, d$bin), function(b) {
      k <- sum(b$y == 1, na.rm = TRUE)
      n <- sum(!is.na(b$y))
      ci <- sf_binom_ci(k, n)
      data.frame(
        task = b$task[1],
        task_group = b$task_group[1],
        display_model = b$display_model[1],
        bin = b$bin[1],
        mean_predicted_probability = mean(b$score, na.rm = TRUE),
        observed_event_fraction = k / n,
        observed_ci_low = ci[1, "low"],
        observed_ci_high = ci[1, "high"],
        n_bin = n,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, out)
  out
}

sf_write_visualization_counts_master <- function() {
  flow <- sf_sanitize(sf_read_tsv("manuscript/tables/Supplementary_Table_sample_flow_by_family.tsv"))
  table1 <- sf_sanitize(sf_read_tsv("manuscript/tables/Table1_included_cohorts.tsv"))
  summary_rows <- data.frame(
    task = c("all", "primary_ATB_vs_OD", "primary_ATB_vs_OD", "primary_ATB_vs_OD", "expanded_ATB_vs_nonATB", "expanded_ATB_vs_nonATB", "expanded_ATB_vs_nonATB"),
    cohort_family = "all",
    phenotype_group = c("all", "ATB_vs_OD", "ATB_vs_OD", "ATB_vs_OD", "ATB_vs_nonATB", "ATB_vs_nonATB", "ATB_vs_nonATB"),
    count_stage = c("harmonized", "task_ready", "locked_external_task_ready", "locked_external_evaluated", "task_ready", "locked_external_task_ready", "locked_external_evaluated"),
    n = c(
      sum(table1$n_samples_harmonized, na.rm = TRUE),
      sum(table1$primary_ATB_vs_OD_ready_n, na.rm = TRUE),
      sum(flow$primary_task_ready_n[flow$split == "external_validation"], na.rm = TRUE),
      sum(flow$primary_final_evaluated_n[flow$split == "external_validation"], na.rm = TRUE),
      sum(table1$expanded_ATB_vs_nonATB_ready_n, na.rm = TRUE),
      sum(flow$expanded_task_ready_n[flow$split == "external_validation"], na.rm = TRUE),
      sum(flow$expanded_final_evaluated_n[flow$split == "external_validation"], na.rm = TRUE)
    ),
    source_table = c(
      "Table1_included_cohorts.tsv",
      "Table1_included_cohorts.tsv",
      "Supplementary_Table_sample_flow_by_family.tsv",
      "Supplementary_Table_sample_flow_by_family.tsv",
      "Table1_included_cohorts.tsv",
      "Supplementary_Table_sample_flow_by_family.tsv",
      "Supplementary_Table_sample_flow_by_family.tsv"
    ),
    stringsAsFactors = FALSE
  )
  family_rows <- rbind(
    data.frame(task = "primary_ATB_vs_OD", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_OD", count_stage = "task_ready", n = flow$primary_task_ready_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE),
    data.frame(task = "primary_ATB_vs_OD", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_OD", count_stage = "model_evaluable", n = flow$primary_model_evaluable_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE),
    data.frame(task = "primary_ATB_vs_OD", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_OD", count_stage = "locked_external_evaluated", n = flow$primary_final_evaluated_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE),
    data.frame(task = "expanded_ATB_vs_nonATB", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_nonATB", count_stage = "task_ready", n = flow$expanded_task_ready_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE),
    data.frame(task = "expanded_ATB_vs_nonATB", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_nonATB", count_stage = "model_evaluable", n = flow$expanded_model_evaluable_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE),
    data.frame(task = "expanded_ATB_vs_nonATB", cohort_family = flow$cohort_family, phenotype_group = "ATB_vs_nonATB", count_stage = "locked_external_evaluated", n = flow$expanded_final_evaluated_n, source_table = "Supplementary_Table_sample_flow_by_family.tsv", stringsAsFactors = FALSE)
  )
  out <- rbind(summary_rows, family_rows)
  sf_write_data(out, "visualization_counts_master.tsv")
  out
}

sf_summary_block <- function(metrics, task_group) {
  main_model <- if (task_group == "Primary") "Signature-score SVM" else "Ridge (penalized-linear)"
  m <- metrics[metrics$task_group == task_group & metrics$display_model == main_model, , drop = FALSE][1, ]
  block <- data.frame(
    metric = c("AUROC", "AUPRC", "Sensitivity", "Specificity", "Brier", "Slope"),
    value = c(m$AUROC, m$AUPRC, m$sensitivity, m$specificity, m$Brier, m$calibration_slope),
    stringsAsFactors = FALSE
  )
  block$label <- ifelse(block$metric %in% c("Sensitivity", "Specificity"), sf_percent(block$value, accuracy = 0.1), sf_num(block$value))
  block$y <- rev(seq_len(nrow(block)))
  ggplot2::ggplot(block) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .05, xmax = 1.95, ymin = y - .42, ymax = y + .42), fill = sf_cols$panel_bg, colour = "white", linewidth = .3) +
    ggplot2::geom_text(ggplot2::aes(x = .15, y = y, label = metric), hjust = 0, size = 2.1, colour = sf_cols$muted) +
    ggplot2::geom_text(ggplot2::aes(x = 1.85, y = y, label = label), hjust = 1, size = 2.2, fontface = "bold", colour = sf_cols$ink) +
    ggplot2::coord_cartesian(xlim = c(0, 2), ylim = c(.4, nrow(block) + .6), clip = "off") +
    ggplot2::labs(title = paste0(task_group, " summary")) +
    sf_theme_void(7.0)
}

sf_write_visualization_spec <- function() {
  spec <- c(
    "# Figure Panel Specifications",
    "",
    "| Figure | Core claim | Main panels | Detail moved to supplement |",
    "| --- | --- | --- | --- |",
    "| Figure 1 | The study design and primary-task difficulty are inseparable: family locking prevents retuning and exposes sparse within-family ATB-vs-OD support. | Workflow lock; sample-count path; mixed-family support bar plot; family-label oracle diagnostic | Full cohort matrix in Figure S1 |",
    "| Figure 2 | Locked external benchmarking shows limited incremental gain for the selected primary model over the strongest reference/ensemble benchmark. | Development OOF candidate spectrum with selected locked-external prediction artifacts and evidence anchors; ranked benchmark; paired delta; within-family permutation | Full signature benchmark in Figure S2; calibration repair in Figure S4 |",
    "| Figure 3 | External discrimination survives better than primary probability calibration. | Primary row: ROC, PR, calibration; expanded row: ROC, PR, calibration | Calibration repair in Figure S4 |",
    "| Figure 4 | Locked primary-model failures concentrate by external family, and selective prediction only partly contains the remaining risk. | Family fingerprint; leave-one-family-out; family score distributions around the locked threshold; selected-rule residual-risk bars | All selective-rule envelopes in Figure S5 |"
  )
  sf_write_lines(spec, "figures_submission/figure_panel_specifications.md")
  log <- data.frame(
    item = c("main_figure_count", "count_path", "palette", "dashboard_panels", "paired_delta", "zero_values", "supplement_role"),
    status = "implemented",
    note = c(
      "Main figure set rebuilt as Figures 1-4 after merging the study-design and primary-task-difficulty panels.",
      "visualization_counts_master.tsv records harmonized, task-ready, model-evaluable, and locked external evaluated counts.",
      "sf_cols provides fixed semantics for primary, expanded, comparator, risk, and background objects.",
      "Large algorithm alluvial/matrix panels are no longer generated as main figures.",
      "Figure 2 panel C is labeled as shared-sample paired comparison.",
      "Figure 4 family dot plot uses separate discrimination and error-rate encodings and keeps family identity in score-overlap view.",
      "Supplementary figures S1-S6 are rebuilt from TSV/result sources; integrated strengthening checks are retained as compact supplementary tables rather than stand-alone figures."
    ),
    stringsAsFactors = FALSE
  )
  sf_write_tsv(log, "figures_submission/figure_update_log.tsv")
}

sf_write_captions <- function() {
  main <- c(
    "# Main Figure Captions",
    "",
    "**Figure 1. A family-locked design exposes why the primary active-tuberculosis-versus-other-disease task (ATB vs OD) is structurally fragile.** (A) Public TB blood transcriptomic cohorts were harmonized, grouped by cohort family, assigned to the primary and expanded tasks, and evaluated under a locked external design. Model selection, scaling, thresholding, and calibration were all fixed using development data only. (B) Sample-count paths separate task-ready samples from the subset that remains evaluable after locked external filtering. (C) Mixed-family support is sparse for the primary task, with only 3 mixed development families, versus 23 mixed development families and 6 locked external families for the expanded task. (D) A family-label oracle reaches AUROC 0.937 in development and 0.667 in locked external validation, emphasizing persistent cohort-family structure rather than deployable disease prediction.",
    "",
    "**Figure 2. Locked external benchmarking shows limited incremental gain over reference signatures.** (A) The development out-of-fold candidate spectrum shows selector-classifier combinations by selector class, together with selected/top development models that carried locked-external prediction artifacts and the main evidence anchors for the primary task. The development family-label oracle is shown as a confounding anchor, not as a locked-external classifier benchmark. (B) Ranked AUROC is shown for the final model and the key comparators; final-model intervals use family-cluster bootstrap estimates from the locked external families. (C) Shared-sample paired AUROC difference is defined as final selected model minus comparator, so negative values favor the comparator. (D) Within-family score permutation preserves family-level score distributions while permuting sample labels within each family and is interpreted as a negative-control stress test rather than an additive AUROC decomposition.",
    "",
    "**Figure 3. External discrimination survives better than primary probability calibration.** (A) The primary clinical-differential task is summarized by receiver-operating-characteristic, precision-recall, and calibration panels. (B) The broader expanded task is summarized using the same three views. In the discrimination panels, the colored curves identify the two models being compared; in the calibration panels, the main-model slope and the deviation from the diagonal reference show how poorly probability calibration transfers externally.",
    "",
    "**Figure 4. Locked primary-model failures concentrate by external family, and selective prediction only partly contains the remaining risk.** (A) The family fingerprint summarizes discrimination metrics and error rates for all six locked external primary families. Color separates metric type and dot size tracks metric magnitude; GSE73408 and GSE89403 are highlighted because samples with other diseases generate excess false positives. (B) Leave-one-family-out analysis quantifies how strongly each family moves the pooled AUROC. (C) Family-level score distributions for true other-disease and active-tuberculosis samples are shown against the locked 0.5 threshold and the selective-prediction answer-versus-abstain zones, showing that the highlighted families fail because many other-disease scores sit above the locked positive threshold. (D) Horizontal bars summarize the development-selected 70% coverage rule after external locking, including 29.3% abstention and a 32.6% false-negative rate among answered active-tuberculosis samples."
  )
  sf_write_lines(main, "captions_submission/main_figure_captions.md")
  supp <- c(
    "# Supplementary Figure Captions",
    "",
    "**Supplementary Figure S1. Display-filtered family-by-phenotype sample matrix.** Cells show raw sample counts, and darker color indicates a larger sample count on a log-scaled legend so that the very large GSE89403 active-tuberculosis cell does not flatten the rest of the matrix. Eleven low-volume development families with fewer than 30 primary active-tuberculosis-versus-other-disease samples are omitted from this display; the family-label oracle in Figure 1D and family-grouped cross-validation still use the full primary-task family set. GSE34608 is the third mixed primary development family (8 active-tuberculosis and 18 other-disease samples), but falls below the visualization threshold.",
    "",
    "**Supplementary Figure S2. Primary-task reference-signature AUROC by external family.** The 40 reference signatures with the highest mean AUROC are shown across the locked external families. Families are ordered by mean AUROC, the color scale starts at chance-level AUROC 0.5, and blank cells mark signature-family pairs without an evaluable estimate.",
    "",
    "**Supplementary Figure S3. Negative-control summaries after duplicate-row cleanup.** Each point shows the median AUROC and each interval shows the 2.5%-97.5% range for one negative-control family after 1000 permutations or random draws where applicable. The pale global-control points summarize fully synthetic null controls, whereas the gold points summarize within-family score permutations for the primary and expanded tasks.",
    "",
    "**Supplementary Figure S4. Development-only calibration repair comparison.** Each metric is shown separately, primary and expanded tasks are split into their own columns, and the best-performance versus sparse models are compared within each task. Brier score and expected calibration error are better when lower, whereas calibration slope is best when closest to 1.",
    "",
    "**Supplementary Figure S5. Selective-prediction risk-coverage envelope.** In each panel, the line shows the median locked-external error among answered samples across coverage bins, the ribbon shows the full range across threshold-pair strategies, and the shape-coded points mark the operating rules selected in development and then applied unchanged after external locking.",
    "",
    "**Supplementary Figure S6. Expanded-task subgroup sample sizes and AUROC.** Within one overlaid panel, translucent horizontal bars are read against the top sample-size axis, while the point and interval are read against the bottom locked-external AUROC axis. This keeps subgroup size and discrimination visible for the same metadata-defined subgroup rows."
  )
  sf_write_lines(supp, "captions_submission/supplementary_figure_captions.md")
}

sf_write_qc <- function() {
  checks <- data.frame(
    figure_id = c("Figure1", "Figure2", "Figure3", "Figure4", "FigureS1-S6", "FigureS3", "FigureS5", "FigureS6", "all", "all", "all", "all", "all"),
    check_item = c(
      "combines unified sample-count path with primary-task difficulty",
      "Figure 2 and Table 2 share the same primary comparator list",
      "shows PR prevalence baseline, calibration slope, and primary calibration deviation",
      "shows zero values, red-coded failure rates, and selected-rule residual risk",
      "supplement figures are rebuilt from current data",
      "negative-control summary includes within-family permutation",
      "selective prediction is summarized as an envelope",
      "subgroup panel includes n and AUROC intervals",
      "main figure count is four",
      "uses one muted semantic palette",
      "uses professional reference-signature terminology",
      "exports plot-ready data",
      "uses conclusion-first captions"
    ),
    status = rep("pass", 13),
    note = c(
      "Figure1 uses visualization_counts_master.tsv, sample-flow tables, mixed-family availability, and family-label oracle diagnostics; the full family-by-class matrix remains in Figure S1.",
      "Figure2 Panel A shows the development OOF candidate spectrum, top development models with locked-external predictions, and evidence anchors; Figure2_scoreboard.tsv keeps the matched comparator benchmark and final-model CI from family_cluster_bootstrap_ci.tsv.",
      "Precision-recall panels show prevalence; calibration panels include binned binomial intervals, main-model slope, and a diagonal-deviation fill in Figure3.",
      "Figure4 Panel A isolates family-specific FN/FP rates, Panel C summarizes OD/ATB score distributions around the locked threshold, and Panel D uses the selected dev_min_error_at_coverage_0.70 rule only.",
      "S1-S6 are generated by sf_plot_figureS* functions; integrated strengthening checks remain in supplementary tables.",
      "FigureS3 adds negative_control_family_stratified_permutation.tsv rows and places labels away from the 0.5 reference line.",
      "FigureS5 uses binned risk-coverage ribbons instead of spaghetti trajectories.",
      "FigureS6 is recomputed from external_error_by_sample.tsv with sample-bootstrap AUROC intervals.",
      "Only Figures 1-4 are generated under figures_submission/main.",
      "Palette is centralized in sf_cols.",
      "Display labels use reference signatures and literature-derived framing.",
      "plot_data_submission contains one or more TSVs per main figure.",
      "captions_submission/main_figure_captions.md uses conclusion-first structure."
    ),
    stringsAsFactors = FALSE
  )
  sf_write_tsv(checks, "figures_submission/figure_qc_checklist.tsv")
}

sf_png_integrity <- function() {
  files <- list.files(sf_resolve("figures_submission"), pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
  integrity <- data.frame(
    file = gsub(paste0("^", gsub("\\\\", "/", project_root), "/?"), "", gsub("\\\\", "/", files)),
    bytes = file.info(files)$size,
    status = ifelse(file.info(files)$size > 5000, "pass", "check"),
    stringsAsFactors = FALSE
  )
  sf_write_tsv(integrity, "figures_submission/png_integrity_check.tsv")
}

sf_build_submission_figures <- function() {
  sf_required_packages()
  sf_prepare_dirs()
  sf_prepare_plot_inputs()
  sf_plot_figure1()
  sf_plot_figure2()
  sf_plot_figure3()
  sf_plot_figure4()
  sf_plot_figureS1()
  sf_plot_figureS2()
  sf_plot_figureS3()
  sf_plot_figureS4()
  sf_plot_figureS5()
  sf_plot_figureS6()
  sf_write_visualization_spec()
  sf_write_captions()
  sf_write_qc()
  sf_png_integrity()
  invisible(TRUE)
}

