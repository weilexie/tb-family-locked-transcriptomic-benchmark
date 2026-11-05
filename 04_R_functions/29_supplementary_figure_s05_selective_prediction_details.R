sf_plot_figureS5 <- function() {
  risk <- sf_sanitize(sf_read_tsv("results/ml_final/selective_prediction_risk_coverage.tsv"))
  risk$task_label <- factor(
    sf_short_task(risk$task),
    levels = c("Primary", "Expanded"),
    labels = c("Primary task", "Expanded task")
  )
  risk$role_label <- factor(
    ifelse(risk$selected_role == "best_performance", "Best performance model", "Sparse model"),
    levels = c("Best performance model", "Sparse model")
  )
  risk$panel_label <- paste(risk$task_label, risk$role_label, sep = " / ")
  risk$.risk_row <- seq_len(nrow(risk))
  risk$coverage_bin <- cut(risk$external_coverage, breaks = seq(0, 1, by = .05), include.lowest = TRUE)
  env <- do.call(rbind, lapply(split(risk, paste(risk$task_label, risk$role_label, risk$coverage_bin, sep = "__")), function(d) {
    data.frame(
      task_label = d$task_label[1],
      role_label = d$role_label[1],
      panel_label = d$panel_label[1],
      coverage = mean(d$external_coverage, na.rm = TRUE),
      error_median = stats::median(d$external_error_rate_among_answered, na.rm = TRUE),
      error_low = min(d$external_error_rate_among_answered, na.rm = TRUE),
      error_high = max(d$external_error_rate_among_answered, na.rm = TRUE),
      n_strategies = nrow(d),
      stringsAsFactors = FALSE
    )
  }))
  selected <- risk[risk$is_selected_strategy %in% TRUE, , drop = FALSE]
  selected$selected_strategy <- ifelse(is.na(selected$selected_strategy), "selected", selected$selected_strategy)
  selected$selected_strategy_label <- dplyr::case_when(
    selected$selected_strategy == "dev_balanced_accuracy_max" ~ "Maximum balanced accuracy",
    selected$selected_strategy == "dev_min_error_at_coverage_0.70" ~ "Minimum error at 70% coverage",
    selected$selected_strategy == "dev_sensitivity_ge_0.90_max_coverage" ~ "Maximum coverage with sensitivity >= 90%",
    selected$selected_strategy == "dev_specificity_ge_0.70_max_coverage" ~ "Maximum coverage with specificity >= 70%",
    TRUE ~ gsub("_", " ", selected$selected_strategy)
  )
  selected$selected_strategy_label <- factor(
    selected$selected_strategy_label,
    levels = c(
      "Maximum balanced accuracy",
      "Minimum error at 70% coverage",
      "Maximum coverage with sensitivity >= 90%",
      "Maximum coverage with specificity >= 70%"
    )
  )
  p <- ggplot2::ggplot(env, ggplot2::aes(x = coverage, y = error_median, colour = task_label, fill = task_label)) +
    ggplot2::geom_vline(xintercept = .70, linetype = 3, linewidth = .22, colour = sf_cols$muted) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = error_low, ymax = error_high), alpha = .13, colour = NA) +
    ggplot2::geom_line(linewidth = .72) +
    ggplot2::geom_point(
      data = selected,
      ggplot2::aes(x = external_coverage, y = external_error_rate_among_answered, shape = selected_strategy_label),
      inherit.aes = FALSE,
      size = 2.75,
      fill = sf_cols$warning_light,
      colour = sf_cols$ink,
      stroke = .25
    ) +
    ggplot2::facet_grid(task_label ~ role_label) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_colour_manual(values = c("Primary task" = sf_cols$primary, "Expanded task" = sf_cols$expanded), guide = "none") +
    ggplot2::scale_fill_manual(values = c("Primary task" = sf_cols$primary_light, "Expanded task" = sf_cols$expanded_light), guide = "none") +
    ggplot2::scale_shape_manual(values = c(21, 22, 23, 24), name = "Locked operating rule") +
    ggplot2::labs(x = "Answered coverage", y = "Error among answered") +
    sf_theme(7.0) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(7, 18, 7, 18),
      legend.position = "bottom",
      legend.title.align = .5
    )
  risk_out <- risk[, setdiff(names(risk), ".risk_row"), drop = FALSE]
  sf_write_data(risk_out, "FigureS5_selective_prediction_operating_point_details.tsv")
  sf_write_data(env, "FigureS5_selective_prediction_envelope.tsv")
  sf_save(p, "figures_submission/supplement/FigureS5_selective_prediction_operating_point_details", 8.8, 6.2)
}
