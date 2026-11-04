sf_plot_figureS4 <- function() {
  cal <- sf_sanitize(sf_read_tsv("results/ml_final/calibration_repair_external_summary.tsv"))
  cal$task_label <- factor(
    sf_short_task(cal$task),
    levels = c("Primary", "Expanded"),
    labels = c("Primary task", "Expanded task")
  )
  cal$role_label <- factor(
    ifelse(cal$selected_role == "best_performance", "Best performance", "Penalized-linear"),
    levels = c("Best performance", "Penalized-linear")
  )
  cal$method_label <- dplyr::case_when(
    cal$calibration_method == "raw_uncalibrated" ~ "Raw score",
    cal$calibration_method == "development_oof_platt" ~ "Platt scaling",
    cal$calibration_method == "development_oof_temperature" ~ "Temperature scaling",
    cal$calibration_method == "development_oof_beta" ~ "Beta calibration",
    cal$calibration_method == "development_oof_isotonic" ~ "Isotonic regression",
    TRUE ~ gsub("_", " ", cal$calibration_method)
  )
  method_levels <- c("Raw score", "Platt scaling", "Temperature scaling", "Beta calibration", "Isotonic regression")
  cal$method_label <- factor(cal$method_label, levels = rev(method_levels))
  long <- tidyr::pivot_longer(cal, c("Brier", "ECE", "calibration_slope"), names_to = "metric", values_to = "value")

  metric_limits <- function(x) {
    rng <- range(x, na.rm = TRUE)
    span <- diff(rng)
    if (!is.finite(span) || span <= 0) span <- max(abs(rng), 0.05)
    c(rng[1] - span * 0.12, rng[2] + span * 0.12)
  }

  metric_breaks <- function(metric_id) {
    if (metric_id == "Brier") return(c(0.15, 0.20, 0.25))
    if (metric_id == "ECE") return(c(0.05, 0.10, 0.15))
    if (metric_id == "calibration_slope") return(c(0, 1, 2))
    waiver()
  }

  metric_labels <- function(metric_id) {
    if (metric_id %in% c("Brier", "ECE")) return(function(x) sprintf("%.2f", x))
    if (metric_id == "calibration_slope") return(function(x) sprintf("%.0f", x))
    waiver()
  }

  plot_metric <- function(metric_id, panel_title, ref_line = NULL) {
    d <- long[long$metric == metric_id, , drop = FALSE]
    p_metric <- ggplot2::ggplot(d, ggplot2::aes(x = value, y = method_label, colour = role_label, shape = role_label))
    if (!is.null(ref_line)) {
      p_metric <- p_metric + ggplot2::geom_vline(xintercept = ref_line, linetype = 2, linewidth = .25, colour = sf_cols$muted)
    }
    p_metric +
      ggplot2::geom_point(position = ggplot2::position_dodge(width = .5), size = 2.15, alpha = .95) +
      ggplot2::facet_wrap(~task_label, nrow = 1) +
      ggplot2::scale_colour_manual(values = c("Best performance" = sf_cols$signature, "Penalized-linear" = sf_cols$warning), name = "Selected model") +
      ggplot2::scale_shape_manual(values = c("Best performance" = 16, "Penalized-linear" = 17), name = "Selected model") +
      ggplot2::scale_x_continuous(breaks = metric_breaks(metric_id), labels = metric_labels(metric_id)) +
      ggplot2::coord_cartesian(xlim = metric_limits(d$value)) +
      ggplot2::labs(title = panel_title, x = NULL, y = NULL) +
      sf_theme(6.8) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 7.4, face = "bold"),
        strip.text = ggplot2::element_text(size = 5.8),
        axis.text.y = ggplot2::element_text(size = 5.2),
        axis.text.x = ggplot2::element_text(size = 4.8),
        panel.grid.major.x = ggplot2::element_line(colour = sf_cols$grid, linewidth = .25),
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "bottom"
      )
  }

  p <- (plot_metric("Brier", "Brier score (lower is better)") |
    plot_metric("ECE", "Expected calibration error (lower is better)") |
    plot_metric("calibration_slope", "Calibration slope (ideal = 1)", ref_line = 1)) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  sf_write_data(long, "FigureS4_calibration_repair_comparison.tsv")
  sf_save(p, "figures_submission/supplement/FigureS4_calibration_repair_comparison", 11.1, 5.1)
}
