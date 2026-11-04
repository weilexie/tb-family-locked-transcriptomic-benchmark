sf_plot_figure3 <- function() {
  primary <- sf_label_curves(sf_read_tsv("plot_inputs_submission/Figure5_primary_curves.tsv"))
  expanded <- sf_label_curves(sf_read_tsv("plot_inputs_submission/Figure5_expanded_curves.tsv"))
  curves <- rbind(primary, expanded)
  metrics <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure5_calibration_summary.tsv"))
  metrics$task_group <- sf_short_task(metrics$task)
  metrics$display_model <- sf_display_model(metrics$display_name)
  metrics <- unique(metrics)
  cal <- sf_calibration_bins()
  main_model <- c(
    "Primary" = "Signature-score SVM",
    "Expanded" = "Ridge (penalized-linear)"
  )
  model_cols <- c(
    "Signature-score SVM" = sf_cols$primary,
    "Penalized-linear primary" = sf_cols$sparse,
    "Ridge (penalized-linear)" = "#6F8EC7",
    "Expanded random forest" = sf_cols$warning,
    "Top-10 signature mean ensemble" = "#82B9AA",
    "Top-5 signature mean ensemble" = "#82B9AA"
  )
  task_legend_models <- function(task_group) {
    if (task_group == "Primary") {
      c("Signature-score SVM", "Top-10 signature mean ensemble")
    } else {
      c("Ridge (penalized-linear)", "Top-10 signature mean ensemble")
    }
  }
  curves$display_model <- factor(curves$display_model, levels = names(model_cols))
  cal$display_model <- factor(cal$display_model, levels = names(model_cols))
  metrics$display_model <- factor(metrics$display_model, levels = names(model_cols))
  metric_value <- function(task_group, model, field) {
    idx <- which(metrics$task_group == task_group & as.character(metrics$display_model) == model)
    if (length(idx) == 0) return(NA_real_)
    metrics[[field]][idx[1]]
  }

  plot_roc <- function(task_group, panel_letter, show_legend = FALSE) {
    dat <- curves[curves$task_group == task_group & curves$curve_type == "ROC", ]
    auc <- metric_value(task_group, main_model[[task_group]], "AUROC")
    ggplot2::ggplot(dat, ggplot2::aes(x = fpr, y = tpr, colour = display_model)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = sf_cols$negative, linewidth = .3) +
      ggplot2::geom_path(linewidth = .72, alpha = .95, show.legend = show_legend) +
      ggplot2::annotate("label", x = .62, y = .10, label = paste0("AUROC ", sf_num(auc)), hjust = 0, size = 2.15, fill = "white", colour = sf_cols$ink, linewidth = .14, label.padding = grid::unit(.08, "lines")) +
      ggplot2::scale_colour_manual(values = model_cols, breaks = task_legend_models(task_group), drop = TRUE, name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(.25, .50, .75, 1.00), labels = scales::number_format(accuracy = .01)) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(x = "False positive rate", y = "True positive rate") +
      sf_theme(7.0) +
      ggplot2::theme(
        legend.position = if (show_legend) c(.68, .18) else "none",
        legend.direction = "vertical",
        legend.background = ggplot2::element_rect(fill = "white", colour = sf_cols$border, linewidth = .18),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.key.width = grid::unit(5.6, "mm"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
        plot.margin = ggplot2::margin(8, 8, 8, 14)
      ) +
      sf_panel_tag(panel_letter)
  }
  plot_pr <- function(task_group, panel_letter, show_legend = FALSE) {
    dat <- curves[curves$task_group == task_group & curves$curve_type == "PR", ]
    prev <- metric_value(task_group, main_model[[task_group]], "prevalence")
    auprc <- metric_value(task_group, main_model[[task_group]], "AUPRC")
    ggplot2::ggplot(dat, ggplot2::aes(x = recall, y = precision, colour = display_model)) +
      ggplot2::geom_hline(yintercept = prev, linetype = 3, colour = sf_cols$muted, linewidth = .35) +
      ggplot2::geom_path(linewidth = .72, alpha = .95, show.legend = show_legend) +
      ggplot2::annotate("label", x = .04, y = .11, label = paste0("AUPRC ", sf_num(auprc), "\nprevalence ", sf_num(prev)), hjust = 0, size = 2.05, fill = "white", colour = sf_cols$ink, linewidth = .14, label.padding = grid::unit(.08, "lines"), lineheight = .9) +
      ggplot2::scale_colour_manual(values = model_cols, breaks = task_legend_models(task_group), drop = TRUE, name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(.25, .50, .75, 1.00), labels = scales::number_format(accuracy = .01)) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
      ggplot2::labs(x = "Recall", y = "Precision") +
      sf_theme(7.0) +
      ggplot2::theme(
        legend.position = if (show_legend) c(.34, .20) else "none",
        legend.direction = "vertical",
        legend.background = ggplot2::element_rect(fill = "white", colour = sf_cols$border, linewidth = .18),
        legend.key.height = grid::unit(2.8, "mm"),
        legend.key.width = grid::unit(5.6, "mm"),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
        plot.margin = ggplot2::margin(8, 8, 8, 10)
      ) +
      sf_panel_tag(panel_letter)
  }
  plot_cal <- function(task_group, panel_letter) {
    dat <- cal[cal$task_group == task_group, ]
    slope <- metric_value(task_group, main_model[[task_group]], "calibration_slope")
    main_dat <- dat[dat$display_model == main_model[[task_group]], , drop = FALSE]
    main_dat <- main_dat[order(main_dat$mean_predicted_probability), , drop = FALSE]
    deviation_fill <- ifelse(task_group == "Primary", sf_cols$failure_light, sf_cols$expanded_light)
    label_fill <- ifelse(task_group == "Primary", sf_cols$failure_light, sf_cols$expanded_light)
    label_colour <- ifelse(task_group == "Primary", sf_cols$failure, sf_cols$expanded)
    ggplot2::ggplot(dat, ggplot2::aes(x = mean_predicted_probability, y = observed_event_fraction, colour = display_model)) +
      ggplot2::geom_ribbon(data = main_dat, ggplot2::aes(x = mean_predicted_probability, ymin = pmin(observed_event_fraction, mean_predicted_probability), ymax = pmax(observed_event_fraction, mean_predicted_probability)), inherit.aes = FALSE, fill = deviation_fill, alpha = .62) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = sf_cols$negative, linewidth = .3) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = observed_ci_low, ymax = observed_ci_high), width = .025, linewidth = .25, alpha = .55, show.legend = FALSE) +
      ggplot2::geom_path(linewidth = .55, alpha = .9, show.legend = FALSE) +
      ggplot2::geom_point(ggplot2::aes(size = n_bin), alpha = .75, show.legend = FALSE) +
      ggplot2::annotate("label", x = .05, y = .94, label = paste0("slope ", sf_num(slope)), hjust = 0, vjust = 1, size = 2.0, linewidth = .15, fill = label_fill, colour = label_colour, fontface = "bold", label.padding = grid::unit(.08, "lines"), lineheight = .9) +
      ggplot2::scale_size(range = c(1.0, 3.2), guide = "none") +
      ggplot2::scale_colour_manual(values = model_cols, breaks = task_legend_models(task_group), drop = TRUE, name = NULL) +
      ggplot2::scale_x_continuous(breaks = c(.25, .50, .75, 1.00), labels = scales::number_format(accuracy = .01)) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      ggplot2::labs(x = "Predicted ATB probability", y = "Observed ATB fraction") +
      sf_theme(7.0) +
      ggplot2::theme(
        legend.position = "none",
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
        plot.margin = ggplot2::margin(8, 8, 8, 10)
      ) +
      sf_panel_tag(panel_letter)
  }
  fig_a <- plot_roc("Primary", "A", show_legend = TRUE)
  fig_b <- plot_pr("Primary", "", show_legend = TRUE)
  fig_c <- plot_cal("Primary", "")
  fig_d <- plot_roc("Expanded", "B", show_legend = TRUE)
  fig_e <- plot_pr("Expanded", "", show_legend = TRUE)
  fig_f <- plot_cal("Expanded", "")
  sf_write_data(curves, "Figure3_curve_data.tsv")
  sf_write_data(metrics, "Figure3_metric_summary.tsv")
  sf_write_data(cal, "Figure3_calibration_bins.tsv")
  fig <- patchwork::wrap_plots(
    A = fig_a,
    B = fig_b,
    C = fig_c,
    D = fig_d,
    E = fig_e,
    F = fig_f,
    design = "
      ABC
      DEF
    ",
    heights = c(1, 1),
    widths = c(1, 1, 1.05),
    guides = "keep"
  ) +
    sf_annotation()
  sf_save(fig, "figures_submission/main/Figure3_external_discrimination_and_calibration", 12.6, 7.35)
}
