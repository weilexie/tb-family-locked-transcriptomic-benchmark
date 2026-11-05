sf_plot_figure4 <- function() {
  perf <- sf_sanitize(sf_read_tsv("results/ml_final/external_per_family_performance.tsv"))
  perf <- perf[perf$task == "primary_ATB_vs_OD" & perf$selected_role == "best_performance", c("cohort_family", "AUROC", "sensitivity", "specificity", "n"), drop = FALSE]
  names(perf)[names(perf) == "n"] <- "n_samples"
  perf_long <- tidyr::pivot_longer(perf, c("AUROC", "sensitivity", "specificity"), names_to = "metric", values_to = "value")
  perf_long$metric <- dplyr::recode(perf_long$metric, sensitivity = "Sensitivity", specificity = "Specificity")
  fail <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure6_failure_rates.tsv"))
  fail <- fail[, c("cohort_family", "error_label", "error_rate"), drop = FALSE]
  names(fail) <- c("cohort_family", "metric", "value")
  finger <- rbind(perf_long[, c("cohort_family", "metric", "value")], fail)
  metric_raw_levels <- c("AUROC", "Sensitivity", "Specificity", "FN among ATB", "FP among OD")
  finger <- finger[finger$metric %in% metric_raw_levels, ]
  finger <- merge(
    expand.grid(cohort_family = unique(perf$cohort_family), metric = metric_raw_levels, stringsAsFactors = FALSE),
    finger,
    by = c("cohort_family", "metric"),
    all.x = TRUE
  )
  finger$value[is.na(finger$value) & finger$metric %in% c("FN among ATB", "FP among OD")] <- 0
  finger$metric <- factor(finger$metric, levels = metric_raw_levels)
  finger$metric_label <- factor(
    as.character(finger$metric),
    levels = metric_raw_levels,
    labels = c("AUROC", "Sensitivity", "Specificity", "FN rate\namong ATB", "FP rate\namong OD")
  )
  finger$family_short <- sf_family_short(finger$cohort_family)
  fp_rank <- finger[as.character(finger$metric) == "FP among OD", c("family_short", "value"), drop = FALSE]
  finger$family_short <- factor(finger$family_short, levels = fp_rank$family_short[order(fp_rank$value, fp_rank$family_short)])
  finger$label <- ifelse(is.na(finger$value), "NA", ifelse(finger$value < .995, sprintf("%.2f", finger$value), "1.00"))
  finger$is_failure_metric <- as.character(finger$metric) %in% c("FN among ATB", "FP among OD")
  finger$dot_fill <- ifelse(
    finger$is_failure_metric,
    sf_mix_col("white", sf_cols$failure, finger$value),
    sf_mix_col("#F4F6F7", sf_cols$expanded, finger$value)
  )
  finger$metric_family <- ifelse(finger$is_failure_metric, "Error rate", "Discrimination")
  finger$value_size <- pmax(finger$value, .03)
  error_profile <- finger[finger$metric %in% c("FN among ATB", "FP among OD"), , drop = FALSE]
  error_profile$error_label <- factor(
    as.character(error_profile$metric),
    levels = c("FN among ATB", "FP among OD"),
    labels = c("FN rate among ATB", "FP rate among OD")
  )
  error_profile$display_value <- scales::percent(error_profile$value, accuracy = 1)
  error_rank <- stats::aggregate(value ~ family_short, error_profile, max, na.rm = TRUE)
  error_profile$family_short <- factor(error_profile$family_short, levels = error_rank$family_short[order(error_rank$value)])
  collapse_rows <- data.frame(
    family_short = factor(c("GSE73408", "GSE89403"), levels = levels(error_profile$family_short)),
    stringsAsFactors = FALSE
  )
  collapse_rows <- collapse_rows[!is.na(collapse_rows$family_short), , drop = FALSE]
  p5a <- ggplot2::ggplot(finger, ggplot2::aes(x = metric_label, y = family_short)) +
    ggplot2::geom_rect(data = collapse_rows, ggplot2::aes(xmin = .5, xmax = length(metric_raw_levels) + .5, ymin = as.numeric(family_short) - .46, ymax = as.numeric(family_short) + .46), inherit.aes = FALSE, fill = sf_cols$failure_light, alpha = .28, colour = NA) +
    ggplot2::geom_point(ggplot2::aes(fill = metric_family, size = value_size, alpha = value_size), shape = 21, colour = sf_cols$ink, stroke = .22) +
    ggplot2::scale_fill_manual(values = c("Discrimination" = sf_cols$expanded, "Error rate" = sf_cols$failure), name = NULL) +
    ggplot2::scale_size_continuous(
      range = c(2.1, 5.8),
      breaks = c(.25, .50, .75, 1.00),
      labels = c("0.25", "0.50", "0.75", "1.00"),
      name = "Dot size"
    ) +
    ggplot2::scale_alpha_continuous(range = c(.30, .92), guide = "none") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, nrow = 1, override.aes = list(size = 3.0, alpha = .9)),
      size = ggplot2::guide_legend(order = 2, nrow = 1, override.aes = list(shape = 21, fill = "white", colour = sf_cols$ink, alpha = 1))
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", lineheight = .88),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = ggplot2::margin(1, 0, 0, 0),
      legend.key.width = grid::unit(4.5, "mm")
    ) +
    sf_panel_tag("A")

  loo <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure6_leave_one_family_out.tsv"))
  loo <- loo[loo$task == "primary_ATB_vs_OD", ]
  loo$family_short <- sf_family_short(loo$removed_family)
  loo$family_order <- forcats::fct_reorder(loo$family_short, loo$delta_AUROC)
  loo$direction <- ifelse(loo$delta_AUROC < 0, "Performance drops", "Performance increases")
  p5b <- ggplot2::ggplot(loo, ggplot2::aes(x = delta_AUROC, y = family_order)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = sf_cols$negative, linewidth = .35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = delta_AUROC, yend = family_order, colour = direction), linewidth = .55) +
    ggplot2::geom_point(ggplot2::aes(fill = direction), shape = 21, size = 2.5, colour = sf_cols$ink, stroke = .25) +
    ggplot2::scale_colour_manual(values = c("Performance drops" = sf_cols$warning, "Performance increases" = sf_cols$expanded), guide = "none") +
    ggplot2::scale_fill_manual(values = c("Performance drops" = sf_cols$warning, "Performance increases" = sf_cols$expanded), name = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, override.aes = list(size = 2.7))) +
    ggplot2::coord_cartesian(xlim = c(-0.052, 0.036), clip = "off") +
    ggplot2::labs(x = "Delta AUROC after removing family", y = NULL) +
    sf_theme(7.2) +
    ggplot2::theme(
      legend.position = c(.80, .20),
      legend.direction = "vertical",
      legend.background = ggplot2::element_rect(fill = "white", colour = sf_cols$border, linewidth = .18),
      legend.key.height = grid::unit(2.8, "mm"),
      legend.key.width = grid::unit(4.4, "mm"),
      legend.text = ggplot2::element_text(size = 6.2)
    ) +
    sf_panel_tag("B")

  score <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure6_score_distribution.tsv"))
  score$error_type <- factor(score$error_type, levels = c("TN", "FP", "FN", "TP"))
  score$outcome_group <- ifelse(score$error_type %in% c("TP", "TN"), "Correct", "Error")
  score$family_short <- sf_family_short(score$family_label)
  score$family_short <- factor(score$family_short, levels = levels(error_profile$family_short))
  score$true_label <- factor(score$true_label, levels = c("OD", "ATB"))
  sel <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure6_selective_summary.tsv"))[1, ]
  collapse_score <- data.frame(family_short = factor(c("GSE73408", "GSE89403"), levels = levels(score$family_short)), stringsAsFactors = FALSE)
  collapse_score <- collapse_score[!is.na(collapse_score$family_short), , drop = FALSE]
  score_xmax <- length(levels(score$family_short)) + .5
  zone_label_x <- score_xmax + .18
  p5c <- ggplot2::ggplot(score, ggplot2::aes(x = family_short, y = predicted_probability)) +
    ggplot2::annotate("rect", xmin = .5, xmax = score_xmax, ymin = 0, ymax = sel$lower_threshold, fill = sf_cols$signature_light, alpha = .55) +
    ggplot2::annotate("rect", xmin = .5, xmax = score_xmax, ymin = sel$lower_threshold, ymax = sel$upper_threshold, fill = sf_cols$warning_light, alpha = .40) +
    ggplot2::annotate("rect", xmin = .5, xmax = score_xmax, ymin = sel$upper_threshold, ymax = 1, fill = sf_cols$primary_light, alpha = .55) +
    ggplot2::geom_rect(data = collapse_score, ggplot2::aes(xmin = as.numeric(family_short) - .5, xmax = as.numeric(family_short) + .5, ymin = 0, ymax = 1), inherit.aes = FALSE, fill = sf_cols$failure, alpha = .10) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = true_label), width = .60, outlier.shape = NA, linewidth = .28, colour = sf_cols$ink, position = ggplot2::position_dodge(width = .70), alpha = .62) +
    ggplot2::geom_point(ggplot2::aes(fill = true_label), shape = 21, size = .54, alpha = .34, colour = sf_cols$ink, stroke = .03, position = ggplot2::position_jitterdodge(jitter.width = .08, jitter.height = 0, dodge.width = .70)) +
    ggplot2::geom_hline(yintercept = unique(score$locked_threshold)[1], linetype = 2, linewidth = .35, colour = sf_cols$ink) +
    ggplot2::geom_hline(yintercept = c(sel$lower_threshold, sel$upper_threshold), linetype = 3, linewidth = .35, colour = sf_cols$warning) +
    ggplot2::annotate("text", x = zone_label_x, y = sel$lower_threshold / 2, label = paste0("Answer OD\n< ", sf_num(sel$lower_threshold)), hjust = 1, vjust = .5, size = 2.0, colour = sf_cols$muted, lineheight = .92) +
    ggplot2::annotate("text", x = zone_label_x, y = mean(c(sel$lower_threshold, sel$upper_threshold)), label = paste0("Abstain\n", sf_num(sel$lower_threshold), "-", sf_num(sel$upper_threshold)), hjust = 1, vjust = .5, size = 2.0, colour = sf_cols$muted, lineheight = .92) +
    ggplot2::annotate("text", x = zone_label_x, y = sel$upper_threshold + (1 - sel$upper_threshold) / 2, label = paste0("Answer ATB\n> ", sf_num(sel$upper_threshold)), hjust = 1, vjust = .5, size = 2.0, colour = sf_cols$muted, lineheight = .92) +
    ggplot2::annotate("text", x = zone_label_x, y = unique(score$locked_threshold)[1], label = paste0("Locked ", sf_num(unique(score$locked_threshold)[1])), hjust = 1, vjust = -0.55, size = 1.95, colour = sf_cols$ink) +
    ggplot2::scale_fill_manual(values = c("OD" = sf_cols$signature_light, "ATB" = sf_cols$primary), name = "True class") +
    ggplot2::coord_cartesian(ylim = c(0, 1.045), clip = "off") +
    ggplot2::labs(x = NULL, y = "Predicted ATB probability") +
    sf_theme(7.2) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3.0, alpha = .62))) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, size = 6.2),
      legend.position = c(.16, .965),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_rect(fill = "white", colour = sf_cols$border, linewidth = .18),
      legend.key.height = grid::unit(2.8, "mm"),
      legend.key.width = grid::unit(5.0, "mm"),
      plot.margin = ggplot2::margin(3, 30, 3, 3)
    ) +
    sf_panel_tag("C")

  asym_long <- data.frame(
    selected_strategy = sel$selected_strategy,
    metric = factor(
      c("Coverage", "Abstention", "FN among answered ATB", "FP among answered OD"),
      levels = rev(c("Coverage", "Abstention", "FN among answered ATB", "FP among answered OD"))
    ),
    value = c(
      sel$answered_n / sel$n,
      sel$uncertain_n / sel$n,
      sel$false_negative_rate_among_answered_ATB,
      sel$false_positive_rate_among_answered_comparator
    ),
    stringsAsFactors = FALSE
  )
  asym_long$display_value <- scales::percent(asym_long$value, accuracy = 0.1)
  asym_long$is_focus <- as.character(asym_long$metric) == "FN among answered ATB"
  asym_long$bar_group <- dplyr::case_when(
    as.character(asym_long$metric) == "FN among answered ATB" ~ "ATB miss rate",
    as.character(asym_long$metric) == "FP among answered OD" ~ "OD false positive rate",
    TRUE ~ "Coverage / abstention"
  )
  p5d <- ggplot2::ggplot(asym_long, ggplot2::aes(x = value, y = metric, fill = bar_group)) +
    ggplot2::geom_col(width = .58, colour = "white", linewidth = .32) +
    ggplot2::geom_text(ggplot2::aes(label = display_value), hjust = -0.12, size = 2.25, fontface = "bold", colour = sf_cols$ink) +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, max(asym_long$value, na.rm = TRUE) * 1.22), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = c("Coverage / abstention" = sf_cols$signature_light, "ATB miss rate" = sf_cols$failure, "OD false positive rate" = sf_cols$warning), guide = "none") +
    ggplot2::labs(x = "Share", y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, colour = sf_cols$border, linewidth = .34)) +
    sf_panel_tag("D")

  sf_write_data(finger, "Figure4_family_fingerprint.tsv")
  sf_write_data(loo, "Figure4_leave_one_family_influence.tsv")
  sf_write_data(score, "Figure4_score_distribution.tsv")
  sf_write_data(asym_long, "Figure4_selective_asymmetry.tsv")
  fig <- patchwork::wrap_plots(
    A = p5a,
    B = p5b,
    C = p5c,
    D = p5d,
    design = "
      AB
      CD
    ",
    heights = c(1.0, 1.02),
    widths = c(1.05, .95),
    guides = "keep"
  ) +
    sf_annotation()
  sf_save(fig, "figures_submission/main/Figure4_family_heterogeneity_and_residual_risk", 12.0, 8.15)
}
