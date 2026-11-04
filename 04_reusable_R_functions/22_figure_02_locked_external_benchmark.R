sf_plot_figure2 <- function() {
  h2h <- sf_sanitize(sf_read_tsv("manuscript/tables/Table2_reference_signature_locked_external_head_to_head.tsv"))
  h2h <- h2h[h2h$analysis_id == "primary_ATB_vs_OD", , drop = FALSE]
  for (nm in c("evaluation_scope", "comparator_reason", "selection_basis")) {
    if (!nm %in% names(h2h)) h2h[[nm]] <- ""
  }
  for (nm in c("AUROC_ci_low", "AUROC_ci_high", "AUPRC_ci_low", "AUPRC_ci_high")) {
    if (!nm %in% names(h2h)) h2h[[nm]] <- NA_real_
  }
  fc <- sf_sanitize(sf_read_tsv("results/ml_final/family_cluster_bootstrap_ci.tsv"))
  for (i in seq_len(nrow(h2h))) {
    if (grepl("^final_", h2h$model_layer[i])) {
      role <- sub("^final_", "", h2h$model_layer[i])
      hit <- fc[fc$task == h2h$analysis_id[i] & fc$selected_role == role, , drop = FALSE]
      if (nrow(hit) > 0) {
        h2h$AUROC_ci_low[i] <- hit$AUROC_ci_low[1]
        h2h$AUROC_ci_high[i] <- hit$AUROC_ci_high[1]
        h2h$AUPRC_ci_low[i] <- hit$AUPRC_ci_low[1]
        h2h$AUPRC_ci_high[i] <- hit$AUPRC_ci_high[1]
      }
    }
  }
  h2h$display_name <- sf_display_model(h2h$signature_or_model)
  h2h$category <- dplyr::case_when(
    h2h$model_layer == "final_best_performance" ~ "Final selected model",
    h2h$model_layer == "final_best_sparse" ~ "Sparse final model",
    h2h$model_layer == "reference_signature_ensemble" ~ "Signature ensemble",
    grepl("top_locked_external_overlap_free", h2h$comparator_reason, fixed = TRUE) ~ "Top overlap-free signature",
    TRUE ~ "Canonical/development signature"
  )
  h2h$category <- factor(h2h$category, levels = c("Final selected model", "Sparse final model", "Signature ensemble", "Top overlap-free signature", "Canonical/development signature"))
  h2h$AUROC_ci_low_plot <- ifelse(is.na(h2h$AUROC_ci_low), h2h$AUROC, h2h$AUROC_ci_low)
  h2h$AUROC_ci_high_plot <- ifelse(is.na(h2h$AUROC_ci_high), h2h$AUROC, h2h$AUROC_ci_high)
  h2h$has_auroc_ci <- is.finite(h2h$AUROC_ci_low) & is.finite(h2h$AUROC_ci_high) & (h2h$AUROC_ci_high > h2h$AUROC_ci_low)
  h2h <- h2h[order(h2h$AUROC), ]
  h2h$axis_label <- sf_benchmark_axis_label(h2h$display_name, as.character(h2h$category), h2h$n_families, h2h$n_samples)
  h2h$display_name <- factor(h2h$display_name, levels = unique(h2h$display_name))
  h2h$axis_label <- factor(h2h$axis_label, levels = unique(h2h$axis_label))
  h2h$row_id <- seq_len(nrow(h2h))
  shared_y_levels <- levels(h2h$axis_label)
  shared_y_expand <- ggplot2::expansion(mult = c(0.02, 0.04))
  cat_vals <- c("Final selected model" = sf_cols$primary, "Sparse final model" = sf_cols$sparse, "Signature ensemble" = sf_cols$ensemble, "Top overlap-free signature" = sf_cols$signature, "Canonical/development signature" = sf_cols$negative)
  dev_grid <- sf_sanitize(sf_read_tsv("results/ml_final/development_combination_performance.tsv"))
  dev_grid <- dev_grid[is.finite(dev_grid$AUROC), , drop = FALSE]
  dev_grid$task_group <- sf_short_task(dev_grid$task)
  dev_grid$selector_group <- sf_selector_group(dev_grid$selector)
  selector_levels <- c("Stepwise GLM selectors", "Model-based selectors", "Stable gene panels", "Stable genes + signatures", "Reference-signature scores")
  dev_grid$selector_group <- factor(dev_grid$selector_group, levels = selector_levels)
  dev_grid$selector_y <- as.numeric(dev_grid$selector_group)
  group_median <- stats::aggregate(AUROC ~ selector_group + selector_y, dev_grid, median, na.rm = TRUE)
  ext_grid <- sf_sanitize(sf_read_tsv("results/ml_final/all_model_development_external_evaluation.tsv"))
  ext_grid <- ext_grid[ext_grid$task == "primary_ATB_vs_OD" & is.finite(ext_grid$external_AUROC), , drop = FALSE]
  wf_src <- sf_sanitize(sf_read_tsv("results/ml_final/negative_control_family_stratified_permutation.tsv"))
  floor_row <- wf_src[wf_src$task == "primary_ATB_vs_OD" & wf_src$selected_role == "best_performance", , drop = FALSE][1, ]
  final_row <- h2h[h2h$model_layer == "final_best_performance", , drop = FALSE][1, ]
  reference_benchmark_row <- h2h[!grepl("^final", h2h$model_layer) & is.finite(h2h$AUROC), , drop = FALSE]
  reference_benchmark_row <- reference_benchmark_row[which.max(reference_benchmark_row$AUROC), , drop = FALSE]
  oracle <- sf_sanitize(sf_read_tsv("results/ml_final/family_label_confounding_check.tsv"))
  oracle_primary_dev <- oracle[oracle$task == "primary_ATB_vs_OD" & oracle$split == "development", , drop = FALSE][1, ]
  spectrum <- data.frame(
    point_type = "anchor",
    anchor_role = c("chance", "within_family_permutation_benchmark", "final_locked_model", "best_reference_or_ensemble_benchmark", "family_label_oracle_bound"),
    label = c("Chance\n0.500", sprintf("Within-family\npermutation %.3f", floor_row$AUROC_median), sprintf("Final model\n%.3f", final_row$AUROC), sprintf("Reference/ensemble\nbenchmark %.3f", reference_benchmark_row$AUROC), sprintf("Development family-label\noracle %.3f", oracle_primary_dev$family_oracle_AUROC)),
    AUROC = c(0.5, floor_row$AUROC_median, final_row$AUROC, reference_benchmark_row$AUROC, oracle_primary_dev$family_oracle_AUROC),
    source_table = c("reference", "negative_control_family_stratified_permutation.tsv", "Table2_reference_signature_locked_external_head_to_head.tsv", "Table2_reference_signature_locked_external_head_to_head.tsv", "family_label_confounding_check.tsv"),
    stringsAsFactors = FALSE
  )
  spectrum$anchor_class <- dplyr::case_when(
    spectrum$anchor_role %in% c("within_family_permutation_benchmark", "best_reference_or_ensemble_benchmark") ~ "negative_control_or_reference_benchmark",
    spectrum$anchor_role %in% c("chance", "family_label_oracle_bound") ~ "boundary",
    TRUE ~ "final_model"
  )
  dev_points <- data.frame(
    point_type = "development_oof_grid",
    anchor_role = "development_candidate_grid",
    label = "",
    AUROC = dev_grid$AUROC,
    source_table = "development_combination_performance.tsv",
    anchor_class = NA_character_,
    selector_group = as.character(dev_grid$selector_group),
    task = dev_grid$task,
    task_group = dev_grid$task_group,
    development_AUROC = dev_grid$AUROC,
    external_AUROC = NA_real_,
    stringsAsFactors = FALSE
  )
  ext_grid$display_model <- sf_display_model(ext_grid$selector_classifier)
  ext_points <- data.frame(
    point_type = "top_development_models_with_external_predictions",
    anchor_role = "primary_locked_external_prediction_available",
    label = "",
    AUROC = ext_grid$external_AUROC,
    source_table = "all_model_development_external_evaluation.tsv",
    anchor_class = NA_character_,
    selector_group = "Locked external predictions",
    task = ext_grid$task,
    task_group = sf_short_task(ext_grid$task),
    development_AUROC = ext_grid$development_AUROC,
    external_AUROC = ext_grid$external_AUROC,
    stringsAsFactors = FALSE
  )
  spectrum$selector_group <- NA_character_
  spectrum$task <- NA_character_
  spectrum$task_group <- "Primary"
  spectrum$development_AUROC <- NA_real_
  spectrum$external_AUROC <- NA_real_
  spectrum_out <- rbind(spectrum, dev_points, ext_points)
  anchor_y <- max(dev_grid$selector_y, na.rm = TRUE) + 1.45
  external_y <- max(dev_grid$selector_y, na.rm = TRUE) + .55
  spectrum$label_y <- 1 + c(.28, -.28, .28, -.28, .28)
  spectrum$y <- 1
  benchmark_anchors <- spectrum[spectrum$anchor_class == "negative_control_or_reference_benchmark", , drop = FALSE]
  boundary <- spectrum[spectrum$anchor_class == "boundary", , drop = FALSE]
  final_anchor <- spectrum[spectrum$anchor_class == "final_model", , drop = FALSE]
  p3a_anchor <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = .5, xend = oracle_primary_dev$family_oracle_AUROC, y = 1, yend = 1), linewidth = .42, colour = sf_cols$ink, linetype = 2, alpha = .78) +
    ggplot2::geom_segment(ggplot2::aes(x = floor_row$AUROC_median, xend = reference_benchmark_row$AUROC, y = 1, yend = 1), linewidth = 1.05, colour = sf_cols$warning) +
    ggplot2::geom_point(data = boundary, ggplot2::aes(x = AUROC, y = y), shape = 21, size = 2.35, fill = "white", colour = sf_cols$ink, stroke = .35) +
    ggplot2::geom_point(data = benchmark_anchors, ggplot2::aes(x = AUROC, y = y), shape = 21, size = 4.25, fill = sf_cols$warning, colour = sf_cols$ink, stroke = .34) +
    ggplot2::geom_point(data = final_anchor, ggplot2::aes(x = AUROC, y = y), shape = 21, size = 3.35, fill = sf_cols$primary, colour = sf_cols$ink, stroke = .30) +
    ggrepel::geom_label_repel(data = spectrum, ggplot2::aes(x = AUROC, y = y, label = label, fill = anchor_class), nudge_y = spectrum$label_y - 1, min.segment.length = 0, segment.size = .18, size = 2.0, label.padding = grid::unit(.08, "lines"), max.overlaps = Inf, colour = sf_cols$ink, show.legend = FALSE) +
    ggplot2::scale_x_continuous(limits = c(.48, .96), breaks = c(.5, .6, .7, .8, .9), labels = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = 1, labels = "Evidence\nanchors") +
    ggplot2::scale_fill_manual(values = c("boundary" = "white", "negative_control_or_reference_benchmark" = sf_cols$warning_light, "final_model" = sf_cols$primary_light), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(.48, .96), ylim = c(.56, 1.44), clip = "off") +
    ggplot2::labs(x = NULL, y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.0, colour = sf_cols$ink),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(4, 5.5, 0, 5.5)
    ) +
    sf_panel_tag("A")
  p3a_points <- ggplot2::ggplot() +
    ggplot2::geom_point(data = dev_grid, ggplot2::aes(x = AUROC, y = selector_y, colour = selector_group, shape = task_group), alpha = .40, size = 1.10, position = ggplot2::position_jitter(height = .12, width = 0)) +
    ggplot2::geom_segment(data = group_median, ggplot2::aes(x = AUROC, xend = AUROC, y = selector_y - .28, yend = selector_y + .28), linewidth = .55, colour = sf_cols$ink) +
    ggplot2::geom_point(data = ext_grid, ggplot2::aes(x = external_AUROC, y = external_y), shape = 23, fill = sf_cols$signature_light, colour = sf_cols$ink, stroke = .25, alpha = .92, size = 2.35, position = ggplot2::position_jitter(height = .055, width = 0)) +
    ggplot2::scale_y_continuous(
      breaks = c(seq_along(selector_levels), external_y),
      labels = c(selector_levels, sprintf("Locked external\npredictions (n=%d)", nrow(ext_grid)))
    ) +
    ggplot2::scale_colour_manual(values = c(
      "Reference-signature scores" = sf_cols$ensemble,
      "Stable genes + signatures" = sf_cols$expanded,
      "Stable gene panels" = sf_cols$signature,
      "Model-based selectors" = sf_cols$primary,
      "Stepwise GLM selectors" = sf_cols$warning
    ), guide = "none") +
    ggplot2::scale_shape_manual(values = c("Primary" = 16, "Expanded" = 1), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(.48, .96), ylim = c(.35, external_y + .45), clip = "on") +
    ggplot2::labs(x = "AUROC", y = NULL) +
    sf_theme(7.0) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(colour = sf_cols$ink, size = 2.3), nrow = 1)) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 6.0, colour = sf_cols$ink),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, colour = sf_cols$border, linewidth = .34),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      plot.margin = ggplot2::margin(0, 5.5, 3, 5.5)
    )
  p3a <- p3a_anchor / p3a_points + patchwork::plot_layout(heights = c(.34, .82))

  auroc_col_x <- 0.975
  auprc_col_x <- 1.040
  p3rank <- ggplot2::ggplot(h2h, ggplot2::aes(y = axis_label, x = AUROC)) +
    ggplot2::geom_segment(
      data = h2h[h2h$has_auroc_ci, , drop = FALSE],
      ggplot2::aes(x = AUROC_ci_low_plot, xend = AUROC_ci_high_plot, yend = axis_label),
      linewidth = .50,
      colour = sf_cols$ink,
      alpha = .70
    ) +
    ggplot2::geom_point(ggplot2::aes(fill = category), shape = 21, size = 2.55, colour = sf_cols$ink, stroke = .22) +
    ggplot2::geom_text(
      ggplot2::aes(x = auroc_col_x, label = sprintf("%0.3f", AUROC)),
      hjust = 1,
      size = 1.88,
      colour = sf_cols$muted,
      family = "mono"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = auprc_col_x, label = sprintf("%0.3f", AUPRC)),
      hjust = 1,
      size = 1.88,
      colour = sf_cols$muted,
      family = "mono"
    ) +
    ggplot2::annotate(
      "text",
      x = auroc_col_x,
      y = Inf,
      label = "AUROC",
      hjust = 1,
      vjust = 1.35,
      size = 1.95,
      fontface = "bold",
      colour = sf_cols$muted,
      family = "mono"
    ) +
    ggplot2::annotate(
      "text",
      x = auprc_col_x,
      y = Inf,
      label = "AUPRC",
      hjust = 1,
      vjust = 1.35,
      size = 1.95,
      fontface = "bold",
      colour = sf_cols$muted,
      family = "mono"
    ) +
    ggplot2::scale_fill_manual(values = cat_vals, guide = "none") +
    ggplot2::scale_y_discrete(limits = shared_y_levels, expand = shared_y_expand) +
    ggplot2::scale_x_continuous(breaks = c(0.6, 0.7, 0.8), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::coord_cartesian(xlim = c(0.55, 1.06), clip = "off") +
    ggplot2::labs(x = "AUROC on locked external primary task", y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.y = ggplot2::element_text(size = 5.15, lineheight = .82),
      plot.tag = ggplot2::element_text(face = "bold", size = 9.6, colour = sf_cols$ink, hjust = 0, vjust = 1),
      plot.margin = ggplot2::margin(6, 7, 4, 5.5)
    ) +
    sf_panel_tag("B", x = 0.00, y = 1.00)

  delta <- sf_sanitize(sf_read_tsv("manuscript/tables/Supplementary_Table_paired_family_bootstrap_model_differences.tsv"))
  delta <- delta[delta$task == "primary_ATB_vs_OD" & delta$model_b %in% h2h$signature_or_model, , drop = FALSE]
  delta$comparator <- sf_display_model(delta$model_b)
  axis_lookup <- unique(h2h[, c("display_name", "axis_label"), drop = FALSE])
  delta <- merge(delta, axis_lookup, by.x = "comparator", by.y = "display_name", all.x = TRUE)
  delta$axis_label <- factor(delta$axis_label, levels = levels(h2h$axis_label))
  delta <- delta[!is.na(delta$axis_label), , drop = FALSE]
  delta_xlim <- range(c(delta$delta_AUROC_ci_low, delta$delta_AUROC_ci_high), na.rm = TRUE)
  delta_xlim[1] <- floor((delta_xlim[1] - 0.01) / 0.05) * 0.05
  delta_xlim[2] <- ceiling((delta_xlim[2] + 0.01) / 0.05) * 0.05
  p3c <- ggplot2::ggplot(delta, ggplot2::aes(y = axis_label, x = delta_AUROC)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = sf_cols$muted, linewidth = .22) +
    ggplot2::geom_segment(ggplot2::aes(x = delta_AUROC_ci_low, xend = delta_AUROC_ci_high, yend = axis_label), linewidth = .50, colour = sf_cols$ink) +
    ggplot2::geom_point(ggplot2::aes(fill = delta_AUROC < 0), shape = 21, size = 2.35, colour = sf_cols$ink, stroke = .22) +
    ggplot2::scale_fill_manual(values = c("TRUE" = sf_cols$warning, "FALSE" = sf_cols$primary_light), guide = "none") +
    ggplot2::scale_y_discrete(limits = shared_y_levels, drop = FALSE, expand = shared_y_expand) +
    ggplot2::coord_cartesian(xlim = delta_xlim, clip = "on") +
    ggplot2::labs(x = "Delta AUROC", y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 9.6, colour = sf_cols$ink, hjust = 0, vjust = 1),
      plot.margin = ggplot2::margin(6, 5.5, 4, 20)
    ) +
    sf_panel_tag("C", x = -0.075, y = 1.00)

  wf <- sf_sanitize(sf_read_tsv("results/ml_final/negative_control_family_stratified_permutation.tsv"))
  wf$task_label <- sf_short_task(wf$task)
  wf$role_label <- paste(wf$task_label, ifelse(wf$selected_role == "best_performance", "best", "sparse"), sep = "\n")
  wf$role_label <- factor(wf$role_label, levels = wf$role_label[order(wf$task_label, wf$AUROC_median)])
  wf$display_value <- sprintf("%.3f", wf$AUROC_median)
  p3d <- ggplot2::ggplot(wf, ggplot2::aes(x = role_label, y = AUROC_median, fill = task_label)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, linewidth = .22, colour = sf_cols$muted) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = AUROC_q025, ymax = AUROC_q975), linewidth = .54, colour = sf_cols$ink) +
    ggplot2::geom_point(shape = 21, size = 3.0, colour = sf_cols$ink, stroke = .24) +
    ggplot2::geom_text(ggplot2::aes(label = display_value), vjust = -1.0, size = 2.05, colour = sf_cols$ink) +
    ggplot2::scale_fill_manual(values = c("Primary" = sf_cols$failure_light, "Expanded" = sf_cols$expanded_light), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0.495, 0.62), clip = "off") +
    ggplot2::labs(x = NULL, y = "AUROC median") +
    sf_theme(7.0) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, lineheight = .86, size = 5.8),
      panel.border = ggplot2::element_rect(fill = NA, colour = sf_cols$border, linewidth = .34),
      plot.tag = ggplot2::element_text(face = "bold", size = 9.6, colour = sf_cols$ink, hjust = 0, vjust = 1),
      plot.margin = ggplot2::margin(6, 5.5, 4, 5.5)
    ) +
    sf_panel_tag("D", x = 0.00, y = 1.00)

  sf_write_data(spectrum_out, "Figure2_discrimination_spectrum.tsv")
  sf_write_data(h2h, "Figure2_scoreboard.tsv")
  sf_write_data(delta, "Figure2_delta.tsv")
  sf_write_data(wf, "Figure2_within_family_negative_control.tsv")
  fig <- patchwork::wrap_plots(
    A = p3a,
    B = p3rank,
    C = p3c,
    D = p3d,
    design = "
      AAA
      BCD
    ",
    heights = c(1.08, 1.15),
    widths = c(1.18, .84, .64),
    guides = "keep"
  ) +
    sf_annotation()
  sf_save(fig, "figures_submission/main/Figure2_locked_external_benchmark", 12.4, 8.8)
}
