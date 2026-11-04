sf_plot_figure1 <- function() {
  counts <- sf_write_visualization_counts_master()
  flow <- data.frame(
    step = 1:5,
    label = c(
      "Public TB blood\ncohorts",
      "Phenotype\nharmonization",
      "Cohort-family\nlocking",
      "Two diagnostic\ntasks",
      "Locked external\nbenchmark"
    ),
    sublabel = c("50 studies\n6074 samples", "ATB / OD / LTBI / HC", "families locked\nas units", "primary and\nexpanded", "reference signatures\n+ ML models"),
    x = 1:5,
    y = 1,
    stringsAsFactors = FALSE
  )
  p1a <- ggplot2::ggplot(flow) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x - .42, xmax = x + .42, ymin = .68, ymax = 1.25), fill = "white", colour = sf_cols$border, linewidth = .35) +
    ggplot2::geom_text(ggplot2::aes(x = x, y = 1.03, label = label), size = 2.62, fontface = "bold", lineheight = .88) +
    ggplot2::geom_text(ggplot2::aes(x = x, y = .78, label = sublabel), size = 1.86, colour = sf_cols$muted, lineheight = .88) +
    ggplot2::geom_segment(data = flow[flow$step < 5, ], ggplot2::aes(x = x + .43, xend = x + .57, y = 1, yend = 1), arrow = grid::arrow(length = grid::unit(2.1, "mm")), linewidth = .35, colour = sf_cols$muted) +
    ggplot2::geom_vline(xintercept = 4.5, linetype = 2, linewidth = .35, colour = sf_cols$failure) +
    ggplot2::annotate("label", x = 3.35, y = 1.39, label = "Development-only:\nfeature selection, scaling,\nthresholds, calibration", fill = sf_cols$warning_light, colour = sf_cols$ink, linewidth = .18, size = 1.82, label.padding = grid::unit(.08, "lines"), lineheight = .86) +
    ggplot2::annotate("label", x = 5.10, y = 1.39, label = "No external\nretuning", fill = sf_cols$failure_light, colour = sf_cols$ink, linewidth = .18, size = 1.82, label.padding = grid::unit(.08, "lines"), lineheight = .86) +
    ggplot2::coord_cartesian(xlim = c(.45, 5.55), ylim = c(.64, 1.45), clip = "off") +
    sf_theme_void(7.7) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 3, 1, 3)) +
    sf_panel_tag("A")

  flow_counts <- data.frame(
    task = c("Primary\nATB vs OD", "Expanded\nATB vs nonATB"),
    task_ready = c(
      counts$n[counts$task == "primary_ATB_vs_OD" & counts$count_stage == "task_ready" & counts$cohort_family == "all"],
      counts$n[counts$task == "expanded_ATB_vs_nonATB" & counts$count_stage == "task_ready" & counts$cohort_family == "all"]
    ),
    external_task_ready = c(
      counts$n[counts$task == "primary_ATB_vs_OD" & counts$count_stage == "locked_external_task_ready" & counts$cohort_family == "all"],
      counts$n[counts$task == "expanded_ATB_vs_nonATB" & counts$count_stage == "locked_external_task_ready" & counts$cohort_family == "all"]
    ),
    final_evaluated = c(
      counts$n[counts$task == "primary_ATB_vs_OD" & counts$count_stage == "locked_external_evaluated" & counts$cohort_family == "all"],
      counts$n[counts$task == "expanded_ATB_vs_nonATB" & counts$count_stage == "locked_external_evaluated" & counts$cohort_family == "all"]
    ),
    colour_group = c("Primary", "Expanded"),
    stringsAsFactors = FALSE
  )
  count_long <- tidyr::pivot_longer(flow_counts, c("task_ready", "external_task_ready", "final_evaluated"), names_to = "stage", values_to = "n")
  count_long$stage <- factor(count_long$stage, levels = c("task_ready", "external_task_ready", "final_evaluated"), labels = c("Task-ready", "External\ntask-ready", "Locked\nevaluated"))
  count_long$stage_x <- as.numeric(count_long$stage)
  end_lab <- count_long[count_long$stage == "Locked\nevaluated", , drop = FALSE]
  end_lab$task_label <- ifelse(end_lab$colour_group == "Primary", "Primary\nATB vs OD", "Expanded\nATB vs nonATB")
  end_lab$label_y <- ifelse(end_lab$colour_group == "Primary", end_lab$n - 170, end_lab$n + 170)
  p1b <- ggplot2::ggplot(count_long, ggplot2::aes(x = stage_x, y = n, colour = colour_group, group = colour_group)) +
    ggplot2::annotate("rect", xmin = .93, xmax = 3.10, ymin = 0, ymax = max(count_long$n) * 1.06, fill = NA, colour = sf_cols$border, linewidth = .34) +
    ggplot2::geom_line(linewidth = .72, alpha = .92) +
    ggplot2::geom_point(ggplot2::aes(fill = colour_group), shape = 21, size = 3.35, colour = sf_cols$ink, stroke = .26) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.82, fontface = "bold", size = 2.20, colour = sf_cols$ink) +
    ggplot2::geom_label(data = end_lab, ggplot2::aes(x = 2.78, y = label_y, label = task_label, colour = colour_group), inherit.aes = FALSE, hjust = 1, size = 1.82, lineheight = .84, label.padding = grid::unit(.08, "lines"), linewidth = .14, fill = "white", show.legend = FALSE) +
    ggplot2::scale_x_continuous(breaks = 1:3, labels = levels(count_long$stage)) +
    ggplot2::scale_y_continuous(labels = scales::comma, breaks = c(0, 2000, 4000, 6000)) +
    ggplot2::scale_colour_manual(values = c("Primary" = sf_cols$primary, "Expanded" = sf_cols$expanded), guide = "none") +
    ggplot2::scale_fill_manual(values = c("Primary" = sf_cols$primary_light, "Expanded" = sf_cols$expanded_light), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(.88, 3.16), ylim = c(0, max(count_long$n) * 1.13), clip = "off") +
    ggplot2::labs(x = NULL, y = "Samples") +
    sf_theme(7.0) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 6.4, lineheight = .88),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 1)),
      legend.position = "none"
    ) +
    sf_panel_tag("B")

  mixed <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure2_mixed_family_availability.tsv"))
  locked_external <- mixed$n_mixed_families[mixed$task == "Primary" & mixed$split == "Locked external"][1]
  constraints <- data.frame(
    label = c("Primary mixed\ndevelopment", "Expanded mixed\ndevelopment", "Locked external\nfamilies"),
    n = c(
      mixed$n_mixed_families[mixed$task == "Primary" & mixed$split == "Development"][1],
      mixed$n_mixed_families[mixed$task == "Expanded" & mixed$split == "Development"][1],
      locked_external
    ),
    group = c("Primary", "Expanded", "External"),
    stringsAsFactors = FALSE
  )
  constraints$label <- factor(constraints$label, levels = constraints$label)
  p1c <- ggplot2::ggplot(constraints, ggplot2::aes(x = label, y = n, fill = group)) +
    ggplot2::geom_col(width = .58, colour = "white", linewidth = .35) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.34, fontface = "bold", size = 3.05, colour = sf_cols$ink) +
    ggplot2::scale_fill_manual(values = c("Primary" = sf_cols$primary_light, "Expanded" = sf_cols$expanded_light, "External" = sf_cols$signature_light), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, max(constraints$n) * 1.20), clip = "off") +
    ggplot2::labs(x = NULL, y = "Mixed cohort families") +
    sf_theme(7.0) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", lineheight = .88), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)), legend.position = "none") +
    sf_panel_tag("C")

  fam <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure2_family_phenotype_matrix.tsv"))
  fam_primary <- fam[fam$phenotype %in% c("ATB", "OD"), , drop = FALSE]
  fam_primary$phenotype <- factor(fam_primary$phenotype, levels = c("ATB", "OD"))
  fam_primary$split_label <- factor(fam_primary$split, levels = c("Development", "Locked external"))
  fam_primary$family_short <- sf_family_short(fam_primary$family_label)
  fam_primary$log_n <- log10(fam_primary$n + 1)
  primary_families <- unique(fam_primary[, c("split_label", "family_short", "family_n"), drop = FALSE])
  primary_grid <- merge(
    primary_families,
    data.frame(phenotype = levels(fam_primary$phenotype), stringsAsFactors = FALSE),
    all = TRUE
  )
  fam_plot <- merge(
    primary_grid,
    fam_primary[, c("split_label", "family_short", "phenotype", "n"), drop = FALSE],
    by = c("split_label", "family_short", "phenotype"),
    all.x = TRUE
  )
  fam_plot$n[is.na(fam_plot$n)] <- 0
  fam_plot$split_label <- factor(fam_plot$split_label, levels = c("Development", "Locked external"))
  fam_plot$phenotype <- factor(fam_plot$phenotype, levels = c("ATB", "OD"))
  fam_order <- primary_families$family_short[order(primary_families$split_label, -primary_families$family_n, primary_families$family_short)]
  fam_plot$family_short <- factor(fam_plot$family_short, levels = unique(fam_order))
  fam_plot$log_n <- log10(fam_plot$n + 1)
  fam_plot$label <- ifelse(fam_plot$n > 0, as.character(fam_plot$n), "")
  p1_matrix <- ggplot2::ggplot(fam_plot, ggplot2::aes(x = phenotype, y = family_short, fill = log_n)) +
    ggplot2::geom_tile(colour = "white", linewidth = .28) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 1.48, colour = sf_cols$ink) +
    ggplot2::facet_grid(split_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(low = "#F5F7F8", high = sf_cols$primary, name = "log10(n+1)") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "left", title.hjust = .5, barwidth = grid::unit(25, "mm"), barheight = grid::unit(2.3, "mm"))) +
    ggplot2::labs(x = NULL, y = NULL) +
    sf_theme(6.25) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 4.85),
      axis.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y = ggplot2::element_text(size = 6.0, angle = 90, face = "bold"),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = ggplot2::margin(1, 0, 0, 0),
      legend.title = ggplot2::element_text(size = 5.8, face = "bold"),
      legend.text = ggplot2::element_text(size = 5.5)
    ) +
    sf_panel_tag("D")

  oracle <- sf_sanitize(sf_read_tsv("results/ml_final/family_label_confounding_check.tsv"))
  oracle <- oracle[oracle$task == "primary_ATB_vs_OD", , drop = FALSE]
  oracle$split_label <- ifelse(oracle$split == "development", "Development", "Locked external")
  oracle$split_label <- factor(oracle$split_label, levels = c("Development", "Locked external"))
  final_primary <- sf_sanitize(sf_read_tsv("results/ml_final/family_cluster_bootstrap_ci.tsv"))
  final_primary <- final_primary[final_primary$task == "primary_ATB_vs_OD" & final_primary$selected_role == "best_performance", , drop = FALSE][1, ]
  p1d <- ggplot2::ggplot(oracle, ggplot2::aes(x = split_label, y = family_oracle_AUROC, fill = split_label)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2, colour = sf_cols$negative, linewidth = .3) +
    ggplot2::geom_hline(yintercept = final_primary$AUROC, linetype = 3, colour = sf_cols$primary, linewidth = .38) +
    ggplot2::geom_segment(ggplot2::aes(xend = split_label, y = .5, yend = family_oracle_AUROC), linewidth = .62, colour = sf_cols$muted) +
    ggplot2::geom_point(shape = 21, size = 3.85, colour = sf_cols$ink, stroke = .28) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", family_oracle_AUROC)), vjust = -0.85, fontface = "bold", size = 2.55, colour = sf_cols$ink) +
    ggplot2::annotate("label", x = 1.50, y = final_primary$AUROC + .034, label = paste0("Final model ", sf_num(final_primary$AUROC)), hjust = 0, size = 1.82, colour = sf_cols$primary, fill = "white", linewidth = .14, label.padding = grid::unit(.08, "lines")) +
    ggplot2::scale_fill_manual(values = c("Development" = sf_cols$failure_light, "Locked external" = sf_cols$signature_light), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0.45, 1.0), clip = "off") +
    ggplot2::labs(x = NULL, y = "Family-label oracle AUROC") +
    sf_theme(6.9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
      plot.margin = ggplot2::margin(5.5, 6, 5.5, 7)
    ) +
    sf_panel_tag("D")

  callout <- data.frame(
    item = c("Task-ready samples", "Mixed development families", "Development family-label oracle"),
    label = c("2538", "3", sprintf("%.3f", oracle$family_oracle_AUROC[oracle$split_label == "Development"][1])),
    interpretation = c("sample count alone looks adequate", "primary contrast support is thin", "family structure can classify labels"),
    stringsAsFactors = FALSE
  )

  sf_write_data(flow, "Figure1_pipeline.tsv")
  sf_write_data(count_long, "Figure1_sample_count_path.tsv")
  sf_write_data(constraints, "Figure1_design_constraints.tsv")
  sf_write_data(fam_plot, "Figure1_primary_phenotype_matrix.tsv")
  sf_write_data(oracle, "Figure1_family_label_oracle.tsv")
  sf_write_data(callout, "Figure1_primary_difficulty_callout.tsv")
  fig <- patchwork::wrap_plots(
    A = p1a,
    B = p1b,
    C = p1c,
    D = p1d,
    design = "
      AA
      BC
      BD
    ",
    widths = c(1.48, 1),
    heights = c(.58, 1, 1),
    guides = "keep"
  ) +
    sf_annotation()
  sf_save(fig, "figures_submission/main/Figure1_study_design_and_task_structure", 11.4, 6.8)
}
