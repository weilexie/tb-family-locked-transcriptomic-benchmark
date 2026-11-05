sf_plot_figureS3 <- function() {
  neg <- sf_sanitize(sf_read_tsv("manuscript/tables/Supplementary_Table_negative_control_summary.tsv"))
  neg <- unique(neg)
  neg$task <- "summary"
  neg$selected_role <- "summary"
  neg$control_family <- "global_negative_control"
  neg$display_label <- dplyr::case_when(
    neg$negative_control == "random_uniform_score" ~ "Random uniform scores",
    neg$negative_control == "selected_external_score_label_permutation" ~ "Locked external label permutation",
    neg$negative_control == "random_reference_signature_benchmark_spread" ~ "Random reference-signature benchmark spread",
    TRUE ~ gsub("_", " ", neg$negative_control)
  )

  wf <- sf_sanitize(sf_read_tsv("results/ml_final/negative_control_family_stratified_permutation.tsv"))
  wf$control_family <- "within_family_score_permutation"
  wf$display_label <- paste0(
    "Within-family permutation: ",
    ifelse(wf$task == "primary_ATB_vs_OD", "primary task", "expanded task"),
    ", ",
    ifelse(wf$selected_role == "best_performance", "best-performance model", "sparse model")
  )
  common <- c("task", "selected_role", "negative_control", "control_family", "n_iterations", "AUROC_median", "AUROC_q025", "AUROC_q975", "interpretation", "display_label")
  neg <- neg[, common, drop = FALSE]
  wf <- wf[, common, drop = FALSE]
  neg_all <- rbind(neg, wf)
  neg_all <- unique(neg_all)
  neg_all$control_type <- ifelse(
    neg_all$control_family == "within_family_score_permutation",
    "Within-family score permutation",
    "Global negative controls"
  )
  neg_all$control_type <- factor(neg_all$control_type, levels = c("Global negative controls", "Within-family score permutation"))
  neg_all$negative_control_label <- stringr::str_wrap(neg_all$display_label, 28)
  neg_all <- neg_all[order(neg_all$AUROC_median), ]
  neg_all$negative_control_label <- factor(neg_all$negative_control_label, levels = rev(neg_all$negative_control_label))
  neg_all$display_value <- sprintf("%.3f", neg_all$AUROC_median)
  p <- ggplot2::ggplot(neg_all, ggplot2::aes(x = AUROC_median, y = negative_control_label)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, linewidth = .25, colour = sf_cols$muted) +
    ggplot2::geom_segment(ggplot2::aes(x = AUROC_q025, xend = AUROC_q975, yend = negative_control_label), linewidth = .58, colour = sf_cols$ink) +
    ggplot2::geom_point(ggplot2::aes(fill = control_type), shape = 21, size = 2.85, colour = sf_cols$ink, stroke = .24) +
    ggplot2::geom_label(ggplot2::aes(x = 0.875, label = display_value), hjust = 1, size = 2.05, colour = sf_cols$ink, fill = "white", linewidth = 0, label.padding = grid::unit(.05, "lines")) +
    ggplot2::scale_fill_manual(
      values = c("Global negative controls" = sf_cols$signature_light, "Within-family score permutation" = sf_cols$warning_light),
      name = NULL
    ) +
    ggplot2::coord_cartesian(xlim = c(0.35, 0.88), clip = "off") +
    ggplot2::labs(x = "Median AUROC across permutations (95% interval)", y = NULL) +
    sf_theme(7.8) +
    ggplot2::theme(legend.position = "bottom")
  sf_write_data(neg_all, "FigureS3_permutation_negative_controls.tsv")
  sf_save(p, "figures_submission/supplement/FigureS3_permutation_negative_controls", 7.2, 4.6)
}
