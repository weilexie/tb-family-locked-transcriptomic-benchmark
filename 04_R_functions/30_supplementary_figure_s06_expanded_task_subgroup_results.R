sf_plot_figureS6 <- function() {
  pred <- sf_sanitize(sf_read_tsv("results/ml_final/external_error_by_sample.tsv"))
  pred <- pred[pred$task == "expanded_ATB_vs_nonATB" & pred$selected_role == "best_sparse", , drop = FALSE]
  pred$score <- pred$predicted_probability
  vars <- c("geographical_region", "tissue", "platform_type", "hiv_status", "age_group")
  rows <- list()
  for (v in vars) {
    lv <- sort(unique(pred[[v]]))
    for (x in lv) {
      d <- pred[pred[[v]] == x, , drop = FALSE]
      if (nrow(d) < 30 || length(unique(d$true_label)) < 2) next
      m <- compute_binary_metrics(d$true_label, d$score, threshold = 0.5)
      ci <- sf_bootstrap_auc_ci(d$true_label, d$score, n_boot = 200, seed = 20260511 + length(rows))
      rows[[length(rows) + 1]] <- data.frame(
        subgroup_variable = v,
        subgroup_level = x,
        subgroup_label = paste0(v, ": ", x),
        n = m$n,
        n_pos = m$n_pos,
        AUROC = m$AUROC,
        AUROC_ci_low = ci[1],
        AUROC_ci_high = ci[2],
        AUPRC = m$AUPRC,
        sensitivity = m$sensitivity,
        specificity = m$specificity,
        stringsAsFactors = FALSE
      )
    }
  }
  sg <- do.call(rbind, rows)
  sg$subgroup_variable <- factor(sg$subgroup_variable, levels = vars)
  sg$subgroup_group <- dplyr::recode(
    as.character(sg$subgroup_variable),
    geographical_region = "Geographic region",
    tissue = "Tissue",
    platform_type = "Platform",
    hiv_status = "HIV status",
    age_group = "Age group"
  )
  pretty_level <- function(x) {
    y <- gsub("_", " ", as.character(x))
    y <- tools::toTitleCase(tolower(y))
    y <- gsub("Hiv", "HIV", y, fixed = TRUE)
    y <- gsub("Rna Seq", "RNA-seq", y, fixed = TRUE)
    y <- gsub("Rna-Seq", "RNA-seq", y, fixed = TRUE)
    y <- gsub("Us", "US", y, fixed = TRUE)
    y <- gsub("Uk", "UK", y, fixed = TRUE)
    y
  }
  sg$plot_label <- paste0(sg$subgroup_group, ": ", pretty_level(sg$subgroup_level))
  label_order <- unique(sg$plot_label[order(match(as.character(sg$subgroup_variable), vars), -sg$n, sg$plot_label)])
  sg$plot_label <- factor(sg$plot_label, levels = rev(label_order))
  sg$y <- match(as.character(sg$plot_label), levels(sg$plot_label))
  auroc_xlim <- c(0.35, 1.02)
  n_axis_max <- ceiling(max(sg$n, na.rm = TRUE) / 100) * 100
  n_to_x <- function(n) auroc_xlim[1] + (as.numeric(n) / n_axis_max) * diff(auroc_xlim)
  x_to_n <- function(x) (as.numeric(x) - auroc_xlim[1]) / diff(auroc_xlim) * n_axis_max
  sg$n_x <- n_to_x(sg$n)
  top_breaks_n <- c(0, 300, 600, 900, n_axis_max)
  top_breaks_n <- unique(top_breaks_n[top_breaks_n <= n_axis_max])
  top_breaks_x <- n_to_x(top_breaks_n)
  p <- ggplot2::ggplot(sg, ggplot2::aes(x = AUROC, y = y)) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = auroc_xlim[1], xmax = n_x, ymin = y - .34, ymax = y + .34),
      fill = sf_cols$expanded_light,
      alpha = .55,
      colour = sf_cols$border,
      linewidth = .18
    ) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = sf_cols$muted, linewidth = .25) +
    ggplot2::geom_segment(ggplot2::aes(x = AUROC_ci_low, xend = AUROC_ci_high, yend = y), colour = sf_cols$ink, linewidth = .48) +
    ggplot2::geom_point(shape = 21, size = 2.35, fill = sf_cols$expanded_light, colour = sf_cols$ink, stroke = .22) +
    ggplot2::scale_x_continuous(
      limits = auroc_xlim,
      breaks = c(0.5, 0.7, 0.9),
      sec.axis = ggplot2::sec_axis(~ x_to_n(.), name = "Sample size (n)", breaks = top_breaks_n)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(levels(sg$plot_label)),
      labels = levels(sg$plot_label),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(x = "External AUROC (95% bootstrap interval)", y = NULL) +
    sf_theme(6.8) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5.5),
      axis.title.x.top = ggplot2::element_text(size = 6.1, colour = sf_cols$ink),
      axis.text.x.top = ggplot2::element_text(size = 5.0, colour = sf_cols$muted),
      axis.ticks.x.top = ggplot2::element_line(colour = sf_cols$border),
      panel.grid.major.x = ggplot2::element_line(colour = sf_cols$grid, linewidth = .25),
      panel.grid.major.y = ggplot2::element_line(colour = "#F3F6F8", linewidth = .25),
      plot.margin = ggplot2::margin(7, 18, 7, 18)
    )
  sf_write_data(sg, "FigureS6_expanded_task_subgroup_results.tsv")
  sf_save(p, "figures_submission/supplement/FigureS6_expanded_task_subgroup_results", 8.8, 7.6)
}
