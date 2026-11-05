sf_plot_figureS1 <- function() {
  fam <- sf_sanitize(sf_read_tsv("plot_inputs_submission/Figure2_family_phenotype_matrix.tsv"))
  fam$phenotype <- factor(fam$phenotype, levels = c("ATB", "OD", "LTBI", "HC"))
  fam$split <- factor(fam$split, levels = c("Development", "Locked external"))
  families <- unique(fam[, c("split", "cohort_family", "family_label", "family_n"), drop = FALSE])
  grid <- merge(families, data.frame(phenotype = levels(fam$phenotype), stringsAsFactors = FALSE), all = TRUE)
  full <- merge(grid, fam[, c("split", "cohort_family", "phenotype", "n"), drop = FALSE], by = c("split", "cohort_family", "phenotype"), all.x = TRUE)
  full$n[is.na(full$n)] <- 0
  full$phenotype <- factor(full$phenotype, levels = c("ATB", "OD", "LTBI", "HC"))
  full$split <- factor(full$split, levels = c("Development", "Locked external"))
  full$family_label <- factor(full$family_label, levels = unique(families$family_label[order(families$split, -families$family_n, families$family_label)]))
  full$log_n <- log10(full$n + 1)
  legend_counts <- c(0, 1, 10, 100, 500)
  p <- ggplot2::ggplot(full, ggplot2::aes(x = phenotype, y = family_label, fill = log_n)) +
    ggplot2::geom_tile(colour = "white", linewidth = .35) +
    ggplot2::facet_grid(split ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(
      low = "#FBFCFD",
      high = "#DABCC7",
      breaks = log10(legend_counts + 1),
      labels = legend_counts,
      name = "Sample count",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(3.5, "mm")
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    sf_theme(7.0) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 5.6),
      legend.position = "bottom",
      legend.title.align = .5
    )
  sf_write_data(full, "FigureS1_full_family_phenotype_matrix.tsv")
  sf_save(p, "figures_submission/supplement/FigureS1_full_family_phenotype_matrix", 8.8, 9.4)
}
