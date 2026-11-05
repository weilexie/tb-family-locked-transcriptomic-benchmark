sf_plot_figureS2 <- function() {
  mat <- sf_sanitize(sf_read_tsv("results/ml_final/published_signature_external_per_family.tsv"))
  mat <- mat[mat$task == "primary_ATB_vs_OD" & mat$evaluation_scope == "all_locked_external", , drop = FALSE]
  mat$signature_label <- sf_display_model(mat$signature_id)
  mat$family_label <- sf_family_short(mat$cohort_family)
  sig_rank <- stats::aggregate(AUROC ~ signature_label, mat, mean, na.rm = TRUE)
  sig_rank <- sig_rank[order(-sig_rank$AUROC), , drop = FALSE]
  keep_sig <- head(sig_rank$signature_label, 40)
  mat <- mat[mat$signature_label %in% keep_sig, , drop = FALSE]
  fam_rank <- stats::aggregate(AUROC ~ family_label, mat, mean, na.rm = TRUE)
  family_levels <- fam_rank$family_label[order(fam_rank$AUROC)]
  mat$family_label <- factor(stringr::str_wrap(mat$family_label, width = 12), levels = stringr::str_wrap(family_levels, width = 12))
  mat$signature_label <- factor(mat$signature_label, levels = rev(sig_rank$signature_label[sig_rank$signature_label %in% keep_sig]))
  p <- ggplot2::ggplot(mat, ggplot2::aes(x = family_label, y = signature_label, fill = AUROC)) +
    ggplot2::geom_tile(colour = "white", linewidth = .18) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FAFBFC", "#EEF3F1", "#DCE8E3", sf_cols$ensemble),
      limits = c(.5, 1),
      oob = scales::squish,
      na.value = "#EEF2F4",
      name = "AUROC",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        barwidth = grid::unit(42, "mm"),
        barheight = grid::unit(3.5, "mm")
      )
    ) +
    ggplot2::labs(x = "External family", y = NULL) +
    sf_theme(6.2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = .5, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 4.8),
      legend.position = "bottom",
      legend.title.align = .5
    )
  sf_write_data(mat, "FigureS2_full_reference_signature_matrix.tsv")
  sf_save(p, "figures_submission/supplement/FigureS2_full_reference_signature_matrix", 8.7, 10.4)
}
