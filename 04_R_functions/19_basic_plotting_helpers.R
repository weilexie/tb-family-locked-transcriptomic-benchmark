plot_pdf_png <- function(path_base, plot_fun, width = 8, height = 5) {
  pdf_path <- project_path(paste0(path_base, ".pdf"))
  png_path <- project_path(paste0(path_base, ".png"))
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(pdf_path, width = width, height = height)
  plot_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = width * 120, height = height * 120, res = 120)
  plot_fun()
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, png = png_path))
}

safe_barplot <- function(values, main, ylab = "", las = 2, col = "steelblue") {
  values <- values[!is.na(values)]
  clean_label <- function(x) {
    x <- as.character(x)
    x <- gsub("[\r\n\t]+", " ", x)
    x <- gsub("[[:cntrl:]]+", " ", x)
    trimws(x)
  }
  names(values) <- clean_label(names(values))
  main <- clean_label(main)
  ylab <- clean_label(ylab)
  if (length(values) == 0) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No data available")
  } else {
    graphics::barplot(values, main = main, ylab = ylab, las = las, col = col, cex.names = 0.75)
  }
}
