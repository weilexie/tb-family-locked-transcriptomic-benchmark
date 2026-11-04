calibration_summary <- function(y_true, prob) {
  y <- as_binary(y_true)
  eps <- 1e-6
  p <- pmin(pmax(prob, eps), 1 - eps)
  fit <- try(stats::glm(y ~ qlogis(p), family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  co <- coef(fit)
  data.frame(calibration_intercept = unname(co[1]), calibration_slope = unname(co[2]))
}

