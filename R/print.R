#' @export
print.splcr_result <- function(x, digits = 4, ...) {
  stopifnot(is.list(x))
  method <- if (!is.null(x$method)) x$method else "spl"
  title <- switch(method,
                  spl_cr = "Smoothed Predictive Likelihood (SPL-CR)",
                  spl_fg = "Smoothed Predictive Likelihood (SPL-FG)",
                  "Smoothed Predictive Likelihood")
  cat(title, "\n", sep = "")
  if (!is.null(x$kernel) || !is.null(x$k_nn) || !is.null(x$folds)) {
    cat("Kernel: ", x$kernel,
        "  kNN: ", x$k_nn,
        "  folds: ", x$folds, "\n", sep = "")
  }
  if (!is.null(x$Lhat)) {
    cat("CV mean log-score (Lhat): ", format(round(x$Lhat, digits), nsmall = digits), "\n", sep = "")
  }
  if (!is.null(x$se_L)) {
    cat("Jackknife SE: ", format(round(x$se_L, digits), nsmall = digits), "\n", sep = "")
  }
  # Optional CI if fold_means available
  if (!is.null(x$fold_means)) {
    fm <- x$fold_means[is.finite(x$fold_means)]
    R <- length(fm)
    if (R >= 2L && !is.null(x$Lhat) && !is.null(x$se_L)) {
      crit <- stats::qt(0.975, df = R - 1L)
      lo <- x$Lhat - crit * x$se_L
      hi <- x$Lhat + crit * x$se_L
      cat("95% CI: [", format(round(lo, digits), nsmall = digits),
          ", ", format(round(hi, digits), nsmall = digits), "]\n", sep = "")
    }
  }
  invisible(x)
}
