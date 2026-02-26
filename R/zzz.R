.onLoad <- function(libname, pkgname) {
  # Ensure Rcpp is loaded so its C-callables used by BEGIN_RCPP are available
  requireNamespace("Rcpp", quietly = TRUE)
  invisible(NULL)
}
