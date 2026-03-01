#' Compare multiple model specifications using SPL scores
#'
#' Convenience wrapper to evaluate a list of formulas under SPL-CR or SPL-FG
#' across one or more kernels using a shared fold assignment.
#'
#' @param models Named list of formulas.
#' @param data A data.frame with time and status columns.
#' @param method Either \code{"csh"} or \code{"fg"}.
#' @param kernels Character vector of kernels.
#' @param ... Passed to \code{spl_csh()} or \code{spl_fg()}.
#'
#' @return A data.frame with columns: model, kernel, Lhat, se_L.
#' @export
spl_compare <- function(models,
                        data,
                        method = c("cr", "fg"),
                        kernels = c("uniform","triangular","epanechnikov","quartic"),
                        ...) {

  method <- match.arg(method)
  if (!is.list(models) || length(models) == 0L) stop("models must be a non-empty list of formulas.")
  if (is.null(names(models))) stop("models must be a named list.")

  dots <- list(...)

  time <- if (!is.null(dots$time)) dots$time else "time"
  status <- if (!is.null(dots$status)) dots$status else "status"
  folds <- if (!is.null(dots$folds)) dots$folds else 10L
  stratify <- if (!is.null(dots$stratify)) dots$stratify else TRUE
  seed <- if (!is.null(dots$seed)) dots$seed else NULL

  st <- as.integer(data[[status]])
  fold_id <- .assign_folds(status = st, n_fold = folds, seed = seed, stratify = stratify)

  res <- list()
  ii <- 1L

  for (mname in names(models)) {
    form <- models[[mname]]
    for (k in kernels) {
      if (method == "cr") {
        out <- do.call(spl_csh, c(list(formula = form, data = data, kernel = k, fold_id = fold_id), dots))
      } else {
        out <- do.call(spl_fg, c(list(formula = form, data = data, kernel = k, fold_id = fold_id), dots))
      }
      res[[ii]] <- data.frame(
        model = mname,
        kernel = k,
        Lhat = out$Lhat,
        se_L = out$se_L,
        stringsAsFactors = FALSE
      )
      ii <- ii + 1L
    }
  }

  do.call(rbind, res)
}
