spl_compare <- function(models,
                        data,
                        method = c("csh", "fg"),
                        kernels = NULL,
                        ...) {

  method <- match.arg(method)
  if (is.null(kernels)) kernels <- c("uniform", "triangular", "epanechnikov", "quartic")
  if (!is.list(models) || length(models) == 0L) stop("models must be a non-empty list of formulas.")
  if (is.null(names(models))) stop("models must be a named list.")

  dots <- list(...)

  time <- if (!is.null(dots$time)) dots$time else "time"
  status <- if (!is.null(dots$status)) dots$status else "status"
  folds <- if (!is.null(dots$folds)) dots$folds else 10L
  stratify <- if (!is.null(dots$stratify)) dots$stratify else TRUE
  seed <- if (!is.null(dots$seed)) dots$seed else NULL

  st <- as.integer(data[[status]])
  fold_id <- dots$fold_id
  dots$fold_id <- NULL
  if (is.null(fold_id)) fold_id <- .assign_folds(status = st, n_fold = folds, seed = seed, stratify = stratify)

  res <- list()
  ii <- 1L

  for (mname in names(models)) {
    form <- models[[mname]]
    for (k in kernels) {
      if (method == "csh") {
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
