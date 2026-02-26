.kernel_id <- function(kernel) {
  if (is.numeric(kernel) && length(kernel) == 1L) return(as.integer(kernel))
  k <- tolower(as.character(kernel)[1])
  switch(k,
         uniform       = 1L,
         triangular    = 2L,
         epanechnikov  = 3L,
         quartic       = 4L,
         stop("Unknown kernel: ", kernel,
              ". Use one of: uniform, triangular, epanechnikov, quartic."))
}

.rhs_only <- function(formula) {
  if (!inherits(formula, "formula")) stop("formula must be a formula object.")
  if (length(formula) == 2L) return(formula)          # ~ x1 + x2
  if (length(formula) == 3L) return(stats::as.formula(paste0("~", deparse(formula[[3L]]))))
  stop("Unsupported formula format.")
}

.assign_folds <- function(status, n_fold = 10L, seed = NULL, stratify = TRUE) {
  n <- length(status)
  if (!is.null(seed)) set.seed(seed)

  if (!stratify) {
    idx <- sample.int(n)
    folds <- rep(seq_len(n_fold), length.out = n)
    out <- integer(n); out[idx] <- folds
    return(out)
  }

  st <- as.integer(status)
  levs <- sort(unique(st))
  out <- integer(n)
  for (lv in levs) {
    id <- which(st == lv)
    k <- length(id)
    if (k == 0) next
    id_shuf <- sample(id, k)
    folds_lv <- rep(seq_len(n_fold), length.out = k)
    out[id_shuf] <- folds_lv
  }
  out
}

.jk_se_from_foldmeans <- function(fold_means, n_total) {
  fold_means <- fold_means[is.finite(fold_means)]
  R <- length(fold_means)
  if (R < 2L) stop("Need at least 2 valid folds for jackknife SE.")
  Lhat <- mean(fold_means)
  sigma2_jk <- (R - 1) / R * sum((fold_means - Lhat)^2)
  se_L <- sqrt(sigma2_jk) / sqrt(n_total)
  list(Lhat = Lhat, sigma2_jk = sigma2_jk, se_L = se_L, R = R)
}

.cumhaz_at <- function(bh_time, bh_cumhaz, t) {
  if (!is.finite(t) || t <= 0) return(0)
  j <- findInterval(t, bh_time)
  if (j <= 0) return(0)
  bh_cumhaz[j]
}

.extract_cox_basehaz <- function(cox_fit) {
  bh <- survival::basehaz(cox_fit, centered = FALSE)
  if ("strata" %in% names(bh)) {
    stop("strata() terms are not supported in this minimal package version.")
  }
  bh <- bh[order(bh$time), , drop = FALSE]
  list(time = bh$time, cumhaz = bh$hazard)
}

.extract_cox_jumps <- function(cox_fit) {
  bh <- survival::basehaz(cox_fit, centered = FALSE)
  if ("strata" %in% names(bh)) {
    stop("strata() terms are not supported in this minimal package version.")
  }
  bh <- bh[order(bh$time), , drop = FALSE]
  dh <- c(bh$hazard[1], diff(bh$hazard))
  dh[dh < 0] <- 0
  cbind(time = bh$time, jump = dh)
}

.make_zero_row <- function(train_df, vars) {
  if (length(vars) == 0L) return(data.frame(`1` = 1)[, FALSE])
  out <- setNames(vector("list", length(vars)), vars)
  for (v in vars) {
    x <- train_df[[v]]
    if (is.factor(x)) {
      out[[v]] <- factor(levels(x)[1], levels = levels(x))
    } else {
      out[[v]] <- 0
    }
  }
  as.data.frame(out, stringsAsFactors = TRUE)
}

.align_test_to_train_levels <- function(train_df, test_df, vars) {
  out <- test_df
  for (v in vars) {
    if (!v %in% names(out) || !v %in% names(train_df)) next
    xt <- train_df[[v]]
    if (is.factor(xt)) {
      lev <- levels(xt)
      out[[v]] <- as.character(out[[v]])
      out[[v]][!(out[[v]] %in% lev)] <- lev[1]
      out[[v]] <- factor(out[[v]], levels = lev)
    } else {
      out[[v]] <- suppressWarnings(as.numeric(out[[v]]))
      out[[v]][!is.finite(out[[v]])] <- 0
    }
  }
  out
}
