spl_csh <- function(formula, data, time = "time", status = "status",
                    causes = NULL, folds = 10L, k_nn = 5L, kernel = "uniform",
                    seed = NULL, stratify = TRUE, fold_id = NULL,
                    return_scores = FALSE, verbose = FALSE) {
  kid <- .kernel_id(kernel)
  if (length(kid) != 1L || is.na(kid) || !kid %in% 1:4) stop("Invalid kernel.")
  kernel <- c("uniform", "triangular", "epanechnikov", "quartic")[[kid]]
  if (!is.data.frame(data) || !all(c(time, status) %in% names(data))) stop("Invalid data or time/status columns.")
  tt <- data[[time]]; st <- data[[status]]; n <- nrow(data)
  if (!is.numeric(tt) || any(!is.finite(tt)) || any(tt <= 0)) stop("Follow-up times must be finite and positive.")
  if (!is.numeric(st) || any(!is.finite(st)) || any(st < 0 | st != floor(st))) stop("Status must contain nonnegative integers.")
  if (is.null(causes)) causes <- sort(unique(st[st > 0]))
  if (!is.numeric(causes) || !length(causes) || any(!is.finite(causes)) ||
      any(causes <= 0 | causes != floor(causes)) || anyDuplicated(causes) ||
      !setequal(causes, unique(st[st > 0]))) stop("Include every observed cause exactly once.")
  causes <- sort(as.integer(causes))
  if (length(k_nn) != 1L || !is.finite(k_nn) || k_nn < 1 || k_nn != floor(k_nn)) stop("Invalid k_nn.")
  rhs <- .rhs_only(formula)
  if (!all(all.vars(rhs) %in% names(data))) stop("Formula variables are missing from data.")
  if (anyNA(data[, unique(c(time, status, all.vars(rhs))), drop = FALSE])) stop("Missing model data.")
  if (is.null(fold_id)) {
    if (length(folds) != 1L || !is.finite(folds) || folds < 2 || folds > n || folds != floor(folds)) stop("Invalid folds.")
    fold_id <- .assign_folds(st, folds, seed, stratify)
  }
  if (!is.numeric(fold_id) || length(fold_id) != n || any(!is.finite(fold_id)) ||
      any(fold_id < 1 | fold_id != floor(fold_id))) stop("Invalid fold_id.")
  folds <- max(fold_id)
  if (folds < 2L || !identical(sort(unique(as.integer(fold_id))), seq_len(folds))) stop("Folds must be consecutive and nonempty.")
  fold_id <- as.integer(fold_id)
  scores <- numeric(n); fold_means <- numeric(folds); fold_sizes <- integer(folds)
  truncated <- integer(folds)
  fit_warnings <- character()
  for (r in seq_len(folds)) {
    test_idx <- which(fold_id == r)
    train <- data[fold_id != r, , drop = FALSE]
    test <- data[test_idx, , drop = FALSE]
    tau <- max(train[[time]])
    outside <- test[[time]] > tau
    test[[time]] <- pmin(test[[time]], tau)
    test[[status]][outside] <- 0L
    truncated[r] <- sum(outside)
    sc <- numeric(nrow(test))
    for (cc in causes) {
      if (!any(train[[status]] == cc)) stop("Training fold ", r, " has no events for cause ", cc, ".")
      fit_data <- train
      event_name <- ".splcr_event"
      while (event_name %in% names(fit_data)) event_name <- paste0(event_name, "_")
      fit_data[[event_name]] <- as.integer(train[[status]] == cc)
      lhs <- substitute(survival::Surv(TIME, EVENT), list(TIME = as.name(time), EVENT = as.name(event_name)))
      form <- stats::as.formula(call("~", lhs, rhs[[2L]]), env = environment(rhs))
      fit <- withCallingHandlers(
        survival::coxph(form, data = fit_data, x = TRUE, y = TRUE, na.action = stats::na.fail),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, paste0("Fold ", r, ", cause ", cc, ": ", conditionMessage(w)))
        })
      if (any(!is.finite(stats::coef(fit)))) stop("Nonfinite Cox coefficient in fold ", r, ".")
      bh <- survival::basehaz(fit, centered = FALSE)
      if ("strata" %in% names(bh)) stop("Stratified Cox baselines are not supported.")
      increments <- c(bh$hazard[1L], diff(bh$hazard))
      if (any(!is.finite(increments)) || any(increments < -1e-10)) stop("Invalid cumulative baseline hazard.")
      baseline <- .csh_lc_baseline(.csh_baseline(bh$time, pmax(increments, 0), tau,
                                                k = k_nn, kernel = kernel, tail = "error"))
      risk <- exp(stats::predict(fit, newdata = test, type = "lp", reference = "zero"))
      if (any(!is.finite(risk)) || any(risk <= 0)) stop("Invalid Cox risk multiplier.")
      contribution <- -.csh_lc_cumulative(baseline, test[[time]]) * risk
      ev <- test[[status]] == cc
      contribution[ev] <- contribution[ev] + log(.csh_lc_hazard(baseline, test[[time]][ev]) * risk[ev])
      sc <- sc + contribution
    }
    if (any(!is.finite(sc))) stop("Nonfinite predictive score in fold ", r, ".")
    scores[test_idx] <- sc
    fold_means[r] <- mean(sc)
    fold_sizes[r] <- length(sc)
    if (isTRUE(verbose)) message(sprintf("[spl_csh] fold %d/%d: n_test=%d, mean=%.6f", r, folds, length(sc), fold_means[r]))
  }
  jk <- .jk_se_from_foldmeans(fold_means, n)
  out <- list(method = "spl_csh", kernel = kernel, k_nn = as.integer(k_nn),
              folds = as.integer(folds), fold_means = fold_means, Lhat = mean(scores),
              sigma2_jk = jk$sigma2_jk, se_L = jk$se_L, fold_id = fold_id,
              algorithm = "local_scaled", fold_sizes = fold_sizes,
              truncated = truncated, fit_warnings = unique(fit_warnings))
  if (isTRUE(return_scores)) out$scores <- scores
  class(out) <- c("splcr_result", "list")
  out
}
