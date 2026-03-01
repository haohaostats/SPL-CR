#' Smoothed Predictive Likelihood for Cause-Specific Hazards (SPL-CSH)
#'
#' Computes cross-validated smoothed predictive log-scores under the cause-specific
#' hazards representation by smoothing baseline cumulative hazard jumps on each
#' training fold using a k-nearest-neighbour kernel smoother.
#'
#' This function is implemented to be numerically consistent with the reference
#' analysis scripts used in the accompanying manuscript.
#'
#' Time-unit note:
#' \itemize{
#'   \item This function is unit-agnostic (days / years / months are all allowed).
#'   \item Absolute log-scores depend on the time unit because event terms involve a hazard density.
#'   \item For fair model comparison, use the same time unit across all competing models.
#' }
#'
#' @param formula A model formula. Recommended format is a right-hand side only formula
#'   such as \code{~ x1 + x2}. If a full formula is supplied, only the right-hand side
#'   is used.
#' @param data A data.frame containing at least \code{time} and \code{status}.
#' @param time Name of the follow-up time column.
#' @param status Name of the event indicator column (0=censoring; positive integers are causes).
#' @param causes Integer vector of causes to include. Default uses all positive values in \code{status}.
#' @param folds Number of CV folds.
#' @param k_nn Number of nearest neighbours for kNN bandwidth.
#' @param kernel Kernel name or id (uniform, triangular, epanechnikov, quartic).
#' @param seed Optional integer seed to generate folds. If NULL, no seed is set.
#' @param stratify Whether to stratify fold assignment by \code{status}.
#' @param fold_id Optional precomputed fold assignment (integer vector length nrow(data)).
#'   If provided, \code{seed} and \code{stratify} are ignored.
#' @param return_scores Whether to return per-observation out-of-fold scores.
#' @param verbose Whether to print fold progress.
#'
#' @return An object of class \code{splcr_result}.
#' @export
spl_csh <- function(formula,
                   data,
                   time = "time",
                   status = "status",
                   causes = NULL,
                   folds = 10L,
                   k_nn = 5L,
                   kernel = "epanechnikov",
                   seed = NULL,
                   stratify = TRUE,
                   fold_id = NULL,
                   return_scores = FALSE,
                   verbose = FALSE) {
  
  stopifnot(is.data.frame(data))
  if (!all(c(time, status) %in% names(data))) stop("data must contain columns: ", time, ", ", status)
  
  rhs <- .rhs_only(formula)
  kid <- .kernel_id(kernel)
  
  n <- nrow(data)
  st <- as.integer(data[[status]])
  
  if (is.null(causes)) causes <- sort(unique(st[st != 0]))
  causes <- sort(as.integer(causes))
  if (length(causes) == 0L) stop("No event causes found.")
  
  if (is.null(fold_id)) {
    fold_id <- .assign_folds(status = st, n_fold = folds, seed = seed, stratify = stratify)
  } else {
    if (length(fold_id) != n) stop("fold_id must have length nrow(data).")
    fold_id <- as.integer(fold_id)
    folds <- max(fold_id)
  }
  
  fold_means <- rep(NA_real_, folds)
  scores_all <- if (isTRUE(return_scores)) rep(NA_real_, n) else NULL
  eps <- 1e-9
  
  # helper consistent with the manuscript scripts:
  # H(t|x) = H0_hat(t) * risk, where H0_hat comes from survfit(coxph)
  .CH_cox_safe_fold <- function(sf_obj, cox_fit, subject_row, t0) {
    # sf_obj is survfit(cox_fit) computed on the training fold
    if (is.null(sf_obj) || is.null(sf_obj$n.event) || sum(sf_obj$n.event) == 0) return(0)
    if (!is.finite(t0) || t0 <= 0) return(0)
    
    risk <- as.numeric(stats::predict(cox_fit, newdata = subject_row, type = "risk"))
    if (!is.finite(risk)) return(0)
    
    H0 <- sf_obj$cumhaz[sf_obj$time <= t0]
    if (length(H0) == 0) return(0)
    as.numeric(tail(H0, 1)) * risk
  }
  
  for (r in seq_len(folds)) {
    test_idx <- which(fold_id == r)
    train_idx <- which(fold_id != r)
    
    d_tr <- data[train_idx, , drop = FALSE]
    d_te <- data[test_idx,  , drop = FALSE]
    
    # Fit Cox per cause and store:
    #  - cox model
    #  - survfit object (for cumulative hazard used in total_H_t)
    #  - baseline jumps from basehaz(centered=FALSE) for smoothing
    fits <- vector("list", length(causes)); names(fits) <- as.character(causes)
    sf_list <- vector("list", length(causes)); names(sf_list) <- as.character(causes)
    jumps_list <- vector("list", length(causes)); names(jumps_list) <- as.character(causes)
    
    for (cc in causes) {
      d_tr$status_k <- as.integer(as.integer(d_tr[[status]]) == cc)
      form_cc <- stats::as.formula(paste0("survival::Surv(", time, ", status_k) ", deparse(rhs)))
      fit_cc <- survival::coxph(form_cc, data = d_tr, x = TRUE)
      
      fits[[as.character(cc)]] <- fit_cc
      sf_list[[as.character(cc)]] <- survival::survfit(fit_cc)
      jumps_list[[as.character(cc)]] <- .extract_cox_jumps(fit_cc)  # basehaz(centered=FALSE) diff
    }
    
    # Score test fold
    sc <- numeric(length(test_idx))
    
    for (j in seq_along(test_idx)) {
      subject <- d_te[j, , drop = FALSE]
      t0 <- as.numeric(subject[[time]])
      d0 <- as.integer(subject[[status]])
      
      # total cumulative hazard across causes at t0, using survfit-based cumhaz (paper code)
      total_H_t <- 0
      for (cc in causes) {
        total_H_t <- total_H_t + .CH_cox_safe_fold(sf_list[[as.character(cc)]],
                                                   fits[[as.character(cc)]],
                                                   subject, t0)
      }
      
      if (d0 == 0L) {
        sc[j] <- -total_H_t
      } else {
        if (!(d0 %in% causes)) {
          sc[j] <- log(eps) - total_H_t
        } else {
          model_j <- fits[[as.character(d0)]]
          jumps_j <- jumps_list[[as.character(d0)]]
          
          smooth_h0_t <- .calculate_smoothed_hazard_cpp(
            Jumps_H0 = as.matrix(jumps_j),
            bw = -1,
            nn = as.integer(k_nn),
            t = t0,
            method_id = as.integer(kid)
          )
          
          # scale by exp(lp) as in manuscript scripts
          lp_j <- as.numeric(stats::predict(model_j, newdata = subject, type = "lp"))
          risk_score_j <- exp(lp_j)
          smooth_h_t <- smooth_h0_t * risk_score_j
          
          log_h_term <- log(smooth_h_t + eps)
          sc[j] <- log_h_term - total_H_t
        }
      }
    }
    
    fold_means[r] <- mean(sc, na.rm = TRUE)
    if (isTRUE(return_scores)) scores_all[test_idx] <- sc
    
    if (isTRUE(verbose)) {
      message(sprintf("[spl_csh] fold %d/%d: n_test=%d, mean=%.6f",
                      r, folds, length(test_idx), fold_means[r]))
    }
  }
  
  jk <- .jk_se_from_foldmeans(fold_means, n_total = n)
  
  out <- list(
    method = "spl_csh",
    kernel = kernel,
    k_nn = as.integer(k_nn),
    folds = as.integer(folds),
    fold_means = fold_means,
    Lhat = jk$Lhat,
    sigma2_jk = jk$sigma2_jk,
    se_L = jk$se_L,
    fold_id = fold_id
  )
  if (isTRUE(return_scores)) out$scores <- scores_all
  class(out) <- c("splcr_result", "list")
  out
}