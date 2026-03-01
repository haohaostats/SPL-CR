
#' Smoothed Predictive Likelihood for Fine--Gray (SPL-FG)
#'
#' Computes cross-validated smoothed predictive log-scores under a Fine--Gray
#' representation. This robust implementation uses \code{predictRisk} and explicit
#' data binding to safely bypass internal centering/scoping issues in cross-validation.
#'
#' Time-unit note:
#' \itemize{
#'   \item This function is unit-agnostic (days / years / months are all allowed).
#'   \item Absolute log-scores depend on the time unit because event terms involve a density.
#'   \item For fair model comparison, use the same time unit across all competing models.
#' }
#'
#' @inheritParams spl_csh
#' @param density_scale Multiplicative scale applied to event densities before taking logs.
#'   Use this only if you intentionally want to re-express densities on another time scale.
#'   For exact replication of manuscript PBC/Melanoma analyses (time in days), keep
#'   \code{density_scale = 1}.
#'
#' @return An object of class \code{splcr_result}.
#' @export
spl_fg <- function(formula,
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
                   verbose = FALSE,
                   density_scale = 1) {
  
  stopifnot(is.data.frame(data))
  if (!all(c(time, status) %in% names(data))) {
    stop("data must contain columns: ", time, ", ", status)
  }
  
  rhs <- .rhs_only(formula)
  vars <- all.vars(rhs)
  kid <- .kernel_id(kernel)
  
  n <- nrow(data)
  st_all <- as.integer(data[[status]])
  if (is.null(causes)) causes <- sort(unique(st_all[st_all != 0]))
  causes <- sort(as.integer(causes))
  if (length(causes) < 1L) stop("No event causes found.")
  
  if (is.null(fold_id)) {
    fold_id <- .assign_folds(status = st_all, n_fold = folds, seed = seed, stratify = stratify)
  } else {
    if (length(fold_id) != n) stop("fold_id must have length nrow(data).")
    fold_id <- as.integer(fold_id)
    folds <- max(fold_id)
  }
  
  if (!is.numeric(density_scale) || length(density_scale) != 1L ||
      !is.finite(density_scale) || density_scale <= 0) {
    stop("density_scale must be a positive finite numeric scalar.")
  }
  
  eps <- 1e-9
  
  # ---- helpers ---------------------------------------------------------------
  
  .norm_status <- function(df, causes_keep) {
    st <- df[[status]]
    if (is.factor(st)) {
      st <- suppressWarnings(as.integer(as.character(st)))
    } else {
      st <- as.integer(st)
    }
    st[is.na(st)] <- 0L
    keep <- c(0L, as.integer(causes_keep))
    st[!(st %in% keep)] <- 0L
    df[[status]] <- st
    df
  }
  
  .align_levels_row <- function(trn, one_row, vars_local) {
    out <- one_row
    for (v in vars_local) {
      if (!v %in% names(out) || !v %in% names(trn)) next
      xv <- out[[v]]
      xt <- trn[[v]]
      
      if (is.factor(xt)) {
        lev <- levels(xt)
        val <- if (is.factor(xv)) as.character(xv)[1] else as.character(xv)[1]
        if (!val %in% lev) val <- lev[1]
        out[[v]] <- factor(val, levels = lev)
      } else {
        if (!is.numeric(xv)) out[[v]] <- suppressWarnings(as.numeric(out[[v]]))
        if (!is.finite(out[[v]])) out[[v]] <- 0
      }
    }
    out
  }
  
  .make_zero_row <- function(trn, vars_local) {
    if (length(vars_local) == 0L) return(data.frame(`1` = 1)[, FALSE])
    out <- setNames(vector("list", length(vars_local)), vars_local)
    for (v in vars_local) {
      x <- trn[[v]]
      if (is.factor(x)) {
        lev <- levels(x)
        if (length(lev) == 0L) lev <- sort(unique(as.character(x)))
        out[[v]] <- factor(lev[1], levels = lev)
      } else {
        out[[v]] <- 0
      }
    }
    as.data.frame(out, stringsAsFactors = TRUE)
  }
  
  .fit_fgr_per_cause <- function(train_df, vars_local, causes_local) {
    keep_cols <- unique(c(time, status, vars_local))
    trn <- train_df[, intersect(colnames(train_df), keep_cols), drop = FALSE]
    trn <- .norm_status(trn, causes_local)
    
    fgr_list <- list()
    avail <- c()
    
    for (cc in causes_local) {
      if (sum(trn[[status]] == cc) == 0) {
        if (isTRUE(verbose)) message(sprintf("[spl_fg] skip cause %s (no events in training fold)", cc))
        next
      }
      
      rhs_str <- if (length(vars_local) == 0L) "1" else paste(vars_local, collapse = " + ")
      
      # IMPORTANT:
      # Use unqualified Hist(...) in the stored formula (manuscript-compatible),
      # while resolving it safely here via local binding.
      Hist <- prodlim::Hist
      form_fg <- stats::as.formula(
        paste0("Hist(", time, ", ", status, ", cens.code = 0) ~ ", rhs_str),
        env = environment()
      )
      
      fit <- tryCatch({
        tmp <- riskRegression::FGR(formula = form_fg, data = trn, cause = as.integer(cc))
        # Bind fold-local objects explicitly for predictRisk() stability in CV
        tmp$call$data <- trn
        tmp$call$formula <- form_fg
        tmp
      }, error = function(e) {
        if (isTRUE(verbose)) {
          message(sprintf("[spl_fg] FGR error for cause %s: %s", cc, conditionMessage(e)))
        }
        NULL
      })
      
      if (is.null(fit)) next
      fgr_list[[as.character(cc)]] <- fit
      avail <- c(avail, cc)
    }
    
    list(models = fgr_list, causes = sort(as.integer(avail)), trn = trn)
  }
  
  # Manuscript-aligned version:
  # Accept explicit times_grid so all causes can share the same all-cause event grid
  # in the training fold, matching the reference scripts.
  .baseline_cif_jumps_from_FGR <- function(fgr_fit, trn, vars_local, time_col, status_col, cc,
                                           times_grid = NULL) {
    if (is.null(times_grid)) {
      times_grid <- sort(unique(trn[[time_col]][trn[[status_col]] == cc]))
      if (length(times_grid) < 2L) {
        times_grid <- sort(unique(trn[[time_col]][trn[[status_col]] != 0]))
      }
      if (length(times_grid) < 2L) {
        tmax <- max(trn[[time_col]], na.rm = TRUE)
        times_grid <- seq(0.05 * tmax, 0.95 * tmax, length.out = 25)
      }
    }
    
    nd0 <- .make_zero_row(trn, vars_local)
    
    F0_vals <- tryCatch(
      as.numeric(riskRegression::predictRisk(fgr_fit, newdata = nd0, times = times_grid)),
      error = function(e) rep(NA_real_, length(times_grid))
    )
    
    F0_vals[!is.finite(F0_vals)] <- 0
    F0_vals <- pmin(pmax(F0_vals, 0), 1)
    dF <- diff(c(0, F0_vals))
    dF[dF < 0] <- 0
    
    cbind(time = times_grid, dF = dF)
  }
  
  .alpha_from_predictRisk <- function(fgr_fit, newdata_row, t0, F0s, eps_local = 1e-9) {
    F0s <- min(max(F0s, eps_local), 1 - eps_local)
    
    Fhat <- tryCatch(
      as.numeric(riskRegression::predictRisk(fgr_fit, newdata = newdata_row, times = t0))[1],
      error = function(e) NA_real_
    )
    if (!is.finite(Fhat)) return(NA_real_)
    
    Fhat <- min(max(Fhat, eps_local), 1 - eps_local)
    
    num <- log1p(-Fhat)
    den <- log1p(-F0s)
    alpha <- num / den
    
    if (!is.finite(alpha)) return(NA_real_)
    alpha <- min(max(alpha, exp(-6)), exp(6))
    alpha
  }
  
  fold_score_FG <- function(train_df, test_df) {
    causes_train <- sort(unique(as.integer(train_df[[status]][train_df[[status]] != 0])))
    if (length(causes_train) < 1L) return(rep(NA_real_, nrow(test_df)))
    
    vars_local <- vars
    fit_out <- .fit_fgr_per_cause(train_df, vars_local, causes_train)
    fgr_list <- fit_out$models
    causes_avail <- fit_out$causes
    trn <- fit_out$trn
    
    if (length(causes_avail) == 0L) return(rep(NA_real_, nrow(test_df)))
    
    # Manuscript-consistent: use all-cause training event times as common grid
    event_times_all <- sort(unique(trn[[time]][trn[[status]] != 0]))
    
    baseline_jumps <- list()
    for (cc in causes_avail) {
      tg <- if (length(event_times_all) >= 2L) event_times_all else NULL
      baseline_jumps[[as.character(cc)]] <- .baseline_cif_jumps_from_FGR(
        fgr_list[[as.character(cc)]],
        trn,
        vars_local,
        time,
        status,
        cc,
        times_grid = tg
      )
    }
    
    keep_cols <- unique(c(time, status, vars_local))
    tst <- test_df[, intersect(colnames(test_df), keep_cols), drop = FALSE]
    tst <- .norm_status(tst, causes_avail)
    
    score <- rep(NA_real_, nrow(tst))
    
    for (i in seq_len(nrow(tst))) {
      ti <- as.numeric(tst[[time]][i])
      di <- as.integer(tst[[status]][i])
      xi <- .align_levels_row(trn, tst[i, , drop = FALSE], vars_local)
      
      F_list <- list()
      f_list <- list()
      
      for (cc in causes_avail) {
        fit_c <- fgr_list[[as.character(cc)]]
        Jc <- baseline_jumps[[as.character(cc)]]
        
        # Smooth baseline level/density (F0s, f0s)
        if (is.null(Jc) || nrow(Jc) == 0L) {
          F0s <- eps
          f0s <- eps
        } else {
          pr <- tryCatch(
            .smooth_level_and_density_cpp(
              Jumps_H = as.matrix(Jc),
              nn = as.integer(k_nn),
              t = ti,
              kernel_id = as.integer(kid)
            ),
            error = function(e) c(0, 0)
          )
          F0s <- min(max(pr[1], eps), 1 - eps)
          f0s <- max(pr[2], eps)
        }
        
        # Reverse-engineer alpha from predictRisk; if invalid, fall back to baseline
        alpha <- .alpha_from_predictRisk(fit_c, xi, ti, F0s, eps_local = eps)
        
        if (!is.finite(alpha) || alpha <= 0) {
          F_list[[as.character(cc)]] <- F0s
          f_list[[as.character(cc)]] <- f0s
        } else {
          one_minus <- 1 - F0s
          F_ind <- 1 - (one_minus ^ alpha)
          f_ind <- alpha * (one_minus ^ (alpha - 1)) * f0s
          
          F_list[[as.character(cc)]] <- min(max(F_ind, eps), 1 - eps)
          f_list[[as.character(cc)]] <- max(f_ind, eps)
        }
      }
      
      # Overall survival S(t|x) = 1 - sum_c F_c(t|x)
      if (length(F_list) == 0L) {
        S_t <- 1.0
      } else {
        F_sum <- sum(vapply(F_list, function(z) ifelse(is.finite(z), z, 0), 0.0))
        F_sum <- min(max(F_sum, 0), 1 - eps)
        S_t <- max(1 - F_sum, eps)
      }
      
      # SPL-FG score
      if (di == 0L) {
        score[i] <- log(S_t)
      } else {
        fj <- f_list[[as.character(di)]]
        if (!length(fj) || !is.finite(fj)) fj <- eps
        score[i] <- log(max(fj * density_scale, eps))
      }
    }
    
    score
  }
  
  # ---- CV loop ---------------------------------------------------------------
  fold_means <- rep(NA_real_, folds)
  scores_all <- if (isTRUE(return_scores)) rep(NA_real_, n) else NULL
  
  for (f in seq_len(folds)) {
    test_idx <- which(fold_id == f)
    train_idx <- which(fold_id != f)
    
    scv <- tryCatch(
      fold_score_FG(train_df = data[train_idx, , drop = FALSE],
                    test_df  = data[test_idx,  , drop = FALSE]),
      error = function(e) {
        if (isTRUE(verbose)) {
          message(sprintf("[spl_fg] fold %d failed: %s", f, conditionMessage(e)))
        }
        rep(NA_real_, length(test_idx))
      }
    )
    
    fold_means[f] <- if (all(!is.finite(scv))) NA_real_ else mean(scv, na.rm = TRUE)
    if (isTRUE(return_scores)) scores_all[test_idx] <- scv
    
    if (isTRUE(verbose)) {
      na_n <- sum(!is.finite(scv))
      message(sprintf("[spl_fg] fold %d/%d: n_test=%d, NA=%d, mean=%.6f",
                      f, folds, length(test_idx), na_n, fold_means[f]))
    }
  }
  
  jk <- .jk_se_from_foldmeans(fold_means, n_total = n)
  
  out <- list(
    method = "spl_fg",
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