spl_compare_wide <- function(x,
                             value_col = "Lhat",
                             digits = 4L,
                             kernel_order = c("Uniform", "Triangular", "Epanechnikov", "Quartic")) {
  if (!is.data.frame(x)) stop("x must be a data.frame (typically the output of spl_compare()).")

  req <- c("model", "kernel", value_col)
  miss <- setdiff(req, names(x))
  if (length(miss) > 0L) {
    stop("Missing required columns in x: ", paste(miss, collapse = ", "))
  }

  dat <- x[, c("model", "kernel", value_col), drop = FALSE]

  k_raw <- tolower(as.character(dat$kernel))
  k_map <- c(
    "uniform"      = "Uniform",
    "triangular"   = "Triangular",
    "epanechnikov" = "Epanechnikov",
    "quartic"      = "Quartic"
  )
  dat$kernel <- ifelse(k_raw %in% names(k_map), unname(k_map[k_raw]), as.character(dat$kernel))

  model_levels <- unique(as.character(dat$model))
  dat$model <- factor(as.character(dat$model), levels = model_levels)

  dat <- dat[!duplicated(dat[, c("model", "kernel")]), , drop = FALSE]

  wide <- stats::reshape(
    dat,
    idvar = "model",
    timevar = "kernel",
    direction = "wide"
  )

  wide <- wide[order(wide$model), , drop = FALSE]

  prefix <- paste0(value_col, ".")
  nm <- names(wide)
  names(wide) <- sub(paste0("^", gsub("\\.", "\\\\.", prefix)), "", nm)
  names(wide)[names(wide) == "model"] <- "Model"

  desired <- c("Model", kernel_order)
  keep <- intersect(desired, names(wide))
  extra <- setdiff(names(wide), keep)
  wide <- wide[, c(keep, extra), drop = FALSE]

  num_cols <- setdiff(names(wide), "Model")
  for (cc in num_cols) {
    if (is.numeric(wide[[cc]])) wide[[cc]] <- round(wide[[cc]], digits = digits)
  }

  rownames(wide) <- NULL

  attr(wide, "value_col") <- value_col
  attr(wide, "digits") <- as.integer(digits)
  class(wide) <- c("splcr_compare_wide", "data.frame")

  wide
}

.pad_right <- function(x, width) {
  x <- as.character(x)
  paste0(x, strrep(" ", pmax(0L, width - nchar(x, type = "width"))))
}

.pad_left <- function(x, width) {
  x <- as.character(x)
  paste0(strrep(" ", pmax(0L, width - nchar(x, type = "width"))), x)
}

print.splcr_compare_wide <- function(x,
                                     digits = attr(x, "digits", exact = TRUE) %||% 4L,
                                     mark_best = TRUE,
                                     title = NULL,
                                     ...) {
  stopifnot(is.data.frame(x))

  if (is.null(title)) {
    metric <- attr(x, "value_col", exact = TRUE)
    metric <- if (is.null(metric)) "Lhat" else metric
    title <- sprintf("\U0001f4ca SPL model comparison  \u2022  metric = %s", metric)
  }

  dat <- x
  cols <- names(dat)

  num_cols <- setdiff(cols, "Model")
  num_cols <- num_cols[vapply(dat[num_cols], is.numeric, logical(1))]

  disp <- dat

  for (cc in num_cols) {
    disp[[cc]] <- formatC(disp[[cc]], format = "f", digits = digits)
  }

  if (isTRUE(mark_best) && length(num_cols) > 0L && nrow(dat) > 0L) {
    for (cc in num_cols) {
      vals <- dat[[cc]]
      if (all(!is.finite(vals))) next
      idx <- which(vals == max(vals, na.rm = TRUE))

      disp[[cc]][idx] <- paste0("\u2605 ", disp[[cc]][idx])
    }
  }

  widths <- integer(length(cols))
  for (j in seq_along(cols)) {
    col_txt <- c(cols[j], as.character(disp[[j]]))
    widths[j] <- max(nchar(col_txt, type = "width"), na.rm = TRUE)
  }

  cat(title, "\n", sep = "")
  cat(strrep("\u2500", max(sum(widths) + 3L * (length(widths) - 1L), nchar(title, type = "width"))), "\n", sep = "")

  header <- character(length(cols))
  for (j in seq_along(cols)) {
    if (j == 1L) {
      header[j] <- .pad_right(cols[j], widths[j])
    } else {
      header[j] <- .pad_left(cols[j], widths[j])
    }
  }
  cat(paste(header, collapse = "   "), "\n", sep = "")
  cat(strrep("\u2500", sum(widths) + 3L * (length(widths) - 1L)), "\n", sep = "")

  for (i in seq_len(nrow(disp))) {
    row_out <- character(length(cols))
    for (j in seq_along(cols)) {
      txt <- as.character(disp[i, j, drop = TRUE])
      if (j == 1L) {
        row_out[j] <- .pad_right(txt, widths[j])
      } else {
        row_out[j] <- .pad_left(txt, widths[j])
      }
    }
    cat(paste(row_out, collapse = "   "), "\n", sep = "")
  }

  cat(strrep("\u2500", sum(widths) + 3L * (length(widths) - 1L)), "\n", sep = "")
  if (isTRUE(mark_best) && length(num_cols) > 0L) {
    cat("Note: \u2605 marks the highest value in each kernel column (larger is better).\n")
  }

  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
