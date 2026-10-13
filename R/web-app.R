run_splcr_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the web app. Please install it first.")
  }

  app <- shiny::shinyApp(ui = .splcr_app_ui(), server = .splcr_app_server())
  shiny::runApp(app, ...)
}

.splcr_get_fn <- function(name) {
  env <- environment()
  while (!identical(env, emptyenv())) {
    if (exists(name, envir = env, inherits = FALSE)) {
      return(get(name, envir = env, inherits = FALSE))
    }
    env <- parent.env(env)
  }
  if (exists(name, envir = .GlobalEnv, inherits = TRUE)) {
    return(get(name, envir = .GlobalEnv, inherits = TRUE))
  }
  if ("package:splcr" %in% search() && exists(name, envir = as.environment("package:splcr"), inherits = FALSE)) {
    return(get(name, envir = as.environment("package:splcr"), inherits = FALSE))
  }
  if ("splcr" %in% loadedNamespaces()) {
    ns <- asNamespace("splcr")
    if (exists(name, envir = ns, inherits = FALSE)) {
      return(get(name, envir = ns, inherits = FALSE))
    }
  }
  stop(sprintf("Required function '%s' is not available in the current app environment.", name))
}

.splcr_demo_pbc_data <- function() {
  data("pbc", package = "survival", envir = environment())
  pbc_data <- subset(pbc, !is.na(trt))
  pbc_data$status <- ifelse(
    pbc_data$status == 2, 1L,
    ifelse(pbc_data$status == 1, 2L, 0L)
  )
  pbc_data$time <- as.numeric(pbc_data$time)
  pbc_data$logbili <- log(pbc_data$bili)
  pbc_data$sex <- factor(pbc_data$sex)
  pbc_data$stage <- factor(pbc_data$stage)
  pbc_data$ascites <- factor(pbc_data$ascites)

  keep <- c("time", "status", "age", "sex", "stage", "logbili", "albumin", "protime", "ascites")
  pbc_data[stats::complete.cases(pbc_data[, keep]), keep, drop = FALSE]
}

.splcr_prepare_data <- function(dat) {
  stopifnot(is.data.frame(dat))
  dat[] <- lapply(dat, function(x) {
    if (is.character(x) || is.logical(x)) factor(x) else x
  })
  dat
}

.splcr_numeric_choices <- function(dat) {
  names(dat)[vapply(dat, is.numeric, logical(1))]
}

.splcr_status_choices <- function(dat) {
  .splcr_numeric_choices(dat)
}

.splcr_predictor_choices <- function(dat, time_col = NULL, status_col = NULL) {
  setdiff(names(dat), c(time_col, status_col))
}

.splcr_formula_from_text <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) stop("Formula text cannot be empty.")
  if (!grepl("^~", txt)) txt <- paste("~", txt)
  stats::as.formula(txt)
}

.splcr_models_from_text <- function(txt) {
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE)))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) stop("Provide at least one model definition.")

  models <- vector("list", length(lines))
  names(models) <- character(length(lines))

  for (i in seq_along(lines)) {
    line <- lines[[i]]
    parts <- strsplit(line, "\\s*(=|:)\\s*", perl = TRUE)[[1]]
    if (length(parts) < 2L) {
      stop("Each model line must use 'name = predictors' or 'name: predictors'.")
    }
    name <- trimws(parts[1])
    rhs <- trimws(paste(parts[-1], collapse = " = "))
    if (!nzchar(name) || !nzchar(rhs)) stop("Model names and formulas must be non-empty.")
    names(models)[i] <- name
    models[[i]] <- .splcr_formula_from_text(rhs)
  }

  models
}

.splcr_round_df <- function(dat, digits = 4L) {
  out <- dat
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  out
}

.splcr_status_label_map <- function(values, source = "upload") {
  vals <- sort(unique(as.integer(values)))
  labels <- vapply(vals, function(v) {
    if (identical(source, "demo")) {
      return(switch(as.character(v), `0` = "Censored", `1` = "Death", `2` = "Transplant", paste("Cause", v)))
    }
    if (v == 0L) return("Censored")
    paste("Cause", v)
  }, character(1))
  stats::setNames(labels, as.character(vals))
}

.splcr_status_table <- function(values, source = "upload") {
  map <- .splcr_status_label_map(values, source = source)
  data.frame(
    Code = as.integer(names(map)),
    Label = unname(map),
    check.names = FALSE
  )
}
.splcr_parse_status_labels <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) return(setNames(character(0), character(0)))
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE)))
  lines <- lines[nzchar(lines)]
  out <- character(0)
  for (line in lines) {
    parts <- strsplit(line, "\\s*=\\s*", perl = TRUE)[[1]]
    if (length(parts) < 2L) next
    key <- trimws(parts[1])
    val <- trimws(paste(parts[-1], collapse = " = "))
    if (!grepl("^-?[0-9]+$", key) || !nzchar(val)) next
    out[key] <- val
  }
  out
}

.splcr_status_table_from_map <- function(map) {
  data.frame(
    Code = as.integer(names(map)),
    Label = unname(map),
    check.names = FALSE
  )
}

.splcr_result_summary <- function(res) {
  ci_low <- NA_real_
  ci_high <- NA_real_
  fm <- res$fold_means[is.finite(res$fold_means)]
  if (length(fm) >= 2L && is.finite(res$Lhat) && is.finite(res$se_L)) {
    crit <- stats::qt(0.975, df = length(fm) - 1L)
    ci_low <- res$Lhat - crit * res$se_L
    ci_high <- res$Lhat + crit * res$se_L
  }

  .splcr_round_df(data.frame(
    Method = res$method,
    Kernel = res$kernel,
    kNN = res$k_nn,
    Folds = res$folds,
    Lhat = res$Lhat,
    SE = res$se_L,
    CI_low = ci_low,
    CI_high = ci_high,
    check.names = FALSE
  ))
}

.splcr_table_html <- function(dat, digits = 4L, highlight_max = FALSE) {
  dat <- .splcr_round_df(dat, digits = digits)
  cell_class <- matrix("", nrow = nrow(dat), ncol = ncol(dat))

  if (highlight_max && nrow(dat) > 0L) {
    num_cols <- which(vapply(dat, is.numeric, logical(1)))
    for (j in num_cols) {
      vals <- dat[[j]]
      if (all(!is.finite(vals))) next
      idx <- which(vals == max(vals, na.rm = TRUE))
      cell_class[idx, j] <- "is-best"
    }
  }

  header <- shiny::tags$tr(lapply(names(dat), function(nm) shiny::tags$th(nm)))
  rows <- lapply(seq_len(nrow(dat)), function(i) {
    shiny::tags$tr(
      lapply(seq_len(ncol(dat)), function(j) {
        value <- dat[i, j, drop = TRUE]
        txt <- if (is.numeric(value)) formatC(value, digits = digits, format = "f") else as.character(value)
        shiny::tags$td(class = cell_class[i, j], txt)
      })
    )
  })

  shiny::tags$table(
    class = "splcr-table",
    shiny::tags$thead(header),
    shiny::tags$tbody(rows)
  )
}

.splcr_metric_card <- function(title, value) {
  shiny::div(
    class = "metric-card",
    shiny::tags$div(class = "metric-label", title),
    shiny::tags$div(class = "metric-value", value)
  )
}

.splcr_app_ui <- function() {
  shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        body { background: linear-gradient(180deg, #f6f1e8 0%, #ffffff 52%, #eef5f2 100%); color: #1f2a2a; }
        .app-shell { max-width: 1440px; margin: 0 auto; }
        .hero { margin: 18px 0 14px; padding: 24px 28px; border-radius: 20px;
                background: linear-gradient(135deg, #113d3c 0%, #285b57 55%, #dca95b 100%);
                color: #f7f1e8; box-shadow: 0 18px 48px rgba(17, 61, 60, 0.18); }
        .hero h1 { margin: 0 0 8px; font-size: 30px; font-weight: 700; letter-spacing: 0.3px; }
        .hero p { margin: 0; max-width: 900px; font-size: 15px; line-height: 1.6; opacity: 0.95; }
        .panel-card { background: rgba(255,255,255,0.94); border: 1px solid rgba(17, 61, 60, 0.08);
                      border-radius: 18px; padding: 18px; box-shadow: 0 10px 30px rgba(17, 61, 60, 0.08); margin-bottom: 16px; }
        .section-title { margin: 0 0 10px; font-size: 16px; font-weight: 700; color: #113d3c; }
        .small-note { color: #516362; font-size: 12px; line-height: 1.5; }
        .shiny-input-container { width: 100% !important; }
        .btn-primary { background: #113d3c; border-color: #113d3c; }
        .btn-primary:hover { background: #0d3130; border-color: #0d3130; }
        .nav-tabs { border-bottom: 1px solid rgba(17, 61, 60, 0.1); }
        .nav-tabs > li > a { color: #113d3c; }
        .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover {
          border: none; border-bottom: 3px solid #d08d38; color: #113d3c; background: transparent; font-weight: 700;
        }
        .metric-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 10px; margin-bottom: 14px; }
        .metric-card { background: linear-gradient(180deg, #f5f8f7 0%, #edf3f1 100%); border: 1px solid #d9e6e1; border-radius: 14px; padding: 12px 14px; }
        .metric-label { font-size: 11px; text-transform: uppercase; letter-spacing: 0.08em; color: #5e7271; margin-bottom: 6px; }
        .metric-value { font-size: 21px; font-weight: 700; color: #113d3c; }
        .plot-card { background: #faf8f4; border: 1px solid #e8dfd3; border-radius: 16px; padding: 10px; }
        .plot-caption { color: #5f6f6d; font-size: 12px; margin-top: 8px; }
        .table-wrap { overflow-x: auto; }
        .splcr-table { width: 100%; border-collapse: separate; border-spacing: 0; font-size: 13px; }
        .splcr-table thead th { background: #123f3d; color: #f7f1e8; padding: 10px 12px; text-align: left; font-weight: 700; }
        .splcr-table thead th:first-child { border-top-left-radius: 12px; }
        .splcr-table thead th:last-child { border-top-right-radius: 12px; }
        .splcr-table tbody td { padding: 10px 12px; border-bottom: 1px solid #ece5da; background: #fffdfa; }
        .splcr-table tbody tr:nth-child(even) td { background: #f8f5ee; }
        .splcr-table tbody td.is-best { font-weight: 700; color: #0e5751; background: #e7f5ef; }
        .table-note { color: #5e7271; font-size: 12px; margin-top: 8px; }
        pre { background: #f7f5f1; border: 1px solid #e6ddd0; border-radius: 12px; }
      "))
    ),
    shiny::div(
      class = "app-shell",
      shiny::div(
        class = "hero",
        shiny::h1("SPL-CR Web Companion"),
        shiny::p("Interactive validation for competing-risks prediction models. Review data quality first, inspect basic event patterns, then run SPL-CSH, SPL-FG, or multi-model comparison from the browser.")
      ),
      shiny::fluidRow(
        shiny::column(
          width = 4,
          shiny::div(
            class = "panel-card",
            shiny::h3(class = "section-title", "Data Overview"),
            shiny::radioButtons(
              "data_source",
              "Data source",
              choices = c("Built-in PBC demo" = "demo", "Upload CSV" = "upload")
            ),
            shiny::conditionalPanel(
              "input.data_source === 'upload'",
              shiny::fileInput("csv_file", "CSV file", accept = c(".csv", "text/csv")),
              shiny::checkboxInput("strings_factor", "Convert text columns to factor", TRUE),
              shiny::textAreaInput(
                "status_labels_text",
                "Custom status labels",
                rows = 4,
                value = "0 = Censored\n1 = Death\n2 = Transplant",
                placeholder = "Example:\n0 = Censored\n1 = Death\n2 = Transplant"
              ),
              shiny::actionButton("generate_status_template", "Generate template", class = "btn-primary"),
              shiny::p(class = "small-note", "Define the English label for each numeric status code using one line per code, or generate a template from the current status column.")
            ),
            shiny::uiOutput("data_metrics")
          ),
          shiny::div(
            class = "panel-card",
            shiny::h3(class = "section-title", "Analysis Setup"),
            shiny::selectInput("time_col", "Time column", choices = NULL),
            shiny::selectInput("status_col", "Status column", choices = NULL),
            shiny::selectizeInput("predictors", "Predictors", choices = NULL, multiple = TRUE),
            shiny::actionButton("use_predictors", "Use selected predictors", class = "btn-primary"),
            shiny::hr(),
            shiny::selectInput("method", "Method", choices = c("SPL-CSH" = "csh", "SPL-FG" = "fg")),
            shiny::selectInput("kernel", "Kernel", choices = c("uniform", "triangular", "epanechnikov", "quartic"), selected = "uniform"),
            shiny::checkboxGroupInput("kernels_compare", "Compare kernels", choices = c("uniform", "triangular", "epanechnikov", "quartic"), selected = c("uniform", "triangular", "epanechnikov", "quartic")),
            shiny::numericInput("folds", "CV folds", value = 10, min = 2, max = 50, step = 1),
            shiny::numericInput("k_nn", "k-NN smoothing", value = 5, min = 1, max = 50, step = 1),
            shiny::numericInput("seed", "Seed", value = 20250101, min = 1, step = 1),
            shiny::checkboxInput("stratify", "Stratify folds by status", TRUE),
            shiny::hr(),
            shiny::p(class = "small-note", "Status coding must follow the package convention: 0 = censoring, positive integers = event causes. All result tables are displayed to four decimal places.")
          )
        ),
        shiny::column(
          width = 8,
          shiny::div(
            class = "panel-card",
            shiny::tabsetPanel(
              id = "main_tabs",
              shiny::tabPanel(
                "Data Preview",
                shiny::br(),
                shiny::h4(class = "section-title", "Status Labels"),
                shiny::uiOutput("status_map_ui"),
                shiny::br(),
                shiny::h4(class = "section-title", "Preview"),
                shiny::uiOutput("data_head_ui")
              ),
              shiny::tabPanel(
                "Single Model",
                shiny::br(),
                shiny::textInput("formula_text", "Formula", value = "~ logbili + albumin + protime + ascites"),
                shiny::actionButton("run_single", "Run analysis", class = "btn-primary"),
                shiny::br(),
                shiny::br(),
                shiny::uiOutput("single_summary_ui"),
                shiny::verbatimTextOutput("single_print")
              ),
              shiny::tabPanel(
                "Model Comparison",
                shiny::br(),
                shiny::textAreaInput(
                  "models_text",
                  "Models (one per line)",
                  rows = 9,
                  value = paste(
                    "M1: age + sex + stage",
                    "M2: logbili + albumin",
                    "M3: logbili + albumin + protime + ascites",
                    "M4: age + sex + stage + logbili + albumin",
                    "M5: age + sex + stage + logbili + albumin + protime + ascites",
                    sep = "\n"
                  )
                ),
                shiny::actionButton("run_compare", "Compare models", class = "btn-primary"),
                shiny::br(),
                shiny::br(),
                shiny::h4(class = "section-title", "Long Format Results"),
                shiny::uiOutput("compare_table_ui"),
                shiny::br(),
                shiny::h4(class = "section-title", "Wide Comparison Table"),
                shiny::uiOutput("compare_wide_ui"),
                shiny::div(class = "table-note", "Highest value in each kernel column is highlighted. All numeric results are shown with four decimal places.")
              )
            )
          )
        )
      )
    )
  )
}

.splcr_app_server <- function() {
  function(input, output, session) {
    shiny::observeEvent(input$method, {
      kernels <- c("uniform", "triangular", "epanechnikov", "quartic")
      selected <- if (identical(input$method, "csh")) "uniform" else "epanechnikov"
      shiny::updateSelectInput(session, "kernel", choices = kernels, selected = selected)
      shiny::updateCheckboxGroupInput(session, "kernels_compare", choices = kernels, selected = kernels)
    }, ignoreInit = FALSE)
    dataset <- shiny::reactive({
      if (identical(input$data_source, "demo")) return(.splcr_demo_pbc_data())
      shiny::req(input$csv_file)
      dat <- utils::read.csv(input$csv_file$datapath, stringsAsFactors = isTRUE(input$strings_factor), check.names = FALSE)
      .splcr_prepare_data(dat)
    })

    status_map <- shiny::reactive({
      dat <- dataset()
      shiny::req(input$status_col)
      base_map <- .splcr_status_label_map(dat[[input$status_col]], source = input$data_source)
      if (!identical(input$data_source, "upload")) return(base_map)
      custom_map <- .splcr_parse_status_labels(input$status_labels_text)
      if (!length(custom_map)) return(base_map)
      base_map[names(custom_map)] <- unname(custom_map)
      base_map
    })

    shiny::observeEvent(dataset(), {
      dat <- dataset()
      time_choices <- .splcr_numeric_choices(dat)
      status_choices <- .splcr_status_choices(dat)
      default_time <- if ("time" %in% names(dat)) "time" else if (length(time_choices)) time_choices[[1]] else names(dat)[[1]]
      default_status <- if ("status" %in% names(dat)) "status" else if (length(status_choices)) status_choices[[1]] else names(dat)[[1]]
      predictor_choices <- .splcr_predictor_choices(dat, default_time, default_status)
      default_predictors <- intersect(c("logbili", "albumin", "protime", "ascites"), predictor_choices)

      shiny::updateSelectInput(session, "time_col", choices = time_choices, selected = default_time)
      shiny::updateSelectInput(session, "status_col", choices = status_choices, selected = default_status)
      shiny::updateSelectizeInput(session, "predictors", choices = predictor_choices, selected = default_predictors, server = TRUE)
    }, ignoreInit = FALSE)

    shiny::observeEvent(list(input$time_col, input$status_col, dataset()), {
      dat <- dataset()
      predictor_choices <- .splcr_predictor_choices(dat, input$time_col, input$status_col)
      current <- intersect(input$predictors, predictor_choices)
      shiny::updateSelectizeInput(session, "predictors", choices = predictor_choices, selected = current, server = TRUE)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$use_predictors, {
      txt <- if (length(input$predictors) == 0L) "~ 1" else paste("~", paste(input$predictors, collapse = " + "))
      shiny::updateTextInput(session, "formula_text", value = txt)
    })

    shiny::observeEvent(input$generate_status_template, {
      shiny::req(identical(input$data_source, "upload"), input$status_col)
      dat <- dataset()
      shiny::req(input$status_col %in% names(dat))
      vals <- sort(unique(as.integer(dat[[input$status_col]])))
      map <- .splcr_status_label_map(vals, source = "upload")
      template <- paste(sprintf("%s = %s", names(map), unname(map)), collapse = "\n")
      shiny::updateTextAreaInput(session, "status_labels_text", value = template)
    })

    analysis_args <- shiny::reactive({
      shiny::req(dataset(), input$time_col, input$status_col)
      dat <- dataset()
      validate_msg <- NULL
      if (!is.numeric(dat[[input$time_col]])) validate_msg <- "The selected time column must be numeric."
      if (!is.numeric(dat[[input$status_col]])) validate_msg <- "The selected status column must be numeric and use 0, 1, 2, ... coding."
      shiny::validate(shiny::need(is.null(validate_msg), validate_msg))

      list(
        data = dat,
        time = input$time_col,
        status = input$status_col,
        folds = as.integer(input$folds),
        k_nn = as.integer(input$k_nn),
        seed = as.integer(input$seed),
        stratify = isTRUE(input$stratify)
      )
    })

    single_result <- shiny::eventReactive(input$run_single, {
      args <- analysis_args()
      form <- .splcr_formula_from_text(input$formula_text)
      shiny::withProgress(message = "Running SPL analysis", value = 0.2, {
        res <- if (identical(input$method, "csh")) {
          do.call(.splcr_get_fn("spl_csh"), c(list(formula = form, kernel = input$kernel), args))
        } else {
          do.call(.splcr_get_fn("spl_fg"), c(list(formula = form, kernel = input$kernel), args))
        }
        shiny::incProgress(0.8)
        res
      })
    })

    compare_result <- shiny::eventReactive(input$run_compare, {
      shiny::req(length(input$kernels_compare) > 0L)
      args <- analysis_args()
      models <- .splcr_models_from_text(input$models_text)
      shiny::withProgress(message = "Comparing model specifications", value = 0.2, {
        res <- do.call(
          .splcr_get_fn("spl_compare"),
          c(list(models = models, method = input$method, kernels = input$kernels_compare), args)
        )
        shiny::incProgress(0.8)
        res
      })
    })

    output$data_metrics <- shiny::renderUI({
      dat <- dataset()
      status_vals <- if (!is.null(input$status_col) && input$status_col %in% names(dat)) dat[[input$status_col]] else 0
      shiny::div(
        class = "metric-grid",
        .splcr_metric_card("Rows", format(nrow(dat), big.mark = ",")),
        .splcr_metric_card("Columns", ncol(dat)),
        .splcr_metric_card("Event Causes", length(setdiff(sort(unique(status_vals)), 0)))
      )
    })

    output$status_map_ui <- shiny::renderUI({
      shiny::req(input$status_col)
      tbl <- .splcr_status_table_from_map(status_map())
      shiny::div(class = "table-wrap", .splcr_table_html(tbl, digits = 4L))
    })

    output$data_head_ui <- shiny::renderUI({
      dat <- utils::head(dataset(), 8)
      shiny::div(class = "table-wrap", .splcr_table_html(dat, digits = 4L))
    })

    output$single_summary_ui <- shiny::renderUI({
      res <- single_result()
      shiny::tagList(
        if (length(res$fit_warnings)) shiny::div(class = "alert alert-warning", "Cox fitting warnings were recorded; inspect model convergence before interpreting the score."),
        shiny::div(class = "table-wrap", .splcr_table_html(.splcr_result_summary(res), digits = 4L))
      )
    })

    output$single_print <- shiny::renderPrint({
      print(single_result())
    })

    output$compare_table_ui <- shiny::renderUI({
      shiny::div(class = "table-wrap", .splcr_table_html(compare_result(), digits = 4L))
    })

    output$compare_wide_ui <- shiny::renderUI({
      wide <- .splcr_get_fn("spl_compare_wide")(compare_result())
      shiny::div(class = "table-wrap", .splcr_table_html(wide, digits = 4L, highlight_max = TRUE))
    })
  }
}
