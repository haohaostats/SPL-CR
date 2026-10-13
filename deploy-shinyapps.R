if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Install 'rsconnect' first: install.packages('rsconnect')")
}

rsconnect::deployApp(
  appDir = ".",
  appFiles = "app.R",
  appPrimaryDoc = "app.R",
  appName = "splcr-web",
  account = "haostats"
)
