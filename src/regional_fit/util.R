read_csv <- function(path, ...) {
  read.csv(path, ..., stringsAsFactors = FALSE, check.names = FALSE)
}


write_png <- function(filename, code, ...) {
  png(filename, ...)
  on.exit(dev.off())
  force(code)
}


version_check <- function(package, version) {
  if (packageVersion(package) < version) {
    stop(sprintf(
      paste("Please update %s with: drat:::add('ncov-ic');",
            "install.packages('%s')"),
      package, package))
  }
}


get_n_threads <- function() {
  n <- as.integer(Sys.getenv("CONTEXT_CORES", Sys.getenv("MC_CORES", 1)))
  message(sprintf("Running on %d threads", n))
  n
}
