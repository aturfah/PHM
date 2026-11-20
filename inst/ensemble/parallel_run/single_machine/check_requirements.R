#!/usr/bin/env Rscript
## check_requirements.R -- preflight checks, run by launch.sh before it starts
## any replicate: the shared scripts exist, the PHM package (and therefore
## every package it imports) loads, and CFG$dataPath exists and is readable.
## Catching a problem once, up front, beats finding it M times over in M
## crashed background replicates.
##
## Usage: Rscript check_requirements.R <config.R path> [shared script dir]
## [shared script dir] is parallel_run/, where run_replicate.R and
## run_combine.R live; launch.sh passes it. Without it those two are not
## checked.
##
## Prints one line and exits 0 if everything checks out; otherwise lists
## every problem found and exits 1. It does not open the data file, so a
## wrongly shaped matrix shows up only when the first replicate runs.

cliArgs <- commandArgs(trailingOnly = TRUE)
if (length(cliArgs) < 1 || length(cliArgs) > 2) {
  stop("usage: Rscript check_requirements.R <config.R path> [shared script dir]")
}
configPath <- cliArgs[1]
sharedDir <- if (length(cliArgs) >= 2) cliArgs[2] else NA

if (!file.exists(configPath)) {
  stop(sprintf("config file '%s' not found", configPath))
}
source(configPath)

problems <- character(0)

if (!is.na(sharedDir)) {
  for (f in c("run_replicate.R", "run_combine.R")) {
    p <- file.path(sharedDir, f)
    if (!file.exists(p)) problems <- c(problems, sprintf("missing script: %s", p))
  }
}

## Loading PHM's namespace also loads every package it imports, so this one
## check covers mclust, mvtnorm and the rest.
phmOk <- tryCatch({ loadNamespace("PHM"); TRUE }, error = function(e) {
  problems <<- c(problems, sprintf(
    "the PHM package could not be loaded (%s) -- install it with remotes::install_github(\"aturfah/PHM\"), or R CMD INSTALL from a local checkout",
    conditionMessage(e)))
  FALSE
})
if (phmOk && !exists("fitDensityReplicate", envir = asNamespace("PHM"))) {
  problems <- c(problems, sprintf(
    "the installed PHM (version %s) has no fitDensityReplicate() -- install a version that includes the ensemble functions",
    as.character(utils::packageVersion("PHM"))))
}

if (is.null(CFG$dataPath) || !file.exists(CFG$dataPath)) {
  problems <- c(problems, sprintf(
    "CFG$dataPath ('%s') does not exist -- set it in %s to the saved data matrix (an N x D numeric matrix saved with saveRDS())",
    if (is.null(CFG$dataPath)) "NULL" else CFG$dataPath, configPath))
} else if (file.access(CFG$dataPath, mode = 4) != 0) {
  problems <- c(problems, sprintf("CFG$dataPath ('%s') exists but isn't readable", CFG$dataPath))
}

if (length(problems) > 0) {
  cat("preflight checks failed:\n", file = stderr())
  for (p in problems) cat(sprintf("  - %s\n", p), file = stderr())
  quit(status = 1, save = "no")
}

cat(sprintf("preflight checks passed (PHM %s, dataPath OK)\n", as.character(utils::packageVersion("PHM"))))
