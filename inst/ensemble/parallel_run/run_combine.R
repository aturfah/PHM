#!/usr/bin/env Rscript
## run_combine.R -- the reduce step: reads the M replicate results written by
## run_replicate.R, pools them with PHM::combineDensityReplicates() (which
## runs PHMv2()), and saves the fit to CFG$outPath. Used by both run styles;
## run it once, after every replicate has finished.
##
## Needs only the M small result files, never the data, so its cost does not
## depend on the data size.
##
## Usage: Rscript run_combine.R [config.R path]
##   config.R defaults to config.R in the current directory.

cliArgs <- commandArgs(trailingOnly = TRUE)
if (length(cliArgs) > 1) {
  stop("usage: Rscript run_combine.R [config.R path]")
}
configPath <- if (length(cliArgs) >= 1) cliArgs[1] else "config.R"

if (!file.exists(configPath)) {
  stop(sprintf("config file '%s' not found -- run this from the directory holding config.R, or pass its path: Rscript run_combine.R <config.R path>", configPath))
}
source(configPath)
suppressPackageStartupMessages(library(PHM))
if (!is.null(CFG$codeBase)) {
  message("note: CFG$codeBase is ignored; the fitting code comes from the installed PHM package")
}

resultFiles <- file.path(CFG$resultsDir, sprintf("replicate_%03d.rds", seq_len(CFG$M)))
missing <- resultFiles[!file.exists(resultFiles)]
if (length(missing) > 0) {
  stop(sprintf(
    "missing %d of %d replicate result file(s) (e.g. %s) -- did every replicate finish? Check the replicate logs, rerun the failed ones, then run this again",
    length(missing), CFG$M, missing[1]))
}

replicateResults <- lapply(resultFiles, readRDS)

fit <- PHM::combineDensityReplicates(
  replicateResults, replicateWeighting = CFG$replicateWeighting,
  numCores = CFG$numCores, mcSamples = CFG$mcSamples, batchSize = CFG$batchSize,
  label = CFG$label
)

saveRDS(fit, CFG$outPath)
cat(sprintf("[run_combine] wrote %s (%d pooled components from %d replicates)\n",
            CFG$outPath, length(fit$paramsList), CFG$M))
