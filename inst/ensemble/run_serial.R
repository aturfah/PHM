#!/usr/bin/env Rscript
## run_serial.R -- fits the whole ensemble in one R session, one replicate
## after another, with PHM::fitPHMEnsemble(). No parallelism and no
## per-replicate result files: only the final fit is written, to
## CFG$outPath.
##
## Usage: Rscript run_serial.R [config.R path]
##   config.R defaults to config.R in the current directory.
##
## The replicates match what parallel_run/ produces from the same config.R
## (same seeds), so the pooled components are identical. DeltaPmc and the merge
## heights can differ very slightly, because PHMv2() estimates DeltaPmc by Monte
## Carlo.

cliArgs <- commandArgs(trailingOnly = TRUE)
if (length(cliArgs) > 1) {
  stop("usage: Rscript run_serial.R [config.R path]")
}
configPath <- if (length(cliArgs) >= 1) cliArgs[1] else "config.R"

if (!file.exists(configPath)) {
  stop(sprintf("config file '%s' not found -- run this from the directory holding config.R, or pass its path: Rscript run_serial.R <config.R path>", configPath))
}
source(configPath)
suppressPackageStartupMessages(library(PHM))
if (!is.null(CFG$codeBase)) {
  message("note: CFG$codeBase is ignored; the fitting code comes from the installed PHM package")
}

X <- readRDS(CFG$dataPath)
cat(sprintf("[run_serial] fitting %d replicates one after another on a %d x %d matrix\n",
            CFG$M, nrow(X), ncol(X)))

fit <- PHM::fitPHMEnsemble(
  X, M = CFG$M, rho = CFG$rho, G = CFG$G, modelNames = CFG$modelNames,
  replicateWeighting = CFG$replicateWeighting, seedBase = CFG$seedBase,
  numCores = CFG$numCores, mcSamples = CFG$mcSamples, batchSize = CFG$batchSize,
  stage2Init = CFG$stage2Init, nRestartTries = CFG$nRestartTries,
  GStep = CFG$GStep, GMax = CFG$GMax, extendVerbose = CFG$extendVerbose,
  label = CFG$label
)

saveRDS(fit, CFG$outPath)
cat(sprintf("[run_serial] wrote %s (%d pooled components from %d replicates)\n",
            CFG$outPath, length(fit$paramsList), CFG$M))
