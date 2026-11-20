#!/usr/bin/env Rscript
## run_replicate.R -- fits ONE ensemble replicate (the "map" step) with
## PHM::fitDensityReplicate() and writes resultsDir/replicate_<idx>.rds.
## Used by both run styles:
##   single_machine/launch.sh   passes the replicate index as an argument;
##   slurm_array/replicate.sbatch  leaves it out, and the index comes from
##                                 $SLURM_ARRAY_TASK_ID.
##
## Usage: Rscript run_replicate.R [idx] [config.R path]
##   idx       integer in 1..M; if omitted, $SLURM_ARRAY_TASK_ID is used.
##   config.R  defaults to config.R in the current directory.
##
## The result file holds only the replicate's fitted components plus its
## log-likelihood and BIC, never the data, so it is small whatever the data
## size. Rerunning a replicate gives the same result (its seed comes from
## PHM::phmSeeds(CFG$seedBase, CFG$M)) and overwrites the same file.

cliArgs <- commandArgs(trailingOnly = TRUE)
usage <- "usage: Rscript run_replicate.R [idx] [config.R path]"
idxArg <- NA
if (length(cliArgs) >= 1 && grepl("^[0-9]+$", cliArgs[1])) {
  idxArg <- as.integer(cliArgs[1])
  cliArgs <- cliArgs[-1]
}
if (length(cliArgs) > 1) stop(usage)
configPath <- if (length(cliArgs) == 1) cliArgs[1] else "config.R"

if (!file.exists(configPath)) {
  stop(sprintf("config file '%s' not found -- run this from the directory holding config.R, or pass its path: %s", configPath, usage))
}
source(configPath)
suppressPackageStartupMessages(library(PHM))
if (!is.null(CFG$codeBase)) {
  message("note: CFG$codeBase is ignored; the fitting code comes from the installed PHM package")
}

idx <- if (!is.na(idxArg)) idxArg else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(idx)) {
  stop(sprintf("no replicate index: pass it as the first argument or set SLURM_ARRAY_TASK_ID -- %s", usage))
}
seeds <- PHM::phmSeeds(CFG$seedBase, CFG$M)
if (idx < 1 || idx > length(seeds)) {
  stop(sprintf("replicate index %d is out of range for CFG$M = %d in %s", idx, CFG$M, configPath))
}

X <- readRDS(CFG$dataPath)

result <- PHM::fitDensityReplicate(
  X, idx = idx, seed = seeds[idx], rho = CFG$rho, G = CFG$G,
  modelNames = CFG$modelNames, stage2Init = CFG$stage2Init,
  nRestartTries = CFG$nRestartTries, GStep = CFG$GStep, GMax = CFG$GMax,
  extendVerbose = CFG$extendVerbose, label = CFG$label
)

dir.create(CFG$resultsDir, showWarnings = FALSE, recursive = TRUE)
outFile <- file.path(CFG$resultsDir, sprintf("replicate_%03d.rds", idx))
saveRDS(result, outFile)

cat(sprintf("[run_replicate] idx=%d done -- G=%d (model=%s), loglik=%.2f -- wrote %s\n",
            idx, result$G, result$modelName, result$loglik, outFile))
