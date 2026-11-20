#!/usr/bin/env Rscript
## local_driver.R -- single-script alternative to launch.sh +
## run_replicate.R + run_combine.R: fits all M replicates in parallel
## within ONE R process, using base R's own parallel::mclapply() (ships
## with R -- no extra packages to install), then
## runs the combine step itself once they're all back. No background job
## files to gather, and no separate worker script.
##
## mclapply() works by forking the current R process (macOS/Linux only --
## NOT Windows; see the check below). Each fork starts as a copy-on-write
## clone of this process's memory, so X doesn't need to be re-read or
## re-serialized per replicate the way it's re-read from disk in
## run_replicate.R, but every replicate's own working set still adds up
## across forks running at once.
##
## Trade-off vs. launch.sh: a crashed replicate is reported inline (see
## the try-error check below) but you don't get launch.sh's one-log-file-
## per-replicate for sifting through output when several things go wrong
## at once.
##
## mc.preschedule = FALSE below is deliberate, not the mclapply() default:
## with the default (TRUE), mclapply() groups several replicates onto
## each forked worker and runs them sequentially inside it, and I
## confirmed (by injecting a failure) that if one replicate in a group
## hard-crashes its worker, mclapply() marks EVERY replicate in that same
## group as failed too -- including ones that already finished
## successfully before the crash. mc.preschedule = FALSE gives one fork
## per replicate instead, so a crash is isolated to the one replicate
## that actually failed, matching launch.sh's per-process isolation. The
## cost is one fork's overhead per replicate rather than per core, which
## is negligible next to an actual density-fitting replicate's own runtime.
##
## Usage: Rscript local_driver.R [config.R path]
##        MAX_CONCURRENT=4 Rscript local_driver.R [config.R path]
## config.R is looked up in the current working directory by default;
## pass its path explicitly as the only argument to use one that lives
## somewhere else.

cliArgs <- commandArgs(trailingOnly = TRUE)
if (length(cliArgs) > 1) {
  stop("usage: Rscript local_driver.R [config.R path]")
}
configPath <- if (length(cliArgs) >= 1) cliArgs[1] else "config.R"

if (!file.exists(configPath)) {
  stop(sprintf("config file '%s' not found -- run this from the directory holding config.R, or pass its path: Rscript local_driver.R <config.R path>", configPath))
}
source(configPath)
suppressPackageStartupMessages(library(PHM))
if (!is.null(CFG$codeBase)) {
  message("note: CFG$codeBase is ignored; the fitting code comes from the installed PHM package")
}

if (.Platform$OS.type == "windows") {
  stop("parallel::mclapply() forks processes and does not work on Windows -- use parallel::makeCluster()/parLapply() instead, or one of the other parallel_run/ styles")
}

X <- readRDS(CFG$dataPath)
seeds <- PHM::phmSeeds(CFG$seedBase, CFG$M)

maxCores <- as.integer(Sys.getenv("MAX_CONCURRENT", parallel::detectCores()))
cat(sprintf("[local_driver] fitting %d replicates, up to %d at a time\n", CFG$M, maxCores))

replicateResults <- parallel::mclapply(seq_len(CFG$M), function(idx) {
  PHM::fitDensityReplicate(X, idx = idx, seed = seeds[idx], rho = CFG$rho, G = CFG$G,
                       modelNames = CFG$modelNames, stage2Init = CFG$stage2Init,
                       nRestartTries = CFG$nRestartTries, GStep = CFG$GStep, GMax = CFG$GMax,
                       extendVerbose = CFG$extendVerbose, label = CFG$label)
}, mc.cores = maxCores, mc.preschedule = FALSE)

## mclapply() returns a try-error object for any replicate that errored,
## and NULL for one whose fork was killed (e.g. by the operating system for
## running out of memory), rather than failing the whole call -- check for
## both before pooling, so a broken replicate never ends up in the ensemble.
failed <- vapply(replicateResults, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("replicate(s) %s failed -- see the error printed above; a replicate with no error message was most likely killed for running out of memory (try a smaller MAX_CONCURRENT)", paste(which(failed), collapse = ", ")))
}

fit <- PHM::combineDensityReplicates(
  replicateResults, replicateWeighting = CFG$replicateWeighting,
  numCores = CFG$numCores, mcSamples = CFG$mcSamples, batchSize = CFG$batchSize,
  label = CFG$label
)

saveRDS(fit, CFG$outPath)
cat(sprintf("[local_driver] wrote %s (%d pooled components from %d replicates)\n",
            CFG$outPath, length(fit$paramsList), CFG$M))
