## k-means++ initial centers (D^2 weighting). Used by fitDensityReplicate()'s
## stage2Init = "kmeans++" Stage 2 random-restart option. Returns a
## k x ncol(X) matrix of chosen centers.
kmeansppInit <- function(X, k, seed) {
  set.seed(seed)
  n <- nrow(X)
  centers <- matrix(NA_real_, k, ncol(X))
  centers[1, ] <- X[sample.int(n, 1), ]
  d2 <- rowSums(sweep(X, 2, centers[1, ], "-")^2)
  for (i in 2:k) {
    probs <- d2 / sum(d2)
    idx <- sample.int(n, 1, prob = probs)
    centers[i, ] <- X[idx, ]
    d2 <- pmin(d2, rowSums(sweep(X, 2, centers[i, ], "-")^2))
  }
  centers
}

## Hard-assigns each row of X to the nearest row of centers (squared
## Euclidean distance), with no iterative refinement. Used by
## fitDensityReplicate()'s stage2Init = "randomFromData" Stage 2 random-restart
## option. Returns an integer vector of length nrow(X).
assignToNearestCenter <- function(X, centers) {
  d2 <- sapply(seq_len(nrow(centers)), function(j) rowSums(sweep(X, 2, centers[j, ], "-")^2))
  max.col(-d2, ties.method = "first")
}

## Automatic BIC-driven (model, G) search with range extension. Runs
## mclust::mclustBIC() over GInit, and if the BIC-optimal G lands on the
## upper boundary of the range tried, extends the search upward in steps
## of GStep (merging the new BIC cells onto the ones already computed via
## mclust::mclustBICupdate()) until an interior optimum is found or GMax
## is reached. Used by fitDensityReplicate() for each replicate's Stage 1 fit.
##
## This is not an approximation of a single wider search -- it computes
## exactly the same result, incrementally: mclustBICupdate() only
## concatenates BIC matrices computed on the same data, and
## Mclust(data, x = bic) just reads off the best cell, so no EM is ever
## re-run for a G value already tried in an earlier round.
##
## Hierarchical-clustering initialization for every Mclust()/mclustBIC()
## call follows mclust's own session-wide mclust.options("hcModelName")
## setting; set it with mclust::mclust.options() before calling this
## function (or fitDensityReplicate()/fitPHMEnsemble()) if a different
## initialization is desired. In a distributed run, this is a process-level
## R option -- set it identically in every worker's session, since it
## will not propagate across worker processes on its own.
##
## Returns an ordinary mclust::Mclust() object for the best (model, G)
## found over the (possibly extended) range.
mclustAutoExtend <- function(data, GInit = 1:25, modelNames = NULL,
                              GStep = 5, GMax = 100,
                              label = "mclust", verbose = TRUE, ...) {
  GHi <- max(GInit)
  bic <- mclust::mclustBIC(data, G = GInit, modelNames = modelNames, verbose = FALSE, ...)

  repeat {
    bestG <- as.integer(sub("^.*,([0-9]+)$", "\\1", names(mclust::pickBIC(bic, 1))))
    if (bestG < GHi || GHi >= GMax) break   # interior optimum, or hit the safety cap

    newHi <- min(GHi + GStep, GMax)
    if (verbose) {
      cat(sprintf("[%s] mclustAutoExtend: BIC-optimal G=%d at upper boundary of G<=%d -- extending search to %d:%d\n",
                  label, bestG, GHi, GHi + 1, newHi))
    }
    ## The extension can fail in three ways, all handled the same: keep the
    ## BIC table built so far and stop extending. (1) mclustBIC() itself
    ## errors. (2) It returns, but every new (model, G) cell is NA -- e.g.
    ## every EEE fit at the larger G went singular. Passing an all-NA table
    ## to mclustBICupdate() makes it compare NA > number and error out
    ## ("missing value where TRUE/FALSE needed"), which used to kill the
    ## whole replicate. (3) mclustBICupdate() errors for any other reason.
    ## Printed regardless of `verbose`, because the result is then capped at
    ## the current boundary and that is worth seeing in the log.
    bicExt <- tryCatch(
      mclust::mclustBIC(data, G = (GHi + 1):newHi, modelNames = modelNames, verbose = FALSE, ...),
      error = function(e) NULL)
    bicNew <- if (!is.null(bicExt) && any(is.finite(bicExt))) {
      tryCatch(mclust::mclustBICupdate(bic, bicExt), error = function(e) NULL)
    } else {
      NULL
    }
    if (is.null(bicNew)) {
      why <- if (is.null(bicExt)) "mclustBIC() failed" else
             if (!any(is.finite(bicExt))) "no model fit successfully" else "merging BIC tables failed"
      cat(sprintf("[%s] mclustAutoExtend: extension to G=%d:%d unusable (%s) -- keeping best fit with G<=%d (G=%d)\n",
                  label, GHi + 1, newHi, why, GHi, bestG))
      break
    }
    bic <- bicNew
    GHi <- newHi
  }

  if (verbose && bestG == GMax) {
    cat(sprintf("[%s] mclustAutoExtend: stopped at hard cap GMax=%d with BIC still best at the boundary -- consider raising GMax\n",
                label, GMax))
  }

  mclust::Mclust(data, x = bic)
}

## Converts the result of mclust::me() (an EM fit, not a full Mclust
## object) into the PHM-compatible per-component parameter list
## (prob/mean/var/class). Used by fitDensityReplicate()'s Stage 2. See
## constructPmcParamsMclust() for the output format.
buildParamsFromEM <- function(emFit) {
  K <- length(emFit$parameters$pro)
  lapply(seq_len(K), function(j) list(
    prob = emFit$parameters$pro[j],
    mean = emFit$parameters$mean[, j, drop = FALSE],
    var  = emFit$parameters$variance$sigma[, , j, drop = FALSE],
    class = j))
}

#' Fit one PHM ensemble replicate (Stage 1 + optional Stage 2)
#'
#' @description
#' The single-replicate "map" step of the PHM ensemble: draws a (sub)sample
#' of `X`, fits a `(model, G)` by BIC (with automatic range extension via
#' `mclustAutoExtend()`), optionally re-fits that `(model, G)` by EM from a
#' fresh random initialization (Stage 2), and returns just the resulting
#' per-component parameter list plus the replicate-level log-likelihood and
#' BIC used for weighting -- never the data itself, and never anything that
#' scales with the (sub)sample size. This is the unit of work
#' [fitPHMEnsemble()] runs `M` times; it has no dependency on any other
#' replicate's data or result, so it is the piece to parallelize across a
#' cluster, with [combineDensityReplicates()] as the matching reduce step
#' that needs all `M` results together.
#'
#' @details
#' Data for this replicate is a subsample of size `round(rho * n)`, drawn
#' from `X` without replacement by default (`replace = TRUE` gives a
#' bootstrap-style resample instead). At `rho == 1` without replacement this
#' draws all `n` rows -- a full permutation of `X`, equivalent to using `X`
#' as-is, since row order never affects any downstream step.
#'
#' If `X` passed here is already this replicate's pre-assigned shard (e.g.
#' a disjoint partition of a larger dataset, materialized once upstream and
#' handed to this worker on its own), call with `rho = 1` to use all of it.
#'
#' Stage 1 (always run): `(model, G)` is chosen via `mclustAutoExtend()`,
#' mclust's own (hc-init) BIC search with automatic range extension.
#'
#' Stage 2 (optional, controlled by `stage2Init`) is where this replicate's
#' contribution to ensemble diversity comes from when it runs: the chosen
#' `(model, G)` is re-fit by EM ([mclust::me()]) starting from a fresh
#' random initialization. Different random starts can converge to
#' different local optima of the same `(model, G)` -- this is why `M > 1`
#' replicates of [fitPHMEnsemble()] can still be meaningful even when
#' `rho == 1` (identical data every replicate): the randomness comes from
#' the EM initialization, not from the data. (With more than 2000 rows,
#' mclust's own Stage 1 initialization already uses a random subset of the
#' rows, `mclust.options("subset")`, so Stage 1 fits vary at `rho == 1` too;
#' with 2000 rows or fewer they are identical.) A degenerate fit (e.g. an
#' empty component) is retried with a fresh random seed, up to
#' `nRestartTries` times, before falling back to the plain Stage 1 fit.
#'
#' @param X Numeric data matrix (\eqn{N \times D}) this replicate draws its
#'   subsample from -- the whole ensemble's data, or (with `rho = 1`) just
#'   this replicate's own pre-assigned shard.
#' @param idx This replicate's index within the ensemble (`1..M`). Used
#'   only to tag component `class` labels so they stay distinct once pooled
#'   by [combineDensityReplicates()], and in progress messages -- never to
#'   derive randomness, so it has no effect on the fit itself.
#' @param seed Seed for this replicate's data subsampling and (if
#'   `stage2Init` is set) its random restart(s). Pass a pre-generated,
#'   deterministic per-replicate seed here -- e.g. `phmSeeds(seedBase, M)[idx]`,
#'   computed the same way by every worker -- rather than letting each worker choose
#'   its own, so a distributed run stays reproducible and matches what a
#'   single-process [fitPHMEnsemble()] call would have produced for the
#'   same `seedBase`.
#' @param rho Subsample fraction, in `(0, 1]`, drawn without replacement
#'   from `X` unless `replace = TRUE`. Use `rho = 1` when `X` already *is*
#'   this replicate's shard.
#' @param G Initial range of number-of-components values to try for this
#'   replicate's Stage 1 BIC search; automatically extended per
#'   `mclustAutoExtend()` if needed.
#' @param modelNames mclust covariance-model family/families to restrict
#'   the Stage 1 `(model, G)` BIC search to, e.g. `"EII"` or `"VVV"` -- or
#'   a vector of several. `NULL` (default) imposes no restriction.
#'   Whatever model mclust settles on is also what Stage 2's EM re-fit
#'   uses.
#' @param stage2Init Stage 2 random-restart method, or `NULL` to skip
#'   Stage 2 entirely (this replicate's contribution is then Stage 1's own
#'   converged fit, used as-is). One of `NULL`, `"kmeans"`, `"kmeans++"`,
#'   `"random"`, `"randomFromData"`:
#'   - `"kmeans"`: base R `kmeans()` with `Gk` random data rows as
#'     starting centers (Forgy init), run to convergence.
#'   - `"kmeans++"`: `kmeansppInit()` centers, run to convergence via
#'     `kmeans()`.
#'   - `"random"`: each point independently assigned a uniformly random
#'     hard label -- no spatial information at all, no clustering run.
#'   - `"randomFromData"`: `Gk` random data rows used directly as
#'     centers, with a single nearest-center assignment pass
#'     (`assignToNearestCenter()`) and no iteration.
#' @param nRestartTries Number of Stage 2 random-restart attempts before
#'   falling back to the Stage 1 fit.
#' @param GStep,GMax Passed to `mclustAutoExtend()`.
#' @param extendVerbose Whether `mclustAutoExtend()` should print progress
#'   messages as it extends its search range.
#' @param label Used only in progress messages.
#' @param replace Whether to draw the subsample with replacement (a bootstrap
#'   resample when `rho = 1`). Default `FALSE`.
#' @param ... Further arguments passed to [mclust::mclustBIC()] in the Stage 1
#'   search.
#'
#' @return A list with elements:
#'   - `params`: this replicate's per-component `prob`/`mean`/`var`/`class`
#'     list (`class` already suffixed with `idx`, so pooling several
#'     replicates' `params` never collides).
#'   - `loglik`, `bic`: this replicate's log-likelihood and BIC, as used by
#'     [combineDensityReplicates()]'s `"likelihood"`/`"bic"` weighting.
#'   - `idx`, `G`, `modelName`: bookkeeping, echoed back for convenience.
#'
#'   Small and serialization-friendly regardless of how large `X`/the
#'   subsample was -- this is what a cluster map task should write out
#'   (e.g. as a single small `.rds` file) for the reduce step to collect.
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 5), rnorm(200, -5)), ncol = 2, byrow = TRUE)
#' rep1 <- fitDensityReplicate(dat, idx = 1, seed = 123, rho = 0.8, G = 1:5)
#'
#' @export
fitDensityReplicate <- function(X, idx, seed, rho = 1.0, G = 1:25, modelNames = NULL,
                                 stage2Init = NULL, nRestartTries = 5,
                                 GStep = 5, GMax = 100, extendVerbose = TRUE,
                                 label = "example", replace = FALSE, ...) {
  stage2Methods <- c("kmeans", "kmeans++", "random", "randomFromData")
  if (!is.null(stage2Init) && !stage2Init %in% stage2Methods) {
    stop(sprintf("stage2Init must be NULL or one of: %s (got '%s')",
                  paste(stage2Methods, collapse = ", "), stage2Init))
  }
  if (!is.null(modelNames)) {
    validModels <- c("E", "V", mclust::mclust.options("emModelNames"))
    bad <- setdiff(modelNames, validModels)
    if (length(bad) > 0) {
      stop(sprintf("modelNames contains unrecognized mclust model name(s): %s (valid: %s)",
                    paste(bad, collapse = ", "), paste(validModels, collapse = ", ")))
    }
  }

  n <- nrow(X)
  set.seed(seed)

  ## ---- data for this replicate ----
  ## A single draw of round(rho*n) rows, without replacement unless
  ## replace = TRUE. Without replacement, at rho == 1
  ## this draws all n rows without replacement -- a full permutation of
  ## X, i.e. every row exactly once -- so it needs no special case: row
  ## order never affects any downstream step (Mclust/kmeans/me all treat
  ## rows as an unordered collection), so "all of X, permuted" behaves
  ## identically to "all of X, as-is".
  subIdx <- sample.int(n, round(rho * n), replace = replace)
  subDat <- X[subIdx, , drop = FALSE]

  ## ---- Stage 1 (always run): (model, G) chosen via mclust's own
  ## default (hc-init) BIC search, auto-extending the G range upward if
  ## the optimum sits at its boundary (see mclustAutoExtend() above) ----
  mcl0 <- mclustAutoExtend(subDat, GInit = G, modelNames = modelNames,
                            GStep = GStep, GMax = GMax,
                            label = paste0(label, " rep", idx), verbose = extendVerbose, ...)
  Gk <- mcl0$G
  modelName <- mcl0$modelName

  if (is.null(stage2Init)) {
    ## ---- no Stage 2: Stage 1's own converged fit, used as-is ----
    params <- constructPmcParamsMclust(mcl0)
    loglik <- mcl0$loglik
    bicVal <- mcl0$bic
    cat(sprintf("[%s] PHM replicate %d: G=%d (model=%s), loglik=%.2f (no Stage 2)\n",
                label, idx, Gk, modelName, mcl0$loglik))
  } else {
    ## ---- Stage 2: re-fit EM from a FRESH random initialization at the
    ## chosen (model, G) -- see stage2Init options above. Retry with a
    ## new random seed if a given start produces a degenerate fit, up to
    ## nRestartTries times, before falling back to the plain mcl0 fit. ----
    emr <- NULL
    for (attempt in seq_len(nRestartTries)) {
      z0 <- switch(stage2Init,
        "kmeans" = {
          km <- stats::kmeans(subDat, centers = Gk, nstart = 1, iter.max = 100)
          mclust::unmap(km$cluster, groups = seq_len(Gk))
        },
        "kmeans++" = {
          set.seed(seed * 100 + attempt)
          init <- kmeansppInit(subDat, Gk, seed * 100 + attempt)
          km <- stats::kmeans(subDat, centers = init, iter.max = 100)
          mclust::unmap(km$cluster, groups = seq_len(Gk))
        },
        "random" = {
          set.seed(seed * 100 + attempt)
          labels <- sample.int(Gk, nrow(subDat), replace = TRUE)
          mclust::unmap(labels, groups = seq_len(Gk))
        },
        "randomFromData" = {
          set.seed(seed * 100 + attempt)
          ctrIdx <- sample.int(nrow(subDat), Gk)
          labels <- assignToNearestCenter(subDat, subDat[ctrIdx, , drop = FALSE])
          mclust::unmap(labels, groups = seq_len(Gk))
        })
      emr <- tryCatch(mclust::me(modelName = modelName, data = subDat, z = z0),
                       error = function(e) NULL)
      if (!is.null(emr) && is.finite(emr$loglik)) break
      emr <- NULL
    }
    if (is.null(emr)) {
      cat(sprintf("[%s] PHM replicate %d: random-restart EM (%s) failed after %d tries, falling back to mclust's own fit\n",
                  label, idx, stage2Init, nRestartTries))
      params <- constructPmcParamsMclust(mcl0)
      loglik <- mcl0$loglik
      bicVal <- mcl0$bic
    } else {
      params <- buildParamsFromEM(emr)
      loglik <- emr$loglik
      bicVal <- mclust::bic(modelName = modelName, loglik = emr$loglik,
                             n = nrow(subDat), d = ncol(subDat), G = Gk,
                             equalPro = FALSE)
    }
    cat(sprintf("[%s] PHM replicate %d: G=%d (model=%s)%s\n", label, idx, Gk, modelName,
                if (!is.null(emr)) sprintf(", loglik=%.2f [%s]", emr$loglik, stage2Init) else " (fallback fit)"))
  }

  for (g in seq_along(params)) params[[g]]$class <- paste(params[[g]]$class, idx, sep = "_")

  list(params = params, loglik = loglik, bic = bicVal, idx = idx, G = Gk, modelName = modelName)
}

## Pools the M replicate results from fitDensityReplicate() into one
## component list, scaling each replicate's (already normalized to 1)
## component probabilities by that replicate's weight, so the pooled mass
## sums to 1. Weights: "equal" (1/M each), "likelihood" (softmax over the
## replicates' log-likelihoods) or "bic" (softmax over BIC/2), both
## log-sum-exp stabilized. Shared by combineDensityReplicates() and
## globalDensityEstimation().
poolDensityReplicates <- function(replicateResults, replicateWeighting = "equal",
                                  label = "example") {
  M <- length(replicateResults)
  replicateLoglik <- vapply(replicateResults, function(r) r$loglik, numeric(1))
  replicateBic    <- vapply(replicateResults, function(r) r$bic, numeric(1))
  allParams       <- lapply(replicateResults, function(r) r$params)

  if (replicateWeighting == "equal") {
    repW <- rep(1 / M, M)
  } else {
    score <- switch(replicateWeighting,
      "likelihood" = replicateLoglik,
      "bic"        = replicateBic / 2)
    shift <- max(score)
    repWUnnorm <- exp(score - shift)
    repW <- repWUnnorm / sum(repWUnnorm)
    cat(sprintf("[%s] PHM: replicate weights (%s-based): %s\n",
                label, replicateWeighting, paste(sprintf("%.4g", repW), collapse = ", ")))
  }

  ensParams <- do.call(c, allParams)
  repIdxOf <- rep(seq_len(M), times = vapply(allParams, length, integer(1)))
  for (i in seq_along(ensParams)) ensParams[[i]]$prob <- ensParams[[i]]$prob * repW[repIdxOf[i]]
  ensParams
}

#' Combine PHM ensemble replicates into one PHMv2 fit
#'
#' @description
#' The "reduce" step of the PHM ensemble: takes the `M` small per-replicate
#' results produced independently by [fitDensityReplicate()] -- one per
#' worker, in a distributed run, collected back together -- computes the
#' replicate weights, pools every replicate's components into a single
#' list, and fits [PHMv2()] on the pooled ensemble. Needs only what
#' [fitDensityReplicate()] returned, never the original data, so this
#' step's cost scales with the total number of pooled components (roughly
#' `M` times the average `G` per replicate) and `mcSamples`, not with the
#' dataset size -- it is a single, comparatively cheap job even when the
#' `M` replicates were fit on a very large dataset in total.
#'
#' Must be run only after all `M` replicates have finished: unlike
#' [fitDensityReplicate()], this step is not itself parallelizable (the
#' `"likelihood"`/`"bic"` weighting needs every replicate's score together,
#' and `PHMv2()` fits a single ensemble), so a cluster job graph should
#' have the `M` [fitDensityReplicate()] tasks as a parallel "map" stage,
#' followed by exactly one call to this function as the "reduce" stage.
#'
#' @param replicateResults A list of length `M`, each element the return
#'   value of a [fitDensityReplicate()] call (i.e. a list with `params`,
#'   `loglik`, `bic`, `idx`, ...). Order in this list does not matter --
#'   each element's own `idx` (baked into its components' `class` labels)
#'   is what distinguishes replicates, not position.
#' @param replicateWeighting How much total probability mass each of the
#'   `M` replicates contributes to the pooled ensemble, before [PHMv2()]
#'   is fit. One of:
#'   - `"equal"`: default. Every replicate contributes exactly `1/M` of
#'     the total mass, regardless of how many components it has or how
#'     well it fits.
#'   - `"likelihood"`: each replicate is weighted by a softmax over the
#'     `M` replicate-level log-likelihoods (log-sum-exp stabilized, since
#'     raw log-likelihoods are joint densities over an entire (sub)sample
#'     and can run into the thousands-negative range). Only a fair
#'     comparison when every replicate's (sub)sample was the same size.
#'   - `"bic"`: like `"likelihood"`, but weights by the BIC-penalized
#'     log-likelihood (`bic / 2`, on the same scale as plain log-
#'     likelihood) instead of the raw log-likelihood, so a replicate
#'     doesn't get more say just because `mclustAutoExtend()` let it
#'     settle on a larger `G`. Same same-sample-size caveat as
#'     `"likelihood"`.
#' @param numCores Number of cores for [PHMv2()]'s Monte Carlo
#'   \eqn{\Delta P_{\rm mc}} computation.
#' @param mcSamples Monte Carlo sample size passed to [PHMv2()].
#' @param batchSize Monte Carlo batch size passed to [PHMv2()].
#' @param label Used only in progress messages.
#'
#' @return The fitted [PHMv2()] object; `$paramsList` holds the pooled
#'   components from all `M` replicates.
#'
#' @export
combineDensityReplicates <- function(replicateResults,
                                      replicateWeighting = c("equal", "likelihood", "bic"),
                                      numCores = 2, mcSamples = 2e5, batchSize = 1e5,
                                      label = "example") {
  replicateWeighting <- match.arg(replicateWeighting)

  M <- length(replicateResults)
  ensParams <- poolDensityReplicates(replicateResults, replicateWeighting, label)

  cat(sprintf("[%s] PHM: fitting PHMv2 on %d pooled components from %d replicates...\n",
              label, length(ensParams), M))
  PHMv2(paramsList = ensParams, scaling = "alpha", mcSamples = mcSamples,
        batchSize = batchSize, numCores = numCores)
}

#' Per-replicate seeds for a PHM ensemble
#'
#' @description
#' The `M` seeds used for the replicates of a PHM ensemble, drawn from a
#' single base seed. [fitPHMEnsemble()] uses this internally, and every script
#' in `inst/ensemble/` calls it with the same `seedBase` and `M`, so replicate
#' `idx` always gets `phmSeeds(seedBase, M)[idx]` no matter which process or
#' machine fits it. Nothing about the seeds has to be written to a file or
#' passed between jobs.
#'
#' @param seedBase Base seed.
#' @param M Number of replicates.
#'
#' @return Integer vector of length `M`.
#'
#' @examples
#' phmSeeds(1000, 5)
#'
#' @export
phmSeeds <- function(seedBase, M) {
  set.seed(seedBase)
  sample.int(1e7, M)
}

#' Fit a PHM ensemble
#'
#' @description
#' Fits an ensemble of `M` independent mixture-model replicates and pools
#' their components into a single [PHMv2()] fit. Each replicate's density
#' is estimated from a (sub)sample of the data, with a `(model, G)`
#' search that automatically extends its range when BIC wants more
#' components than were tried (see `mclustAutoExtend()`); optionally, a
#' second stage re-fits that `(model, G)` by EM from a fresh random
#' initialization, which is where most of the ensemble's diversity comes
#' from when the subsample fraction `rho` is 1.
#'
#' @details
#' This function is now a thin, serial driver: it generates the `M`
#' per-replicate seeds, calls [fitDensityReplicate()] once per replicate in
#' an ordinary `lapply()`, and passes all `M` results to
#' [combineDensityReplicates()]. It is kept for convenience and backward
#' compatibility -- for example, from an interactive session, or when `M`
#' replicates comfortably fit in one process. To run the `M` replicates in
#' parallel (multicore, or across a cluster), call [fitDensityReplicate()]
#' directly for each replicate -- each call is fully independent and needs
#' nothing from any other replicate -- then call
#' [combineDensityReplicates()] once, directly, on the collected results,
#' instead of calling this function.
#'
#' @param X Numeric data matrix (\eqn{N \times D})
#' @param M Number of ensemble replicates
#' @param rho Subsample fraction per replicate, in `(0, 1]`, drawn without
#'   replacement
#' @param G Initial range of number-of-components values to try for each
#'   replicate's Stage 1 BIC search; automatically extended per
#'   `mclustAutoExtend()` if needed
#' @param modelNames mclust covariance-model family/families to restrict
#'   the Stage 1 `(model, G)` BIC search to, e.g. `"EII"` or `"VVV"` -- or
#'   a vector of several. `NULL` (default) imposes no restriction.
#'   Whatever model mclust settles on for a replicate is also what Stage
#'   2's EM re-fit uses for that replicate.
#' @param replicateWeighting See [combineDensityReplicates()].
#' @param seedBase Seed used to generate the `M` per-replicate seeds
#'   (passed to [fitDensityReplicate()] as its `seed` argument)
#' @param numCores Number of cores for [PHMv2()]'s Monte Carlo
#'   \eqn{\Delta P_{\rm mc}} computation
#' @param mcSamples Monte Carlo sample size passed to [PHMv2()]
#' @param batchSize Monte Carlo batch size passed to [PHMv2()]
#' @param stage2Init Stage 2 random-restart method, or `NULL` to skip
#'   Stage 2 entirely. See [fitDensityReplicate()].
#' @param nRestartTries Number of Stage 2 random-restart attempts before
#'   falling back to the Stage 1 fit for a replicate
#' @param GStep,GMax Passed to `mclustAutoExtend()`
#' @param extendVerbose Whether `mclustAutoExtend()` should print progress
#'   messages as it extends its search range
#' @param label Used only in progress messages
#'
#' @return The fitted [PHMv2()] object; `$paramsList` holds the pooled
#'   components from all `M` replicates.
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 5), rnorm(200, -5)), ncol = 2, byrow = TRUE)
#' phm <- fitPHMEnsemble(dat, M = 2, rho = 0.8, G = 1:5, numCores = 1)
#'
#' @export
fitPHMEnsemble <- function(X, M = 5, rho = 1.0, G = 1:25, modelNames = NULL,
                            replicateWeighting = c("equal", "likelihood", "bic"),
                            seedBase = 1000,
                            numCores = 2, mcSamples = 2e5, batchSize = 1e5,
                            stage2Init = NULL, nRestartTries = 5,
                            GStep = 5, GMax = 100, extendVerbose = TRUE,
                            label = "example") {
  replicateWeighting <- match.arg(replicateWeighting)

  ## ---- M per-replicate seeds, generated once, centrally -- this is the
  ## only piece of state every replicate needs that isn't already in its
  ## own arguments, which is exactly why fitDensityReplicate() takes a
  ## single pre-generated seed rather than (seedBase, idx): a distributed
  ## caller (e.g. the scripts in inst/ensemble/) calls phmSeeds() with the
  ## same seedBase and M and hands seeds[idx] to each worker. ----
  seeds <- phmSeeds(seedBase, M)

  ## ---- "map": M independent replicate fits. Serial here (plain lapply);
  ## a distributed caller runs this loop's body as M separate cluster/
  ## multicore tasks instead and collects their results the same way. ----
  replicateResults <- lapply(seq_len(M), function(idx) {
    fitDensityReplicate(X, idx = idx, seed = seeds[idx], rho = rho, G = G,
                         modelNames = modelNames, stage2Init = stage2Init,
                         nRestartTries = nRestartTries, GStep = GStep, GMax = GMax,
                         extendVerbose = extendVerbose, label = label)
  })

  ## ---- "reduce": one combine step over all M results ----
  combineDensityReplicates(replicateResults, replicateWeighting = replicateWeighting,
                            numCores = numCores, mcSamples = mcSamples,
                            batchSize = batchSize, label = label)
}
