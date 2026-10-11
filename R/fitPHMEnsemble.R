## k-means++ initial centers (D^2 weighting). Used by fitPHMEnsemble()'s
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
## fitPHMEnsemble()'s stage2Init = "randomFromData" Stage 2 random-restart
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
## is reached. Used by fitPHMEnsemble() for each replicate's Stage 1 fit.
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
## function (or fitPHMEnsemble()) if a different initialization is
## desired.
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
    bicExt <- tryCatch(
      mclust::mclustBIC(data, G = (GHi + 1):newHi, modelNames = modelNames, verbose = FALSE, ...),
      error = function(e) NULL)
    if (is.null(bicExt)) {
      if (verbose) {
        cat(sprintf("[%s] mclustAutoExtend: extension to %d:%d failed -- stopping with G<=%d\n",
                    label, GHi + 1, newHi, GHi))
      }
      break
    }
    bic <- mclust::mclustBICupdate(bic, bicExt)
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
## (prob/mean/var/class). Used by fitPHMEnsemble()'s Stage 2. See
## constructPmcParamsMclust() for the output format.
buildParamsFromEM <- function(emFit) {
  K <- length(emFit$parameters$pro)
  lapply(seq_len(K), function(j) list(
    prob = emFit$parameters$pro[j],
    mean = emFit$parameters$mean[, j, drop = FALSE],
    var  = emFit$parameters$variance$sigma[, , j, drop = FALSE],
    class = j))
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
#' Data per replicate is a subsample of size `round(rho * n)`, drawn
#' *without* replacement (not a bootstrap resample). At `rho == 1` this
#' draws all `n` rows without replacement -- a full permutation of `X`,
#' equivalent to using `X` as-is, since row order never affects any
#' downstream step.
#'
#' Stage 1 (always run): `(model, G)` is chosen via `mclustAutoExtend()`,
#' mclust's own (hc-init) BIC search with automatic range extension.
#'
#' Stage 2 (optional, controlled by `stage2Init`) is where ensemble
#' diversity comes from when it runs: the chosen `(model, G)` is re-fit
#' by EM ([mclust::me()]) starting from a fresh, replicate-specific
#' random initialization. Different random starts can converge to
#' different local optima of the same `(model, G)` -- this is why `M >
#' 1` can still be meaningful even when `rho == 1` (identical data every
#' replicate): the randomness comes from the EM initialization, not from
#' the data. A degenerate fit (e.g. an empty component) is retried with a
#' fresh random seed, up to `nRestartTries` times, before falling back to
#' the plain Stage 1 fit for that replicate.
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
#' @param replicateWeighting How much total probability mass each of the
#'   `M` replicates contributes to the pooled ensemble, before [PHMv2()]
#'   is fit. One of:
#'   - `"equal"`: default. Every replicate contributes exactly `1/M` of
#'     the total mass, regardless of how many components it has or how
#'     well it fits.
#'   - `"likelihood"`: each replicate is weighted by a softmax over the
#'     `M` replicate-level log-likelihoods (log-sum-exp stabilized, since
#'     raw log-likelihoods are joint densities over an entire (sub)sample
#'     and can run into the thousands-negative range).
#'   - `"bic"`: like `"likelihood"`, but weights by the BIC-penalized
#'     log-likelihood (`bic / 2`, on the same scale as plain log-
#'     likelihood) instead of the raw log-likelihood, so a replicate
#'     doesn't get more say just because `mclustAutoExtend()` let it
#'     settle on a larger `G`.
#' @param seedBase Seed used to generate the `M` per-replicate seeds
#' @param numCores Number of cores for [PHMv2()]'s Monte Carlo
#'   \eqn{\Delta P_{\rm mc}} computation
#' @param mcSamples Monte Carlo sample size passed to [PHMv2()]
#' @param batchSize Monte Carlo batch size passed to [PHMv2()]
#' @param stage2Init Stage 2 random-restart method, or `NULL` to skip
#'   Stage 2 entirely (each replicate's contribution is then Stage 1's
#'   own converged fit, used as-is). One of `NULL`, `"kmeans"`,
#'   `"kmeans++"`, `"random"`, `"randomFromData"`:
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
  set.seed(seedBase)
  seeds <- sample.int(1e7, M)

  allParams <- vector("list", M)
  replicateLoglik <- numeric(M)
  replicateBic <- numeric(M)
  for (idx in seq_len(M)) {
    set.seed(seeds[idx])

    ## ---- data for this replicate ----
    ## A single without-replacement draw of round(rho*n) rows. At rho == 1
    ## this draws all n rows without replacement -- a full permutation of
    ## X, i.e. every row exactly once -- so it needs no special case: row
    ## order never affects any downstream step (Mclust/kmeans/me all treat
    ## rows as an unordered collection), so "all of X, permuted" behaves
    ## identically to "all of X, as-is".
    subIdx <- sample.int(n, round(rho * n), replace = FALSE)
    subDat <- X[subIdx, ]

    ## ---- Stage 1 (always run): (model, G) chosen via mclust's own
    ## default (hc-init) BIC search, auto-extending the G range upward if
    ## the optimum sits at its boundary (see mclustAutoExtend() above) ----
    mcl0 <- mclustAutoExtend(subDat, GInit = G, modelNames = modelNames,
                              GStep = GStep, GMax = GMax,
                              label = paste0(label, " rep", idx), verbose = extendVerbose)
    Gk <- mcl0$G
    modelName <- mcl0$modelName

    if (is.null(stage2Init)) {
      ## ---- no Stage 2: Stage 1's own converged fit, used as-is ----
      params <- constructPmcParamsMclust(mcl0)
      replicateLoglik[idx] <- mcl0$loglik
      replicateBic[idx] <- mcl0$bic
      cat(sprintf("[%s] PHM replicate %d/%d: G=%d (model=%s), loglik=%.2f (no Stage 2)\n",
                  label, idx, M, Gk, modelName, mcl0$loglik))
    } else {
      ## ---- Stage 2: re-fit EM from a FRESH random initialization at the
      ## chosen (model, G) -- see stage2Init options above. Retry with a
      ## new random seed if a given start produces a degenerate fit, up to
      ## nRestartTries times, before falling back to the plain mcl0 fit. ----
      emr <- NULL
      for (attempt in seq_len(nRestartTries)) {
        z0 <- switch(stage2Init,
          "kmeans" = {
            km <- kmeans(subDat, centers = Gk, nstart = 1, iter.max = 100)
            mclust::unmap(km$cluster, groups = seq_len(Gk))
          },
          "kmeans++" = {
            set.seed(seeds[idx] * 100 + attempt)
            init <- kmeansppInit(subDat, Gk, seeds[idx] * 100 + attempt)
            km <- kmeans(subDat, centers = init, iter.max = 100)
            mclust::unmap(km$cluster, groups = seq_len(Gk))
          },
          "random" = {
            set.seed(seeds[idx] * 100 + attempt)
            labels <- sample.int(Gk, nrow(subDat), replace = TRUE)
            mclust::unmap(labels, groups = seq_len(Gk))
          },
          "randomFromData" = {
            set.seed(seeds[idx] * 100 + attempt)
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
        cat(sprintf("[%s] PHM replicate %d/%d: random-restart EM (%s) failed after %d tries, falling back to mclust's own fit\n",
                    label, idx, M, stage2Init, nRestartTries))
        params <- constructPmcParamsMclust(mcl0)
        replicateLoglik[idx] <- mcl0$loglik
        replicateBic[idx] <- mcl0$bic
      } else {
        params <- buildParamsFromEM(emr)
        replicateLoglik[idx] <- emr$loglik
        replicateBic[idx] <- mclust::bic(modelName = modelName, loglik = emr$loglik,
                                          n = nrow(subDat), d = ncol(subDat), G = Gk,
                                          equalPro = FALSE)
      }
      cat(sprintf("[%s] PHM replicate %d/%d: G=%d (model=%s)%s\n", label, idx, M, Gk, modelName,
                  if (!is.null(emr)) sprintf(", loglik=%.2f [%s]", emr$loglik, stage2Init) else " (fallback fit)"))
    }
    for (g in seq_along(params)) params[[g]]$class <- paste(params[[g]]$class, idx, sep = "_")
    allParams[[idx]] <- params
  }

  ## ---- replicate weights: "equal" (1/M each, the original behavior),
  ## "likelihood" (softmax over the M replicate-level log-likelihoods), or
  ## "bic" (softmax over the M replicate-level BIC/2 values). Both
  ## non-equal schemes use the log-sum-exp shift for numerical stability
  ## -- see @param replicateWeighting above. Each weight replaces the flat
  ## 1/M scale factor applied to that replicate's own (already-normalized-
  ## to-1) component probabilities, so the M replicates' total pooled mass
  ## still sums to 1 under every scheme. ----
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

  cat(sprintf("[%s] PHM: fitting PHMv2 on %d pooled components...\n", label, length(ensParams)))
  PHMv2(paramsList = ensParams, scaling = "alpha", mcSamples = mcSamples,
        batchSize = batchSize, numCores = numCores)
}
