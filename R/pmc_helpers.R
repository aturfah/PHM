sampleDistbn <- function(distbn_params, num_samples) {
  mean <- distbn_params$mean
  var <- distbn_params$var
  prob <- distbn_params$prob

  if (is.null(nrow(mean))) mean <- matrix(mean, nrow=1, ncol=1)

  D <- nrow(mean)
  K <- length(prob)

  if (K > 1) {
    ## Split the samples across each sub-distribution
    Nk_vec <- c(rmultinom(1, num_samples, prob))

    ## With the counts sample from the subdivision
    output <- lapply(1:K, function(idx) {
      if (Nk_vec[idx] == 0) return(NULL)
      sampleDistbn(list(
        mean=mean[, idx, drop=F],
        var=var[, , idx, drop=F],
        prob=1
      ),
      Nk_vec[idx]
      )
    })
    return(do.call(rbind, output))
  }

  mvtnorm::rmvnorm(num_samples, mean, matrix(var[, , 1], ncol=D, nrow=D))
}

sampleMixture <- function(paramsList, num_samples) {
  prob_vec <- sapply(paramsList, function(x) sum(x$prob))
  Nk_vec <- c(stats::rmultinom(1, num_samples, prob_vec))

  paramsList[which(Nk_vec == 0)] <- NULL ## Don't sample components with 0 samples
  Nk_vec <- Nk_vec[which(Nk_vec != 0)]

  output <- lapply(1:length(paramsList), function(idx) {
    distbn_params <- paramsList[[idx]]
    sampleDistbn(distbn_params, Nk_vec[idx])
  })

  do.call(rbind, output)
}

generateDistbnFunc <- function(distbn_params) {
  ## NOTE: This takes one element of the list paramsList that everything else takes

  ## Normalize the probabilities of each component
  M <-  length(distbn_params$prob)
  rel_prop <- distbn_params$prob / sum(distbn_params$prob)

  ## Output function
  function(x) {
    D <- ncol(x)
    res <- sapply(1:M, function(idx) {
      rel_prop[idx] * mvtnorm::dmvnorm(x,
                              distbn_params$mean[, idx, drop=F],
                              matrix(
                                distbn_params$var[, , idx, drop=F],
                                nrow=D, ncol=D
                              ))
    })

    if (is.null(dim(res))) return(sum(res))

    return(apply(res, 1, sum))
  }
}

generatePosteriorProbFunc <- function(paramsList, j) {
  function(x) {
    if (is.null(dim(x))) stop("x must be an NxD matrix")

    num <- sum(paramsList[[j]]$prob) * generateDistbnFunc(paramsList[[j]])(x)
    den <- 0
    K <- length(paramsList)
    for (k in 1:K) {
      den <- den +
        sum(paramsList[[k]]$prob) * generateDistbnFunc(paramsList[[k]])(x)
    }
    den <- ifelse(den == 0, 1, den) ## Should probably be able to do this in a MC way

    return(num / den)
  }
}

#' Compute Posterior Matrix for a set of parameteres
#'
#' @param paramsList List of mixture parameteres
#' @param data Data matrix for which to compute posterior probabilities. Rows are observations
#'
#' @description Compute the posterior probability matrix
#'
#' @return An \eqn{N\times K} matrix
#' @export
computePosteriorProbMatrix <- function(paramsList, data) {
  if (is.null(dim(data))) data <- matrix(data, ncol=1)
  K <- length(paramsList)
  distbn.mat <- do.call(cbind,
                       lapply(1:K, function(j) {
                         densJ <- generateDistbnFunc(paramsList[[j]])
                         sum(paramsList[[j]]$prob) * densJ(data)
                       }))
  t(apply(distbn.mat, 1, function(x) if (sum(x) != 0) x / sum(x) else x + 1 / length(x)))
}


## parallel::mclapply() does not stop when a worker fails. A task that
## throws an error comes back as a "try-error" object, and a worker killed by
## the operating system (most often for running out of memory) comes back as
## NULL, with only a warning. Code that then combines the results with
## rbind() or c() silently drops that part of the work -- for PHMv2() this
## means Monte Carlo batches go missing and every DeltaPmc comes out too
## small. parallelApply() runs mclapply() and stops with an error instead.
parallelApply <- function(X, FUN, numCores, ...) {
  res <- parallel::mclapply(X, FUN, mc.cores = numCores, ...)
  failed <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
  if (any(failed)) {
    errs <- vapply(res[failed], function(r) {
      if (inherits(r, "try-error")) trimws(conditionMessage(attr(r, "condition"))) else "no result returned"
    }, character(1))
    killed <- vapply(res, is.null, logical(1))
    stop(sprintf("%d of %d parallel task(s) failed (first: %s).%s",
                 sum(failed), length(res), errs[1],
                 if (any(killed)) paste0(" A task with no result was most likely killed by the ",
                                         "operating system for running out of memory; try a ",
                                         "smaller numCores or batchSize.") else ""),
         call. = FALSE)
  }
  res
}

posteriorMatrixMCPmc <- function(paramsList, mcSamples, batchSize, numCores, verbose) {
  num_batches <- ceiling(mcSamples / batchSize)
  applyFunc <- lapply
  if (numCores > 1) applyFunc <- function(x, func) parallelApply(x, func, numCores)

  if (verbose) cat("Drawing MC Samples\n")
  mc_obs <- sampleMixture(paramsList, mcSamples)
  if (is.null(dim(mc_obs))) mc_obs <- matrix(mc_obs, ncol=1)

  if (verbose) cat("Computing MC Posterior Matrix\n")
  if (num_batches == 1) {
    post_mat <- computePosteriorProbMatrix(paramsList, mc_obs)
  } else {
    batch_idx <- t(sapply(1:num_batches, function(idx) {
      lb <- 1 + batchSize * (idx - 1)
      ub <- min(batchSize * idx, mcSamples)
      c(lb, ub)
    }))
    chunked_matrix <- lapply(1:nrow(batch_idx), function(idx) {
      v <- batch_idx[idx, ]
      mc_obs[v[1]:v[2], , drop=F]
    })
    temp <- applyFunc(chunked_matrix, function(submat) {
      computePosteriorProbMatrix(paramsList, submat)
    })
    if (verbose) cat("Consolidating chunked matrix\n")
    post_mat <- do.call(rbind, temp)
  }

  return(post_mat)
}


weightedMclust <- function(data, weights,
                           init_data=NULL, init_idx=NULL,
                           G=NULL, modelNames=NULL,
                           verbose=F,
                           ...) {

  if (is.null(dim(data))) data <- matrix(data, ncol=1)

  if (is.null(init_data)) {
    if (is.null(init_idx)) {
      init_data <- data
    } else {
      init_data <- data[init_idx, , drop=F]
    }
  }

  hc_init <- mclust::hc(data = init_data,
                modelName = mclust::mclust.options("hcModelName"),
                use = mclust::mclust.options("hcUse"))

  ## Prepare the for loop
  if (is.null(G)) G <- 1:10
  if (is.null(modelNames)) {
      if (ncol(data) == 1) {
        if (G == 1) {
          modelNames <- "X"
        } else {
          modelNames <- c("E", "V")
        }
      } else {
        modelNames <- mclust::mclust.options("emModelNames")
      }
  }


  ## For provided G, fit the baseline GMM and then tune with weights
  total_weights  <- sum(weights)
  d <- ncol(data)
  if (verbose) {
    cat("......Cluster Size:", nrow(init_data), "\n")
    cat("......Observations:", nrow(data), "\n")
    cat("......Total Weights:", total_weights, "\n")
  }

  grid <- expand.grid(G = G, mn = modelNames, stringsAsFactors = FALSE)

  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid$G[i]
    mn <- grid$mn[i]

    ## Don't allow models where number of params is > # of samples
    mcl_params <- mclust::nMclustParams(mn, d, g)
    if (verbose) {
      cat("........Model:", mn, "| G:", g, "| Params:", mcl_params, "\n")
    }

    if (mcl_params > total_weights) {
      out <- list(bic=-Inf)
      if (verbose) {
        cat(paste("..........Skipping; insufficient samples\n"))
      }
      attributes(out) <- list(returnCode=-342)
      return(out)
    }

    mcl <- mclust::Mclust(init_data,
                      G=g, modelNames=mn,
                      initialization=list(hcPairs=hc_init),
                      verbose=F, ...)

    if (is.null(mcl$parameters)) {
      cat(paste("..........Skipping; Mclust did not fit properly\n"))
      mcl$bic <- -Inf
      return(mcl)
    }

    ## For G = 1, Mclust() reports its single-component model names ("X",
    ## "XII", "XXI", "XXX"), which estep()/mstep() don't accept; use the
    ## equivalent "V" names for both the E-step below and me.weighted().
    if (g == 1) {
      mn <- gsub("X", "V", mcl$modelName)
      mcl$modelName <- mn
    }

    ## Get the Z matrix for data from mcl
    mcl$data <- data
    mcl$z <- mclust::estep(data, mn, mcl$parameters)$z
    mcl$z <- mcl$z + .Machine$double.eps^2
    mcl$z <- mcl$z / rowSums(mcl$z)

    ## The weighted EM can break down (e.g. a component's covariance becomes
    ## singular and its parameters turn NA); skip that candidate, as above.
    fit <- tryCatch(do.call(mclust::me.weighted, c(list(weights=weights), mcl)),
                    error = function(e) NULL)
    if (is.null(fit) || anyNA(unlist(fit$parameters[c("pro", "mean", "variance")]))) {
      if (verbose) cat("..........Skipping; weighted EM did not converge\n")
      out <- list(bic=-Inf)
      attributes(out) <- list(returnCode=-1)
      return(out)
    }

    ## me.weighted() returns a log-likelihood but no BIC, which the model
    ## selection below needs. Use the total weight as the effective sample
    ## size, as the parameter-count check above does.
    fit$bic <- if (is.finite(fit$loglik)) {
      2 * fit$loglik - mcl_params * log(total_weights)
    } else {
      -Inf
    }
    fit
  })

  maxBIC <- -Inf
  model <- NULL
  for (mcl in res) {
    if (is.null(attributes(mcl)$returnCode) || attributes(mcl)$returnCode < 0) next

    if (mcl$bic > maxBIC) model <- mcl
    maxBIC <- max(maxBIC, mcl$bic)
  }

  model
}

reweightedMclust <- function(data, weights, reweights,
                             init_data=NULL, init_idx=NULL,
                             G=NULL, modelNames=NULL, ...) {

  if (is.null(dim(data))) data <- matrix(data, ncol=1)

  if (is.null(init_data)) {
    if (is.null(init_idx)) {
      init_data <- data
    } else {
      init_data <- data[init_idx, , drop=F]
    }
  }

  hc_init <- mclust::hc(data = init_data,
                modelName = mclust::mclust.options("hcModelName"),
                use = mclust::mclust.options("hcUse"))

  ## Prepare the for loop
  if (is.null(G)) G <- 1:10
  if (is.null(modelNames)) {
    if (ncol(data) == 1) {
      if (G == 1) {
        modelNames <- "X"
      } else {
        modelNames <- c("E", "V")
      }
    } else {
      modelNames <- mclust::mclust.options("emModelNames")
    }
  }

  ## For provided G, fit the baseline GMM and then tune with weights
  res <- expand.grid(G=G, mn=modelNames) %>%
    data.frame() %>%
    apply(1, function (id_vec) {
      g <- as.numeric(id_vec["G"])
      mn <- id_vec["mn"]

      ## Don't allow models where number of params is > # of samples
      mcl_params <- mclust::nMclustParams(mn, ncol(data), g)
      if (mcl_params > sum(weights)) {
        out <- list(bic=-Inf)
        cat(paste("\tSkipping", g, mn, "insufficient samples",
                  round(sum(weights), 4), "for params", mcl_params), "\n")
        attributes(out) <- list(returnCode=-342)
        return(out)
      }

      mcl <- mclust::Mclust(init_data,
                            G=g, modelNames=mn,
                            initialization=list(hcPairs=hc_init),
                            verbose=F, ...)

      ## Model fails to fit, automatically fail
      if (is.null(mcl$parameters)) {
        mcl$bic <- -Inf
        return(mcl)
      }

      ## Model names don't align with estep
      mn_old <- mn
      if (g == 1) {
        mn <- gsub("X", "V", mcl$modelName)
      }
      # print(paste(g, mn_old, mn, mcl$modelName))

      ## Get the Z matrix for data from mcl
      mcl$data <- data
      mcl$z <- mclust::estep(data, mn, mcl$parameters)$z
      mcl$z <- mcl$z + .Machine$double.eps^2
      mcl$z <- mcl$z / rowSums(mcl$z)

      do.call(mclust::me.weighted, c(list(weights=reweights), mcl))
    })

  maxBIC <- -Inf
  model <- NULL
  for (mcl in res) {
    if (is.null(attributes(mcl)$returnCode) || attributes(mcl)$returnCode < 0) next

    if (mcl$bic > maxBIC) model <- mcl
    maxBIC <- max(maxBIC, mcl$bic)
  }

  model
}
