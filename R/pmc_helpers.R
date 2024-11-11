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
    return(Reduce(rbind, output))
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

  Reduce(rbind, output)
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


computePosteriorProbMatrix <- function(paramsList, data) {
  if (is.null(dim(data))) data <- matrix(data, ncol=1)
  K <- length(paramsList)
  distbn.mat <- Reduce(cbind,
                       lapply(1:K, function(j) {
                         densJ <- generateDistbnFunc(paramsList[[j]])
                         sum(paramsList[[j]]$prob) * densJ(data)
                       }))
  t(apply(distbn.mat, 1, function(x) x / sum(x)))
}


posteriorMatrixMCPmc <- function(paramsList, mcSamples, batchSize, numCores, verbose) {
  num_batches <- ceiling(mcSamples / batchSize)
  applyFunc <- lapply
  if (numCores > 1) applyFunc <- function(x, func) parallel::mclapply(x, func, mc.cores=numCores)

  mc_obs <- sampleMixture(paramsList, mcSamples)
  if (is.null(dim(mc_obs))) mc_obs <- matrix(mc_obs, ncol=1)

  if (verbose) cat("Computing MC Posterior Matrix\n")
  if (num_batches == 1) {
    post_mat <- computePosteriorProbMatrix(paramsList, mc_obs)
  } else {
    post_mat <- Reduce(rbind, applyFunc(1:num_batches, function(idx) {
      lb <- 1 + batchSize * (idx - 1)
      ub <- min(batchSize * idx, mcSamples)
      computePosteriorProbMatrix(paramsList, mc_obs[lb:ub, , drop=F])
    }))
  }

  return(post_mat)
}
