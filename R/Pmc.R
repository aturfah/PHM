#' Construct \eqn{P_{\rm{mc}}} parameter list from `Mclust` output
#'
#' @description
#' Take the output of `Mclust()` and format it for consumption by `compute*Pmc` and [PHM()] functions
#'
#' @details
#' This function takes the parameters object from the output of `Mclust()` and
#' transforms it to define a GMM for use in the `compute*Pmc` and `PHM` functions. The output
#' is a list of lists, where each sublist corresponds ot the parameters of a single
#' cluster distribution (default is a single Gaussian component, K = 1).
#'
#' If `singleElement = TRUE`, then all components will be combined into a single
#' list.
#'
#' The parameters in each sublist are
#' - `mean` should be a `DxK` matrix where each column corresponds to a component mean
#' - `var` should be a `DxDxK` array where each slice corresponds to a covariance matrix
#' - `prob` should be a `K` dimensional vector for the proportions within the parent mixture
#'
#' @param mclustObj Output of `Mclust()` function call
#' @param singleElement Boolean, whether to combine into a single list element
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' mcl <- Mclust(dat)
#' constructPmcParamsMclust(mcl)
#'
#' @returns List of lists where each sublist contains the parameters for the
#' mixture component or a single list if `singleEleemnt` is `TRUE`.
#' @export
constructPmcParamsMclust <- function(mclustObj, singleElement=F) {
  params <- mclustObj$parameters

  K <- length(params$pro)
  if (mclustObj$d == 1) {
    params$mean <- matrix(params$mean, nrow=1)
    params$variance$sigma <- array(params$variance$sigmasq, dim=c(1, 1, K))
  }

  output <- lapply(1:K, function(j) {
    list(
      prob=params$pro[j],
      mean=params$mean[, j, drop=F],
      var=params$variance$sigma[, , j, drop=F],
      class=j
    )
  })

  if (singleElement) {
    g <- function(a, b) {
      list(
        mean=cbind(a$mean, b$mean),
        prob=c(a$prob, b$prob),
        var=abind::abind(a$var, b$var, along=3),
        class=NULL
      )
    }
    output <- Reduce(g, output)
  }
  return(output)
}

#' Construct \eqn{P_{\rm{mc}}} parameter list from a partition
#'
#' @description
#' Estimates the mixture model density for a partition of the data using `Mclust`
#'
#' @details
#' Performs a naive density estimation, fitting a GMM to each partition separately.
#' Only observations in a cluster are considered to fit the GMM.
#'
#' See [constructPmcParamsMclust()] for a description of the output.
#'
#' @param partition Vector of labels for observations
#' @param data Numeric matrix
#' @param ... Parameters passed to [mclust::Mclust()]
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#'
#' @returns List of lists where each sublist contains the
#' proportion, mean, covariance matrix estimates
#' for each cluster.
#' @export
constructPmcParamsPartition <- function(partition, data, ...) {
  if (is.null(dim(data))) data <- matrix(data, ncol=1)

  label_ids <- unique(partition)
  lapply(label_ids, function(k) {
    dat <- data[which(partition == k), , drop=F]
    mcl <- mclust::Mclust(dat, verbose=F, ...)
    pars <- constructPmcParamsMclust(mcl, T)

    pars$class <- k
    pars$prob <- pars$prob * mean(partition == k)
    return(pars)
  })
}


#' Construct \eqn{P_{\rm{mc}}} parameter list from a partition with weights
#'
#' @description
#' Estimates the weighted mixture model density for a partition of the data using `Mclust`.
#'
#' @details
#' This procedure attempts to account for clustering uncertainty when estimating the cluster densities.
#' For a given observation \eqn{x_i} and cluster \eqn{C_j} we estimate its cluster-specific weight based on an observation-cluster distance \eqn{d(x_i, C_j)}.
#' By default we take the observation-cluster distance to be the smallest Euclidean distance to any member in that cluster, and compute the weight for cluster \eqn{j} as:
#' \deqn{w_{ij} = \frac{e^{-d(x_i, C_j)}}{\sum_{k=1}^K e^{-d(x_i, C_k)}}}
#' [constructPmcParamsPartition()] is a special case of this procedure with weights of 1 if an observation is in a cluster and 0 otherwise.
#' The weights are then passed to a weighted EM procedure to estimate the cluster-specific density via GMM.
#' 
#' See [constructPmcParamsMclust()] for a description of the output.
#'
#' @param partition Vector of labels for observations
#' @param data Numeric matrix for data (\eqn{N \times p})
#' @param weights Optional precomputed weight matrix (\eqn{N \times K})
#' @param threshold Threshold past which to not include weights (for computational efficiency)
#' @param linkFunc Function to use to combine all distances within a cluster for weight calculation. Default is \code{min}
#' @param verbose Whether to print out debug messages. Default is \code{FALSE}
#' @param ... Parameters passed to [mclust::Mclust()]
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params_w <- constructPmcParamsWeightedPartition(partition, dat, G=1:5)
#'
#' @returns List of lists where each sublist contains the
#' proportion, mean, covariance matrix estimates
#' for each cluster.
#' @export
constructPmcParamsWeightedPartition <- function(partition, data, weights=NULL, threshold=1e-4, linkFunc=min, verbose=F, ...) {
  if (!is.factor(partition)) partition <- as.factor(partition)
  K <- length(unique(partition))

  if (is.null(dim(data)))
    data <- matrix(data, ncol=1)
  else if (!is.matrix(data)) {
    stop("data must be a matrix")
  }

  if (nrow(data) != length(partition))
    stop("data and partition must have same number of observations")

  ## Validate weights
  if (is.null(weights)) {
    # dist_mat <- as.matrix(dist(data, "euc"))
    dist_to_clust <- sapply(unique(partition), function(k) {
      clust_idx <- which(partition == k)
      clust_mat <- data[clust_idx, ]
      sapply(1:nrow(data), function(idx) {
        clust_dist <- colSums((t(clust_mat) - data[idx, ])^2)
        linkFunc(clust_dist[which(clust_dist > 0)])
      })
    })
    mindist <- apply(dist_to_clust, 1, min)
    weights <- exp(-1 * (dist_to_clust - mindist))
    weights <- weights / rowSums(weights)
  } else if (is.matrix(weights)) {
    if (ncol(weights) != K)
      stop("weights must have 1 column per cluster")
    if (nrow(weights) != nrow(data))
      stop("weights must have 1 row per observation")
  } else {
    stop("Weights must be provided as an {Observations}x{Clusters} matrix")
  }

  label_ids <- unique(partition)
  lapply(label_ids, function(k) {
    idx <- which(label_ids == k)
    dat <- data[which(partition == k), , drop=F]
    w <- weights[, idx]

    ## Filter ignore observations with w < threshold
    valid <- which(w >= threshold)
    data_valid <- data[valid, , drop=F]
    w <- w[valid]

    if (verbose) {
      cat("Estimating Cluster", which(label_ids == k), 
          "of", length(label_ids), "Nk =", nrow(dat), "\n")
      cat("\tWeights:", round(sum(w), 4), "\n")
    }

    wmcl <- weightedMclust(data_valid, weights = w, init_data = dat, ...)
    pars <- constructPmcParamsMclust(wmcl, T)
    pars$class <- k
    pars$prob <- pars$prob * nrow(dat) / nrow(data)
    return(pars)
  })
}


#' Monte Carlo \eqn{P_{\rm{mc}}} computation
#'
#' @description
#' Compute Monte Carlo estimate of \eqn{P_{\rm mc}} for a given cluster configuration based on estimated GMM densities.
#'
#' @details
#' \eqn{P_{\rm{mc}}} can be difficult to evaluate as standard cubature methods tend to perform poorly in higher dimensions.
#' We can approximate it for a \eqn{K}-cluster configuration using a Monte Carlo integral of the form
#' \deqn{\hat P_{{\rm mc}} = \frac{1}{M} \sum_{i=1}^{M} \sum_{j=1}^K 2 \left(1 - \pi_j(x_i) \right) \, \pi_j(x_i)}
#' Where the \eqn{M} observations are sampled from the overall data density \eqn{P(x)}
#'
#' @param paramsList List containing lists with each component GMM parameters. See [constructPmcParamsMclust] for format of components.
#' @param mcSamples Numeric for number of MC samples to use to approximate the integral.
#' @param batchSize Numeric for the observations to assign to each core. Helps with memory concerns. Default `mcSamples`.
#' @param numCores Number of cores to use in parallel::mclapply call. Default is 1.
#' @param verbose Boolean whether to print output messages
#'
#' @examples 
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computeMonteCarloPmc(params, 1e5, verbose=T)
#' 
#' @return Monte Carlo estimate of \eqn{P_{\rm mc}}
#' @export
computeMonteCarloPmc <- function(paramsList, mcSamples=1e5, batchSize=mcSamples, numCores=1, verbose=F) {
  K <- length(paramsList)
  if (K == 1) return(0)

  probs <- lapply(paramsList, function(x) x$prob)
  total_prob <- Reduce(sum, probs)

  if (total_prob != 1) {
    warning("Probabilities in paramsList do not sum to 1. Re-scaling...")
    for (idx in 1:K) {
      paramsList[[idx]]$prob <- paramsList[[idx]]$prob / total_prob
    }
  }

  post_mat <- posteriorMatrixMCPmc(paramsList, mcSamples, batchSize, numCores, verbose)

  return(sum(post_mat * (1 - post_mat)) / mcSamples)
}


#' Monte Carlo \eqn{\Delta P_{\rm {mc}}} Matrix computation
#'
#' @description
#' Compute the \eqn{\Delta P_{\rm mc}} matrix for a set of clusters, where the \eqn{ij^{th}} element is \eqn{\Delta P_{\rm mc}^{(i, j)}}.
#'
#' @details
#' Each step of the PHM algorithm reduces the overall \eqn{P_{\rm mc}} by the \eqn{\Delta P_{\rm mc}} value of the merged clusters.
#' For each pair of clusters \eqn{j, k} we estimate their \eqn{\Delta P_{\rm mc}^{(j, k)}} value.
#' \deqn{\Delta \hat P_{\rm mc}^{(j, k)} = \frac{1}{M} \sum_{i=1}^M \pi_j(x_i) \, \pi_k(x_i) }
#' Where the \eqn{M} observations are sampled from the overall data density \eqn{P(x)}.
#' Note that \eqn{\Delta \hat P_{\rm mc}^{(j, k)} = \Delta \hat P_{\rm mc}^{(k, j)}}
#' 
#' @inheritParams computeMonteCarloPmc
#'
#' @examples 
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computeMonteCarloDeltaPmcMatrix(params, 1e5, verbose=T)
#' 
#' @return \eqn{K \times K} matrix with each pair of clusters' \eqn{\Delta P_{\rm{mc}}} value.
#' @export
computeMonteCarloDeltaPmcMatrix <- function(paramsList, mcSamples=1e6, batchSize=mcSamples, numCores=1, verbose=F) {
  K <- length(paramsList)

  output <- matrix(0, K, K)
  post_mat <- posteriorMatrixMCPmc(paramsList, mcSamples, batchSize, numCores, verbose)

  if (verbose) cat("Evaluating Integral\n")
  for (i in 1:(K-1)) {
    for (j in (i+1):K) {
      post_i <- post_mat[, i]
      post_j <- post_mat[, j]

      output[i, j] <- sum(post_i * post_j) / mcSamples
      output[j, i] <- output[i, j]
    }
  }

  output
}


#' Cubature \eqn{P_{\rm{mc}}} computation
#'
#' Compute Pmc for a given Gaussian mixture distribution based on `cubature` package
#'
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param integralControl List specifying arguments to pass to [cubature::cubintegrate()]. See details.
#' 
#' @details
#' For a given cluster configuration, the overall misclassification probability \eqn{P_{\rm mc}} can be evaluated as
#' \deqn{P_{\rm mc} = \sum_{j=1}^K \int 2 \left(1 - \pi_j(x)\right) \, \pi_j(x) P(x) dx}
#' 
#' This integral is implemented using the `cubature` function. the `integralControl` variable accepts arguments to the [cubature::cubintegrate()] function
#' The defaults for this function are:
#' \itemize{
#'  \item \code{method}: Which integration method is to be used. Default is \code{hcubature}
#'  \item \code{lowerLimit} and \code{upperLimit}: The bounds of integration. Default is \eqn{\pm \infty}.
#'  \item \code{maxEval}: Sets the maximum number of integral evaluations. Default is 1e6
#'  \item \code{relTol}: Sets the convergence tolerance for the integration. Default is 1e-5
#' }
#'
#' @examples 
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computePmc(params)
#' 
#' @return Output from the `cubature::cubintegrate()` function.
#' @export
computePmc <- function(paramsList, integralControl=list()) {

  ## Default Integral Parameters
  intCont <- list(
    method="hcubature",
    lowerLimit=-Inf,
    upperLimit=Inf,
    maxEval=1e6,
    relTol=1e-5,
    nVec=1024L
  )
  intCont[names(integralControl)] <- integralControl

  ## Get the parameters
  K <- length(paramsList)
  D <- nrow(paramsList[[1]]$mean)

  ## Pmc trivially 0 for a single cluster
  if (K == 1) return(list(
    integral=0,
    error=NULL,
    neval=NULL,
    returnCode=0
  ))

  ## Input validation: Make sure component probabilities all scale to 1
  probs <- lapply(paramsList, function(x) x$prob)
  total_prob <- Reduce(sum, probs)

  if (total_prob != 1) {
    warning("Probabilities in paramsList do not sum to 1. Re-scaling...")
    for (idx in 1:K) {
      paramsList[[idx]]$prob <- paramsList[[idx]]$prob / total_prob
    }
  }

  ## Input validation: Make sure dimensions line up for means/covariances
  for (j in 1:K) {
    M <- length(paramsList[[j]]$prob)
    ## Means checking
    stopifnot(D == nrow(paramsList[[j]]$mean))
    stopifnot(M == ncol(paramsList[[j]]$mean))

    ## Covariance Matrix checking
    stopifnot(D == dim(paramsList[[j]]$var[1]))
    stopifnot(D == dim(paramsList[[j]]$var[2]))
    stopifnot(M == dim(paramsList[[j]]$var[3]))
  }
  rm(j, M)

  ## Functions to perform integral
  pmcIntegralFunc <- function(x, cubatureFunc=T) {
    if (is.null(dim(x))) x <- matrix(x, nrow=1)

    if (cubatureFunc == T) x <- t(x)

    distbn.mat <- sapply(1:K, function(j) {
      comp_prob <- sum(paramsList[[j]]$prob)
      comp_prob * generateDistbnFunc(paramsList[[j]])(x)
    })
    fX <- apply(distbn.mat, 1, sum)

    post.mat <- t(apply(distbn.mat, 1, function(v) {
      denom <- sum(v)
      if (denom == 0) denom <- 1 ## If this is 0 then fX is 0 so it doesn't matter

      v / denom
    }))

    res <- sapply(1:K, function(j) {
      postJ <- post.mat[, j]
      postJ * (1 - postJ) * fX
    })

    if (is.null(dim(res))) return(sum(res))

    output <- apply(res, 1, sum)
    if (cubatureFunc == T) return(matrix(output, ncol=nrow(x)))

    output
  }

  res <- cubature::cubintegrate(
    f=pmcIntegralFunc,
    lower=rep(intCont$lowerLimit, D),
    upper=rep(intCont$upperLimit, D),
    nVec=intCont$nVec,
    relTol=intCont$relTol,
    maxEval=intCont$maxEval,
    method=intCont$method
  )

  return(res)
}


#'  \eqn{\Delta P_{\rm {mc}}} Matrix computation
#'
#' Compute the \eqn{\Delta P_{\rm mc}} matrix for a set of clusters, where the \eqn{ij^{th}} element is \eqn{\Delta P_{\rm mc}^{(i, j)}}.
#'
#' @inheritParams computePmc
#'
#' @details 
#' #' Each step of the PHM algorithm reduces the overall \eqn{P_{\rm mc}} by the \eqn{\Delta P_{\rm mc}} value of the merged clusters.
#' For each pair of clusters \eqn{j, k}, \eqn{\Delta P_{\rm mc}} is 
#' \deqn{\Delta P_{\rm mc}^{(j, k)} = \int \pi_j(x) \, \pi_k(x) P(x) dx}
#' Where the relationship \eqn{P_{\rm mc} = \sum_{i < j} 2\Delta P_{\rm mc}^{(i, j)}}
#' See [computePmc] for description of `integralControl` parameters.
#' 
#' @return \eqn{K \times K} matrix with each pair of clusters' contribution to \eqn{P_{\rm{mc}}}
#' 
#' @examples 
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computeDeltaPmcMatrix(params)
#' 
#' @export
computeDeltaPmcMatrix <- function(paramsList, integralControl=list()) {
  K <- length(paramsList)
  D <- nrow(paramsList[[1]]$mean)
  output <- matrix(0, nrow=K, ncol=K)

  intCont <- list(
    method="cuhre",
    lowerLimit=-Inf,
    upperLimit=Inf,
    maxEval=1e6,
    relTol=1e-6,
    nVec=1024L
  )
  intCont[names(integralControl)] <- integralControl


  generatePmcReductionFunc <- function(i, j) {
    post_i <- generatePosteriorProbFunc(paramsList, i)
    post_j <- generatePosteriorProbFunc(paramsList, j)

    function (x) {
      if (is.null(dim(x))) x <- matrix(x, nrow=1)
      x <- t(x) ## Cubature passes in transposed to what we expect

      fX <- sapply(1:K, function(k) {
        comp_prob <- sum(paramsList[[k]]$prob)
        comp_prob * generateDistbnFunc(paramsList[[k]])(x)
      })
      fX <- apply(fX, 1, sum)

      output <- post_i(x) * post_j(x) * fX
      return(matrix(output, ncol=nrow(x)))
    }
  }

  for (i in 1:(K-1)) {
    for (j in ((i+1):K)) {
      res <- cubature::cubintegrate(
        f=generatePmcReductionFunc(i, j),
        lower=rep(intCont$lowerLimit, D),
        upper=rep(intCont$upperLimit, D),
        nVec=intCont$nVec,
        method=intCont$method,
        maxEval=intCont$maxEval,
        relTol=intCont$relTol
      )

      output[i, j] <- res$integral
      output[j, i] <- output[i, j]
    }
  }

  output
}


#' Pairwise \eqn{P_{\rm{mc}}} Matrix computation
#'
#' @description TODO: FILL ME IN
#'
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param mc Boolean whether to compute \eqn{P_{\rm mc}} with Monte Carlo integration or cubature integration
#' @param ... Additional parameters passed to [computePmc()] or [computeMonteCarloPmc()]
#'
#' @return \eqn{K \times K} matrix with Pairwise \eqn{P_{\rm{mc}}} values for each pair of clusters
#'
#' export
computePairwisePmcMatrix <- function(paramsList, mc=T, ...) {
  K <- length(paramsList)
  output <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i <= j) next

      params_i <- paramsList[[i]]
      params_j <- paramsList[[j]]

      params_i$prob <- 0.5 * params_i$prob / sum(params_i$prob)
      params_j$prob <- 0.5 * params_j$prob / sum(params_j$prob)

      if (mc) {
        ppmc <- computeMonteCarloPmc(list(params_i, params_j), verbose=F, ...)
      } else {
        ppmc <- computePmc(list(params_i, params_j), verbose=F, ...)$integral
      }
      output[i, j] <- ppmc
      output[j, i] <- ppmc
    }
  }
  return(output)
}
