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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computeMonteCarloPmc(params, 1e5, verbose=TRUE)
#'
#' @return Monte Carlo estimate of \eqn{P_{\rm mc}}
#' @export
computeMonteCarloPmc <- function(paramsList, mcSamples=1e5, batchSize=mcSamples, numCores=1, verbose=F) {
  K <- length(paramsList)
  if (K == 1) return(0)

  probs <- lapply(paramsList, function(x) x$prob)
  total_prob <- do.call(sum, probs)

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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' computeMonteCarloDeltaPmcMatrix(params, 1e5, verbose=TRUE)
#'
#' @return \eqn{K \times K} matrix with each pair of clusters' \eqn{\Delta P_{\rm{mc}}} value.
#' @export
computeMonteCarloDeltaPmcMatrix <- function(paramsList, mcSamples=1e6, batchSize=mcSamples, numCores=1, verbose=F) {
  K <- length(paramsList)

  post_mat <- posteriorMatrixMCPmc(paramsList, mcSamples, batchSize, numCores, verbose)

  if (verbose) cat("Evaluating Integral\n")
  output <- crossprod(post_mat) / mcSamples
  diag(output) <- 0

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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
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
  total_prob <- do.call(sum, probs)

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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
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
computePairwisePmcMatrixOLD <- function(paramsList, mc=T, ...) {
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
        ppmc <- computePmc(list(params_i, params_j), ...)$integral
      }
      output[i, j] <- ppmc
      output[j, i] <- ppmc
    }
  }
  return(output)
}

#' Decompose cluster solutions
#'
#' @description TODO: FILL ME IN
#'
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#'
#' @return List of parameters with redundant components consolidated inso a single component.
#'
#' @export
decomposeParams <- function(paramsList) {
  params_decomp <- lapply(paramsList, function(par) {
    lapply(1:length(par$prob), function(idx){
      list(
        prob=par$prob[idx],
        mean=par$mean[, idx, drop=F],
        var=par$var[, , idx, drop=F],
        class=paste0(par$class, "_", idx)
      )
    })
  })

  do.call(c, params_decomp)
}

