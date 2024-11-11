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
#' # dat <- c(rnorm(100), rnorm(100, 3))
#' # mcl <- Mclust(dat)
#' # constructPmcParamsMclust(mcl)
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

  } else {
    return(output)
  }
}

#' Construct \eqn{P_{\rm{mc}}} parameter list from a partition
#'
#' @description
#' Estimates the mixture model density for a partition of the data using `Mclust`
#'
#' @details
#' TODO: Fill me in.
#' 
#' See [constructPmcParamsMclust()] for a description of the output.
#'
#' @param partition Vector of labels for observations
#' @param data Numeric matrix
#'
#' @examples
#' # dat <- c(rnorm(100), rnorm(100, 3))
#' # partition <- c(rep(1, 100), rep(2, 100))
#' # constructPmcParamsPartition(partition, dat)
#'
#' @returns List of lists where each sublist contains the GMM density estimates
#' for each group in the partition.
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


#' Monte Carlo \eqn{P_{\rm{mc}}} computation
#'
#' @description 
#' Compute Monte Carlo estimate of Pmc associated with given cluster parameters based on Gaussian cluster assumption.
#' 
#' @details TODO: Fill me in.
#' 
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param mcSamples Numeric for number of MC samples to use to approximate the integral. Default 1e5
#' @param batchSize Numeric for the observations to assign to each core. Helps with memory concerns. Default 1e3.
#' @param numCores Number of cores to use in parallel::mclapply call. Default is `floor(detectCores() / 2)`.
#' @param verbose Boolean whether to print output messages
#'
#' @return Value of Pmc from the Monte Carlo integral
#' @export
computeMonteCarloPmc <- function(paramsList, mcSamples=1e5, batchSize=mcSamples, numCores=1, verbose=T) {
  K <- length(paramsList)
  if (K == 1) return(0)

  post_mat <- posteriorMatrixMCPmc(paramsList, mcSamples, batchSize, numCores, verbose)

  return(sum(post_mat * (1 - post_mat)) / mcSamples)
}


#' Monte Carlo $\Delta P_{\rm {mc}}$ Matrix computation
#'
#' Compute the merging Pmc reduction matrix estimated via Monte Carlo integration
#'
#' See [computeMonteCarloPmc()] for a breakdown of the parameters
#' @export
computeMonteCarloDeltaPmcMatrix <- function(paramsList, mcSamples=1e6, batchSize=mcSamples, numCores=1, verbose=T) {
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
#' @param integralControl Specifies arguments for integration methods. See details.
#'
#' @details
#' TODO: Clear this up
#' For `pcubature` the control variables are `method` for the integration
#' For `pcubature` the control variables are `lowerLimit` and `upperLimit` for the integral (defaults to \pm Inf)
#' For `pcubature` the control variables are `maxEval` which sets the number of integral evaluations (defaults to 1e6)
#' For `pcubature` the control variables are `relTol` which sets the maximum tolerance (defaults to 1e-5)
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


#'  $\Delta P_{\rm {mc}}$ Matrix computation
#'
#' Compute the merging Pmc reduction matrix calculated using `cubature` package
#'
#' See [computePmc()] for a breakdown of the parameters
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
      res <- cubintegrate(
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
