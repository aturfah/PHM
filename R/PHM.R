#' Find a clustering for a given \eqn{P_{\rm{mc}}} threshold
#'
#' @param phm Output of [PHM()]
#' @param threshold \eqn{P_{\rm{mc}}} threshold, default is 0.01
#'
#' @description Each step of the PHM algorithm reduces \eqn{P_{\rm{mc}}}. This gives the results from the PHM algorithm terminated when \eqn{P_{\rm{mc}}} falls below some specified threshold.
#'
#' @return Result of the PHM merging procedure terminated when the \eqn{P_{\rm{mc}}} threshold is satisfied
#' export
thresholdPHM <- function(phm, threshold=0.01) {
  kappa <- length(phm)
  for (k in kappa:1) {
    if (phm[[k]]$pmc < threshold) {
      break
    }
  }

  return(phm[[k]])
}

mergeParams <- function(par1, par2) {
  out <- list(
    prob=c(par1$prob, par2$prob),
    mean=cbind(par1$mean, par2$mean),
    var=abind::abind(par1$var, par2$var, along=3)
  )

  if (!is.null(par1$class) & !is.null(par2$class)) {
    out$class <- paste(par1$class, par2$class, sep="|")
  }

  return(out)
}


#' PHM Algorithm
#'
#' @description 
#' Implements the PHM algorithm which constructs a clustering hierarchy by successively merging clusters with the largest \eqn{\Delta P_{\rm mc}} values.
#' 
#' @param mclustObj Output from [mclust::Mclust()]
#' @param paramsList A list generated from [constructPmcParamsMclust()], [constructPmcParamsPartition()], [constructPmcParamsPartition()] providing the initial cluster parameter estimates
#' @param partition A vector providing obseration partition memberships for the initial state
#' @param data An \eqn{N \times D} matrix of observations
#' @param verbose Boolean whether to suppress debug statements
#' @param computePosterior Boolean whether to compute the Posterior matrix for the merges
#' @param partitionModel If no partition provided, the covariance structure to estimate the density for each partition using [constructPmcParamsPartition()]
#' @param partitionMaxComponents If specifying `partition`, the maximum number of components to estimate the density for each partition
#' @param mc Boolean whether to use Monte Carlo integration to evaluate the \eqn{\Delta P_{\rm{mc}}} matrix
#' @param storeDeltaPmc Boolean whether to store the \eqn{\Delta P_{\rm{mc}}} matrix at each merging step. Initial \eqn{\Delta P_{\rm{mc}}} will be stored.
#' @param storeParams Boolean whether to store the intermediate mixture parameters. Initial parameters will be stored.
#' @param ... Parameters pased to either [computeDeltaPmcMatrix()] or [computeMonteCarloDeltaPmcMatrix()] to evaluate the \eqn{\Delta P_{\rm{mc}}} matrix
#' 
#' @examples 
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=T)
#' partition <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#' params <- constructPmcParamsPartition(partition, dat, G=1:5)
#' phm <- PHM(paramsList=params, data=dat, partition=partition)
#' 
#' @return
#' A list of lists for each step of the PHM algorithm. Each sublist contains
#' \itemize{
#'  \item \code{clusters}: Number of clusters \eqn{K} at this merge
#'  \item \code{posterior_matrix}: \eqn{N\times K} matrix of posterior cluster probabilities
#'  \item \code{labels}: Partition of the observations
#'  \item \code{pmc_change}: Value of \eqn{\Delta P_{\rm mc}} leading to this value of \eqn{K}
#'  \item \code{params}: Cluster-specific densities
#'  \item \code{pmc_components}: Number of original clusters involved in this merge
#'  \item \code{pmc_accum}: Accumulated \eqn{\Delta P_{\rm mc}} in this subtree (unused)
#'  \item \code{min_merge_pmc}: Minimum value of \eqn{\Delta P_{\rm mc}} for all merges in this subtree
#'  \item \code{merge_components}: Index of components merged in this step
#'  \item \code{pmc}: Overall \eqn{P_{\rm mc}} remaining in the cluster configuration
#'  \item \code{pmc_matrix}: \eqn{\Delta P_{\rm mc}} matrix for the remaining clusters
#' }
#' 
#' @export
PHM <- function(mclustObj=NULL, 
                paramsList=NULL, partition=NULL, data=NULL,
                deltaPmc=NULL,
                verbose=T,
                computePosterior=T,
                partitionWeightedDensity=T,
                storeDeltaPmc=T,
                storeParams=T,
                partitionModel="VVI",
                partitionMaxComponents=10,
                mc=T, ...) {

  ## Validate input parameter values
  posterior_matrix <- NULL
  data_labels <- NULL
  if (!is.null(mclustObj)) {
    if (is.null(data) && computePosterior) {
      ## Data/Posterior Matrix assumed from the object
      data <- mclustObj$data
      posterior_matrix <- mclustObj$z
      data_labels <- mclustObj$classification
    }
    paramsList <- constructPmcParamsMclust(mclustObj)
  } else if (!is.null(partition)) {
    if (is.null(data)) {
      stop("Must provide data")
    }
    if (is.null(paramsList)) {
      if (verbose) cat("Estimating partition densities\n")
      partDensFunc <- if (partitionWeightedDensity) {
        constructPmcParamsWeightedPartition
      } else {
        constructPmcParamsPartition
      }
      paramsList <- partDensFunc(partition, data,
                                modelNames=partitionModel,
                                G=1:partitionMaxComponents)
    }
    part_levels <- sapply(paramsList, function(x) x$class)
    data_labels <- as.numeric(factor(partition, levels=part_levels))
  } else if (!is.null(paramsList)) {
    if (is.null(data)) {
      stop("Must provide data")
    }
    if (!is.null(partition)) {
      part_levels <- sapply(paramsList, function(x) x$class)
      data_labels <- as.numeric(factor(partition, levels=part_levels))
    }
  } else {
    stop("Must provide one of mclustObj, paramsList, or partition")
  }

  if (computePosterior & is.null(data)) {
    stop("Need data to compute posterior")
  }

  ## Prepare the Posterior Matrix and Class labels
  if (computePosterior == T) {
    if (is.null(posterior_matrix)) {
      if (verbose) cat("PHM Constructing Posterior Data Matrix\n")
      posterior_matrix <- computePosteriorProbMatrix(paramsList, data)
      if (is.null(partition)) {
        data_labels <- apply(posterior_matrix, 1, which.max)
      }
    }
  }
  N <- nrow(data)
  K <- length(paramsList)

  ## Construct Delta Pmc Matrix
  if (verbose) cat("Constructing Delta Pmc Matrix\n")
  delta_pmc <- if (mc) {
    computeMonteCarloDeltaPmcMatrix(paramsList, verbose=verbose, ...)
  } else {
    computeDeltaPmcMatrix(paramsList, ...)
  }
  diag(delta_pmc) <- 0
  if (verbose) cat("Complete; Running PHM algorithm\n")

  ## Parameters/data to store
  output <- lapply(1:K, function(k) list(
    clusters=k,
    posterior_matrix=if(computePosterior) {matrix(1, nrow=N)} else {NULL},
    labels=if(computePosterior) {rep(1, N)} else {NULL},
    pmc_change=NA,
    params=NULL,
    pmc_components=1,
    pmc_accum=0,
    min_merge_pmc=NA,
    merge_components=c(-1, -1),
    pmc=0))

  pmc <- sum(delta_pmc)
  tmp_delta <- delta_pmc
  tmp_params <- paramsList
  components <- rep(1, K)
  tmp_comp <- 2
  pmc_accum <- rep(0, K)
  component_map <- lapply(1:K, function(x) x)

  for (idx in K:2) {
    ## Store the results
    class_labels <- sapply(tmp_params, function(x) x$class)
    output[[idx]]$labels <- class_labels[data_labels]
    output[[idx]]$pmc <- sum(tmp_delta)
    if (storeParams || idx == K) {
      output[[idx]]$params <- tmp_params
    }
    output[[idx]]$pmc_components <- tmp_comp
    if (idx < K) {
      output[[idx]]$pmc_change <- output[[idx+1]]$pmc - output[[idx]]$pmc
      output[[idx]]$pmc_accum <- pmc_accum[i]
      pmc_accum[i] <- pmc_accum[i] + output[[idx]]$pmc_change
      output[[idx]]$min_merge_pmc <- min_merge_pmc
    }
    if (storeDeltaPmc || idx == K) {
      output[[idx]]$pmc_matrix <- tmp_delta
    }
    output[[idx]]$posterior_matrix <- posterior_matrix

    ## Identify the components to merge
    maxval <- max(tmp_delta, na.rm=T)
    cand_rows <- which(tmp_delta == maxval, arr.ind=T)

    cand_rows <- cand_rows[which(cand_rows[, 1] < cand_rows[, 2]), , drop=F]
    i <- cand_rows[1, 1]
    j <- cand_rows[1, 2]
    output[[idx]]$merge_components <- c(i, j) ## store this

    ## New row for DeltaPmc Matrix
    row_i <- tmp_delta[i, ]
    row_j <- tmp_delta[j, ]
    new_row_delta <- row_i + row_j
    new_row_delta[i] <- 0
    new_row_delta <- new_row_delta[-j]

    ## Track minimum Pmc between the components
    component_map[[i]] <- c(component_map[[i]], component_map[[j]])
    candidate_posns <- as.matrix(expand.grid(component_map[[i]], component_map[[i]]))
    candidate_posns <- candidate_posns[which(candidate_posns[, 1] - candidate_posns[, 2] != 0), ]
    min_merge_pmc <- delta_pmc[candidate_posns]
    min_merge_pmc <- min(min_merge_pmc[which(min_merge_pmc > 0)]) ## Ignore 0 values for this case
    component_map[[j]] <- NULL


    ## Keep track of the components and accumulated Pmc for height
    pmc_accum[i] <- pmc_accum[i] + pmc_accum[j]
    pmc_accum[i] <- pmc_accum[i]
    pmc_accum <- pmc_accum[-j]

    # tmp_comp <- max(components[i], components[j])
    components[i] <- components[i] + components[j]
    tmp_comp <- components[i]
    components <- components[-j]

    ## Construct new Merge Parameter
    merged_params <- mergeParams(tmp_params[[i]], tmp_params[[j]])

    ## Shrink matrix my removing j row/column
    tmp_delta <- tmp_delta[-j, , drop=F]
    tmp_delta <- tmp_delta[, -j, drop=F]

    ## Replace i terms with i+j
    tmp_delta[i, ] <- new_row_delta
    tmp_delta[, i] <- new_row_delta

    ## Update the parameters
    tmp_params[[i]] <- merged_params
    tmp_params[[j]] <- NULL

    ## Update the cluster labels; make sure always 1:K
    if (computePosterior) {
      posterior_matrix[, i] <- posterior_matrix[, i] + posterior_matrix[, j]
      posterior_matrix <- posterior_matrix[, -j, drop=F]
      data_labels <- apply(posterior_matrix, 1, which.max)
    }

  }

  output[[1]]$pmc_change <- output[[2]]$pmc
  output[[1]]$pmc_components <- tmp_comp
  output[[1]]$pmc_accum <- pmc - output[[1]]$pmc_change
  output[[1]]$min_merge_pmc <- min_merge_pmc

  output
}

#' Compute Posterior Matrix based on PHM merging
#'
#' @details TODO: Fill me in
#'
#' @param phm Output from [PHM()]
#' @param data \eqn{N \times D} matrix of observations
#' @param initK Number of clusters from which to start computing the posterior matrix
#' 
#' @return List of the same structure as from [PHM()] with `posterior_matrix` and `labels` fields calculated for the specified elements.
#' @export
addPosteriorMatrix <- function(phm, data, initK=length(phm)) {
  posterior_matrix <- computePosteriorProbMatrix(phm[[initK]]$params, data)
  data_labels <- apply(posterior_matrix, 1, which.max)

  for (j in initK:2) {
    class_labels <- sapply(phm[[j]]$params, function(x) x$class)
    phm[[j]]$posterior_matrix <- posterior_matrix
    phm[[j]]$labels <- class_labels[data_labels]

    i <- phm[[j]]$merge_components[1]
    j <- phm[[j]]$merge_components[2]

    posterior_matrix[, i] <- posterior_matrix[, i] + posterior_matrix[, j]
    posterior_matrix <- posterior_matrix[, -j, drop=F]
    data_labels <- apply(posterior_matrix, 1, which.max)
  }

  phm[[1]]$posterior_matrix <- posterior_matrix
  phm[[1]]$labels <- data_labels

  return(phm)
}

#' Recover the posterior matrix based on PHM merging
#'
#' @details TODO: Fill me in
#'
#' @param phm Output from [PHM()]
#' @param K Number of clusters for which to compute the deltaPmc Matrix
#' 
#' @return TODO
#' @export
recoverDeltaPmcMatrix <- function(phm, K) {
  stopifnot(K < length(phm))
  
  ## Get the closest value in PHM from which we can construct K
  valid_idx <- which(sapply(phm, function(x) !is.null(x$pmc_matrix)))
  min_idx <- min(valid_idx[which(valid_idx > K)])
  
  if (!is.finite(min_idx)) stop("No valid deltaPmc matrix identified.")
  

  tmp_delta <- phm[[min_idx]]$pmc_matrix

  for (idx in min_idx:K) {
    if (idx %% 50 == 0) print(dim(tmp_delta))
    i <- phm[[idx]]$merge_components[1]
    j <- phm[[idx]]$merge_components[2]
    
    row_i <- tmp_delta[i, ]
    row_j <- tmp_delta[j, ]
    new_row_delta <- row_i + row_j
    new_row_delta[i] <- 0
    new_row_delta <- new_row_delta[-j]
    
    tmp_delta <- tmp_delta[-j, , drop=F]
    tmp_delta <- tmp_delta[, -j, drop=F]
    
    ## Replace i terms with i+j
    tmp_delta[i, ] <- new_row_delta
    tmp_delta[, i] <- new_row_delta
  }
  
  tmp_delta
}

#' Recover the parameters based on the PHM merging
#'
#' @details TODO: Fill me in
#'
#' @param phm Output from [PHM()]
#' @param K Number of clusters for which to compute the deltaPmc Matrix
#' @param paramsToKeep Parameters of interest
#' 
#' @return TODO
#' @export
recoverPHMParams <- function(phm, K, paramsToKeep=c("prob", "mean", "var", "class")) {
  stopifnot(K < length(phm))
  
  ## Get the closest value in PHM from which we can construct K
  valid_idx <- which(sapply(phm, function(x) !is.null(x$params)))
  min_idx <- min(valid_idx[which(valid_idx > K)])
  
  if (!is.finite(min_idx)) stop("No valid parameters identified.")
  
  ## Clear out names that we won't be using
  tmp_params <- phm[[min_idx]]$params
  tmp_params <- lapply(tmp_params, function(l) {
    for (par_name in names(l)) {
      if (!(par_name %in% paramsToKeep)) {
        l[[par_name]] <- NULL
      }
    }
    l
  })

  for (idx in min_idx:K) {
    if (idx %% 100 == 0) print(length(tmp_params))
    i <- phm[[idx]]$merge_components[1]
    j <- phm[[idx]]$merge_components[2]

    tmp_params[[i]] <- mergeParams(tmp_params[[i]], tmp_params[[j]])
    tmp_params[[j]] <- NULL
  }
  
  ## Once again clear out names introduced by mergeParams
  tmp_params <- lapply(tmp_params, function(l) {
    for (par_name in names(l)) {
      if (!(par_name %in% paramsToKeep)) {
        l[[par_name]] <- NULL
      }
    }
    l
  })

  tmp_params
}


#' Recover the posterior matrix based on the PHM merging
#'
#' @details TODO: Fill me in
#'
#' @param phm Output from [PHM()]
#' @param K Number of clusters for which to compute the deltaPmc Matrix
#' @param data Matrix of observations for which to compute the posterior matrix
#' @param computePosterior Whether to recompute the posterior matrix if not available
#' 
#' @return TODO
#' @export
recoverPosterior <- function(phm, K, data=NULL, computePosterior=F) {
  stopifnot(K < length(phm))

  ## Get the closest value in PHM from which we can construct K
  valid_idx <- which(sapply(phm, function(x) !is.null(x$posterior_matrix)))
  min_idx <- min(valid_idx[which(valid_idx > K)])

  if (!is.finite(min_idx)) {
    if (computePosterior) {
      ## Get the combined parameter and just recompute the posterior matrix
      params <- phm[[min_idx]]$params
      if (is.null(params)) {
        warning("No parameters computed, reconstructing from PHM process to compute posterior")
        params <- recoverPHMParams(phm, K)
        if (is.null(data)) {
          stop("To compute matrix, must provide data")
        }
        posterior_matrix <- computePosteriorProbMatrix(params, data)
        return(posterior_matrix)
      }
    } else {
      stop("No valid posterior matrix identified. Rerun with computePosterior=T to compute the posterior matrix.")
    }
  } else {
    ## If we have the posterior matrix, grab it
    posterior_matrix <- phm[[min_idx]]$posterior_matrix
  }

  for (idx in min_idx:K) {
    if (idx %% 100 == 0) print(length(tmp_params))
    i <- phm[[idx]]$merge_components[1]
    j <- phm[[idx]]$merge_components[2]
    
    posterior_matrix[, i] <- posterior_matrix[, i] + posterior_matrix[, j]
    posterior_matrix <- posterior_matrix[, -j, drop=F]
  }

  data_labels <- apply(posterior_matrix, 1, which.max)

  list(
    matrix=posterior_matrix,
    labels=data_labels
  )
}
