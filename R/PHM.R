#' Find the solution for PHM at a specified threshold
#'
#' @param phm_output Output of `PHM`
#' @param threshold Pmc threshold, default is 0.01
#'
#' TODO: Fill in Returns + More detailed description
#'
#' @return Cluster results based on GMM satisfying specified threshold
#' @export
thresholdPHM <- function(phm_output, threshold=0.01) {
  kappa <- length(phm_output)
  for (k in kappa:1) {
    if (phm_output[[k]]$pmc < threshold) {
      break
    }
  }

  return(phm_output[[k]])
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
#' TODO: FILL ME IN
#'
#' @return FILL ME IN
#' @export
PHM <- function(mclustObj=NULL, paramsList=NULL, partition=NULL, data=NULL,
                verbose=T,
                computePosterior=T,
                partitionModel="VVI",
                partitionMaxComponents=10,
                mc=T, ...) {

  ## Validate input parameter values
  posterior_matrix <- NULL
  data_labels <- NULL
  if (!is.null(mclustObj)) {
    if (is.null(data) & computePosterior) {
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
    if (verbose) cat("Estimating partition densities\n")
    paramsList <- constructPmcParamsPartition(partition, data,
                                           modelNames=partitionModel,
                                           G=1:partitionMaxComponents)
  } else if (!is.null(paramsList)) {
    if (is.null(data)) {
      stop("Must provide data")
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
      data_labels <- apply(posterior_matrix, 1, which.max)
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

  ## Parameters/data to store
  output <- lapply(1:K, function(k) list(
    clusters=k,
    posterior_matrixrix=if(computePosterior) {matrix(1, nrow=N)} else {NULL},
    labels=if(computePosterior) {rep(1, N)} else {NULL},
    pmc_change=NA,
    params=NULL,
    merge_components=c(-1, -1),
    pmc=0))

  pmc <- sum(delta_pmc)
  tmp_delta <- delta_pmc
  tmp_params <- paramsList

  for (idx in K:2) {
    ## Store the results
    class_labels <- sapply(tmp_params, function(x) x$class)
    output[[idx]]$labels <- class_labels[data_labels]
    output[[idx]]$pmc <- sum(tmp_delta)
    output[[idx]]$params <- tmp_params
    if (idx < K) {
      output[[idx]]$pmc_change <- output[[idx+1]]$pmc - output[[idx]]$pmc
    }
    output[[idx]]$pmc_matrix <- tmp_delta
    output[[idx]]$posterior_matrixrix <- posterior_matrix

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

    ## New row for Pairwise Pmc Matrix
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

  output
}

#' Compute Posterior Matrix based on PHM merging
#'
#' @details TODO: Fill me in
#'
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

