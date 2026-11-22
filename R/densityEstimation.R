#' Global (Ensemble) Density Estimation
#'
#' @description
#' Fits `M` Gaussian mixture models, each to a (sub)sample of `X`, and pools
#' their components into one mixture density, returned as a component list
#' ready for [PHMv2()]. Each replicate is fit by [fitDensityReplicate()], the
#' same function the PHM ensemble pipeline uses (see [fitPHMEnsemble()]), so
#' the same `M`, subsample size and seeds give the same components here as
#' in [fitPHMEnsemble()] or the scripts in `inst/ensemble/`.
#'
#' @details
#' Each replicate draws `round(rho * nrow(X))` rows (or `subsampleSize`
#' rows), without replacement unless `replace = TRUE`, and selects a model
#' and number of components by BIC. If the best `G` is the largest value
#' tried, the range is widened automatically (see `GStep`/`GMax`, passed
#' through `...`). The `M` replicates are pooled with equal weights, so
#' each contributes `1/M` of the total probability.
#'
#' With `saveDir` set, each replicate's result is saved there and reused on
#' later calls with the same `prefix`, subsample size and replicate index,
#' unless `overwrite = TRUE`.
#'
#' @param X Numeric data matrix (\eqn{N \times D}).
#' @param M Number of replicates. Default 1, a single density estimate.
#' @param rho Fraction of rows used per replicate. Default `NULL`, all rows.
#' @param subsampleSize Number of rows used per replicate, as an alternative
#'   to `rho`. Default `NULL`.
#' @param replace Whether to draw each replicate's rows with replacement.
#' @param seeds Seeds for the `M` replicates. The default,
#'   `phmSeeds(1000, M)`, matches [fitPHMEnsemble()]'s default `seedBase`.
#' @param numCores Number of replicates fit in parallel (forked processes;
#'   not available on Windows).
#' @param saveDir Folder in which to save each replicate's result for reuse.
#'   `NULL` (default) saves nothing.
#' @param prefix File name prefix for saved replicates. Ignored if `saveDir`
#'   is `NULL`.
#' @param overwrite Whether to refit replicates that are already saved.
#'   Ignored if `saveDir` is `NULL`.
#' @param verbose Whether to print progress messages.
#' @param ... Further arguments to [fitDensityReplicate()], such as `G`,
#'   `modelNames`, `GStep`, `GMax` and `stage2Init`.
#'
#' @examples
#' set.seed(1)
#' dat <- matrix(c(rnorm(200), rnorm(200, 5), rnorm(200, -5)), ncol = 2, byrow = TRUE)
#' params <- globalDensityEstimation(dat, M = 2, rho = 0.8, G = 1:5)
#' length(params)
#'
#' @returns List of lists, one per pooled component, each with the
#' component's `prob`, `mean`, `var` and `class`.
#'
#' @export
globalDensityEstimation <- function(X,
                                    M=1,
                                    rho=NULL,
                                    subsampleSize=NULL,
                                    replace=F,
                                    seeds=phmSeeds(1000, M),
                                    numCores=1,
                                    saveDir=NULL,
                                    prefix="global",
                                    overwrite=F,
                                    verbose=F,
                                    ...) {
    ## Validate inputs
    if (!is.null(rho) && !is.null(subsampleSize)) {
        stop("Cannot specify both rho and subsampleSize parameters")
    }
    if (length(seeds) != M) {
        stop("Length of seeds must be equal to M")
    }

    ## Subsample size, expressed as rho for fitDensityReplicate()
    N <- nrow(X)
    if (!is.null(subsampleSize)) {
        rho <- subsampleSize / N
    } else if (is.null(rho)) {
        rho <- 1
    }
    subsampleSize <- round(rho * N)
    verbosityPrint("Subsample size per replicate:", subsampleSize, vValue=verbose, vLevel=1)

    ## Prepare save directory
    if (!is.null(saveDir) && !dir.exists(saveDir)) dir.create(saveDir, recursive = TRUE)

    ## Fit replicates, in parallel if requested
    apply_func <- lapply
    if (numCores > 1 && M > 1) {
        apply_func <- function(X, FUN) {
            parallelApply(X, FUN, numCores, mc.preschedule=FALSE)
        }
        verbosityPrint("Running in Parallel;", numCores, "cores", vValue=verbose, vLevel=1)
    }

    results <- apply_func(1:M, function(idx) {
        filename <- if (is.null(saveDir)) NA else
            file.path(saveDir, paste(prefix, "ssize", subsampleSize, "Repl", idx, "replicate.rds", sep="_"))

        if (!is.na(filename) && file.exists(filename) && !overwrite) {
            verbosityPrint("Loading Replicate", idx, vValue=verbose, vLevel=2)
            return(readRDS(filename))
        }

        verbosityPrint("Running Replicate", idx, "| Seed:", seeds[idx], vValue=verbose, vLevel=2)
        args <- list(X, idx = idx, seed = seeds[idx], rho = rho, replace = replace, ...)
        if (is.null(args$extendVerbose)) args$extendVerbose <- isTRUE(verbose) | verbose > 0
        if (is.null(args$label)) args$label <- prefix
        res <- do.call(fitDensityReplicate, args)
        if (!is.na(filename)) saveRDS(res, filename)
        res
    })

    failed <- vapply(results, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
    if (any(failed)) {
        stop(sprintf("replicate(s) %s failed", paste(which(failed), collapse = ", ")))
    }

    poolDensityReplicates(results, replicateWeighting = "equal", label = prefix)
}


#' Localized (Ensemble) Density Estimation
#'
#' @description EEP
#'
#' @details FILL ME IN
#' 
#' @inheritParams globalDensityEstimation
#' @param partitionFunc Function taking input matrix X as input and returning a vector of partition labels.
#' @param partitionFuncParams List of parameters to be passed to partitionFunc.
#' @param weightFunc Function that acceps input matrix X and cluster vector as inputs, returns an NxK weight matrix.
#' @param weightFuncParams List of parameters to be passed to weightFunc.
#' @param weightThreshold For weighted GMM, ignore observations with weight below this value.
#' @param preserveClusters Whether to keep the components for a partition together or to treat them separately
#' 
#' @export 
localizedDensityEstimation <- function(X,
                                       partitionFunc=NULL,
                                       partitionFuncParams=list(),
                                       weightFunc=NULL,
                                       weightFuncParams=list(),
                                       weightThreshold=5e-3,
                                       preserveClusters=F,
                                       M=1, 
                                       rho=NULL,
                                       subsampleSize=NULL,
                                       replace=F,
                                       seeds=1:M,
                                       numCores=1,
                                       saveDir=NULL,
                                       prefix="global",
                                       overwrite=F,
                                       verbose=F,
                                       ...) {
    ## Validate inputs
    if (!is.null(rho) && !is.null(subsampleSize)) {
        stop("Cannot specify both rho and subsampleSize parameters")
    }
    if (length(seeds) != M) {
        stop("Length of seeds must be equal to M")
    }
    if (!is.double(weightThreshold)) {
        stop("Must use double for weightThreshold")
    }

    ## Logic for subsample size
    N <- nrow(X)
    subsamp <- F
    if (!is.null(subsampleSize)) {
        verbosityPrint("Subsample Size Provided =", subsampleSize, vValue=verbose, vLevel=1)
        subsamp <- T
    } else if (!is.null(rho)) {
        ## rho provided, use that
        subsampleSize = ceiling(N * rho)
        verbosityPrint("Rho Provided; Subsample Size =", subsampleSize, vValue=verbose, vLevel=1)
        subsamp <- T
    } else {
        subsampleSize <- N
        verbosityPrint("Using Full Sample; N =", subsampleSize, vValue=verbose, vLevel=1)
    }

    ## Prepare save directory
    backup_replicates <- F
    if (!is.null(saveDir)) {
        backup_replicates <- T
        if (!dir.exists(saveDir)) dir.create(saveDir)
        verbosityPrint("Saving to:", saveDir, vValue=verbose, vLevel=1)
    } else {
        saveDir <- ""
    }
    ## If multicore then use mclapply (via parallelApply(), which stops if a worker fails)
    apply_func <- lapply
    if (numCores > 1 && M > 1) {
        apply_func <- function(X, FUN) {
            parallelApply(X, FUN, numCores)
        }
        verbosityPrint("Running in Parallel;", numCores, "cores", vValue=verbose, vLevel=1)
    } else {
        verbosityPrint("Running with Single Process", vValue=verbose, vLevel=1)
    }

    results <- apply_func(1:M, function(idx) {
        seed <- seeds[idx]
        set.seed(seed)

        ## Whether this replicate has already been run; check in saveDir if specified
        filename <- paste(prefix, 
            "ssize", subsampleSize,
            "Repl", idx,
            "f.RDS", sep="_"
        )
        filename <- file.path(saveDir, filename)
        need_to_run <- !backup_replicates || overwrite || !file.exists(filename)

        if (need_to_run) {
            verbosityPrint("Running Replicate", idx, vValue=verbose, vLevel=2)

            ## Subset the observations OR use entire dataset
            if (subsamp) {
                sub_idx <- sample.int(N, size = subsampleSize, replace=replace)
                sub_X <- X[sub_idx, ]
            } else {
                sub_idx <- 1:N
                sub_X <- X
            }

            ## Partition Data + Construct Weight Matrix
            verbosityPrint("Running Clustering", vValue=verbose, vLevel=3)
            partitionFuncParams[["seed"]] <- seed
            partitionFuncParams[["X"]] <- sub_X
            partition <- do.call(partitionFunc, partitionFuncParams)


            verbosityPrint("Constructing Weights", vValue=verbose, vLevel=3)
            weightFuncParams[["X"]] <- sub_X
            weightFuncParams[["partition"]] <- partition
            weights <- do.call(weightFunc, weightFuncParams)


            ## Construct Weighted Density
            label_ids <- sort(unique(partition))
            params <- lapply(label_ids, function(k) {
                ## Get subset corresponding to this cluster
                clust_idx_k <- which(partition == k)
                sub_X_k <- sub_X[clust_idx_k, , drop=F]
                w <- weights[, k, drop=F]
                prob_reweight <- length(clust_idx_k) / nrow(sub_X)

                ## Apply Weight Threshold
                valid <- which(w >= weightThreshold)
                sub_X_k_valid <- sub_X[valid, , drop=F]
                w <- w[valid]

                verbosityPrint("Density Estimation: Cluster", k, vValue=verbose, vLevel=3)
                wmcl <- weightedMclust(sub_X_k_valid, 
                                       init_data=sub_X_k,
                                       weights=w, 
                                       verbose=verbose, ...)
                pars <- constructPmcParamsMclust(wmcl, T)
                pars$class <- paste(idx, k, sep="_")
                if (preserveClusters) {
                    pars$prob <- pars$prob * prob_reweight
                } else {
                    pars <- decomposeParams(list(pars))
                    for (idx in 1:length(pars)) {
                        pars[[idx]]$prob <- pars[[idx]]$prob * prob_reweight
                    }
                }

                pars
            })

            ## Write output
            if (backup_replicates) saveRDS(params, filename)
        } else {
            verbosityPrint("Loading Replicate", idx, vValue=verbose, vLevel=2)
            params <- readRDS(filename)
        }

        if (!preserveClusters) params <- do.call(c, params)

        params
    })
    results <- do.call(c, results)

    ## Normalize the Probabilities
    for (idx in 1:length(results)) {
        results[[idx]]$prob <- results[[idx]]$prob / M
    }

    results
}

###########################################
#### Sample Functions for Localized DE ####
###########################################



#' Naive Cluster Weight Calculation
#'
#' @description Constructs a one-hot weight matrix based on cluster assignments.
#'
#' @details An observation is assigned weight 1 for cluster k if it is a member of cluster k, 0 otherwise.
#' 
#' @param X Data matrix
#' @param partition Partition of observations in X into clusters
#' 
#' @export 
naiveWeightFunc <- function(X, partition) {
  one_hot <- matrix(0, nrow = length(partition), ncol = max(partition))
  one_hot[cbind(seq_along(partition), partition)] <- 1
  
  colnames(one_hot) <- 1:max(partition)
  one_hot
}


#' Distance-based Cluster Weight Calculation
#'
#' @description Constructs a weighting matrix based on the smallest distance to observations in a cluster.
#'
#' @details FILL ME IN
#' 
#' @inheritParams naiveWeightFunc
#' @param removeSelfDistance Whether to count distance to self (i.e., distance to assigned cluster is 0)
#' @param scaling Factor by which to multiply distances. Larger values places more weight on the assigned cluster.
#' @param aggFunc Function to aggregate distances from the observation to the cluster. Default is min (minimum distance used).
#' 
#' @export 
distanceWeightFunc <- function(X, partition, removeSelfDistance=F, scaling=1, aggFunc=min) {
  X <- as.matrix(X)
  N <- nrow(X)
  unique_clusters <- sort(unique(partition))
  K <- length(unique_clusters)
  row_sq <- rowSums(X^2)
  
  ## Get minimum cluster distance
  dist_to_clust <- matrix(0, nrow = N, ncol = K)
  for (k_idx in seq_along(unique_clusters)) {
    k <- unique_clusters[k_idx]
    clust_idx <- which(partition == k)
    
    ## Identity: ||a - b||^2 = ||a||^2 + ||b||^2 - 2<a, b>
    dist_to_k <- outer(row_sq, row_sq[clust_idx], "+")
    dist_to_k <- dist_to_k - 2 * (X %*% t(X[clust_idx, , drop = FALSE]))

    ## Allow Self-Distance?
    if (removeSelfDistance) {
      for (j in seq_along(clust_idx)) {
        dist_to_k[clust_idx[j], j] <- Inf
      }
    }

    ## Negative Values?
    dist_to_k[dist_to_k < 0] <- 0
    
    ## Minimum Distance to observation in cluster
    ## Alternative: matrixStats::rowMins(D2_k))
    dist_to_clust[, k_idx] <- apply(dist_to_k, 1, aggFunc)
  }

  ## Euclidean Distance; above is squared distance
  dist_to_clust <- scaling * sqrt(dist_to_clust)
  
  ## Softmax Distance for Weight
  mindist <- apply(dist_to_clust, 1, min)
  weights <- exp(-(dist_to_clust - mindist))
  weights <- weights / rowSums(weights)
  
  colnames(weights) <- unique_clusters
  
  weights
}


plotDensity2D <- function(paramsList,
                          X=NULL,
                          partition=NULL,
                          colors=RColorBrewer::brewer.pal(12, "Paired"),
                          xlim=NULL,
                          ylim=NULL,
                          gridResolution=200,
                          densityLevels=c(5e-2, 1e-1),
                          densityLevelWidth=0.5,
                          colorDensity=F,
                          textSize=8,
                          legendPosition="none") {

  ## Get component names
  density_classes <- sapply(paramsList, function(x) x$class)

  ## Verify data is 2D / Prepare X
  if (!is.null(X)) {
    if (ncol(X) != 2) stop("X must be 2D")
    X <- data.frame(X)
    colnames(X) <- c("X1", "X2")

    if (!is.null(partition)) {
      X$part <- factor(partition)
    } else {
      X$part <- rep(1, nrow(X))
    }
  }
  
  ## Set limits properly; default is 10% margin
  if (!is.null(X)) {
    if (is.null(xlim)) {
      xlim <- c( min(X[, 1]), max(X[, 1]) )
      xlim[1] <- xlim[1] - abs(xlim[1] * 0.1)
      xlim[2] <- xlim[2] + abs(xlim[2] * 0.1)
    }
    if (is.null(ylim)) {
      ylim <- c( min(X[, 2]), max(X[, 2]) )
      ylim[1] <- ylim[1] - abs(ylim[1] * 0.1)
      ylim[2] <- ylim[2] + abs(ylim[2] * 0.1)
    }    
  } else {
    if (is.null(xlim) || is.null(ylim))
      stop("If X not specified then xlim and ylim must be specified!")
  }

  ## Construct Grid for Evaluation
  mat <- expand.grid(X=seq(min(xlim), max(xlim), length.out=gridResolution),
                     Y=seq(min(ylim), max(ylim), length.out=gridResolution)) %>%
    as.matrix()
  
  ## Density Matrix
  density_mat <- sapply(paramsList, function(x) {
    K <- length(x$prob)
    
    tmp <- sapply(1:K, function(idx) {
      x$prob[idx] * mvtnorm::dmvnorm(mat, x$mean[, idx], x$var[, , idx]) # / sum(x$prob)
    })
    rowSums(tmp)
  })
  colnames(density_mat) <- density_classes

  dens_df <- data.frame(mat, dens=density_mat) %>%
    tidyr::pivot_longer(cols=dplyr::starts_with("dens")) %>%
    dplyr::mutate(name=stringr::str_remove(name, "dens."),
                  name=factor(name, levels=density_classes))

  if (colorDensity) {
    if (length(density_classes) > length(colors)) {
        cat("Suppressing, not enough colors provided.\n")
        dens_df <- dplyr::mutate(dens_df, col="One")
        density_colors <- rep("#000", length(density_classes))
    } else {
        dens_df <- dplyr::mutate(dens_df, col=name)
        density_colors <- colors
    }
  } else {
    dens_df <- dplyr::mutate(dens_df, col="One")
    density_colors <- rep("#000", length(density_classes))
  }

  ## Generate Plot
  plt <- ggplot2::ggplot()
  
  ## Draw observations
  if (!is.null(X)) {
    plt <- plt + ggplot2::geom_point(
      ggplot2::aes(x=X1,
          y=X2,
          color=part),
      alpha=0.5,
      data=X
    ) +
    ggplot2::scale_color_manual(values=colors) + 
    ggnewscale::new_scale_color()
  }
  
  ## Draw Density
  plt <- plt + ggplot2::geom_contour(
    ggplot2::aes(x=X, y=Y, group=name, z=value, color=col), 
    breaks = densityLevels,
    linewidth=densityLevelWidth, 
    data=dens_df) +
    ggplot2::scale_color_manual(values=density_colors) + 
    ggnewscale::new_scale_color()
  
  ## Formatting
  plt <- plt +
    
    ggplot2::scale_x_continuous(limits=xlim) +
    ggplot2::scale_y_continuous(limits=ylim) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position=legendPosition,
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   text=ggplot2::element_text(size=textSize))

  plt
}


##########################################################
### Old Density Estimation Functions; to be Deprecated ###
##########################################################

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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
#' library(mclust)
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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
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
#' dat <- matrix(c(rnorm(200), rnorm(200, 3), rnorm(200, -3)), ncol=2, byrow=TRUE)
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


#' Ensemble density estimation via Mclust
#'
#' @description EEP
#'
#' @details FILL ME IN
#'
#' @param data Numeric matrix for data
#' @param G Number of components for each ensembled model
#' @param replicates Number of replicates over which to ensemble
#' @param subsampSize Number of observations to draw in each replicate. If NULL use entire dataset in each replicate
#' @param saveDir Where to save intermediate density estimates (for rerunning analyses). NULL is not saving
#' @param prefix Filename prefix for saving intermediate density estimates
#' @param verbose Whether to print out logging messages
#' @param numCores Number of cores over which to parallelize
#' @param seeds List of seeds to use across replicates
#' @param ... Parameters to provide to Mclust
#'
#' @examples
#' set.seed(1)
#'
#' @returns List of lists where each sublist contains the
#' proportion, mean, covariance matrix estimates
#' for each component in the ensembled model.
#'
#' @export
constructPmcParamsSubAggMclust <- function(data,
                         replicates,
                         subsampSize=NULL,
                         G=NULL,
                         saveDir=NULL,
                         prefix="subAggMclust_",
                         verbose=F, numCores=1, 
                         seeds=NULL,
                         ...) {

  if (!is.null(seeds)) stopifnot(length(seeds) == replicates)

  ## If we want to save the results, set that up
  if (!is.null(saveDir)) {
    if (!dir.exists(saveDir)) dir.create(saveDir)
  }

  ## Default 10 components
  if (is.null(G)) G <- 1:10

  ## Multicore Processing if specified
  apply_func <- lapply
  if (numCores > 1) apply_func <- function(X, FUN) {
    parallelApply(X, FUN, numCores)
  }

  ## Break down into replicates
  if (is.null(seeds)) {
    seeds <- rep(NULL, replicates)
  }
  if (is.null(subsampSize)) {
    subsamp_dat <- lapply((1:replicates), function(idx) {
      list(
        idx=idx,
        dat=data,
        seed=seeds[idx]
      )
    })
  } else {
    subsamp_idx <- lapply(1:replicates, function(idx) {
      set.seed(seeds[idx])
      sample.int(nrow(data), subsampSize, replace=F)
    })
    subsamp_dat <- lapply(1:replicates, function(idx) {
      list(
        idx=idx,
        dat=data[subsamp_idx[[idx]], ],
        seed=seeds[idx]
      )
    })
  }


  ## Run the models of interest
  gmm_res <- apply_func(subsamp_dat, function(sub_obj) {
    idx <- sub_obj$idx
    subsamp_dat <- sub_obj$dat
    seed <- sub_obj$seed
    if (!is.null(seed)) set.seed(seed)

    ## Generate filename
    filename <- paste(prefix, "N", subsampSize, "G", max(G),
                      "repl", idx, "f.RData", sep="_")
    filename <- file.path(saveDir, filename)

    ## Either this file exists or filename is empty because saveDir is NULL
    need_to_run <- !file.exists(filename) || is.null(saveDir)

    if (need_to_run) {
      mcl <- mclust::Mclust(subsamp_dat, G=G,
                            verbose=verbose, ...)
      params <- constructPmcParamsMclust(mcl)
      for (g in 1:mcl$G) {
        params[[g]]$class <- paste(params[[g]]$class, idx, sep="_")
      }

      if (!is.null(saveDir)) save(params, subsamp_dat, file=filename)
    } else {
      load(filename)
    }

    params
  })

  ## Average over all estimated densities
  gmm_res <- do.call(c, gmm_res)
  for (idx in 1:length(gmm_res)) {
    gmm_res[[idx]]$prob <- gmm_res[[idx]]$prob / replicates
  }

  gmm_res
}

## Older localized ensemble; superseded by localizedDensityEstimation().
#' @noRd
constructPmcParamsLocalizedEnsemble <- function(data,
                                                clustFunc,
                                                replicates,
                                                subsampSize=NULL,
                                                G=NULL,
                                                saveDir=NULL,
                                                prefix="localizedEnsemble_",
                                                weighted=T,
                                                verbose=F, numCores=1, 
                                                seeds=NULL,
                                                ...) {
  if (is.null(seeds)) {

    seeds = 1:replicates
  } else {
    stopifnot(length(seeds) == replicates)
  }

  ## If we want to save the results, set that up
  if (!is.null(saveDir)) {
    if (!dir.exists(saveDir)) dir.create(saveDir)
  }

  ## Default 10 components
  if (is.null(G)) G <- 1:10

  ## Multicore Processing if specified
  apply_func <- lapply
  if (numCores > 1) apply_func <- function(X, FUN) {
    parallelApply(X, FUN, numCores)
  }

  ## Whether to do weighted or naive density
  density_func <- if (weighted) {
    constructPmcParamsWeightedPartition
  } else {
    constructPmcParamsPartition
  }

  ## Break down into replicates
  if (is.null(seeds)) {
    seeds <- rep(NULL, replicates)
  }
  if (is.null(subsampSize)) {
    subsamp_dat <- lapply((1:replicates), function(idx) {
      list(
        idx=idx,
        dat=data,
        seed=seeds[idx]
      )
    })
  } else {
    subsamp_idx <- lapply(1:replicates, function(idx) {
      set.seed(seeds[idx])
      sample.int(nrow(data), subsampSize, replace=F)
    })
    subsamp_dat <- lapply(1:replicates, function(idx) {
      list(
        idx=idx,
        dat=data[subsamp_idx[[idx]], ],
        seed=seeds[idx]
      )
    })
  }

  ## From here do the thing
  ensemble <- apply_func(subsamp_dat, function(sub_obj) {
    idx <- sub_obj$idx
    subsamp_dat <- sub_obj$dat
    seed <- sub_obj$seed
    if (!is.null(seed)) set.seed(seed)

    ## Generate filename
    filename <- paste(prefix, ifelse(weighted, "_Weighted_", "_Naive_"),
                      "N", subsampSize, "G", max(G),
                      "repl", idx, "f.RData", sep="_")
    filename <- file.path(saveDir, filename)

    ## Either this file exists or filename is empty because saveDir is NULL
    need_to_run <- !file.exists(filename) || is.null(saveDir)

    if (need_to_run) {
      partition = clustFunc(subsamp_dat, seed=seed)

      params <- density_func(partition, subsamp_dat, G=G, ...)
      for (g in seq_along(params)) {
        params[[g]]$class <- paste(idx, params[[g]]$class, sep="_")
      }

      if (!is.null(saveDir)) save(params, subsamp_dat, file=filename)
    } else {
      load(filename)
    }

    params
  })

  ## Average over all estimated densities
  ensemble <- do.call(c, ensemble)
  for (idx in 1:length(ensemble)) {
    ensemble[[idx]]$prob <- ensemble[[idx]]$prob / replicates
  }

  decomposeParams(ensemble)
}



#' Ensemble weighted density estimation for partitions
#'
#' @description EEP
#'
#' @details PEW
#'
#' @inheritParams constructPmcParamsWeightedPartition
#' @param replicates Number of replicates over which to be ensembled
#' @param subsampSize Maximum number of observations used per cluster in each replicate
#' @param saveDir Where to save intermediate density estimates (for rerunning analyses). NULL is not saving
#' @param prefix Filename prefix for saving intermediate density estimates
#' @param verbose Whether to print out logging messages
#' @param numCores Number of cores over which to parallelize
#'
#' @examples
#' set.seed(1)
#'
#' @returns List of lists where each sublist contains the
#' proportion, mean, covariance matrix estimates
#' for each cluster.
#'
#' @export
constructPmcParamsSubAggPartition <- function(data,
                            partition,
                            replicates,
                            subsampSize,
                            linkFunc=min,
                            threshold=1e-4,
                            saveDir=NULL,
                            prefix="subAggPartition_",
                            numCores=1, 
                            verbose=F,
                            ...) {
  stop("This method does not work; use constructPmcParamsLocalizedEnsemble")
  ## Prepare for use across all cores
  label_ids <- sort(unique(partition))
  K <- length(label_ids)
  clust_idx <- lapply(label_ids, function(x) {which(partition == x)})
  clust_size <- sapply(clust_idx, length)
  clust_prob <- clust_size / sum(clust_size)

 apply_func <- if(numCores == 1) lapply else function(x, f) parallelApply(x, f, numCores)

  params_raw <- lapply(1:replicates, function(repl) {
    if (verbose) cat("Replicate", repl, "\n")
    ## Prepare subsamples from each cluster
    subsamp_idx <- lapply(clust_idx, function(v) {
      if (length(v) <= subsampSize) return(v)
      sample(v, subsampSize)
    })

    subsamp_size <- sapply(subsamp_idx, length)

    subsamp_reweights_K <- (clust_size) / (subsamp_size)
    subsamp_reweights <- lapply(1:K, function(idx) rep(subsamp_reweights_K[idx], subsamp_size[idx]))
    subsamp_reweights <- do.call(c, subsamp_reweights)
    subsamp_labels <- lapply(1:K, function(idx) rep(label_ids[idx], subsamp_size[idx]))
    subsamp_labels <- do.call(c, subsamp_labels)
    subsamp_mat <- lapply(subsamp_idx, function(v) {
      data[v, ]
    })
    subsamp_mat <- do.call(rbind, subsamp_mat)
    if (verbose) cat("\tSubsampling Complete", repl, "\n")


    ## Generate filename
    filename_w <- paste(prefix, "N", subsampSize, "repl", repl, "weights.RData", sep="_")
    filename_w <- file.path(saveDir, filename_w)

    ## Either this file exists or filename is empty because saveDir is NULL
    need_to_run <- !file.exists(filename_w) || is.null(saveDir)
    if (need_to_run) {
      ## Compute Weights for subsampled matrix
      dist_to_clust <- apply_func(label_ids, function(k) {
        clust_idx <- which(subsamp_labels == k)
        clust_mat <- subsamp_mat[clust_idx, ]

        sapply(1:nrow(subsamp_mat), function(idx) {
          clust_dist <- colSums((t(clust_mat) - subsamp_mat[idx, ])^2)
          linkFunc(clust_dist[which(clust_dist > 0)])
        })
      })
      dist_to_clust <- do.call(cbind, dist_to_clust)
      mindist <- apply(dist_to_clust, 1, min)
      weights <- exp(-1 * (dist_to_clust - mindist))
      weights <- weights/rowSums(weights)
      reweights <- apply(weights, 2, function(v) v * subsamp_reweights)

      if (!is.null(saveDir)) save(weights, reweights, file=filename_w)
    } else {
      load(filename_w)
    }
    if (verbose) cat("\tWeighting Complete\n")

    ## Now we estimate the cluster-specific densities
    if (verbose) cat("\tParam Estimation Beginning\n")
    pars <- apply_func(label_ids, function(k) {
      filename_k <- paste(prefix, "N", subsampSize, "repl", repl, "cluster", k, "density.RData", sep="_")
      filename_k <- file.path(saveDir, filename_k)
      need_to_run <- !file.exists(filename_k) || is.null(saveDir)

      idx <- which(label_ids == k)
      if (need_to_run) {
        ## Cluster-specific weights and observations
        dat <- subsamp_mat[which(subsamp_labels == k), , drop=F]
        w <- weights[, idx]
        rw <- reweights[, idx]

        ## Threshold observations with too small weights
        valid <- which(w >= threshold)
        data_valid <- subsamp_mat[valid, , drop=F]
        w <- w[valid]
        rw <- rw[valid]

        wmcl <- reweightedMclust(data_valid, weights = w, reweights=rw,
                                init_data = dat, ...)
        if (!is.null(saveDir)) save(wmcl, file=filename_k)
      } else {
        load(filename_k)
      }

      if (verbose) cat("\t\tParam Estimation Cluster", k, "Complete\n")
      ## Combine the parameters
      pars_ <- constructPmcParamsMclust(wmcl, T)
      pars_$class <- k
      pars_$prob <- pars_$prob * clust_prob[idx] / replicates
      pars_
    })

    pars
  })
  if (verbose) cat("\tParam Estimation Complete\n")

  lapply(1:K, function(k) {
    pars <- lapply(params_raw, function(x) x[[k]])

    res <- Reduce(mergeParams, pars)
    res$class <- k

    res
  })
}
