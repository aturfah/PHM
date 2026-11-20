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
    if (verbose) cat("Subsample size per replicate:", subsampleSize, "\n")

    ## Prepare save directory
    if (!is.null(saveDir) && !dir.exists(saveDir)) dir.create(saveDir, recursive = TRUE)

    ## Fit replicates, in parallel if requested
    apply_func <- lapply
    if (numCores > 1 && M > 1) {
        apply_func <- function(X, FUN) {
            parallel::mclapply(X, FUN, mc.cores=numCores, mc.preschedule=FALSE)
        }
        if (verbose) cat("Running in Parallel;", numCores, "cores\n")
    }

    results <- apply_func(1:M, function(idx) {
        filename <- if (is.null(saveDir)) NA else
            file.path(saveDir, paste(prefix, "ssize", subsampleSize, "Repl", idx, "replicate.rds", sep="_"))

        if (!is.na(filename) && file.exists(filename) && !overwrite) {
            if (verbose) cat("..Loading Replicate", idx, "\n")
            return(readRDS(filename))
        }

        if (verbose) cat("..Running Replicate", idx, "| Seed:", seeds[idx], "\n")
        args <- list(X, idx = idx, seed = seeds[idx], rho = rho, replace = replace, ...)
        if (is.null(args$extendVerbose)) args$extendVerbose <- verbose
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
        if (verbose) cat("Subsample Size Provided =", subsampleSize, "\n")
        subsamp <- T
    } else if (!is.null(rho)) {
        ## rho provided, use that
        subsampleSize = ceiling(N * rho)
        if (verbose) cat("Rho Provided; Subsample Size =", subsampleSize, "\n")
        subsamp <- T
    } else {
        subsampleSize <- N
        if (verbose) cat("Using Full Sample; N =", subsampleSize, "\n")
    }

    ## Prepare save directory
    backup_replicates <- F
    if (!is.null(saveDir)) {
        backup_replicates <- T
        if (!dir.exists(saveDir)) dir.create(saveDir)
    } else {
        saveDir <- ""
    }

    ## If multicore then use mclapply
    apply_func <- lapply
    if (numCores > 1 && M > 1) {
        apply_func <- function(X, FUN) {
            parallel::mclapply(X, FUN, mc.cores=numCores)
        }
        if (verbose) cat("Running in Parallel;", numCores, "cores\n")
    } else {
        if (verbose) cat("Running with Single Process\n")
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
            if (verbose) cat("..Running Replicate", idx, "\n")

            ## Subset the observations OR use entire dataset
            if (subsamp) {
                sub_idx <- sample.int(N, size = subsampleSize, replace=replace)
                sub_X <- X[sub_idx, ]
            } else {
                sub_idx <- 1:N
                sub_X <- X
            }

            ## Partition Data + Construct Weight Matrix
            if (verbose) cat("....Running Clustering\n")
            partitionFuncParams[["seed"]] <- seed
            partitionFuncParams[["X"]] <- sub_X
            partition <- do.call(partitionFunc, partitionFuncParams)


            if (verbose) cat("....Constructing Weights\n")
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

                if (verbose) {
                    cat("....Density Estimation: Cluster", k, "\n")
                }

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
            if (verbose) cat("..Loading Replicate", idx, "\n")
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