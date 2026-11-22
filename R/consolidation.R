## Files for density component consolidation

#' Pairwise \eqn{P_{\rm{mc}}} Matrix computation
#'
#' @description TODO: FILL ME IN
#'
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param mcSamples Monte Carlo Samples for Pmc calculation
#' @param numCores Number of cores for computations
#' @param equalProbs Whether to ignore prior weights for each component (i.e., Pmc based solely on mean/variance)
#' @param threshold Pmc threshold for Pairwise Pmc, all values below threshold are treated as 0
#' @param verbose Whether to print progress messages
#' @param thresholdTol Tolerance for the search over pairs: the search for the pair whose \eqn{P_{\rm{mc}}} is closest to `threshold` stops once it is within this distance
#'
#' @return \eqn{K \times K} matrix with Pairwise \eqn{P_{\rm{mc}}} values for each pair of clusters
#'
#' @export
computePairwisePmcMatrix <- function(paramsList, mcSamples, numCores=1, threshold=0, thresholdTol=1e-2, equalProbs=F, verbose=T) {
  K <- length(paramsList)
  mat <- expand.grid(i=1:K, j=1:K)
  mat <- mat[which(mat[, 1] > mat[, 2]), , drop=F]
  display_digits <- 1+floor(-log10(thresholdTol))

  ## Multicore functionality
  apply_func <- lapply
  if (numCores > 1) {
    apply_func <- function(x, func) parallelApply(x, func, numCores)
  }

  ## Equalize probabilities across all components
  if (equalProbs) {
    tgt_probs <- 1/K
    for (idx in 1:length(paramsList)) {
      paramsList[[idx]]$prob <- paramsList[[idx]]$prob * tgt_probs / sum(paramsList[[idx]]$prob)
    }
  }

  if (verbose) cat("Mahalanobis Beginning\r")
  elapsed_dist <- system.time({
    dist_vals <- apply_func(1:nrow(mat), function(idx) {
      v <- mat[idx, ]
      i <- as.numeric(v[1])
      j <- as.numeric(v[2])

      par_i <- paramsList[[i]]
      par_j <- paramsList[[j]]

      if (equalProbs) {
        cov_pool <- par_i$var + par_j$var
        cov_pool <- cov_pool[, , 1] / 2
      } else {
        total_prob <- par_i$prob + par_j$prob

        cov_pool <- par_i$prob * par_i$var + par_j$prob * par_j$var
        cov_pool <- cov_pool[, , 1] / total_prob
      }

      mah_dist <- mahalanobis(as.vector(par_i$mean), 
                              as.vector(par_j$mean), 
                              cov_pool)
      mah_dist
    })
    dist_vals <- do.call(c, dist_vals)
    dist_order <- order(dist_vals)
    mat$dist <- dist_vals
    mat <- mat[dist_order, ]
  })
  if (verbose) cat("Mahalanobis Complete:", elapsed_dist[3], "secs\n")


  if (verbose) cat("Search Beginning\n")
  elapsed_search <- system.time({
    pmc <- 0
    lbound <- nrow(mat)
    ubound <- 1
    index <- lbound
    delta <- threshold

    while(abs(delta) > thresholdTol && threshold != 0) {
      prev_index <- index
      index <- ceiling((lbound + ubound) / 2)
      vals <- mat[index, ]

      par_i <- paramsList[[as.numeric(vals[1])]]
      par_j <- paramsList[[as.numeric(vals[2])]]

      total_prob <- par_i$prob + par_j$prob
      par_i$prob <- par_i$prob / total_prob
      par_j$prob <- par_j$prob / total_prob

      pmc <- computeMonteCarloPmc(list(par_i, par_j), mcSamples)
      delta <- threshold - pmc
      if (verbose) cat("\tIndex:", 
                        stringr::str_pad(index, ceiling(log10(nrow(mat))), "left"),
                        "| Pmc (delta):", round(pmc, display_digits),
                        "(",round(delta, display_digits), ")", "\n")

      if (index == 1) break

      if (delta > 0) {
        lbound <- index
      } else (
        ubound <- index
      )
    }
  })
  prev_index <- index

  if (verbose) cat("Search complete:", elapsed_search[3], "secs\n")
  if (verbose) cat("\tSearch Index:", prev_index, "/", nrow(mat), "\n")

  elapsed_sample <- system.time({
    mc_draws <- apply_func(paramsList, function(pars) {
      mvtnorm::rmvnorm(mcSamples,
              mean=pars$mean,
              sigma=pars$var[, , 1])
    })
  })
  if (verbose) cat("Sampling Complete:", elapsed_sample[3], "secs\n")

  if (verbose) cat("Evaluation Beginning\r")
  pmcMat <- matrix(0, nrow=K, ncol=K)
  elapsed_eval <- system.time({
    tmp <- apply_func(1:prev_index, function(index) {
      vals <- mat[index, ]
      i <- as.numeric(vals[1])
      j <- as.numeric(vals[2])
      par_i <- paramsList[[i]]
      par_j <- paramsList[[j]]

      group_split <- rmultinom(1, mcSamples, 
                               c(par_i$prob, par_j$prob))

      samp_i <- mc_draws[[i]]
      samp_i <- samp_i[sample(1:mcSamples, group_split[1], T), ]
      samp_j <- mc_draws[[j]]
      samp_j <- samp_j[sample(1:mcSamples, group_split[2], T), ]
      samp <- rbind(samp_i, samp_j)
      post_mat <- computePosteriorProbMatrix(list(par_i, par_j), samp)

      pmc <- sum(post_mat * (1 - post_mat)) / nrow(samp)

      list(i=i, j=j, pmc=pmc)
    })

    for (idx in 1:length(tmp)) {
      pmcMat[tmp[[idx]]$i, tmp[[idx]]$j] <- tmp[[idx]]$pmc
      pmcMat[tmp[[idx]]$j, tmp[[idx]]$i] <- tmp[[idx]]$pmc
    }
  })
  if (verbose) cat("Evaluation Complete:", elapsed_eval[3], "secs\n")
  elapsed_total <- elapsed_sample[3] + elapsed_eval[3] + elapsed_dist[3] + elapsed_search[3]
  if (verbose) cat("Mission Complete:", elapsed_total, "secs\n")

  list(
    mat=pmcMat,
    elapsed_sample=elapsed_sample[3],
    elapsed_eval=elapsed_eval[3],
    elapsed_search=elapsed_search[3],
    elapsed_mahal=elapsed_dist[3],
    elapsed=elapsed_total
  )
}

#' Consolidate redundant components
#'
#' @description TODO: FILL ME IN
#'
#' @param paramsList List containing lists with each component GMM parameters. See `generateDistbnFunc` for format of components.
#' @param pmcMatrix Pairwise Pmc matrix
#' @param threshold Pmc threshold for component consolidation. Any pairs of components with Pmc > threshold will be combined.
#'
#' @return List of parameters with redundant components consolidated inso a single component.
#'
#' @export
consolidateParams <- function(paramsList, pmcMatrix, threshold) {
  merge_components <- list()
  for (i in 1:nrow(pmcMatrix)) {
    for (j in i:nrow(pmcMatrix)) {
      if (i == j) next
      if (pmcMatrix[i, j] > threshold) {
        merge_components[[length(merge_components) + 1]] <- c(i, j)
      }
    }
  }

  merge_clusters <- list()
  for (posn in merge_components) {
    i <- posn[1]
    j <- posn[2]

    if (length(merge_clusters) > 0) {
      valid_clust <- sapply(merge_clusters, function(mc) {
        i %in% mc || j %in% mc
      })
      valid_idx <- which(valid_clust)
    } else {
      valid_idx <- numeric()
    }


    if (length(valid_idx) == 0) {
      merge_clusters[[length(merge_clusters) + 1]] <- posn
    } else if (length(valid_idx) == 1) {
      merge_clusters[[valid_idx]] <- c(merge_clusters[[which(valid_clust)]], posn)
      merge_clusters[[valid_idx]] <- unique(merge_clusters[[valid_idx]])
    } else {
      base_idx <- valid_idx[1]
      merge_clusters[[base_idx]] <- do.call(c, merge_clusters[valid_idx])
      for (idx in length(valid_idx):2) merge_clusters[[valid_idx[idx]]] <- NULL

      merge_clusters[[base_idx]] <- c(merge_clusters[[base_idx]], posn)
      merge_clusters[[base_idx]] <- unique(merge_clusters[[base_idx]])
    }
  }

  ## Grab parameters to be consolidated
  group_params <- lapply(merge_clusters, function(v) {
    paramsList[v]
  })

  for (idx in 1:length(paramsList)) {
    to_be_merged <- any(sapply(merge_clusters, function(x) idx %in% x))
    if (!to_be_merged) {
      group_params[[length(group_params) + 1]] <- list(paramsList[[idx]])
    }
  }

  ## Consolidate parameters
  consolidated_params <- lapply(group_params, function(v) {
    if (length(v) == 1) return(v[[1]])

    probs <- sapply(v, function(vv) vv$prob)
    means <- lapply(v, function(vv) vv$mean)
    vars <- lapply(v, function(vv) vv$var)

    ## Combine means
    mu <- lapply(1:length(probs), function(idx) {
      probs[idx] * means[[idx]] / sum(probs)
    })
    mu <- Reduce(`+`, mu)
    mu <- rowSums(mu)
    mu <- matrix(mu, ncol=1)

    ## Combine variances
    sigma <- lapply(1:length(probs), function(idx) {
      probs[idx] * vars[[idx]] / sum(probs)
    })
    sigma <- Reduce('+', sigma)

    list(
      prob=sum(probs),
      mean=mu,
      var=sigma,
      class=paste(sapply(v, function(vv) vv$class), collapse="|")
    )
  })
}