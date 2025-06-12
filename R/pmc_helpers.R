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
  t(apply(distbn.mat, 1, function(x) if (sum(x) != 0) x / sum(x) else x + 1 / length(x)))
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


weightedMclust <- function(data, weights, 
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
  
  hc_init <- hc(data = init_data, 
                modelName = mclust::mclust.options("hcModelName"), 
                use = mclust::mclust.options("hcUse"))
  
  ## Prepare the for loop
  if (is.null(G)) G <- 1:10
  if (is.null(modelNames)) modelNames <- mclust::mclust.options("emModelNames")
  
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

      do.call("me.weighted", c(list(weights=weights), mcl))
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
