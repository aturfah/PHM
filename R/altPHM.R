#' PHM Algorithm (supporting alternate merging strategies)
#' 
#' @description TODO: Fill me in 
#' 
#' @param par1
#' 
#' @examples 
#' set.seed(1)
#' 
#' @return
#' What is the return object
#' 
#' @export
PHMv2 <- function(paramsList=NULL,
                  deltaPmc=NULL,
                  scaling=c("unscaled", "average", "alpha"),
                  monteCarlo=T,
                  numCores=1,
                  verbose=F,
                  ...) {

  scaling <- match.arg(scaling)

  if (is.null(paramsList) && is.null(deltaPmc)) stop("One of paramsList or deltaPmc must be specified")
  
  ## Verify paramsList has prob, mean, and var
  if (!is.null(paramsList)) {
    tmp <- lapply(paramsList, function(x) x$prob)
    if (any(sapply(tmp, is.null))) stop("Missing $prob in at least one element of paramsList")
    tmp <- lapply(paramsList, function(x) x$mean)
    if (any(sapply(tmp, is.null))) stop("Missing $mean in at least one element of paramsList")
    tmp <- lapply(paramsList, function(x) x$var)
    if (any(sapply(tmp, is.null))) stop("Missing $var in at least one element of paramsList")    
  }
  
  ## Compute deltaPmc if not provided
  if (is.null(deltaPmc)) {
    if (verbose) cat("Computing deltaPmc matrix...\n")
    if (monteCarlo) {
      deltaPmc <- computeMonteCarloDeltaPmcMatrix(
        paramsList,
        numCores=numCores,
        ...,
        verbose=verbose
      )
    } else {
      deltaPmc <- computeDeltaPmcMatrix(paramsList, ...)
    }
  } else {
    ## If provided do some validation
    if (ncol(deltaPmc) != nrow(deltaPmc)) stop("deltaPmc must be a square matrix")
    
    if (!is.null(paramsList)) {
      if (length(paramsList) != ncol(deltaPmc)) stop("paramsList mismatch with deltaPmc")
    }
  }
  
  ## Set up merging procedure
  K <- ncol(deltaPmc)
  component_tracker <- lapply(1:K, function(x) x)
  merge_components <- list(NA)
  merge_deltaPmc <- numeric(K)
  merge_value <- numeric(K)
  alpha_vec <- sapply(params, function(x) sum(x$prob))
  
  orig_deltaPmc <- deltaPmc

  ## Results
  if (verbose) cat("Commencing merging procedure\n")
  while (K > 1) {
    ## Matrix for merging rule
    value_matrix <- 2 * deltaPmc
    if (scaling != "unscaled") {
        if (scaling == "average") {
            ## Get the Average Component Pmc
            comp_size <- sapply(component_tracker, length)
            sz <- outer(comp_size, comp_size)
            value_matrix <- value_matrix / sz
        } else if (scaling == "alpha") {
            ## Scale by alpha_i * alpha_j
            alph <- outer(alpha_vec, alpha_vec)
            value_matrix <- value_matrix / alph
        } else {
            stop("Unsopported scaling provided")
        }
    }
    value_matrix[lower.tri(value_matrix)] <- 0

    ## Get index of maximal element
    max_val <- max(value_matrix)
    max_ind <- which(value_matrix == max_val, arr.ind=T)
    i <- max_ind[1, 1]
    j <- max_ind[1, 2]
    max_pmc <- 2 * deltaPmc[i, j]
    
    ## Update Pmc Matrix
    row_i <- deltaPmc[i, ]
    row_j <- deltaPmc[j, ]

    new_row <- row_i + row_j
    new_row[i] <- 0
    new_row <- new_row[-j]
    
    deltaPmc <- deltaPmc[, -j, drop=F]
    deltaPmc <- deltaPmc[-j, , drop=F]
    deltaPmc[i, ] <- new_row
    deltaPmc[, i] <- new_row
    
    ## Track which components are part of which component
    component_tracker[[i]] <- c(
      component_tracker[[i]],
      component_tracker[[j]]
    )
    component_tracker[[j]] <- NULL  

    ## Update component probability matrix
    alpha_vec[i] <- alpha_vec[i] + alpha_vec[j]
    alpha_vec <- alpha_vec[-j]
    
    ## Construct output
    merge_components[[K]] <- c(i=i, j=j)
    merge_deltaPmc[K] <- max_pmc
    merge_value[K] <- max_val

    K <- K-1
  }
  
  if (verbose) cat("Complete!\n")
  list(
    deltaPmc=orig_deltaPmc,
    paramsList=paramsList,
    mergeComps=merge_components,
    mergeValues=merge_value,
    mergeDeltaPmc=merge_deltaPmc,
    mergeCriterion=scaling
  )
}


#' Convert an original PHM object to the new PHM format
#' 
#' @param phm Output from [PHM()] function
#' 
#' @export 
convertToPHMv2 <- function(phm) {
  ## Convert an original PHM to the new format
  res <- list()
  
  res$deltaPmc <- phm[[length(phm)]]$pmc_matrix
  res$paramsList <- phm[[length(phm)]]$params
  res$mergeComps <- lapply(phm, function(x) x$merge_components)
  res$mergeComps[[1]] <- NA
  res$mergeValues <- sapply(phm, function(x) x$pmc_change)[1:(length(phm)-1)]
  res$mergeValues <- c(NA, res$mergeVals)
  res$mergeCriterion <- "unscaled"
  res$mergeDeltaPmc <- res$mergeVals

  res
}

#############################
#### Recover X Functions ####
#############################

#' Recover merged groups at given value of \eqn{K}
#' z
#' @export 
recoverGroupsv2 <- function(phm, k) {
  K <- ncol(phm$deltaPmc)
  ct <- lapply(1:K, function(x) x)

  for (idx in K:(k+1)) {
    res <- phm$mergeComp[[idx]]
    i <- res[1]
    j <- res[2]
    ct[[i]] <- c(
      ct[[i]],
      ct[[j]]
    )
    ct[[j]] <- NULL
  }
  
  ct
}

#' Recover posterior at given value of \eqn{K}
#' 
#' @export 
recoverPosteriorv2 <- function (phm, k, posterior) {
  grps <- recoverGroups(phm, k)

  post <- lapply(grps, function (idx) {
    rowSums(posterior[, idx, drop=F])
  })
  
  do.call(cbind, post)
}

#' Recover deltaPmc matrix at given value of \eqn{K}
#' 
#' @export 
recoverDeltaPmcv2 <- function(phm, k) {
  grps <- recoverGroups(phm, k)
  deltaPmc <- phm$deltaPmc

  tmp_mat <- lapply(grps, function(idx) {
    tmp <- rowSums(deltaPmc[, idx, drop=F])
    tmp[idx] <- 0
  })
  tmp_mat <- do.call(cbind, tmp_mat)

  output <- lapply(grps, function(idx) {
    colSums(tmp_mat[idx, , drop=F])
  })

  do.call(rbind, output)
}

#####################################
#### Visualize New PHM Functions ####
#####################################

constructVisData <- function(phmObj,
                             scaleHeights="unscaled",
                             groupProbs=NULL) {
  K <- ncol(phmObj$deltaPmc)
  alpha_vec <- sapply(phmObj$paramsList, function(x) sum(x$prob))

  ## Whether or not we track groupProbs
  gprobs_unspec <- FALSE
  if (is.null(groupProbs)) {
    gprobs_unspec <- TRUE
    groupProbs <- rep(1, K)
  }

  ## Which components are merged
  merge_components <- phmObj$mergeComps[K:2]


  ## For the original PHM scale the heights based 
  ## Alternate scalings shouldn't need this to occur
  height <- phmObj$mergeValues[K:2]
  if (phmObj$mergeCriterion == "unscaled") {
    clust_sizes <- rep(list(1), K)
    sizes <- numeric(K-1)

    for (idx in 1:(K-1)) {
      mc <- merge_components[[idx]]
      mc1 <- mc[1]
      mc2 <- mc[2]

      clust_sizes[[mc1]] <- clust_sizes[[mc1]] + clust_sizes[[mc2]]
      clust_sizes[[mc2]] <- NULL
      sizes[idx] <- clust_sizes[[mc1]]
    }
    print(height)
    print(choose(sizes, 2))
    height <- height / choose(sizes, 2)
    print(height)
  }

  if (scaleHeights == "unscaled") {
    height <- 1 / height
  } else if (scaleHeights == "log10") {
    height <- -log10(height)
  } else if (scaleHeights %in% c("pmcdist")) {
    load(system.file("extdata", "pmc_scale_function.RData", 
                     package = "PHM"))
    height <- inv_log10(log10(height))
    height <- height^2
    height <- height + (1:length(height)) * 1e-6 ## Slight height offset
  }
  print(height)

  ## Track components merging
  output <- data.frame()
  groupProbs_new <- groupProbs
  merge_tree <- as.list(1:K)
  height_tracker <- rep(list(list(base=0, height=0)), K)
  component_id_map <- 1:K
  base_height <- 0
  for (idx in 1:(K-1)) {
    mcs <- unname(merge_components[[idx]])
    hgt <- height[idx]
    ## Verify that height doesn't decrease
    for (posn in mcs) {
      hgt <- max(hgt, height_tracker[[posn]]$base + 1e-6)
    }

    ## Set new minimum height
    for (posn in 1:length(height_tracker)) height_tracker[[posn]]$height <- hgt
    
    ## Get the new rows
    new_rows <- rbind(
      c(
        ID=component_id_map[mcs[1]],
        y=height_tracker[[mcs[1]]]$base,
        yend=height_tracker[[mcs[1]]]$height,
        gprob=groupProbs_new[mcs[[1]]]
      ),
      c(
        ID=component_id_map[mcs[2]],
        y=height_tracker[[mcs[2]]]$base,
        yend=height_tracker[[mcs[2]]]$height,
        gprob=groupProbs_new[mcs[[2]]]
      )
    )
    
    if (nrow(output) == 0) {
      output <- data.frame(new_rows)
    } else {
      output <- rbind(output, new_rows)
    }
    
    ## Update groupProbs post-merge
    if (gprobs_unspec) {
      groupProbs_new <- groupProbs_new[-1] ## These are all 1 so just remove first element
    } else {
      probs1 <- alpha_vec[mcs[[1]]]
      probs2 <- alpha_vec[mcs[[2]]]
      
      ## Update alpha and groupProbs
      groupProbs_new[mcs[1]] <- (
        groupProbs_new[mcs[1]] * probs1 + groupProbs_new[mcs[2]] * probs2
      ) / (probs1 + probs2)
      groupProbs_new <- groupProbs_new[-mcs[2]]
      
      alpha_vec[mcs[[1]]] <- alpha_vec[mcs[[1]]] + alpha_vec[mcs[[2]]]
      alpha_vec <- alpha_vec[-mcs[[2]]]
    }
    
    ## Make note of the combined merge nodes so we have the merges stored somewhere
    merge_tree[[mcs[1]]] <- list(
      left=merge_tree[[mcs[1]]],
      right=merge_tree[[mcs[2]]],
      order=idx
    )
    merge_tree[[mcs[2]]] <- NULL
    
    ## Remove the component and re-map
    base_height <-  hgt
    height_tracker[[mcs[2]]] <- NULL
    height_tracker[[mcs[1]]]$base <- base_height
    
    
    component_id_map <- component_id_map[-mcs[2]]
  }
  
  order_x <- function(merge_res) {
    if (typeof(merge_res) == "list") {
      if (length(merge_res) == 1) {
        return(order_x(merge_res[[1]]))
      } else {
        return(c(order_x(merge_res[[1]]), order_x(merge_res[[2]])))
      }
    } else {
      return(merge_res)
    }
  }
  x_posns <- order_x(merge_tree)
  
  map_xposns <- function(vec) {
    sapply(vec, function(x) which(x_posns == x))
  }
  output <- dplyr::mutate(output, 
                          x=ifelse(y==0, map_xposns(ID), NA))
  
  while(any(is.na(output))) {
    output <- output %>%
      dplyr::left_join(output, by=c("y"="yend")) %>%
      dplyr::rename(x=x.x,
                    ID=ID.x,
                    # pmc=pmc.x,
                    # pmc_pct=pmc_pct.x,
                    # pmc_change=pmc_change.x,
                    gprob=gprob.x) %>%
      dplyr::group_by(ID, y, yend, x, gprob, 
                      #pmc, pmc_change, pmc_pct
                      ) %>%
      dplyr::summarize(x=ifelse(all(is.na(x)), mean(x.y), mean(x)),
                       xend=x, .groups="keep") %>%
      dplyr::ungroup() %>%
      dplyr::arrange(-yend) %>%
      dplyr::select(-dplyr::ends_with(".y"))
  }
  
  ## Add horizontal components
  horiz_comps <- output %>%
    dplyr::group_by(yend, #pmc, pmc_change, pmc_pct
                    ) %>%
    dplyr::summarize(
      xend=max(x),
      x=min(x),
      gprob=mean(gprob),
      .groups="keep"
    ) %>%
    dplyr::mutate(y=yend)
  output <- dplyr::bind_rows(
    output,
    horiz_comps
  )

  list(
    df=output,
    xlab=x_posns,
    display_names=sapply(phmObj$paramsList, function(x) x$class)[x_posns]
  )
}

#' Plot PHM Dendrogram (for v2)
#' 
#' @export 
plotPHMv2Dendrogram <- function(phmObj,
                             scaleHeights=c("log10", "unscaled", "pmcdist"),
                             colors=NULL,
                             displayAxis=c("box", "label", "index", "none"),
                             displayAxisSize=NULL,
                             colorAxis=NULL,
                             groupProbs=NULL,
                             groupColorMax="black",
                             groupColorMin="lightgray") {
  
  displayAxis <- match.arg(displayAxis)
  scaleHeights <- match.arg(scaleHeights)
  K <- ncol(phmObj$deltaPmc)

  phm_dendro_data <- constructVisData(phmObj,
                                      scaleHeights=scaleHeights,
                                      groupProbs=groupProbs)
  
  
  ## Visualization params
  if (is.null(colorAxis)) {
    colorAxis <- displayAxis == "box" || !is.null(colors)
  }
  
  ## Default is the paired pallette; only if K < 12
  if (colorAxis && is.null(colors)) {
    if (K > 12) {
      warning("More clusters than in default pallette. Suppressing coloring.")
      colors <- rep("#000000", K)
    } else {
      colors <- RColorBrewer::brewer.pal(K, "Paired")
    }
  } else if (colorAxis && K > length(colors)) {
    warning("Not enough colors provided. Suppressing coloring.")
    colors <- rep("#000000", K)
  } else if (!colorAxis) {
    colors <- rep("#000000", K)
  }
  
  ## Set axis theme
  if (is.null(displayAxisSize)) {
    if (displayAxis == "box") {
      displayAxisSize <- 8
    } else {
      displayAxisSize <- 10
    }
  }
  displayAxisFmt <- ggtext::element_markdown(
    color=unname(colors[phm_dendro_data$xlab]),
    size=displayAxisSize)
  displayAxisLabels <- if (displayAxis == "box") {
    rep("\U25A0", K)
  } else if (displayAxis == "label") {
    phm_dendro_data$display_names
  } else if (displayAxis == "index") {
    phm_dendro_data$xlab
  } else {
    NULL
  }
  
  
  scale_func <- ggplot2::scale_y_continuous
  
  plt <- ggplot2::ggplot(phm_dendro_data$df,
                         ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    # ggplot2::geom_segment(ggplot2::aes(color=gprob),
    #                       data=dplyr::filter(phm_dendro_data$df, linetype=="dashed"),
    #                       linetype="dashed") +
    # ggplot2::geom_segment(ggplot2::aes(color=gprob),
    #                       data=dplyr::filter(phm_dendro_data$df, linetype=="solid"),
    #                       linetype="solid") +
    ggplot2::geom_segment(ggplot2::aes(color=gprob)) +
    # xlab("Mixture Component ID") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(breaks=1:K,
                                labels=displayAxisLabels) +
    ggplot2::scale_color_gradient(low=groupColorMin, high=groupColorMax) +
    scale_func(expand=ggplot2::expansion(mult=c(0, 0.05))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x=displayAxisFmt,
                   axis.text.y=ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid=ggplot2::element_blank(),
                   panel.spacing=ggplot2::unit(0, "lines"),
                   panel.border=ggplot2::element_blank(),
                   legend.position="none"
    )
  
  plt
}

#' Plot PHM Matrix (for v2)
#' 
#' @export 
plotPHMv2Heatmap <- function(phmObj,
                          colors=NULL,
                          displayAxis=c("box", "label", "index", "none"),
                          displayAxisSize=NULL,
                          colorAxis=NULL,
                          gridColor="black",
                          fillLimits=NULL,
                          fillScale=c("log10", "pmcdist"),
                          legendPosition="none") {
  K <- ncol(phmObj$deltaPmc)
  displayAxis <- match.arg(displayAxis)
  fillScale <- match.arg(fillScale)
  
  if (is.null(colorAxis)) {
    colorAxis <- displayAxis == "box" || !is.null(colors)
  }
  
  ## Default is the paired pallette; only if K < 12
  if (colorAxis && is.null(colors)) {
    if (K > 12) {
      warning("More clusters than in default pallette. Suppressing coloring.")
      colors <- rep("#000000", K)
    } else {
      colors <- RColorBrewer::brewer.pal(K, "Paired")
    }
  } else if (colorAxis && K > length(colors)) {
    warning("Not enough colors provided. Suppressing coloring.")
    colors <- rep("#000000", K)
  } else if (!colorAxis) {
    colors <- rep("#000000", K)
  }
  
  ## Set axis theme
  if (is.null(displayAxisSize)) {
    if (displayAxis == "box") {
      displayAxisSize <- 8
    } else {
      displayAxisSize <- 10
    }
  }
  phm_dendro_data <- constructVisData(phmObj)
  displayAxisFmt <- ggtext::element_markdown(
    color=unname(colors[phm_dendro_data$xlab]),
    size=displayAxisSize)
  displayAxisLabels <- if (displayAxis == "box") {
    rep("\U25A0", K)
  } else if (displayAxis == "label") {
    phm_dendro_data$display_names
  } else if (displayAxis == "index") {
    phm_dendro_data$xlab
  } else {
    NULL
  }
  
  ## Construct the merging matrix
  label_map <- as.list(1:K)
  merge_matrix <- matrix(NA, nrow=K, ncol=K)
  for (x in K:2) {
    mc <- phmObj$mergeComps[[x]]
    grid <- expand.grid(label_map[[mc[1]]],
                        label_map[[mc[2]]])
    
    merge_matrix[grid[, 1], grid[, 2]] <- phmObj$mergeValues[x]
    merge_matrix[grid[, 2], grid[, 1]] <- phmObj$mergeValues[x]
    
    label_map[[mc[1]]] <- c(label_map[[mc[1]]], label_map[[mc[2]]])
    label_map[[mc[2]]] <- NULL
  }
  colnames(merge_matrix) <- paste0("V", 1:K)
  
  fillScaleFunc <- log10
  if (fillScale == "pmcdist") {
    load(system.file("extdata", "pmc_scale_function.RData", 
                     package = "PHM"))
    fillScaleFunc <- function(x) -inv_log10(log10(x))
  }
  
  matrix_long <- merge_matrix %>%
    dplyr::as_tibble(.name_repair = "unique") %>%
    tibble::rowid_to_column(var = "X") %>%
    tidyr::gather(key = "Y", value = "Z", -1) %>%
    dplyr::mutate(Y = as.numeric(gsub("V", "", Y)),
                  Z.old=Z,
                  Z=fillScaleFunc(Z))
  
  plot_lims <- fillLimits
  if (is.null(fillLimits)) {
    plot_lims <- c(min(matrix_long$Z.old, na.rm=T),
                   max(matrix_long$Z.old, na.rm=T))
  } else if (is.na(fillLimits[1])) {
    plot_lims[1] <- min(matrix_long$Z.old, na.rm=T)
  } else if (is.na(fillLimits[2])) {
    plot_lims[2] <- max(matrix_long$Z.old, na.rm=T)
  }
  
  if (fillScale == "log10") {
    plot_lims <- log10(plot_lims)
  } else {
    plot_lims <- -inv_log10(log10(plot_lims))
  }
  
  mid_point <- mean(plot_lims)
  
  matrix_long %>%
    dplyr::mutate(Z=Z,
                  Z.mod = Z,
                  Z.mod = ifelse(X == Y, "--", Z.mod)) %>%
    dplyr::mutate(X=factor(X, levels=phm_dendro_data$xlab, ordered=T),
                  Y=factor(Y, levels=phm_dendro_data$xlab, ordered=T)) %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, fill = Z)) +
    ggplot2::geom_tile(color = gridColor) +
    ggplot2::scale_fill_gradient2(limits = plot_lims,
                                  low = "blue",
                                  # mid="white",
                                  midpoint = mid_point,
                                  high = "red",
                                  aesthetics = "fill",
                                  breaks=c(plot_lims[1],
                                           plot_lims[2]),
                                  labels=c("Low", "High"),
                                  name="Similarity") +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::scale_x_discrete(labels=displayAxisLabels) +
    ggplot2::scale_y_discrete(labels=displayAxisLabels) +
    ggplot2::theme(legend.position = legendPosition,
                   panel.grid = ggplot2::element_blank(),
                   axis.text.x = displayAxisFmt,
                   axis.text.y = displayAxisFmt)
}
