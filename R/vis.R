constructPHMDendrogramData <- function(phm,
                                       scaleHeights="pmcdist",
                                       heightValue="merge",
                                       mergeLabels="delta",
                                       threshold=1e-3,
                                       groupProbs=NULL) {
  K <- length(phm)
  gprobs_unspec <- FALSE
  if (is.null(groupProbs)) {
    gprobs_unspec <- TRUE
    groupProbs <- rep(1, K)
  }

  pmc <- phm[[K]]$pmc
  pmc_remains <- sapply(K:2, function(k) phm[[k]]$pmc)
  pmc_change <- sapply((K-1):1, function(k) phm[[k]]$pmc_change)
  pmc_components <- sapply((K-1):1, function(k) phm[[k]]$pmc_components)
  pmc_min <- sapply((K-1):1, function(k) phm[[k]]$min_merge_pmc)

  height <- if (heightValue == "min") {
    pmc_min
  } else {
    pmc_change / choose(pmc_components, 2)
    # pmc_change / 2^pmc_components 
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
  merge_components <- t(sapply(K:2, function(k) phm[[k]]$merge_components))

  ## Figure out the components that merge together
  groupProbs_new <- groupProbs
  output <- data.frame()
  merge_tree <- lapply(1:K, function(k) k)
  height_tracker <- lapply(1:K, function(k) list(base=0, height=0))
  component_id_map <- 1:K
  base_height <- 0
  for (idx in 1:(K-1)) {
    ## Figure out which components merge
    mcs <- unname(merge_components[idx, ])
    hgt <- height[idx]

    ## Verify that height doesn't go down
    for (posn in mcs) {
      hgt <- max(hgt, height_tracker[[posn]]$base + 1e-6)
    }

    ## Add this height to all components in height_tracker
    for (posn in 1:length(height_tracker)) height_tracker[[posn]]$height <- hgt

    ## For the components that were merged, add a vertical bar for them
    new_rows <- rbind(c(ID=component_id_map[mcs[1]],
                        y=height_tracker[[mcs[1]]]$base,
                        yend=height_tracker[[mcs[1]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_pct=pmc_change[idx]/pmc,
                        pmc_change=pmc_change[idx],
                        gprob=groupProbs_new[mcs[1]] ),
                      c(ID=component_id_map[mcs[2]],
                        y=height_tracker[[mcs[2]]]$base,
                        yend=height_tracker[[mcs[2]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_pct=pmc_change[idx]/pmc,
                        pmc_change=pmc_change[idx],
                        gprob=groupProbs_new[mcs[2]]))
    if (nrow(output) == 0)  {
      output <- data.frame(new_rows)
    } else {
      output <- rbind(new_rows, output)
    }

    ## Combine the group probs based on parameters
    if (gprobs_unspec) {
      groupProbs_new <- rep(1, length(groupProbs_new) - 1)
    } else {
      if (is.null(phm[[K - idx + 1]]$params)) {
        phm[[K - idx + 1]]$params <- recoverPHMParams(phm, K - idx + 1, "prob")
      }
      probs1 <- sum(phm[[K - idx + 1]]$params[[mcs[1]]]$prob)
      probs2 <- sum(phm[[K - idx + 1]]$params[[mcs[2]]]$prob)
      groupProbs_new[mcs[1]] <- (groupProbs_new[mcs[1]] * probs1 +
          groupProbs_new[mcs[2]] * probs2) / (probs1 + probs2)
      groupProbs_new <- groupProbs_new[-mcs[2]]
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

  ## Figure out where everything goes relative to base node
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
  output <- dplyr::mutate(output, x=ifelse(y==0, map_xposns(ID), NA))

  while(any(is.na(output))) {
    output <- output %>%
      dplyr::left_join(output, by=c("y"="yend")) %>%
      dplyr::rename(x=x.x,
                    ID=ID.x,
                    pmc=pmc.x,
                    pmc_pct=pmc_pct.x,
                    pmc_change=pmc_change.x,
                    gprob=gprob.x) %>%
      dplyr::group_by(ID, y, yend, x, pmc, pmc_change, pmc_pct, gprob) %>%
      dplyr::summarize(x=ifelse(all(is.na(x)), mean(x.y), mean(x)),
                xend=x, .groups="keep") %>%
      dplyr::ungroup() %>%
      dplyr::arrange(-yend) %>%
      dplyr::select(-dplyr::ends_with(".y"))
  }

  ## Add horizontal components
  horiz_comps <- output %>%
    dplyr::group_by(yend, pmc, pmc_change, pmc_pct) %>%
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

  output <-  dplyr::mutate(output, linetype=ifelse(pmc_change < threshold, "dashed", "solid"))

  output$label_column <- if (mergeLabels == "delta") {
    output$pmc_change
  } else if (mergeLabels == "percent") {
    output$pmc_pct
  } else {
    output$pmc
  }

  labels <- output %>%
    dplyr::filter(is.na(ID)) %>%
    dplyr::mutate(xposn=(x+xend)/2,
           lab=ifelse(
             pmc_change < 1e-3,
             formatC(label_column, format = "e", digits = 2),
             round(label_column, 4)
           ))

  list(
    df=output,
    xlab=x_posns,
    labels=labels,
    display_names=sapply(phm[[K]]$params, function(x) x$class)[x_posns]
  )
}

#' Visualize PHM merging procedure via Dendrogram
#'
#' @description Visualize the PHM merging procedure using a dendrogram. Visualization options, such as displaying the merge \eqn{\Delta P_{\rm mc}} value and tracking group membership across merges is included.
#'
#' @param phm Output from [PHM()]
#' @param colors Vector of \eqn{K} hex codes to color the leaf node labels
#' @param scaleHeights String specifying how to set the heights in the dendrogram. See Details for more information.
#' @param heightValue Whether to use the \eqn{\Delta P_{\rm mc}} or \eqn{\min \Delta P_{\rm mc}} value to determine branch height for a merge. See Details for more information.
#' @param threshold Error threshold for the integral past which to represent the merges as dashed lines
#' @param suppressLabels Boolean whether or not to display \eqn{P_{\rm{mc}}} reduction labels on the dendrogram or not
#' @param mergeLabels String indicating what value to display in the labels on the dendrogram
#' @param mergeLabelsSize Text size for the numeric value in the labels on the dendrogram
#' @param mergeLabelsPadding Padding for the merge label text
#' @param mergeLabelsR Radius for the rounded edges of the merge label box
#' @param displayAxis String indicating what label to place on the leaf nodes of the dendrogram. See Details for more information.
#' @param displayAxisSize Text size for the leaf node labels
#' @param colorAxis Whether or not to color the labels on the leaf nodes
#' @param groupProbs Vector of class probability conditional on base group membership
#' @param groupColorMax Color of lines corresponding to high group probability
#' @param groupColorMin Color of lines corresponding to low group probability
#'
#' @details 
#' 
#' There are two options for the value for the \eqn{P_{\rm mc}} height for merging subtrees \eqn{\mathcal{S}_1, \mathcal{S}_2} (\eqn{m = |\mathcal{S}_1| + |\mathcal{S}_2|}):
#' \itemize{
#'  \item \code{"merge"}: \eqn{\binom{m}{2}^{-1} \Delta P_{\rm mc}^{(\mathcal{S}_1, \mathcal{S}2)}}
#'  \item \code{"min"}: \eqn{\min_{i \in \mathcal{S}_1,\, j \in \mathcal{S}_2} \Delta P_{\rm mc}^{(i, j)}}
#' }
#' 
#' Once the value for the tree height is obtained, there are three options to scale.
#' \itemize{
#'  \item \code{"unscaled"}: The height can be left unscaled
#'  \item \code{"log10"}: A \eqn{\log_{10}} scaling can be applied to the height to better reveal different clustering resolutions.
#'  \item \code{"pmcdist"}: A spline-based scaling where we map the value to a linear distance between two Gaussian clusters.
#' }
#' 
#' The \code{displayAxis} parameter controls what to display on the axes of the heatmap. \code{box} displays a standard box, which is most useful for color-coded axes. \code{label} displays the cluster label, as specified in the \code{phm} parameters. \code{index} displays the numeric index of each cluster in the \code{phm} parameter list. \code{none} suppresses the cluster label entirely, and is mose useful for when there are a large number of clusters.
#'
#' @return A ggplot object
#' @export
plotPHMDendrogram <- function(phm, colors=NULL,
                              scaleHeights=c("log10", "unscaled", "pmcdist"),
                              heightValue=c("merge", "min"),
                              threshold=0,
                              suppressLabels=F,
                              mergeLabels=c("delta", "pmc", "percent"),
                              mergeLabelsSize=2,
                              mergeLabelsBorderSize=0.15,
                              mergeLabelsPadding=0.15,
                              mergeLabelsR=0.1,
                              displayAxis=c("box", "label", "index", "none"),
                              displayAxisSize=NULL,
                              colorAxis=NULL,
                              groupProbs=NULL,
                              groupColorMax="black",
                              groupColorMin="lightgray") {
  scaleHeights = match.arg(scaleHeights)
  displayAxis <- match.arg(displayAxis)
  mergeLabels <- match.arg(mergeLabels)
  heightValue <- match.arg(heightValue)

  K <- length(phm)

  if (!is.null(groupProbs) && length(groupProbs) != K) {
    stop("If specified, groupProbs must have same length as phm")
  }
  if (!is.null(groupProbs)) {
    groupProbs <- unname(groupProbs)
  }

  pmc_dendro_data <- constructPHMDendrogramData(phm,
                                                scaleHeights = scaleHeights,
                                                heightValue = heightValue,
                                                mergeLabels=mergeLabels,
                                                threshold=threshold,
                                                groupProbs=groupProbs)

  ## By default, only color with "box" axis display; otherwise no color
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
        color=unname(colors[pmc_dendro_data$xlab]),
        size=displayAxisSize)
  displayAxisLabels <- if (displayAxis == "box") {
    rep("\U25A0", K)
  } else if (displayAxis == "label") {
    pmc_dendro_data$display_names
  } else if (displayAxis == "index") {
    pmc_dendro_data$xlab
  } else {
    NULL
  }


  scale_func <- ggplot2::scale_y_continuous # ggplot2::scale_y_log10
  plt <- ggplot2::ggplot(pmc_dendro_data$df,
                         ggplot2::aes(x=x, y=y, xend=xend, yend=yend)) +
    ggplot2::geom_segment(ggplot2::aes(color=gprob),
                          data=dplyr::filter(pmc_dendro_data$df, linetype=="dashed"),
                          linetype="dashed") +
    ggplot2::geom_segment(ggplot2::aes(color=gprob),
                          data=dplyr::filter(pmc_dendro_data$df, linetype=="solid"),
                          linetype="solid") +
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

  if (!suppressLabels) {
    plt <- plt + ggplot2::geom_label(data=pmc_dendro_data$labels,
                            ggplot2::aes(x=xposn, y=y, label=lab, alpha=NULL),
                            size=mergeLabelsSize,
                            label.size=mergeLabelsBorderSize,
                            label.padding = ggplot2::unit(mergeLabelsPadding, "lines"),
                            label.r = ggplot2::unit(mergeLabelsR, "lines"))
  }

  plt
}

#' Generate the distruct plot from the posterior matrix
#'
#' @description Visualize a distruct plot based on either a partition or the posterior cluster probabilities.
#'
#' @param phm Output from the [PHM()] function
#' @param K Number of clusters for which to generate the distruct plot
#' @param labels Ground truth class labels for the observations (ordered factor vector)
#' @param colors Optinal vector with colors for the mixture components
#' @param axisTextSize Size for axis labels
#' @param partition Whether to visualize from the posterior matrix or partition labels
#'
#' @details 
#' In the case of visualizing for a partition, the posterior probabilities are set to 1 if it is the cluster the obesrvation is assigned to and 0 otherwise.
#'
#' @return A ggplot object
#' @export
plotPHMDistruct <- function(phm, K=length(phm),
                            colors=NULL,
                            labels=NULL,
                            axisTextSize=6,
                            partition=F) {

  ## Validation for colors
  if (is.null(colors) && K > 12) {
    stop("Too many clusters for default coloring scheme")
  } else if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(K, "Set1")
  }

  ## Validation for labels
  if (partition && (K == length(phm))) {
    ## Onehot encode the cluster labels to get posterior
    clust_labels <- phm[[K]]$labels
    clust_levels <- sapply(phm[[K]]$params, function(x) x$class)
    clust_labels <- factor(clust_labels, levels=clust_levels)
    post_mat <- sapply(clust_labels, function(idx) {
      v <- numeric(K)
      v[idx] <- 1
      v
    })
    post_mat <- t(post_mat)
  } else {
    post_mat <- phm[[K]]$posterior_matrix
    if (is.null(post_mat)) stop("Does phm have a posterior matrix computed? Run addPosteriorMatrix() to compute.")
  }

  if (is.null(labels)) {
    if (is.null(phm[[K]]$labels)) stop("Must either provide labels or have them computed. Run addPosteriorMatrix() to compute.")
    labels <- phm[[K]]$labels
  }
  if (!is.factor(labels)) labels <- factor(labels)
  labels_order <- levels(labels)

  ## Define break points for the distruct plot
  group_counts <- table(labels)
  label_positions <- sapply(1:length(group_counts), function(idx) {
    ifelse(idx != 1, sum(group_counts[1:(idx-1)]), 0) + group_counts[idx] / 2
  })

  ## Form data for plotting
  df <- data.frame(labels=labels, posterior=post_mat)

  group_columns_df <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(max_column = which.max(dplyr::c_across(dplyr::where(is.numeric)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::starts_with("posterior")) %>%
    dplyr::group_by(labels, max_column) %>%
    dplyr::summarize(count=dplyr::n()) %>%
    dplyr::group_by(labels) %>%
    dplyr::filter(count == max(count))

  group_columns <- dplyr::pull(group_columns_df, "max_column")
  names(group_columns) <- dplyr::pull(group_columns_df, "labels")

  plt_data <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sort_val = dplyr::c_across(dplyr::where(is.numeric))[group_columns[labels]]) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(labels, -sort_val) %>%
    dplyr::mutate(row=1:dplyr::n()) %>%
    tidyr::pivot_longer(cols=dplyr::starts_with("posterior")) %>%
    dplyr::arrange(row, value)

  ## Positions for vertical lines
  vert_lines <- plt_data %>%
    dplyr::group_by(labels) %>%
    dplyr::summarize(minX=min(row),
              maxX=max(row)) %>%
    dplyr::mutate(minX=minX-0.5, maxX=maxX+0.5,
           y=0, yend=1) %>%
    tidyr::pivot_longer(cols=c(minX, maxX))

  ## Distruct plot
  plt <- plt_data %>%
    ggplot2::ggplot(ggplot2::aes(x=row, y=value, fill=name)) +
    ggplot2::geom_bar(position="stack", stat="identity") +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::geom_segment(data=vert_lines,
                 aes(x=value, xend=value, y=y, yend=yend, fill=NULL),
                 linetype="dashed", linewidth=0.3) +
    ggplot2::scale_fill_manual(values=colors) +
    ggplot2::scale_x_continuous(expand=c(0, 0),
                       labels=labels_order,
                       breaks=label_positions,
                       minor_breaks = NULL) +
    ggplot2::scale_y_continuous(expand=c(0, 0),
                       breaks=NULL, minor_breaks = NULL, n.breaks=0) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x=ggplot2::element_text(hjust=0.5, size=axisTextSize),
        axis.ticks.y=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        panel.grid=ggplot2::element_blank(),
        legend.position="none",
        plot.margin = ggplot2::unit(c(0.01, 0.03, 0, 0.01), "inches"))

  return(plt)
}

#' Plot \eqn{\Delta P_{\rm{mc}}} matrix
#'
#' @description Visualize the matrix of \eqn{\Delta P_{\rm mc}} values
#'
#' @param phm Output from [PHM()]
#' @param K Number of clusters for which to visualize the heatmap
#' @param colors Vector of \eqn{K} hex codes to color the leaf node labels
#' @param displayAxis String indicating what label to place along the axis (corresponding to clusters)
#' @param displayAxisSize Text size for the axis labels
#' @param colorAxis Whether or not to color the axis labels
#' @param visScale Whether to display the raw \eqn{\Delta P_{\rm{mc}}} values or scale them to be percent of total \eqn{P_{\rm{mc}}}
#' @param visSize Text size for the values inside the heatmap
#' @param visThreshold At what value suppress the value and show "< (visThreshold)"
#' @param visDigits Number of digits to round the displayed values
#'
#' @return A ggplot object
#' @export
plotPmcMatrix <- function(phm, K=length(phm), colors=NULL,
                          displayAxis=c("box", "label", "index", "none"),
                          displayAxisSize=NULL,
                          colorAxis=NULL,
                          visScale=c("absolute", "percent"),
                          visSize=2,
                          visThreshold=1e-3,
                          visDigits=3) {
  visScale <- match.arg(visScale)
  displayAxis = match.arg(displayAxis)

  if (is.null(colorAxis)) {
    colorAxis <- displayAxis == "box" || !is.null(colors)
  }

  ## Format Colors
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

  ## Format axis
  if (is.null(displayAxisSize)) {
    if (displayAxis == "box") {
        displayAxisSize <- 8
    } else {
        displayAxisSize <- 10
    }
  }
  if (displayAxis == "box") {
    displayAxisLabels <- rep("\U25A0", K)
  } else if (displayAxis == "index") {
    displayAxisLabels <- 1:K
  } else if (displayAxis == "label") {
    displayAxisLabels <- sapply(phm[[K]]$params, function(x) x$class)
  } else {
    displayAxisLabels <- NULL
  }

  displayAxisFmt <- ggtext::element_markdown(
        color=unname(colors),
        size=displayAxisSize)


  ## Format Threshold
  fmt_func <- function(x) round(x, visDigits)

  pmc_matrix <- phm[[K]]$pmc_matrix
  pmc <- sum(pmc_matrix)
  plotLimits <- if (visScale == "percent") 1 else pmc

  colnames(pmc_matrix) <- paste0("V", 1:K)

  pmc_matrix_long <- (2 * pmc_matrix) %>%
    # Data wrangling
    dplyr::as_tibble(.name_repair="unique") %>%
    tibble::rowid_to_column(var="X") %>%
    tidyr::gather(key="Y", value="Z", -1) %>%
    dplyr::mutate(Y=as.numeric(gsub("V","",Y)))

  if (visScale == "percent") pmc_matrix_long$Z <- pmc_matrix_long$Z / pmc

  htmp <- pmc_matrix_long %>%
    dplyr::mutate(Z.mod = ifelse(
        Z < visThreshold,
        paste0("< ", visThreshold),
        fmt_func(Z)
      ),
      Z.mod = ifelse(X == Y, "--", Z.mod)) %>%
    dplyr::filter(X <= Y) %>%
    # Viz
    ggplot2::ggplot(aes(X, Y, fill=Z)) +
    ggplot2::geom_tile(color="black") +
    ggplot2::geom_text(aes(label=Z.mod, fill=NULL),
              size=visSize) +
    ggplot2::scale_fill_gradient(limits=c(0, plotLimits),
                        low="white", high="red",
                        aesthetics="fill") +
    ## Format the axis labels
    ggplot2::scale_x_continuous(breaks=1:K,
                       minor_breaks = NULL,
                       labels=displayAxisLabels) +
    ggplot2::scale_y_continuous(breaks=1:K,
                       minor_breaks = NULL,
                       labels=displayAxisLabels) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme(
          legend.position="NULL",
          panel.grid=ggplot2::element_blank(),
          axis.text.x=displayAxisFmt,
          axis.text.y=displayAxisFmt)

  return(htmp)
}


#' Visualize PHM dendrogram structure with a heatmap
#'
#' @description 
#' For a pair of clusters \eqn{i, j}, the heatmap position \eqn{i, j} is the value of \eqn{\Delta P_{\rm mc}} where the clusters are first merged. This allows for a more straightforward visualization of the multi-resolution structure of heatmaps when combined with the dendrogram.
#'
#' @param phm Output from [PHM()]
#' @param colors Vector of \eqn{K} hex codes to color the leaf node labels
#' @param displayAxis String indicating what label to place along the axis (corresponding to clusters). See Details for more information.
#' @param displayAxisSize Text size for the axis labels
#' @param colorAxis Whether or not to color the axis labels
#' @param gridColor What color to make the heatmap grid
#' @param fillLimits Optional vector to manually set limits of the fill scaling. Default is to use the min and max \eqn{\Delta P_{\rm mc}} values from \code{phm}
#' @param fillScale Whether to use \eqn{\log_{10} \Delta P_{\rm mc}} or the spline scaling for the heatmap color
#' @param legendPosition Where to put the legend for the heatmap colors. Default is to suppress.
#'
#' @details
#' Consider a pair of initial clusters \eqn{i, j}, and let \eqn{\mathcal{S}} be the smallest subtree containing both \eqn{i, j} (\eqn{\mathcal{S}_i} contains cluster \eqn{i} and \eqn{\mathcal{S}_{j}} contains cluster \eqn{j}). The pairwise value for \eqn{i, j}, \eqn{\rho(i, j)} is based on \eqn{\Delta P_{\rm mc}^{ \left(\mathcal{S}(i), \mathcal{S}(j)\right)}}, which is the \eqn{\Delta P_{\rm mc}} value at which \eqn{i, j} are merged into a single cluster. 
#' 
#' As with the dendrogram heights ([plotPHMDendrogram()]) there are multiple options for scaling \eqn{\rho(i, j)}. The first is \code{log10} scaling, where \eqn{\rho(i, j) = \log_{10} \Delta P_{\rm mc}}. The second is a spline-based scaling, \code{pmcdist}, where we map the \eqn{\Delta P_{\rm mc}} value to a linear distance between two Gaussian clusters.
#' 
#' The \code{displayAxis} parameter controls what to display on the axes of the heatmap. \code{box} displays a standard box, which is most useful for color-coded axes. \code{label} displays the cluster label, as specified in the \code{phm} parameters. \code{index} displays the numeric index of each cluster in the \code{phm} parameter list. \code{none} suppresses the cluster label entirely, and is mose useful for when there are a large number of clusters.
#'
#' @return A ggplot object
#' @export
plotPHMMatrix <- function(phm, colors=NULL,
                          displayAxis=c("box", "label", "index", "none"),
                          displayAxisSize=NULL,
                          colorAxis=NULL,
                          gridColor="black",
                          fillLimits=NULL,
                          # fillValue=c("deltaPmc", "minDeltaPmc", "scaledDeltaPmc"),
                          fillScale=c("log10", "pmcdist"),
                          legendPosition="none") {
  displayAxis <- match.arg(displayAxis)
  fillScale <- match.arg(fillScale)
  K <- length(phm)

  ## Same preprocessing as plotPmcMatrix
  ## By default, only color with "box" axis display; otherwise no color
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
  pmc_dendro_data <- constructPHMDendrogramData(phm)
  displayAxisFmt <- ggtext::element_markdown(
        color=unname(colors[pmc_dendro_data$xlab]),
        size=displayAxisSize)
  displayAxisLabels <- if (displayAxis == "box") {
    rep("\U25A0", K)
  } else if (displayAxis == "label") {
    pmc_dendro_data$display_names
  } else if (displayAxis == "index") {
    pmc_dendro_data$xlab
  } else {
    NULL
  }

  ## Construct Merging Matrix
  label_map <- as.list(seq_along(phm))
  merge_matrix <- matrix(NA, nrow=length(phm), ncol=length(phm))
  for (x in rev(seq_along(phm))) {
    mc <- phm[[x]]$merge_components

    if (sum(mc) < 0 || x == 1) next

    grid <- expand.grid(label_map[[mc["row"]]],
                        label_map[[mc["col"]]])

    merge_matrix[grid[, 1], grid[, 2]] <- phm[[x-1]]$pmc_change
    merge_matrix[grid[, 2], grid[, 1]] <- phm[[x-1]]$pmc_change

    ## Update
    label_map[[mc["row"]]] <- c(label_map[[mc["row"]]],
                                label_map[[mc["col"]]])
    label_map[[mc["col"]]] <- NULL
  }
  colnames(merge_matrix) <- paste0("V", seq_along(phm))

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

  ## Can be cleaned up, vestigial code from Pmc Matrix
  matrix_long %>%
    dplyr::mutate(Z=Z,
                  Z.mod = Z,
                  Z.mod = ifelse(X == Y, "--", Z.mod)) %>%
    dplyr::mutate(X=factor(X, levels=pmc_dendro_data$xlab, ordered=T),
                  Y=factor(Y, levels=pmc_dendro_data$xlab, ordered=T)) %>%
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


#' Visualize regions contributing to \eqn{P_{\rm mc}} in a 2D plot
#'
#' @description Visualize the point-specific \eqn{P_{\rm mc}} over a grid of points to visually inspect cluster contributions to \eqn{P_{\rm mc}}
#' 
#'
#' @param paramsList A list generated from [constructPmcParamsMclust()], [constructPmcParamsPartition()], [constructPmcParamsPartition()] providing the initial cluster parameter estimates. Parameters should be 2D
#' @param data A \eqn{N \times 2} matrix containing observations
#' @param partition A vector containing class labels
#' @param colors Vector of \eqn{K} hex codes to color the leaf node labels
#' @param xlim Vector with min/max values for the x-axis `NULL` sets this based on the maximum and minimum x-values for observations
#' @param ylim Vector with min/max values for the y-axis `NULL` sets this based on the maximum and minimum y-values for observations
#' @param suppressPmc Flag whether to display the regions contributing to \eqn{P_{\rm mc}} in the plot
#' @param numPmcPatches Number of points along each axis in the \eqn{P_{mc}} grid to evaluate. Higher values gives a smoother visualization
#' @param logPmcThreshold Threshold for \eqn{P_{\rm mc}} values to display. All coordinates with \eqn{\log_{10} P_{\rm mc}} below this value will be ignored.
#' @param logPmcMidpoint Midpoint in the \eqn{P_{\rm mc}} color gradient
#' @param PmcColor Hex value for the maximum value in the \eqn{P_{\rm mc}} color gradient
#' @param suppressDensity Flag whether to overlay the cluster-specific densities
#' @param densityLevels Values of the density at which to display the level curves
#' @param densityLevelWidth Line width for the density curves
#' @param suppressObservations Flag whether to display the observations from \code{data}
#' @param pointSize Size of the data scatterplot points
#' @param labelSize Text size for the cluster labels
#' @param textSize Text size for the plots
#' @param legendPosition Where to put the legend for the heatmap colors. Default is to suppress
#'
#' @details
#' For an arbitrary point \eqn{\mathbf x}, using the Monte Carlo evaluation of \eqn{P_{\rm mc}} its point-specific \eqn{P_{\rm mc}} can be calculated as
#' \deqn{\sum_{k=1}^K \pi_k(\mathbf{x}) (1 - \pi_k(\mathbf{x})) \times f(x=\mathbf{x}) }
#' Where \eqn{f(\mathbf{x})} is the overall data density evaluated at a point. Note that for a sample of \eqn{M} points, the average of these point-specific \eqn{P_{\rm mc}} values produce the Monte Carlo estimate of \eqn{P_{\rm mc}} described in Turfah and Wen (2025).
#' 
#' The point-specific \eqn{P_{\rm mc}} of each point in a grid is evaluated and visualized. The cluster-specific density level sets and/or the partitioned observations can be overlaid over this to better understand how the \eqn{P_{\rm mc}} value was obtained.
#' 
#' @return A ggplot object
#' @export
plotPmc2D <- function(paramsList,
                      data,
                      partition,
                      colors=RColorBrewer::brewer.pal(length(paramsList), "Paired"),
                      xlim=NULL,
                      ylim=NULL,
                      suppressPmc=F,
                      numPmcPatches=200,
                      logPmcThreshold=-4,
                      logPmcMidpoint=0.6 * logPmcThreshold,
                      PmcColor="#222",
                      suppressDensity=F,
                      densityLevels=c(1e-2, 1e-1),
                      densityLevelWidth=0.3,
                      suppressObservations=F,
                      pointSize=0.5,
                      labelSize=2,
                      textSize=9,
                      legendPosition="none"
                      ) {
  ## Validate parameter dimensions
  if (length(partition) != nrow(data)) stop("Observations mismatch between partition and data")
  if (ncol(data) != 2) stop("Data must be of dimension 2")
  if (!is.null(xlim) && length(xlim) != 2) stop("xlim must be of length 2 or NULL")
  if (!is.null(ylim) && length(ylim) != 2) stop("ylim must be of length 2 or NULL")

  ## No color specified => black
  if (is.null(colors)) colors <- rep("#000000", length(paramsList))

  ## Define plot limits
  if (is.null(xlim)) {
    xlim <- c(min(data[, 1]), max(data[, 1]))
    if (xlim[1] < 0) {
      xlim[1] <- xlim[1] * 1.1
    } else {
      xlim[1] <- xlim[1] / 1.1
    }
    if (xlim[2] < 0) {
      xlim[2] <- xlim[2] / 1.1
    } else {
      xlim[2] <- xlim[2] * 1.1
    }
  }
  if (is.null(ylim)) {
    ylim <- c(min(data[, 2]), max(data[, 2]))
    if (ylim[1] < 0) {
      ylim[1] <- ylim[1] * 1.1
    } else {
      ylim[1] <- ylim[1] / 1.1
    }
    if (ylim[2] < 0) {
      ylim[2] <- ylim[2] / 1.1
    } else {
      ylim[2] <- ylim[2] * 1.1
    }
  }

  data_df <- data.frame(
    data, 
    g=as.factor(partition)
  )

  ## Construct Grid of points at which to evaluate Pmc
  mat <- expand.grid(X=seq(min(xlim), max(xlim), length.out=numPmcPatches),
                     Y=seq(min(ylim), max(ylim), length.out=numPmcPatches)) %>%
    as.matrix()

  ## Evaluate Pmc
  density_mat <- sapply(paramsList, function(x) {
    K <- length(x$prob)

    tmp <- sapply(1:K, function(idx) {
      x$prob[idx] * mvtnorm::dmvnorm(mat, x$mean[, idx], x$var[, , idx]) / sum(x$prob)
    })
    rowSums(tmp)
  })
  posterior <- density_mat / rowSums(density_mat)
  pmc <- rowSums(posterior * (1 - posterior)) * rowSums(density_mat)

  ## Cluster Density Matrix for visualization
  class_labels <- sapply(paramsList, function(x) x$class)
  names(colors) <- class_labels
  dens_df <- data.frame(mat, dens=density_mat) %>%
    tidyr::pivot_longer(cols=starts_with("dens")) %>%
    dplyr::mutate(name=stringr::str_remove(name, "dens."),
           name=as.numeric(name),
           name=class_labels[name],
           name=factor(name, levels=class_labels))

  ## Generate plot
  plt <- data.frame(mat, pmc=log10(pmc)) %>%
    dplyr::filter(pmc > logPmcThreshold) %>%
    ggplot2::ggplot()

  if (!suppressPmc) {
    labels_func <- function(x) {
      x_fmt <- format(round(x, 1), trim = TRUE)
      x_fmt[x == logPmcThreshold] <- paste0("≤ ", logPmcThreshold)
      x_fmt
    }
    plt <- plt + ggplot2::geom_point(aes(x=X, y=Y, color=pmc)) +
      ggplot2::scale_color_gradient2(midpoint=logPmcMidpoint,
                            high=PmcColor, low="white",
                            limits=c(logPmcThreshold,
                                    max(log10(pmc))),
                            na.value="white",
                            name="log(Pmc)",
                            labels=labels_func) +
      ggnewscale::new_scale_color()
  }
  if (!suppressObservations) {
    plt <- plt + ggplot2::geom_point(
        aes(x=X1,
            y=X2,
            color=g),
        size=pointSize,
        alpha=0.5,
        data = data_df
      ) +
      ggplot2::scale_color_manual(values=colors,
                                  guide="none") +
      ggnewscale::new_scale_color()
  }
  if (!suppressDensity) {
    plt <- plt + ggplot2::geom_contour(aes(x=X, y=Y, color=name, group=name, z=value), 
                 breaks = densityLevels,
                 linewidth=densityLevelWidth, data=dens_df) +
      ggplot2::scale_color_manual(values=colors,
                                  guide="none")
  }
  plt <- plt + 
    ## Cluster labels
    ggplot2::geom_label(
      aes(x=X1, y=X2, label=g),
      size=labelSize,
      data=data_df %>% group_by(g) %>% dplyr::summarize(X1=mean(X1), X2=mean(X2))
    ) +
    ## Formatting
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

  return(plt)
}