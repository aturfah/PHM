constructPHMDendrogramData <- function(phm, uniformHeights=F, mergeLabels="delta", threshold=1e-3) {
  K <- length(phm)
  pmc <- phm[[K]]$pmc
  pmc_remains <- sapply(K:2, function(k) phm[[k]]$pmc)
  pmc_change <- sapply((K-1):1, function(k) phm[[k]]$pmc_change)
  height <- if (uniformHeights) {1:(K-1)} else {
    pmc / (pmc_remains)
  }
  merge_components <- t(sapply(K:2, function(k) phm[[k]]$merge_components))

  ## Figure out the components that merge together
  output <- data.frame()
  merge_tree <- lapply(1:K, function(k) k)
  height_tracker <- lapply(1:K, function(k) list(base=0, height=0))
  component_id_map <- 1:K
  base_height <- 0
  for (idx in 1:(K-1)) {
    ## Figure out which components merge
    mcs <- unname(merge_components[idx, ])
    hgt <- height[idx]

    ## Add this height to all components in height_tracker
    for (posn in 1:length(height_tracker)) height_tracker[[posn]]$height <- hgt # height_tracker[[posn]]$height + hgt

    ## For the components that were merged, add a vertical bar for them
    new_rows <- rbind(c(ID=component_id_map[mcs[1]], 
                        y=height_tracker[[mcs[1]]]$base,
                        yend=height_tracker[[mcs[1]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_pct=pmc_change[idx]/pmc,
                        pmc_change=pmc_change[idx]),
                      c(ID=component_id_map[mcs[2]],
                        y=height_tracker[[mcs[2]]]$base,
                        yend=height_tracker[[mcs[2]]]$height,
                        pmc=pmc_remains[idx],
                        pmc_pct=pmc_change[idx]/pmc,
                        pmc_change=pmc_change[idx]))
    if (nrow(output) == 0)  {
      output <- data.frame(new_rows)
    } else {
      output <- rbind(new_rows, output)
    }

    ## Make note of the combined merge nodes so we have the merges stored somewhere
    merge_tree[[mcs[1]]] <- list(
      merge_tree[[mcs[1]]],
      merge_tree[[mcs[2]]]
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
                    pmc_change=pmc_change.x) %>%
      dplyr::group_by(ID, y, yend, x, pmc, pmc_change, pmc_pct) %>%
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
  } else if (mergeLabels == "pct") {
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
#' @description TODO: Fill me in
#'
#' @param phm Output from [PHM()]
#' @param colors Vector of \eqn{K} hex codes to color the leaf node labels
#' @param uniformHeights Boolean whether the difference in heights of the merges should be constant or dependent on the reduction of \eqn{P_{\rm{mc}}}
#' @param threshold Error threshold for the integral past which to represent the merges as dashed lines
#' @param suppressLabels Boolean whether or not to display \eqn{P_{\rm{mc}}} reduction labels on the dendrogram or not
#' @param mergeLabels String indicating what value to display in the labels on the dendrogram
#' @param mergeLabelsSize Text size for the numeric value in the labels on the dendrogram
#' @param mergeLabelsPadding Padding for the merge label text
#' @param mergeLabelsR Radius for the rounded edges of the merge label box
#' @param displayAxis String indicating what label to place on the leaf nodes of the dendrogram
#' @param displayAxisSize Text size for the leaf node labels
#' @param colorAxis Whether or not to color the labels on the leaf nodes
#' 
#' @details TODO: Fill me in
#' 
#' @export 
plotPHMDendrogram <- function(phm, colors=NULL, 
                              uniformHeights=F, 
                              threshold=0,
                              suppressLabels=F,
                              mergeLabels=c("delta", "pmc", "percent"),
                              mergeLabelsSize=2,
                              mergeLabelsBorderSize=0.15,
                              mergeLabelsPadding=0.15,
                              mergeLabelsR=0.1,
                              displayAxis=c("box", "label", "index", "none"),
                              displayAxisSize=NULL,
                              colorAxis=NULL) {
  displayAxis <- match.arg(displayAxis)
  mergeLabels <- match.arg(mergeLabels)
  K <- length(phm)
  pmc_dendro_data <- constructPHMDendrogramData(phm, 
                                                uniformHeights = uniformHeights,
                                                mergeLabels=mergeLabels,
                                                threshold=threshold)

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

  offset <- 1
  scale_func <- ggplot2::scale_y_log10
  if (uniformHeights) {
    offset <- 0
    scale_func <- ggplot2::scale_y_continuous
  }
  
  plt <- ggplot2::ggplot(
      pmc_dendro_data$df, 
      ggplot2::aes(x=x, y=y+offset, xend=xend, yend=yend+offset)) +
    ggplot2::geom_segment(data=dplyr::filter(pmc_dendro_data$df, linetype=="dashed"), linetype="dashed") +
    ggplot2::geom_segment(data=dplyr::filter(pmc_dendro_data$df, linetype=="solid"), linetype="solid") +
    # xlab("Mixture Component ID") +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(breaks=1:K,
                                    labels=displayAxisLabels) +
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
                            ggplot2::aes(x=xposn, y=y+offset, label=lab),
                            size=mergeLabelsSize,
                            label.size=mergeLabelsBorderSize,
                            label.padding = ggplot2::unit(mergeLabelsPadding, "lines"),
                            label.r = ggplot2::unit(mergeLabelsR, "lines"))
  }

  plt
}

#' Generate the distruct plot from the posterior matrix
#'
#' @description TODO: Fill me in
#'
#' @param phm Output from the [PHM()] function
#' @param K Number of clusters for which to generate the distruct plot
#' @param labels Ground truth class labels for the observations (ordered factor vector)
#' @param colors Optinal vector with colors for the mixture components
#' @param axisTextSize Size for axis labels
#'
#' @details TODO: Fill me in
#' 
#' @export
plotPHMDistruct <- function(phm, K=length(phm),
                            colors=NULL,
                            labels=NULL,
                            axisTextSize=6) {

  ## Validation for colors
  if (is.null(colors) && K > 12) {
    stop("Too many clusters for default coloring scheme")
  } else if (is.null(colors)) {
    colors <- RColorBrewer::brewer.pal(K, "Set1")
  }

  ## Validation for labels
  post_mat <- phm[[K]]$posterior_matrix
  if (is.null(post_mat)) stop("Does phm have a posterior matrix computed? Run addPosteriorMatrix() to compute.")

  if (is.null(labels)) {
    if (is.null(phm[[K]]$labels)) stop("Must either provide labels or have them computed. Run addPosteriorMatrix() to compute.")
    labels <- phm[[K]]$labels
  }
  if (!is.factor(labels)) labels <- factor(labels)
  labels_order=levels(labels)

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

#' Plot \eqn{\Delta \P_{\rm{mc}}} matrix
#'
#' @description TODO: Fill  me in
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
#' @details TODO: Fill me in
#' 
#' @return pew
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