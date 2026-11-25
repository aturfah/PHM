## export_phm_json.R
##
## Generic exporter: turns ANY PHMv2() result into the JSON schema consumed
## by the "PHM Explorer" interactive app. Always exports the FULL fidelity
## merge tree (every fitted component as a leaf, every merge) -- the app
## itself handles collapsing that down to a smaller number of leaves/blocks
## for display, live, so this script never needs to know about a display
## threshold.
##
## Usage:
##   source("export_phm_json.R")
##
##   ## fit is already 2D -- X doubles as both the density input and the plot
##   export_phm_json(phm2, X = X, out_path = "my_result.json", label = "my analysis")
##
##   ## fit is higher-dimensional -- X (full-dimensional) still drives the
##   ## posterior shading, but the scatter is laid out with a separate 2D
##   ## embedding (UMAP, t-SNE, PCA, ...) that has nothing to do with the fit
##   embedding <- umap::umap(X)$layout
##   export_phm_json(phm2, X = X, xy = embedding, out_path = "my_result.json")
##
##   ## no access to the original data at all, only someone else's posteriors
##   ## and a layout -- both supplied directly, X isn't needed
##   export_phm_json(phm2, xy = embedding, z = my_posteriors, out_path = "my_result.json")
##
## Arguments:
##   phm2      A PHMv2() result (required). Works for a single-fit result or
##             a pooled/ensemble paramsList (e.g. from
##             constructPmcParamsSubAggMclust) -- nothing here assumes an
##             mclust object exists.
##   X         Optional original-space data, matrix/data.frame, used ONLY to
##             compute posterior probabilities when `z` isn't supplied. Its
##             column count must match the fit's own dimensionality (every
##             column PHMv2 was fit on) -- this is deliberately NOT truncated
##             to 2 columns for the posteriors. When `xy` isn't given, X[, 1:2]
##             also serves as the scatter layout (see `xy`).
##   xy        Optional n x 2 coordinates for the scatter panel specifically
##             -- what the point actually LOOKS like on screen. Any 2D
##             layout works: the fit's own two dimensions, a UMAP or t-SNE
##             embedding computed separately, a PCA projection, and so on.
##             It's never used for anything but placing the dots. If omitted
##             and `X` is given, the first two columns of `X` are used
##             (xy = X[, 1:2]) -- for a 2D fit that is the fitted space itself.
##   z         Optional n x K posterior-probability matrix (K = number of
##             fitted components), e.g. mcl$z for a single Mclust fit. If X
##             is supplied but z is not, it is computed generically from
##             phm2$paramsList (Gaussian mixture posterior via Cholesky
##             log-density + softmax), so this works even when there is no
##             mclust object at all (pooled/ensemble fits included). Supply
##             this directly when you have posteriors from elsewhere (or a
##             non-Gaussian fit) -- X is then not needed at all.
##   annotation  Optional length-n factor/character/logical/integer vector of
##             categorical labels, one per row of the layout (`xy`, or `X` when
##             it is already 2D). Adds a second, static scatter panel in which
##             every category gets its own color. It is display-only: it plays
##             no part in the fit or in the posterior shading, and it works even
##             without `X`/`z` (then only that static panel is shown). Factor
##             level order is kept as the legend order; NA = "unlabeled".
##   annotation_name  Optional short name for the annotation (e.g. "cell type"),
##             shown as the static panel's title.
##   out_path  Output JSON path.
##   label     Short free-text string identifying this fit/run, shown in
##             the page footer (e.g. "compound benchmark, n=11700, seed=202").
##   digits    Rounding for point coordinates and heights (default 4).
##   prob_eps  Posterior probabilities below this are treated as zero and
##             dropped from the (sparse) exported representation (default
##             0.005) -- mixture posteriors are typically concentrated on a
##             handful of components per point, so this keeps the export
##             small regardless of how many components were fitted.
##
## Alignment: row i of `xy`, `z` and `annotation` must describe the same
## observation as row i of `X` -- this function can only check the lengths,
## and stops with an error if they disagree; verifying that the rows are in
## the same ORDER is the caller's job.
##
## The scatter panel appears only when BOTH a 2D layout (`xy`, or X[, 1:2]
## when `xy` is omitted) and a source of posterior probabilities (`z`, or `X` in
## the fit's native dimensionality) are available. A layout with no way to
## shade it is not enough on its own -- the app has nothing else to color
## the points by.

suppressPackageStartupMessages({
  library(PHM)
  library(jsonlite)
})

export_phm_json <- function(phm2, X = NULL, xy = NULL, z = NULL, annotation = NULL,
                             annotation_name = NULL, out_path = "phm_export.json",
                             label = NULL, digits = 4, prob_eps = 0.005) {

  stopifnot(!is.null(phm2), !is.null(phm2$paramsList))
  K <- length(phm2$paramsList)
  d <- length(phm2$paramsList[[1]]$mean)

  if (!is.null(X)) X <- as.matrix(X)
  if (!is.null(xy)) {
    xy <- as.matrix(xy)
    if (ncol(xy) != 2) stop("xy must have exactly 2 columns (it's the scatter panel's plot layout).")
  }

  ## No separate layout given: plot the first two columns of X. For a 2D fit
  ## that is the fitted space itself; for higher-dimensional X it is just
  ## X[, 1:2] (X itself, all columns, still drives the posteriors).
  if (is.null(xy) && !is.null(X) && ncol(X) >= 2) {
    if (ncol(X) > 2) {
      message("xy not supplied: plotting the first two columns of X (X[, 1:2]) in the scatter panel. ",
              "Pass `xy` to use a different 2D layout (e.g. UMAP/t-SNE/PCA).")
    }
    xy <- X[, 1:2, drop = FALSE]
  }

  has_layout <- !is.null(xy)
  n <- if (has_layout) nrow(xy) else if (!is.null(X)) nrow(X) else NA_integer_

  ## ---- alignment checks: sizes only -- row ORDER can't be verified here ----
  if (has_layout && !is.null(X) && nrow(X) != nrow(xy)) {
    stop(sprintf("xy has %d row(s) but X has %d -- they must describe the same observations, in the same order.",
                 nrow(xy), nrow(X)))
  }
  if (has_layout && !is.null(z) && nrow(z) != nrow(xy)) {
    stop(sprintf("z has %d row(s) but the layout has %d.", nrow(z), nrow(xy)))
  }
  has_annot <- !is.null(annotation)
  if (has_annot) {
    if (!has_layout) stop("annotation needs a 2D layout: pass `xy` (or a 2-column X) as well.")
    if (length(annotation) != n) {
      stop(sprintf("annotation has %d element(s) but the layout has %d row(s).", length(annotation), n))
    }
    if (!is.factor(annotation)) annotation <- factor(annotation)
    annotation <- droplevels(annotation)
  }

  ## ---- posterior probabilities (n x K), only needed if there's a layout to shade ----
  if (has_layout && is.null(z)) {
    if (is.null(X)) {
      ## With an annotation the static category panel still has something to
      ## show; only the posterior-shaded panel is dropped (no z to shade it).
      if (!has_annot) {
        warning("xy was supplied without X or z -- no posterior probabilities to shade the ",
                "scatter with, so the scatter panel will be omitted. Pass `z` directly, or ",
                "pass `X` (the fit's own original-space data) so it can be computed.")
        has_layout <- FALSE
      }
    } else if (ncol(X) != d) {
      stop(sprintf(
        paste("X has %d column(s) but the fit has %d dimension(s) -- pass the fit's FULL",
              "original-space data as X (used only to compute posteriors; never truncated",
              "or plotted), and use `xy` for whatever 2D layout you actually want plotted",
              "(e.g. a UMAP/t-SNE embedding). If you don't have the original data, compute",
              "posteriors yourself and pass them as `z` instead."),
        ncol(X), d
      ))
    } else {
      log_dens <- matrix(-Inf, n, K)
      for (k in seq_len(K)) {
        mu_k <- as.vector(phm2$paramsList[[k]]$mean)
        Sigma_k <- matrix(phm2$paramsList[[k]]$var, d, d)
        prob_k <- phm2$paramsList[[k]]$prob
        ch <- tryCatch(chol(Sigma_k), error = function(e) NULL)
        if (is.null(ch)) next
        logdet <- 2 * sum(log(diag(ch)))
        Xc <- sweep(X, 2, mu_k, "-")
        zz <- t(backsolve(ch, t(Xc), transpose = TRUE))
        maha <- rowSums(zz^2)
        log_dens[, k] <- log(prob_k) - 0.5 * d * log(2 * pi) - 0.5 * logdet - 0.5 * maha
      }
      m <- apply(log_dens, 1, max)
      unnorm <- exp(log_dens - m)
      z <- unnorm / rowSums(unnorm)
    }
  }

  has_scatter <- has_layout && !is.null(z)
  has_points <- has_scatter || (has_layout && has_annot)

  ## ---- authoritative leaf order + dendrogram line-segment geometry ----
  vis <- PHM:::constructVisData(phm2, K = K, scaleHeights = "d*")
  xlab_order <- as.integer(vis$xlab)          # position (1..K) -> component label
  pos_of <- setNames(seq_along(xlab_order), xlab_order)

  ## ---- reconstruct the exact merge tree via recoverGroupsv2 ----
  partitions <- list()
  for (k in K:1) partitions[[as.character(k)]] <- recoverGroupsv2(phm2, k)

  find_merge <- function(part_hi, part_lo) {
    for (lo in part_lo) {
      matches <- Filter(function(hi) all(hi %in% lo), part_hi)
      if (length(matches) == 2 && identical(sort(unlist(matches)), sort(lo))) return(matches)
    }
    NULL
  }

  merges <- list()
  for (k in K:2) {
    hi <- partitions[[as.character(k)]]
    lo <- partitions[[as.character(k - 1)]]
    m <- find_merge(hi, lo)
    merges[[length(merges) + 1]] <- list(left = m[[1]], right = m[[2]])
  }
  stopifnot(length(merges) == K - 1)
  lab <- vis$labels   # one row per merge, same chronological order as `merges`

  ## Heatmap tile colors are NOT precomputed/exported here. plotPHMv2Heatmap()
  ## colors cell (i,j) purely as a function of the d*-scaled height of the
  ## merge that first joins i and j -- exactly the `height` every internal
  ## node already carries below -- so the app derives tile colors from that
  ## directly, live, for whatever leaf count is currently on screen. That
  ## avoids ever materializing a K x K color matrix (it would dominate the
  ## export size for any fit with more than a few dozen components) and
  ## means the heatmap recolors correctly at any display threshold without
  ## needing a second export.

  ## ---- sparse per-point posterior probabilities (only components with
  ## prob > prob_eps are kept; this is what lets the export stay small even
  ## for a pooled/ensemble fit with hundreds of components) ----
  sparse_z <- NULL
  if (has_scatter) {
    sparse_z <- vector("list", n)
    for (i in seq_len(n)) {
      row <- z[i, ]
      idx <- which(row > prob_eps)
      ord <- idx[order(-row[idx])]
      sparse_z[[i]] <- lapply(ord, function(k) list(unname(k), round(row[k], 3)))
    }
  }

  ## ---- build 2K-1 "nodes": K leaves + (K-1) internal merges, in the SAME
  ## full-fidelity form the interactive app already knows how to draw; the
  ## app collapses this down to fewer leaves for display, at render time ----
  nodes <- vector("list", K + (K - 1))

  for (pos in seq_len(K)) {
    comp <- xlab_order[pos]
    nodes[[pos]] <- list(
      type = "leaf", id = paste0("leaf", comp), members = comp,
      pos_min = pos, pos_max = pos, height = 0,
      x = pos, xend = pos, y = 0, yend = 0
    )
  }
  for (i in seq_len(K - 1)) {
    mrg <- merges[[i]]
    all_members <- sort(c(mrg$left, mrg$right))
    positions <- pos_of[as.character(all_members)]
    nodes[[K + i]] <- list(
      type = "internal", id = paste0("node", i),
      members = all_members,
      left_members = I(mrg$left), right_members = I(mrg$right),
      pos_min = min(positions), pos_max = max(positions),
      height = lab$y[i],
      x = lab$x[i], xend = lab$xend[i], y = lab$y[i], yend = lab$yend[i]
    )
  }

  segs <- as.data.frame(vis$df)[, c("x", "y", "xend", "yend")]

  export <- list(
    meta = list(
      K = K, mergeCriterion = phm2$mergeCriterion,
      n = if (has_points) n else NULL,
      label = label,
      generatedAt = format(Sys.time(), "%Y-%m-%d %H:%M %Z")
    ),
    xlabOrder = xlab_order,
    dendroSegments = segs,
    nodes = lapply(nodes, function(nd) { nd$members <- I(nd$members); nd })
  )
  if (has_points) {
    export$points <- list(x1 = round(xy[, 1], digits), x2 = round(xy[, 2], digits))
  }
  if (has_scatter) export$z <- sparse_z
  if (has_points && has_annot) {
    export$annotation <- list(
      name = if (is.null(annotation_name)) "annotation" else annotation_name,
      levels = I(levels(annotation)),
      codes = I(as.integer(annotation))      # 1-based index into levels; NA -> null
    )
  }

  write_json(export, out_path, auto_unbox = TRUE, digits = digits, na = "null")
  cat("wrote", out_path, ",", file.size(out_path), "bytes  (K =", K,
      if (has_points) paste0(", n = ", n) else ", no scatter data",
      if (has_points && has_annot) paste0(", annotation: ", nlevels(annotation), " categories") else "", ")\n")
  invisible(export)
}
