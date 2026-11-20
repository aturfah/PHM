#' @keywords internal
"_PACKAGE"

## Several mclust functions rebuild a call with a bare function name and
## evaluate it in the caller's frame: Mclust() calls `mclustBIC`, hc() calls
## `hcVVV`/`hcEII`/..., and so on. When the caller is PHM code, that lookup
## goes through PHM's namespace, so mclust is imported whole. Without this,
## those functions fail with "could not find function" unless the user has
## attached mclust with library(mclust).
#' @import mclust
#' @importFrom stats mahalanobis rmultinom
NULL

## Column names used inside ggplot2::aes() and dplyr verbs in the plotting
## code. Declared here so R CMD check does not report them as undefined
## global variables.
utils::globalVariables(c(
  "ID", "ID.x", "X", "X1", "X2", "Xcenter", "Xraw", "Xwidth",
  "Y", "Ycenter", "Yraw", "Ywidth", "Z", "Z.mod",
  "count", "g", "gprob", "gprob.x", "lab", "label_column",
  "linetype", "linetype.x", "maxX", "max_column", "minX", "name", "part",
  "pmc.x", "pmc_change.x", "pmc_pct", "pmc_pct.x", "sort_val", "value",
  "x", "x.x", "x.y", "xend", "xposn", "y", "yend"
))

## The spline that maps log10(Pmc) to the "pmcdist" height scale used by the
## dendrogram and heatmap plots, stored in inst/extdata.
pmcScaleFunction <- function() {
  env <- new.env()
  load(system.file("extdata", "pmc_scale_function.RData", package = "PHM"), envir = env)
  env$inv_log10
}
