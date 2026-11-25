# PHM workflows

Worked recipes for the situations that come up most, using `PHMv2()` (the default in this skill). Section 11 covers the original `PHM()`. Every snippet assumes

```r
library(PHM); library(mclust)    # mclust only because some recipes call Mclust() directly
set.seed(1)
```

and a numeric matrix `X` (N × D). Recipes marked "tested" are exercised by `tests/smoke_test.R` (in the folder that holds this skill). The ensemble pipeline has its own file, `ensemble.md`, and PHM Explorer is in `explorer.md`.

## Contents

1. [GMM to hierarchy (tested)](#1-gmm-to-hierarchy)
2. [Any clustering method to hierarchy (tested)](#2-any-clustering-method-to-hierarchy)
3. [Comparing several clusterers on one dataset (tested)](#3-comparing-several-clusterers-on-one-dataset)
4. [Soft memberships as weights (tested)](#4-soft-memberships-as-weights)
5. [Choosing a resolution and extracting labels (tested)](#5-choosing-a-resolution-and-extracting-labels)
6. [Merge criteria, reusing ΔP_mc, zooming (tested)](#6-merge-criteria-reusing-δp_mc-zooming)
7. [Ensemble: consensus over subsamples (tested)](#7-ensemble-consensus-over-subsamples)
7b. [Localized density estimation (tested)](#7b-localized-density-estimation)
8. [High-dimensional and large data](#8-high-dimensional-and-large-data)
9. [Single-cell embeddings (Seurat)](#9-single-cell-embeddings-seurat)
10. [Where does P_mc come from? (2-D diagnostic)](#10-where-does-p_mc-come-from-2-d-diagnostic)
11. [The original `PHM()` workflow (tested)](#11-the-original-phm-workflow)

## 1. GMM to hierarchy

Overfit a GMM on purpose (many components), then let PHM show which components belong together. This is the setting of the nested-diamonds example, where a 64-component fit reveals three layers of structure.

```r
data("moons", package = "PHM")
X <- moons
mcl    <- Mclust(X, G = 1:15)
params <- constructPmcParamsMclust(mcl)
phm2   <- PHMv2(paramsList = params, batchSize = 1e5)

dend <- plotPHMv2Dendrogram(phm2, displayAxis = "label")
heat <- plotPHMv2Heatmap(phm2, displayAxis = "label")
ggpubr::ggarrange(dend, heat, ncol = 1)
```

For many components (say more than 12) use `displayAxis = "none"` in both plots. In 2-D, `monteCarlo = FALSE` (cubature) is affordable and exact up to tolerance, but Monte Carlo is the safer default.

## 2. Any clustering method to hierarchy

PHM only needs labels. The weighted constructor turns them into per-cluster densities.

```r
hcl    <- hclust(dist(X), method = "ward.D2")
labels <- cutree(hcl, k = 10)                     # or kmeans(X, 10, nstart = 5)$cluster, Leiden, DensityCut, ...

params <- constructPmcParamsWeightedPartition(labels, X)          # slow for N in the tens of thousands
phm2   <- PHMv2(paramsList = params, batchSize = 1e5)
```

Notes:

- Start from a deliberately fine partition (more clusters than you believe exist); PHM merges upward, so it cannot split an under-segmented cluster.
- Use integer or character labels. Factors also work, but `as.integer(factor(x))` is the safe habit.
- Extra arguments go to `Mclust`: `constructPmcParamsWeightedPartition(labels, X, G = 1:3, modelNames = c("VVI", "VVV"))` limits the per-cluster mixture size and covariance structure, which also speeds things up.
- To compare against the naive estimator, swap in `constructPmcParamsPartition(labels, X)`. It typically reports lower P_mc between neighbouring clusters (see the weighted-density tutorial).

## 3. Comparing several clusterers on one dataset

Run PHMv2 on each solution with the same `scaling`, and put the heatmaps on a common color scale.

```r
cands <- list(
  gmm    = mcl$classification,
  hclust = cutree(hcl, 8),
  kmeans = kmeans(X, 8, nstart = 5)$cluster
)
phms <- lapply(cands, function(lab) {
  PHMv2(paramsList = constructPmcParamsWeightedPartition(lab, X), batchSize = 1e5)
})

rng <- range(unlist(lapply(phms, function(p) p$mergeValues[-1])))   # raw merge-criterion scale
plots <- lapply(names(phms), function(nm) {
  ggpubr::ggarrange(
    plotPHMv2Dendrogram(phms[[nm]], displayAxis = "label") + ggplot2::ggtitle(nm),
    plotPHMv2Heatmap(phms[[nm]], displayAxis = "label", fillLimits = rng),
    ncol = 1)
})
ggpubr::ggarrange(plotlist = plots, nrow = 1)
```

In the v2 heatmap `fillLimits` is on the raw `mergeValues` scale (index 1 of `mergeValues` is a placeholder zero, hence `[-1]`). What to look for: clusters that merge early, at large values, in every method are the same "cluster" seen at different resolutions; a split that survives to the top of the tree in one method but not another is method-specific. The tutorials choose K for k-means and hierarchical clustering by silhouette score, and let the GMM choose G by BIC.

## 4. Soft memberships as weights

Any N × K matrix of non-negative memberships can replace the default distance weights. Its **columns must follow `unique(labels)` (order of first appearance)**, not sorted label order.

```r
lab <- mcl$classification
W   <- mcl$z[, unique(lab), drop = FALSE]        # columns reordered to match unique(lab)
params <- constructPmcParamsWeightedPartition(lab, X, weights = W, G = 1:3)
phm2   <- PHMv2(paramsList = params, batchSize = 1e5)
```

The weighted-density tutorial does the same with fuzzy c-means (`ppclust::fcm(X, centers = K)`, using `fcl$u` as weights and `fcl$cluster` as labels). Reorder `fcl$u`'s columns with `unique(fcl$cluster)` whenever the first observation is not in cluster 1. In the tutorial's 1-D simulation, c-means weights recovered the true overlap best, the unweighted estimator under-estimated it and the default distance weights over-estimated it.

## 5. Choosing a resolution and extracting labels

```r
pmc_by_k <- cumsum(phm2$mergeDeltaPmc)           # remaining P_mc at each K (index = K)
path <- data.frame(K = seq_along(pmc_by_k), pmc = pmc_by_k, drop = phm2$mergeDeltaPmc)
path                                             # `drop` = P_mc removed by the merge from K to K-1; a big value = overlapping clusters merged

K_star <- max(which(pmc_by_k < 0.01))            # finest K with remaining P_mc < 0.01

union_labels <- function(phm2, params, orig_labels, k) {
  grp <- recoverGroupsv2(phm2, k)                # positions in `params` for each merged cluster
  cls <- sapply(params, `[[`, "class")           # position -> original label
  out <- integer(length(orig_labels))
  for (g in seq_along(grp)) out[orig_labels %in% cls[grp[[g]]]] <- g
  out
}
lab <- union_labels(phm2, params, labels, K_star)     # exact unions of the input clusters

table(truth, lab)                                # against known classes, if available
adjustedRandIndex(truth, lab)                    # mclust

# Soft alternative: posterior probabilities of the merged clusters
post     <- recoverPosteriorv2(phm2, K_star, computePosteriorProbMatrix(params, X))
lab_soft <- apply(post, 1, which.max)
```

`pmc_by_k` needs no extra function because remaining P_mc at K is the sum of the ΔP_mc matrix at K, and each merge removes `mergeDeltaPmc[k]` from it. The threshold is a tuning choice: report it, and show the dendrogram so the reader sees the alternatives. For a GMM input use `mcl$classification` as `labels`. `union_labels()` gives hard labels that are exact unions of the input clusters; `lab_soft` follows the fitted densities and can put a few boundary points elsewhere, and its integers use the same cluster order as `recoverGroupsv2()`.

## 6. Merge criteria, reusing ΔP_mc, zooming

The expensive part is the ΔP_mc matrix, and it does not depend on the merge criterion. Compute it once, then try criteria:

```r
phm_alpha <- PHMv2(paramsList = params, scaling = "alpha", batchSize = 1e5)
D <- phm_alpha$deltaPmc
phm_raw   <- PHMv2(paramsList = params, deltaPmc = D, scaling = "unscaled")
phm_avg   <- PHMv2(paramsList = params, deltaPmc = D, scaling = "average")

plotPHMv2Dendrogram(phm_alpha, scaleHeights = "d*")        # d* is for alpha runs
plotPHMv2Dendrogram(phm_raw,   scaleHeights = "log10")     # non-alpha: log10, unscaled or pmcdist
plotPHMv2Heatmap(phm_raw,      fillScale = "pmcdist")
```

Passing both `paramsList` and `deltaPmc` skips the Monte Carlo step (the two are only checked for matching size). `"unscaled"` reproduces the original `PHM()` merges; `"alpha"` (default) divides by the product of the subtrees' total mixture weights; `"average"` divides by the product of the subtree sizes in original clusters. To zoom in:

```r
grp <- recoverGroupsv2(phm_alpha, k = 4)
plotPHMv2Heatmap(phm_alpha, subsetClusters = grp[[1]], displayAxis = "label")   # one branch of the K = 4 solution
plotPHMv2Dendrogram(phm_alpha, initK = 6)                                       # start the picture at K = 6
```

`convertToPHMv2(phm)` turns an existing `PHM()` result into a v2-style object for these plots.

## 7. Ensemble: consensus over subsamples

Fit many GMMs on random subsamples, pool all their components, and let PHM merge the redundant ones. Structure that appears in most replicates shows up as components merging early, at large ΔP_mc. This is the package's main route from raw data.

```r
fit <- fitPHMEnsemble(X, M = 10, rho = 0.5, G = 1:25, numCores = 2)   # returns a PHMv2 object
saveRDS(fit, "phm_fit.rds")

pmc_by_k <- cumsum(fit$mergeDeltaPmc)
plotPHMv2Dendrogram(fit, initK = 20, displayAxis = "none")

K <- 2
post <- recoverPosteriorv2(fit, K, computePosteriorProbMatrix(fit$paramsList, X))
lab  <- apply(post, 1, which.max)
```

Pooled components are labelled `"<component>_<replicate>"` (e.g. `"3_7"`). `globalDensityEstimation(X, M = 10, rho = 0.5)` returns the same pooled `paramsList` without running PHMv2 (and can fit replicates in parallel with `numCores` and cache them with `saveDir`), so `PHMv2(paramsList = globalDensityEstimation(...))` with a different `scaling` or `mcSamples` is the way to vary the reduce step. For runs that take hours, or to fit replicates on several cores or a cluster with per-replicate reruns, use the command-line pipeline in `ensemble.md`.

Older code uses `constructPmcParamsSubAggMclust(X, replicates = 20, subsampSize = 1000, G = 1:10, saveDir = "subagg_cache")`, which pools plain `Mclust` fits (labels `"<component>_<replicate>"`) without G extension or per-replicate seeds. It still works and reproduces the published analyses.

For a clustering method instead of a GMM, use `localizedDensityEstimation()` with `M` and `rho` (section 7b), which reruns the clusterer on each subsample. `constructPmcParamsSubAggPartition()` is disabled in the current package (it stops with "This method does not work").

## 7b. Localized density estimation

When the user has a clustering *method* rather than fixed labels, `localizedDensityEstimation()` runs it (on the whole data or on each of M subsamples), weights points by distance, and fits a weighted GMM per cluster.

```r
km <- function(X, K, seed = 1, ...) {            # receives X and seed; must return labels 1..K
  set.seed(seed); kmeans(X, K, nstart = 5)$cluster
}
params <- localizedDensityEstimation(X, partitionFunc = km, partitionFuncParams = list(K = 8),
                                     weightFunc = distanceWeightFunc, preserveClusters = TRUE,
                                     G = 1:3)
fit <- PHMv2(paramsList = params, batchSize = 1e5)
```

`preserveClusters = TRUE` keeps one mixture per cluster (leaves = clusters); `FALSE` splits them into single Gaussians (leaves = components). `M` and `rho` make it an ensemble of partitions. Labels must be consecutive integers starting at 1 (see `api.md`). `weightFunc = naiveWeightFunc` gives the unweighted version.

## 8. High-dimensional and large data

- Reduce first. The tutorials run PHM on the top principal components (HGDP: the elbow of the PC variances; PBMC: 9 PCs from an elbow computed on 600 PCs). PHM's Gaussian densities and Monte Carlo integral degrade quickly with dimension.

```r
pc  <- prcomp(expr, center = TRUE)                       # rows = observations
npc <- PCAtools::findElbowPoint(pc$sdev[1:200]^2)        # or choose by variance explained / domain knowledge
X   <- pc$x[, 1:npc]
```

- Cost drivers: the number of starting clusters K (the ΔP_mc matrix is K × K), N (Monte Carlo posteriors, weighted-partition distances), and D. Levers, in order of usefulness: fewer PCs, restricted `modelNames` for the GMMs, `batchSize = 1e5`, `numCores` (each worker holds its own batch, so more cores also means more memory), `localizedDensityEstimation(..., rho = )` on subsamples (or supplying `weights`) for the density step, lower `mcSamples` for a first look (raise it for final numbers).
- Save expensive objects: `saveRDS(params, "params.rds")`, `saveRDS(phm2, "phm2.rds")`. The tutorials recommend this because density estimation and ΔP_mc dominate the runtime.

## 9. Single-cell embeddings (Seurat)

The PBMC tutorial runs PHM on the leading PCs of a Seurat object with several clusterers (GMM, DensityCut, Leiden through monocle3, k-means, hierarchical) and compares the hierarchies.

```r
library(Seurat)
embed  <- Embeddings(obj, "pca")[, 1:npc]                    # matrix, cells x PCs
labels <- as.integer(Idents(obj))                            # e.g. Leiden clusters from FindClusters(..., algorithm = 4)

params <- constructPmcParamsWeightedPartition(labels, embed)
phm2   <- PHMv2(paramsList = params, batchSize = 1e5)
plotPHMv2Dendrogram(phm2, displayAxis = "label")

lab4 <- union_labels(phm2, params, labels, k = 4)            # helper from section 5
table(obj$seurat_annotations, lab4)                          # relate to annotations
```

## 10. Where does P_mc come from? (2-D diagnostic)

For 2-D data, `plotPmc2D` shades the regions that contribute most to P_mc and overlays each cluster's density contours, which explains why two clusters have a large ΔP_mc. It works on `params`, so it is independent of `PHM()` vs `PHMv2()`.

```r
plotPmc2D(params, X, labels, colors = RColorBrewer::brewer.pal(length(params), "Paired"),
          PmcColor = "#222222",                              # the default "#222" is rejected by ggplot2
          numPmcPatches = 150, densityLevels = c(1e-2, 1e-1))
```

`suppressPmc = TRUE`, `suppressDensity = TRUE` and `suppressObservations = TRUE` switch layers off.

## 11. The original `PHM()` workflow

Use it to reproduce the paper or tutorials, or when you want per-K labels and posteriors stored on the object and the `plotPmcMatrix()` / `plotPHMDistruct()` figures.

```r
# from a GMM (each component is a cluster)
phm <- PHM(mclustObj = mcl, batchSize = 1e5)

# from a partition
params <- constructPmcParamsWeightedPartition(labels, X)
phm    <- PHM(paramsList = params, data = X, partition = labels, batchSize = 1e5)
# shortcut: PHM(partition = labels, data = X)   (weighted, VVI, up to 10 components per cluster)

dend <- plotPHMDendrogram(phm, scaleHeights = "log10", displayAxis = "label")
heat <- plotPHMMatrix(phm, fillScale = "pmcdist", displayAxis = "label")
ggpubr::ggarrange(dend, heat, ncol = 1)
plotPmcMatrix(phm, K = 5)                                    # ΔP_mc matrix at K = 5, with numbers
plotPHMDistruct(phm, K = 5, labels = factor(truth), colors = RColorBrewer::brewer.pal(5, "Set1"))
```

`phm[[K]]` is the K-cluster solution (`labels`, `posterior_matrix`, `pmc`, `pmc_change`, `merge_components`, `params`, `pmc_matrix`; see `api.md`). Extracting a solution:

```r
path <- data.frame(K = seq_along(phm), pmc = sapply(phm, `[[`, "pmc"), drop = sapply(phm, `[[`, "pmc_change"))
sol    <- thresholdPHM(phm, threshold = 0.01)               # finest level with remaining P_mc < 0.01
labels <- sol$labels                                         # merged clusters look like "2|5" (character once any merge happened)
K      <- sol$clusters
phm[[4]]$posterior_matrix                                    # do not use recoverPosterior(): it crashes on default output
```

- Comparing several `PHM()` runs: give `plotPHMMatrix` a common `fillLimits` on the raw `pmc_change` scale, `range(na.omit(unlist(lapply(phms, function(p) sapply(p, "[[", "pmc_change")))))`.
- Memory: `computePosterior = FALSE`, `storeParams = FALSE`, `storeDeltaPmc = FALSE`; fill posteriors later with `addPosteriorMatrix(phm, X)` and rebuild a level with `recoverPHMParams()` / `recoverDeltaPmcMatrix()`.
- Coarser-level `labels` come from the merged posterior argmax, so they are not guaranteed to be exact unions of the input clusters (in a test with 8 clusters, about 15 of 1,000 points were relabeled at K = 7). For exact unions from a `PHM()` result, convert and use the v2 helper from section 5: `union_labels(convertToPHMv2(phm), params, labels, k)`.
- `PHM(deltaPmc = )` is ignored; to reuse a precomputed matrix use `PHMv2(deltaPmc = D, scaling = "unscaled")`, which gives the same merges as `PHM()`.
