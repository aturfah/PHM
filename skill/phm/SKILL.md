---
name: phm
description: Use the PHM R package (github.com/aturfah/PHM; Turfah & Wen 2025, "Multiresolution Clustering of Genomic Data") and its companion PHM Explorer for multiresolution and consensus clustering. PHM estimates a density (a Gaussian mixture, from raw data or from any existing clustering), measures cluster overlap with the misclassification probability P_mc / ΔP_mc (the "distinguishability criterion"), and merges clusters into a hierarchy shown as a dendrogram and heatmap. Covers density estimation (global, localized, weighted partitions, subsampled ensembles), PHMv2 and the original PHM(), cutting the tree, plots, the PHM ensemble pipeline (fitPHMEnsemble, fitDensityReplicate, combineDensityReplicates, the serial / single-machine / Slurm scripts in inst/ensemble), and exporting a fit to the interactive PHM Explorer page. Use this skill whenever the user mentions PHM, PHMv2, PHM ensemble, PHM Explorer, P_mc, Pmc, ΔPmc, deltaPmc, misclassification probability, distinguishability, merging GMM or Mclust components into a hierarchy, multiresolution or consensus clustering, reconciling k-means/hclust/Leiden/GMM solutions, "which clusters are really distinct", fitPHMEnsemble, globalDensityEstimation, constructPmcParams*, plotPHMv2Dendrogram, or view_in_explorer, even if they do not name the package.
---

# PHM: multiresolution clustering with the distinguishability criterion

PHM is an R package. It takes a density made of many Gaussian pieces, decides how strongly each pair of pieces overlaps, and merges them step by step into a hierarchy from K clusters down to 1. The hierarchy shows which clusters are clearly distinct and which are refinements of coarser structure, so different clusterers (or different K) can be reconciled instead of arguing about one "best" K.

There are two ways in:

1. **From raw data** (no clustering yet): fit a deliberately rich Gaussian mixture, usually as an ensemble over random subsamples (`fitPHMEnsemble()`), and let PHM merge its components. This is the package's main route for real data and the one the ensemble pipeline scales up.
2. **From an existing clustering** (labels from k-means, hclust, Leiden, DensityCut, ...): turn the labels into per-cluster densities (`constructPmcParamsWeightedPartition()`), then run `PHMv2()`.

Both end in a `PHMv2` object, which the `recover*v2()` helpers, the `plotPHMv2*()` plots and PHM Explorer all read. Default to `PHMv2()` and its helpers. The original `PHM()` is covered in "When to use the original `PHM()`" below.

## The idea in one paragraph

Each cluster has a fitted density (a Gaussian mixture), which gives a posterior cluster probability π_j(x) at every point. The **misclassification probability** is P_mc = E[ Σ_j π_j(x)(1 − π_j(x)) ], the overlap between clusters under the fitted overall density. It splits over pairs: P_mc = Σ_{i<j} 2·ΔP_mc(i,j) with ΔP_mc(i,j) = E[π_i(x)·π_j(x)], estimated by Monte Carlo. Merging clusters i and j removes 2·ΔP_mc(i,j) from P_mc, and the merged cluster's ΔP_mc row is the sum of the two old rows. PHM repeatedly merges the best pair according to a merge criterion built from ΔP_mc. A tiny ΔP_mc at a merge means the two sides were well separated (they only join near the top of the tree); a large one means they were nearly indistinguishable.

## Setup

- Install: `remotes::install_github("aturfah/PHM")` (or `devtools::`). The bundled `scripts/setup_phm.R` checks and installs everything: `Rscript scripts/setup_phm.R --check-only` reports what is missing and installs nothing; without the flag it installs the CRAN dependencies and PHM, then runs a small end-to-end check. `--explorer` adds jsonlite and checks for python3 (PHM Explorer), `--optional` adds ggpubr and ppclust. Run `--check-only` first, and get the user's go-ahead before installing packages on their machine. Source builds (Homebrew R on macOS, Linux) need C and Fortran compilers and take several minutes.
- Load with `library(PHM)`. The package imports mclust and ggplot2 itself, so its own functions work without attaching anything else. Attach `library(mclust)` when the user's own code calls `Mclust()`, `hc()` or `adjustedRandIndex()` unqualified, because `Mclust()` looks up `mclustBIC` by bare name in the caller's environment.
- Data must be a numeric **matrix** (N × D), not a data.frame (`as.matrix()` first), with finite values.
- Gaussian densities and the Monte Carlo integral behave best in low dimension. For high-dimensional data (expression, genotypes) run PHM on a PCA embedding (the package's tutorials use 5 to 30 PCs).
- ΔP_mc is Monte Carlo. Call `set.seed()` before `PHMv2()` when results must be reproducible. Ensemble replicates get their own seeds from `phmSeeds()`.
- Parallel code uses forked `parallel::mclapply` (`numCores`), so it runs serially on Windows. If a forked worker dies (most often killed for running out of memory) the package stops with "k of n parallel task(s) failed ... try a smaller numCores or batchSize" instead of returning partial results.

## Pick the entry point

| What the user has | Build the mixture | Then |
|---|---|---|
| Raw data, no clustering (the default for discovery) | `fit <- fitPHMEnsemble(X, M = 10, rho = 0.5)` | done: `fit` is already a `PHMv2` object |
| Raw data, large N or long runs, several cores or a Slurm cluster | the command-line pipeline in `system.file("ensemble", package = "PHM")` (see `references/ensemble.md`) | `readRDS(outPath)` is a `PHMv2` object |
| Raw data, a single GMM | `params <- globalDensityEstimation(X)` or `constructPmcParamsMclust(Mclust(X, G = 1:30))` | `PHMv2(paramsList = params)` |
| Cluster labels from k-means / hclust / Leiden / DensityCut / ... | `constructPmcParamsWeightedPartition(labels, X)` (recommended) or `constructPmcParamsPartition(labels, X)` (naive, fast) | same |
| Soft memberships (fuzzy c-means, `mcl$z`) | `constructPmcParamsWeightedPartition(labels, X, weights = W)` | same |
| A clustering *function* to apply inside the pipeline, optionally over subsamples | `localizedDensityEstimation(X, partitionFunc = f, weightFunc = distanceWeightFunc, ...)` | same |
| Their own Gaussian mixtures | hand-build the list described in `references/api.md` | same |
| Only a precomputed ΔP_mc matrix | none | `PHMv2(deltaPmc = D, scaling = "unscaled")` (alpha scaling needs `paramsList`) |

Weighted vs naive partition densities: the naive estimator fits each cluster only on its own points, so it under-estimates tails and under-states overlap. The weighted one gives every point a soft weight for every cluster (from its distance to the cluster), which accounts for uncertainty in the partition. Prefer weighted unless N is large (see Scale).

Ensemble vs single GMM: one Mclust fit picks a single G and starting point; pooling M fits on random subsamples (`rho` of the rows each, without replacement) smooths over that choice, and PHM then merges the redundant components across replicates. Components that describe the same structure in many replicates merge early, at large ΔP_mc. Labels of pooled components look like `"3_7"` (component 3 of replicate 7).

## Core recipe

```r
library(PHM)

## route 1: raw data -> ensemble (pooled mixture + PHMv2 in one call)
fit <- fitPHMEnsemble(X, M = 10, rho = 0.5, G = 1:25, numCores = 2)
params <- fit$paramsList

## route 2: labels from any clusterer -> per-cluster densities -> PHMv2
params <- constructPmcParamsWeightedPartition(labels, X)   # G, modelNames pass through to Mclust
set.seed(1)
fit <- PHMv2(paramsList = params, mcSamples = 1e6, batchSize = 1e5, numCores = 1)
```

A `PHMv2` result is a small list: `deltaPmc` (the K × K ΔP_mc matrix), `paramsList`, `mergeComps` (`mergeComps[[k]]` is `c(i, j)`, the positions merged to go from k clusters to k − 1), `mergeValues` (merge criterion at each merge), `mergeDeltaPmc` (unscaled 2·ΔP_mc of each merge) and `mergeCriterion`. Index 1 of the merge vectors is a placeholder. Everything else is derived with the helpers below. `fitPHMEnsemble()` and `combineDensityReplicates()` call `PHMv2()` with its defaults (`scaling = "alpha"`) and `mcSamples = 2e5`.

### Read out a resolution

```r
pmc_by_k <- cumsum(fit$mergeDeltaPmc)           # remaining P_mc at each K (index = K); 0 at K = 1
K_star   <- max(which(pmc_by_k < 0.01))         # finest K whose remaining P_mc is below 0.01
                                                # (thresholdPHM() does this for PHM() results)
groups <- recoverGroupsv2(fit, k = K_star)      # list: positions in paramsList making up each cluster
recoverDeltaPmcv2(fit, K_star)                  # K_star x K_star ΔP_mc matrix
recoverParamsv2(fit, K_star)                    # merged cluster densities

# Soft labels for every observation: posterior of the merged clusters, then argmax
post <- recoverPosteriorv2(fit, K_star, computePosteriorProbMatrix(fit$paramsList, X))   # N x K_star
lab  <- apply(post, 1, which.max)

# Labels that are exact unions of the INPUT clusters (partition route only)
union_labels <- function(fit, params, orig_labels, k) {
  grp <- recoverGroupsv2(fit, k)
  cls <- sapply(params, `[[`, "class")          # position in params -> original label
  out <- integer(length(orig_labels))
  for (g in seq_along(grp)) out[orig_labels %in% cls[grp[[g]]]] <- g
  out
}
```

The 0.01 threshold is an example, not a rule. Show the user the P_mc trajectory (`pmc_by_k`) and the dendrogram and let the scientific question pick the resolution; the merge heights also show natural gaps. Positions in `paramsList` follow `unique(partition)` (order of first appearance) for the partition constructors, component order for Mclust, and replicate-then-component order for ensembles. Group order in `groups`, columns of `post` and integers in `lab` all use the same order. For ensembles and other pooled mixtures there are no input labels, so use the posterior route.

`computePosteriorProbMatrix()` builds an N × K matrix. With thousands of pooled components and a large N, compute it in row blocks, or on a subsample, and keep only the merged columns.

### Merge criterion (`scaling`)

- `"alpha"` (default): ranks pairs by ΔP_mc divided by the product of the two subtrees' total mixture weights, i.e. overlap relative to cluster mass.
- `"unscaled"`: ranks by raw ΔP_mc. This is the rule of the original `PHM()` and gives the same merge order and the same `pmc_by_k`.
- `"average"`: ΔP_mc divided by the product of the subtree sizes counted in original clusters.

The ΔP_mc matrix does not depend on the criterion, so `PHMv2(paramsList = params, deltaPmc = fit$deltaPmc, scaling = "unscaled")` re-merges without new Monte Carlo. `pmc_by_k` is the sum of the current ΔP_mc matrix whichever criterion chose the merges, so it compares across criteria only where the clusters at that K agree.

## Visualize

```r
dend <- plotPHMv2Dendrogram(fit, scaleHeights = "d*", displayAxis = "label")
heat <- plotPHMv2Heatmap(fit, fillScale = "d*", displayAxis = "label")
ggpubr::ggarrange(dend, heat, ncol = 1)          # they share the leaf order

plotPHMv2Heatmap(fit, subsetClusters = groups[[1]])   # zoom into one branch
plotPHMv2Dendrogram(fit, initK = 30)                  # start the picture at K = 30 (essential for big ensembles)
plotPmc2D(params, X, labels, PmcColor = "#222222")    # 2-D data only: where P_mc comes from
```

- `scaleHeights` / `fillScale`: `"d*"` (default, alpha runs only), `"log10"`, `"unscaled"` (dendrogram only) and `"pmcdist"` (not for alpha runs). For a non-alpha run use `"log10"` or `"pmcdist"`.
- Default axis boxes use the Paired palette (K ≤ 12). For more leaves pass `colors =` or `displayAxis = "none"`, and for pooled ensembles with hundreds of components use `initK` so the plot starts at a readable resolution.
- `plotPmc2D()`'s default `PmcColor = "#222"` is rejected by current ggplot2 ("Malformed colour string"); always pass a six-digit color such as `"#222222"`.
- To compare runs side by side, give the heatmaps common `fillLimits` on the raw `mergeValues` scale (`range(unlist(lapply(runs, function(p) p$mergeValues[-1])))`), with the same `scaling` in every run.

## Interactive viewing: PHM Explorer

PHM Explorer is a self-contained HTML page (dendrogram, ΔP_mc heatmap, and a scatter shaded by posterior probability) for browsing a `PHMv2` fit, collapsing it live to any number of leaves, and hovering branches to see which points they claim. The skill bundles it in `scripts/explorer/`:

```sh
Rscript scripts/explorer/r/view_in_explorer.R phm_fit.rds X.rds --scatter=scatter.rds
```

writes `phm_fit_export.json` and a standalone `phm_fit_explorer.html` next to the fit and opens it. `X.rds` adds the scatter panel, shaded by posterior probability and laid out on `X[, 1:2]`; `--scatter=` gives better 2-D coordinates (UMAP, t-SNE) and an optional `annotation` column, which matters when the data have more than two dimensions. From R, `source("scripts/explorer/r/export_phm_json.R")` and call `export_phm_json(fit, X = X, xy = embedding, annotation = cell_type, out_path = "fit.json")`. The exporter expects an alpha-scaled fit whose `paramsList` entries are single Gaussians (ensembles, `globalDensityEstimation()`, Mclust). For partition fits with multi-component clusters, pass `z = computePosteriorProbMatrix(params, X)` instead of `X`. Details, the scatter-data format and troubleshooting are in `references/explorer.md`. Requires jsonlite and python3.

## Scale and speed

- The cost drivers are the number of starting components K (the ΔP_mc matrix is K × K; pooled ensembles reach thousands), N (posteriors, weighted-partition distances) and D.
- ΔP_mc uses `mcSamples` draws (`PHMv2()` default 1e6, ensemble default 2e5), processed `batchSize` at a time. Lower `batchSize` caps memory; `numCores` multiplies it, since each forked worker holds its own batch. A failure such as "1 of 4 parallel task(s) failed ... no result returned" means a worker was killed: lower `numCores` or `batchSize`.
- `monteCarlo = FALSE` (cubature) is only practical in very low dimension.
- `constructPmcParamsWeightedPartition` is roughly O(N²·D). For tens of thousands of points supply your own `weights`, use the naive constructor, fit on a subsample, or run `localizedDensityEstimation()` with `rho < 1` (the older `constructPmcParamsSubAggPartition()` is disabled).
- Ensembles: each replicate fits `round(rho * N)` rows. Restrict `modelNames` (e.g. `c("VII", "VEI", "EII")`) in high dimension to keep parameter counts below the subsample size, and cap `GMax`. For long runs use the command-line pipeline, which fits replicates as independent jobs and reruns a failed one alone with the same seed.
- Save expensive objects (`saveRDS(fit, ...)`); density estimation and ΔP_mc dominate runtime.

## Pitfalls that cost people time

1. **Positions vs labels.** `recoverGroupsv2()` and `mergeComps` refer to positions in `paramsList`, not label values. Map through `sapply(params, "[[", "class")`.
2. **`weights` column order.** A user-supplied `weights` matrix must have columns in the order of first appearance of the labels (`unique(partition)`), not sorted order. Otherwise the cluster densities are silently swapped.
3. **Tiny clusters.** A cluster with only a handful of points makes `constructPmcParamsWeightedPartition` fail ("argument is of length zero" / "a dimension is zero"), even with restricted `modelNames`. Drop it, fold it into its nearest cluster, or use `constructPmcParamsPartition(labels, X, modelNames = "VVI", G = 1)`.
4. **Start fine.** PHM only merges. Start from more clusters (or components) than you believe exist; it cannot split an under-segmented cluster.
5. **Identical replicates.** With N ≤ 2000, `rho = 1` and no `stage2Init`, every ensemble replicate fits the same rows from the same deterministic start, so the ensemble adds nothing. Use `rho < 1` or set `stage2Init`. Above 2000 rows mclust initializes from a random 2000-row subset (`mclust.options("subset")`), so replicates differ.
6. **Tiny top-level merges.** In ensembles a few outlier components often join only at the very top with negligible ΔP_mc and tiny mass. Check the mixture weight of each group (`sapply(groups, function(g) sum(sapply(params[g], function(p) sum(p$prob)))))`) before calling a top split a real cluster.
7. **`recoverPosterior()` (v1) crashes** ("subscript out of bounds") on default `PHM()` output. Read `phm[[K]]$posterior_matrix`, call `addPosteriorMatrix()`, or use `recoverPosteriorv2(convertToPHMv2(phm), ...)`.
8. **`PHM(deltaPmc = )` is ignored.** To reuse a matrix, call `PHMv2(paramsList, deltaPmc = D, scaling = "unscaled")`.
9. **Argument names** are camel-case and differ by function: `subsampSize` in `constructPmcParamsSubAggMclust`, `subsampleSize` in `globalDensityEstimation`, `rho` in the ensemble functions; `mcSamples`, `batchSize`, `numCores` throughout.

## When to use the original `PHM()`

`PHM()` merges by raw ΔP_mc (like `scaling = "unscaled"`) and returns a list indexed by K, with `labels`, `posterior_matrix`, `pmc`, `pmc_change`, `merge_components`, `params` and `pmc_matrix` stored at every level. Reach for it when the user reproduces the paper or the older tutorials, or wants per-K labels and posteriors on the object or the v1 figures (`plotPHMDendrogram`, `plotPHMMatrix`, `plotPmcMatrix`, `plotPHMDistruct`). `thresholdPHM(phm, 0.01)` returns the finest level below a remaining-P_mc threshold. v1 plots take only `PHM()` results and v2 plots only `PHMv2` results; `convertToPHMv2(phm)` converts one way. `phm[[K]]$labels` at coarser K come from the posterior argmax, so they are not exact unions of the input clusters. Details in `references/workflows.md` (section 11) and `references/api.md`.

## Reading results for the user

State what the hierarchy says in plain terms. Early merges (large ΔP_mc, low in the tree) are heavily overlapping clusters that are only distinguishable at fine resolution, candidates for over-splitting. Clusters that stay separate until the top of the tree (tiny ΔP_mc at the merge) are robustly distinct. Say how different clusterers or resolutions agree at coarse levels. Mention that results depend on the density model fitted to each cluster and that Monte Carlo noise matters for very small ΔP_mc (raise `mcSamples` if a conclusion hinges on one). If ground-truth labels exist, `table(truth, lab)` and `mclust::adjustedRandIndex(truth, lab)` are the standard checks.

## Bundled files

- `references/api.md`: every exported function with its signature, defaults, the `paramsList` format, output structures and bundled datasets. Read it when writing or debugging calls.
- `references/workflows.md`: worked recipes (GMM, comparing clusterers, soft weights, choosing a resolution, merge criteria, consensus via subsampling, localized density estimation, large data, Seurat embeddings, the original `PHM()`).
- `references/ensemble.md`: the PHM ensemble in R and from the command line (serial, one machine, Slurm), `config.R` fields, seeds, failures and reruns, scaling advice.
- `references/explorer.md`: exporting to PHM Explorer, scatter-data format, reading the display, troubleshooting.
- `scripts/setup_phm.R`: checks for and installs PHM and its dependencies.
- `scripts/explorer/`: the PHM Explorer exporter (`r/export_phm_json.R`), one-command viewer (`r/view_in_explorer.R`), page template and builder (`src/`).

Cite as: Turfah A, Wen X (2025). Multiresolution clustering of genomic data. bioRxiv, doi:10.1101/2025.06.13.659529. The package is GPL-3.0; questions go to its authors (aturfah/PHM on GitHub).
