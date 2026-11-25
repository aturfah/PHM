# PHM API reference

Signatures and defaults below were taken from the installed package (`args()`) at PHM 1.1.0, GitHub main after the ensemble pipeline was merged (November 2026). Where roxygen text and code disagree, this file follows the code. Check `args(fn)` if the installed version is different.

## Contents

- [Parameter list format (`paramsList`)](#parameter-list-format-paramslist)
- [Density estimation from raw data](#density-estimation-from-raw-data)
- [Building `paramsList` from a clustering](#building-paramslist-from-a-clustering)
- [PHMv2 (recommended default)](#phmv2-recommended-default)
- [The original PHM()](#the-original-phm)
- [Working with PHM() results](#working-with-phm-results)
- [P_mc primitives](#p_mc-primitives)
- [Visualization](#visualization)
- [Under-documented utilities](#under-documented-utilities)
- [Bundled datasets](#bundled-datasets)
- [Known bugs and doc gaps](#known-bugs-and-doc-gaps)

## Parameter list format (`paramsList`)

Everything downstream of density estimation consumes a list with one element per **cluster**. Each cluster's density is itself a Gaussian mixture with M components (M = 1 is a single Gaussian):

```r
list(
  list(prob = <numeric, length M>,        # mixture weights; sum over ALL clusters must be 1
       mean = <D x M matrix>,             # one column per component
       var  = <D x D x M array>,          # one covariance slice per component
       class = <label>),                  # cluster label; shown on plot axes and used to name merged clusters
  ...
)
```

- `prob` is the cluster's share of the population times each component's share inside the cluster, so `sum(unlist(lapply(params, `[[`, "prob")))` is 1. PHM rescales with a warning if not.
- For a 1-D component keep `mean` as a 1 × M matrix and `var` as a 1 × 1 × M array (`drop = FALSE` when subsetting).
- Ensemble, `globalDensityEstimation()` and Mclust lists have M = 1 per element; partition constructors and `localizedDensityEstimation(preserveClusters = TRUE)` produce M ≥ 1.

## Density estimation from raw data

### The ensemble functions

```r
fitPHMEnsemble(X, M = 5, rho = 1, G = 1:25, modelNames = NULL,
               replicateWeighting = c("equal", "likelihood", "bic"), seedBase = 1000,
               numCores = 2, mcSamples = 2e5, batchSize = 1e5,
               stage2Init = NULL, nRestartTries = 5, GStep = 5, GMax = 100,
               extendVerbose = TRUE, label = "example")
```

Runs `fitDensityReplicate()` for `idx = 1..M` with seeds `phmSeeds(seedBase, M)`, one after another, then `combineDensityReplicates()`. Returns a **`PHMv2` object** (fields below), with the pooled components in `$paramsList`. `numCores` is only for the PHMv2 Monte Carlo step, not for the replicates.

```r
fitDensityReplicate(X, idx, seed, rho = 1, G = 1:25, modelNames = NULL, stage2Init = NULL,
                    nRestartTries = 5, GStep = 5, GMax = 100, extendVerbose = TRUE,
                    label = "example", replace = FALSE, ...)
```

The map step. `set.seed(seed)`, draws `round(rho * N)` rows (without replacement unless `replace = TRUE`), then Stage 1: `mclustBIC` over `G` and `modelNames`; when the best G is the top of the range it extends by `GStep` up to `GMax` (an extension that fails keeps the best fit so far and says so). Optional Stage 2 (`stage2Init` one of `"kmeans"`, `"kmeans++"`, `"random"`, `"randomFromData"`): re-runs EM for the chosen model and G from a random start, retrying up to `nRestartTries` times, falling back to Stage 1. `...` goes to `mclustBIC`. Returns `list(params, loglik, bic, idx, G, modelName)`; `params` is a paramsList whose `class` values are component indices.

```r
combineDensityReplicates(replicateResults, replicateWeighting = c("equal", "likelihood", "bic"),
                         numCores = 2, mcSamples = 2e5, batchSize = 1e5, label = "example")
```

The reduce step. `replicateResults` is the list of `fitDensityReplicate()` returns (M = its length). Weights each replicate (`equal` = 1/M; `likelihood` = softmax of log-likelihoods; `bic` = softmax of BIC/2), multiplies its component probabilities, relabels classes `"<component>_<idx>"`, pools in replicate order, and runs `PHMv2()` (alpha scaling). Returns the `PHMv2` object.

```r
phmSeeds(seedBase, M)    # set.seed(seedBase); sample.int(1e7, M)
```

The replicate seeds. Replicate `i` always gets `phmSeeds(seedBase, M)[i]`, so it can be refit alone with the same result. Changing M keeps the first seeds only if seedBase is unchanged (the draw is sequential).

### Global and localized density estimation

```r
globalDensityEstimation(X, M = 1, rho = NULL, subsampleSize = NULL, replace = FALSE,
                        seeds = phmSeeds(1000, M), numCores = 1, saveDir = NULL,
                        prefix = "global", overwrite = FALSE, verbose = FALSE, ...)
```

Returns a **paramsList** (not a PHMv2 object): the pooled components of M `fitDensityReplicate()` fits, weighted equally. `rho` or `subsampleSize` set rows per replicate (default all rows). `...` goes to `fitDensityReplicate()` (`G`, `modelNames`, `GStep`, `GMax`, `stage2Init`). `numCores` fits replicates in parallel (forked). `saveDir` caches each replicate as `<prefix>_ssize_<n>_Repl_<idx>_replicate.rds` and reuses it unless `overwrite`. With default seeds, `globalDensityEstimation(X, M, rho)` equals `fitPHMEnsemble(X, M, rho)$paramsList`.

```r
localizedDensityEstimation(X, partitionFunc = NULL, partitionFuncParams = list(),
                           weightFunc = NULL, weightFuncParams = list(), weightThreshold = 0.005,
                           preserveClusters = FALSE, M = 1, rho = NULL, subsampleSize = NULL,
                           replace = FALSE, seeds = 1:M, numCores = 1, saveDir = NULL,
                           prefix = "global", overwrite = FALSE, verbose = FALSE, ...)
```

Per replicate: `labels <- partitionFunc(X = Xsub, seed = seeds[idx], <partitionFuncParams>)`, `weights <- weightFunc(X = Xsub, partition = labels, <weightFuncParams>)` (N × K), then a weighted Mclust per cluster, ignoring points whose weight is below `weightThreshold`. `...` goes to Mclust (`G`, `modelNames`). Returns a paramsList: one element per cluster (a mixture, class `"<replicate>_<label>"`) with `preserveClusters = TRUE`, or one element per Gaussian component (class `"<replicate>_<label>_<m>"`) with `FALSE`. Probabilities are renormalized over the whole list.

Requirements that are easy to miss: the partition function is always called with a `seed` argument, so give it a `seed` (or `...`) parameter; and labels must be the integers `1..K` with no gaps, because weights are looked up as `weights[, k]` by label value (recode with `as.integer(factor(labels))` inside the function; HDBSCAN noise label 0 must be turned into real clusters).

Weight functions:

- `naiveWeightFunc(X, partition)`: one-hot.
- `distanceWeightFunc(X, partition, removeSelfDistance = FALSE, scaling = 1, aggFunc = min)`: `w_ij ∝ exp(-scaling * d(x_i, C_j))`, with d the Euclidean distance from x_i to the nearest member of cluster j (`aggFunc = mean` etc. changes the linkage). Columns follow `sort(unique(partition))`. O(N²) memory per cluster block.

## Building `paramsList` from a clustering

| Function | Use |
|---|---|
| `constructPmcParamsMclust(mclustObj, singleElement = FALSE)` | One cluster per `Mclust` component (`class` = component index). `singleElement = TRUE` returns all components as one element. |
| `constructPmcParamsPartition(partition, data, ...)` | Naive: for each label, `Mclust` on that cluster's points only. `...` goes to `Mclust` (`G`, `modelNames`). Fast, tolerant of small clusters, under-estimates tails. |
| `constructPmcParamsWeightedPartition(partition, data, weights = NULL, threshold = 1e-4, linkFunc = min, verbose = FALSE, ...)` | Weighted EM per cluster. Default weights `w_ij = exp(-d_ij) / Σ_k exp(-d_ik)`, `d_ij` the squared Euclidean distance from point i to its nearest point in cluster j (`linkFunc`). `threshold` drops points with smaller weight for a cluster. `weights`: optional N × K matrix with columns in `unique(partition)` order. `...` goes to `Mclust`. Errors on data.frames and on clusters with only a few points. |
| `constructPmcParamsSubAggMclust(data, replicates, subsampSize = NULL, G = NULL, saveDir = NULL, prefix = "subAggMclust_", verbose = FALSE, numCores = 1, seeds = NULL, ...)` | Older subsample-and-aggregate GMM: `Mclust` on `replicates` subsamples of `subsampSize` rows (without replacement; whole data if NULL; `G` defaults to 1:10), components pooled with probabilities divided by `replicates`, labels `"<j>_<replicate>"`. Kept for the published analyses (e.g. HIHA tutorial); prefer `fitPHMEnsemble()`/`globalDensityEstimation()` for new work, which add automatic G extension, seeds per replicate and the command-line pipeline. |
| `constructPmcParamsSubAggPartition(...)` | Disabled: the current code stops with "This method does not work; use constructPmcParamsLocalizedEnsemble", and that function is internal. Use `localizedDensityEstimation()` with `M`/`rho`, or the weighted constructor on a subsample. |

## PHMv2 (recommended default)

```r
PHMv2(paramsList = NULL, deltaPmc = NULL, scaling = c("alpha", "unscaled", "average"),
      monteCarlo = TRUE, numCores = 1, verbose = FALSE, ...)
```

Provide `paramsList`, a precomputed square `deltaPmc`, or both (alpha needs `paramsList`; both skips the Monte Carlo). `...` goes to `computeMonteCarloDeltaPmcMatrix()` (`mcSamples` default 1e6, `batchSize` default `mcSamples`). Merge rule: pick the pair maximizing 2·ΔP_mc divided by 1 (`unscaled`), the product of the subtrees' numbers of original clusters (`average`), or the product of their total mixture weights (`alpha`).

| Field | Meaning |
|---|---|
| `deltaPmc` | Original K × K ΔP_mc matrix (zero diagonal; sums to P_mc) |
| `paramsList` | The input parameters (if given) |
| `mergeComps` | List of length K; `mergeComps[[k]]` is `c(i = , j = )`, the positions merged to go from k to k − 1 (`NA` at index 1) |
| `mergeValues` | Value of the (possibly scaled) merge criterion for each merge |
| `mergeDeltaPmc` | Unscaled 2·ΔP_mc of each merge; `cumsum(mergeDeltaPmc)[k]` is the remaining P_mc at k |
| `mergeCriterion` | The `scaling` used |

Readers:

- `recoverGroupsv2(phm, k)`: list of original-cluster positions per cluster at k.
- `recoverPosteriorv2(phm, k, posterior)`: sums the columns of an N × K posterior matrix (e.g. `computePosteriorProbMatrix(phm$paramsList, X)`) within each group; N × k.
- `recoverParamsv2(phm, k, paramsToKeep = c("prob", "mean", "var", "class"))`: merged cluster parameters at k (mixtures).
- `recoverDeltaPmcv2(phm, k)`: k × k ΔP_mc matrix.

`pmc_by_k[K] == sum(recoverDeltaPmcv2(phm, K))` for every scaling, and for `"unscaled"` it equals `sapply(phm1, "[[", "pmc")` of the corresponding `PHM()` run.

## The original PHM()

```r
PHM(mclustObj = NULL, paramsList = NULL, partition = NULL, data = NULL,
    deltaPmc = NULL,            # not used
    verbose = TRUE, computePosterior = TRUE,
    partitionWeightedDensity = TRUE, storeDeltaPmc = TRUE, storeParams = TRUE,
    partitionModel = "VVI", partitionMaxComponents = 10,
    mc = TRUE, ...)             # ... -> computeMonteCarloDeltaPmcMatrix (mc=TRUE) or computeDeltaPmcMatrix (mc=FALSE)
```

Input resolution: `mclustObj` wins; otherwise `partition` (builds the weighted paramsList with `modelNames = partitionModel`, `G = 1:partitionMaxComponents` unless `partitionWeightedDensity = FALSE`); otherwise `paramsList`. `data` is required except with `mclustObj`.

Returns a list of length K. `phm[[k]]`:

| Field | Meaning |
|---|---|
| `clusters` | k |
| `labels` | Observation labels at k. Input labels at the starting K; merged clusters joined with `"|"`; coarser levels from the merged posterior argmax. |
| `posterior_matrix` | N × k posteriors (NULL if `computePosterior = FALSE`) |
| `pmc` | Remaining P_mc at k (0 at k = 1) |
| `pmc_change` | Drop in P_mc from the merge that produced k from k + 1 (NA at the start); 2 × the merged pair's ΔP_mc |
| `merge_components` | (i, j), positions at level k + 1 that were merged |
| `pmc_components` | Number of original clusters in the newly merged cluster |
| `min_merge_pmc` | Minimum pairwise ΔP_mc between original clusters inside the merged subtree |
| `params` | Merged cluster parameters (stored when `storeParams = TRUE`; always at the starting K) |
| `pmc_matrix` | ΔP_mc matrix at k (stored when `storeDeltaPmc = TRUE`; always at the starting K) |

## Working with PHM() results

| Function | Notes |
|---|---|
| `thresholdPHM(phm, threshold = 0.01)` | Returns `phm[[k]]` for the largest k whose `pmc` is below `threshold`. |
| `addPosteriorMatrix(phm, data, initK = length(phm))` | Fills `posterior_matrix` and `labels` at every level. Use after `computePosterior = FALSE`. |
| `recoverPHMParams(phm, K, paramsToKeep = ...)` | Merged parameters at K, rebuilt from the nearest stored level. Requires K < length(phm). |
| `recoverDeltaPmcMatrix(phm, K)` | ΔP_mc matrix at K. Requires K < length(phm). |
| `recoverPosterior(phm, K, data = NULL, computePosterior = FALSE)` | Broken on default output (see Known bugs). |
| `computePosteriorProbMatrix(paramsList, data)` | N × K posterior cluster probabilities under the fitted densities. |
| `convertToPHMv2(phm)` | Converts a `PHM()` result so the v2 readers and plots accept it. |

## P_mc primitives

| Function | Notes |
|---|---|
| `computeMonteCarloPmc(paramsList, mcSamples = 1e5, batchSize = mcSamples, numCores = 1, verbose = FALSE)` | Monte Carlo P_mc = mean over draws of Σ_j π_j(1 − π_j). |
| `computeMonteCarloDeltaPmcMatrix(paramsList, mcSamples = 1e6, batchSize = mcSamples, numCores = 1, verbose = FALSE)` | K × K ΔP_mc matrix as `crossprod(posterior) / mcSamples`, zero diagonal. Batches are spread over `numCores` forked workers; a failed worker stops the call with an error. |
| `computePmc(paramsList, integralControl = list())` | Cubature P_mc (`$integral`). Low dimension only. |
| `computeDeltaPmcMatrix(paramsList, integralControl = list())` | Cubature ΔP_mc matrix (default `method = "cuhre"`, `relTol = 1e-6`, `maxEval = 1e6`). One integral per pair; slow. |

## Visualization

All return ggplot objects.

```r
plotPHMv2Dendrogram(phmObj, initK = NULL, scaleHeights = c("d*", "log10", "unscaled", "pmcdist"),
    threshold = 0, colors = NULL, displayAxis = c("box", "label", "index", "none"),
    displayAxisSize = NULL, displayValue = FALSE, displayValueSize = NULL, colorAxis = NULL,
    groupProbs = NULL, groupColorMax = "black", groupColorMin = "lightgray")

plotPHMv2Heatmap(phmObj, initK = NULL, colors = NULL, displayAxis = c("box", "label", "index", "none"),
    displayAxisSize = NULL, colorAxis = NULL, gridColor = "black", fillLimits = NULL,
    fillScale = c("d*", "log10", "pmcdist"), legendPosition = "none", subsetClusters = NULL, scaleAxis = FALSE)
```

`"d*"` requires an alpha run; `"pmcdist"` is not available for alpha runs. `initK` starts the picture at a coarser level (leaves = the clusters at initK). `subsetClusters` (positions of original clusters) zooms into part of the tree. `groupProbs` (length = number of leaves) colors branches by a per-cluster value, e.g. the fraction of cells from a known group; merged branches show the weighted average. `threshold` draws merges below it dashed.

v1 plots (take a `PHM()` result):

```r
plotPHMDendrogram(phm, colors = NULL, scaleHeights = c("log10", "unscaled", "pmcdist"),
    heightValue = c("merge", "min"), threshold = 0, suppressLabels = FALSE,
    mergeLabels = c("delta", "pmc", "percent"), ..., displayAxis = c("box", "label", "index", "none"), ...)
plotPHMMatrix(phm, colors = NULL, displayAxis = ..., gridColor = "black", fillLimits = NULL,
    fillScale = c("log10", "pmcdist"), legendPosition = "none")
plotPmcMatrix(phm, K = length(phm), colors = NULL, displayAxis = ..., visScale = c("absolute", "percent"),
    visSize = 2, visThreshold = 1e-3, visDigits = 3)
plotPHMDistruct(phm, K = length(phm), colors = NULL, labels = NULL, axisTextSize = 6, partition = FALSE)
```

Takes a paramsList:

```r
plotPmc2D(paramsList, data, partition, colors = brewer.pal(length(paramsList), "Paired"),
    xlim = NULL, ylim = NULL, suppressPmc = FALSE, numPmcPatches = 200, logPmcThreshold = -4,
    logPmcMidpoint = 0.6 * logPmcThreshold, PmcColor = "#222", suppressDensity = FALSE,
    densityLevels = c(1e-2, 1e-1), densityLevelWidth = 0.3, suppressObservations = FALSE,
    pointSize = 0.5, labelSize = 2, textSize = 9, legendPosition = "none")
```

2-column data only. Shades where point-wise P_mc is large and overlays density contours. Pass `PmcColor = "#222222"` (see Known bugs). Lower `numPmcPatches` for speed.

- Palettes: axis boxes use Paired (K ≤ 12); `plotPHMDistruct` defaults to Set1 (9 colors). Pass `colors` for more.
- `plotPHMMatrix` `fillLimits` are on the raw `pmc_change` scale; `plotPHMv2Heatmap` `fillLimits` on the raw `mergeValues` scale.

## Under-documented utilities

Roxygen text still says "TODO"; treat as experimental.

- `computePairwisePmcMatrix(paramsList, mcSamples, numCores = 1, threshold = 0, thresholdTol = 0.01, equalProbs = FALSE, verbose = TRUE)`: pairwise P_mc between single-Gaussian clusters as two-cluster problems, skipping pairs below `threshold`. Returns `list(mat, ...)`.
- `consolidateParams(paramsList, pmcMatrix, threshold)`: merges clusters whose pairwise P_mc exceeds `threshold` into single Gaussians (moment matching). Single-component clusters only.
- `decomposeParams(paramsList)`: splits each cluster's mixture into one element per Gaussian (labels `"<class>_<m>"`).

## Bundled datasets

`data("name", package = "PHM")`, or just use the name after `library(PHM)` (lazy data). All N × 2 numeric matrices without labels:

- `crater` (1750): a dense disk inside a ring.
- `density_example` (3000): two crossing elongated Gaussians, plus a ring around a round Gaussian.
- `double_ring` (5000): two concentric rings, rows 1–2500 inner (radius 2–3), 2501–5000 outer (6–7).
- `moons` (1000): two interlocking half-rings, rows 1–500 and 501–1000.
- `nested_diamonds` (6400): 64 Gaussian clusters in three nested levels (4 × 4 × 4), 100 rows each.

## Known bugs and doc gaps

- `recoverPosterior(phm, K)` crashes ("subscript out of bounds") on default `PHM()` output: its loop `min_idx:(K+1)` runs upward when the posterior is stored at K. The same loop shape is in `recoverDeltaPmcMatrix` and `recoverPHMParams`, which are guarded by `stopifnot(K < length(phm))`.
- `plotPmc2D()`'s default `PmcColor = "#222"` is rejected by ggplot2 ("Malformed colour string ... Must contain either 6 or 8 hex values"). Pass a six-digit color.
- `PHM(deltaPmc = )` is accepted but not used (now documented as such).
- `computeMonteCarloPmc`/`computePmc` roxygen formulas show a factor of 2 per term that the code does not use. Code: P_mc = Σ_j π_j(1 − π_j); the ΔP_mc matrix sums to P_mc.
- `localizedDensityEstimation()` defaults `seeds = 1:M` while `globalDensityEstimation()` defaults `seeds = phmSeeds(1000, M)`, and both default `prefix = "global"`; pass `prefix` explicitly when caching both in one `saveDir`.
- `constructPmcParamsSubAggPartition()` is exported but stops immediately; its error message names `constructPmcParamsLocalizedEnsemble`, which is not exported. The supported replacement is `localizedDensityEstimation()`.
- Several older help pages (`recoverPosterior`, `recoverDeltaPmcMatrix`, the utilities above) have "TODO" details.
