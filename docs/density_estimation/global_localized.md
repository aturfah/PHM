# Density Estimation

The PHM package implements two methods for density estimation:

- **Global** density estimation involves using the full sample and estimating the density using a Gaussian Mixture Model (GMM).
- **Localized** density estimation involves first partitioning the data into clusters and then estimating each cluster-specific density (via a GMM) separately. 

## Global Density Estimation

The PHM package implements the global density estimation via the `globalDensityEstimation()` function. Each replicate is fit by `fitDensityReplicate()`, the same function used by the PHM ensemble pipeline (see [Running a PHM Ensemble](ensemble.md)): it selects a covariance model and number of components `G` by BIC with `Mclust`, and widens the `G` range automatically when the best `G` is the largest one tried. Arguments such as `G`, `modelNames`, `GStep`, `GMax` and `stage2Init` are passed through to `fitDensityReplicate()`.

This function also implements the ensemble density estimation procedure (via subagging). The number of ensemble replicates is controlled by the parameter `M`. Specifying either the `rho` or `subsampleSize` parameters control the number of samples used to estimate the density in the ensemble replicates; each replicate uses `round(rho * N)` rows. The replicates are pooled with equal weights. With the default seeds, `globalDensityEstimation(X, M, rho)` returns the same components as `fitPHMEnsemble(X, M, rho)$paramsList`.

Some example function calls are provided below. For a more detailed description of the function and its arguments, please see the package documentation.

```
## A single mixture fit to all of example_data (G = 1:25 to start, widened if needed)
globalDensityEstimation(X=example_data)

## Control the number of components and only consider spherical covariance structure
globalDensityEstimation(X=example_data, 
                        G=1:20, 
                        modelNames=c("EEE"))

## Fit multiple GMMs based on the entire dataset
globalDensityEstimation(X=example_data, 
                        M=5)

## Ensemble density estimation via subagging (no replacement)
globalDensityEstimation(X=example_data, 
                        M=5,
                        rho=0.5)
globalDensityEstimation(X=example_data, 
                        M=5,
                        subsampleSize=1500)

## Resample the entire dataset (with replacement)
globalDensityEstimation(X=example_data, 
                        M=5,
                        rho=1,
                        replace=T)
```

## Localized Density Estimation

The PHM package implements a localized density estimation procedure via a three-step process. First, the observations are partitioned according to some clustering method (e.g., *k*-means). Second, each observation is assigned a weight of belonging to the inferred clusters. Finally, these weights are used to estimate the cluster-specific density (in the form of a GMM) via a weighted EM procedure. 


This overall procedure is implemented in the `localizedDensityEstimation()` function. This function takes the partitoning procedure (`partitionFunc`) as well as the procedure to determine the cluster-specific weights (`weightFunc`) as arguments as well. The PHM package implements two methods to compute the cluster-specific weights: `naiveWeightFunc()` which gives the assigned cluster a weight of 1 with all other clusters a weight of 0, as well as the distance-based weighting scheme described in [Turfah and Wen (2026)](#). For computational efficiency, a `threshold` parameter is provided to ignore observations with weights falling below the specified value when performing density estimation for that cluster.

Additional parameters to these functions can be provided in the form of a list using the corresponding `...FuncParams` argument. The density estimation is performed using the `Mclust` function; as such any valid arguments for that function can also be provided here.

The cluster-specific distributions are estimated as a mixture density, leading to the representation as a mixture of mixtures. If the overall mixture distribution (i.e., ignoring the cluster labels and treating all mixture components separately) is of interest then this can be specified via the `preserveClusters` parameter.

This function implements an ensemble density estimation procedure (with subagging) with the same parameterization as `globalDensityEstimation()` described above. In this case, the procedure described above is applied across the ensemble replicates, each producing its own partition and corresponding mixture density. 

Some example calls for the localized density estimation are provided below.

```
library(dbscan) ## For HDBSCAN
library(ClusterR) ## For k-means

## Clustering Functions
kmeansClusterFunc <- function(X, K, num_init=1, seed=1, ...) {
  clust <- KMeans_rcpp(X, K, num_init=num_init, seed = seed)
  clust$clusters
}

hdbscanClusterFunc <- function(X, minPts=5, singletonNoise=T, ...) {
  clust <-hdbscan(X, minPts = minPts)$cluster

  ## Turn noise points into singleton clusters
  if (singletonNoise) {
    max_K <- max(clust)
    clust[clust == 0] <- max_K + 1 + seq_along(which(clust==0))
  }
  
  clust
}

## Use HDBSCAN to partition the entire sample with naive weighting; preserve clusters for density
localizedDensityEstimation(example_data,
                           partitionFunc=hdbscanClusterFunc, 
                           weightFunc=naiveWeightFunc,
                           preserveClusters = T)

## Distance-based weighted density estimation based on k-means; ignore cluster labels
localizedDensityEstimation(example_data,
                           partitionFunc=kmeansClusterFunc, 
                           partitionFuncParams=list(K=7)
                           weightFunc=distanceWeightFunc, 
                           preserveClusters = F)

## Ensemble of k-means based on 50% of the data
localizedDensityEstimation(example_data,
                           partitionFunc=kmeansClusterFunc, 
                           partitionFuncParams=list(K=7)
                           weightFunc=distanceWeightFunc, 
                           M=5,
                           rho=0.5,
                           preserveClusters = F)
```
