# Localized Density Estimation and Weight Construction

This page describes how to set up the localized density estimation, and constructing custom `weightFunc` methods.
Broadly speaking, localized density estimation adopts a divide-and-conquer strategy for density estimation, constructing density estimates for well-separated regions (as defined by a partition) separately. We incorporate a weighted density estimation scheme to introduce ambiguitry in cluster assignments for observations near cluster boundaries.

For more details on the `localizedDensityEstimation` function, please see the [Global and Localized Density Estimation Page](density_estimation/global_localized.md).

## Implementing Functions for Localized Density Estimation

The localized density estimation function `localizedDensityEstimation` has two mandatory function arguments: `partitionFunc` and `weightFunc`.

### `partitionFunc` and `partitionFuncParams`

Like the name suggests, the `partitionFunc` is the function which produces the partition, defining the distinct regions for the localized density estimation procedure. This function must take two required arguments: `X` (for the data matrix) and `seed` for setting the seed. It must return a vector with partition labels, where cluster labels are integers (1, 2, 3, ..., *K*). An example for `k`-means using the Silhouette coefficient to select the number of clusters from a range of *K* values (default is `2:10`) is provided below.

```
library(ClusterR)

## k-means selecting K via silhouette coefficient
silhouetteClusterFunc <- function(X, K_range=2:10, num_init=1, seed=1) {
  max_sil <- -Inf
  max_soln <- NULL
  for (K in K_range) {
    clust <- KMeans_rcpp(X, K, num_init = num_init, seed = seed)$clusters
    sil <- silhouette_of_clusters(X, clust)$silhouette_global_average
    
    if (sil > max_sil) {
      max_sil <- sil
      max_soln <- clust
    }
  }

  max_soln
}

## Function call specifying this function
dens <- localizedDensityEstimation(
    X, 
    partitionFunc=silhouetteClusterFunc,
    <Remainder of arguments...>
)
```

The `localizedDensityEstimation` method allows for named parameters to be passed and specified via the `partitionFuncParams` argument, which accepts a list of arguments to be passed to the `partitionFunc`. If in the instance below we wanted to have *k*-means five different initializations for the clustering, and change the $K$ range from `2:20`, this can be accomplished with the following function call.

```
dens <- localizedDensityEstimation(
    X, 
    partitionFunc=silhouetteClusterFunc,
    partitionFuncParams=list(num_init=5, K_range=2:20),
    <Remainder of arguments...>
)
```

### `weightFunc` and `weightFuncParams`

Once the partition has been constructed, it is then used to construct the weights. The `weightFunc` must accept two required arguments: `X` (for the data matrix) as well as `partition` (for the output of `partitionFunc`). It must return an $N \times K$ weight matrix (i.e., rows sum to 1) where entry *(i, j)* is the proportion of weight for observation *i* assigned to cluster *j*.

The PHM package implements two weight functions: `naiveWeightFunc` and `distanceWeightFunc`. `naiveWeightFunc` assigns an observation a weight of 1 for the cluster to which it is assigned in the partition and 0 otherwise. The `distanceWeightFunc` assigns weights based on the minimum distance (i.e., nearest neighbor) to observations within clusters. Please see the function documentation for additional details; the implementation for `naiveWeightFunc` is shown below to demonstrate the basic structure.

```
naiveWeightFunc <- function(X, partition) {
  one_hot <- matrix(0, nrow = length(partition), ncol = max(partition))
  one_hot[cbind(seq_along(partition), partition)] <- 1
  
  colnames(one_hot) <- 1:max(partition)
  one_hot
}
```

In the same vein as `partitionFunc`, additional arguments can be provided to `weightFunc` via the `weightFuncParams` argument. In the case of the `distanceWeightFunc`, if we wanted to use the $10^{th}$ nearest neighbor instead of the minimum distance we could do this by specifying the `aggFunc` parameter in the manner shown below.

```
kth_smallest <- function(values, k=10) {
  if (k > length(values)) return(max(values))

  sort.int(values, partial=k)[k]
}

## Use the k-means with silhouette clustering function
dens <- localizedDensityEstimation(
    X, 
    partitionFunc=silhouetteClusterFunc,
    partitionFuncParams=list(num_init=5),
    weightFunc=distanceWeightFunc,
    weightFuncParams=list(aggFunc=kth_smallest)
)
```

## Example 

Here we will illustrate the different weighting schemes and partition functions on a simulated 2D dataset with the following structure:

- One third of the observations belong to two Gassian clusters arranged in an "X"
- One third of the observations are arranged in a ring
- One third of the observations belong to a Gaussian cluster inside of the ring

<center>
<b>IMAGE GOES HERE</b>
</center>

Here we will use a simple $k$-means partitioning function with pre-specified $K=10$

```
library(ClusterR)

kmeansClusterFunc <- function(X, K=10, seed=1, ...) {
  clust <- KMeans_rcpp(X, K, seed = seed)
  clust$clusters
}
```

We will also compare the estimated densities using three different weighting schemes: Naive (`naiveWeightFunc`), Distance-based weighting (`distanceWeightFunc`), as well as the distance-based weighting using the $10^{th}$ nearest neighbor distance (as described above).

For the cluster-specific densities we will allow up to 5 components per cluster. The cluster grouping of the components.

```
library(PHM)
library(mclust)
data("density_example", package="PHM")

density_naive <- localizedDensityEstimation(
  density_example, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = naiveWeightFunc,
  preserveClusters = T,
  G=1:6
)

density_mindist <- localizedDensityEstimation(
  density_example, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = distanceWeightFunc,
  preserveClusters = T,
  G=1:6
)

kth_smallest <- function(values, k=10) {
  if (k > length(values)) return(max(values))
  
  sort.int(values, partial=k)[k]
}

density_kdist <- localizedDensityEstimation(
  density_example, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = distanceWeightFunc,
  weightFuncParams=list(aggFunc=kth_smallest),
  preserveClusters = T,
  G=1:6
)
```

Below we visualize the cluster-estimated densities for the three density estimation procedures (top row). We also visualize the mixture density discarding the cluster labels (bottom row), using the `decomposeParams` function.

```
library(ggpubr)
library(ggplot2)

density_levels <- 0.005
plt_naive <- plotDensity2D(density_naive, 
                           density_example, 
                           kmeans_partition, 
                           densityLevels = density_levels,
                           densityLevelWidth=0.8,
                           colorDensity = T)
plt_mindist <- plotDensity2D(density_mindist, 
                              density_example, 
                              kmeans_partition, 
                              densityLevels = density_levels,
                             densityLevelWidth=0.8,
                             colorDensity = T)
plt_kdist <- plotDensity2D(density_kdist, 
                               density_example, 
                               kmeans_partition, 
                               densityLevels = density_levels,
                           densityLevelWidth=0.8,
                           colorDensity = T)


plt_naive_decomp <- plotDensity2D(decomposeParams(density_naive), 
                                   density_example, 
                                   kmeans_partition, 
                                   densityLevels = density_levels,
                                   colorDensity = F)
plt_mindist_decomp <- plotDensity2D(decomposeParams(density_mindist), 
                             density_example, 
                             kmeans_partition, 
                             densityLevels = density_levels,
                             colorDensity = F)
plt_kdist_decomp <- plotDensity2D(decomposeParams(density_kdist), 
                           density_example, 
                           kmeans_partition, 
                           densityLevels = density_levels,
                           colorDensity = F)

ggarrange(
  plt_naive + ggtitle("Cluster Density (Onehot)"),
  plt_mindist + ggtitle("Cluster Density (Minimum Dist.)"),
  plt_kdist + ggtitle("Cluster Density (10-neighbor Dist.)"),
  plt_naive_decomp + ggtitle("Decomposed (Onehot)"),
  plt_mindist_decomp + ggtitle("Decomposed (Minimum Dist.)"),
  plt_kdist_decomp + ggtitle("Decomposed (10-neighbor Dist.)"),
  ncol=3,
  nrow=2
)
```

<center>
<b>IMAGE GOES HERE</b>
</center>

