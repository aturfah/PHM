## Example to show the different weighting schemes that can be used
library(dplyr)
library(PHM)
library(mclust)
library(ClusterR)
devtools::load_all("~/projects/dissertation/PHM")

# rm(list=ls())

kmeansClusterFunc <- function(X, K=10, seed=1, ...) {
  clust <- KMeans_rcpp(X, K, seed = seed)
  clust$clusters
}

density_naive <- localizedDensityEstimation(
  example_data, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = naiveWeightFunc,
  G=1:6,
  preserveClusters = T,
  verbose=T
)

density_mindist <- localizedDensityEstimation(
  example_data, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = distanceWeightFunc,
  preserveClusters = T,
  G=1:6,
  verbose=T
)

kth_smallest <- function(values, k=10) {
  if (k > length(values)) return(max(values))
  
  sort.int(values, partial=k)[k]
}


density_kdist <- localizedDensityEstimation(
  example_data, 
  partitionFunc = kmeansClusterFunc,
  weightFunc = distanceWeightFunc,
  weightFuncParams=list(aggFunc=kth_smallest),
  preserveClusters = T,
  G=1:6,
  verbose=T
)

kmeans_partition <- kmeansClusterFunc(example_data)

library(ggpubr)
library(ggplot2)

density_levels <- 0.005
plt_naive <- plotDensity2D(density_naive, 
                           example_data, 
                           kmeans_partition, 
                           densityLevels = density_levels,
                           densityLevelWidth=0.8,
                           colorDensity = T)
plt_mindist <- plotDensity2D(density_mindist, 
                              example_data, 
                              kmeans_partition, 
                              densityLevels = density_levels,
                             densityLevelWidth=0.8,
                             colorDensity = T)
plt_kdist <- plotDensity2D(density_kdist, 
                               example_data, 
                               kmeans_partition, 
                               densityLevels = density_levels,
                           densityLevelWidth=0.8,
                           colorDensity = T)


plt_naive_decomp <- plotDensity2D(decomposeParams(density_naive), 
                                   example_data, 
                                   kmeans_partition, 
                                   densityLevels = density_levels,
                                   colorDensity = F)
plt_mindist_decomp <- plotDensity2D(decomposeParams(density_mindist), 
                             example_data, 
                             kmeans_partition, 
                             densityLevels = density_levels,
                             colorDensity = F)
plt_kdist_decomp <- plotDensity2D(decomposeParams(density_kdist), 
                           example_data, 
                           kmeans_partition, 
                           densityLevels = density_levels,
                           colorDensity = F)

ggarrange(
  plt_naive + ggtitle("Cluster Density (Onehot)"),
  plt_mindist + ggtitle("Cluster Density (Minimum Dist.)"),
  plt_kdist + ggtitle("Cluster Density (10-neighbor Dist.)"),
  plt_naive_decomp + ggtitle("Decomposed Density (Onehot)"),
  plt_mindist_decomp + ggtitle("Decomposed Density (Minimum Dist.)"),
  plt_kdist_decomp + ggtitle("Decomposed Density (10-neighbor Dist.)"),
  ncol=3,
  nrow=2
)

