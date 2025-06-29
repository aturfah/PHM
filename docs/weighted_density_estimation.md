# Weighted Density Estimation

We consider the hierarchical clustering partition for k = 10 for the double ring example. We first load the necessary packages and cluster the data 

```
library(PHM)
library(mclust)
library(RColorBrewer)


## Start from hierarchical clustering partition with k = 10
data("double_ring", package = "PHM")
hcl <- hclust(dist(double_ring), method = "ward.D2")
hcl_labels <- cutree(hcl, 10)
```

The data and partition is visualized below

<center>
<img src="figures/weighted_density/data_plot.png" alt="Alt Text" width="600" height="300">
</center>


Based on this partition, we compare the results for estimating the cluster-specific densities via GMM in two ways

1. *Unweighted*: Use only the observations in a given cluster to estimate its cluster-specific density
2. *Weighted*: Each observation is assigned a cluster-specific weight based on its minimum distance to a given cluster. Specifically, the "distance" for an observation $i$ to a given cluster $k$ is defined as
$$
    d_{x_i \mapsto C_k} = \min_{j \in C_k} d(x_i, x_j)
$$
and the weight for observation $i$ to cluster $k$ is
$$
    w_{ik} = \frac{e^{-d_{x_k \mapsto C_k}}}{ \sum_{l=1}^K e^{-d_{x_i \mapsto C_l}} }
$$
