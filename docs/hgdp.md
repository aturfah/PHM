# Human Genome Diversity Project

Load the necessary packages

```
library(PHM)
library(mclust)
library(ggplot2)
library(PCAtools)
library(stringr)
library(dplyr)
library(densitycut)
library(ClusterR)
library(Seurat)
library(dbscan)
```

Prepare and preprocess the data

```
data_raw <- read.table("data_files/hgdp.txt") %>%
    data.frame()

allele_dosage_matrix <- apply(data_raw[, 4:ncol(data_raw)], 2, function(vec) {
ref_allele <- substr(vec[1], 1, 1) 

str_count(vec, ref_allele)
})
sampling_labels <- data_raw[, 2]
geographic_labels <- factor(data_raw[, 3],
                            levels=c("EUROPE",
                                    "CENTRAL_SOUTH_ASIA",
                                    "MIDDLE_EAST", 
                                    "EAST_ASIA",
                                    "OCEANIA",
                                    "AMERICA",
                                    "AFRICA"), ordered=T)
geographic_levels <- levels(geographic_labels)
geographic_levels <- case_when(
geographic_levels == "EUROPE" ~ "Europe",
geographic_levels == "CENTRAL_SOUTH_ASIA" ~ "\nC/S Asia",
geographic_levels == "MIDDLE_EAST" ~ "Middle East",
geographic_levels == "EAST_ASIA" ~ "\nEast Asia",
geographic_levels == "OCEANIA" ~ "Oceania",
geographic_levels == "AMERICA" ~ "\nAmerica",
geographic_levels == "AFRICA" ~ "Africa",
TRUE ~ geographic_levels
)
levels(geographic_labels) <- geographic_levels

PC_RES <- prcomp(allele_dosage_matrix)
PC_MATRIX <- PC_RES$x  
```

PCs chosen based on the elbow point

```
sum(PC_RES$sdev[1:200]^2) / sum(PC_RES$sdev^2) ## ~85%
NUM_PCS <- findElbowPoint(PC_RES$sdev[1:200]^2)
DATA <- PC_MATRIX[, 1:NUM_PCS]
```

Fit models and run PHM

```
set.seed(1)
## GMM
mcl <- Mclust(DATA, warn=T, G=1:20)
phm_mcl <- PHM(mclustObj=mcl,
                data=DATA,
                batchSize=1e5)

## DensityCut
dccl <- DensityCut(DATA, show.plot = F)
dc_labels <- match(dccl$cluster, sort(unique(dccl$cluster)))
dc_params <- constructPmcParamsWeightedPartition(dc_labels, DATA)
phm_dc <- PHM(paramsList=dc_params, data=DATA, partition=dc_labels)

## Calculate the Silhouette Scores for k-means and hierarchical clustering
hcl <- hclust(dist(DATA), method="ward.D2")
cluster_silhouettes <- sapply(2:20, function(k) {
    hcl_labels <- cutree(hcl, k)
    kcl_labels <- KMeans_rcpp(DATA, k, num_init = 5)$clusters
    c(K=k,
    sil_hcl=silhouette_of_clusters(DATA, hcl_labels)$silhouette_global_average,
    sil_kcl=silhouette_of_clusters(DATA, kcl_labels)$silhouette_global_average)
})

## k-means
kcl_k_sil <- cluster_silhouettes[1, which.max(cluster_silhouettes[3, ])]
kcl_sil <- KMeans_rcpp(DATA, kcl_k_sil, num_init = 5)

kcl_sil_params <- constructPmcParamsWeightedPartition(kcl_sil$clusters, DATA)
phm_ksil <- PHM(paramsList=kcl_sil_params, data=DATA, partition=kcl_sil$clustkers)

## Hierarchical Clustering
hcl_k_sil <- cluster_silhouettes[1, which.max(cluster_silhouettes[2, ])]
hcl_sil_params <- constructPmcParamsWeightedPartition(cutree(hcl, hcl_k_sil), DATA)
phm_hsil <- PHM(paramsList=hcl_sil_params, data=DATA, partition=cutree(hcl, hcl_k_sil))

## HDBSCAN
hdbcl <- hdbscan(DATA, 5)
plot(hdbcl, show_flat=T)

## Leiden
hgdp_seurat <- prepare_seurat(DATA)
hgdp_seurat <- FindNeighbors(hgdp_seurat, dims = 1:ncol(DATA))
hgdp_leiden <- FindClusters(hgdp_seurat, algorithm = 4) ## Leiden

leiden_labels <- Idents(hgdp_leiden)
leiden_K <- length(unique(leiden_labels))

leiden_params <- constructPmcParamsWeightedPartition(
    leiden_labels,
    DATA
)

phm_leiden <- PHM(paramsList = leiden_params,
                    data=DATA,
                    batchSize=1e5)

```

Visualize the results

```
## Get everything on the same scale
matrix_range <- range(c(3.742256e-06, 6.311475e-02, 6.335561e-06, 3.553335e-02))
matrix_range[1] <- matrix_range[1] * 0.9
matrix_range[2] <- matrix_range[2] * 1.1

## GMM
colors <- RColorBrewer::brewer.pal(mcl$G, "Paired")
dendro <- plotPHMDendrogram(phm_mcl, scaleHeights = "log10", colors=colors, suppressLabels = T)
distruct <- plotPHMDistruct(phm_mcl, labels=geographic_labels, colors=colors)
phm_matrix <- plotPHMMatrix(phm_mcl, colors=colors,
                            fillScale = "Pmc",
                            fillLimits = matrix_range,
                            legendPosition = "right")

hgdp_combined <- gridExtra::arrangeGrob(
    dendro + plot_theme + ggtitle("(a)"),
    phm_matrix + plot_theme + ggtitle("(b)"),
    distruct + plot_theme + ggtitle("(c)"),
    layout_matrix = matrix(c(1, 1, 3, 2, 2, 3), ncol=2)
)

## DensityCut
colors <- RColorBrewer::brewer.pal(length(dc_params), "Paired")
dc_distruct <- plotPHMDistruct(phm_dc, 
                                labels=geographic_labels, 
                                colors=colors, 
                                partition=T)
dc_distruct2 <- plotPHMDistruct(phm_dc, 
                                labels=geographic_labels, 
                                colors=colors)
dc_dendro <- plotPHMDendrogram(phm_dc, colors=colors,
                                scaleHeights = "log10",
                                suppressLabels = T)
dc_phm_matrix <- plotPHMMatrix(phm_dc, colors=colors,
                            fillScale = "Pmc",
                            legendPosition = "right")

  dc_hgdp_combined <- gridExtra::arrangeGrob(
    dc_dendro + plot_theme + ggtitle("(a)"),
    dc_phm_matrix + plot_theme + ggtitle("(b)"),
    dc_distruct + plot_theme + ggtitle("(c)"),
    dc_distruct2 + plot_theme + ggtitle("(d)"),
    layout_matrix = matrix(c(1, 1, 3, 4, 2, 2, 3, 4), ncol=2)
  )

## k-means
colors <- RColorBrewer::brewer.pal(length(kcl_sil_params), "Paired")
ksil_distruct <- plotPHMDistruct(phm_ksil,
                                labels=geographic_labels,
                                colors=colors,
                                partition=T)
ksil_distruct2 <- plotPHMDistruct(phm_ksil,
                                labels=geographic_labels,
                                colors=colors)
ksil_dendro <- plotPHMDendrogram(phm_ksil,
                                scaleHeights = "log10", 
                                colors=colors,
                                suppressLabels = T)
ksil_phm_matrix <- plotPHMMatrix(phm_ksil,
                                colors=colors,
                                fillScale = "Pmc",
                                legendPosition = "right")

ksil_hgdp_combined <- gridExtra::arrangeGrob(
    ksil_dendro + plot_theme + ggtitle("(a)"),
    ksil_phm_matrix + plot_theme + ggtitle("(b)"),
    ksil_distruct + plot_theme + ggtitle("(c)"),
    ksil_distruct2 + plot_theme + ggtitle("(d)"),
    layout_matrix = matrix(c(1, 1, 3, 4, 2, 2, 3, 4), ncol=2)
  )

## Hierarchical Clustering
colors <- RColorBrewer::brewer.pal(length(hcl_sil_params), "Paired")
hsil_distruct <- plotPHMDistruct(phm_hsil, 
                                labels=geographic_labels, 
                                colors=colors,
                                partition=T)
hsil_distruct2 <- plotPHMDistruct(phm_hsil, 
                                labels=geographic_labels, 
                                colors=colors)
hsil_dendro <- plotPHMDendrogram(phm_hsil,
                                scaleHeights = "log10", 
                                colors=colors,
                                suppressLabels = T)
hsil_phm_matrix <- plotPHMMatrix(phm_hsil,
                                colors=colors,
                                fillScale = "Pmc",
                                legendPosition = "right")

hsil_hgdp_combined <- gridExtra::arrangeGrob(
    hsil_dendro + plot_theme + ggtitle("(a)"),
    hsil_phm_matrix + plot_theme + ggtitle("(b)"),
    hsil_distruct + plot_theme + ggtitle("(c)"),
    hsil_distruct2 + plot_theme + ggtitle("(c)"),
    layout_matrix = matrix(c(1, 1, 3, 4, 2, 2, 3, 4), ncol=2)
  )

## Leiden
colors <- RColorBrewer::brewer.pal(length(leiden_params), "Paired")
leiden_distruct <- plotPHMDistruct(phm_leiden, 
                                labels=geographic_labels, 
                                colors=colors,
                                partition = T)
leiden_distruct2 <- plotPHMDistruct(phm_leiden, 
                                    labels=geographic_labels, 
                                    colors=colors)
leiden_dendro <- plotPHMDendrogram(phm_leiden,
                                scaleHeights = "log10", 
                                colors=colors,
                                suppressLabels = T)
leiden_phm_matrix <- plotPHMMatrix(phm_leiden,
                                colors=colors,
                                fillScale = "Pmc",
                                legendPosition = "right")
leiden_hgdp_combined <- gridExtra::arrangeGrob(
leiden_dendro + plot_theme + ggtitle("(a)"),
leiden_phm_matrix + plot_theme + ggtitle("(b)"),
leiden_distruct + plot_theme + ggtitle("(c)"),
leiden_distruct2 + plot_theme + ggtitle("(d)"),
layout_matrix = matrix(c(1, 1, 3, 4, 2, 2, 3, 4), ncol=2)
)
```

GMM Results

DensityCut results

Hierarchical Clustering Results

k-means results

Leiden results

Helper Functions

```
prepare_seurat <- function(data) {
  colnames(data) <- paste0("X", 1:ncol(data))
  rownames(data) <- paste0("R", 1:nrow(data))  
  seurat_obj <- SeuratObject::CreateSeuratObject(round(t(data)), assay="doot")
  seurat_embed <- SeuratObject::CreateDimReducObject(
    embeddings = data,
    key = "PC_",
    assay = "doot"
  )
  seurat_obj[["pca"]] <- seurat_embed
  
  seurat_obj
}
```