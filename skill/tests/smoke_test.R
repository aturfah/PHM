# Smoke test for the claims and snippets in phm/SKILL.md and phm/references/*.md
#
# Run:  Rscript tests/smoke_test.R
# Needs the PHM package and its dependencies installed
# (remotes::install_github("aturfah/PHM")) and jsonlite + python3 for the Explorer checks.
# Written against PHM 1.1.0 at GitHub main after the ensemble pipeline (Nov 2026).
#
# Each check() prints PASS/FAIL; the script exits non-zero if any check fails.

suppressPackageStartupMessages({
  library(PHM)
  library(mclust)
})

failures <- 0
check <- function(desc, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    cat("   error:", conditionMessage(e), "\n"); FALSE
  })
  cat(if (ok) "PASS " else "FAIL ", desc, "\n", sep = "")
  if (!ok) failures <<- failures + 1
}
raises <- function(expr) inherits(try(expr, silent = TRUE), "try-error")

## ---- toy data: 4 clusters, two of them close together -----------------------
set.seed(1)
N <- 150
X <- rbind(
  matrix(rnorm(N * 2, c(0, 0)),  ncol = 2, byrow = TRUE),
  matrix(rnorm(N * 2, c(2.5, 0)), ncol = 2, byrow = TRUE),
  matrix(rnorm(N * 2, c(0, 8)),  ncol = 2, byrow = TRUE),
  matrix(rnorm(N * 2, c(9, 8)),  ncol = 2, byrow = TRUE)
)
truth <- rep(1:4, each = N)

cat("\n== Entry points ==\n")
params <- constructPmcParamsWeightedPartition(truth, X, G = 1:2)
check("params list has prob/mean/var/class per cluster",
      length(params) == 4 && all(sapply(params, function(p) all(c("prob", "mean", "var", "class") %in% names(p)))))
check("cluster probs sum to 1", abs(sum(sapply(params, function(p) sum(p$prob))) - 1) < 1e-8)

set.seed(2)
phm <- PHM(paramsList = params, data = X, partition = truth,
           mcSamples = 2e5, batchSize = 1e5, verbose = FALSE)
check("phm is indexed by K (length == starting K)", length(phm) == 4)
check("phm[[K]]$clusters == K", all(sapply(seq_along(phm), function(k) phm[[k]]$clusters == k)))
check("level-K labels reproduce the input partition", all(as.character(phm[[4]]$labels) == as.character(truth)))
check("merged labels joined with '|'", any(grepl("|", phm[[2]]$labels, fixed = TRUE)))
check("posterior matrix is N x K", all(dim(phm[[3]]$posterior_matrix) == c(nrow(X), 3)))

cat("\n== P_mc bookkeeping ==\n")
pm <- phm[[4]]$pmc_matrix
check("pmc_matrix is symmetric with zero diagonal", isSymmetric(pm) && all(diag(pm) == 0))
check("phm[[K]]$pmc == sum(pmc_matrix)", abs(phm[[4]]$pmc - sum(pm)) < 1e-12)
pmc_mc <- computeMonteCarloPmc(params, mcSamples = 2e5)
check("sum(deltaPmc) matches computeMonteCarloPmc (no extra factor of 2)",
      abs(sum(pm) - pmc_mc) < 0.02)
check("pmc decreases monotonically to 0", {
  v <- sapply(phm, function(x) x$pmc); all(diff(v) >= -1e-12) && v[1] == 0
})
m <- phm[[4]]$merge_components
check("pmc_change at K=3 equals 2 * deltaPmc[i,j] of the merged pair",
      abs(phm[[3]]$pmc_change - 2 * pm[m[1], m[2]]) < 1e-12)
check("the two overlapping clusters (1,2) are merged first", all(sort(m) == c(1, 2)))

cat("\n== Extracting a solution ==\n")
sol <- thresholdPHM(phm, threshold = 0.01)
check("thresholdPHM is exported and returns a level", is.list(sol) && !is.null(sol$clusters))
check("thresholdPHM picks the finest K with remaining pmc < threshold", {
  v <- sapply(phm, function(x) x$pmc); sol$clusters == max(which(v < 0.01))
})
check("recoverPHMParams(phm, 2) returns 2 clusters", length(recoverPHMParams(phm, 2)) == 2)
check("recoverDeltaPmcMatrix(phm, 2) is 2 x 2", all(dim(recoverDeltaPmcMatrix(phm, 2)) == c(2, 2)))
check("recover* refuses K == length(phm)", raises(recoverPHMParams(phm, 4)))
# Known upstream bug (PHM 1.1.0): recoverPosterior() loops min_idx:(K+1), which runs *upward*
# when the posterior is already stored at K, and crashes. If this prints FIXED, update SKILL.md pitfall 6.
if (raises(recoverPosterior(phm, 2))) {
  cat("KNOWN BUG still present: recoverPosterior(phm, K) errors on default output\n")
} else {
  cat("FIXED upstream? recoverPosterior now works: update phm/SKILL.md pitfall 6\n")
}
check("phm[[K]]$posterior_matrix is the supported way to get posteriors",
      all(dim(phm[[2]]$posterior_matrix) == c(nrow(X), 2)))
phm_np <- PHM(paramsList = params, data = X, partition = truth, computePosterior = FALSE,
              mcSamples = 1e5, verbose = FALSE)
phm_ap <- addPosteriorMatrix(phm_np, X)
check("addPosteriorMatrix fills posterior_matrix and labels at every level",
      all(dim(phm_ap[[2]]$posterior_matrix) == c(nrow(X), 2)) && length(phm_ap[[2]]$labels) == nrow(X))

cat("\n== PHM(deltaPmc=) is ignored ==\n")
bogus <- matrix(0, 4, 4)
set.seed(2)
phm_b <- PHM(paramsList = params, data = X, partition = truth, deltaPmc = bogus,
             mcSamples = 2e5, batchSize = 1e5, verbose = FALSE)
check("PHM() ignores a supplied deltaPmc", sum(phm_b[[4]]$pmc_matrix) > 0)

cat("\n== Shorthand + mclustObj entry points ==\n")
set.seed(3)
phm_short <- PHM(partition = truth, data = X, partitionMaxComponents = 2,
                 mcSamples = 1e5, verbose = FALSE)
check("PHM(partition=, data=) shorthand runs", length(phm_short) == 4)
mcl <- Mclust(X, G = 4, verbose = FALSE)
set.seed(3)
phm_mcl <- PHM(mclustObj = mcl, mcSamples = 1e5, batchSize = 1e5, verbose = FALSE)
check("PHM(mclustObj=) runs and uses G components", length(phm_mcl) == mcl$G)

cat("\n== Data validation ==\n")
check("data.frame is rejected by weighted partition",
      raises(constructPmcParamsWeightedPartition(truth, as.data.frame(X))))

cat("\n== Tiny clusters ==\n")
X_tiny   <- rbind(X[1:100, ], X[N + 1:100, ], X[3 * N + 1:3, ])
lab_tiny <- c(rep(1, 100), rep(2, 100), rep(3, 3))
check("weighted partition errors on a 3-point cluster even with modelNames='VVI', G=1",
      raises(constructPmcParamsWeightedPartition(lab_tiny, X_tiny, modelNames = "VVI", G = 1)))
check("naive partition with VVI, G=1 tolerates a 3-point cluster",
      !raises(constructPmcParamsPartition(lab_tiny, X_tiny, modelNames = "VVI", G = 1)))

cat("\n== Weights column order ==\n")
# Relabel so that label 2 is the cluster centred at (0,0) and label 1 the one at (2.5,0),
# and move a label-2 row first so that unique(partition) = c(2, 1, 3, 4) (not sorted).
lab_swapped <- ifelse(truth == 1, 2, ifelse(truth == 2, 1, truth))
first_rows  <- c(1, N + 1, 2 * N + 1, 3 * N + 1)
row_order   <- c(first_rows, setdiff(seq_along(lab_swapped), first_rows))
lab_s <- lab_swapped[row_order]
X_s   <- X[row_order, ]
first_seen <- as.character(unique(factor(lab_s)))
check("unique(partition) is order-of-appearance (2,1,3,4), not sorted",
      identical(first_seen, c("2", "1", "3", "4")))
W_seen   <- sapply(first_seen,           function(k) as.numeric(lab_s == k))
W_sorted <- sapply(sort(unique(lab_s)),  function(k) as.numeric(lab_s == k))
p_seen   <- constructPmcParamsWeightedPartition(lab_s, X_s, weights = W_seen,   G = 1)
p_sorted <- constructPmcParamsWeightedPartition(lab_s, X_s, weights = W_sorted, G = 1)
mean_x <- function(p, k) p[[which(sapply(p, function(z) z$class) == k)]]$mean[1, 1]
check("first-appearance column order: label 1 -> x ~ 2.5 and label 2 -> x ~ 0",
      abs(mean_x(p_seen, 1) - 2.5) < 0.5 && abs(mean_x(p_seen, 2) - 0) < 0.5)
check("sorted column order silently swaps the two clusters",
      abs(mean_x(p_sorted, 1) - 0) < 0.5 && abs(mean_x(p_sorted, 2) - 2.5) < 0.5)

cat("\n== Nothing beyond PHM needs attaching ==\n")
run_r <- function(code) system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)), stdout = TRUE, stderr = FALSE)
out <- run_r(paste0(
  "suppressPackageStartupMessages(library(PHM)); set.seed(1); ",
  "x <- matrix(rnorm(400), ncol = 2); l <- rep(1:2, each = 100); ",
  "r <- try(constructPmcParamsWeightedPartition(l, x, G = 1), silent = TRUE); cat(!inherits(r, 'try-error'))"))
check("weighted partition works with only PHM attached (mclust is imported)", identical(tail(out, 1), "TRUE"))
out <- run_r(paste0(
  "suppressPackageStartupMessages(library(PHM)); set.seed(1); ",
  "x <- matrix(rnorm(400), ncol = 2); r <- try(mclust::Mclust(x, G = 1:2, verbose = FALSE), silent = TRUE); ",
  "cat(inherits(r, 'try-error'))"))
check("user code calling mclust::Mclust() without library(mclust) fails (needs mclustBIC in scope)",
      identical(tail(out, 1), "TRUE"))
check("plotPmcMatrix works without attaching ggplot2", !raises(plotPmcMatrix(phm, K = 3)))

cat("\n== Visualization ==\n")
check("plotPHMDendrogram returns ggplot", inherits(plotPHMDendrogram(phm, displayAxis = "label"), "ggplot"))
check("plotPHMMatrix returns ggplot", inherits(plotPHMMatrix(phm, fillScale = "pmcdist"), "ggplot"))
check("plotPmcMatrix returns ggplot", inherits(plotPmcMatrix(phm, K = 3), "ggplot"))
check("plotPHMDistruct returns ggplot",
      inherits(plotPHMDistruct(phm, K = 4, labels = factor(truth)), "ggplot"))
cols4 <- RColorBrewer::brewer.pal(4, "Paired")
check("plotPmc2D with its default PmcColor '#222' errors (KNOWN BUG; if this FAILs, update SKILL.md)",
      raises(plotPmc2D(params, X, truth, colors = cols4, numPmcPatches = 40)))
check("plotPmc2D works with PmcColor = '#222222'",
      !raises(plotPmc2D(params, X, truth, colors = cols4, numPmcPatches = 40, PmcColor = "#222222")))
check("plotPHMDistruct default palette fails for K = 10 (Set1 has 9 colors)", {
  # build a fake phm with K=10 by reusing structure is heavy; test the palette claim directly
  suppressWarnings(length(RColorBrewer::brewer.pal(10, "Set1"))) == 9
})

cat("\n== PHMv2 ==\n")
set.seed(4)
phm2 <- PHMv2(paramsList = params, mcSamples = 2e5, batchSize = 1e5)
check("PHMv2 default scaling is 'alpha'", phm2$mergeCriterion == "alpha")
set.seed(4)
phm2u <- PHMv2(paramsList = params, scaling = "unscaled", mcSamples = 2e5, batchSize = 1e5)
check("PHMv2 unscaled merges the same first pair as PHM()",
      all(sort(unname(phm2u$mergeComps[[4]])) == sort(phm[[4]]$merge_components)))
check("alpha scaling requires paramsList",
      raises(PHMv2(deltaPmc = phm2$deltaPmc, scaling = "alpha")))
check("PHMv2 accepts a precomputed deltaPmc", {
  r <- PHMv2(deltaPmc = phm2$deltaPmc, scaling = "unscaled"); r$mergeCriterion == "unscaled"
})
groups <- recoverGroupsv2(phm2, k = 2)
check("recoverGroupsv2 returns 2 groups covering all 4 components", length(groups) == 2 && setequal(unlist(groups), 1:4))
post <- recoverPosteriorv2(phm2, 2, computePosteriorProbMatrix(params, X))
check("recoverPosteriorv2 gives N x 2 rows summing to 1",
      all(dim(post) == c(nrow(X), 2)) && all(abs(rowSums(post) - 1) < 1e-8))
check("recoverParamsv2 returns 2 clusters", length(recoverParamsv2(phm2, 2)) == 2)
check("recoverDeltaPmcv2 is 2 x 2", all(dim(recoverDeltaPmcv2(phm2, 2)) == c(2, 2)))
check("plotPHMv2Dendrogram (default d*) works with alpha", inherits(plotPHMv2Dendrogram(phm2), "ggplot"))
check("plotPHMv2Dendrogram 'pmcdist' errors with alpha", raises(plotPHMv2Dendrogram(phm2, scaleHeights = "pmcdist")))
check("plotPHMv2Dendrogram 'd*' errors without alpha", raises(plotPHMv2Dendrogram(phm2u, scaleHeights = "d*")))
check("plotPHMv2Dendrogram log10 works for unscaled", inherits(plotPHMv2Dendrogram(phm2u, scaleHeights = "log10"), "ggplot"))
check("plotPHMv2Heatmap works, with subsetClusters",
      inherits(plotPHMv2Heatmap(phm2, subsetClusters = groups[[1]]), "ggplot"))
conv <- convertToPHMv2(phm)
check("convertToPHMv2 result plots with v2 dendrogram",
      inherits(plotPHMv2Dendrogram(conv, scaleHeights = "log10"), "ggplot"))

cat("\n== Utilities ==\n")
check("decomposeParams splits mixtures into single components",
      length(decomposeParams(params)) == sum(sapply(params, function(p) length(p$prob))))
sub <- constructPmcParamsSubAggMclust(X, replicates = 2, G = 1:2)
check("constructPmcParamsSubAggMclust (arg 'replicates') runs", is.list(sub) && length(sub) >= 2)
check("constructPmcParamsSubAggMclust takes 'subsampSize'",
      "subsampSize" %in% names(formals(constructPmcParamsSubAggMclust)))
check("bundled datasets load", {
  for (d in c("crater", "double_ring", "moons", "nested_diamonds")) data(list = d, package = "PHM")
  is.matrix(double_ring) && ncol(double_ring) == 2
})

cat("\n== Workflow recipes (workflows.md) on the bundled moons data ==\n")
data("moons", package = "PHM")
Xm <- moons
set.seed(10)

## 1. GMM -> hierarchy
mclm <- Mclust(Xm, G = 1:15, verbose = FALSE)
phm_g <- PHM(mclustObj = mclm, mcSamples = 2e5, batchSize = 1e5, verbose = FALSE)
check("recipe 1: GMM hierarchy has length G and monotone pmc",
      length(phm_g) == mclm$G && all(diff(sapply(phm_g, function(x) x$pmc)) >= -1e-12))
check("recipe 1: dendrogram + heatmap with axis 'none' work for many components", {
  inherits(plotPHMDendrogram(phm_g, suppressLabels = TRUE, displayAxis = "none"), "ggplot") &&
    inherits(plotPHMMatrix(phm_g, displayAxis = "none", fillScale = "pmcdist"), "ggplot")
})

## 2. any clustering -> hierarchy (integer, character and factor labels)
hcl <- hclust(dist(Xm), method = "ward.D2")
labs <- cutree(hcl, k = 8)
params_m <- constructPmcParamsWeightedPartition(labs, Xm, G = 1:2, modelNames = c("VVI", "VVV"))
phm_h <- PHM(paramsList = params_m, data = Xm, partition = labs, mcSamples = 2e5, batchSize = 1e5, verbose = FALSE)
check("recipe 2: hclust partition -> PHM with K = 8", length(phm_h) == 8)
lab_chr <- paste0("c", labs)
p_chr <- constructPmcParamsWeightedPartition(lab_chr, Xm, G = 1)
phm_chr <- PHM(paramsList = p_chr, data = Xm, partition = lab_chr, mcSamples = 1e5, verbose = FALSE)
check("recipe 2: character labels work", length(phm_chr) == 8 && any(grepl("|", phm_chr[[2]]$labels, fixed = TRUE)))
lab_fac <- factor(labs)
p_fac <- constructPmcParamsWeightedPartition(lab_fac, Xm, G = 1)
phm_fac <- PHM(paramsList = p_fac, data = Xm, partition = lab_fac, mcSamples = 1e5, verbose = FALSE)
check("recipe 2: factor labels work", length(phm_fac) == 8)

## 3. compare clusterers with a common fill scale
cands <- list(gmm = mclm$classification, hclust = labs, kmeans = kmeans(Xm, 8, nstart = 5)$cluster)
phms <- lapply(cands, function(lab) {
  p <- constructPmcParamsWeightedPartition(lab, Xm, G = 1)
  PHM(paramsList = p, data = Xm, partition = lab, mcSamples = 1e5, batchSize = 1e5, verbose = FALSE)
})
rng <- range(na.omit(unlist(lapply(phms, function(p) sapply(p, `[[`, "pmc_change")))))
check("recipe 3: heatmaps accept shared fillLimits on the raw pmc_change scale",
      all(sapply(phms, function(p) inherits(plotPHMMatrix(p, displayAxis = "label", fillScale = "pmcdist", fillLimits = rng), "ggplot"))))

## 4. soft memberships as weights
labz <- mclm$classification
Wz <- mclm$z[, unique(labz), drop = FALSE]
pz <- constructPmcParamsWeightedPartition(labz, Xm, weights = Wz, G = 1:2)
phm_z <- PHM(paramsList = pz, data = Xm, partition = labz, mcSamples = 1e5, batchSize = 1e5, verbose = FALSE)
check("recipe 4: soft-membership weights (columns in unique(labels) order) run", length(phm_z) == length(unique(labz)))

## 5. choose a resolution, extract labels
path <- data.frame(K = seq_along(phm_h), pmc = sapply(phm_h, `[[`, "pmc"), drop = sapply(phm_h, `[[`, "pmc_change"))
sol <- thresholdPHM(phm_h, threshold = 0.01)
labels_sol <- sol$labels
truth_m <- labs
check("recipe 5: path data.frame, thresholdPHM labels, ARI and table run", {
  nrow(path) == 8 && length(labels_sol) == nrow(Xm) &&
    is.numeric(adjustedRandIndex(truth_m, labels_sol)) && is.table(table(truth_m, labels_sol))
})
check("recipe 5: strsplit on '|' recovers original labels",
      all(unlist(strsplit(as.character(phm_h[[2]]$labels[1]), "|", fixed = TRUE)) %in% as.character(1:8)))
check("recipe 5: plotPHMDistruct with Set1 colors", {
  Ks <- sol$clusters
  inherits(plotPHMDistruct(phm_h, K = Ks, labels = factor(truth_m), colors = RColorBrewer::brewer.pal(max(Ks, 3), "Set1")), "ggplot")
})

## 5b. exact unions of the input clusters
union_labels <- function(phm, params, orig_labels, k) {
  grp <- recoverGroupsv2(convertToPHMv2(phm), k = k)
  cls <- sapply(params, `[[`, "class")
  out <- integer(length(orig_labels))
  for (g in seq_along(grp)) out[orig_labels %in% cls[grp[[g]]]] <- g
  out
}
split_count <- function(u) sum(sapply(unique(labs), function(o) length(unique(u[labs == o])) > 1))
check("recipe 5: union_labels gives exact unions (no input cluster split) at K = 7, 4, 2",
      all(sapply(c(7, 4, 2), function(k) {
        u <- union_labels(phm_h, params_m, labs, k); length(unique(u)) == k && split_count(u) == 0
      })))
lk7 <- as.character(phm_h[[7]]$labels)
n_split7 <- sum(sapply(unique(labs), function(o) length(unique(lk7[labs == o])) > 1))
cat("INFO pitfall 5: at K = 7, ", n_split7, " of 8 input clusters have points relabeled by the posterior argmax",
    " (SKILL.md says 3 to 4; it varies with the Monte Carlo draw)\n", sep = "")

## 6. PHMv2 alpha, zoom, initK
phm2m <- PHMv2(paramsList = params_m, scaling = "alpha", mcSamples = 2e5, batchSize = 1e5, numCores = 1)
grp <- recoverGroupsv2(phm2m, k = 4)
postm <- recoverPosteriorv2(phm2m, 4, computePosteriorProbMatrix(params_m, Xm))
check("recipe 6: v2 groups/posterior give 4 clusters", length(grp) == 4 && ncol(postm) == 4)
check("recipe 6: v2 heatmap subset + initK dendrogram + d* heatmap", {
  inherits(plotPHMv2Heatmap(phm2m, subsetClusters = grp[[1]], displayAxis = "label"), "ggplot") &&
    inherits(plotPHMv2Dendrogram(phm2m, initK = 6), "ggplot") &&
    inherits(plotPHMv2Heatmap(phm2m, fillScale = "d*", displayAxis = "label"), "ggplot")
})

## 7. consensus over subsamples
cache <- tempfile("subagg_cache")
p_sub <- constructPmcParamsSubAggMclust(Xm, replicates = 3, subsampSize = 300, G = 1:3, saveDir = cache)
phm_sub <- PHM(paramsList = p_sub, data = Xm, computePosterior = FALSE, mcSamples = 1e5, batchSize = 1e5, verbose = FALSE)
check("recipe 7: SubAggMclust -> PHM(computePosterior = FALSE) -> dendrogram without labels",
      length(phm_sub) == length(p_sub) && inherits(plotPHMDendrogram(phm_sub, displayAxis = "none", suppressLabels = TRUE), "ggplot"))
check("recipe 7: saveDir caches replicate fits", length(list.files(cache)) >= 3)
check("constructPmcParamsSubAggPartition is disabled (stops), as api.md says",
      raises(constructPmcParamsSubAggPartition(Xm, labs, replicates = 2, subsampSize = 150, G = 1)))

## 10. 2-D diagnostic
check("recipe 10: plotPmc2D with 8 clusters",
      !raises(plotPmc2D(PmcColor = "#222222", params_m, Xm, labs, colors = RColorBrewer::brewer.pal(length(params_m), "Paired"),
                         numPmcPatches = 60, densityLevels = c(1e-2, 1e-1))))

cat("\n== PHMv2 as the default workflow (SKILL.md core recipe, workflows.md 1, 3, 5, 6, 7) ==\n")
remaining_pmc <- function(p2) cumsum(p2$mergeDeltaPmc)

## remaining P_mc from the merge history, for every scaling
p2s <- lapply(c("alpha", "unscaled", "average"), function(sc) {
  set.seed(1); PHMv2(paramsList = params_m, scaling = sc, mcSamples = 2e5, batchSize = 1e5)
})
names(p2s) <- c("alpha", "unscaled", "average")
check("cumsum(mergeDeltaPmc) equals sum(recoverDeltaPmcv2) at every K, for all three scalings",
      all(sapply(p2s, function(p2) {
        direct <- sapply(1:8, function(k) if (k == 1) 0 else sum(recoverDeltaPmcv2(p2, k)))
        max(abs(remaining_pmc(p2) - direct)) < 1e-10
      })))
check("mergeDeltaPmc[1] is a placeholder zero and remaining P_mc is monotone",
      all(sapply(p2s, function(p2) p2$mergeDeltaPmc[1] == 0 && all(diff(remaining_pmc(p2)) >= -1e-12))))
set.seed(1)
p1_ref <- PHM(paramsList = params_m, data = Xm, partition = labs, mcSamples = 2e5, batchSize = 1e5, verbose = FALSE)
p2_ref <- PHMv2(deltaPmc = p1_ref[[8]]$pmc_matrix, scaling = "unscaled")
check("PHMv2(deltaPmc = <v1 matrix>, 'unscaled') reproduces PHM() merges and remaining P_mc exactly",
      all(sapply(2:8, function(k) all(unname(p2_ref$mergeComps[[k]]) == p1_ref[[k]]$merge_components))) &&
        max(abs(remaining_pmc(p2_ref) - sapply(p1_ref, `[[`, "pmc"))) < 1e-10)
K_star <- max(which(remaining_pmc(p2s$unscaled) < 0.01))
check("K_star from the threshold is a valid K", K_star >= 1 && K_star <= 8)

## mclust route
p_mcl <- constructPmcParamsMclust(mclm)
set.seed(3); phm2_g <- PHMv2(paramsList = p_mcl, mcSamples = 2e5, batchSize = 1e5)
check("recipe 1 (v2): Mclust -> constructPmcParamsMclust -> PHMv2", ncol(phm2_g$deltaPmc) == mclm$G)
check("recipe 1 (v2): v2 dendrogram and heatmap with axis 'label' and 'none'",
      all(sapply(c("label", "none"), function(ax)
        inherits(plotPHMv2Dendrogram(phm2_g, displayAxis = ax), "ggplot") &&
          inherits(plotPHMv2Heatmap(phm2_g, displayAxis = ax), "ggplot"))))

## exact-union labels straight from a v2 object, and the soft alternative
union_labels_v2 <- function(phm2, params, orig_labels, k) {
  grp <- recoverGroupsv2(phm2, k)
  cls <- sapply(params, `[[`, "class")
  out <- integer(length(orig_labels))
  for (g in seq_along(grp)) out[orig_labels %in% cls[grp[[g]]]] <- g
  out
}
p2a <- p2s$alpha
check("recipe 5 (v2): union_labels on a v2 object gives exact unions at K = 7, 4, 2",
      all(sapply(c(7, 4, 2), function(k) {
        u <- union_labels_v2(p2a, params_m, labs, k); length(unique(u)) == k && split_count(u) == 0
      })))
u4 <- union_labels_v2(p2a, params_m, labs, 4)
post4 <- recoverPosteriorv2(p2a, 4, computePosteriorProbMatrix(params_m, Xm))
lab_soft4 <- apply(post4, 1, which.max)
agree <- mean(ave(u4, lab_soft4, FUN = function(z) as.integer(names(sort(table(z), decreasing = TRUE))[1])) == u4)
check("recipe 5 (v2): posterior argmax labels mostly agree with the union labels (cluster order matches recoverGroupsv2)",
      all(dim(post4) == c(nrow(Xm), 4)) && agree > 0.9 && mean(lab_soft4 == u4) > 0.9)
gmm_union <- union_labels_v2(phm2_g, p_mcl, mclm$classification, 3)
check("recipe 5 (v2): union_labels for a GMM with mcl$classification", length(unique(gmm_union)) == 3)

## reuse of the delta matrix and merge criteria
D <- p2s$alpha$deltaPmc
r_raw <- PHMv2(paramsList = params_m, deltaPmc = D, scaling = "unscaled")
r_avg <- PHMv2(paramsList = params_m, deltaPmc = D, scaling = "average")
check("recipe 6: paramsList + deltaPmc skips Monte Carlo (deltaPmc unchanged) and honours the criterion",
      identical(r_raw$deltaPmc, D) && r_raw$mergeCriterion == "unscaled" && r_avg$mergeCriterion == "average")
check("recipe 6: deltaPmc alone works for 'average' but not 'alpha'",
      !raises(PHMv2(deltaPmc = D, scaling = "average")) && raises(PHMv2(deltaPmc = D, scaling = "alpha")))
check("recipe 6: d* for alpha, log10/pmcdist/unscaled for non-alpha plots",
      inherits(plotPHMv2Dendrogram(p2a, scaleHeights = "d*"), "ggplot") &&
        inherits(plotPHMv2Dendrogram(r_raw, scaleHeights = "log10"), "ggplot") &&
        inherits(plotPHMv2Dendrogram(r_raw, scaleHeights = "unscaled"), "ggplot") &&
        inherits(plotPHMv2Heatmap(r_raw, fillScale = "pmcdist"), "ggplot"))
check("recipe 6: alpha dendrogram accepts 'log10' and 'unscaled' too",
      inherits(plotPHMv2Dendrogram(p2a, scaleHeights = "log10"), "ggplot") &&
        inherits(plotPHMv2Dendrogram(p2a, scaleHeights = "unscaled"), "ggplot"))

## comparing runs (v2 fillLimits on the raw mergeValues scale)
phms2 <- lapply(cands, function(lab) {
  PHMv2(paramsList = constructPmcParamsWeightedPartition(lab, Xm, G = 1), mcSamples = 1e5, batchSize = 1e5)
})
rng2 <- range(unlist(lapply(phms2, function(p) p$mergeValues[-1])))
check("recipe 3 (v2): heatmaps share fillLimits = range(mergeValues[-1]); dendrograms accept ggtitle",
      all(sapply(names(phms2), function(nm)
        inherits(plotPHMv2Heatmap(phms2[[nm]], displayAxis = "label", fillLimits = rng2), "ggplot") &&
          inherits(plotPHMv2Dendrogram(phms2[[nm]], displayAxis = "label") + ggplot2::ggtitle(nm), "ggplot"))))

## consensus with v2
phm2_sub <- PHMv2(paramsList = p_sub, mcSamples = 1e5, batchSize = 1e5, numCores = 1)
check("recipe 7 (v2): SubAggMclust -> PHMv2 -> dendrogram with axis 'none'",
      ncol(phm2_sub$deltaPmc) == length(p_sub) && inherits(plotPHMv2Dendrogram(phm2_sub, displayAxis = "none"), "ggplot"))

## v1 and v2 families do not mix
phm_v2_obj <- p2s$unscaled
check("v1 plot functions error on a PHMv2 result", raises(plotPHMDendrogram(phm_v2_obj)) && raises(plotPmcMatrix(phm_v2_obj)))
check("v2 plot functions error on a raw PHM() result", raises(plotPHMv2Dendrogram(phm_h)))

## what must be attached for the v2 path
out2 <- run_r(paste0(
  "suppressPackageStartupMessages({library(PHM); library(mclust)}); set.seed(1); ",
  "x <- rbind(matrix(rnorm(200), ncol = 2), matrix(rnorm(200, 4), ncol = 2)); l <- rep(1:2, each = 100); ",
  "p <- constructPmcParamsWeightedPartition(l, x, G = 1); m <- PHMv2(paramsList = p, mcSamples = 1e5); ",
  "a <- plotPHMv2Dendrogram(m); b <- plotPHMv2Heatmap(m); cat(inherits(a, 'ggplot') && inherits(b, 'ggplot'))"))
check("v2 pipeline and plots run with only PHM and mclust attached (no ggplot2/dplyr)", identical(tail(out2, 1), "TRUE"))


cat("\n== Ensemble and density estimation (SKILL.md entry points, ensemble.md) ==\n")
data("double_ring", package = "PHM")
set.seed(20)
Xr <- double_ring[sample(nrow(double_ring), 1200), ]
check("phmSeeds(seedBase, M) is set.seed(seedBase); sample.int(1e7, M)",
      identical(phmSeeds(1000, 4), { set.seed(1000); sample.int(1e7, 4) }))
rep1 <- fitDensityReplicate(Xr, idx = 1, seed = phmSeeds(1000, 3)[1], rho = 0.5, G = 1:6, extendVerbose = FALSE)
check("fitDensityReplicate returns params, loglik, bic, idx, G, modelName",
      all(c("params", "loglik", "bic", "idx", "G", "modelName") %in% names(rep1)))
rep1b <- fitDensityReplicate(Xr, idx = 1, seed = phmSeeds(1000, 3)[1], rho = 0.5, G = 1:6, extendVerbose = FALSE)
check("a replicate refit with the same seed is identical", identical(rep1, rep1b))
check("G is extended beyond the starting range when BIC wants more (G = 1:6 start)", rep1$G > 6)
ens <- fitPHMEnsemble(Xr, M = 3, rho = 0.5, G = 1:6, numCores = 1, mcSamples = 5e4, extendVerbose = FALSE)
check("fitPHMEnsemble returns a PHMv2 object with alpha scaling",
      all(c("deltaPmc", "paramsList", "mergeComps", "mergeValues", "mergeDeltaPmc", "mergeCriterion") %in% names(ens)) &&
        ens$mergeCriterion == "alpha")
check("pooled components are labelled '<component>_<replicate>', single Gaussians, probs sum to 1",
      all(grepl("^[0-9]+_[1-3]$", sapply(ens$paramsList, `[[`, "class"))) &&
        all(sapply(ens$paramsList, function(p) length(p$prob)) == 1) &&
        abs(sum(sapply(ens$paramsList, `[[`, "prob")) - 1) < 1e-8)
check("pooled order is replicate 1 first", all(grepl("_1$", sapply(ens$paramsList, `[[`, "class")[seq_len(rep1$G)])))
reps <- lapply(1:3, function(i) fitDensityReplicate(Xr, idx = i, seed = phmSeeds(1000, 3)[i], rho = 0.5, G = 1:6, extendVerbose = FALSE))
set.seed(5); c1 <- combineDensityReplicates(reps, numCores = 1, mcSamples = 5e4)
check("map + reduce by hand pools the same components as fitPHMEnsemble",
      isTRUE(all.equal(c1$paramsList, ens$paramsList)))
gde <- globalDensityEstimation(Xr, M = 3, rho = 0.5, G = 1:6, extendVerbose = FALSE)
check("globalDensityEstimation(X, M, rho) equals fitPHMEnsemble(X, M, rho)$paramsList", isTRUE(all.equal(gde, ens$paramsList)))
cache <- tempfile("gde_cache")
gde2 <- globalDensityEstimation(Xr, M = 2, rho = 0.5, G = 1:6, saveDir = cache, extendVerbose = FALSE)
check("globalDensityEstimation saveDir caches one file per replicate", length(list.files(cache)) == 2)
check("rho = 1, N <= 2000, no Stage 2: replicates are identical (pitfall 5)", {
  small <- Xr[1:600, ]
  a <- fitDensityReplicate(small, 1, phmSeeds(1000, 2)[1], rho = 1, G = 1:6, extendVerbose = FALSE)
  b <- fitDensityReplicate(small, 2, phmSeeds(1000, 2)[2], rho = 1, G = 1:6, extendVerbose = FALSE)
  isTRUE(all.equal(a$loglik, b$loglik)) && a$G == b$G
})
check("stage2Init = 'kmeans' makes those replicates differ", {
  small <- Xr[1:600, ]
  a <- fitDensityReplicate(small, 1, phmSeeds(1000, 2)[1], rho = 1, G = 1:6, stage2Init = "kmeans", extendVerbose = FALSE)
  b <- fitDensityReplicate(small, 2, phmSeeds(1000, 2)[2], rho = 1, G = 1:6, stage2Init = "kmeans", extendVerbose = FALSE)
  !isTRUE(all.equal(a$params, b$params))
})
K2 <- 2
post_e <- recoverPosteriorv2(ens, K2, computePosteriorProbMatrix(ens$paramsList, Xr))
lab_e <- apply(post_e, 1, which.max)
truth_r <- as.integer(sqrt(rowSums(Xr^2)) > 4.5) + 1
check("core recipe route 1: ensemble on double_ring separates the two rings at K = 2 (ARI > 0.95)",
      mclust::adjustedRandIndex(truth_r, lab_e) > 0.95)
grp_e <- recoverGroupsv2(ens, K2)
mass <- sapply(grp_e, function(g) sum(sapply(ens$paramsList[g], function(p) sum(p$prob))))
check("pitfall 6 snippet: group masses sum to 1", abs(sum(mass) - 1) < 1e-8)
check("ensemble dendrogram with initK and displayAxis = 'none'",
      inherits(plotPHMv2Dendrogram(ens, initK = 10, displayAxis = "none"), "ggplot"))
check("reduce step can be redone from paramsList + deltaPmc without Monte Carlo",
      PHMv2(paramsList = ens$paramsList, deltaPmc = ens$deltaPmc, scaling = "unscaled")$mergeCriterion == "unscaled")
check("ensemble command-line scripts are installed",
      all(file.exists(file.path(system.file("ensemble", package = "PHM"),
          c("run_serial.R", "config/config_template.R", "parallel_run/run_replicate.R", "parallel_run/run_combine.R",
            "parallel_run/single_machine/launch.sh", "parallel_run/slurm_array/submit.sh",
            "example/double_ring/config.R")))))

## localized density estimation (workflows.md 7b)
km <- function(X, K, seed = 1, ...) { set.seed(seed); kmeans(X, K, nstart = 5)$cluster }
loc <- localizedDensityEstimation(Xm, partitionFunc = km, partitionFuncParams = list(K = 6),
                                  weightFunc = distanceWeightFunc, preserveClusters = TRUE, G = 1:3)
check("localizedDensityEstimation preserveClusters = TRUE: one mixture per cluster, class '<rep>_<label>'",
      length(loc) == 6 && all(grepl("^1_[1-6]$", sapply(loc, `[[`, "class"))))
loc2 <- localizedDensityEstimation(Xm, partitionFunc = km, partitionFuncParams = list(K = 6),
                                   weightFunc = distanceWeightFunc, preserveClusters = FALSE, G = 1:3)
check("preserveClusters = FALSE: single Gaussians, probs sum to 1",
      all(sapply(loc2, function(p) length(p$prob)) == 1) && abs(sum(sapply(loc2, `[[`, "prob")) - 1) < 1e-8)
check("workflows 7b: PHMv2 on localized params runs", length(PHMv2(paramsList = loc, mcSamples = 5e4)$mergeComps) == 6)
km_noseed <- function(X, K) kmeans(X, K)$cluster
check("a partition function without a seed (or ...) argument fails (it is always passed seed)",
      raises(localizedDensityEstimation(Xm, partitionFunc = km_noseed, partitionFuncParams = list(K = 4),
                                        weightFunc = naiveWeightFunc)))
km_gap <- function(X, K, seed = 1, ...) { set.seed(seed); c(1, 3, 5, 7)[kmeans(X, K, nstart = 3)$cluster] }
check("labels that are not 1..K fail (weights looked up by label value)",
      raises(localizedDensityEstimation(Xm, partitionFunc = km_gap, partitionFuncParams = list(K = 4),
                                        weightFunc = distanceWeightFunc, G = 1)))
w <- distanceWeightFunc(Xm, km(Xm, 4))
check("distanceWeightFunc gives an N x K matrix with rows summing to 1", all(dim(w) == c(nrow(Xm), 4)) && all(abs(rowSums(w) - 1) < 1e-8))

cat("\n== PHM Explorer (explorer.md) ==\n")
explorer_dir <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                                        "..", "phm", "scripts", "explorer"), mustWork = FALSE)
check("bundled Explorer files exist",
      all(file.exists(file.path(explorer_dir, c("r/export_phm_json.R", "r/view_in_explorer.R", "src/template.html", "src/build.py")))))
if (requireNamespace("jsonlite", quietly = TRUE) && nzchar(Sys.which("python3"))) {
  source(file.path(explorer_dir, "r", "export_phm_json.R"))
  tdir <- tempfile("explorer"); dir.create(tdir)
  invisible(capture.output(export_phm_json(ens, X = Xr, out_path = file.path(tdir, "ens.json"), label = "test")))
  js <- jsonlite::fromJSON(file.path(tdir, "ens.json"), simplifyVector = FALSE)
  check("export of an ensemble fit has meta, nodes, xlabOrder and points",
        all(c("meta", "nodes", "xlabOrder", "points") %in% names(js)) && js$meta$K == length(ens$paramsList))
  params_mc <- constructPmcParamsWeightedPartition(labs, Xm, G = 2, modelNames = "VVV")   # two Gaussians per cluster
  fit_part <- PHMv2(paramsList = params_mc, mcSamples = 5e4)
  check("multi-component partition clusters: export with X fails (explorer.md)",
        raises(capture.output(export_phm_json(fit_part, X = Xm, out_path = file.path(tdir, "p.json")))))
  check("... and works with z = computePosteriorProbMatrix() and xy",
        !raises(capture.output(export_phm_json(fit_part, xy = Xm, z = computePosteriorProbMatrix(params_mc, Xm),
                                               out_path = file.path(tdir, "p.json")))))
  fit_un <- PHMv2(paramsList = params_mc, deltaPmc = fit_part$deltaPmc, scaling = "unscaled")
  check("non-alpha fits are rejected by the exporter",
        raises(capture.output(export_phm_json(fit_un, xy = Xm, z = computePosteriorProbMatrix(params_mc, Xm),
                                              out_path = file.path(tdir, "u.json")))))
  saveRDS(ens, file.path(tdir, "fit.rds")); saveRDS(Xr, file.path(tdir, "X.rds"))
  saveRDS(data.frame(x = Xr[, 1], y = Xr[, 2], annotation = factor(truth_r)), file.path(tdir, "plot.rds"))
  st <- system2(file.path(R.home("bin"), "Rscript"),
                c(shQuote(file.path(explorer_dir, "r", "view_in_explorer.R")), shQuote(file.path(tdir, "fit.rds")),
                  shQuote(file.path(tdir, "X.rds")), paste0("--plot=", shQuote(file.path(tdir, "plot.rds")))),
                stdout = FALSE, stderr = FALSE, env = "BROWSER=true")
  html <- file.path(tdir, "fit_explorer.html")
  check("view_in_explorer.R writes <fit>_export.json and <fit>_explorer.html with the data embedded",
        file.exists(file.path(tdir, "fit_export.json")) && file.exists(html) &&
          !any(grepl("__DATA_JSON__", readLines(html, warn = FALSE), fixed = TRUE)))
  saveRDS(data.frame(x = Xr[-1, 1], y = Xr[-1, 2]), file.path(tdir, "plot_bad.rds"))
  st_bad <- system2(file.path(R.home("bin"), "Rscript"),
                    c(shQuote(file.path(explorer_dir, "r", "view_in_explorer.R")), shQuote(file.path(tdir, "fit.rds")),
                      shQuote(file.path(tdir, "X.rds")), paste0("--plot=", shQuote(file.path(tdir, "plot_bad.rds")))),
                    stdout = FALSE, stderr = FALSE)
  check("view_in_explorer.R refuses plot data whose row count differs from X", !identical(st_bad, 0L))
  out_b <- system2("python3", c(shQuote(file.path(explorer_dir, "src", "build.py")), "--data", shQuote(file.path(tdir, "ens.json")),
                                 "--out", shQuote(file.path(tdir, "ens.html"))), stdout = TRUE, stderr = TRUE)
  check("headless route: build.py turns an export into a page", file.exists(file.path(tdir, "ens.html")))
} else cat("SKIP Explorer checks: jsonlite or python3 missing\n")

cat("\n== Setup script (phm/scripts/setup_phm.R) ==\n")
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
setup_script <- normalizePath(file.path(dirname(this_file), "..", "phm", "scripts", "setup_phm.R"), mustWork = FALSE)
check("setup_phm.R exists", file.exists(setup_script))
status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(setup_script), "--check-only"), stdout = FALSE, stderr = FALSE)
check("setup_phm.R --check-only exits 0 when everything is installed", identical(status, 0L))
status_bad <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(setup_script), "--bogus"), stdout = FALSE, stderr = FALSE)
check("setup_phm.R rejects unknown arguments", !identical(status_bad, 0L))
src <- readLines(setup_script)
start <- grep("^required <- c\\(", src)
end <- start + which(src[(start + 1):length(src)] == ")")[1]
script_pkgs <- sub("^\\s+([A-Za-z0-9.]+)\\s+=.*$", "\\1", src[(start + 1):(end - 1)])
imports <- trimws(strsplit(gsub("\n", " ", packageDescription("PHM")$Imports), ",")[[1]])
imports <- sub("\\s*\\(.*\\)$", "", imports)
check("setup_phm.R required-package list matches PHM's DESCRIPTION Imports (minus base 'parallel')",
      setequal(script_pkgs, setdiff(imports, "parallel")))

cat("\n")
if (failures > 0) {
  cat(failures, "check(s) FAILED\n"); quit(status = 1)
} else cat("All checks passed\n")
