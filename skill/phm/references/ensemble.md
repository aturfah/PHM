# The PHM ensemble

A PHM ensemble fits M Gaussian mixtures, each on a random subsample of the data, pools their components into one mixture, and builds the PHM hierarchy on the pool with `PHMv2()`. The replicates are independent, so the work splits into a **map** step (one replicate, any number in parallel) and a **reduce** step (pool + PHMv2, once).

## Contents

- [Functions](#functions)
- [In an R session](#in-an-r-session)
- [From the command line](#from-the-command-line)
- [config.R fields](#configr-fields)
- [Failures, reruns, reproducibility](#failures-reruns-reproducibility)
- [Choosing settings](#choosing-settings)
- [After the run](#after-the-run)

## Functions

| Function | Role |
|---|---|
| `fitDensityReplicate(X, idx, seed, rho, G, ...)` | Map: one replicate. Draws `round(rho * N)` rows without replacement, picks a covariance model and G by BIC, widening G by `GStep` up to `GMax` when the best G is the largest tried, optionally refits from a random start (`stage2Init`). Returns `list(params, loglik, bic, idx, G, modelName)`. |
| `combineDensityReplicates(replicateResults, replicateWeighting = "equal", ...)` | Reduce: weights replicates (`equal`, `likelihood` or `bic`), pools components (labels `"<component>_<replicate>"`), runs `PHMv2()`. Returns a `PHMv2` object. |
| `fitPHMEnsemble(X, M, rho, ...)` | Both, in one session, replicates one after another. |
| `phmSeeds(seedBase, M)` | The M replicate seeds from one base seed. |
| `globalDensityEstimation(X, M, rho, ...)` | Same replicates, returns the pooled `paramsList` only (parallel `numCores`, on-disk cache `saveDir`). |

## In an R session

```r
library(PHM)
fit <- fitPHMEnsemble(X, M = 10, rho = 0.5, G = 1:25, numCores = 2)
saveRDS(fit, "phm_fit.rds")
```

Equivalent by hand, which is what the command-line scripts do:

```r
seeds <- phmSeeds(1000, 10)
reps  <- lapply(1:10, function(i) fitDensityReplicate(X, idx = i, seed = seeds[i], rho = 0.5, G = 1:25))
fit   <- combineDensityReplicates(reps, replicateWeighting = "equal", numCores = 2)
```

To fit replicates in parallel inside R, use `parallel::mclapply` over `idx`, or `globalDensityEstimation(X, M = 10, rho = 0.5, numCores = 4, saveDir = "reps")` followed by `PHMv2(paramsList = params, mcSamples = 2e5, batchSize = 1e5)`; with default seeds the pooled components are identical to `fitPHMEnsemble()`.

## From the command line

The package installs scripts that read one `config.R` per run. Find them with

```sh
ENS=$(Rscript -e 'cat(system.file("ensemble", package = "PHM"))')
```

(in a clone of the repo they are in `inst/ensemble/`). Each folder there has a README.

```sh
mkdir my_run && cd my_run
Rscript -e 'saveRDS(as.matrix(read.csv("data.csv")), "X.rds")'     # input: one .rds holding an N x D numeric matrix
cp $ENS/config/config_template.R config.R                          # edit dataPath, resultsDir, outPath, M, rho, ...

Rscript $ENS/run_serial.R                                          # 1. serial, one R session (fitPHMEnsemble)
MAX_CONCURRENT=5 bash $ENS/parallel_run/single_machine/launch.sh   # 2. one machine, 5 replicates at a time, then combine
bash $ENS/parallel_run/slurm_array/submit.sh                       # 3. Slurm: array of M tasks + combine job (afterok)
MAX_CONCURRENT=10 bash $ENS/parallel_run/slurm_array/submit.sh     #    ... at most 10 tasks running at once

Rscript $ENS/parallel_run/run_replicate.R 3                        # fit replicate 3 only (test, or rerun)
Rscript $ENS/parallel_run/run_combine.R                            # pool whatever replicate_*.rds exist for 1..M
```

Every script takes the config path as its last optional argument (default `./config.R`); relative paths inside `config.R` are relative to the directory the scripts run from. `launch.sh` first runs `check_requirements.R`, which reports every setup problem at once. On Slurm, `run_replicate.R` reads its index from `SLURM_ARRAY_TASK_ID`; edit `--mem`, `--time` and module lines in `replicate.sbatch` / `combine.sbatch` for the cluster, and make `--cpus-per-task` in `combine.sbatch` match `numCores`.

Outputs: `resultsDir/replicate_001.rds`, ... (one small file per replicate, never the data), `outPath` (the `PHMv2` object), and `logs/` next to `config.R`.

The worked example `$ENS/example/double_ring/` (README + `config.R` with M = 10, rho = 0.5) runs the bundled `double_ring` data all three ways; start there when a user is new to the pipeline.

## config.R fields

`config.R` defines one list, `CFG`:

| Field | Default in template | Meaning |
|---|---|---|
| `dataPath` | placeholder | `.rds` of the N × D numeric matrix (finite, at least 2 columns; names ignored). Every replicate reads it whole, then subsamples. |
| `resultsDir` | placeholder | Folder for `replicate_<idx>.rds`. |
| `outPath` | placeholder | Where the final `PHMv2` object is saved. |
| `M` | 40 | Number of replicates. |
| `rho` | 0.25 | Fraction of rows per replicate, without replacement, in (0, 1]. |
| `G` | `1:25` | Starting range of mixture sizes (extended automatically). |
| `modelNames` | NULL (all) | mclust covariance families allowed, e.g. `c("VII", "VEI", "EII")` in high dimension. Validated up front. |
| `GStep`, `GMax` | 5, 100 | Extension step and hard cap for G. Hitting `GMax` prints a message suggesting a larger cap. |
| `extendVerbose` | FALSE | Print G-extension progress; turn on for a new dataset. |
| `seedBase` | 1000 | Base for `phmSeeds()`. |
| `stage2Init` | NULL | Stage 2 EM restart: NULL (off), `"kmeans"`, `"kmeans++"`, `"random"`, `"randomFromData"`. |
| `nRestartTries` | 5 | Stage 2 retries before falling back to the Stage 1 fit. |
| `replicateWeighting` | `"equal"` | `"equal"`, `"likelihood"`, or `"bic"` (softmax of log-likelihood or BIC/2 across replicates). |
| `mcSamples`, `batchSize` | 2e5, 1e5 | Monte Carlo settings for the reduce-step `PHMv2()`. |
| `numCores` | 2 | Cores for the reduce-step Monte Carlo only. Replicate concurrency is `MAX_CONCURRENT` (one machine) or the array (Slurm). |
| `label` | `"cluster_run"` | Tag on log lines. |

## Failures, reruns, reproducibility

- Replicate `i` always uses `phmSeeds(seedBase, M)[i]`, so rerunning `run_replicate.R i` gives the same result, and serial, one-machine and Slurm runs of the same `config.R` give identical pooled components.
- If any replicate fails, `launch.sh` and the Slurm dependency stop before the reduce step. Look at `logs/`, fix the cause, rerun only the failed indices with `run_replicate.R`, then `run_combine.R`. `run_combine.R` stops if any `replicate_<idx>.rds` for 1..M is missing.
- A replicate whose G extension fails (an all-NA BIC block) keeps its best fit and logs "extension to G=a:b unusable ... keeping best fit"; that is not an error.
- In R, forked parallel steps (`numCores > 1`) stop with "k of n parallel task(s) failed (first: ...)" when a worker errors or is killed. "no result returned" usually means the OS killed it for memory: lower `numCores` or `batchSize`.
- Changing `M`, `seedBase`, `rho`, `G` or `modelNames` changes the replicates: use a fresh `resultsDir`, or old `replicate_*.rds` files will be pooled with new ones.

## Choosing settings

- **rho and M.** Each replicate must see enough points to fit its largest model: `round(rho * N)` well above the parameter count (for VVV, G·(D + D(D+1)/2 + 1)). Typical: rho 0.25–0.5 with M 10–50. Smaller rho is faster per replicate and gives more diversity; larger M makes the pooled density smoother and the ΔP_mc matrix bigger.
- **N ≤ 2000 with rho = 1 and no Stage 2** gives identical replicates (mclust's hierarchical initialization is deterministic on ≤ 2000 rows). Use rho < 1 or `stage2Init`. Above 2000 rows mclust initializes on a random 2000-row subset, so replicates differ even at rho = 1.
- **Dimension.** Run on PCs (5–30). In 20+ dimensions restrict `modelNames` to diagonal or spherical families, which keeps fits stable and fast.
- **Reduce-step cost.** The pool has K = Σ G_i components (often hundreds to a few thousand). ΔP_mc costs about `mcSamples` × K density evaluations and K² memory; each forked worker holds a `batchSize` × K posterior block. With K in the thousands use `batchSize = 2e4`–`5e4`, few cores, and enough memory on the combine job.
- **Weighting.** `equal` is the default and the safest; `likelihood`/`bic` favor better-fitting replicates and are comparable only because every replicate uses the same subsample size.

## After the run

```r
fit <- readRDS("phm_fit.rds")
length(fit$paramsList)                         # pooled components
pmc_by_k <- cumsum(fit$mergeDeltaPmc)
plotPHMv2Dendrogram(fit, initK = 30, displayAxis = "none")

grp  <- recoverGroupsv2(fit, 5)
mass <- sapply(grp, function(g) sum(sapply(fit$paramsList[g], function(p) sum(p$prob))))
mass                                            # tiny-mass groups are usually outlier components, not clusters

post <- recoverPosteriorv2(fit, 5, computePosteriorProbMatrix(fit$paramsList, X))
lab  <- apply(post, 1, which.max)
```

For interactive browsing: `Rscript <skill>/scripts/explorer/r/view_in_explorer.R phm_fit.rds X.rds` (see `explorer.md`).
