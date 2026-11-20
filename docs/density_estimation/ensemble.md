# Running a PHM Ensemble

A PHM ensemble fits `M` Gaussian mixture models, each to a random subsample of
the data, pools their components into one mixture, and builds the PHM
hierarchy on the pooled mixture with `PHMv2()`. Because the replicates are
independent, the fitting splits into a *map* step, run once per replicate and
in parallel if you like, and a *reduce* step, run once at the end.

| Function | Role |
|---|---|
| `fitDensityReplicate()` | Map step: fits one replicate. Chooses a covariance model and `G` by BIC, widening the `G` range when the best `G` is the largest tried, optionally refits from a random start (`stage2Init`), and returns the replicate's components, log-likelihood and BIC. |
| `combineDensityReplicates()` | Reduce step: weights the replicates (`"equal"` by default, or `"likelihood"` / `"bic"`), pools their components and runs `PHMv2()`. |
| `fitPHMEnsemble()` | Both steps in one R session, one replicate after another. |
| `phmSeeds()` | The `M` replicate seeds, drawn from one base seed, so any replicate can be refit alone with the same result. |

`globalDensityEstimation()` fits its replicates with the same
`fitDensityReplicate()`; see [Global and Localized Density
Estimation](global_localized.md).

## In an R session

```r
library(PHM)

fit <- fitPHMEnsemble(double_ring, M = 10, rho = 0.5, G = 1:25)

length(fit$paramsList)              # pooled components from all 10 replicates
groups <- recoverGroupsv2(fit, 2)   # cut the tree at two clusters
plotPHMv2Dendrogram(fit)
```

Each replicate fits `round(rho * N)` rows drawn without replacement, with seed
`phmSeeds(seedBase, M)[idx]` (default `seedBase = 1000`). The defaults are no
Stage 2 refit (`stage2Init = NULL`) and equal replicate weights
(`replicateWeighting = "equal"`).

## From the command line

For long runs, or to use several cores or a cluster, the package installs a
set of scripts that read their settings from one `config.R` file:

```sh
ENS=$(Rscript -e 'cat(system.file("ensemble", package = "PHM"))')

mkdir my_run && cd my_run
cp $ENS/config/config_template.R config.R     # then edit dataPath, resultsDir, outPath, M, ...

Rscript $ENS/run_serial.R                                        # one session, serial
MAX_CONCURRENT=5 bash $ENS/parallel_run/single_machine/launch.sh # one machine, 5 at a time
bash $ENS/parallel_run/slurm_array/submit.sh                     # Slurm job array of M tasks
Rscript $ENS/parallel_run/run_replicate.R 1                      # test one replicate
```

The input is one `.rds` file holding an N x D numeric matrix
(`saveRDS(X, "X.rds")`). Every style gives the same pooled components from the
same `config.R`. If a replicate fails, the pooling step does not run; rerun the
failed replicate (it reuses its seed), then run
`Rscript $ENS/parallel_run/run_combine.R`.

The scripts are in [`inst/ensemble/`](https://github.com/aturfah/PHM/tree/main/inst/ensemble)
of the repository, with a README in each folder. The
[`example/double_ring/`](https://github.com/aturfah/PHM/tree/main/inst/ensemble/example/double_ring)
folder walks through a complete run on the `double_ring` dataset in all three
styles.

## Viewing the result

The fit saved to `outPath` is a `PHMv2` object. Besides the package's own
functions (`recoverGroupsv2()`, `recoverPosteriorv2()`,
`plotPHMv2Dendrogram()`, `plotPHMv2Heatmap()`), PHM Explorer opens it in the
browser with

```sh
Rscript PHM_Explorer/r/view_in_explorer.R my_run/phm_fit.rds my_run/X.rds
```
