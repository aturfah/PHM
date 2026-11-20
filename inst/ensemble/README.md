# PHM ensemble: run scripts

Scripts for running a PHM ensemble analysis from the command line, in serial,
in parallel on one machine, or on a Slurm cluster. The fitting itself is done
by the PHM package's ensemble functions:

| Function | Role |
|---|---|
| `fitDensityReplicate()` | Map step: fits one replicate (a Gaussian mixture on a random subsample) and returns its components, log-likelihood and BIC. |
| `combineDensityReplicates()` | Reduce step: pools the replicates' components and runs `PHMv2()` on them. |
| `fitPHMEnsemble()` | Both steps in one R session, one replicate after another. |
| `phmSeeds()` | The replicate seeds, drawn from one base seed. |

`globalDensityEstimation()` uses the same `fitDensityReplicate()` for each
replicate, so the same settings give the same pooled components there too.
See each function's help page (for example `?fitDensityReplicate`) for the
details.

## Where the scripts are

They are installed with the package. Set `ENS` to their folder once per shell:

```sh
ENS=$(Rscript -e 'cat(system.file("ensemble", package = "PHM"))')
```

In a clone of the PHM repo they are also at `inst/ensemble/`.

| Path under `$ENS` | What it holds |
|---|---|
| [`config/config_template.R`](config/) | The one file of settings every script reads. Copy it to `config.R` in a working folder. |
| `run_serial.R` | Fits the whole ensemble in one R session with `fitPHMEnsemble()`. |
| [`parallel_run/`](parallel_run/) | `run_replicate.R` and `run_combine.R`, shared by the two parallel styles: [`single_machine/`](parallel_run/single_machine/) and [`slurm_array/`](parallel_run/slurm_array/). |
| [`example/`](example/) | A worked example on the package's `double_ring` dataset, run all three ways. |

## Quick start

Keep each analysis in its own working folder:

```sh
mkdir my_run && cd my_run
cp $ENS/config/config_template.R config.R
# edit config.R: at least dataPath, resultsDir, outPath and M
bash $ENS/parallel_run/single_machine/launch.sh
```

`launch.sh` checks the setup, fits the `M` replicates a few at a time, and then
pools them; the fit is written to `outPath`. On a Slurm cluster use
`bash $ENS/parallel_run/slurm_array/submit.sh` instead, and to fit everything
in one R session, `Rscript $ENS/run_serial.R`. To test a single replicate:

```sh
Rscript $ENS/parallel_run/run_replicate.R 1
```

## Input data

One `.rds` file holding a plain numeric matrix, one row per observation and one
column per variable:

```r
saveRDS(X, "X.rds")    # X is an N x D numeric matrix
```

- At least two columns, and every value finite (no `NA`, `NaN` or `Inf`).
- No labels, factors or extra columns. Row and column names are ignored.
- Scaling is up to you. Spherical models (`EII`, `VII`) are sensitive to it;
  `VVV` and `EEE` are not.
- Each replicate fits on `round(rho * N)` rows, so that subsample has to be
  well above the parameter count of the largest model allowed by `G` and
  `modelNames`.

Every replicate reads the whole matrix and then subsamples it, so the matrix
has to fit in memory wherever a replicate runs.

## How a run works

1. **Seeds.** `phmSeeds(seedBase, M)` gives the `M` seeds. Every script
   recomputes the same vector, so replicate `i` always gets the same seed and
   no seed file is passed around.
2. **Map: one replicate per job.** `fitDensityReplicate()` draws
   `round(rho * N)` rows without replacement and chooses a covariance model
   and number of components `G` by BIC. If the best `G` is at the top of the
   range tried, it widens the range in steps of `GStep` up to `GMax`.
   Optionally (`stage2Init`), it refits that model by EM from a random start.
   Each replicate writes a small `replicate_<idx>.rds` with its components,
   log-likelihood and BIC, never the data.
3. **Reduce: once.** `combineDensityReplicates()` reads all `M` results,
   weights the replicates (`equal` by default, or `likelihood` / `bic`),
   pools their components into one mixture and runs `PHMv2()` on it. The
   result is saved to `outPath`.

If any replicate fails, the reduce step does not run. Fix the cause, rerun
only the failed replicates (each reuses its own seed), then run
`run_combine.R`.

## Output

- `resultsDir/replicate_001.rds`, ...: one per replicate, a list with
  `params`, `loglik`, `bic`, `idx`, `G` and `modelName`.
- `outPath`: a `PHMv2` object with `deltaPmc`, `paramsList` (the pooled
  components), `mergeComps`, `mergeValues`, `mergeDeltaPmc` and
  `mergeCriterion`. Cut the tree at `k` clusters with `recoverGroupsv2(fit, k)`
  or `recoverPosteriorv2()`, plot it with `plotPHMv2Dendrogram()`, or open it
  in PHM Explorer with `PHM_Explorer/r/view_in_explorer.R`.
- `logs/`: one log per replicate (single machine) or per array task (Slurm),
  next to the `config.R` the run used.
