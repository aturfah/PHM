# config/

`config_template.R` defines one list, `CFG`, holding every setting a run
needs. Copy it to `config.R` in a working folder for each analysis and edit
the copy:

```sh
cp $ENS/config/config_template.R my_run/config.R
```

(`$ENS` is the scripts folder; see [`../README.md`](../README.md).) Every
script reads the same `config.R`, so the replicate jobs and the combine job
cannot disagree about `M`, `seedBase`, `rho` and the rest. The template
comments each field in detail; this page is the summary.

## Fields

**Paths**

| Field | Meaning | Template value |
|---|---|---|
| `dataPath` | The input `.rds` (an N x D numeric matrix). | placeholder |
| `resultsDir` | Where each replicate writes `replicate_<idx>.rds`. Created if missing. | placeholder |
| `outPath` | Where the final `PHMv2` fit is saved. | placeholder |

**Replicates (map step)**

| Field | Meaning | Template value |
|---|---|---|
| `M` | Number of replicates. | `40` |
| `rho` | Fraction of rows each replicate draws, without replacement, in (0, 1]. | `0.25` |
| `G` | Starting range for the number of components. | `1:25` |
| `modelNames` | mclust covariance models to search, e.g. `"EEE"`. `NULL` searches mclust's default set. | `NULL` |
| `GStep`, `GMax` | How far the `G` range widens when the best `G` sits at its top: `GStep` values at a time, never past `GMax`. | `5`, `100` |
| `extendVerbose` | Log each widening step. | `FALSE` |
| `seedBase` | Seed from which all `M` replicate seeds are drawn. | `1000` |
| `stage2Init` | Optional EM refit from a random start: `NULL` (skip), `"kmeans"`, `"kmeans++"`, `"random"` or `"randomFromData"`. | `NULL` |
| `nRestartTries` | Stage 2 attempts before falling back to the plain BIC fit. | `5` |

**Pooling (reduce step)**

| Field | Meaning | Template value |
|---|---|---|
| `replicateWeighting` | `"equal"` (1/M each), `"likelihood"` or `"bic"` (softmax over the replicates' log-likelihood or BIC/2). | `"equal"` |
| `mcSamples`, `batchSize` | Monte Carlo sample size and batch size for the ΔPmc matrix. | `2e5`, `1e5` |
| `numCores` | Cores for the ΔPmc computation in the combine step only. It has no effect on how many replicates run at once. | `2` |
| `label` | Tag printed on every log line. | `"cluster_run"` |

## Things to know

- **Relative paths** in `dataPath`, `resultsDir` and `outPath` are resolved
  from the directory you run the scripts in. Run from the folder that holds
  `config.R`, or use absolute paths.
- **Subsample size.** With large `N`, a small `rho` keeps replicates fast;
  with small `N`, make sure `round(rho * N)` is well above the parameter
  count of the largest model `G` and `modelNames` allow.
- **`"likelihood"` and `"bic"` weighting** compare replicates fairly only
  because every replicate uses a subsample of the same size. With
  `"likelihood"` or `"bic"` on large subsamples, one replicate often takes
  nearly all the weight; `"equal"` keeps every replicate's say.
- **Slurm:** `parallel_run/slurm_array/submit.sh` reads `M` from this file
  and sizes the job array to match, so there is nothing to keep in sync.
- **Older config files** from the standalone phm_ensemble_code scripts still
  work: a `codeBase` field is ignored (with a note), because the fitting code
  now comes from the installed PHM package, and a `phmSeeds()` defined in the
  file is not used.
