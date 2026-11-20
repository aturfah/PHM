# parallel_run/

Scripts that run the replicates in parallel and then pool them. Both styles
read the same `config.R` and produce the same result from the same settings;
they differ only in where the replicate jobs run. Both use the two shared
scripts in this folder:

| Script | Role |
|---|---|
| `run_replicate.R` | Fits one replicate: `Rscript run_replicate.R [idx] [config.R]`. Without `idx`, it takes the index from `SLURM_ARRAY_TASK_ID`. |
| `run_combine.R` | Pools the `M` results into the final fit: `Rscript run_combine.R [config.R]`. |

| Folder | Use it when |
|---|---|
| [`single_machine/`](single_machine/) | You have one computer and want to use its cores. No scheduler needed. |
| [`slurm_array/`](slurm_array/) | You have a cluster that runs Slurm. Each replicate is one task of a job array. |

## What both styles share

- **Where `config.R` comes from.** Every script looks for `config.R` in the
  current directory, or takes its path as an argument. The usual pattern is
  to keep a run's `config.R`, data and results together in one working
  folder and call the scripts by their full path from there.
- **Where the code comes from.** Every script runs `library(PHM)`, so the
  fitting code is the installed PHM package. The scripts find each other by
  their own location, so call them by full path (`$ENS/...`) from any working
  directory.
- **Seeds.** Replicate `i` always uses the same seed, recomputed from
  `seedBase` and `M`, so any replicate can be rerun alone and gives the same
  result.
- **Files.** Each replicate writes `resultsDir/replicate_<idx>.rds` (three
  digits, e.g. `replicate_007.rds`). The combine script reads exactly
  `replicate_001.rds` to `replicate_<M>.rds` and writes `outPath`. Logs go to
  `logs/` next to `config.R`.
- **Failures.** If any replicate fails, the combine step does not run, so a
  partial ensemble is never pooled. Fix the cause, rerun the failed
  replicates, then run the combine script by hand.

One caution: the combine script checks only that each expected file exists.
If `resultsDir` still holds files from an earlier run with different
settings, they would be pooled too. Use a fresh `resultsDir` for each run.
