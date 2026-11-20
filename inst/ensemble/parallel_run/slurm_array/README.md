# parallel_run/slurm_array/

Runs each replicate as one task of a Slurm job array, then pools them in a
second job that Slurm starts only after every array task has succeeded.
Nothing of yours has to stay running while the jobs wait or run.

`$ENS` below is the scripts folder; see [`../../README.md`](../../README.md). The
repo (or installed package) must be on a filesystem the compute nodes can read,
and PHM must be installed where the nodes' R can load it.

## Before the first run

1. **Set up R on the nodes.** Both `.sbatch` files have commented lines such
   as `# module load R`. Replace them with whatever your cluster uses so that
   `Rscript` and the packages from the top-level README are available.
2. **Set resources.** The defaults are 1 CPU, 4 GB and 2 hours per replicate,
   and 2 CPUs, 8 GB and 30 minutes for the combine job. Raise memory and time
   for large data. `--cpus-per-task` in `combine.sbatch` should match
   `numCores` in `config.R`.

## Running

Run `submit.sh` on a login node with R loaded, because it reads `M` from
`config.R`:

```sh
module load R                                 # or however your cluster provides R
cd my_run                                     # holds config.R
bash $ENS/parallel_run/slurm_array/submit.sh
# or: bash $ENS/parallel_run/slurm_array/submit.sh /elsewhere/config.R
MAX_CONCURRENT=5 bash $ENS/parallel_run/slurm_array/submit.sh   # at most 5 tasks at once
```

`submit.sh` submits a job array of `M` tasks (`--array=1-M`, taken from
`config.R`, so the `--array` line inside `replicate.sbatch` is only a fallback
for submitting that file by hand). With `MAX_CONCURRENT=N` set, the array is
`1-M%N`. It then submits the combine job with
`--dependency=afterok:<array job id>`, prints both job ids and exits. Logs go
to `logs/replicate_<task id>.log` and `logs/combine.log` next to `config.R`.

## When a task fails

The combine job then never starts; Slurm leaves it pending with reason
`DependencyNeverSatisfied`. Cancel it with `scancel <combine job id>`, fix
the cause, and resubmit only the failed tasks. Each task recomputes its own
seed, so rerunning task 11 gives exactly what task 11 would have given:

```sh
cd my_run
SD=$ENS/parallel_run/slurm_array
sbatch --array=11,24 --output=logs/replicate_%a.log $SD/replicate.sbatch $SD $PWD/config.R
# once those succeed:
sbatch --output=logs/combine.log $SD/combine.sbatch $SD $PWD/config.R
```

## Files

| File | Role |
|---|---|
| `submit.sh` | Reads `M` from `config.R`, submits the `M`-task array job and the dependent combine job. |
| `replicate.sbatch` | Array job script. Takes this folder and the `config.R` path as arguments (`submit.sh` passes both). |
| `combine.sbatch` | Combine job script, same arguments. |

Both job scripts call the shared `../run_replicate.R` and `../run_combine.R`.
Without an index argument, `run_replicate.R` takes the replicate number from
`SLURM_ARRAY_TASK_ID`.

These scripts were tested with a stand-in for `sbatch` that runs each array
task with its `SLURM_ARRAY_TASK_ID` set and starts the combine job only after
every task succeeds. They have not yet been run on a real Slurm cluster.
