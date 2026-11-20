# parallel_run/single_machine/

Runs the replicates in parallel on one computer, with no scheduler. There are
two ways to do it; pick one.

| Script | How replicates run |
|---|---|
| `launch.sh` | Each replicate is its own background `Rscript` process with its own log file. Uses the shared `../run_replicate.R` and `../run_combine.R`. |
| `local_driver.R` | All replicates run inside one R session through `parallel::mclapply()`. |

Both give the same result from the same `config.R`.

## launch.sh

```sh
cd my_run                                                   # holds config.R
bash $ENS/parallel_run/single_machine/launch.sh
MAX_CONCURRENT=4 bash $ENS/parallel_run/single_machine/launch.sh   # at most 4 at once
bash $ENS/parallel_run/single_machine/launch.sh /elsewhere/config.R
```

What it does, in order:

1. Runs `check_requirements.R`, which checks that the shared scripts exist,
   that the PHM package (with everything it imports) loads, and that
   `dataPath` exists and is readable. It reports every problem at once and
   stops before starting any replicate. It does not open the data file, so a wrongly shaped matrix
   shows up only when the first replicate fails.
2. Starts replicates in batches of `MAX_CONCURRENT` (default: the number of
   cores) and waits for each whole batch before starting the next. Batching
   keeps the script compatible with the bash 3.2 that ships with macOS.
3. Records any replicate that exits with an error and keeps going, so every
   replicate gets its turn.
4. If all succeeded, runs `run_combine.R`. If any failed, prints a message,
   skips the combine step and exits with status 1.

Logs: `logs/replicate_<i>.log`, next to `config.R`. The combine step prints
to the terminal.

### Running pieces by hand

```sh
Rscript $ENS/parallel_run/run_replicate.R <idx> [config.R]   # one replicate, idx in 1..M
Rscript $ENS/parallel_run/run_combine.R [config.R]           # pool all M results
```

This is how to finish a run after a failure: rerun the failed indices (each
gets its original seed), then run the combine script. `run_replicate.R`
overwrites that one replicate's result file. `run_combine.R` stops and names
the first missing file if any of the `M` results are absent.

## local_driver.R

```sh
Rscript $ENS/parallel_run/single_machine/local_driver.R [config.R]
MAX_CONCURRENT=4 Rscript $ENS/parallel_run/single_machine/local_driver.R
```

Fits every replicate in one command, then pools them and writes `outPath`.
Each replicate gets its own forked process (`mc.preschedule = FALSE`), so one
replicate's error does not take down the others. If any replicate fails,
including a fork killed by the operating system (usually for running out of
memory), the script stops before pooling and lists which ones failed.

Limits compared with `launch.sh`:

- No per-replicate log files, and nothing is saved until every replicate is
  done, so a failure means rerunning the whole set.
- Does not work on Windows, which cannot fork processes.

## Memory

Every replicate running at the same time needs its own working memory, and
under `launch.sh` each one also loads the full data matrix. For large data,
set `MAX_CONCURRENT` by what the machine's RAM can hold, not by its core
count.

## Files

| File | Role |
|---|---|
| `launch.sh` | Orchestrates a full run, as above. |
| `check_requirements.R` | Preflight checks. `Rscript check_requirements.R <config.R>` runs it alone. |
| `local_driver.R` | The one-session alternative to all of the above. |
