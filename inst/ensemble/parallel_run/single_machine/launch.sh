#!/usr/bin/env bash
## launch.sh -- runs all M replicates as background OS processes on this
## machine, in batches of MAX_CONCURRENT at a time (default: this
## machine's core count), then runs the combine step once they've all
## finished successfully. This is the only "orchestration" here -- it's a
## plain shell script that runs in the foreground of your own terminal
## session for as long as fitting takes; there's no scheduler tracking
## anything for you the way Slurm does in slurm_array/.
##
## config.R is looked up relative to wherever YOU run this from (your
## current working directory) by default, not relative to this script --
## so a run's own config.R, data, and results can all live together in a
## working directory of your choosing, with this script (and the shared
## ../run_replicate.R and ../run_combine.R) invoked by path:
##   cd /wherever/your/run/lives        # config.R is here
##   bash $ENS/parallel_run/single_machine/launch.sh
## where $ENS is inst/ensemble in a PHM checkout, or
## system.file("ensemble", package = "PHM") in an installed copy.
## or pass a config.R that lives somewhere else entirely as an argument,
## from any directory:
##   ./launch.sh /path/to/your/run/config.R
##   MAX_CONCURRENT=4 ./launch.sh /path/to/your/run/config.R
##
## logs/ is created next to whichever config.R this run actually used
## (not next to this script), so a run's log files stay with that run's
## own working directory rather than piling up inside the code checkout.
##
## Before launching anything, runs check_requirements.R: confirms the PHM
## package (and so everything it imports) is installed and CFG$dataPath
## exists -- once, up front, with one clear message, rather than finding out
## the same way M times over via M separately crashed background replicates
## and log files.
##
## Batches (rather than immediately refilling a slot the instant any one
## job finishes) on purpose: that finer-grained scheduling needs `wait -n`,
## which stock macOS ships without (its default /bin/bash is 3.2; `wait -n`
## needs bash >= 4.3). This batched version runs correctly on that default
## shell. If you've installed a newer bash (e.g. via Homebrew) and want
## jobs to start the moment a slot frees up rather than waiting for the
## whole batch, swap the inner `wait` loop for `wait -n` in a counting
## loop, or just use GNU parallel / `xargs -P` if you have either.

set -euo pipefail

## Where this script lives; the shared run_replicate.R/run_combine.R are one
## level up. Resolved so this can be invoked from any working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

CONFIG_PATH="${1:-config.R}"
if [ ! -f "$CONFIG_PATH" ]; then
  echo "config file '$CONFIG_PATH' not found -- run this from the directory holding config.R, or pass its path: ./launch.sh <config.R path>" >&2
  exit 1
fi
CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"

Rscript "$SCRIPT_DIR/check_requirements.R" "$CONFIG_PATH" "$SHARED_DIR"

LOG_DIR="$(dirname "$CONFIG_PATH")/logs"
mkdir -p "$LOG_DIR"

M=$(Rscript -e "source(\"$CONFIG_PATH\"); cat(CFG\$M)")
maxJobs=${MAX_CONCURRENT:-$(Rscript -e 'cat(parallel::detectCores())')}
echo "running $M replicates, up to $maxJobs at a time (set MAX_CONCURRENT to change this)"

idx=1
anyFailed=0
while [ "$idx" -le "$M" ]; do
  batchEnd=$(( idx + maxJobs - 1 ))
  if [ "$batchEnd" -gt "$M" ]; then batchEnd=$M; fi

  pids=()
  for i in $(seq "$idx" "$batchEnd"); do
    Rscript "$SHARED_DIR/run_replicate.R" "$i" "$CONFIG_PATH" > "$LOG_DIR/replicate_${i}.log" 2>&1 &
    pid=$!
    pids+=("$pid")
    echo "  launched replicate $i (pid $pid) -> $LOG_DIR/replicate_${i}.log"
  done

  for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
      anyFailed=1
    fi
  done

  idx=$(( batchEnd + 1 ))
done

if [ "$anyFailed" -ne 0 ]; then
  echo "one or more replicates failed -- check $LOG_DIR/replicate_*.log before running the combine step" >&2
  exit 1
fi

echo "all $M replicates finished -- running the combine step"
Rscript "$SHARED_DIR/run_combine.R" "$CONFIG_PATH"
