#!/usr/bin/env bash
## submit.sh -- the only "orchestration" piece in this style: submits the
## map array job, then the reduce job with a dependency on it. Runs once,
## briefly, and exits -- Slurm itself tracks which array tasks finished
## and enforces the dependency, so nothing here needs to keep running.
##
## config.R is looked up relative to wherever YOU run this from (your
## current working directory) by default, not relative to this script --
## so a run's own config.R, data, and results can live together in a
## working directory of your choosing, with this script invoked by path:
##   cd /wherever/your/run/lives        # config.R is here
##   bash $ENS/parallel_run/slurm_array/submit.sh   ($ENS: see ../README.md)
## or pass a config.R that lives somewhere else entirely as an argument,
## from any directory:
##   bash $ENS/parallel_run/slurm_array/submit.sh /path/to/your/run/config.R
##
## logs/ is created next to whichever config.R this run actually used
## (not next to this script).

set -euo pipefail

## Where this script (and replicate.sbatch/combine.sbatch, its siblings)
## actually live -- resolved so submit.sh can be invoked from any working
## directory, passed through to the sbatch job scripts below so THEY can
## find the shared ../run_replicate.R and ../run_combine.R regardless of
## Slurm's own working directory for the job.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_PATH="${1:-config.R}"
if [ ! -f "$CONFIG_PATH" ]; then
  echo "config file '$CONFIG_PATH' not found -- run this from the directory holding config.R, or pass its path: ./submit.sh <config.R path>" >&2
  exit 1
fi
CONFIG_PATH="$(cd "$(dirname "$CONFIG_PATH")" && pwd)/$(basename "$CONFIG_PATH")"

LOG_DIR="$(dirname "$CONFIG_PATH")/logs"
mkdir -p "$LOG_DIR"

## The array size comes from CFG$M, so the number of tasks always matches
## config.R; this overrides the #SBATCH --array line in replicate.sbatch.
## Needs Rscript on the machine you submit from (load R there first if your
## cluster uses modules). MAX_CONCURRENT=N, if set, caps how many tasks run
## at once (Slurm's --array=1-M%N); unset, Slurm runs as many as it can.
if ! command -v Rscript > /dev/null; then
  echo "Rscript not found -- load R on this machine first (e.g. module load R); submit.sh reads M from config.R" >&2
  exit 1
fi
M=$(Rscript -e "source(\"$CONFIG_PATH\"); cat(CFG\$M)")
if ! [[ "$M" =~ ^[0-9]+$ ]] || [ "$M" -lt 1 ]; then
  echo "could not read a positive CFG\$M from $CONFIG_PATH (got '$M')" >&2
  exit 1
fi
ARRAY="1-$M"
if [ -n "${MAX_CONCURRENT:-}" ]; then ARRAY="$ARRAY%$MAX_CONCURRENT"; fi

mapJobId=$(sbatch --parsable \
  --array="$ARRAY" \
  --output="$LOG_DIR/replicate_%a.log" \
  "$SCRIPT_DIR/replicate.sbatch" "$SCRIPT_DIR" "$CONFIG_PATH")
combineJobId=$(sbatch --parsable --dependency=afterok:"$mapJobId" \
  --output="$LOG_DIR/combine.log" \
  "$SCRIPT_DIR/combine.sbatch" "$SCRIPT_DIR" "$CONFIG_PATH")

echo "map array job:     $mapJobId  ($M tasks, --array=$ARRAY; see $LOG_DIR/replicate_*.log)"
echo "combine job:        $combineJobId  (queued to run after $mapJobId; see $LOG_DIR/combine.log)"
echo "final fit will land at the outPath set in $CONFIG_PATH once $combineJobId finishes"
