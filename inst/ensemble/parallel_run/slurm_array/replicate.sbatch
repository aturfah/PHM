#!/bin/bash
#SBATCH --job-name=phm_replicate
#SBATCH --array=1-40            # submit.sh replaces this with 1-M from config.R
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/replicate_%a.log

## Adjust to how this cluster provides R -- e.g.:
# module load R
# source activate phm-env

## submit.sh passes this script's own directory and the config.R path it
## resolved as $1/$2, so run_replicate.R is found regardless of Slurm's
## working directory for the job and uses the same config.R submit.sh
## used (its --output above is overridden by submit.sh too, to put logs
## next to that config.R rather than here). Defaults let this still work
## if you `sbatch replicate.sbatch` directly from inside this directory,
## with a config.R in that same directory.
SCRIPT_DIR="${1:-.}"
CONFIG_PATH="${2:-config.R}"

Rscript "$SCRIPT_DIR/../run_replicate.R" "$CONFIG_PATH"
