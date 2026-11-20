#!/bin/bash
#SBATCH --job-name=phm_combine
#SBATCH --cpus-per-task=2       # should match CFG$numCores in config.R
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/combine.log

## Adjust to how this cluster provides R -- e.g.:
# module load R
# source activate phm-env

## submit.sh passes this script's own directory and the config.R path it
## resolved as $1/$2 -- see replicate.sbatch's comment for why, and the
## same defaults for running this directly without submit.sh.
SCRIPT_DIR="${1:-.}"
CONFIG_PATH="${2:-config.R}"

Rscript "$SCRIPT_DIR/../run_combine.R" "$CONFIG_PATH"
