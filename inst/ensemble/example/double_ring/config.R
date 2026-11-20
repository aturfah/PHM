## config.R for the double-ring example: M = 10 replicates, each fit to a
## random half of PHM's double_ring dataset (rho = 0.5), with the template's
## defaults for everything else (no Stage 2 refit, equal replicate weights).
## Every field is explained in config/config_template.R; only the choices
## specific to this example are commented here.
##
## Paths are relative, so run every command from the folder holding this
## file. The data file is written in the first step of README.md:
##   Rscript -e 'saveRDS(PHM::double_ring, "double_ring_X.rds")'

CFG <- list(
  dataPath   = "double_ring_X.rds",   # 5000 x 2 matrix, two concentric rings
  resultsDir = "results",
  outPath    = "double_ring_phm_fit.rds",

  M   = 10,
  rho = 0.5,     # each replicate fits a random half: 2500 of the 5000 points

  ## A uniform ring is not Gaussian, so each ring needs many components to
  ## tile it. The search starts at 1:25 and widens on its own if BIC wants more.
  G = 1:25,
  modelNames = NULL,
  GStep = 5,
  GMax = 100,
  extendVerbose = TRUE,

  seedBase = 1000,

  ## No Stage 2 refit: each replicate already differs because it fits its own
  ## random half of the data.
  stage2Init = NULL,
  nRestartTries = 5,

  replicateWeighting = "equal",
  mcSamples = 2e5,
  batchSize = 1e5,
  numCores = 2,

  label = "double_ring"
)
