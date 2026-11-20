## config_template.R -- copy this to config.R in a working directory for
## one analysis run, then edit it. Every script in PHM's inst/ensemble/
## (run_serial.R, run_replicate.R, run_combine.R, local_driver.R, launch.sh,
## submit.sh) reads the same config.R, so they can never disagree about M,
## seedBase, rho, etc. -- there is exactly one place these numbers live.
## The fitting code itself comes from the installed PHM package.
##
## dataPath/resultsDir/outPath are placeholders that need editing for every
## new dataset or run. Relative paths are taken relative to the directory
## you run the scripts from, which is normally the one holding config.R.

CFG <- list(

  ## ============================================================
  ## PATHS -- consumed by the scripts, never by PHM's
  ## fitDensityReplicate()/combineDensityReplicates() themselves: those two
  ## functions take X and a plain argument list directly, not a CFG object.
  ## These three paths are the scripts' own bookkeeping, not part of the
  ## package API.
  ## ============================================================

  # Where the dataset lives: a saved N x D numeric matrix (no labels, no
  # wrapping list -- just the matrix itself, e.g. saveRDS(X, "X.rds")).
  # Read once per replicate by run_replicate.R / local_driver.R,
  # each of which hands it to fitDensityReplicate()'s `X` argument. The
  # replicate-level subsampling (see `rho` below) happens AFTER this read,
  # so every replicate currently reads the full matrix even though it then
  # only uses a fraction of it -- fine for a dataset that fits in memory on
  # one machine, worth knowing if it doesn't.
  dataPath   = "/path/to/your/X.rds",

  # Where each replicate's small result lands: run_replicate.R (or the map
  # step of local_driver.R) writes one replicate_<idx>.rds file
  # per replicate here; run_combine.R reads every file matching
  # replicate_%03d.rds for idx in 1..M back out of this same directory
  # before pooling them. Never read or written by fitDensityReplicate()/
  # combineDensityReplicates() directly -- it's the entry-point scripts'
  # own file-based hand-off between the map and reduce phases.
  resultsDir = "/path/to/your/results",

  # Where the final fitted ensemble is written: run_combine.R (or the
  # reduce step of local_driver.R) saveRDS()s combineDensityReplicates()'s
  # return value here once -- the fitted PHMv2 object, whose `$paramsList`
  # holds the pooled, merged components from all M replicates.
  outPath    = "/path/to/your/phm_fit.rds",

  ## ============================================================
  ## STAGE 1 -- always run, inside fitDensityReplicate() for every
  ## replicate. Chooses a (model, G) by BIC via mclustAutoExtend(), on a
  ## (sub)sample of X.
  ## ============================================================

  # How many ensemble replicates to fit. NOT itself a parameter of
  # fitDensityReplicate() or combineDensityReplicates() -- neither function
  # takes an M argument (fitDensityReplicate() fits exactly one replicate
  # given its own `idx`/`seed`; combineDensityReplicates() infers M as
  # length(replicateResults)). M lives here purely as the shared loop bound
  # the orchestration scripts need: PHM::phmSeeds() uses it to draw M
  # seeds, run_combine.R uses it to know how many replicate_*.rds files to
  # expect, and local_driver.R uses it as its seq_len(M) loop
  # range. (In the original single-process fitPHMEnsemble(), M WAS a real
  # function argument, for the same reason: it was the loop bound inside
  # that one function. Splitting fitting from pooling moved that bound out
  # to here instead.)
  M = 40,

  # Fraction of X each replicate samples WITHOUT replacement, in (0, 1].
  # Used directly inside fitDensityReplicate(): `subIdx <- sample.int(n,
  # round(rho * n), replace = FALSE)` is the very first thing it does with
  # X, before Stage 1 or Stage 2 ever sees the data. rho = 1 draws all n
  # rows without replacement -- a full permutation, statistically identical
  # to using X as-is, since nothing downstream cares about row order.
  # Smaller rho means a faster, less data-hungry replicate, at the cost of
  # each individual replicate seeing less of the structure (M of them
  # pooled together is what's supposed to make up for that).
  rho = 0.25,

  # The initial candidate range for the number of mixture components G,
  # passed straight through to mclustAutoExtend()'s `GInit` argument (which
  # hands it to mclust::mclustBIC() as its own `G=`). This is a STARTING
  # range, not a hard limit -- see GStep/GMax below for what happens if
  # BIC's optimum lands right at the top of this range.
  G = 1:25,

  # Restricts which mclust covariance-model families (e.g. "EII", "VVV")
  # Stage 1's BIC search is allowed to consider; NULL imposes no
  # restriction and lets mclust search its full default set. Passed to
  # mclustAutoExtend()'s `modelNames` argument (and on into every
  # mclust::mclustBIC() call it makes); also validated up front inside
  # fitDensityReplicate() against mclust::mclust.options("emModelNames")
  # before any fitting starts, so a typo here fails fast rather than deep
  # inside a long replicate fit. Whatever model Stage 1 settles on for a
  # replicate is also what that replicate's Stage 2 EM re-fit (if any)
  # uses -- modelNames has no separate effect on Stage 2.
  modelNames = NULL,

  # If mclustAutoExtend()'s BIC-optimal G lands exactly at the current
  # range's upper boundary (suggesting BIC would prefer even more
  # components than were tried), it extends the search upward by GStep
  # components at a time -- re-using every BIC value already computed via
  # mclust::mclustBICupdate() rather than re-fitting from scratch -- until
  # either an interior optimum is found or GMax (next) is hit.
  GStep = 5,

  # The hard cap GStep-driven extension will not go past, regardless of
  # where BIC still wants to go. If a replicate's search hits this cap, a
  # warning-style message is printed (when extendVerbose is TRUE)
  # suggesting GMax be raised -- fitting still completes, capped at GMax
  # components, rather than extending forever.
  GMax = 100,

  # Whether mclustAutoExtend() prints its own progress messages (range
  # extended, extension failed, hit the GMax cap) while it works. Purely
  # about console/log verbosity -- has no effect on what gets fit. Worth
  # turning on (TRUE) the first time you run a new dataset through this,
  # to see whether your G/GStep/GMax choices are actually adequate for it;
  # the double-rings example turns it off once it's been checked once.
  extendVerbose = FALSE,

  # Seeds fitDensityReplicate()'s own data subsampling, and (if stage2Init
  # is set) its Stage 2 random restart(s) -- NOT seedBase directly: each
  # replicate gets its own seed, deterministically derived from seedBase
  # and that replicate's idx by PHM::phmSeeds(). Changing seedBase
  # reproducibly changes every replicate's random draws at once; keeping
  # it fixed is what lets run_replicate.R and run_combine.R (run as
  # separate processes, possibly on separate machines) agree on which seed
  # belongs to which replicate without any file passing a seed around.
  seedBase = 1000,

  ## ============================================================
  ## STAGE 2 -- optional, inside fitDensityReplicate(), after Stage 1. A
  ## fresh EM re-fit of Stage 1's chosen (model, G) from a random start.
  ## Off by default. It matters most when rho is close to 1 and the data has
  ## 2000 rows or fewer: then every replicate fits the same rows from the
  ## same deterministic starting point and the Stage 1 fits come out
  ## identical. (With more than 2000 rows, mclust starts each fit from a
  ## random 2000-row subset, mclust.options("subset"), drawn with the
  ## replicate's own seed, so replicates differ even at rho = 1.)
  ## ============================================================

  # The random-restart method for Stage 2, or NULL (the default) to skip
  # Stage 2 entirely (a replicate's contribution is then just Stage 1's own
  # converged fit, used as-is). One of NULL, "kmeans", "kmeans++",
  # "random", "randomFromData" -- see fitDensityReplicate()'s own roxygen
  # help page (?PHM::fitDensityReplicate) for exactly what each does to initialize the
  # EM re-fit (mclust::me()).
  stage2Init = NULL,

  # How many times Stage 2 retries with a fresh random seed if a given
  # start produces a degenerate fit (e.g. an empty component), before
  # giving up and falling back to the plain Stage 1 fit for that
  # replicate. Has no effect at all when stage2Init is NULL.
  nRestartTries = 5,

  ## ============================================================
  ## REDUCE STEP -- used only inside combineDensityReplicates(), after
  ## every replicate's fitDensityReplicate() result has been collected.
  ## Pools the M replicates' components and fits PHMv2() on the pooled
  ## ensemble. None of these affect Stage 1 or Stage 2 -- they're read only
  ## once, by the one process that runs the reduce step.
  ## ============================================================

  # How much of the pooled ensemble's total probability mass each
  # replicate contributes, before PHMv2() ever sees it. One of:
  #   "equal"      -- the default: every replicate gets exactly 1/M,
  #                   regardless of fit quality or how many components it has.
  #   "likelihood" -- a softmax over the M replicates' own log-likelihoods
  #                   (log-sum-exp stabilized).
  #   "bic"        -- the same softmax, but over BIC/2 instead of raw
  #                   log-likelihood, so a replicate doesn't get more say
  #                   just because mclustAutoExtend() let it settle on a
  #                   larger G.
  # "likelihood"/"bic" are only a fair comparison across replicates whose
  # (sub)samples were the same size -- true by construction here, since
  # every replicate uses the same rho against the same n.
  replicateWeighting = "equal",

  # Monte Carlo sample size and batch size for PHMv2()'s pairwise
  # distinguishability-criterion (Delta P_mc) estimation between pooled
  # components -- passed straight through to PHMv2()'s own `mcSamples`/
  # `batchSize` arguments. Bigger mcSamples means a more precise merge
  # decision at the cost of more compute in the one process that runs the
  # reduce step; batchSize only controls how that work is chunked in
  # memory, not the estimate itself.
  mcSamples = 2e5,
  batchSize = 1e5,

  # Cores PHMv2() uses internally for that same Monte Carlo computation --
  # passed to its `numCores` argument. This is the ONLY place numCores
  # does anything: it has no effect on fitDensityReplicate() or Stage 1/
  # Stage 2 at all, and does not control how many replicates run at once
  # (that's MAX_CONCURRENT in single_machine/launch.sh /local_driver.R, or
  # cpus-per-task in slurm_array/*.sbatch -- a different kind of
  # parallelism, at a different stage of the pipeline).
  numCores = 2,

  ## ============================================================
  ## BOOKKEEPING -- read by both fitDensityReplicate() and
  ## combineDensityReplicates(), but only for the text they print to the
  ## console/log; never affects what gets fit.
  ## ============================================================

  # A free-text tag stamped onto every progress message either function
  # prints (e.g. "[cluster_run] PHM replicate 3: G=7 ..."). Purely
  # cosmetic -- useful for telling one run's log lines apart from
  # another's when several runs' output ends up interleaved or
  # concatenated somewhere.
  label = "cluster_run"
)
