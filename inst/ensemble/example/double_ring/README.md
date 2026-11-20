# example/double_ring/

A complete PHM ensemble run on the package's `double_ring` dataset, three
ways: in serial, in parallel on one machine, and on a Slurm cluster. All three
use the same `config.R` (M = 10 replicates, rho = 0.5) and give the same
pooled components.

| File | What it is |
|---|---|
| `config.R` | Settings for this example: M = 10 and rho = 0.5, with the defaults for everything else (no Stage 2 refit, equal replicate weights). |

## The data

`PHM::double_ring`: 5000 points on two concentric rings centered at the
origin, with points spread uniformly over each ring's area:

- inner ring: radius between 2 and 3, 2500 points
- outer ring: radius between 6 and 7, 2500 points

The matrix has no ring labels. A uniform ring is not Gaussian, so the mixture
model has to tile each ring with many components (20 to 24 per replicate
here). PHM's job is to merge those pieces back into two rings.

## Before you start

Install PHM (`remotes::install_github("aturfah/PHM")`). Then set two shell
variables and work in a copy of this folder, so that results never land inside
the package or repo:

```sh
ENS=$(Rscript -e 'cat(system.file("ensemble", package = "PHM"))')
EX=$ENS/example/double_ring
mkdir -p ~/phm_runs
```

Each method below copies the example folder to its own working folder, writes
the data there, and runs from there, because `config.R` uses paths relative to
the current directory.

Run times below were measured on a small 2-core Linux machine. One replicate
takes about 45 seconds of CPU time, so expect less on a recent laptop.

## 1. In serial

All 10 replicates run one after another in a single R session.

```sh
cp -r $EX ~/phm_runs/dr_serial
cd ~/phm_runs/dr_serial
Rscript -e 'saveRDS(PHM::double_ring, "double_ring_X.rds")'
Rscript $ENS/run_serial.R
```

The log prints one line per replicate:

```
[run_serial] fitting 10 replicates one after another on a 5000 x 2 matrix
[double_ring] PHM replicate 1: G=22 (model=EEV), loglik=-10281.76 (no Stage 2)
[double_ring] PHM replicate 2: G=22 (model=EEV), loglik=-10299.72 (no Stage 2)
...
[double_ring] PHM replicate 10: G=23 (model=EEV), loglik=-10275.26 (no Stage 2)
[double_ring] PHM: fitting PHMv2 on 223 pooled components from 10 replicates...
[run_serial] wrote double_ring_phm_fit.rds (223 pooled components from 10 replicates)
```

Had a replicate's best `G` been 25, the top of the starting range, the log
would also show a line saying the search was widened to `26:30`.

Output: `double_ring_phm_fit.rds`. No per-replicate files are written in this
mode. Time: about 7.5 minutes.

## 2. On one machine, 5 replicates at a time

Each replicate runs as its own background process, at most 5 at once.

```sh
cp -r $EX ~/phm_runs/dr_local
cd ~/phm_runs/dr_local
Rscript -e 'saveRDS(PHM::double_ring, "double_ring_X.rds")'
MAX_CONCURRENT=5 bash $ENS/parallel_run/single_machine/launch.sh
```

`launch.sh` first checks the setup, then starts replicates 1 to 5, waits for
all five to finish, starts 6 to 10, and finally pools them:

```
preflight checks passed (PHM 1.1.0, dataPath OK)
running 10 replicates, up to 5 at a time (set MAX_CONCURRENT to change this)
  launched replicate 1 (pid 300) -> /home/you/phm_runs/dr_local/logs/replicate_1.log
  ...
  launched replicate 10 (pid 375) -> /home/you/phm_runs/dr_local/logs/replicate_10.log
all 10 replicates finished -- running the combine step
[double_ring] PHM: fitting PHMv2 on 223 pooled components from 10 replicates...
[run_combine] wrote double_ring_phm_fit.rds (223 pooled components from 10 replicates)
```

Output: `results/replicate_001.rds` to `replicate_010.rds`, one log per
replicate in `logs/`, and `double_ring_phm_fit.rds`. Time: about 4.5 minutes
on 2 cores; with 5 or more free cores, roughly the time of two replicates plus
the pooling step.

The one-session alternative gives the same result:

```sh
MAX_CONCURRENT=5 Rscript $ENS/parallel_run/single_machine/local_driver.R
```

To test a single replicate before a full run:

```sh
Rscript $ENS/parallel_run/run_replicate.R 1
```

## 3. On a Slurm cluster, all 10 at once

Each replicate is one task of a 10-task job array; Slurm runs as many of them
at the same time as the cluster allows. A second job pools the results once
every task has succeeded.

**One-time setup.** PHM must be installed where the compute nodes' R can load
it, and the working folder must be on a filesystem the nodes can read, such as
your home directory. In both `$ENS/parallel_run/slurm_array/replicate.sbatch`
and `combine.sbatch`, replace the commented `# module load R` line with
whatever makes `Rscript` available on your cluster (edit copies of the two
files if the installed package folder is read-only for you). The default
resources (1 CPU, 4 GB, 2 hours per replicate) are plenty for this example.

**Submit**, from a login node with R loaded, because `submit.sh` reads `M` from
`config.R`:

```sh
module load R            # or however your cluster provides R
cp -r $EX ~/phm_runs/dr_slurm
cd ~/phm_runs/dr_slurm
Rscript -e 'saveRDS(PHM::double_ring, "double_ring_X.rds")'
bash $ENS/parallel_run/slurm_array/submit.sh
```

`submit.sh` sets the array to `1-10` from `M` in `config.R`, submits it, then
submits the combine job with `--dependency=afterok`, prints both job ids and
exits:

```
map array job:     123456  (10 tasks, --array=1-10; see /home/you/phm_runs/dr_slurm/logs/replicate_*.log)
combine job:        123457  (queued to run after 123456; see /home/you/phm_runs/dr_slurm/logs/combine.log)
final fit will land at the outPath set in /home/you/phm_runs/dr_slurm/config.R once 123457 finishes
```

To cap how many tasks run at once, set `MAX_CONCURRENT` (for example
`MAX_CONCURRENT=5 bash $ENS/parallel_run/slurm_array/submit.sh` submits
`--array=1-10%5`). Leave it unset to submit all 10 at once.

**Watch it** with `squeue -u $USER`. Each task writes
`logs/replicate_<task>.log` and `results/replicate_<task>.rds`; the combine
job writes `logs/combine.log` and `double_ring_phm_fit.rds`. With all 10 tasks
running together, the whole run takes about as long as one replicate plus the
pooling step.

**If a task fails,** the combine job stays pending (`DependencyNeverSatisfied`).
Cancel it, fix the cause, and resubmit only the failed tasks, for example
task 4:

```sh
scancel 123457
SD=$ENS/parallel_run/slurm_array
sbatch --array=4 --output=logs/replicate_%a.log $SD/replicate.sbatch $SD $PWD/config.R
# after it succeeds:
sbatch --output=logs/combine.log $SD/combine.sbatch $SD $PWD/config.R
```

## 4. Check the result

From the working folder of any of the three runs:

```r
library(PHM)
fit <- readRDS("double_ring_phm_fit.rds")
X   <- readRDS("double_ring_X.rds")

length(fit$paramsList)                        # 223 pooled components
fit$mergeDeltaPmc[2:5]                        # heights of the top merges

post <- computePosteriorProbMatrix(fit$paramsList, X)
cl2  <- apply(recoverPosteriorv2(fit, 2, post), 1, which.max)
ring <- ifelse(sqrt(rowSums(X^2)) < 4.5, "inner", "outer")
table(cluster = cl2, ring)
```

Cut at two clusters, the tree recovers the two rings exactly (which ring gets
cluster number 1 can differ between runs):

```
       ring
cluster inner outer
      1     0  2500
      2  2500     0
```

The final merge, which joins the two rings, has a ΔPmc of about 1e-7, while
the merges just below it are about 0.005 to 0.03, so two clusters is a clear
cut. Below that level the outer ring (154 of the 223 components, against 69
for the inner one) splits into arcs.

The three methods use the same seeds, so `fit$paramsList` is identical across
them. `fit$deltaPmc` and the merge heights differ very slightly between runs,
because `PHMv2()` estimates ΔPmc by Monte Carlo; the two-ring cut is the same.

To view the fit in PHM Explorer (from a clone of the PHM repo):

```sh
Rscript /path/to/PHM/PHM_Explorer/r/view_in_explorer.R double_ring_phm_fit.rds double_ring_X.rds
```

## Why these settings

- **rho = 0.5** means each replicate fits its own random 2500 of the 5000
  points. Different halves give different fits (here G ranges from 20 to 24),
  and across 10 replicates every point is used about 5 times. With rho = 1
  every replicate would see the same rows; on data with 2000 rows or fewer,
  that gives the same fit in every replicate unless `stage2Init` is set.
  (Above 2000 rows, mclust starts each fit from a random 2000-row subset,
  `mclust.options("subset")`, so even rho = 1 gives some variation.)
- **`stage2Init = NULL` and `replicateWeighting = "equal"`** are the defaults:
  no EM refit, and every replicate contributes 1/10 of the pooled mixture.
- **`G = 1:25` with `GStep = 5`.** About 20 to 25 components per replicate are
  needed to tile the two rings; if BIC wants more, the range widens.
