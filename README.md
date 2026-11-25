# PHM: software package for multi-resolution clustering

**Version 1.0**

This package implements the PHM algorithm for multiresolution and consensus clustering analysis, as well as its associated visualizations. 
The method and computational procedures are described in [Turfah and Wen (2025)](https://www.biorxiv.org/content/10.1101/2025.06.13.659529v1) and [Turfah and Wen (2026)](https://arxiv.org/abs/2611.32945).

For questions/commens regarding the package, please contact Ali Turfah (aturfah at umich dot edu) and William Wen (xwen at umich dot edu).

Please see our [user guide](https://aturfah.github.io/PHM) for a detailed guide and use case examples.

## License

Software distributed under the terms of the GNU General Public License version 3.0. See [LICENSE](LICENSE) for more details.

## Installation and Loading

To install please run `devtools::install_github("aturfah/PHM")`

To load in an `R` session, run the following:
```
library(PHM)
library(mclust)
```

## Running a PHM ensemble

`fitPHMEnsemble()` fits an ensemble of Gaussian mixture models on random subsamples of the data and builds the PHM hierarchy on their pooled components. Its two halves, `fitDensityReplicate()` (one replicate) and `combineDensityReplicates()` (pooling plus `PHMv2()`), can also be run separately, so the replicates can be fit in parallel.

The package installs command-line scripts for this in `inst/ensemble/` (find them with `system.file("ensemble", package = "PHM")`): a serial runner, a single-machine parallel runner, a Slurm job-array runner, and a worked example on the `double_ring` dataset. See [`inst/ensemble/README.md`](inst/ensemble/README.md) and the [user guide page](https://aturfah.github.io/PHM/density_estimation/ensemble.html).

