# PHM Explorer

PHM Explorer is a single self-contained HTML page for browsing a `PHMv2` fit: a dendrogram, the ΔP_mc heatmap aligned under it, and (when data are given) a scatter of the points shaded by posterior probability of the hovered branch. It needs no server. A control collapses the full tree live to any number of leaves, so one export serves every resolution.

The skill bundles the tool in `scripts/explorer/` (copies of `PHM_Explorer/r/` and `PHM_Explorer/src/` from the PHM repository):

```
scripts/explorer/r/export_phm_json.R    PHMv2 fit -> JSON export
scripts/explorer/r/view_in_explorer.R   one command: fit.rds -> export + standalone page, then opens it
scripts/explorer/src/template.html      the app (styles, markup, logic) with a __DATA_JSON__ placeholder
scripts/explorer/src/build.py           template + export -> standalone page
```

Requirements: the PHM package, jsonlite, and python3 (for `build.py`).

## One command

```sh
Rscript scripts/explorer/r/view_in_explorer.R <fit.rds> [X.rds] [--scatter=<scatter.rds>]
```

| Argument | Meaning |
|---|---|
| `fit.rds` | A saved `PHMv2` object, e.g. `fitPHMEnsemble()` output or the ensemble pipeline's `outPath`. |
| `X.rds` | Optional. An N × D matrix in the fit's own D dimensions: the fitting data, a subsample of it, or new data. Used to compute posteriors, and without `--scatter=` its first two columns `X[, 1:2]` are the scatter layout. Without it you get the tree and heatmap only. |
| `--scatter=scatter.rds` | Optional (the flag was `--plot=` before Nov 25, 2026). A data frame with numeric `x`, `y` (the scatter layout, replacing `X[, 1:2]`; use it whenever D > 2), optional `annotation` (categorical; adds a second, static panel colored by category; `NA` = unlabeled) and optional `id` (checked against `rownames(X)`). Row i must be observation i of X. |

It writes `<fit>_export.json` and `<fit>_explorer.html` next to the fit and opens the page (`open` on macOS, `xdg-open` on Linux; otherwise it prints the path). It checks for PHM, jsonlite, python3 and its own helper files first and lists every problem at once. It stops, using none of the scatter data, if its rows do not match X in count or `id`. Keep the script inside `scripts/explorer/r/`: it locates `export_phm_json.R` and `../src/build.py` relative to itself.

When running headless (an agent with no browser), skip the opening step by calling the two pieces directly:

```sh
Rscript -e 'source("scripts/explorer/r/export_phm_json.R"); export_phm_json(readRDS("phm_fit.rds"), X = readRDS("X.rds"), out_path = "phm_fit_export.json", label = "my run")'
python3 scripts/explorer/src/build.py --data phm_fit_export.json --out phm_fit_explorer.html
```

then hand the user the `.html` file. The page also has a **load a PHM export...** button and accepts a dropped `.json`, so one built page can view any export.

## From R

```r
library(PHM)
source("scripts/explorer/r/export_phm_json.R")

export_phm_json(
  fit,                         # PHMv2 result (required)
  X          = X,              # optional: data in the fit's own dimension, for posteriors (fitting data, a subsample, or new data)
  xy         = embedding,      # optional: n x 2 layout (UMAP, t-SNE); defaults to X[, 1:2]
  z          = NULL,           # optional: n x K posteriors, instead of X
  annotation = cell_type,      # optional: length-n categories for the static panel
  annotation_name = "cell type",
  out_path   = "fit_export.json",
  label      = "PBMC, M = 20, rho = 0.25",   # shown in the page footer
  digits     = 4,              # rounding of coordinates and heights
  prob_eps   = 0.005           # posteriors below this are dropped from the sparse export
)
```

The scatter appears only when the export has both a layout (`xy`, or `X[, 1:2]` when `xy` is omitted) and posteriors (`z`, or `X` in the fit's dimension). When X has more than two columns and `xy` is omitted, the exporter says it is plotting `X[, 1:2]`; for PCs that is PC1 vs PC2, often a poor picture of the clusters, so prefer a UMAP or t-SNE `xy`.

## What the exporter can and cannot take

- **Alpha-scaled fits only.** Heights are computed on the `"d*"` scale, so a fit made with `scaling = "unscaled"` or `"average"` fails with "d scaling only makes sense with alpha criterion". Re-run `PHMv2(paramsList = fit$paramsList, deltaPmc = fit$deltaPmc, scaling = "alpha")`; it reuses the ΔP_mc matrix, so no Monte Carlo.
- **Single-Gaussian leaves for automatic posteriors.** With `X`, posteriors are computed assuming each `paramsList` element is one Gaussian, which holds for ensembles, `globalDensityEstimation()`, `constructPmcParamsMclust()` and `localizedDensityEstimation(preserveClusters = FALSE)`. For partition fits whose clusters are mixtures (`constructPmcParamsWeightedPartition()`, `preserveClusters = TRUE`), passing `X` fails with "X has 2 column(s) but the fit has 6 dimension(s)". Pass `z = computePosteriorProbMatrix(fit$paramsList, X)` and `xy` instead.
- **`X` must have the fit's dimension**, never truncated to 2 columns; its rows need not be the fitting data (a subsample or new data works). For D > 2 pass the full X and, for a useful picture, a separate 2-D `xy`.
- **Size.** Exports scale to thousands of components and tens of thousands of points (posteriors stored sparsely). For millions of points, export a subsample, taking the same rows from X, `xy` and `annotation`.
- `PHM()` (v1) results: convert first with `convertToPHMv2(phm)`; they then need `paramsList` and alpha scaling, so in practice re-run `PHMv2()` on the v1 `params`.

## Reading the display

- Hover a leaf or branch in the dendrogram: the whole subtree highlights (gold), the matching block of the heatmap lights up, and the scatter shades every point by its posterior probability of belonging to that branch.
- Heatmap cell (i, j) is colored by the height at which leaves i and j first merge: blue = merge late (well separated), red = merge early (weakly distinguishable).
- "show top N of K components" collapses the tree to N leaves by absorbing the finest merges first, exactly as cutting the PHM tree at N clusters. Start big ensembles at 10–30 leaves.
- The branch-detail panel lists the components in the hovered branch, its d*-scaled merge height and, with a scatter, how many points it claims.
- No ground-truth labels are used anywhere; the only labels shown are the optional `annotation`, in its own static panel.

## Troubleshooting

| Message or symptom | Cause and fix |
|---|---|
| "That file doesn't look like a PHM Explorer export" | The JSON lacks top-level `meta`, `nodes` or `xlabOrder`; it was not made by `export_phm_json()` or was edited. Re-export. |
| No scatter panel | The export lacks posteriors (no `X` or `z`) or a layout. Give `X`, or `xy` plus `z`. |
| Scatter looks like an uninformative blob | `xy` was omitted for D > 2 data, so `X[, 1:2]` was plotted. Pass a UMAP/t-SNE `xy` (or `--scatter=`). |
| "X has N column(s) but the fit has M dimension(s)" | X is not the data the fit used, or the leaves are mixtures. See above. |
| "d scaling only makes sense with alpha criterion" | Non-alpha fit; re-run with `scaling = "alpha"` reusing `deltaPmc`. |
| "scatter data has N rows but X has M" / id mismatch | Scatter data not aligned with X; nothing from it was used. |
| "usage: ... [--scatter=<plot.rds>]" when passing `--plot=` | The flag was renamed; use `--scatter=`. |
| `python3 not found` | Install Python 3, or build elsewhere: the `.json` alone can be loaded into any already-built page. |
