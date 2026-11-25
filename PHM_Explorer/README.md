# PHM Explorer

PHM explorer is an interactive dendrogram $+$ heatmap (and optional 2D scatter) for exploring [PHM](https://github.com/aturfah/PHM) clustering results.

## Quick start

Open [`index.html`](index.html) directly in a browser to try it with a preload example (the two-moons benchmark, $n = 1000, K = 10$). 

To explore your own PHM result, export it to JSON from R. Next, load the resulting `.json` file by either clicking on **load a PHM export...** button in the page heade or directly dragging and dropping the JSON file anywhere on the page. 

Please see [`docs/guide.html`](docs/guide.html) for the full walkthrough.

## Directory layout

```
index.html            standalone app, built with the two-moons example as its default
src/
  template.html        source template (styles + markup + app logic), with a
                        __DATA_JSON__ placeholder in place of an embedded dataset
  build.py              regenerates index.html from template.html + a chosen export
r/
  export_phm_json.R    generic exporter: any PHMv2() result -> the JSON the app reads
  view_in_explorer.R   one command: saved fit .rds -> export + standalone page, opened
examples/
  two_moons_export.json            K = 10, n = 1000, with scatter data
  concentric_rings_export.json     K = 20, n = 1500, three nested rings, with scatter data
  compound_11700_export.json       K = 20, n = 11700, with scatter data
  compound_11700_no_scatter.json   same fit, tree only (no X/scatter panel)
docs/
  guide.html            field guide: install, export, load, read the display
```

## Export PHM result and plotting data to JSON file 


A sample R code is shown below

```r
install.packages(c("remotes", "jsonlite"))
remotes::install_github("aturfah/PHM")

library(PHM)
source("r/export_phm_json.R")

fit <- PHMv2(X, ...)   # your fit, X is 2D here

export_phm_json(
  fit,                # required 
  X        = X,       # optional for the scatter panel
  out_path = "my_result_export.json", # required
  label    = "my analysis"  # optional
)
```
The output file `my_result_export.json` is ready to be loaded into `index.html`.

### Optional information for scatter panel

The following information is for plotting the optional scatter panel, which helps link dendrogram and heatmap to the original clustering data. 

+ `X`: an $n \times d$ matrix of $d$-dimensional data points for clustering. It is used for computing the posterior probabilities. It does not have to be the original clustering data used for fitting PHM - a subsample or an indepdent sample is also suitable.

+ `xy`: an $n \times 2$ matrix containing the 2D coordinates for plotting the data points. Commonly used 2D embeddings/coordinates can be generated from principal component analysis (PCA), UMAP, and t-SNE. **If `xy` is not specified, the first 2 columns of `X` matrix is used**.

+ `annotation`: a vector of reference labels for the plotting data. If specified, a second,  **static** scatter panel is added for reference, with one distinct color per category.    

Each row of `X`, `xy` and `annotation` must represent the same data point. The exporter only checks sizes and stops with an error if they disagree. Users are responsible for ensuring the consistent orders across data matrices. 

An example to plot scatter points using the UMAP embedding:
```r
fit <- PHMv2(X, ...)           
embedding <- umap::umap(X)$layout  
export_phm_json(
  fit,
  X        = X,          
  xy       = embedding,  
  out_path = "my_result_export.json"
)
```

An example using reference ground truth labels:
```r
export_phm_json(
  fit,
  X          = X,
  xy         = embedding,                # n x 2, row i = row i of X
  annotation = cell_type,                # length-n factor/character; NA = "unlabeled"
  annotation_name = "cell type",         # optional panel title
  out_path   = "my_result_export.json"
)
```


## Launch explorer from command line

Run `r/view_in_explorer.R` to launch the PHM explorer directly from the command line:

```
Rscript PHM_Explorer/r/view_in_explorer.R <fit.rds> [X.rds] [--scatter=<plot.rds>]
```

| Argument | Meaning |
|---|---|
| `<fit.rds>` | Required. A saved `PHMv2()` result |
| `X.rds` | Optional. The data matrix for the scatter panel |
| `--scatter=<plot.rds>` | Optional. A data frame of 2D coordinates for the scatter panel|

The script requires the `PHM` package, `jsonlite` and `python3` to run. It generates a html file, which can be opened with with `open` on macOS or `xdg-open` on Linux.



## Rebuilding `index.html` with a different default dataset

`index.html` is pre-generated. To create a different default example, run

```
python3 src/build.py --data examples/compound_11700_export.json --out index.html
```
