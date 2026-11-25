#!/usr/bin/env Rscript
## view_in_explorer.R -- converts a saved PHMv2() fit (for example the .rds
## that fitPHMEnsemble() or the scripts in PHM's inst/ensemble/ write to
## CFG$outPath) into a PHM Explorer JSON export, builds a standalone HTML
## page with that export embedded, and opens it -- one command instead of
## exporting and loading the file by hand.
##
## Lives in PHM_Explorer/r/ and finds export_phm_json.R (next to it) and
## src/build.py (one level up) from its own location, so it can be run from
## any directory.
##
## Usage:
##   Rscript view_in_explorer.R <fit.rds> [X.rds] [--scatter=<plot.rds>]
##
## <fit.rds>  Required. A saved PHMv2() result, e.g. from
##            combineDensityReplicates() (whatever CFG$outPath pointed at).
## [X.rds]    Optional. The same data matrix that fit was run on (CFG$dataPath).
##            Adds the scatter panel. Without it you still get a complete,
##            valid export -- just the dendrogram/heatmap, no scatter (see
##            README.md: "A tree-only export is entirely valid."). Without
##            --scatter=, the first two columns of X (X[, 1:2]) are plotted;
##            if X has more than 2 columns, use --scatter= below to plot a
##            better 2D layout (UMAP, t-SNE, PCA, ...) instead.
##
## [--scatter=]  Optional. A data frame for 2D plotting, used only by the Explorer
##            (see README.md, "Launch explorer from command line"). Columns:
##              x, y         required, numeric -- the 2D coordinates to plot;
##                           they replace X as the scatter layout, so this is
##                           how to get a scatter for data with D > 2
##              annotation   optional, categorical -- adds a static scatter
##                           panel, one color per category (NA = unlabeled)
##              id           optional -- if present AND X has rownames, the two
##                           must match exactly, in order
##            Row i must be the same observation as row i of X. That is the
##            caller's responsibility; this script stops with an error only on
##            clear evidence of misalignment (different row counts, or ids that
##            disagree with X's rownames) and then uses NONE of the scatter data.
##            The scatter data never touches the fit.
##
## Output: <fit.rds>'s own directory gets two new files next to it, named
## after <fit.rds> itself -- e.g. phm_fit.rds -> phm_fit_export.json and
## phm_fit_explorer.html -- so results from different runs never collide
## and nothing is written inside PHM_Explorer itself (in particular,
## index.html -- the app's own default-on-open page -- is never touched).

## PHM_Explorer/ is the parent of the r/ folder this script lives in.
PHM_EXPLORER_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) != 1) stop("run this with Rscript: Rscript view_in_explorer.R <fit.rds> [X.rds] [--scatter=<plot.rds>]")
  dirname(dirname(normalizePath(f)))
})

usage <- "usage: Rscript view_in_explorer.R <fit.rds> [X.rds] [--scatter=<plot.rds>]"
cliArgs <- commandArgs(trailingOnly = TRUE)
scatterFlag <- grepl("^--scatter=", cliArgs)
scatterPath <- if (any(scatterFlag)) sub("^--scatter=", "", cliArgs[scatterFlag][1]) else NA
if (sum(scatterFlag) > 1) stop(usage)
posArgs <- cliArgs[!scatterFlag]
if (any(grepl("^--", posArgs)) || length(posArgs) < 1 || length(posArgs) > 2 ||
    (!is.na(scatterPath) && !nzchar(scatterPath))) {
  stop(usage)
}
rdsPath <- posArgs[1]
xPath <- if (length(posArgs) >= 2) posArgs[2] else NA

## --- preflight: the same "every problem at once" approach as
## single_machine/check_requirements.R, rather than failing on the first
## thing and making you fix-and-rerun repeatedly ---
problems <- character(0)

if (!file.exists(rdsPath)) {
  problems <- c(problems, sprintf("fit file '%s' not found", rdsPath))
}
if (!is.na(xPath) && !file.exists(xPath)) {
  problems <- c(problems, sprintf("data file '%s' not found", xPath))
}

if (!is.na(scatterPath) && !file.exists(scatterPath)) {
  problems <- c(problems, sprintf("scatter data file '%s' not found", scatterPath))
}

exporterPath <- file.path(PHM_EXPLORER_DIR, "r", "export_phm_json.R")
buildPyPath <- file.path(PHM_EXPLORER_DIR, "src", "build.py")
if (!dir.exists(PHM_EXPLORER_DIR)) {
  problems <- c(problems, sprintf(
    "PHM_EXPLORER_DIR ('%s') does not exist -- keep this script in PHM_Explorer/r/",
    PHM_EXPLORER_DIR))
} else {
  if (!file.exists(exporterPath)) problems <- c(problems, sprintf("missing %s -- keep this script in PHM_Explorer/r/", exporterPath))
  if (!file.exists(buildPyPath)) problems <- c(problems, sprintf("missing %s -- keep this script in PHM_Explorer/r/", buildPyPath))
}

if (!requireNamespace("PHM", quietly = TRUE)) {
  problems <- c(problems, "R package 'PHM' is not installed -- export_phm_json() needs it: install.packages(c(\"remotes\", \"jsonlite\")); remotes::install_github(\"aturfah/PHM\")")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  problems <- c(problems, "R package 'jsonlite' is not installed: install.packages(\"jsonlite\")")
}
if (nchar(Sys.which("python3")) == 0) {
  problems <- c(problems, "python3 not found on PATH -- needed to run PHM_Explorer/src/build.py")
}

if (length(problems) > 0) {
  cat("view_in_explorer.R: cannot proceed:\n", file = stderr())
  for (p in problems) cat(sprintf("  - %s\n", p), file = stderr())
  quit(status = 1, save = "no")
}

## --- export the fit to JSON ---
source(exporterPath)

fit <- readRDS(rdsPath)
X <- if (!is.na(xPath)) readRDS(xPath) else NULL

## --- optional scatter data: validate before anything is written; on any
## problem, stop and use none of it ---
xy <- NULL
annotation <- NULL
if (!is.na(scatterPath)) {
  scatterDat <- readRDS(scatterPath)
  if (!is.data.frame(scatterDat) && !is.matrix(scatterDat)) {
    stop(sprintf("scatter data '%s' must be a data frame (or matrix) with columns x and y", scatterPath))
  }
  scatterDat <- as.data.frame(scatterDat, stringsAsFactors = FALSE)
  if (!all(c("x", "y") %in% names(scatterDat))) {
    stop(sprintf("scatter data '%s' must have columns named x and y (found: %s)",
                 scatterPath, paste(names(scatterDat), collapse = ", ")))
  }
  if (!is.numeric(scatterDat$x) || !is.numeric(scatterDat$y) ||
      !all(is.finite(scatterDat$x)) || !all(is.finite(scatterDat$y))) {
    stop("scatter data columns x and y must be numeric with no NA/NaN/Inf")
  }
  if (!is.null(X)) {
    if (nrow(scatterDat) != nrow(X)) {
      stop(sprintf(paste0("scatter data has %d rows but X has %d -- they must be the same observations ",
                          "in the same order. Scatter data NOT used."), nrow(scatterDat), nrow(X)))
    }
    if ("id" %in% names(scatterDat) && !is.null(rownames(X)) &&
        !identical(as.character(scatterDat$id), rownames(X))) {
      stop("scatter data's id column does not match X's rownames (values or order). Scatter data NOT used.")
    }
  } else {
    cat("note: no X given, so scatter data alignment with the data cannot be checked (and posterior shading is unavailable)\n")
  }
  xy <- cbind(x = scatterDat$x, y = scatterDat$y)
  if ("annotation" %in% names(scatterDat)) annotation <- scatterDat$annotation
}

rdsDir <- dirname(normalizePath(rdsPath))
baseName <- sub("\\.rds$", "", basename(rdsPath))
jsonPath <- file.path(rdsDir, sprintf("%s_export.json", baseName))
htmlPath <- file.path(rdsDir, sprintf("%s_explorer.html", baseName))

export_phm_json(fit, X = X, xy = xy, annotation = annotation,
                out_path = jsonPath, label = baseName)
cat(sprintf("wrote %s\n", jsonPath))

## --- build a standalone page with that export as its default dataset ---
buildResult <- system2("python3", c(shQuote(buildPyPath), "--data", shQuote(jsonPath), "--out", shQuote(htmlPath)),
                        stdout = TRUE, stderr = TRUE)
buildStatus <- attr(buildResult, "status")
if (!is.null(buildStatus) && buildStatus != 0) {
  cat(paste(buildResult, collapse = "\n"), "\n", file = stderr())
  stop("build.py failed -- see output above")
}
cat(sprintf("wrote %s\n", htmlPath))

## --- open it ---
sysName <- Sys.info()[["sysname"]]
opened <- FALSE
if (sysName == "Darwin") {
  opened <- (system2("open", shQuote(htmlPath)) == 0)
} else if (sysName == "Linux" && nchar(Sys.which("xdg-open")) > 0) {
  opened <- (system2("xdg-open", shQuote(htmlPath)) == 0)
}

if (opened) {
  cat(sprintf("opened %s in your default browser\n", htmlPath))
} else {
  cat(sprintf("open this file in a browser to view it: %s\n", htmlPath))
}
