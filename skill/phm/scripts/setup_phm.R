#!/usr/bin/env Rscript
# Check for, and install, everything needed to run the PHM R package.
#
#   Rscript setup_phm.R                 check, then install whatever is missing (CRAN deps + PHM from GitHub)
#   Rscript setup_phm.R --check-only    report status only; install nothing, create nothing
#   Rscript setup_phm.R --optional      also cover the optional packages used in the skill's recipes
#   Rscript setup_phm.R --explorer      also cover PHM Explorer (jsonlite; checks for python3)
#   Rscript setup_phm.R --update        reinstall PHM from GitHub even if it is already installed
#   Rscript setup_phm.R --lib=/path     install into this library (created if needed)
#   Rscript setup_phm.R --no-verify     skip the final end-to-end smoke run
#
# Exit status is 0 when everything required is present (and the smoke run passes), 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat(paste(readLines(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), n = 13)[2:12],
            collapse = "\n"), "\n")
  quit(status = 0)
}
check_only <- "--check-only" %in% args
want_optional <- "--optional" %in% args
want_explorer <- "--explorer" %in% args
update_phm <- "--update" %in% args
verify <- !("--no-verify" %in% args)
lib_arg <- sub("^--lib=", "", grep("^--lib=", args, value = TRUE))
unknown <- setdiff(args, c("--check-only", "--optional", "--explorer", "--update", "--no-verify", "-h", "--help", grep("^--lib=", args, value = TRUE)))
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = " "), " (try --help)")

## ---- what PHM needs (from its DESCRIPTION: Depends R >= 3.5, Imports below) ----------------
required <- c(
  mvtnorm      = "Gaussian densities and sampling",
  abind        = "stacking covariance arrays",
  mclust       = "GMM fitting (attach it yourself only when your own code calls Mclust())",
  cubature     = "exact integration (monteCarlo = FALSE / mc = FALSE)",
  dplyr        = "plot data preparation",
  tidyr        = "plot data preparation",
  ggplot2      = "all plots",
  RColorBrewer = "default cluster colors",
  ggtext       = "colored axis labels",
  magrittr     = "pipe",
  ggnewscale   = "plots with several fill scales",
  stringr      = "text formatting",
  tibble       = "plot data preparation"
)
optional <- c(
  ggpubr  = "ggarrange() to combine dendrogram + heatmap, as in the skill's recipes",
  ppclust = "fuzzy c-means memberships as weights (weighted-density tutorial)"
)
explorer <- c(
  jsonlite = "writing PHM Explorer exports (scripts/explorer)"
)
github_repo <- "aturfah/PHM"

## ---- helpers --------------------------------------------------------------------------------
ok <- function(x) cat(sprintf("  [ok]      %s\n", x))
miss <- function(x) cat(sprintf("  [MISSING] %s\n", x))
info <- function(x) cat(sprintf("  %s\n", x))
installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)
ver <- function(pkg) as.character(packageVersion(pkg))

cat("== PHM setup check ==\n")
cat("R:", R.version.string, "\n")
if (getRversion() < "3.5.0") stop("PHM needs R >= 3.5.0")

## ---- library location ---------------------------------------------------------------------
if (length(lib_arg)) {
  lib <- normalizePath(lib_arg, mustWork = FALSE)
  if (!dir.exists(lib) && !check_only) dir.create(lib, recursive = TRUE)
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
} else {
  writable <- Filter(function(p) dir.exists(p) && file.access(p, 2) == 0, .libPaths())
  if (length(writable)) {
    lib <- writable[1]
  } else {
    lib <- Sys.getenv("R_LIBS_USER")
    if (!nzchar(lib)) lib <- file.path(Sys.getenv("HOME"), "R", "library")
    if (!check_only) {
      dir.create(lib, recursive = TRUE, showWarnings = FALSE)
      .libPaths(c(lib, .libPaths()))
    }
  }
}
cat("Install library:", lib, if (check_only) "(check-only: nothing will be installed)" else "", "\n")

## CRAN mirror
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

## ---- report status --------------------------------------------------------------------------
report <- function(pkgs, label) {
  cat("\n", label, "\n", sep = "")
  status <- vapply(names(pkgs), installed, logical(1))
  for (p in names(pkgs)) {
    if (status[[p]]) ok(sprintf("%-13s %s", p, ver(p))) else miss(sprintf("%-13s %s", p, pkgs[[p]]))
  }
  names(pkgs)[!status]
}
if (!installed("parallel")) stop("Base package 'parallel' is missing; this R installation is incomplete.")
missing_req <- report(required, "Required CRAN packages:")
missing_opt <- if (want_optional) report(optional, "Optional CRAN packages (--optional):") else character()
if (want_explorer) {
  missing_opt <- c(missing_opt, report(explorer, "PHM Explorer (--explorer):"))
  if (nzchar(Sys.which("python3"))) ok("python3       found (builds the Explorer page)") else
    miss("python3       needed by scripts/explorer/src/build.py; install Python 3")
}

cat("\nPHM itself:\n")
phm_installed <- installed("PHM")
phm_current <- phm_installed && "fitPHMEnsemble" %in% getNamespaceExports("PHM")
if (phm_installed) ok(sprintf("PHM           %s", ver("PHM"))) else miss(sprintf("PHM           github.com/%s", github_repo))
if (phm_installed && !phm_current) {
  info("[OLD]     this PHM predates the ensemble functions (fitPHMEnsemble etc.); rerun with --update")
}

need_cran <- c(missing_req, missing_opt)
need_phm <- !phm_installed || update_phm

if (check_only) {
  cat("\n")
  if (length(missing_req) == 0 && phm_installed) {
    cat("Everything required for PHM is installed.\n")
    if (length(missing_opt)) cat("Optional packages missing:", paste(missing_opt, collapse = ", "), "\n")
    quit(status = 0)
  }
  cat("Not ready. Run without --check-only to install what is missing.\n")
  quit(status = 1)
}

## ---- install ----------------------------------------------------------------------------------
if (length(need_cran)) {
  cat("\nInstalling from CRAN:", paste(need_cran, collapse = ", "), "\n")
  if (identical(.Platform$pkgType, "source")) {
    tools_found <- Sys.which(c("make", "gcc", "gfortran"))
    info("This R installs packages from source, which needs compilers (C and Fortran) and can take several minutes.")
    if (any(tools_found == "")) {
      info(paste0("Not found on PATH: ", paste(names(tools_found)[tools_found == ""], collapse = ", "),
                  ". On macOS install the Xcode command line tools (xcode-select --install) and gfortran;",
                  " on Debian/Ubuntu: sudo apt-get install build-essential gfortran."))
    }
  }
  install.packages(need_cran, lib = lib, Ncpus = max(1L, min(4L, parallel::detectCores(), na.rm = TRUE)))
}

if (need_phm) {
  cat("\nInstalling PHM from GitHub (", github_repo, ")\n", sep = "")
  if (!installed("remotes")) install.packages("remotes", lib = lib)
  if (!installed("remotes")) stop("Could not install 'remotes', which is needed to install PHM from GitHub.")
  still_missing <- Filter(Negate(installed), names(required))
  if (length(still_missing)) {
    cat("Cannot install PHM yet: these dependencies failed to install:", paste(still_missing, collapse = ", "), "\n")
  } else {
    tryCatch(
      remotes::install_github(github_repo, lib = lib, upgrade = "never", dependencies = FALSE, force = update_phm),
      error = function(e) cat("PHM installation failed:", conditionMessage(e), "\n")
    )
  }
}

## ---- final status ---------------------------------------------------------------------------
cat("\n== Result ==\n")
final_missing <- Filter(Negate(installed), c(names(required), "PHM"))
if (length(final_missing)) {
  cat("Still missing:", paste(final_missing, collapse = ", "), "\n")
  cat("Scroll up for the installer's error messages (compiler and system-library problems are the usual cause).\n")
  quit(status = 1)
}
cat("All required packages are installed. PHM version:", ver("PHM"), "\n")
failed_opt <- Filter(Negate(installed), missing_opt)
if (length(failed_opt)) cat("Optional packages that failed to install:", paste(failed_opt, collapse = ", "), "\n")

## ---- end-to-end smoke run ---------------------------------------------------------------------
if (verify) {
  cat("\nRunning a tiny end-to-end check (three clusters -> PHMv2 -> plots)...\n")
  res <- tryCatch({
    suppressPackageStartupMessages({ library(PHM); library(mclust) })
    set.seed(1)
    x <- rbind(matrix(rnorm(300), ncol = 2), matrix(rnorm(300, 3), ncol = 2), matrix(rnorm(300, -3), ncol = 2))
    labels <- rep(1:3, each = 150)
    params <- constructPmcParamsWeightedPartition(labels, x, G = 1)
    phm2 <- PHMv2(paramsList = params, mcSamples = 2e4)
    if ("fitPHMEnsemble" %in% getNamespaceExports("PHM")) {
      ens <- fitPHMEnsemble(x, M = 2, rho = 0.5, G = 1:4, numCores = 1, mcSamples = 2e4, extendVerbose = FALSE)
      stopifnot(length(ens$paramsList) >= 2)
    }
    stopifnot(inherits(plotPHMv2Dendrogram(phm2), "ggplot"), inherits(plotPHMv2Heatmap(phm2), "ggplot"))
    stopifnot(all(is.finite(cumsum(phm2$mergeDeltaPmc))))
    TRUE
  }, error = function(e) { cat("Smoke run failed:", conditionMessage(e), "\n"); FALSE })
  if (!res) quit(status = 1)
  cat("PHM works. Use library(PHM) at the top of your R script (plus library(mclust) if you call Mclust() yourself).\n")
}
quit(status = 0)
