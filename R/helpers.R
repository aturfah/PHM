## Misc Functions

#' Print at specified verbosity level
#'
#' @description Prints the provided text (using \code{cat()}) with a new line if \code{vValue} is not smaller than
#'   \code{vLevel}. \code{TRUE} always prints; \code{FALSE} never prints.
#'
#' @param ... Arguments passed to [cat()]. The indent is prepended directly to
#'   the first argument with no separator.
#' @param vValue Verbosity value provided by the user. \code{TRUE} always
#'   prints; \code{FALSE} never prints; a positive integer prints when
#'   \code{vValue >= vLevel}.
#' @param vLevel Verbosity level required for this print to fire.
#' @param indentLevel Number of times to replicate \code{indentChar} before
#'   the printed text. Defaults to \code{vLevel-1}.
#' @param indentChar Character(s) used to indent. Default is \code{".."}.
#' @param sep Separator passed to [cat()] between values in \code{...}.
#'   Default is \code{" "}.
#'
verbosityPrint <- function(..., vValue, vLevel, indentLevel=vLevel-1, indentChar="..", sep=" ") {
    should_print <- isTRUE(vValue) || (is.numeric(vValue) && vValue >= vLevel)
    if (should_print) {
        leading_text <- paste0(rep(indentChar, indentLevel), collapse="")
        args <- list(...)
        args[[1]] <- paste0(leading_text, args[[1]])
        do.call(cat, c(args, list(sep=sep)))
        cat("\n")
    }
}