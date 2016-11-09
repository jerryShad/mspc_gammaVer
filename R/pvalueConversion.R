## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .pvalueCoversion
##' @param x GRanges objects
##' @param pvalueBase parameter for user defined pvalue range
##' @return GRanges
##' @export
##' @importFrom rtracklayer score
##' @importFrom rtracklayer mcols
##' @author Julaiti Shayiding
##' @example
## .pvalAttr <- .pvalueConversion("", pvalueBase = 1L)

.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  # input param checking
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$pvalue)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  res <- x
  return(res)
}
