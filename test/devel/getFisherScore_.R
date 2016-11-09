## MSPC Project -  Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .get.fisherScore
##' @param pvalueList
##' @param verbose
##' @return IntegerList
##' @export
##' @importFrom metap sumlog
##' @author Julaiti Shayiding
##' @example

get.fisherScore <- function(allPeask.pval, ...) {
  Fisher.score <- suppressWarnings(
    out <- apply(allPeask.pval[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  return(Fisher.score)
}
