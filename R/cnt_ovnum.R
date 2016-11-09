## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .cnt.ovnum
##' @param hitlist list of overlap hit
##' @return Integer vector
##' @export
##' @importFrom S4Vectors Reduce
##' @importFrom S4Vectors lengths
##' @author Julaiti Shayiding
##' @example

.cnt.ovnum <- function(ovHit, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  tot.num <- Reduce('+', lapply(ovHit, lengths))
  res <- tot.num
  return(res)
}
