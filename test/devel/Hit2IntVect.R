## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title
##' @param
##' @param
##' @return
##' @export
##' @importFrom IRanges drop
##' @importFrom rtracklayer as.data.frame
##' @importFrom dplyr bind_rows
##' @importFrom IRanges as.matrix
##' @author  Julaiti Shayiding
##' @example

Hit2IntVec <- function(ovHit, verbose=FALSE, ...) {
  # input param cheching
  if(verbose) {
    cat(">> loading your hitlist...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  res <- lapply(ovHit, function(ele_) {
    out <- drop(ele_)
    out[is.na(out)] <- 0L
    out
  })
  result <- do.call("cbind", res)
  return(result)
}
