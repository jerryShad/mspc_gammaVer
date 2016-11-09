# MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title filtByIndex
#' @param peakset
#' @param ovHit
#' @param replicate.type
#' @return logical
#' @export
#' @author Julaiti Shayiding
#' @example


.filtByIndx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  ovNum <- .cnt.ovnum(ovHit)
  keepMe <- min.c >= ovNum
  return(keepMe)
}


