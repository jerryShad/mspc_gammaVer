## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .keep.sufficientPeaks
##' @param ovHit list of overlap hit
##' @param peakset all peak files after background noise eleminated
##' @param repllicate.type Parameter that indicate whether peakfiles as Biological or Technical
##' @return IntegerList list of overlap hit that meet the filtering condition
##' @export
##' @author Julaliti Shayiding
##' @example

.keep.peaks <- function(ovHit, peakset, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  if(!is(peakset[[1]], "GRanges")) {
    stop("input must conform with GRanges object")
  }
  #stopifnot(length(peakset)>=2)
  if(!is(ovHit[[1]], "CompressedIntegerList")) {
    stop("invalid input")
  }
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  nn <- .cnt.ovnum(ovHit)
  .keep <- nn >= min.c
  out <- lapply(ovHit, function(ele_) {
    ans <- ele_[.keep]
    ans
  })
  res <- out
  return(res)
}
