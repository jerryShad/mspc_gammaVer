## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .peakOverlapping
##' @param peakset
##' @param cur.idx current index of peakFile used as chosen sample
##' @param FUN function selected by user needs when keeping most stringent peak when multiple overlap occur
##' @return IntegerList overlap-hit IntegerList across multiple GRanges simulanously
##' @export
##' @importFrom GenomicRanges findOverlaps
##' @importFrom XVector extractList
##' @importfrom rtracklayer score
##' @importFrom IRanges which.min
##' @author Julaiti Shayiding
##' @example

.peakOverlapping <- function(peakset, cur.idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(cur.idx))
  FUN = match.arg(FUN)   #FIXME
  # set up the entry
  chosen <- peakset[[cur.idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(set[- cur.idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[cur.idx]),names(peakset[- cur.idx]))
  return(res)
}
