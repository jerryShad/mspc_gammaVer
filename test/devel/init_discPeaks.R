.init.discPeaks <- function(peakset, ovHit, replicate.type=c("Biological","Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  len <- length(peakset)-1,
                  len <- length(peakset))
  cnt.ovHit <- as.matrix(Reduce('+', lapply(ovHit, lengths)))
  keepMe <- cnt.ovHit >= min.c
  dropList <- lapply(Hit, function(ele_) ele_[!keep_me])
  init.discardPeaks <- Map(unlist,
                         mapply(extractList, peakset, dropList))
  return(init.discardPeaks)
}
