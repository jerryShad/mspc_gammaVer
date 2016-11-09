# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .setPurification
#' @param peakList list of all confirmed and discarded peaks for each replicates
#' @param replicate.type indicate type of current sample
#' @return GRanges
#' @export
#' @importFrom dplyr setdiff
#' @author Julaiti Shayiding

.setPurification <- function(grs, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(inherits(grs[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  DF <- lapply(grs, function(elm) {
    out <- as(elm, "data.frame")
  })
  if(replicate.type=="Biological") {
    res <- DF[[1]]
  } else {
    res <- setdiff(DF[[1]], DF[[2]])
  }
  out <- as(res, "GRanges")
  return(out)
}

#' @example
myList <- GRangesList(
  saved <- GRanges( seqnames=Rle("chr1", 4),ranges=IRanges(c(3,33,54,91), c(23,42,71,107)),
                    rangeName=c("a1", "a4", "a7", "a11"), p.value=c(1e-22, 1e-6,1e-13, 1e-7)),
  droped <- GRanges(seqnames=Rle("chr1", 5),ranges=IRanges(c(25,33,47,74,91), c(29,42,51,81,107)),
                    rangeName=c("a2", "a4", "a6", "a8", "a11"),
                    p.value=c(1e-3, 1e-6, 1e-4, 1e-5,1e-7))
)

# testme
clean.Set <- .setPurification(myList, "Biological")



