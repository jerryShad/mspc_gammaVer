## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .denoise_peakFiles
##' @param peakset
##' @param denoise_threshold threshold for filtering out all noise features from wholse replicates
##' @return GRangesList list of GRanges without any background noise
##' @export
##' @importFrom rtracklayer score
##' @author Julaiti Shayiding
##' @example
## .denoise_peakFiles(peakset = myDat, denoise_threshold = 1E-4)

.denoise_peakFiles <- function(peakset, denoise_threshold=1E-4, verbose=FALSE, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(denoise_threshold))
  if (verbose) {
    cat(">> filter out background noise peaks from all replicates simultanously...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  all.peaks <- lapply(peakset, function(ele_) {
    if(is.null(ele_$p.value)) {
      peaks <- .pvalueConversion(ele_, pvalueBase = 1L)
    } else {
      keep <- subset(ele_, ele_$p.value < denoise_threshold)
      keep
    }
  })
  res <- all.peaks
  return(res)
}
