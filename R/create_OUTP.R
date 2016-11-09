# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .create_OUTP
#' @param peaks
#' @param pAdjustMethod
#' @param alpha threhold
#' @return GRanges
#' @export
#' @importFrom stats, p.adjust
#' @author Julaiti Shayiding
#' @example

.create_OUTP <- function(peaks, pAdjustMethod="BH", alpha=0.05) {
  # input param checking
  stopifnot(class(peaks)=="GRanges")
  stopifnot(is.numeric(alpha))
  pAdjustMethod = match.arg(pAdjustMethod)
  if(is.null(peaks$p.value)) {
    stop("required slot is missing")
  } else {
    p <- peaks$p.value
    p.adj <- p.adjust(p, method = pAdjustMethod)
    peaks$p.adj <- p.adj
  }
  keepMe <- sapply(peaks, function(elm) {
    res <- elm$p.adj < alpha
  })
  ans <- list(
    keep=peaks[keepMe],
    droped=peaks[!keepMe])
  return(ans)
}

#' example
gr <- GRanges( seqnames=Rle("chr1", 4),ranges=IRanges(c(3,33,54,91), c(23,42,71,107)),
               rangeName=c("a1", "a4", "a7", "a11"), p.value=c(1e-22, 1e-6,1e-13, 1e-7))

.FDR(peaks = gr, pAdjustMethod = "BH", alpha = 0.05)
