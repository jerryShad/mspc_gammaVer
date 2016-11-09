## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title .ovHitAsGRanges
##' @param hitlist
##' @param peakset
##' @param pvalList
##' @return GRangesList
##' @export
##' @importFrom XVector extractList
##' @importFrom IRanges cbind
##' @importFrom BioGenerics Map
##' @author Julaiti Shayiding
##' @example

.Hit2GRanges.expand <- function(ov.hit, ov.sample, idx, allPeak.val) {
  ov.expand <- ov.sample[unlist(extractList(seq_along(ov.sample), ov.hit))]
  all.expand <- mapply(ov.expand, ov.hit, ov.sample[c(idx, seq_len(length(ov.sample))[-idx])])
  corr.VecDim <- lapply(allPeak.val, function(ele_) {
    tmp_pval.df <- data.frame(cbind("p"=ele_, "comb.pvalue"=Fisher.score))
    out <- tmp_pval.df[tmp_pval.df$p != 0.000000e+00, ]
    out$p <- NULL
    out
  })
  res <- Map(cbind, all.expand, corr.VecDim)
  return(res)
}
