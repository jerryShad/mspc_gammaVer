##' MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##'
##' @title .Fisher.stats
##' @param
##' @return
##' @export
##' @importFrom
##' @author
##' @example

.Fisher.stats <- function(hitTB, peakset, stringentThreshold = 1.0E-08, verbose=FALSE, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(stringentThreshold))
  pval_List <- mapply(.get.pvalue, hitTB, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
  }
  pval.TB <- Map(.helper.PV, pval_List)
  pval.TB <- data.frame(pval.TB)
  #pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  #Npresent <- rowSums( !is.na(pval.TB) )
  #comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  comb.pval <- suppressWarnings(
    .globSC <- apply(pval.TB, 1, function(ro) {
      res <- sumlog(ro)$p
    })
  )
  comb.pval <- DataFrame(comb.pval)
  confirmed_idx <- sapply(comb.pval, function(ele_) {
    res <- ele_ <= stringentThreshold
    res
  })
  return(confirmed_idx)
}

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}
