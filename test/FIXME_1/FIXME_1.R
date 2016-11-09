## FIXME

.get.pvalue <- function(ovHit, gr, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(gr)=="GRanges")
  if(!is(ovHit[[1]], "CompressedIntegerList")) {
    stop("invalid input, entry must be IntegerList")
  }
  out <- extractList(gr$p.value, ovHit)
  out <- unlist(out)
  return(out)
}

mapply(.get.pvalue, res.initFilt, all.peakFile)

.helper.getPval <- function(ovHit, grs, ...) {
  # input param checking
  pvlist <- mapply(.get.pvalue, grs)
  res <- sapply(pvlist, function(ele_) {
    out <- ifelse(length(x)>0, x, 0)
    out
  })
  return(res)
}


get.fisherScore <- function(list.pvalue, verbose=FALSE, ...) {
  # input param checking
  Fisher.score <- suppressWarnings(
    out <- apply(list.pvalue[,], 1, function(ro) {
      ans <- sumlog(ro)$p
      ans
    })
  )
  Fisher.score <- unlist(as.integer(Fisher.score))
  return(Fisher.score)
}

corr.VecDim <- lapply(allPeak.val, function(ele_) {
  tmp_pval.df <- data.frame(cbind("p"=ele_, "comb.pvalue"=Fisher.score))
  out <- tmp_pval.df[tmp_pval.df$p != 0.000000e+00, ]
  out$p <- NULL
  out
})

ov.expand <- ov.sample[unlist(extractList(seq_along(ov.sample), ov.hit))]
all.expand <- mapply(ov.expand, ov.hit, ov.sample)
res <- Map(cbind, all.expand, corr.VecDim)


##=================================================================================
###================================================================================
# MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
#
#
#' @title getPvalue
#' @param hit.list
#' @param peakSamples
#' @param
#' @return
#' @export
#' @importFrom
#' @importFrom
#' @author
#' @example

extract.pvalue <- function(hit.list, peakSamples, ...) {
  # input param checking
  stopifnot(inherits(peakSamples[[1]], "GRanges"))
  if (verbose) {
    cat(">> retrieving pvalue for all peaks that sufficiently meet filtering criterian...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  .keep <- .keep.peaks(hit.list, peakSamples, ...)
  # FIX ME : find good way
  res <- mapply(helper.getP, .keep, peakSamples)
  res <- data.frame(mapply(.pval.helper, res))
  return(res)
}

helper.getP <- function(hit, gr, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(gr)=="GRanges")
  stopifnot(missingArg(hit))
  .getP <- extractList(gr$p.value, hit)
  res <- unlist(.getP)
  return(res)
}
