## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @name mspc.core
##' @param
##' @param
##' @return
##' @export
##' @import dplyr
##' @author
##'
.loadPeak <- function(peakFolder, tau.w=1.0E-04, ...) {
  # check input param
  require(dplyr)
  if(missing(peakFolder)) {
    stop("input parameter is missing")
  }
  stopifnot(length(peakFolder)>0)
  if(is(peakFolder[[1]], "GRanges")) {
    peak.grs <- peakFolder
  } else {
    files <- list.files(peakFolder, full.names = TRUE)
    peak.grs <- setNames(
      lapply(files, function(ele_) {
        out <- as(import.bed(ele_), "GRanges")
      }), tools::file_path_sans_ext(basename(files))
    )
  }
  peak.grs <- sort(unique(peak.grs))
  stopifnot(is.numeric(tau.w))
  #noi.DR <- dir.create("data/noise")
  resl <- lapply(peak.grs, function(ele_) {
    if(is.null(ele_$p.value)) {
      ele_ <- .pvalueConversion(ele_, pvalueBase = 1L)
      .splitMe <- function(ele_, tau.w) {
        require(dplyr)
        ww <- as(ele_, "data.frame")
        ww %>%
          filter(p.value > tau.w) %>%
          write.table(., sprintf("noise.%s.bed",ele_), quote=F, sep="\t")
        total.ERs <- as(filter(ww, p.value <= tau.w),"GRanges")
        return(total.ERs)
      }
      res <- splitMe(ele_, tau.w)
    }
  })
  names(resl) <- names(peak.grs)
  return(resl)
}

##-----------------------------------------------------------------------

.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$pvalue)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  return(x)
}


.splitMe <- function(mList, ele_, tau.w) {
  require(dplyr)
  ww <- mList[[ele_]]
  ww %>%
    filter(p.value > tau.w) %>%
    write.table(., sprintf("dropped.%s.csv", ele_), row.names = FALSE)
  total.ERs <- filter(ww, p.value <= tau.w)
  return(total.ERs)
}

#====================================================================================

mylist2 <- lapply(seq_along(myData), function(i) {
  if(is.null(i$p.value)) {
    i <- .pvalueConversion(i, 1L)
  } else {
    splitter <- function(i, tau.w) {
      require(dplyr)
      DF <- myData[[i]]
      DF %>%
        filter(p.value >= tau.w) %>%
        write.csv(., sprintf("dropped.%s.csv", i), row.names = FALSE)
      total.ERs <- filter(DF, p.value <= tau.w)
      return(total.ERs)
    }
    res <- splitter(i, 1.0E-04)
  }
})

