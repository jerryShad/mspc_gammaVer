##' MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##'
##' MSPC Package pipelie - Updated version
##' @description
##' In this script, I intend to double check every piece of code to make sure it is efficient;
##' Also optimize function structure to make complete overall workflow of MSPC Packages ;

##===========================================================================================================
##' @title readPeakFile
##' @description Read bed files from file directory as GRanges objects
##' @return GRanges objects

readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  stopifnot(length(peakFolder)>=1)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      out <- as(import.bed(ele_), "GRanges")
    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}

##' @examples
myData <- readPeakFiles(peakFolder = "peakDat/")

##===========================================================================================================
##' @title .pvalueConversion
##' @description add pvalue as new metadata column for each GRanges
##' @return GRanges objects

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

##===========================================================================================================
##' @title
##' @description
##' @return
##' @example

.denoise_peakFiles <- function(peakFolder, tau.w=1.0E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(tau.w))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  peaks_rmNoi <- lapply(peakFolder, function(ele_) {
    ele_ <- .pvalueConversion(ele_, pvalueBase = 1L)
    out <- subset(ele_, ele_$p.value <= tau.w)
    out
  })
  return(peaks_rmNoi)
}

##' @example
total.ERs <- .denoise_peakFiles(peakFolder = myData, tau.w = 1.0E-04)

##===========================================================================================================
##' @title .peakOverlapping
##' @description finding overlapped regions of current sample which supportedd with rest of replicate by pair-wise
##' @return Integer hit index
##' @example

.peakOverlapping <- function(peakset, idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(idx))
  # set up the entry
  chosen <- peakset[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(peakset[- idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- DataFrame(c(list(que.hit),sup.hit))
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}

##' @example
.hit_1 <- .peakOverlapping(peakset = total.ERs, idx = 1L, FUN = which.max)
.hit_2 <- .peakOverlapping(peakset = total.ERs, idx = 2L, FUN = which.max)
.hit_3 <- .peakOverlapping(peakset = total.ERs, idx = 3L, FUN = which.max)

##----------------------------------------------------------------------------------------------------------
#' @description make each hit table has same pattern

idx <- names(.hit_1)
.hit_2 <- .hit_2[idx]
.hit_3 <- .hit_3[idx]

##----------------------------------------------------------------------------------------------------------
#' @description let hit table as matrix representation by using DataFrame
idx <- names(.hit_1)
.hit_1 <- DataFrame(.hit_1)
.hit_2 <- DataFrame(.hit_2[idx])
.hit_3 <- DataFrame(.hit_3[idx])

Hit <- unique(mapply(rbind, DataFrame(.hit_1), DataFrame(.hit_2), DataFrame(.hit_3)))
Hit <- DataFrame(Hit)

Hit <- DataFrame(lapply(Hit, function(ele_) {
  res <- as(ele_, "CompressedIntegerList")
  res
}))

##==========================================================================================================
##' @title filtByHitIdx
##' @description
##' @return

.filtByHitIdx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  len <- length(peakset)-1,
                  len <- length(peakset))
  min.c <- as.integer(min.c)
  ovNum <- Reduce('+', lapply(ovHit, lengths))
  keepMe <- sapply(ovNum, function(elm) {
    res <- elm >= min.c
    res
  })
  return(keepMe)
}

#' @example
keep_me <- .filtByHitIdx(peakset = total.ERs, Hit,
                       replicate.type = "Biological")

lapply(Hit, `[`, keep_me, drop=FALSE)

keepList <- lapply(Hit, function(ele_) {
  res <- ele_[keep_me]
  res
})

dropList <- lapply(Hit, function(ele_) {
  res <- ele_[!keep_me]
  res
})

#' @example
#library(purrr)
#keepList <- map(Hit, ~.[keep_me])
#dropList <- map(Hit, ~.[!keep_me])

#-------------------------------------------------------------------------------------------------------
#' @description  Alternative solution for filtering overlap hit index given condition that we proposed

func <- function(peakset, ovHit,
                 replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  require(purrr)
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  map(ovHit, lengths) %>%
    reduce(`+`) %>%
    map_lgl(`>=`, min.c) -> keep_me
  return(keep_me)
}

keep_ <- func(total.ERs, .hit_1, "Technical")
keepL <- map(.hit_1, ~.[keep_])
dropL <- map(.hit_1, ~.[!keep_])

#-------------------------------------------------------------------------------------------------------
##' @description peaks that did not meet the sufficient overlap condition won't be proceed to next,
##' instead we keep them as GRanges to get overall discarded peaks for statistical result at the end

Discard.peaks_0 <- Map(unlist, mapply(extractList, total.ERs, dropList))
##======================================================================================================
##' @title Fisher.stats
##' @description get combined pvalue for overlapped peaks by pair wise
##' @return GRanges
##' @example

.Fisher.stats <- function(hitTB, peakset, verbose=FALSE, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
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
  #comb.pval <- apply(pval.TB[1:3], 1,  function(df) sumlog(df)$p)
  pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  Npresent <- rowSums( !is.na(pval.TB) )
  comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  #res <- DataFrame(comb.pval)
  return(comb.pval)
}

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

.fish_stat <- .Fisher.stats(.hit_1, total.ERs)

##====================================================================================================
##' @title
##' @description
##' @return
##' @example



##====================================================================================================
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
  if(!replicate.type %in% c("Biological", "Technical")) {
    stop("wrong type")
  }
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

##===========================================================================================
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
