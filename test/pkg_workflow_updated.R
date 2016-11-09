## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##' @description
##' all updated function code and workflow of MSPC Packages
##'
##==============================================================================================
##' @details read bed files as GRanges objects from file directory
##' function implementation

readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  require
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  stopifnot(length(peakFolder)>=1)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      res <- as(import.bed(ele_), "GRanges")

    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}

##' @examples
myData <- readPeakFiles(peakFolder = "data/")

##==============================================================================================
##' @details splitting GRanges
##' function implementation

.denoise_peakFiles <- function(peakset, tau.w=1.0E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  stopifnot(is.numeric(tau.w))
  #dir.create(file.path(outDir), showWarnings = FALSE)
  #setwd(file.path(outDir))
  res <- lapply(peakset, function(ele_) {
    if(is.null(ele_$p.value)) {
      ele_ <- .pvalueConversion(ele_, pvalueBase = 1L)
    }
    DF <- as.data.frame(ele_)
    require(dplyr)
    DF %>%
      filter(p.value > tau.w) %>%
      export.bed(., sprintf("BackgroundNoise%s.bed", DF), row.names = FALSE)
    total.ERs <- filter(DF, p.value <= tau.w)
    total.ERs <- as(total.ERs, "GRanges")
    return(total.ERs)
  })
}

##========================================================

.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$p.value)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  return(x)
}

##' @example
total.ERs <- .denoise_peakFiles(peakset = myData, tau.w = 1.0E-04)

##==============================================================================================
##' @details peakOverlapping
##' function implementation

.peakOverlapping <- function(peakset, idx=1L, FUN=which.max, ...) {
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
    out.idx0 <- out.idx0[!is.na(out.idx0),]
    ans <- ans[out.idx0]
    ans
  })
  res <- DataFrame(c(list(que.hit), sup.hit))
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}

##' @example
.hit_1 <- .peakOverlapping(peakset = total.ERs, idx = 1L, FUN = which.max)
.hit_2 <- .peakOverlapping(peakset = total.ERs, idx = 2L, FUN = which.max)
.hit_3 <- .peakOverlapping(peakset = total.ERs, idx = 3L, FUN = which.max)

#-----------------------------------------------------------------------------------------------
#' @description integrate hit table to remove duplicated hits

t1 <- as.matrix(DataFrame(.hit_1[names(.hit_3)]))
t2 <- as.matrix(DataFrame(.hit_2[names(.hit_3)]))
t3 <- as.matrix(DataFrame(.hit_3))

#testMe <- unique(mapply(rbind, list(t1,t2,t3)))

#-------------------------------------------------------
Hit <- unique(DataFrame(rbind(t3,t2,t1)))
Hit <- lapply(Hit, function(ele_) {
  res <- as(ele_, "CompressedIntegerList")
  res
})

##=============================================================================================

MSPC.Analyzer <- function(peakset, ovHit, replicate.type=c("Biological","Technical"), tau.s=1.0E-08, ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(tau.s))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  length(peakset)-1,
                  length(peakset))
  cnt.ovHit <- as.matrix(Reduce('+', lapply(ovHit, lengths)))
  keepMe <- cnt.ovHit >= min.c

  dropList <- lapply(ovHit, function(ele_) ele_[!keep_me])
  init.discardPeaks <- Map(unlist,
                           mapply(extractList, peakset, dropList))

  keepList <- lapply(ovHit, function(x) x[keepMe])

  pval_List <- mapply(.get.pvalue, keepList, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
  }
  pval.TB <- as.data.frame(Map(.helper.PV, pval_List))
  #pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  #Npresent <- rowSums( !is.na(pval.TB) )
  #comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  comb.pval <- suppressWarnings(
    .globSC <- apply(pval.TB, 1, function(ro) {
      require(metap)
      res <- sumlog(ro)$p
    })
  )
  comb.pval <- as.matrix(comb.pval)
  ## FIXME: to make more compatible
  Confirmed_idx <- lapply(keepList, function(elm) {
    saved <- sapply(comb.pval, function(x) x <= tau.s)
    res <- elm[saved]
  })
  confirmed <- Map(extractList, peakset, Confirmed_idx)
  .Confirmed.ERs <- lapply(confirmed, function(ele_) {
    res <- unlist(ele_)
    out <- res[!duplicated(res),]
    out
  })
  ## TODO
  Discarded_idx <- lapply(keepList, function(elm) {
    droped <- sapply(comb.pval, function(x) x > tau.s)
    res <- elm[droped]
  })

  .Fisher.discPeaks <- Map(extractList, peakset, Discarded_idx)
  .Fisher.discPeaks <- lapply(.Fisher.discPeaks, function(elm) {
    res <- unlist(elm)
    ans <- res[!duplicated(res),]
    ans
  })
  .Discarded.ERs <- suppressWarnings(mapply(c, .init.discPeaks, .Fisher.discPeaks))
  ## TODO END !
  ##---------------------------------------
  .setPurification <- ifelse(replicate.type=="Biological",
                             res <- .Confirmed.ERs,
                             res <- Map(anti_join, .Confirmed.ERs, .Discarded.ERs))

  ##---------------------------------------
  rslt <- c(.Confirmed.ERs, .Discarded.ERs)
  ans.Reslt <- split(rslt, sub("_.*", "", names(rslt)))[sub("_.*", "", names(.Confirmed.ERs))]
  return(ans.Reslt)
}

.get.pvalue <- function(hit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, hit)
  return(res)
}

##================================================================================================
##' @description
##' This script is used for how to get stringent peakset, weakpeaks set like Musera

library(data.table)
confirmed_dt <- rbindlist(all_Confirmed,idcol=TRUE)
discarded_dt <- rbindlist(all_Discarded,idcol=TRUE)

DT <- rbind(confirmed_dt, discarded_dt)
DT[,gr := ifelse(p.value <= tau.S,"Stringent","Weak"),.id]
DT[,.id:= gsub("[.](confirmed|discarded)","",.id)]
res <- by(DT,DT$.id,FUN = function(x) split(x,x$gr))


##================================================================================================
##' @description
##' This scripts how to manipulate list of confirmed peaks, discarded peaks for BH correction test
x <- c(all_Confirmed, all_Discarded)
ans.Reslt <- split(x, sub("_.*", "", names(x)))[sub("_.*", "", names(all_Confirmed))]

##================================================================================================

func <- function(L1, L2, replicate.type=c("Biological","Technical"))
{
  # input param checking
  if(!is(L1[[1]], "GRanges") || !is(L1[1], "GRanges")) {
    stop("invalid input, entry must be GRanges")
  }
  replicate.type = match.arg(replicate.type)
  if(!type %in% c("Biological","Technical")){
    stop("wrong type, please ")
  }
  if(replicate.type=="Biological") {
    res <- L1
  } else {
    res <- Map(function(x,y) anti_join(x,y), L1, L2)
  }
  return(res)
}

L1 <- Map("data.frame", all_Confirmed)
L2 <- Map("data.frame", all_Discarded)

# testme:
for_BH <- func(L1, L2, "Biological")

#---------------------------------------------------------------------------------------------------
.create_OUTP <- function(peaks, pAdjustMethod="BH", alpha=0.05, ...) {
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

##==================================================================================================


