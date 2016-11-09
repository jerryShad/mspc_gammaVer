## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
##' @description
##' Most updated work flow of MSPC Packages; Must do double check

#-----------------------------------------------------------------------------------------
#' @description read peak files from given file directory as GRanges objects

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
myData <- readPeakFiles(peakFolder = "data/")

#-----------------------------------------------------------------------------------------
#' @description pvalueConversion
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

#-----------------------------------------------------------------------------------------
#' @description remove all background noise peaks from GRanges objects

.denoise_peakFiles <- function(peakset, tau.w=1.0E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(tau.w))
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  reslt <- lapply(peakset, function(ele_) {
    if(is.null(ele_$p.value)) {
      ele_ <- .pvalueConversion(ele_, pvalueBase = 1L)
      #res <- split(ele_, ifelse(ele_$p.value <= tau.w, "total.ERs", "noise"))
      res <- subset(ele_, ele_$p.value <= tau.w)
      res
    }
  })
  return(reslt)
}

##' @example
total.ERs <- .denoise_peakFiles(peakset = myData, tau.w = 1.0E-04)

#--------------------------------------------------------------------------------------
#' @description peakOverlapping over all peak files simulatanously

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
    sort(ans)
  })
  res <- DataFrame(c(list(que=que.hit), sup=sup.hit))
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}

##' @example
.hit_1 <- .peakOverlapping(peakset = total.ERs, idx = 1L, FUN = which.max)
.hit_2 <- .peakOverlapping(peakset = total.ERs, idx = 2L, FUN = which.max)
.hit_3 <- .peakOverlapping(peakset = total.ERs, idx = 3L, FUN = which.max)

#---------------------------------------------------------------------------------------
#' @description integrate hit table to remove duplicated hits

t1 <- as.matrix(DataFrame(.hit_1))
t2 <- as.matrix(DataFrame(.hit_2[names(.hit_1)]))
t3 <- as.matrix(DataFrame(.hit_3[names(.hit_1)]))

#testMe <- unique(mapply(rbind, list(t1,t2,t3)))

#-------------------------------------------------------
Hit <- order(unique(DataFrame(rbind(t1, t2, t3))))
Hit <- lapply(Hit, function(ele_) as(ele_, "CompressedIntegerList"))

#Hit <- DataFrame(mapply(rbind, DataFrame(.hit_1),
#                               DataFrame(.hit_2), DataFrame(.hit_3)))
#Hit <- as.matrix(DataFrame(lapply(Hit, function(ele_) {
#  res <- as(ele_, "CompressedIntegerList")
#  res
#})))

#---------------------------------------------------------------------------------------
#' @description filter out overlap hit by condition that we proposed

.filtByHitIdx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  len <- length(peakset)-1,
                  len <- length(peakset))
  #min.c <- as.integer(min.c)
  ovNum <- as.matrix(Reduce('+', lapply(ovHit, lengths)))
  keepMe <- sapply(ovNum, function(elm) {
    res <- elm >= min.c
    res
  })
  return(keepMe)
}

#' @example
keep_me <- .filtByHitIdx(peakset = total.ERs, Hit,
                         replicate.type = "Biological")

keepList <- lapply(Hit, function(ele_) ele_[keep_me])
dropList <- lapply(Hit, function(ele_) ele_[!keep_me])

#-------------------------------------------------------------------------
Discard.peaks_0 <- Map(unlist, mapply(extractList, total.ERs, dropList))

#--------------------------------------------------------------------------------------
#' @description perform Fisher method to find combined pvalue

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  if(!is(ovHit[1L], "CompressedIntegerList")) {
    stop("entry hit list was not valid")
  }
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

#testme
pvl.list <- Map(.get.pvalue, keepList, total.ERs)

#---------------------------------------------------------------
helper.getPvalue <- function(pvlist, ...) {
  # input param checking
  res <- sapply(pvlist, function(x) {
    out <- ifelse(length(x)>0,
                  x,
                  0)
  })
  return(res)
}

list.pval <- as.matrix(DataFrame(Map(helper.getPvalue, pvl.list)))

get.fisherScore <- function(pv.list,...) {
  # input param checking
  require(metap)
  Fisher.score <- suppressWarnings(
    out <- apply(pv.list[,], 1, function(ro) {
      ans <- sumlog(ro)$p
      ans
    })
  )
  res <- as.matrix(Fisher.score)
  return(res)
}

# testme:
comb.pval <- get.fisherScore(list.pval)

#-------------------------------------------------------------------------------------------

Confirmed_idx <- lapply(keepList, function(ele_) {
  saveMe <- sapply(comb.pval, function(x) x <= 1.0E-08)
  res <- ele_[saveMe]
})
Discarded_idx <- lapply(keepList, function(x) {
  discard <- sapply(comb.pval, function(tt) tt > 1.0E-08)
  x[discard]
})
#rm(comb.pval)

confirmed <- Map(extractList, total.ERs, Confirmed_idx)

Discard.peaks_1 <- mapply(extractList, total.ERs, Discarded_idx)
Discard.peaks_1 <- mapply(unlist, Discard.peaks_1)

#------------------------------------------------------------------
all_Confirmed <- lapply(confirmed, function(ele_) {
  res <- unlist(ele_)
  out <- res[!duplicated(res)]
})

Discard.peaks_1 <- lapply(Discard.peaks_1, function(elm) {
  ans <- elm[!duplicated(elm),]
  ans
})

#------------------------------------------------------------------------------------
all_Discarded <- suppressWarnings(mapply(c, Discard.peaks_0, Discard.peaks_1))


##=========================================================
## Check stringent/weak confirmed or discarded peaks

library(data.table)
confirmed_dt <- rbindlist(all_Confirmed,idcol=TRUE)
discarded_dt <- rbindlist(all_Discarded,idcol=TRUE)

DT <- rbind(confirmed_dt, discarded_dt)
DT[,gr := ifelse(p.value <= tau.S,"Stringent","Weak"),.id]
DT[,.id:= gsub("[.](confirmed|discarded)","",.id)]
res <- by(DT,DT$.id,FUN = function(x) split(x,x$gr))

##===============================================================================

confirm_stats <- lapply(all_Confirmed, function(x) {
  res <- split(x, ifelse(x$p.value <= 1.0E-08,
                         "stringent",
                         "weak"))
})

discard_stats <- lapply(all_Discarded, function(x) {
  res <- split(x, ifelse(x$p.value <= 1.0E-08,
                         "stringent",
                         "weak"))
  res
})



bed_1_con <- all_Confirmed[[1]]
bed_1_dis <- all_Discarded[[1]]

res.stringent <- split(bed_1_con, ifelse(bed_1_con$p.value <= 1e-08, "stringent", "weak"))
res.weak <- split(bed_1_dis, ifelse(bed_1_dis$p.value <= 1e-08, "stringent", "weak"))

#-------------------------------------------------------------------------------------
x <- c(all_Confirmed, all_Discarded)
ans.Reslt <- split(x, sub("_.*", "", names(x)))[sub("_.*", "", names(all_Confirmed))]

#-------------------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------------
