## MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
##
## executble workflow of MSPC Package

##===========================================================================
## readPeakFiles

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

# example
myData <- readPeakFiles(peakFolder = "data/")

##============================================================================
## .pvalueConversion

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

##============================================================================
## .denoise_peakFiles

.denoise_peakFiles <- function(peakFolder, tau.w=1E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(tau.w))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  res <- lapply(peakFolder, function(ele_) {
    ans <- .pvalueConversion(ele_, pvalueBase=1L)
    ans
  })
  ## TODO: FIXME & Beautify me as an elegant version
  filt <- lapply(res, function(ele_) {
    out <- subset(ele_, ele_$p.value < tau.w)
    out
  })
  ## TO END;
  res <- filt
  return(res)
}

# example
all.peakFile <- .denoise_peakFiles(peakFolder = myData, tau.w = 1e-04)

##===========================================================================================
## peakOverlapping

.peakOverlapping <- function(peakset, idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(idx))
  #FUN <- match.arg(FUN, which.min)

  # set up the entry
  chosen <- peakset[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(peakset[- idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}

# example (a.k.a, all needed test must be done):
test_1 <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.min)
test_2 <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.min)
test_3 <- .peakOverlapping(peakset = all.peakFile, idx = 3L, FUN = which.min)

#hitlist <- list(test1 = test_1, test2 = test_2, test3 = test_3)

##=================================================================================
## manipulate all overlap hits from different test to one final hit list

Hit2IntVec <- function(ovHit, verbose=FALSE, ...) {
  # input param cheching
  if(verbose) {
    cat(">> loading your hitlist...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  res <- lapply(ovHit, function(ele_) {
    out <- drop(ele_)
    out[is.na(out)] <- 0L
    out
  })
  result <- do.call("cbind", res)
  return(result)
}

##==================================================================================
## FIXME : make below sub function more compatible
#-----------------------------------------------------------------------------------
.hitVT <- Map(Hit2IntVec, hitlist)
.hitDF <- lapply(.hitVT, as.data.frame)

library(dplyr)
.Hit <- bind_rows(.hitDF) %>% distinct %>%
.hitMT <- as.matrix(.Hit)

.hitOVLP <- apply(.hitMT, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

#------------------------------------------------------------------------------------

##===================================================================================
## .cnt.ovnum : count overlapped peaks

.cnt.ovnum <- function(ovHit, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  if(!is(ovHit[[1]], "CompressedIntegerList")) {
    stop("invalid input, entry hitlist must be IntegerList")
  }
  tot.num <- Reduce('+', lapply(ovHit, lengths))
  res <- tot.num
  return(res)
}

##====================================================================================
## .keep.peaks

.keep.peaks <- function(ovHit, peakset, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  if(!is(peakset[[1]], "GRanges")) {
    stop("input must conform with GRanges object")
  }
  #stopifnot(length(peakset)>=2)
  if(!is(ovHit[[1]], "CompressedIntegerList")) {
    stop("invalid input")
  }
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  nn <- .cnt.ovnum(ovHit)
  .keep <- nn >= min.c
  out <- lapply(ovHit, function(ele_) {
    ans <- ele_[.keep]
    ans
  })
  return(out)
}

# example
keep_list <- .keep.peaks(ovHit = .hitOVLP, peakset = all.peakFile, replicate.type = "Biological")

##====================================================================================
## .get.pvalue : retrieve pvalue for preserved peaks from sufficient overlap condition
.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

##--------------------------------------------------------------
# testme:
pvl.list <- mapply(.get.pvalue, keep_list, all.peakFile)
##====================================================================================
.helper.getPvalue <- function(pvlist, ...) {
  # input param checking
  res <- sapply(pvlist, function(x) {
    out <- ifelse(length(x)>0,
                  x,
                  0)
    out
  })
  return(res)
}

##-----------------------------------------------------------------------------------
# testme:
list.pval <- data.frame(mapply(.helper.getPvalue, pvl.list))
##=====================================================================================
## Fisher method : getting combined pvalue for overlapped peaks

get.fisherScore <- function(pv.list, verbose=FALSE, ...) {
  # input param checking
  require(metap)
  Fisher.score <- suppressWarnings(
    out <- apply(pv.list[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  res <- Fisher.score
  return(res)
}

# testme:
comb.pval <- as(get.fisherScore(list.pval), "List")

#--------------------------------------------------------------------------------------
confirm <- sapply(comb.pval, function(x) x <= 1E-08)
Confirmed_idx <- lapply(keep_list, function(x) x[confirm])
Discarded_idx.1 <- lapply(keep_list, function(x) x[!confirm])

#--------------------------------------------------------------------------------------
Discard.peaks_1 <- mapply(extractList, all.peakFile, Discarded_idx.1)
Discard.peaks_1 <- lapply(Discard.peaks_1, function(ele_) {
  ans <- unlist(ele_)
  ans
})

Discard.peaks_1 <- lapply(Discard.peaks_1, function(elm) {
  ans <- elm[!duplicated(elm),]
  ans
})

##=====================================================================================
confirmed.peaks <- mapply(extractList, all.peakFile, Confirmed_idx)

confirmed.peaks <- lapply(confirmed.peaks, function(ele_) {
  out <- unlist(ele_)
  out
})

confirmed.peaks <- lapply(confirmed.peaks, function(ele_) {
  out <- ele_[!duplicated(ele_),]
  out
})
#========================================================================================














#-----------------------------------------------------------------------------------#
#####################################################################################
#-----------------------------------------------------------------------------------#
myList <- list(f1=IntegerList(1,2,3,4,1,1,1,integer(0),1,2,4),
               f2=IntegerList(1,5,integer(0),integer(0),2,3,4,6,1,5,6),
               f3=IntegerList(1,4,6,7,2,3,3,7,2,5,7))

len <- Reduce('+', lapply(myList, lengths))

keepMe <- len >= length(myList)

res.filt <- lapply(myList, function(elm) {
  ans <- list(
    keep_=elm[keepMe],
    droped_=elm[!keepMe])
  ans
})

library(purrr)

map(myList, lengths) %>%
  reduce(`+`) %>%
  map_lgl(`>=`, length(myList)) -> keep_me

keep_list <- map(myList, ~.[keep_me])
drop_list <- map(myList, ~.[!keep_me])

discarded_peaks <- mapply(extractList, peak.list, drop_list)
allDiscarded_0 <- lapply(discarded_peaks, function(ele_) {
  ans <- unlist(ele_)
  ans
})

