## MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
##
##' @title whole Pipeline of MSPC Packages
##' @description
##' This R scripts served as blue print of MSPC Packages, it render user
##' all pipeline of MSPC Packages to analyze custom data set

#========================================================================================================
#' @description: read bed files as GRanges objects
#'
#' @details function implementation

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

#========================================================================================================
#' @description : attach p.value as new metadata all peaks
#' @details function implementation

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

#========================================================================================================
#' TODO : remove all background noise from all peak files simultanously
#' function implementation

.denoise_peakFiles <- function(peakFolder, tau.w=1.0E-04, verbose=FALSE, ...) {
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
  peaks_rmNoi <- lapply(res, function(ele_) {
    out <- subset(ele_, ele_$p.value <= tau.w)
    out
  })
  ## TO END;
  res <- peaks_rmNoi
  return(res)
}

# testme:
all.peakFile <- .denoise_peakFiles(peakFolder = myData, tau.w = 1.0E-04)

#========================================================================================================
#' TODO : peak overlapping all peak files simultanously
#' function implementation

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
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}


# testme:
test_1 <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.max)
test_2 <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.max)
test_3 <- .peakOverlapping(peakset = all.peakFile , idx = 3L, FUN = which.max)

hitTB <- list(test.1 = test_1, test.2 = test_2, test.3 = test_3)

#=======================================================================================================
#' TODO: to achieve final overlap hit table where all possible overlap pair are included
#' function implementation

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

# testme:
.hitVT <- Map(Hit2IntVec, hitTB)
.hitDF <- lapply(.hitVT, as.data.frame)

library(dplyr)
.Hit <- bind_rows(.hitDF) %>% distinct
.hitMT <- as.matrix(.Hit)

.hitTB <- apply(.hitMT, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

#=======================================================================================================
#' TODO :
#' function implementation
.filtByIndx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  ovNum <- .cnt.ovnum(ovHit)
  keepMe <- min.c >= ovNum
  return(keepMe)
}

# testme :
 keep_me <- .filtByIndx(peakset = all.peakFile, .vHit, replicate.type = "Biological")
#
# library(purrr)
#
# #keep_list <- map(.final_Hit, ~.[keep_me])
drop_list <- map(.vHit, ~.[!keep_me])

#------------------------------------------------------------------------------------------------------
.cnt.ovnum <- function(ovHit, verbose=FALSE, ...) {
  tot.num <- Reduce('+', lapply(ovHit, lengths))
  res <- tot.num
  return(res)
}
#------------------------------------------------------------------------------------------------------

.keep_hit <- function(ovHit, peakset, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  replicate.type = match.arg(replicate.type)
  if(replicate.type=="Biological")
    min.c <- length(peakset)-1
  else
    min.c <- length(peakset)
  #min.c <- as.integer(min.c)
  nn <- .cnt.ovnum(ovHit)
  .keep <- sapply(nn, function(elm) {
    out <- elm >= min.c
  })
  res <- lapply(ovHit, function(ele_) {
    ans <- ele_[.keep]
    ans
  })
  return(res)
}

#------------------------------------------------------------------------------------------------------
.drop_hit <- function(ovHit, peakset, replicate.type=c("Biological", "Technical"), ...) {
  replicate.type = match.arg(replicate.type)
  if(replicate.type=="Biological")
    min.c <- length(peakset)-1
  else
    min.c <- length(peakset)
  nn <- .cnt.ovnum(ovHit)
  .drop <- nn < min.c
  res <- lapply(ovHit, function(ele_) {
    res <- ele_[.drop]
    res
  })
  return(.drop)
}
# testme:
keep_list <- .keep_hit(.vHit, all.peakFile, "Biological")
drop_list <- .drop_hit(.vHit, all.peakFile, "Biological")

#=======================================================================================================
#' @description:
#' discarded peaks because of insufficient overlapped are must be expanded as GRanges, becaues I need to
#' do statistical summary for all discarded peaks when running MSPC Packages
#'
#' @details This is function implementation

droped_list <- lapply(drop_list, function(ele_) {
  ele_ <- drop(ele_)
  ele_ <- ele_[!is.na(ele_)]
  #res <- as(ele_, "IntegerList")
})

Discard.peaks_0 <- mapply(extractList, all.peakFile, droped_list)
Discard.peaks_0 <- mapply(unlist, Discard.peaks_0)

#========================================================================================================
#' @description :
#' peaks with sufficient overlapped hit are not immediatly expanded,
#' instead I'll take its pvalue table based on its overlap hit for finding
#' its global score by Fisher method
#'
#' @details This is the implementation to retrieve corresponding p.value of peaks

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

# testme:
pvl.list <- mapply(.get.pvalue, keep_list, all.peakFile)

#-------------------------------------------------------------------------------------------------------
#' @description :
#'  pvalue table contains zero value, becaues in my overlap hit table, some overlap hit pairs possibly
#' contains zero hit index but it still valid for involving in Fisher method;
#' Here I replace zero index with special index that can be accepted for proceeding Fisher Method
#'
#' @details This is small utility functon that fix this issue efficiently

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

# testme:
list.pval <- data.frame(mapply(.helper.getPvalue, pvl.list))

#=======================================================================================================
#' @description
#' after retrieving pvalue of peaks who has sufficient overlap hit, it is ready to be involved in Fisher method
#' to get global fisher score
#'
#' @details This is the implementation how this done

get.fisherScore <- function(p.list, ...) {
  # input param checking
  require(metap)
  Fisher.score <- suppressWarnings(
    out <- apply(p.list[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  return(Fisher.score)
}

# testme:
comb.pval <- get.fisherScore(list.pval)

#=======================================================================================================
#' @description
#' based on acquired global fisher score,
#' I need to do further filtering for comparing with each of Fisher score with given threshold value, and
#' generate corresponding Confirmed, Discarded peakset respectively
#'
#' @details This is the implementation (still trivial, I can make it)

saveMe <- sapply(comb.pval, function(x) x <= 1.0E-08)
Confirmed_idx <- lapply(keep_list, function(x) x[saveMe])
Discarded_idx <- lapply(keep_list, function(x) x[!saveMe])

#-------------------------------------------------------------------------------------------------------
#' @description
#' according to all confirmed_idx, discarde_idx, I intend to expand them as GRanges objects
#'
#' @details This is the implementation

confirmed.peaks <- mapply(extractList, all.peakFile, Confirmed_idx)
confirmed.peaks <- mapply(unlist, confirmed.peaks)
confirmed.peaks <- lapply(confirmed.peaks, function(ele_) {
  out <- ele_[!duplicated(ele_),]
  out
})


#-------------------------------------------------------------------------------------------------------
Discard.peaks_1 <- mapply(extractList, all.peakFile, Discarded_idx)
Discard.peaks_1 <- mapply(unlist, Discard.peaks_1)
Discard.peaks_1 <- lapply(Discard.peaks_1, function(elm) {
  ans <- elm[!duplicated(elm),]
  ans
})

#=======================================================================================================
#' @description
#' so far, I did two different filtering with different threshold;
#' I need to combine discarded peaks that failed of sufficient overlap condition and one with failed of Fishermethod
#'
#' @details
Discarde_peakList <- suppressWarnings(mapply(c, Discard.peaks_0, Discard.peaks_1))

#=======================================================================================================
#' @description
#' I need to set purification before proceeding BH correction test
#'
#' @details This is the implementation

.setPurification <- function(grs, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(inherits(grs[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
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

#=======================================================================================================
#' @description
#'
#' @details

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

#=======================================================================================================

