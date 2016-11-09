# MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
#
# all workflow of MSPC Package
#'
#'
#========================================================================================
#' TODO : read bed files as GRanges objects
#' function implementation:

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

#======================================================================================
#' TODO : pvalue conversion for score of all peaks
#' function implementation:

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
#========================================================================================
#' TODO : remove all background noise from GRanges objects
#' Function implementation

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
    out <- subset(ele_, ele_$p.value <= tau.w)
    out
  })
  ## TO END;
  res <- filt
  return(res)
}

# testme:
all.peakFile <- .denoise_peakFiles(peakFolder = myData, tau.w = 1E-04)

#=======================================================================================
#' TODO : peakOverlapping
#' Function implementation

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

# testme:
test_1 <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.min)
test_2 <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.min)
test_3 <- .peakOverlapping(peakset = all.peakFile, idx = 3L, FUN = which.min)

hitlist <- list(test.1 = test_1, test.2 = test_2, test.3 = test_3)

#=====================================================================================
#' TODO: manipulate hitlist for solution of data.frame and to obtain final overlap hit table
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
.hitVT <- Map(Hit2IntVec, hitlist)
.hitDF <- lapply(.hitVT, as.data.frame)

library(dplyr)
.Hit <- bind_rows(.hitDF) %>% distinct
.hitMT <- as.matrix(.Hit)

.final_Hit <- apply(.hitMT, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

#====================================================================================
#' TODO : count all overlapped peaks and get vector sum by parallel
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

#====================================================================================
#' TODO : whether which overlapHit is deserved to go next phase by filtering
#' funcion implementation

min.c <- 2L

map(.final_Hit, lengths) %>%
  reduce(`+`) %>%
  map_lgl(`>=`, min.c) -> keep_me

keep_list <- map(.final_Hit, ~.[keep_me])
drop_list <- map(.final_Hit, ~.[!keep_me])

#======================================================================================
#' in order to expand all discarded peaks,
#' I need to reduce the dimension of hit table of dicarded peaks respectively
droped_list <- lapply(drop_list, function(ele_) {
  ele_ <- unique(drop(ele_))
  ele_ <- ele_[!is.na(ele_)]
  res <- as(ele_, "IntegerList")
})

discard.0_peaks <- mapply(extractList, all.peakFile, droped_list)
discard.0_peaks <- mapply(unlist, discard.0_peaks)

#======================================================================================
#' TODO : expand discarded_hitIndex because of insufficient overlap
#' function implementation (FIXME : come up more compatible function)

Discard.peak_0 <- mapply(extractList, all.peakFile, drop_list)
message("peaks are discarded because of insufficient overlap")
Discard.peak_0 <- lapply(Discard.peak_0, function(elm) {
  out <- unlist(elm)
  out
})

Discard.peak_0 <- lapply(Discard.peak_0, function(ele_) {
  res <- ele_[!duplicated(ele_),]
  res
})

#======================================================================================
#' TODO: retrieve pvalue of Kept peaks
#' function implementation

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

# testme:
pvl.list <- mapply(.get.pvalue, keep_list, all.peakFile)

#======================================================================================
#' TODO :
#' function implementation

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
#=====================================================================================
#' TODO : Fisher method
#' function implementation

get.fisherScore <- function(p.list, verbose=FALSE, ...) {
  # input param checking
  Fisher.score <- suppressWarnings(
    out <- apply(p.list[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  return(Fisher.score)
}

# testme:
comb.pval <- get.fisherScore(list.pval)

#=====================================================================================
saveMe <- sapply(comb.pval, function(x) x <= 1E-08)
confirmed_idx <- lapply(keep_list, function(x) x[saveMe])
Discarded_idx <- lapply(keep_list, function(x) x[!saveMe])

#=====================================================================================
#' TODO : to get hit index for all confirmed peaks
#' function implementation (still trivial, FIXME)

confirmed.peaks <- mapply(extractList, all.peakFile, confirmed_idx)
confirmed.peaks <- mapply(unlist, confirmed.peaks)
confirmed.peaks <- lapply(confirmed.peaks, function(ele_) {
  out <- ele_[!duplicated(ele_),]
  out
})

#======================================================================================
#' TODO : to get hit index for second discarded set because of failed Fisher method
#' function implementation (FIXME : implement formal function)
Discard.peaks_1 <- mapply(extractList, all.peakFile, Discarded_idx)
Discard.peaks_1 <- lapply(Discard.peaks_1, function(ele_) {
  ans <- unlist(ele_)
  ans
})

Discard.peaks_1 <- lapply(Discard.peaks_1, function(elm) {
  ans <- elm[!duplicated(elm),]
  ans
})

#======================================================================================
#' TODO : combine initial discarded peaks and second discarded peaks
#' function implementation



