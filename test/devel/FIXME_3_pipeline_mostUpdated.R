# MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
#
#' @description Brandnew pipeline for MSPC Package
#' @details
#' FIX the bug and to get identical result with Musera tools

#-----------------------------------------------------------------------------------------------
#' @description read bed files as GRanges objects
#' This is the implementation
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

#' @example
myData <- readPeakFiles(peakFolder = "data/")

#----------------------------------------------------------------------------------------------
#' @description add pvalue as new metadata of GRanges objects
#' This is function implementation

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

#---------------------------------------------------------------------------------------------
#' @description remove background noise from all GRanges simultanously
#' This is the implementation

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
  res <- peaks_rmNoi
  return(res)
}

#
total.ERs <- .denoise_peakFiles(peakFolder = myData, tau.w = 1.0E-04)

#--------------------------------------------------------------------------------------------
#' @description peak overlapping over GRanges object simultanously
#' This is function implementation

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

#' @example

.hit_1 <- .peakOverlapping(peakset = total.ERs, idx = 1L, FUN = which.max)
.hit_2 <- .peakOverlapping(peakset = total.ERs, idx = 2L, FUN = which.max)
.hit_3 <- .peakOverlapping(peakset = total.ERs, idx = 3L, FUN = which.max)


#------------------------------------------------------------------------------------------
idx <- sort(names(.hit_1))
.hit_2 <- .hit_2[idx]
.hit_3 <- .hit_3[idx]

#'-----------------------------------------------------------------------------------------
#' TODO : try this manipulation for .hit2, .hit3
idx <- sort(names(.hit_1))
.hit_1 <- DataFrame(.hit_1)
.hit_2 <- DataFrame(.hit_2[idx])
.hit_3 <- DataFrame(.hit_3[idx])

Hit <- Map(rbind, .hit_1, .hit_2, .hit_3)
#------------------------------------------------------------------------------------------
#' @description parallel filtering
#' This is list of function

.filtByIndx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  ovNum <- Reduce('+', lapply(ovHit, lengths))
  keepMe <- sapply(ovNum, function(elm) {
    res <- elm >= min.c
    res
  })
  return(keepMe)
}

#' @example
keep_me <- .filtByIndx(peakset = total.ERs, .hit_3,
                       replicate.type = "Biological")

library(purrr)
keepList <- map(.hit_1, ~.[keep_me])
dropList <- map(.hit_1, ~.[!keep_me])


#keep_list <- .keep_hit(.hit_1, total.ERs, "Biological")

Discard.peaks_0 <- mapply(extractList, total.ERs, dropList)
Discard.peaks_0 <- Map(unlist, Discard.peaks_0)
rm(dropList,keep_me)

#---------------------------------------------------------------------------------------
.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

# testme:
pvl.list <- mapply(.get.pvalue, keepList, total.ERs)

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

list.pval <- data.frame(mapply(.helper.getPvalue, pvl.list))
rm(pvl.list)

# get.fisherScore <- function(p.list, ...) {
#   # input param checking
#   require(metap)
#   Fisher.score <- suppressWarnings(
#     out <- apply(p.list[,], 1, function(ro) {
#       ans <- sumlog(ro)$p
#     })
#   )
#   return(Fisher.score)
# }

#--------------------------------------------------------------------------------------
#' @description This is new method for getting combined pvalue

stats    <- rowSums( -2*log( list.pval ), na.rm=T )
Npresent <- rowSums( !is.na(list.pval) )
comb.pval <- pchisq( stats, df=2*Npresent, lower=FALSE )

#--------------------------------------------------------------------------------------

#comb.pval <- get.fisherScore(list.pval)
rm(list.pval)

saveMe <- sapply(comb.pval, function(x) x <= 1.0E-08)
Confirmed_idx <- lapply(keepList, function(x) x[saveMe])
Discarded_idx <- lapply(keepList, function(x) x[!saveMe])
rm(comb.pval)

confirmed <- mapply(extractList, total.ERs, Confirmed_idx)
all_Confirmed.3 <- mapply(unlist, confirmed)
all_Confirmed.3 <- lapply(all_Confirmed.3, function(ele_) {
  out <- ele_[!duplicated(ele_),]
  out
})
rm(confirmed, Confirmed_idx,keepList)

Discard.peaks_1 <- mapply(extractList, total.ERs, Discarded_idx)
Discard.peaks_1 <- mapply(unlist, Discard.peaks_1)
Discard.peaks_1 <- lapply(Discard.peaks_1, function(elm) {
  ans <- elm[!duplicated(elm),]
  ans
})
rm(Discarded_idx)

all_Discarded.3 <- suppressWarnings(mapply(c, Discard.peaks_0, Discard.peaks_1))
#====================================================================================
#' running first test, I get following results:

all_Confirmed.1
all_Discarded.1

all_Confirmed.2
all_Discarded.2

all_Confirmed.3
all_Discarded.3

#------------------------------------------------------------------------------------
.Confirmed.1 <- lapply(all_Confirmed.1, "data.frame")
.Confirmed.2 <- lapply(all_Confirmed.2, "data.frame")
.Confirmed.3 <- lapply(all_Confirmed.3, "data.frame")


nm1 <- sort(names(.Confirmed.1))
confirmed_ERs <- base::Map(function(x,y,z) unique(rbind(x,y,z)), .Confirmed.1, .Confirmed.2[nm1], .Confirmed.3[nm1])
confirmed_ERs <- lapply(confirmed_ERs, "GRanges")
confirmed_ERs <- lapply(confirmed_ERs, function(elm) {
  res <- elm[!duplicated(elm),]
  res
})

.Discard.1 <- lapply(all_Discarded.1, "data.frame")
.Discard.2 <- lapply(all_Discarded.2, "data.frame")
.Discard.3 <- lapply(all_Discarded.3, "data.frame")

nm1 <- sort(names(.Discard.1))
Discard_ERs <- base::Map(function(x,y,z) unique(rbind(x,y,z)), .Discard.1, .Discard.2[nm1], .Discard.3[nm1])
Discard_ERs <- lapply(Discard_ERs, "GRanges")

Discard_ERs <- lapply(Discard_ERs, function(elm) {
  res <- elm[!duplicated(elm),]
  res
})

#confirmed_ERs <- Map(rbind, all_Confirmed.1, all_Confirmed.2, all_Confirmed.3)
#------------------------------------------------------------------------------------
nm1 <- sort(names(myList_1))
BiocGenerics::Map(function(x,y,z)
  unique(rbind(all_Confirmed.1,all_Confirmed.2,all_Confirmed.3)), all_Confirmed.1, all_Confirmed.2[nm1], all_Confirmed.3[nm1])


#------------------------------------------------------------------------------------

func <- function(list, ...) {
  res <- lapply(list, function(elm) {
    out <- as(elm, "data.frame")
    out
  })
  return(res)
}

conf.1 <- func(all_Confirmed.1)
conf.2 <- func(all_Confirmed.2)
conf.3 <- func(all_Confirmed.3)

nm1 <- sort(names(conf.1))
confirmed_peaks <- Map(function(x,y,z) unique(rbind(x,y,z)), conf.1, conf.2[nm1], conf.3[nm1])

confirmed_peaks <- lapply(confirmed_peaks, function(ele_) {
  res <- as(ele_, "GRanges")
  res
})
