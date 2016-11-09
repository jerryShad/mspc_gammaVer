## MSPC Project -  Bioconductor Package for Multiple Sample Peak Calling
##
## code review of whole workflow before building MSPC Packages
##====================================================================
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

## test function
myData <- readPeakFiles(peakFolder = "data/")

##====================================================================
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

##====================================================================

.denoise_peakFiles <- function(peakFolder, denoise_threshold=1E-4, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks from all replicates simultanously...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(denoise_threshold))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  res <- lapply(peakFolder, function(ele_) {
    ans <- .pvalueConversion(ele_, pvalueBase=1L)
    ans
  })
  filt <- lapply(res, function(ele_) {
    out <- subset(ele_, ele_$p.value < denoise_threshold)
    out
  })
  res <- filt
  return(res)
}

## test function

all.peakFile <- .denoise_peakFiles(peakFolder = myData, denoise_threshold = 1e-4)

##======================================================================
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

##======================================================================
## NEW METHOD : all test must be completed to get final hit list before to go next analysing phase

## finish all needed test

first_ov <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.min)
second_ov <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.min)
third_ov <- .peakOverlapping(peakset = all.peakFile, idx = 3L, FUN = which.min)

hitlist <- list(test1 = first_ov, test2 = second_ov, test3 = third_ov)

##======================================================================
Hit2IntVec <- function(ovHit, verbose=FALSE, ...) {
  # input param cheching
  if(!is(ovHit[[1]], "CompressedIntegerList")) {
    stop("entry hitlist must be IntegerList")
  }
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

## test function

hitAsVect <- Map(Hit2IntVec, hitlist)
.res <- lapply(hitAsVect, as.data.frame)

library(dplyr)
final.hit <- bind_rows(.res) %>% distinct
res <- as.matrix(final.hit)

final_hit <- apply(res, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

##======================================================================
.cnt.ovnum <- function(ovHit, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  tot.num <- Reduce('+', lapply(ovHit, lengths))
  res <- tot.num
  return(res)
}

# test function

num <- .cnt.ovnum(ovHit = final_hit)

##======================================================================
.getMinovParam <- function(peakset, replicate.type=c("Biological", "Technical"), verbose=FALSE) {
  # input param check
  if(!is(peakset[[1]], "GRanges")) {
    stop("input must conform with GRanges object")
  }
  stopifnot(length(peakset)>=2)
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  return(min.c)
}

# test function
min.C <- .getMinovParam(peakset = all.peakFile, "Technical")

##============================================================================================
.keep.peaks <- function(ovHit, peakset, replicate.type=c("Biological", "Technical")) {
  nn <- .cnt.ovnum(ovHit)
  min.c <- .getMinovParam(peakset, "Biological")
  keep <- nn >= min.c
  res <- lapply(ovHit, function(ele_) {
    ans <- ele_[keep]
  })
  #message("keep peaks that sufficiently overlapped")
  # call get.Pvalue helper function to continue next ...
  return(res)
}

res.initFilt <- .keep.peaks(ovHit = final_hit, peakset = all.peakFile, replicate.type = "Technical")

##================================================================

func <- function(ov.hit, obj) {
  ans <- extractList(obj$p.value, ov.hit)
  return(ans)
}

# test
pval.list <- mapply(func, res.initFilt, all.peakFile)

.pval.helper <- function(pv.li) {
  ans <- sapply(pv.li, function(x) {
    out <- ifelse(length(x)>0, x, 0)
  })
  return(ans)
}

# test
outp <- data.frame(mapply(.pval.helper, pval.list))

##========================================================================
get.fisherScore <- function(all.pvals, ...) {
  Fisher.score <- suppressWarnings(
    out <- apply(all.pvals[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  return(Fisher.score)
}

# test function

my.fishScore <- get.fisherScore(outp)

##========================================================================
.Hit2GRanges.expand <- function(ov.hit, ov.sample, allPeak.val, ...) {
  ov.expand <- ov.sample[unlist(extractList(seq_along(ov.sample), ov.hit))]
  all.expand <- mapply(ov.expand, ov.hit, ov.sample)

  Fisher.score <- get.fisherScore(all.pvals = allPeak.val)

  corr.VecDim <- lapply(allPeak.val, function(ele_) {
    tmp_pval.df <- data.frame(cbind("p"=ele_, "comb.pvalue"=Fisher.score))
    out <- tmp_pval.df[tmp_pval.df$p != 0.000000e+00, ]
    out$p <- NULL
    out
  })
  res <- Map(cbind, all.expand, corr.VecDim)
  return(res)
}

# test function

expnad.res <- .Hit2GRanges.expand(res.initFilt, all.peakFile)

##============================================================================================

.discard.peaks_0 <- function(final_ovHit, all.peakFiles, replicate.type=c("Biological", "Technical"), idx=1L) {
  #check input param
  nn <- tot.ovnum(final_ovHit)
  min.c <- get.minOv.param(all.peakFiles, "Biological")
  disc.idx <- nn < min.c
  ans <- lapply(ov, function(ele_) {
    out <- ele_[disc.idx]
    out <- Filter(length, out)
  })
  message("peak are discarded because of insufficient overlap")
  res <- mapply(exp_fun, ans, all.peakFiles[c(idx, seq_len(length(allFiles_denoise))[-idx])])
  return(res)
}

##======================================================================

.get.pvalue <- function(final_ovHit, ov.sample, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  ov.pass <- .keep.peaks(final_ovHit, ov.sample, "Biological")
  .get.pval_helper <- function(ov.pass, ov.sample) {
    res <- extractList(ov.sample$p.value, ov.pass)
    res
  }
  allPeaks.pval <- mapply(.get.pval_helper, ov.pass,
                          ov.sample[c(idx, seq_len(length(ov.sample))[-idx])])
  .pval.helper <- function(pv.li) {
    ans <- sapply(pv.li, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
    return(ans)
  }
  res <- data.frame(mapply(.pval.helper, allPeaks.pval))
  return(res)
}

##=========================================================================

## waiting implementation : facilitate output for each replicates and perform BH correction test

func <- function(final_ovHit, peakset, pAdjustMethod="BH",
                 replicate.type=c("Biological", "Technical"),output=NULL, verbose=FALSE, ...) {
  # check input pararm

  padj <- p.adjust(p, method = pAdjustMethod, n = length(p))
}
