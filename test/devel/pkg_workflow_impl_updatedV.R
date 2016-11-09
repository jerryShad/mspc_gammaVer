## MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
##
## Workflow of MSPC Package for unit testing all custom function

#-------------------------------------------------------------------------
# TODO : read peak files
readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  if (verbose)
    cat(">> reading all peakfiles from given folder...\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
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

# testme:
myDat <- readPeakFiles(peakFolder = "data/")

#------------------------------------------------------------------------
# TODO : remove all background noise from every peakfiles simultenously

.denoise_peakFiles <- function(peakset, denoise_threshold=1E-4, verbose=FALSE, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  stopifnot(is.numeric(denoise_threshold))
  if (verbose) {
    cat(">> filter out background noise peaks from all replicates simultanously...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  all.peaks <- lapply(peakset, function(ele_) {
    if(is.null(ele_$p.value)) {
      peaks <- .pvalueConversion(ele_, pvalueBase = 1L)
    } else {
      keep <- subset(ele_, ele_$p.value < denoise_threshold)
      keep
    }
  })
  res <- all.peaks
  return(res)
}

# testme:
all.peakFile <- .denoise_peakFiles(peakFolder = myData, denoise_threshold = 1e-4)

#-----------------------------------------------------------------------------------------
# TODO : convert pvalue and add new metadata column to each purified peakfile

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

#-----------------------------------------------------------------------------------------
# TODO : peakOverlapping

.peakOverlapping <- function(peakset, cur.idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(cur.idx))
  FUN = match.arg(FUN)   #FIXME
  # set up the entry
  chosen <- peakset[[cur.idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(set[- cur.idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[cur.idx]),names(peakset[- cur.idx]))
  return(res)
}

# testme:
test_1 <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.min)
test_2 <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.min)
test_3 <- .peakOverlapping(peakset = all.peakFile, idx = 3L, FUN = which.min)

#-----------------------------------------------------------------------------------
# TODO : combine all test result into one list for further procedure;
# FIXME: implement small utility function that take care of this step elegantly


#-----------------------------------------------------------------------------------
# TODO : get final hitlist for all peakfiles

.get.finalHitList <- function(hitlist, verbose=FALSE, ...) {
  # input param checking
  get.vectList <- Map(.hit2Vect, hitlist)
  hitDF <- lapply(get.vectList, function(ele_) {
    out <- as.data.frame(ele_)
    out
  })
  message("combine all overlap hitlist into one")
  combine.hitList <- as.matrix(bind_rows(hitDF) %>% distinct)
  final.hitList <- apply(combine.hitList, 2, function(ele_) {
    ele_ <- as(ele_, "IntegerList")
    ele_[all(ele_==0L)] <- IntegerList(integer(0))
    ele_
  })
  res <- final.hitList
  return(res)
}

#' @title .hit2Vect

.hit2Vect <- function(ovHit, ...) {
  # input param checking
  if(verbose) {
    cat(">> loading your hitlist...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  message("reverse CompressedList as Integer vector")
  res <- lapply(ovHit, function(ele_) {
    out <- drop(ele_)
    out[is.na(out)] <- 0L
    out
  })
  result <- do.call("cbind", res)
  return(result)
}

# testme:

#--------------------------------------------------------------------------------------
# TODO : get kardinality of overlapped peaks along all peakfiles

.cnt.ovnum <- function(ovHit, verbose=FALSE) {
  if (verbose) {
    cat(">> getting kardinality for all overlapped peaks from replicates...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  tot.num <- Reduce('+', lapply(ovHit, lengths))
  res <- tot.num
  return(res)
}

#--------------------------------------------------------------------------------------
# TODO: keep peaks that meet initial filtering criteria

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
  res <- out
  return(res)
}

#testme :

##==========================================================================================
#-------------------------------------------------------------------------------------------
# TODO : perform fisher method and get Fisher global score
.get.pvalue <- function(ovHit, gr, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(gr)=="GRanges")
  # if(!is(ovHit[[1]], "CompressedIntegerList")) {
  #   stop("invalid input, entry must be IntegerList")
  # }
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


#==============================================================================================
#----------------------------------------------------------------------------------------------

.confirmed.peaks <- function(peak.list, comb.stringThreshold=1E-8, verbose=TRUE) {
  # check input param
  stopifnot(is.numeric(comb.stringThreshold))
  stopifnot(inherits(expanded.peaks[[1]], c("GRanges", "data.frame")))
  if(!"comb.pvalue" %in% colnames(mcols(expanded.peaks[[1]]))) {
    stop("combined pvalue of peaks are not found, unable to proceed")
  }
  res <- lapply(expanded.peaks, function(ele_) {
    ans <- subset(ele_, ele_$comb.pvalue <= comb.stringThreshold)
    ans <- ans[!duplicated(ans),]
    ans
  })
  return(res)
}







