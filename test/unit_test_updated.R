# MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
#
#' TODO : start with mini example

#--------------------------------------------------------------------------
a <- GRanges( seqnames=Rle("chr1", 4),ranges=IRanges(c(3,33,54,91), c(26,42,71,107)),
               rangeName=c("a1", "a2", "a3", "a4"), score=c(22, 6,13, 7))

b <- GRanges(seqnames=Rle("chr1", 6),ranges=IRanges(c(5,12,16,21,37,78), c(9,14,19,29,45,84)),
              rangeName=c("b1", "b2", "b3", "b4", "b5", "b6"),
             score=c(3, 5, 3, 9, 4, 3))

c <- GRanges(seqnames=Rle("chr1", 7),ranges=IRanges(c(1,8,18,35,42,59,81), c(6,13,27,40,46,63,114)),
              rangeName=c("c1", "c2", "c3", "c4","c5","c6","c7"),
             score= c(2.1, 3, 5.1, 3.5, 7, 2, 10))


#--------------------------------------------------------------------------
peak.list <- list(f1=a, f2=b,f3=c)
peak.list <- lapply(peak.list, function(ele_) {
  ans <- .pvalueConversion(ele_, 1L)
  ans
})

# test peakOverlapping function

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
    out.idx0 <- as(FUN(extractList(ele_$p.value, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}
#-----------------------------------------------------------------------------

test_1 <- .peakOverlapping(peakset = peak.list, idx = 1L, FUN = which.min)
test_2 <- .peakOverlapping(peakset = peak.list, idx = 2L, FUN = which.min)
test_3 <- .peakOverlapping(peakset = peak.list, idx = 3L, FUN = which.min)

hit.list <- list(test.1=test_1, test.2=test_2, test.3=test_3)

#---------------------------------------------------------------------
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

res.vec <- Map(Hit2IntVec, hit.list)
.hitDF <- lapply(res.vec, as.data.frame)


library(dplyr)
hit.tmp <- bind_rows(.hitDF) %>% distinct
hitMT <- as.matrix(hit.tmp)

final_hit <- apply(hitMT, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

#------------------------------------------------------------------------
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

#-------------------------------------------------------------------------
.filtBy_suffOvHit <- function(ovHit, peakset, replicate.type=c("Biological", "Technical"), ...) {
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
    ans <- list(keep=ele_[.keep], discard=ele_[!.keep])
    ans
  })
  return(out)
}

#---------------------------------------------------------------------
filt.res <- .filtBy_suffOvHit(final_hit, peak.list, replicate.type = "Technical")


####=================================================================================


#####################################################################################
#-----------------------------------------------------------------------------------#
myList <- list(f1=IntegerList(1,2,3,4,1,1,1,integer(0),1,2,4),
               f2=IntegerList(1,5,integer(0),integer(0),2,3,4,6,1,5,6),
               f3=IntegerList(1,4,6,7,2,3,3,7,2,5,7))

len <- Reduce('+', lapply(final_hit, lengths))

keepMe <- len >= length(peak.list)

res.filt <- lapply(final_hit, function(elm) {
  ans <- list(
    keep_=elm[keepMe],
    droped_=elm[!keepMe])
  ans
})

library(purrr)
map(final_hit, lengths) %>%
  reduce(`+`) %>%
  map_lgl(`>=`, length(final_hit)) -> keep_me

keep_list <- map(final_hit, ~.[keep_me])
drop_list <- map(final_hit, ~.[!keep_me])

discarded_peaks <- mapply(extractList, peak.list, drop_list)
allDiscarded_0 <- lapply(discarded_peaks, function(ele_) {
  ans <- unlist(ele_)
  ans
})

