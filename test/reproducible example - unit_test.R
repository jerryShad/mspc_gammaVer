## code robut function based on following workflow:

aa <- GRanges( seqnames=Rle("chr1", 3),ranges=IRanges(c(2,7,16), c(5,14,20)),
               rangeName=c("a1", "a2", "a3"), score=c(4, 6,9))

bb <- GRanges(seqnames=Rle("chr1", 3),ranges=IRanges(c(4,13,26), c(11,17,28)),
              rangeName=c("b1", "b2", "b3"), score=c(11, 7, 8))

cc <- GRanges(seqnames=Rle("chr1", 4),ranges=IRanges(c(1,4,10, 23), c(3,8,14, 29)),
              rangeName=c("c1", "c2", "c3", "c4"), score= c(4, 6, 3, 8))


ovFunc <- function(set, idx=1L, FUN=which.min) {
  chosen <- set[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(set[-idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(set[idx]),names(set[-idx]))
  return(res)
}

ss <- list(a=aa,b=bb,c=cc) # each replicate must be assigned to named variable.
v1 <- ovFunc(set = ss, idx = 1L, FUN = which.min)
v2 <- ovFunc(set = ss, idx = 2L, FUN = which.min)
v3 <- ovFunc(set = ss, idx = 3L, FUN = which.min)

##========================================================

getVect <- function(myList) {
  # input param checking
  res <- lapply(myList, function(ele_) {
    ans <- drop(ele_)
    ans[is.na(ans)] <- 0L
    ans
  })
  return(res)
}

resVect <- Map(getVect, list(v1, v2, v3))

##========================================================
ans <- lapply(resVect, as.data.frame)

library(dplyr)
myHit <- bind_rows(ans) %>% distinct %>% arrange(a)
HitMatt <- as.matrix(myHit)

myFinalHit <- apply(HitMatt, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

## compute total overlapped number
kd <- Reduce('+', lapply(myFinalHit, lengths))

keep_ <- kd >= 3L

init.pass <- lapply(myFinalHit, function(ele_) {
  pass <- ele_[keep_]
  pass
})

.get.pval <- function(ov.hit, obj) {
  ans <- extractList(obj$score, ov.hit)
  return(ans)
}

pval.list <- mapply(.get.pval, init.pass, ss)
