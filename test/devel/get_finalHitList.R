# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling
#
#' @title .get.finalHitList
#' @param hitlist list of overlap hit
#' @return IntegerList
#' @export
#' @importFrom BioGenerics Map
#' @importFrom dplyr bind_rows
#' @author Julaiti Shayiding
#' @example

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


#'
#' @title .hit2Vect
#' @param ovHit list of overlap hit
#' @return Integer vector
#' @export
#' @importFrom IRanges drop
#' @author Julaiti Shayiding
#' @example

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



