# Bioconductor Package for Multiple Sample Peak Calling
#
#' @title readPeakFiles
#' @param peakFolder set of bed files to be read as GRanges objects
#' @return GRanges object
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom stats setNames
#' @importFrom tools file_path_sans_ext
#' @description
#' This function is served as bootstraper to load set of enriched regions in BED files format,
#' to be available as valid S4 object for efficiently using Bioconductor Software
#' @author  Julaiti Shayiding
#' @example
## myInput <- readPeakFiles(peakFolder = "data/", verbose = FALSE)

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
