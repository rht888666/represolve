#' List frozen benchmark datasets and SHA-256 checksums
#'
#' Reads the bundled \code{DATA_MANIFEST.csv} listing the nine frozen
#' replicate-half datasets archived at Zenodo.
#'
#' @return Data frame of dataset identifiers, file names, byte sizes,
#'   and SHA-256 checksums.
list_datasets <- function() {
    path <- system.file("extdata", "DATA_MANIFEST.csv", package = "represolve")
    if (!nzchar(path) || !file.exists(path)) {
        stop("DATA_MANIFEST.csv not found in installed represolve package")
    }
    utils::read.csv(path, stringsAsFactors = FALSE)
}
