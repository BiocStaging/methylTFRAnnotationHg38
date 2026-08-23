#' @title getGCfreq
#' @description Load the motif GC frequency table for a motif set.
#' Downloaded from AnnotationHub on first use and cached locally thereafter.
#' @param motifSet Motif set to load. One of: altius, cisbpv2, jaspar2020,
#' jaspar2020_distal.
#' @return A named \code{list} of five-row matrices, one per motif.
#' @examples
#' # 1. Create a dummy object and save it to a temporary directory
#' # to bypass AnnotationHub during package checks.
#' mock_dir <- tempdir()
#' mock_file <- file.path(mock_dir, "altius_motif_gcfreq.rds")
#'
#' # (Mocking a 5-row matrix as expected by the function output)
#' mock_data <- list(mock_motif_1 = matrix(0.25, nrow = 5, ncol = 4))
#' saveRDS(mock_data, mock_file)
#'
#' # 2. Temporarily point the package's local directory to tempdir()
#' old_opt <- getOption("methylTFRAnnotationHg38.datadir")
#' options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' # 3. Run the function (it will now read the mock file)
#' x <- getGCfreq(motifSet = "altius")
#'
#' # 4. Restore options and clean up the temporary file
#' options(methylTFRAnnotationHg38.datadir = old_opt)
#' file.remove(mock_file)
#'
#' @export
getGCfreq <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    .resolve_resource(paste0(motifSet, "_motif_gcfreq.rds"))
}
