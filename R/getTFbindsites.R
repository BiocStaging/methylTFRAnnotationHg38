#' @title getTFbindsites
#' @description Retrieve transcription factor binding sites for a motif set.
#' Downloaded from AnnotationHub on first use and cached locally thereafter.
#' @param motifSet Motif set to load. One of: altius, cisbpv2, jaspar2020,
#' jaspar2020_distal.
#' @return A \code{GRangesList}, one element per motif.
#' @import GenomicRanges
#' @examples
#' # 1. Create a dummy object and save it to a temp directory
#' mock_dir <- tempdir()
#' mock_file <- file.path(mock_dir, "altius_tf_bindsites.rds")
#'
#' # Mocking a GRangesList as expected by the function output
#' mock_grl <- GenomicRanges::GRangesList(
#'     mock_motif_1 = GenomicRanges::GRanges("chr1:100-200")
#' )
#' saveRDS(mock_grl, mock_file)
#'
#' # 2. Point the local directory option to the temp directory
#' old_opt <- getOption("methylTFRAnnotationHg38.datadir")
#' options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' # 3. Run the function
#' x <- getTFbindsites()
#'
#' # 4. Restore options and clean up
#' options(methylTFRAnnotationHg38.datadir = old_opt)
#' file.remove(mock_file)
#'
#' @export
getTFbindsites <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    .resolve_resource(paste0(motifSet, "_tf_bindsites.rds"))
}
