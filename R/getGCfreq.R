#' @title Load motif GC frequency tables
#' @description Load the motif GC frequency tables for a motif set.
#' Downloaded from AnnotationHub on first use and cached locally
#' thereafter.
#' @details Each table describes the GC composition around one motif's
#' binding sites. Column \eqn{j} corresponds to the \eqn{j}-th 30 nt
#' window along the (extended) binding site; row \eqn{i} gives the
#' fraction of binding sites whose window at that position falls into
#' genome-wide GC quintile \eqn{i} (row 1 lowest, row 5 highest GC).
#' Columns therefore sum to one. \pkg{methylTFR} combines these
#' fractions with the observed methylation per GC bin to obtain the
#' methylation a motif would be expected to show from sequence
#' composition alone. See \code{\link{motif_gcfreq}} for the data
#' format.
#' @param motifSet Motif set to load. One of
#' \code{availableMotifSets()}: \code{altius}, \code{cisbpv2},
#' \code{jaspar2020}, \code{jaspar2020_distal}. Case-insensitive.
#' @return A named \code{list} of numeric matrices with five rows (GC
#' quintiles), one matrix per motif.
#' @seealso \code{\link{motif_gcfreq}}, \code{\link{getTFbindsites}},
#' \code{\link{getGenomeGC}}
#' @examples
#' # Point the package at a temporary directory holding a small mock
#' # table so the example runs without downloading from AnnotationHub.
#' mock_dir <- file.path(tempdir(), "gcfreq_example")
#' dir.create(mock_dir, showWarnings = FALSE)
#' freq <- matrix(c(0.1, 0.2, 0.4, 0.2, 0.1), nrow = 5, ncol = 3)
#' saveRDS(list(MA0001.1 = freq), file.path(mock_dir,
#' "altius_motif_gcfreq.rds"))
#' old_opt <- options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' gcfreqs <- getGCfreq("altius")
#' gcfreqs$MA0001.1
#' colSums(gcfreqs$MA0001.1)
#'
#' options(old_opt)
#' unlink(mock_dir, recursive = TRUE)
#' @export
getGCfreq <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    .resolve_resource(paste0(motifSet, "_motif_gcfreq.rds"))
}
