#' @title Load transcription factor binding sites
#' @description Load genome-wide transcription factor binding sites for a
#' motif set. Downloaded from AnnotationHub on first use and cached
#' locally thereafter.
#' @details Each element of the returned list holds the binding sites of
#' one motif. Every range is the motif match extended by 200 bp on either
#' side, the footprint window across which \pkg{methylTFR} reads
#' methylation. See \code{\link{tf_bindsites}} for the data format.
#'
#' Sets ending in \code{_distal} only have their own GC frequency
#' table; they share the binding sites of the base set. Passing a
#' \code{_distal} name here returns the base set's binding sites.
#'
#' @param motifSet Motif set to load. One of
#' \code{availableMotifSets()}: \code{altius}, \code{cisbpv2},
#' \code{jaspar2020}, \code{jaspar2020_distal}. Case-insensitive.
#' @return A \code{GRangesList} with one \code{GRanges} per motif,
#' named by motif ID.
#' @seealso \code{\link{tf_bindsites}}, \code{\link{getGCfreq}}
#' @examples
#' # Point the package at a temporary directory holding a small mock
#' # object so the example runs without downloading from AnnotationHub.
#' mock_dir <- file.path(tempdir(), "tfbs_example")
#' dir.create(mock_dir, showWarnings = FALSE)
#' mock_grl <- GenomicRanges::GRangesList(
#'     MA0001.1 = GenomicRanges::GRanges(
#'         c("chr1:1001-1411", "chr2:5001-5411"),
#'         strand = c("+", "-")
#'     )
#' )
#' saveRDS(mock_grl, file.path(mock_dir, "altius_tf_bindsites.rds"))
#' old_opt <- options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' getTFbindsites("altius")
#'
#' options(old_opt)
#' unlink(mock_dir, recursive = TRUE)
#' @export
getTFbindsites <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    motifSet <- sub("_distal$", "", motifSet)
    .resolve_resource(paste0(motifSet, "_tf_bindsites.rds"))
}
