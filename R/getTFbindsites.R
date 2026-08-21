#' @title getTFbindsites
#' @description Retrieve transcription factor binding sites for a motif set. Downloaded from AnnotationHub on first use and cached locally thereafter.
#' @param motifSet Motif set to load. One of: altius, cisbpv2, jaspar2020, jaspar2020_distal.
#' @return A \code{GRangesList}, one element per motif.
#' @examples
#' \donttest{
#' x <- getTFbindsites()
#' }
#' @export
getTFbindsites <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    .resolve_resource(paste0(motifSet, "_tf_bindsites.rds"))
}
