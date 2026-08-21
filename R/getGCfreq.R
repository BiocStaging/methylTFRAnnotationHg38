#' @title getGCfreq
#' @description Load the motif GC frequency table for a motif set. Downloaded from AnnotationHub on first use and cached locally thereafter.
#' @param motifSet Motif set to load. One of: altius, cisbpv2, jaspar2020, jaspar2020_distal.
#' @return A named \code{list} of five-row matrices, one per motif.
#' @examples
#' \donttest{
#' x <- getGCfreq()
#' }
#' @export
getGCfreq <- function(motifSet = "altius") {
    motifSet <- .check_motif_set(motifSet)
    .resolve_resource(paste0(motifSet, "_motif_gcfreq.rds"))
}
