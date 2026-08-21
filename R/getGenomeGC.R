#' @title getGenomeGC
#' @description Load the genome-wide GC distribution.
#' Downloaded from AnnotationHub on first use and cached locally
#' thereafter.
#' @param assembly Genome assembly. Defaults to the assembly this
#' package was built for.
#' @return A \code{GRanges} with GC_bias and GC_bin metadata columns.
#' @import GenomicRanges
#' @importFrom utils packageVersion
#' @examples
#' \dontrun{
#' gc <- getGenomeGC()
#' }
#' @export
getGenomeGC <- function(assembly = .ASSEMBLY) {
    assembly <- tolower(assembly)
    .resolve_resource(paste0("genomewide_GC_", assembly, ".rds"))
}
