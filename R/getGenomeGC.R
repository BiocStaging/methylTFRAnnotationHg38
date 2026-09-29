#' @title Load the genome-wide GC distribution
#' @description Load the genome-wide GC distribution for hg38.
#' Downloaded from AnnotationHub on first use and cached locally
#' thereafter.
#' @details The genome is tiled into non-overlapping 30 nt windows. Each
#' window carries its GC fraction and its genome-wide GC quintile.
#' \pkg{methylTFR} uses the table to assign every methylation call to
#' a GC bin. See \code{\link{genomewide_GC}} for the data format.
#'
#' The package serves a single assembly, so the argument can be left
#' out. It is kept so that code written as
#' \code{getGenomeGC("hg38")} works; any other assembly is an error.
#' @param assembly Genome assembly. Only \code{"hg38"} (the default)
#' is available in this package; case-insensitive.
#' @return A \code{GRanges} with metadata columns \code{GC_bias} (GC
#' fraction of the window) and \code{GC_bin} (quintile, 1 to 5).
#' @seealso \code{\link{genomewide_GC}}, \code{\link{getGCfreq}}
#' @examples
#' # Point the package at a temporary directory holding a small mock
#' # object so the example runs without downloading from AnnotationHub.
#' mock_dir <- file.path(tempdir(), "genomegc_example")
#' dir.create(mock_dir, showWarnings = FALSE)
#' gr <- GenomicRanges::GRanges(c("chr1:1-30", "chr1:31-60", "chr1:61-90"))
#' gr$GC_bias <- c(0.30, 0.43, 0.63)
#' gr$GC_bin <- c(1L, 3L, 5L)
#' saveRDS(gr, file.path(mock_dir, "genomewide_GC_hg38.rds"))
#' old_opt <- options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' getGenomeGC()
#'
#' options(old_opt)
#' unlink(mock_dir, recursive = TRUE)
#' @export
getGenomeGC <- function(assembly = "hg38") {
    if (!is.character(assembly) || length(assembly) != 1L ||
        is.na(assembly) || tolower(assembly) != .ASSEMBLY) {
        stop(
            .PKG_NAME, " provides the GC distribution for ", .ASSEMBLY,
            " only. Call getGenomeGC() without arguments, or load the ",
            "methylTFR annotation package for the assembly you need."
        )
    }
    .resolve_resource(paste0("genomewide_GC_", .ASSEMBLY, ".rds"))
}
