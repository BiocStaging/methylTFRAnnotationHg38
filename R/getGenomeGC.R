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
#' # 1. Create a dummy object and save it to a temp directory
#' mock_dir <- tempdir()
#' mock_file <- file.path(mock_dir, "genomewide_GC_hg38.rds")
#'
#' # Mocking a GRanges object with the expected metadata columns
#' mock_gr <- GenomicRanges::GRanges("chr1:1-100")
#' mock_gr$GC_bias <- 0.5
#' mock_gr$GC_bin <- 1L
#' saveRDS(mock_gr, mock_file)
#'
#' # 2. Point the local directory option to the temp directory
#' old_opt <- getOption("methylTFRAnnotationHg38.datadir")
#' options(methylTFRAnnotationHg38.datadir = mock_dir)
#'
#' # 3. Run the function
#' gc <- getGenomeGC()
#'
#' # 4. Restore options and clean up
#' options(methylTFRAnnotationHg38.datadir = old_opt)
#' file.remove(mock_file)
#'
#' @export
getGenomeGC <- function(assembly = .ASSEMBLY) {
    assembly <- tolower(assembly)
    .resolve_resource(paste0("genomewide_GC_", assembly, ".rds"))
}
