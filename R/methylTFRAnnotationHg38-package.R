#' methylTFRAnnotationHg38: methylTFR annotations for hg38
#'
#' @description
#' Genome annotations that \pkg{methylTFR} needs to compute
#' GC-bias-corrected transcription factor (TF) deviation scores from DNA
#' methylation data on the Homo sapiens hg38 assembly.
#'
#' @details
#' The package provides three kinds of resources, each with an accessor:
#' \describe{
#'   \item{\code{\link{getTFbindsites}}}{Genome-wide TF binding sites,
#'     one \code{GRanges} per motif, extended to the footprint window
#'     (see \code{\link{tf_bindsites}}).}
#'   \item{\code{\link{getGCfreq}}}{Motif GC frequency tables: how
#'     each motif's binding sites distribute over genome-wide GC quintiles
#'     (see \code{\link{motif_gcfreq}}).}
#'   \item{\code{\link{getGenomeGC}}}{The genome-wide GC distribution
#'     in 30 nt windows (see \code{\link{genomewide_GC}}).}
#' }
#' \code{\link{availableMotifSets}} lists the motif sets that can be
#' requested. The three objects are passed to
#' \code{methylTFR::run_methyltfr()}.
#'
#' The data are too large to ship with the package. They are hosted on
#' AnnotationHub, downloaded on first use and served from the local
#' AnnotationHub cache afterwards. The metadata of every resource is
#' listed in \code{system.file("extdata", "metadata.csv", package =
#' "methylTFRAnnotationHg38")}, and the scripts that produced them are in
#' \code{system.file("scripts", package = "methylTFRAnnotationHg38")}.
#'
#' To read the \code{.rds} files from a local directory instead of
#' AnnotationHub (for example on a machine without internet access), set
#' \code{options(methylTFRAnnotationHg38.datadir = "/path/to/files")} or the
#' environment variable \code{METHYL_TFRANNOTATION_HG38_DIR}.
#'
#' @seealso \code{vignette("methylTFRAnnotationHg38")}
#' @importFrom GenomicRanges GRanges
"_PACKAGE"
