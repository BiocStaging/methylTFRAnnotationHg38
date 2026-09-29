#' @title Motif GC frequency tables
#' @name motif_gcfreq
#' @aliases altius_motif_gcfreq cisbpv2_motif_gcfreq jaspar2020_motif_gcfreq
#' jaspar2020_distal_motif_gcfreq
#' @description Per-motif tables describing the GC composition around
#' each motif's binding sites on hg38. Returned by
#' \code{\link{getGCfreq}}; the AnnotationHub record titles are
#' \code{<motifSet>_motif_gcfreq.rds}.
#' @format A named \code{list} with one numeric matrix per motif (names
#' are motif IDs, matching the names of the corresponding
#' \code{\link{tf_bindsites}} list). Each matrix has
#' \describe{
#'   \item{5 rows}{Genome-wide GC quintiles, from lowest (row 1) to
#'     highest (row 5) GC. Quintile boundaries are the ones stored with
#'     \code{\link{genomewide_GC}}, so both objects share one scale.}
#'   \item{one column per window position}{A 30 nt window slides one
#'     base at a time along the binding site (motif plus 200 bp on either
#'     side, padded by a further 65 bp so edge windows are complete).
#'     Column \eqn{j} is the \eqn{j}-th window.}
#' }
#' Entries are the fraction of the motif's binding sites whose window at
#' that position falls into each quintile; every column sums to one.
#' Windows containing \code{N} are ignored.
#' @details The tables are what make the methylTFR deviation score
#' bias-corrected. A motif whose sites sit in GC-rich sequence overlaps
#' hypomethylated CpG islands more often; combining these fractions with
#' the observed methylation per GC bin gives the methylation expected from
#' sequence composition alone, which is then subtracted.
#'
#' Available motif sets:
#' \describe{
#'   \item{\code{jaspar2020}}{JASPAR2020 CORE matrices for
#'     \emph{Homo sapiens} (629 PWMs), scanned genome-wide.}
#'   \item{\code{cisbpv2}}{CIS-BP v2 human PWMs
#'     (\code{human_pwms_v2} from \pkg{chromVARmotifs}), scanned
#'     genome-wide.}
#'   \item{\code{altius}}{Vierstra motif archetypes v1.0 (286
#'     archetype clusters), taken from the published hg38 archetype
#'     scan and de-duplicated by motif score.}
#'   \item{\code{jaspar2020_distal}}{A GC frequency table only: as
#'     \code{jaspar2020}, but restricted to binding sites overlapping
#'     distal regulatory regions (Ensembl Regulatory Build v104). Use it
#'     together with the unrestricted \code{jaspar2020} binding sites.}
#' }
#' @source Built with methylTFRAnnotationBuilder
#' (\url{https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder})
#' from BSgenome sequences; see
#' \code{system.file("scripts", "make-data.R", package =
#' "methylTFRAnnotationHg38")}
#' and \code{system.file("extdata", "metadata.csv", package =
#' "methylTFRAnnotationHg38")} for the exact sources and versions.
#' @seealso \code{\link{getGCfreq}}, \code{\link{methylTFRAnnotationHg38}}
#' @examples
#' ## On first use this downloads the resource from AnnotationHub:
#' \dontrun{
#' gcfreqs <- getGCfreq("jaspar2020")
#' length(gcfreqs)          # number of motifs
#' dim(gcfreqs[[1]])        # 5 x number of window positions
#' colSums(gcfreqs[[1]])    # all 1
#' }
#' @keywords datasets
NULL

#' @title Transcription factor binding sites
#' @name tf_bindsites
#' @aliases altius_tf_bindsites cisbpv2_tf_bindsites jaspar2020_tf_bindsites
#' @description Genome-wide predicted transcription factor binding sites
#' on hg38. Returned by \code{\link{getTFbindsites}}; the
#' AnnotationHub record titles are \code{<motifSet>_tf_bindsites.rds}.
#' @format A \code{GRangesList} with one \code{GRanges} per motif,
#' named by motif ID. Each range is a motif match on the primary
#' chromosomes (chr1-22, X, Y; chrM excluded) extended by 200 bp
#' on either side of the match, so all ranges of a motif have the same
#' width (motif width + 400 bp). The strand of the match is kept; there
#' are no metadata columns. Coordinates are 1-based.
#' @details Available motif sets:
#' \describe{
#'   \item{\code{jaspar2020}}{JASPAR2020 CORE matrices for
#'     \emph{Homo sapiens} (629 PWMs), scanned genome-wide.}
#'   \item{\code{cisbpv2}}{CIS-BP v2 human PWMs
#'     (\code{human_pwms_v2} from \pkg{chromVARmotifs}), scanned
#'     genome-wide.}
#'   \item{\code{altius}}{Vierstra motif archetypes v1.0 (286
#'     archetype clusters), taken from the published hg38 archetype
#'     scan and de-duplicated by motif score.}
#'   \item{\code{jaspar2020_distal}}{A GC frequency table only: as
#'     \code{jaspar2020}, but restricted to binding sites overlapping
#'     distal regulatory regions (Ensembl Regulatory Build v104). Use it
#'     together with the unrestricted \code{jaspar2020} binding sites.}
#' }
#' @source See \code{\link{motif_gcfreq}}.
#' @seealso \code{\link{getTFbindsites}}, \code{\link{methylTFRAnnotationHg38}}
#' @examples
#' ## On first use this downloads the resource from AnnotationHub:
#' \dontrun{
#' tfbs <- getTFbindsites("jaspar2020")
#' length(tfbs)                 # number of motifs
#' head(lengths(tfbs))          # binding sites per motif
#' tfbs[[1]]
#' }
#' @keywords datasets
NULL

#' @title Genome-wide GC distribution
#' @name genomewide_GC
#' @aliases genomewide_GC_hg38
#' @description GC content of the hg38 genome in non-overlapping 30 nt
#' windows. Returned by \code{\link{getGenomeGC}}; the AnnotationHub
#' record title is \code{genomewide_GC_hg38.rds}.
#' @format A \code{GRanges}, one range per 30 nt window on the primary
#' chromosomes, with metadata columns
#' \describe{
#'   \item{\code{GC_bias}}{Numeric, GC fraction of the window (0 to 1);
#'     \code{NaN} for windows made of \code{N}.}
#'   \item{\code{GC_bin}}{Integer, genome-wide GC quintile of the
#'     window, 1 (lowest) to 5 (highest).}
#' }
#' and \code{S4Vectors::metadata()} entries \code{gc_breaks} (the six
#' quintile boundaries), \code{bin_scope} (\code{"genome"}: quantiles
#' taken genome-wide) and \code{step} (window step, 30).
#' @details \pkg{methylTFR} overlaps each methylation call with these
#' windows to assign it a GC bin; the motif tables in
#' \code{\link{motif_gcfreq}} use the same boundaries.
#' @source Computed from the UCSC hg38 BSgenome package; see
#' \code{system.file("scripts", "make-data.R", package =
#' "methylTFRAnnotationHg38")}.
#' @seealso \code{\link{getGenomeGC}}, \code{\link{methylTFRAnnotationHg38}}
#' @examples
#' ## On first use this downloads the resource from AnnotationHub:
#' \dontrun{
#' gc <- getGenomeGC()
#' gc
#' table(gc$GC_bin)
#' S4Vectors::metadata(gc)$gc_breaks
#' }
#' @keywords datasets
NULL
