.PKG_NAME <- "methylTFRAnnotationHg38"
.MOTIF_SETS <- c("altius", "cisbpv2", "jaspar2020", "jaspar2020_distal")
.ASSEMBLY <- "hg38"

#' @keywords internal
#' @noRd
.local_dir <- function() {
    d <- getOption(
        "methylTFRAnnotationHg38.datadir",
        Sys.getenv("METHYL_TFRANNOTATION_HG38_DIR", "")
    )
    if (is.null(d) || identical(d, "")) {
        return(NULL)
    }
    if (!is.character(d) || length(d) != 1L || is.na(d)) {
        stop(
            "The option methylTFRAnnotationHg38.datadir must be a single ",
            "character string (a directory path)."
        )
    }
    if (!dir.exists(d)) {
        stop(
            "Local annotation directory does not exist: ", d,
            "\nUnset options(methylTFRAnnotationHg38.datadir) and the ",
            "METHYL_TFRANNOTATION_HG38_DIR environment variable ",
            "to use AnnotationHub."
        )
    }
    d
}

#' @keywords internal
#' @noRd
#' @description AnnotationHub records belonging to this package, as a
#' character vector of titles named by AnnotationHub ID. Kept separate
#' from .resolve_resource() so the lookup logic can be unit-tested
#' without network access.
.hub_titles <- function() {
    hub <- AnnotationHub::AnnotationHub()
    hits <- AnnotationHub::query(hub, .PKG_NAME)
    stats::setNames(hits$title, names(hits))
}

#' @keywords internal
#' @noRd
.hub_get <- function(id) {
    AnnotationHub::AnnotationHub()[[id]]
}

#' @keywords internal
#' @noRd
#' @description Resolve one annotation resource by file name.
#' Reads from a local directory when one is configured, otherwise
#' from AnnotationHub. The local path exists so the package can be
#' exercised against a freshly built annotation before the data are
#' on the hub, and so users who have downloaded the files by hand
#' can point at them.
.resolve_resource <- function(file) {
    dir <- .local_dir()
    if (!is.null(dir)) {
        path <- file.path(dir, file)
        if (!file.exists(path)) {
            stop(
                "Annotation file not found: ", path,
                "\n(reading from a local directory because ",
                .PKG_NAME, ".datadir is set)"
            )
        }
        return(readRDS(path))
    }
    titles <- .hub_titles()
    idx <- match(file, titles)
    if (is.na(idx)) {
        stop(
            "Resource not found on AnnotationHub: ", file,
            "\nAvailable: ", paste(titles, collapse = ", ")
        )
    }
    .hub_get(names(titles)[idx])
}

#' @keywords internal
#' @noRd
.check_motif_set <- function(motifSet) {
    if (!is.character(motifSet) || length(motifSet) != 1L ||
        is.na(motifSet) || !tolower(motifSet) %in% .MOTIF_SETS) {
        stop(
            "Invalid motif set. Available: ",
            paste(.MOTIF_SETS, collapse = ", ")
        )
    }
    tolower(motifSet)
}

#' @title List the motif sets provided by this package
#' @description Returns the names accepted by the \code{motifSet}
#' argument of \code{\link{getTFbindsites}} and
#' \code{\link{getGCfreq}}.
#' @return A character vector of motif set names.
#' @seealso \code{\link{methylTFRAnnotationHg38}} for an overview of the
#' package.
#' @examples
#' availableMotifSets()
#' @export
availableMotifSets <- function() {
    .MOTIF_SETS
}
