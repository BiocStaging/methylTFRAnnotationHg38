.PKG_NAME <- "methylTFRAnnotationHg38"
.MOTIF_SETS <- c("altius", "cisbpv2", "jaspar2020", "jaspar2020_distal")
.ASSEMBLY <- "hg38"

#' @keywords internal
.local_dir <- function() {
    d <- getOption(
        "methylTFRAnnotationHg38.datadir",
        Sys.getenv("METHYL_TFRANNOTATION_HG38_DIR", "")
    )
    if (nzchar(d)) d else NULL
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
    hub <- AnnotationHub::AnnotationHub()
    hits <- AnnotationHub::query(hub, .PKG_NAME)
    idx <- match(file, hits$title)
    if (is.na(idx)) {
        stop(
            "Resource not found on AnnotationHub: ", file,
            "\nAvailable: ", paste(hits$title, collapse = ", ")
        )
    }
    hits[[names(hits)[idx]]]
}

#' @keywords internal
#' @noRd
.check_motif_set <- function(motifSet) {
    motifSet <- tolower(motifSet)
    if (length(motifSet) != 1 || !motifSet %in% .MOTIF_SETS) {
        stop(
            "Invalid motif set. Available: ",
            paste(.MOTIF_SETS, collapse = ", ")
        )
    }
    motifSet
}
