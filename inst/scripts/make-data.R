## How the data in this package were produced
##
## Every object was built by methylTFRAnnotationBuilder from
## BSgenome.Hsapiens.UCSC.hg38 and the JASPAR2020 motif set.
##
## The build script is reproduced below. It is the script that
## generated the files submitted to AnnotationHub.

#!/usr/bin/env Rscript

#####################################################################
# 02_build_jaspar2020.R
#
# Builds methylTFRAnnotationHg38 for the JASPAR2020 motif set in
# /scratch/icbb/igunduz/mtfr_annotation_test
#
# The expensive step, motif matching across the genome, is skipped if a
# binding-site file is already present. Copy the existing one into
# place first (see PREREQUISITES) and the script only computes the GC
# annotations.
#
# The script is resumable. Every artefact is skipped if it already
# exists, so it can be re-run after an interruption.
#
# Run:
#     cd /scratch/icbb/igunduz/mtfr_annotation_test
#     Rscript 02_build_jaspar2020.R 2>&1 | tee build.log
#
#####################################################################
#
# PREREQUISITES
#
# Copy the existing binding sites into the package's extdata folder,
# under a LOWER CASE name. The builder now normalises file names to
# lower case on both the writing and the reading side, so a file named
# JASPAR2020_tf_bindsites.rds will not be found:
#
#     mkdir -p /scratch/icbb/igunduz/mtfr_annotation_test/methylTFRAnnotationHg38/inst/extdata
#     cp /path/to/JASPAR2020_tf_bindsites.rds \
#        /scratch/icbb/igunduz/mtfr_annotation_test/methylTFRAnnotationHg38/inst/extdata/jaspar2020_tf_bindsites.rds
#
# Do NOT copy the old genomewide_GC_hg38.rds. It was built with
# per-chromosome GC bins and at every genomic position; this run
# produces genome-wide bins restricted to CpG sites, and mixing the two
# would leave the observed and expected sides of a deviation score
# binned on different scales.
#
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(GenomicRanges)
    library(logger)
    library(methylTFRAnnotationBuilder)
})

## ------------------------------------------------------------------
## Configuration
## ------------------------------------------------------------------

base_dir <- "/scratch/icbb/igunduz/mtfr_annotation_test"
assembly <- "Hg38"
motif_set <- "jaspar2020"
cores <- 24 # capped by chromosome count; more than 24 gains nothing
chunk_size <- 15 # motifs per chunk in the GC frequency step

# Optional. Restricting the GC table to windows that overlap a CpG is
# exactly equivalent to the full table for CpG methylomes, and shrinks
# it by roughly threefold. It costs an extra pass over the genome, so it
# is off by default; turn it on only if the GC table turns out to be a
# problem.
restrict_to_cpg <- FALSE

# Also build the distal-restricted GC frequency table. Requires a
# GRanges of distal regions; set to NULL to skip.
distal_rds <- "/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS"

genome <- BSgenome.Hsapiens.UCSC.hg38
pkg_dir <- file.path(base_dir, paste0("methylTFRAnnotation", assembly))
extdata <- file.path(pkg_dir, "inst", "extdata")
cpg_cache <- file.path(base_dir, "cpg_sites_hg38.rds")

## ------------------------------------------------------------------
## Pre-flight
## ------------------------------------------------------------------

if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
}
setwd(base_dir)

log_info("Base directory: {base_dir}")
free_gb <- tryCatch({
    df <- system2("df", c("-BG", shQuote(base_dir)), stdout = TRUE)
    as.numeric(sub("G", "", strsplit(trimws(df[2]), "\\s+")[[1]][4]))
}, error = function(e) NA_real_)
if (!is.na(free_gb)) {
    log_info("Free space: {free_gb} GB")
    if (free_gb < 20) {
        log_warn("Less than 20 GB free. The temp chunk files and the ",
            "binding-site object both need room.")
    }
}

build_distal <- !is.null(distal_rds) && file.exists(distal_rds)
if (!is.null(distal_rds) && !build_distal) {
    log_warn("Distal regions not found at {distal_rds}; skipping the distal table.")
}

motif_sets_declared <- if (build_distal) {
    c(motif_set, paste0(motif_set, "_distal"))
} else {
    motif_set
}

## ------------------------------------------------------------------
## 1. Package scaffold
## ------------------------------------------------------------------

if (!dir.exists(pkg_dir)) {
    log_info("Creating package scaffold for {assembly} ...")
    createMethylTFRPackageScaffold(
        assembly = assembly,
        dest = base_dir,
        motifSets = motif_sets_declared
    )
} else {
    log_info("Scaffold already exists at {pkg_dir}; reusing it.")
}
dir.create(extdata, showWarnings = FALSE, recursive = TRUE)

tf_file <- file.path(extdata, paste0(motif_set, "_tf_bindsites.rds"))
if (file.exists(tf_file)) {
    log_success("Found existing binding sites: {basename(tf_file)} ",
        "({round(file.size(tf_file) / 1e9, 2)} GB). Motif matching will be skipped.")
} else {
    log_warn("No binding-site file at {tf_file}.")
    log_warn("Motif matching will run from scratch, which takes hours.")
    log_warn("To reuse an existing file, copy it there under that exact name.")
}

## ------------------------------------------------------------------
## 2. CpG positions
## ------------------------------------------------------------------

gc_sites <- NULL
if (restrict_to_cpg) {
    if (file.exists(cpg_cache)) {
        log_info("Loading cached CpG positions ...")
        gc_sites <- readRDS(cpg_cache)
    } else {
        log_info("Locating CpG positions genome-wide (one pass, ~10-20 min) ...")
        t0 <- Sys.time()
        gc_sites <- cpgSites(genome)
        log_success("Found {format(length(gc_sites), big.mark = ',')} CpGs ",
            "in {round(difftime(Sys.time(), t0, units = 'mins'), 1)} min.")
        saveRDS(gc_sites, cpg_cache)
    }
    log_info("GC table will be restricted to {format(length(gc_sites), big.mark = ',')} positions.")
} else {
    log_warn("Building the GC table at every genomic position. ",
        "Expect a very large object.")
}

## ------------------------------------------------------------------
## 3. Genome-wide annotations and the motif GC frequency table
## ------------------------------------------------------------------

log_info("Building annotations for {motif_set} ...")
t0 <- Sys.time()
build_annotations(
    annotations = motif_set,
    pkg.base.dir = pkg_dir,
    chunk_size = chunk_size,
    genome = genome,
    cores = cores,
    enhancer = NULL,
    keep_score = FALSE, # methylTFR never reads it
    gc_sites = gc_sites,
    bin_scope = "genome" # matches the quantiles used for the GC freq tables
)
log_success("Done in {round(difftime(Sys.time(), t0, units = 'mins'), 1)} min.")

## ------------------------------------------------------------------
## 4. Distal-restricted GC frequency table
## ------------------------------------------------------------------
## Only the GC frequency table differs; the binding sites are the same
## object. That is why getTFbindsites() has no "jaspar2020_distal"
## entry and analysis code passes "jaspar2020" for the binding sites
## while passing "jaspar2020_distal" for the GC frequencies.

if (build_distal) {
    log_info("Building the distal-restricted GC frequency table ...")
    distal <- readRDS(distal_rds)
    t0 <- Sys.time()
    build_annotations(
        annotations = motif_set,
        pkg.base.dir = pkg_dir,
        chunk_size = chunk_size,
        genome = genome,
        cores = cores,
        enhancer = distal,
        keep_score = FALSE,
        gc_sites = gc_sites,
        bin_scope = "genome"
    )
    log_success("Done in {round(difftime(Sys.time(), t0, units = 'mins'), 1)} min.")
}

## ------------------------------------------------------------------
## 5. Report
## ------------------------------------------------------------------

files <- list.files(extdata, full.names = TRUE)
sizes <- vapply(files, file.size, numeric(1))
report <- data.frame(
    file = basename(files),
    MB = round(sizes / 1e6, 1),
    row.names = NULL
)
report <- report[order(-report$MB), ]

cat("\n---------------- inst/extdata ----------------\n")
print(report, row.names = FALSE)
cat(sprintf("%-44s %8.1f\n", "TOTAL", sum(report$MB)))
cat("----------------------------------------------\n\n")

log_info("Install with:")
log_info("    R CMD INSTALL {pkg_dir}")
log_info("Then verify with:  Rscript 03_verify.R")
