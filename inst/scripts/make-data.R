# How the methylTFRAnnotationHg38 resources were produced
#
# Every resource is built by methylTFRAnnotationBuilder
# (https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder)
# from a BSgenome sequence and a collection of position weight
# matrices. This script records that build and regenerates the files.
#
# Four kinds of resource:
#
#   <set>_tf_bindsites.rds          GRangesList, one GRanges per motif
#   <set>_motif_gcfreq.rds          list of 5 x n matrices, one per motif
#   genomewide_GC_hg38.rds          GRanges of 30 nt windows with GC bins
#   <set>_distal_motif_gcfreq.rds   as above, distal regions only
#
# The genome GC scan takes about a minute. Motif matching and the GC
# frequency tables take several hours per motif set on 24 and 8 cores
# respectively.

library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(GenomicRanges)
library(BiocParallel)
library(methylTFRAnnotationBuilder)

genome <- BSgenome.Hsapiens.UCSC.hg38
pkg_dir <- "methylTFRAnnotationHg38"
extdata <- file.path(pkg_dir, "inst", "extdata")

# chr1-22, X and Y. chrM is excluded: its GC content and its
# methylation are both atypical, and 16 kb contributes nothing to
# genome-wide quantiles.
chromosomes <- standardChrs(genome)


# ------------------------------------------------------------------
# 1. Motif collections
# ------------------------------------------------------------------

# --- jaspar2020 ---------------------------------------------------
# JASPAR2020 CORE filtered to Homo sapiens (NCBI taxonomy 9606): 629
# matrices, converted from counts to PWMs. Keyed on the JASPAR TF
# name.

jaspar2020 <- TFBSTools::toPWM(
    TFBSTools::getMatrixSet(
        JASPAR2020::JASPAR2020,
        list(species = 9606, collection = "CORE")
    )
)

# --- cisbpv2 ------------------------------------------------------
# CIS-BP version 2, from the chromVARmotifs package
# (https://github.com/GreenleafLab/chromVARmotifs), distributed on
# GitHub only. human_pwms_v2 here; mm10 uses mouse_pwms_v2, 884
# matrices. chromVARmotifs' own names are kept unchanged.

utils::data("human_pwms_v2", package = "chromVARmotifs")
cisbpv2 <- chromVARmotifs::human_pwms_v2

# --- altius -------------------------------------------------------
# Vierstra motif archetypes v1.0
# (https://www.vierstra.org/resources/motif_clustering): 286 clusters
# derived from 2,179 motif models drawn from JASPAR, HOCOMOCO, JOLMA,
# TRANSFAC and others. Resources are keyed on the cluster identifier.
##
# These are supplied as genomic positions, not as PWMs, so they do not
# go through findTFBindSites(). The hg38 sites were obtained with
# ChrAccR::getMotifClusterAnnot_altius(), which downloads a
# precomputed tfMotifClusters_hg38.rds. 
##
# Both routes apply the same reduction. Occurrences from every cluster
# are pooled, overlaps are resolved in favour of the higher MOODS
# score, and the survivors are split back out by cluster. That
# competition is global, so one archetype can take a locus from
# another. Resizing to the footprint width happens after the
# reduction; doing it first would let widened neighbours exclude each
# other and change which matches survive.

altius <- readRDS("altius_tf_bindsites.rds") # see above for provenance


# ------------------------------------------------------------------
# 2. Binding sites
# ------------------------------------------------------------------

# matchMotifs(out = "positions") runs once per chromosome, matching
# every motif of a set in the same pass.
##
# Each match is widened to motif width + 400. That width matters:
# methylTFR::computeDeviation() widens the stored range by a further
# 130 bases and reads methylation across offsets -250 to +250 from the
# midpoint, so a narrower stored range truncates every footprint
# without any error being raised.
##
# The motif match score is dropped. methylTFR never reads it, and it
# costs eight bytes on each of several hundred million sites.

for (set in c("jaspar2020", "cisbpv2")) {
    out <- file.path(extdata, paste0(set, "_tf_bindsites.rds"))
    if (file.exists(out)) next
    saveRDS(
        findTFBindSites(
            genome = genome,
            motifs = get(set),
            BPPARAM = MulticoreParam(workers = 24),
            keep_score = FALSE,
            flank = 400,
            chromosomes = chromosomes
        ),
        out
    )
}


# ------------------------------------------------------------------
# 3. Genome-wide GC distribution
# ------------------------------------------------------------------

# GC content in non-overlapping 30 nt windows: 102,942,317 windows for
# hg38, one per 30 bases of the primary chromosomes.
##
# GC is the count of C and G divided by the window width. Dividing by
# the number of called bases instead makes an all-N window 0/0, and
# any step that then drops those windows leaves every later value
# misaligned from the coordinates it is stored against. Dividing by
# the width scores a gap window 0 and keeps it in place.
##
# Windows are assigned to five bins delimited by quantiles of the
# genome-wide GC distribution. Those boundaries are stored on the
# object as metadata(gc)$gc_breaks, and build_annotations() reads them
# from there when binning the motif frequency tables in section 4, so
# the observed and expected sides of a deviation score are guaranteed
# to use the same scale.

gc_file <- file.path(extdata, "genomewide_GC_hg38.rds")
if (!file.exists(gc_file)) {
    saveRDS(
        computeGCgenome(
            genome = genome,
            cores = 12,
            step = 30L,
            bin_scope = "genome",
            chromosomes = chromosomes
        ),
        gc_file
    )
}


# ------------------------------------------------------------------
# 4. Motif GC frequency tables
# ------------------------------------------------------------------

# For each motif, and each offset along its footprint window, the
# proportion of that motif's binding sites whose local 30 nt GC
# content falls in each of the five bins. Columns sum to one.
##
# This is the expected side of the deviation score. Multiplying the
# table by a sample's per-GC-bin mean methylation predicts the
# methylation profile from sequence composition alone, which is then
# subtracted from the observed profile.
##
# build_annotations() reads the binding sites and the genome GC table
# from inst/extdata and skips any resource already present, so it is
# safe to re-run after an interruption.

for (set in c("jaspar2020", "cisbpv2", "altius")) {
    build_annotations(
        annotations = set,
        pkg.base.dir = pkg_dir,
        chunk_size = 15,
        genome = genome,
        cores = 8,
        enhancer = NULL,
        keep_score = FALSE,
        step = 30L,
        bin_scope = "genome",
        chromosomes = chromosomes
    )
}


# ------------------------------------------------------------------
# 5. Distal-restricted frequency table
# ------------------------------------------------------------------

# jaspar2020_distal is section 4 restricted to binding sites
# overlapping distal regulatory regions. Only the frequency table
# differs; the binding sites are the unrestricted ones. That is why
# getTFbindsites() has no jaspar2020_distal entry, and why analysis
# code passes "jaspar2020" for the sites and "jaspar2020_distal" for
# the frequencies.
##
# distal_regions.RDS: Taken from Ensembl Regulatory Build v104, filtered 
# to hg38 and to distal Which annotation the transcription start sites came 
# from, what distance threshold separates distal from proximal, and which build.
# Without it this resource cannot be reproduced.

distal <- readRDS("distal_regions.RDS")

build_annotations(
    annotations = "jaspar2020",
    pkg.base.dir = pkg_dir,
    chunk_size = 15,
    genome = genome,
    cores = 8,
    enhancer = distal,
    keep_score = FALSE,
    step = 30L,
    bin_scope = "genome",
    chromosomes = chromosomes
)

sessionInfo()