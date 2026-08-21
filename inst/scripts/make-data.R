# #################################################################################
# How the methylTFRAnnotationHg38 resources were produced
#
# Every resource is built by methylTFRAnnotationBuilder
# (https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder)
# from a BSgenome sequence and a collection of position weight
# matrices. This script records that build and regenerates the files.
#
# Four kinds of resource produced in inst/extdata:
#   1. <set>_tf_bindsites.rds          GRangesList, one GRanges per motif
#   2. <set>_motif_gcfreq.rds          list of 5 x n matrices, one per motif
#   3. genomewide_GC_hg38.rds          GRanges of 30 nt windows with GC bins
#   4. <set>_distal_motif_gcfreq.rds   as above, distal regions only
#
# Motif sets covered:
#   - jaspar2020: JASPAR2020 CORE (Homo sapiens)
#   - cisbpv2:    CIS-BP v2 (human_pwms_v2 from chromVARmotifs)
#   - altius:     Vierstra motif archetypes v1.0 (from hg38 archetype BED scan)
# #################################################################################

library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(GenomicRanges)
library(BiocParallel)
library(data.table)
library(methylTFRAnnotationBuilder)

genome <- BSgenome.Hsapiens.UCSC.hg38
pkg_dir <- "methylTFRAnnotationHg38"
extdata <- file.path(pkg_dir, "inst", "extdata")
dir.create(extdata, recursive = TRUE, showWarnings = FALSE)

# chr1-22, X and Y. chrM is excluded.
chromosomes <- standardChrs(genome)


# ##############################################################################
# 1. Altius Archetypes v1.0 Binding Sites
# ##############################################################################
# Source: https://www.vierstra.org/resources/motif_clustering
# File: hg38.archetype_motifs.v1.0.bed.gz (286 archetype clusters)
#
# Provenance & Processing:
#   1. Streamed per chromosome via tabix in 20 Mb tiles with 10 kb overlap margin.
#   2. De-duplicated globally across all pooled cluster matches by highest MOODS 
#      score.
#   3. Grouped into a GRangesList by archetype cluster ID (cid).
#   4. Resized symmetrically to (match_width + 400) centered on the motif.

altius_out <- file.path(extdata, "altius_tf_bindsites.rds")

if (!file.exists(altius_out)) {
    release_url <- "https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0"
    bed_url <- paste0(release_url, "/hg38.archetype_motifs.v1.0.bed.gz")
    bed_file <- "hg38.archetype_motifs.v1.0.bed.gz"
    tbi_file <- paste0(bed_file, ".tbi")
    
    if (!file.exists(bed_file)) download.file(bed_url, bed_file, method = "wget", extra = "-c")
    if (!file.exists(tbi_file)) download.file(paste0(bed_url, ".tbi"), tbi_file, method = "wget", extra = "-c")
    
    nonOverlappingByScore <- function(gr) {
        if (length(gr) < 2) return(gr)
        rk <- integer(length(gr))
        rk[order(-gr$score, start(gr), seq_along(gr))] <- seq_along(gr)
        kept <- vector("list", 0)
        idx <- seq_along(gr)
        while (length(idx)) {
            g <- gr[idx]
            r <- rk[idx]
            hits <- findOverlaps(g, ignore.strand = TRUE, drop.self = TRUE)
            best_nb <- rep(.Machine$integer.max, length(g))
            if (length(hits)) {
                agg <- data.table(q = queryHits(hits), r = r[subjectHits(hits)])[, .(m = min(r)), by = q]
                best_nb[agg$q] <- agg$m
            }
            win <- which(r < best_nb)
            if (length(win) == 0L) break
            kept[[length(kept) + 1L]] <- g[win]
            loser <- unique(subjectHits(findOverlaps(g[win], g, ignore.strand = TRUE)))
            idx <- idx[-loser]
        }
        if (length(kept) == 0L) return(gr[0])
        do.call(c, kept)
    }
    
    tile_width <- 20e6
    tile_margin <- 10e3
    chr_len <- GenomeInfoDb::seqlengths(genome)
    
    all_tiles <- list()
    for (chr in chromosomes) {
        len <- as.integer(chr_len[chr])
        starts <- seq(1, len, by = tile_width)
        
        pieces <- lapply(starts, function(core_from) {
            core_to <- min(core_from + tile_width - 1L, len)
            region <- sprintf("%s:%d-%d", chr, max(1L, core_from - tile_margin), min(len, core_to + tile_margin))
            txt <- system2("tabix", c(shQuote(bed_file), region), stdout = TRUE)
            if (length(txt) == 0) return(NULL)
            
            dt <- data.table::fread(text = paste(txt, collapse = "\n"), header = FALSE, sep = "\t",
                                    select = 1:6, col.names = c("chr", "start", "end", "cid", "score", "strand"))
            if (nrow(dt) == 0) return(NULL)
            
            gr <- GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$start + 1L, end = dt$end),
                          strand = dt$strand, cid = dt$cid, score = dt$score)
            gr <- nonOverlappingByScore(gr)
            gr[start(gr) >= core_from & start(gr) <= core_to]
        })
        pieces <- pieces[!vapply(pieces, is.null, logical(1))]
        if (length(pieces)) all_tiles[[chr]] <- do.call(c, pieces)
    }
    
    all_gr <- do.call(c, all_tiles)
    altius_grl <- split(all_gr, all_gr$cid)
    altius_grl <- endoapply(altius_grl, function(g) {
        g <- resize(g, width = width(g)[1] + 400L, fix = "center")
        mcols(g) <- NULL
        g
    })
    
    saveRDS(altius_grl, altius_out)
}


# ##############################################################################
# 2. PWM-based Binding Sites (JASPAR2020 & CIS-BP v2)
# ##############################################################################

# --- JASPAR2020 ---
# CORE filtered to Homo sapiens (NCBI taxonomy 9606): 629 matrices.
jaspar2020 <- TFBSTools::toPWM(
    TFBSTools::getMatrixSet(
        JASPAR2020::JASPAR2020,
        list(species = 9606, collection = "CORE")
    )
)

# --- CIS-BP v2 ---
# Human PWMs from chromVARmotifs.
utils::data("human_pwms_v2", package = "chromVARmotifs")
cisbpv2 <- chromVARmotifs::human_pwms_v2

# Scan genome chromosomes & widen matches by +400 bp
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


# ##############################################################################
# 3. Genome-wide GC Distribution
# ##############################################################################
# Non-overlapping 30 nt windows across primary hg38 chromosomes, binned into
# genome-wide GC quintiles.

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


# ##############################################################################
# 4. Motif GC Frequency Tables
# ##############################################################################

for (set in c("jaspar2020", "cisbpv2", "altius")) {
    out <- file.path(extdata, paste0(set, "_motif_gcfreq.rds"))
    if (file.exists(out)) next
    
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
        chromosomes = chromosomes,
        jaspar_opts = list(species = 9606, collection = "CORE")
    )
}


# ##############################################################################
# 5. Distal-restricted Frequency Table
# ##############################################################################
# Section 4 restricted to binding sites overlapping distal regulatory regions. 
# Only the frequency table differs; the binding sites are the unrestricted ones.
#
# distal_regions.RDS provenance: 
# Taken from Ensembl Regulatory Build v104, filtered to hg38 and to distal.

distal <- readRDS("distal_regions.RDS")
distal_out <- file.path(extdata, "jaspar2020_distal_motif_gcfreq.rds")

if (!file.exists(distal_out)) {
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
        chromosomes = chromosomes,
        jaspar_opts = list(species = 9606, collection = "CORE")
    )
    # The builder will output jaspar2020_distal_motif_gcfreq.rds due to the enhancer param
}

sessionInfo()