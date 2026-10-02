# The accessors read .rds files from a local directory when the
# methylTFRAnnotationHg38.datadir option (or the
# METHYL_TFRANNOTATION_HG38_DIR environment variable)
# is set. These tests write small objects with the same structure as the
# AnnotationHub resources and check they come back unchanged.

make_local_resources <- function(dir) {
    tfbs <- GenomicRanges::GRangesList(
        MA0001.1 = GenomicRanges::GRanges(
            c("chr1:1001-1411", "chr2:5001-5411"),
            strand = c("+", "-")
        )
    )
    freq <- matrix(c(0.1, 0.2, 0.4, 0.2, 0.1), nrow = 5, ncol = 3)
    gc <- GenomicRanges::GRanges(c("chr1:1-30", "chr1:31-60", "chr1:61-90"))
    gc$GC_bias <- c(0.30, 0.43, 0.63)
    gc$GC_bin <- c(1L, 3L, 5L)
    for (s in availableMotifSets()) {
        saveRDS(
            list(MA0001.1 = freq),
            file.path(dir, paste0(s, "_motif_gcfreq.rds"))
        )
        base <- sub("_distal$", "", s)
        saveRDS(tfbs, file.path(dir, paste0(base, "_tf_bindsites.rds")))
    }
    saveRDS(gc, file.path(dir, "genomewide_GC_hg38.rds"))
    list(tfbs = tfbs, freq = freq, gc = gc)
}

with_local_dir <- function(code) {
    dir <- tempfile("annot_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    old <- options(methylTFRAnnotationHg38.datadir = dir)
    on.exit(options(old), add = TRUE)
    ref <- make_local_resources(dir)
    code(dir, ref)
}

test_that("accessors read every resource from a local directory", {
    with_local_dir(function(dir, ref) {
        for (s in availableMotifSets()) {
            freq <- getGCfreq(s)
            expect_type(freq, "list")
            expect_identical(freq$MA0001.1, ref$freq)
            expect_equal(unname(colSums(freq$MA0001.1)), rep(1, 3))
        }
        tfbs <- getTFbindsites(toupper(availableMotifSets()[1]))
        expect_s4_class(tfbs, "GRangesList")
        expect_identical(tfbs, ref$tfbs)

        gc <- getGenomeGC()
        expect_s4_class(gc, "GRanges")
        cols <- names(GenomicRanges::mcols(gc))
        expect_true(all(c("GC_bias", "GC_bin") %in% cols))
        expect_identical(gc, ref$gc)
        expect_identical(getGenomeGC(toupper("hg38")), gc)
    })
})

test_that("a missing local file gives an informative error", {
    with_local_dir(function(dir, ref) {
        file.remove(file.path(dir, "genomewide_GC_hg38.rds"))
        expect_error(getGenomeGC(), "Annotation file not found")
    })
})

test_that("the environment variable is used when the option is unset", {
    dir <- tempfile("annot_env_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    ref <- make_local_resources(dir)
    old_opt <- options(methylTFRAnnotationHg38.datadir = NULL)
    on.exit(options(old_opt), add = TRUE)
    old_env <- Sys.getenv("METHYL_TFRANNOTATION_HG38_DIR", unset = NA)
    Sys.setenv(METHYL_TFRANNOTATION_HG38_DIR = dir)
    on.exit(
        if (is.na(old_env)) Sys.unsetenv("METHYL_TFRANNOTATION_HG38_DIR")
        else Sys.setenv(METHYL_TFRANNOTATION_HG38_DIR = old_env),
        add = TRUE
    )
    expect_identical(methylTFRAnnotationHg38:::.local_dir(), dir)
    expect_identical(getGenomeGC(), ref$gc)
})

test_that("no local directory is used when neither option nor env var is set", {
    old_opt <- options(methylTFRAnnotationHg38.datadir = NULL)
    on.exit(options(old_opt), add = TRUE)
    old_env <- Sys.getenv("METHYL_TFRANNOTATION_HG38_DIR", unset = NA)
    Sys.unsetenv("METHYL_TFRANNOTATION_HG38_DIR")
    on.exit(
        if (!is.na(old_env)) Sys.setenv(METHYL_TFRANNOTATION_HG38_DIR = old_env),
        add = TRUE
    )
    expect_null(methylTFRAnnotationHg38:::.local_dir())
})

test_that("an invalid local directory setting gives an informative error", {
    old <- options(methylTFRAnnotationHg38.datadir = c("a", "b"))
    on.exit(options(old), add = TRUE)
    expect_error(getGenomeGC(), "single character string")

    options(methylTFRAnnotationHg38.datadir = 1)
    expect_error(getGCfreq("altius"), "single character string")

    missing_dir <- file.path(tempdir(), "does_not_exist")
    options(methylTFRAnnotationHg38.datadir = missing_dir)
    expect_error(getTFbindsites("altius"), "does not exist")
})
