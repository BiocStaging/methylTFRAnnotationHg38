# The AnnotationHub branch of the resolver, tested without network
# access by replacing the two small functions that talk to the hub.

test_that("resources are looked up on AnnotationHub by title", {
    old <- options(methylTFRAnnotationHg38.datadir = NULL)
    on.exit(options(old), add = TRUE)
    old_env <- Sys.getenv("METHYL_TFRANNOTATION_HG38_DIR", unset = NA)
    Sys.unsetenv("METHYL_TFRANNOTATION_HG38_DIR")
    on.exit(
        if (!is.na(old_env)) Sys.setenv(METHYL_TFRANNOTATION_HG38_DIR = old_env),
        add = TRUE
    )

    fetched <- character(0)
    local_mocked_bindings(
        .hub_titles = function() c(
            AH1 = "altius_motif_gcfreq.rds",
            AH2 = "altius_tf_bindsites.rds",
            AH3 = "genomewide_GC_hg38.rds"
        ),
        .hub_get = function(id) {
            fetched <<- c(fetched, id)
            paste("object", id)
        }
    )
    expect_identical(getGCfreq("altius"), "object AH1")
    expect_identical(getTFbindsites("altius"), "object AH2")
    expect_identical(getGenomeGC(), "object AH3")
    expect_identical(fetched, c("AH1", "AH2", "AH3"))
    expect_error(getGCfreq("cisbpv2"), "Resource not found on AnnotationHub")
})

test_that("the real AnnotationHub records are reachable", {
    skip_on_cran()
    skip_on_ci()
    skip_if_offline()
    skip_if_not(nzchar(Sys.getenv("METHYLTFR_TEST_HUB")),
        "set METHYLTFR_TEST_HUB=1 to query AnnotationHub"
    )
    titles <- methylTFRAnnotationHg38:::.hub_titles()
    md <- utils::read.csv(system.file("extdata", "metadata.csv",
        package = "methylTFRAnnotationHg38"
    ))
    expect_setequal(unname(titles), md$Title)
})
