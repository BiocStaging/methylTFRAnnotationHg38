test_that("motif set validation rejects unknown sets", {
    expect_error(getGCfreq("not_a_motif_set"), "Invalid motif set")
    expect_error(getGCfreq(c("a", "b")), "Invalid motif set")
})

test_that("motif set validation is case insensitive", {
    expect_identical(methylTFRAnnotationHg38:::.check_motif_set("ALTIUS"), "altius")
})

test_that("metadata.csv covers every declared motif set", {
    md <- utils::read.csv(system.file("extdata", "metadata.csv", package = "methylTFRAnnotationHg38"))
    for (s in methylTFRAnnotationHg38:::.MOTIF_SETS) {
        expect_true(any(grepl(paste0("^", s, "_"), md$Title)),
            info = s
        )
    }
})
