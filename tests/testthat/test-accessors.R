test_that("motif set validation rejects unknown sets", {
    expect_error(getGCfreq("not_a_motif_set"), "Invalid motif set")
    expect_error(getGCfreq(c("a", "b")), "Invalid motif set")
    expect_error(getTFbindsites(NA_character_), "Invalid motif set")
    expect_error(getTFbindsites(1), "Invalid motif set")
    expect_error(getTFbindsites(character(0)), "Invalid motif set")
})

test_that("motif set validation is case insensitive", {
    expect_identical(
        methylTFRAnnotationHg38:::.check_motif_set("ALTIUS"), "altius"
    )
})

test_that("availableMotifSets() returns the declared sets", {
    sets <- availableMotifSets()
    expect_type(sets, "character")
    expect_identical(
        sets, c("altius", "cisbpv2", "jaspar2020", "jaspar2020_distal")
    )
    for (s in sets) {
        expect_identical(methylTFRAnnotationHg38:::.check_motif_set(s), s)
    }
})

test_that("metadata.csv covers every declared motif set and the genome GC", {
    md <- utils::read.csv(system.file("extdata", "metadata.csv",
        package = "methylTFRAnnotationHg38"
    ))
    for (s in availableMotifSets()) {
        expect_true(paste0(s, "_motif_gcfreq.rds") %in% md$Title, info = s)
        base <- sub("_distal$", "", s)
        expect_true(paste0(base, "_tf_bindsites.rds") %in% md$Title, info = s)
    }
    expect_true("genomewide_GC_hg38.rds" %in% md$Title)
    expect_true(all(md$Genome == "hg38"))
    expect_true(all(md$Species == "Homo sapiens"))
})

test_that("getGenomeGC() only accepts this package's assembly", {
    expect_error(getGenomeGC("mm10"), "only")
    expect_error(getGenomeGC(c("hg38", "hg38")), "only")
    expect_error(getGenomeGC(NA_character_), "only")
})
