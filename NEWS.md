# methylTFRAnnotationHg38 0.99.11

Changes in response to the second review:

* Depends on R (>= 4.6.0).
* Input checking: an invalid `methylTFRAnnotationHg38.datadir` option or
  `METHYL_TFRANNOTATION_HG38_DIR` value (not a single string, or a directory
  that does not exist) now gives an informative error. Tests
  added.
* Added `inst/extdata/README.md` describing `metadata.csv` and
  each AnnotationHub resource.
* The vignette's Installation section also installs methylTFR.

# methylTFRAnnotationHg38 0.99.10

* BiocCheck fixes: the data help pages (`?motif_gcfreq`,
  `?tf_bindsites`, `?genomewide_GC`) have a `\value` section, and
  their examples run (reading `metadata.csv`) instead of using
  `\dontrun`.
* Shortened vignette lines to at most 80 characters.
* Added `CITATION.cff`.

# methylTFRAnnotationHg38 0.99.9

Changes in response to the Bioconductor review:

* `Authors@R` lists the funders (ERA-NET Transcan-Neu III - EPILUNAR,
  grant 01KT2409; Saarland University NanoBioMed Young Investigator
  Grant) with the `fnd` role.

* New package help page `?methylTFRAnnotationHg38` summarising the resources and
  linking the accessors.

* New help pages describing the structure of each AnnotationHub resource
  listed in `inst/extdata/metadata.csv`: `?motif_gcfreq`
  (`?altius_motif_gcfreq`, ...), `?tf_bindsites`
  (`?altius_tf_bindsites`, ...) and `?genomewide_GC`
  (`?genomewide_GC_hg38`).

* New exported `availableMotifSets()` listing the accepted motif sets.

* `getGenomeGC()` accepts `"hg38"` (or no argument) and gives an
  informative error for any other assembly.

* `getTFbindsites()` accepts `_distal` set names and returns the
  binding sites of the base set.

* Unit tests now cover the local-directory and AnnotationHub code paths
  (the hub is mocked, so the tests run offline).

* The vignette no longer hides code: the synthetic example data are
  built in a visible chunk, every chunk is evaluated, and each accessor's
  output is printed and explained.

# methylTFRAnnotationHg38 0.99.8

* Initial submission to Bioconductor.

* Provides transcription factor binding sites, motif GC frequency
  tables and a genome-wide GC distribution for hg38, for the
  jaspar2020, cisbpv2 and altius motif sets, plus a distal-restricted
  GC frequency table for jaspar2020.

* Resources are hosted on AnnotationHub and downloaded on first use.

* The genome-wide GC table uses non-overlapping 30 nt windows binned
  into five genome-wide GC quintiles. Bin boundaries are recorded with
  the object so that the motif frequency tables and the genome table
  are always binned on the same scale.

* Fixes on NOTES from biocchecks

* Added 4 spaces as suggested
