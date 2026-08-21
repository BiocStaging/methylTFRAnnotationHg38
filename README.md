# methylTFRAnnotationHg38
[![Test R-universe](https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38/actions/workflows/r-universe-test.yml/badge.svg)](https://github.com/EpigenomeInformatics/methylTFRAnnotationHg38/actions/workflows/r-universe.yml)

`methylTFRAnnotationHg38` provides the genome annotations `methylTFR` needs to compute bias-corrected transcription factor deviation scores on hg38. 

The package supplies three main resources:
*   **Binding sites**: One `GRanges` object per motif, with each range extended so that methylation can be read across the footprint window.
*   **Motif GC frequency tables**: Tables recording how each motif's binding sites distribute across genome-wide GC quintiles.
*   **Genome-wide GC distribution**: Data that assigns each methylation call to a GC bin.

Because the data are too large to ship inside the package, they are hosted on `AnnotationHub`. The data are downloaded dynamically on first use, and subsequent calls are served seamlessly from the local `AnnotationHub` cache.

## Installation

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("AnnotationHub","GenomicRanges", "methylTFRAnnotationHg38"))
```

## Usage

You do not need to interact with `AnnotationHub` directly; the package's accessors handle resolution automatically. The retrieved objects are passed directly to `methylTFR::run_methyltfr()`.

```R
library(methylTFRAnnotationHg38)

tf_bindsites <- getTFbindsites("altius")
gcfreqs      <- getGCfreq("altius")
gc_dist      <- getGenomeGC()
```

### Available Motif Sets

The currently supported motif sets are: `altius`, `cisbpv2`, `jaspar2020`, and `jaspar2020_distal`.

Sets ending in `_distal` provide GC frequency tables computed over distal regulatory regions only. They share unrestricted binding sites with their base set. When using a distal set, pass the base set name to `getTFbindsites()` and the `_distal` name to `getGCfreq()`.

## Local Directory Usage

By default, the accessors will query `AnnotationHub` for the required files. However, if you have downloaded the `.rds` files manually or need to run tests in an environment without internet access, you can bypass the hub. 

To read resources from a local directory instead of `AnnotationHub`, configure the local path using either an R option or an environment variable:
*   **Option**: `options(methylTFRAnnotationHg38.datadir = "/path/to/data")`
*   **Environment Variable**: `Sys.setenv(METHYL_TFRANNOTATION_HG38_DIR = "/path/to/data")`