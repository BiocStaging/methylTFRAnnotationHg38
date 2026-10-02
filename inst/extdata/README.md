# methylTFRAnnotationHg38: `inst/extdata`

This directory contains `metadata.csv`, the AnnotationHub metadata for
the resources served by methylTFRAnnotationHg38. The data files themselves
are hosted on AnnotationHub (they are too large to ship with the
package) and are downloaded on first use by the accessor functions.

## Resources

| Title | R class | Accessor | Help page |
|---|---|---|---|
| `altius_motif_gcfreq.rds` | `list` | `getGCfreq("altius")` | `?motif_gcfreq` |
| `altius_tf_bindsites.rds` | `GRangesList` | `getTFbindsites("altius")` | `?tf_bindsites` |
| `cisbpv2_motif_gcfreq.rds` | `list` | `getGCfreq("cisbpv2")` | `?motif_gcfreq` |
| `cisbpv2_tf_bindsites.rds` | `GRangesList` | `getTFbindsites("cisbpv2")` | `?tf_bindsites` |
| `genomewide_GC_hg38.rds` | `GRanges` | `getGenomeGC()` | `?genomewide_GC` |
| `jaspar2020_distal_motif_gcfreq.rds` | `list` | `getGCfreq("jaspar2020_distal")` | `?motif_gcfreq` |
| `jaspar2020_motif_gcfreq.rds` | `list` | `getGCfreq("jaspar2020")` | `?motif_gcfreq` |
| `jaspar2020_tf_bindsites.rds` | `GRangesList` | `getTFbindsites("jaspar2020")` | `?tf_bindsites` |

* `<set>_tf_bindsites.rds`: genome-wide motif matches for one motif set,
  one `GRanges` per motif, each range extended by 200 bp on either side
  of the match.
* `<set>_motif_gcfreq.rds`: one 5 x n numeric matrix per motif. Rows are
  genome-wide GC quintiles (lowest first), columns are positions of a
  30 nt window along the binding site; each column gives the fraction of
  the motif's sites in each quintile and sums to one.
* `genomewide_GC_hg38.rds`: non-overlapping 30 nt windows across the
  primary chromosomes with their GC fraction (`GC_bias`) and quintile
  (`GC_bin`); the quintile boundaries are stored in the object's metadata.
* `<set>_distal_motif_gcfreq.rds`: as `<set>_motif_gcfreq.rds`, but
  computed only from binding sites overlapping distal regulatory regions
  (Ensembl Regulatory Build v104). It shares the binding sites of the
  base set.

## Columns of `metadata.csv`

| Column | Meaning |
|---|---|
| `Title` | Resource name; also the file name on AnnotationHub |
| `Description` | What the resource contains |
| `BiocVersion` | Bioconductor version the resource was added in |
| `Genome` | Genome assembly (`hg38`) |
| `SourceType` | Format of the source file (`RDS`) |
| `SourceUrl` | Where the input data (motifs, genome) come from |
| `SourceVersion` | Version of the motif collection or genome |
| `Species`, `TaxonomyId` | Organism and NCBI taxonomy ID |
| `Coordinate_1_based` | Coordinates are 1-based (`TRUE`) |
| `DataProvider` | Provider of the motif collection or genome |
| `Maintainer` | Maintainer of the resource |
| `RDataClass` | R class of the object (`GRangesList`, `list`, `GRanges`) |
| `DispatchClass` | How AnnotationHub loads the file (`Rds`) |
| `Location_Prefix`, `RDataPath` | Where the file is stored; joined, they give the download URL |
| `Tags` | Search tags for `AnnotationHub::query()` |

## How the data were made

* `inst/scripts/make-data.R` builds every resource with
  [methylTFRAnnotationBuilder](https://github.com/EpigenomeInformatics/methylTFRAnnotationBuilder).
* `inst/scripts/make-metadata.R` writes this `metadata.csv`.
