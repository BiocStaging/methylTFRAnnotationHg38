# 1. Define the shared (constant) variables for cleaner code
bioc_version    <- "3.23"
genome          <- "hg38"
source_type     <- "RDS"
species         <- "Homo sapiens"
tax_id          <- 9606
coord_1_based   <- TRUE
maintainer      <- "Irem B. Gunduz <irembgunduz@gmail.com>"
dispatch_class  <- "Rds"
location_prefix <- "https://bioconductorhubs.blob.core.windows.net/annotationhub/"

# 2. Build the data frame
metadata <- data.frame(
  Title = c(
    "altius_motif_gcfreq.rds",
    "altius_tf_bindsites.rds",
    "cisbpv2_motif_gcfreq.rds",
    "cisbpv2_tf_bindsites.rds",
    "genomewide_GC_hg38.rds",
    "jaspar2020_distal_motif_gcfreq.rds",
    "jaspar2020_motif_gcfreq.rds",
    "jaspar2020_tf_bindsites.rds"
  ),
  
  Description = c(
    "GC bin frequency tables for ALTIUS motifs on hg38. One five-row matrix per motif giving the number of binding sites falling in each genome-wide GC quintile, used by methylTFR to compute the expected methylation a motif would show from GC content alone.",
    "Genome-wide ALTIUS transcription factor binding site predictions for hg38, one GRanges per motif, each range extended by 200 bases on either side of the motif match so that methylTFR can read methylation across the footprint window.",
    "GC bin frequency tables for CISBPV2 motifs on hg38. One five-row matrix per motif giving the number of binding sites falling in each genome-wide GC quintile, used by methylTFR to compute the expected methylation a motif would show from GC content alone.",
    "Genome-wide CISBPV2 transcription factor binding site predictions for hg38, one GRanges per motif, each range extended by 200 bases on either side of the motif match so that methylTFR can read methylation across the footprint window.",
    "Genome-wide GC content distribution for hg38. A GRanges of tiled windows carrying GC_bias and a GC_bin assignment into genome-wide quintiles. methylTFR uses it to assign each methylation call to a GC bin.",
    "GC bin frequency tables for JASPAR2020_DISTAL motifs on hg38. One five-row matrix per motif giving the number of binding sites falling in each genome-wide GC quintile, used by methylTFR to compute the expected methylation a motif would show from GC content alone.",
    "GC bin frequency tables for JASPAR2020 motifs on hg38. One five-row matrix per motif giving the number of binding sites falling in each genome-wide GC quintile, used by methylTFR to compute the expected methylation a motif would show from GC content alone.",
    "Genome-wide JASPAR2020 transcription factor binding site predictions for hg38, one GRanges per motif, each range extended by 200 bases on either side of the motif match so that methylTFR can read methylation across the footprint window."
  ),
  
  BiocVersion = rep(bioc_version, 8),
  Genome = rep(genome, 8),
  SourceType = rep(source_type, 8),
  
  SourceUrl = c(
    "https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/",
    "https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/",
    "https://github.com/GreenleafLab/chromVARmotifs",
    "https://github.com/GreenleafLab/chromVARmotifs",
    "https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg38/",
    "https://jaspar.elixir.no/",
    "https://jaspar.elixir.no/",
    "https://jaspar.elixir.no/"
  ),
  
  SourceVersion = c(
    "Vierstra motif archetypes v1.0",
    "Vierstra motif archetypes v1.0",
    "CIS-BP v2 (chromVARmotifs pwms_v2)",
    "CIS-BP v2 (chromVARmotifs pwms_v2)",
    "hg38",
    "JASPAR2020 CORE",
    "JASPAR2020 CORE",
    "JASPAR2020 CORE"
  ),
  
  Species = rep(species, 8),
  TaxonomyId = rep(tax_id, 8),
  Coordinate_1_based = rep(coord_1_based, 8),
  
  DataProvider = c(
    "Altius Institute",
    "Altius Institute",
    "CIS-BP",
    "CIS-BP",
    "UCSC",
    "JASPAR",
    "JASPAR",
    "JASPAR"
  ),
  
  Maintainer = rep(maintainer, 8),
  
  RDataClass = c(
    "list",
    "GRangesList",
    "list",
    "GRangesList",
    "GRanges",
    "list",
    "list",
    "GRangesList"
  ),
  
  DispatchClass = rep(dispatch_class, 8),
  Location_Prefix = rep(location_prefix, 8),
  
  RDataPath = c(
    "methylTFRAnnotationHg38/altius_motif_gcfreq.rds",
    "methylTFRAnnotationHg38/altius_tf_bindsites.rds",
    "methylTFRAnnotationHg38/cisbpv2_motif_gcfreq.rds",
    "methylTFRAnnotationHg38/cisbpv2_tf_bindsites.rds",
    "methylTFRAnnotationHg38/genomewide_GC_hg38.rds",
    "methylTFRAnnotationHg38/jaspar2020_distal_motif_gcfreq.rds",
    "methylTFRAnnotationHg38/jaspar2020_motif_gcfreq.rds",
    "methylTFRAnnotationHg38/jaspar2020_tf_bindsites.rds"
  ),
  
  Tags = c(
    "methylTFRAnnotationHg38:GCcontent:MotifAnnotation:ALTIUS",
    "methylTFRAnnotationHg38:TFBS:MotifAnnotation:ALTIUS",
    "methylTFRAnnotationHg38:GCcontent:MotifAnnotation:CISBPV2",
    "methylTFRAnnotationHg38:TFBS:MotifAnnotation:CISBPV2",
    "methylTFRAnnotationHg38:GCcontent:Genome",
    "methylTFRAnnotationHg38:GCcontent:MotifAnnotation:JASPAR2020_DISTAL",
    "methylTFRAnnotationHg38:GCcontent:MotifAnnotation:JASPAR2020",
    "methylTFRAnnotationHg38:TFBS:MotifAnnotation:JASPAR2020"
  ),
  
  stringsAsFactors = FALSE
)

# 3. Write out the metadata.csv file
# Ensure the extdata directory exists before attempting to write
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

write.csv(
  x = metadata, 
  file = "inst/extdata/metadata.csv", 
  row.names = FALSE, 
  quote = TRUE  # This ensures strings containing spaces/commas are properly enclosed
)

# Validate using AnnotationHubData
if (!requireNamespace("AnnotationHubData", quietly = TRUE)) {
    message("Please install AnnotationHubData to validate your metadata.csv")
} else {
    AnnotationHubData::makeAnnotationHubMetadata("inst/extdata")
    message("Metadata validation completed successfully!")
}