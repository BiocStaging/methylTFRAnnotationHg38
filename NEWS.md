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