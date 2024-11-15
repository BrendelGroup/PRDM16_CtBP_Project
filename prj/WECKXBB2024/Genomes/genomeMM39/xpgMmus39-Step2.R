# 4. Build the bowtie2 index, required for bowtie2-aligning reads to the mm39 genome:
#
library(Rbowtie2)
bowtie2_build(references = "BSgenome.Mmusculus.UCSC.mm39chrs.fa",
              bt2Index = file.path("BSgenome.Mmusculus.UCSC.mm39chrs"),
	      "--threads 16 --quiet"
             )
