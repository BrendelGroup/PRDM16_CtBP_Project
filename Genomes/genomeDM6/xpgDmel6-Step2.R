# 4. Build the bowtie2 index, required for bowtie2-aligning reads to the dm6 genome:
#
library(Rbowtie2)
bowtie2_build(references = "BSgenome.Dmelanogaster.UCSC.dm6chrs.fa",
              bt2Index = file.path("BSgenome.Dmelanogaster.UCSC.dm6chrs"),
	      "--threads 16 --quiet"
             )
