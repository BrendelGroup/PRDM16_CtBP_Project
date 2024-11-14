# 1. Load the BSgenome object for Drosophila (assembly dm6):
#
library(BSgenome.Dmelanogaster.UCSC.dm6)

write("\n\nLOG: The downloaded BSgenome object BSgenome.Dmelanogaster.UCSC.dm6 is as follows:\n", stdout())
BSgenome.Dmelanogaster.UCSC.dm6


# 2. We create an object for a Drosophila chromosome-only reference genome:
#
Chromosomes <- paste0("chr", c("2L", "2R", "3L", "3R", "4", "X", "Y", "M"))
ChrSeq <- lapply(Chromosomes, function(x) BSgenome.Dmelanogaster.UCSC.dm6[[x]])
names(ChrSeq) <- Chromosomes
ChrSeqSet <- DNAStringSet(ChrSeq)
write("\n\nLOG: Extracting only the full-length chromosomes gives:\n", stdout())
ChrSeqSet


# 3. Save the sequences to a genome file on disk:
#
write("\n\nLOG: The FASTA-formatted genome file is saved as BSgenome.Dmelanogaster.UCSC.dm6chrs.fa.", stdout())
writeXStringSet(ChrSeqSet, "BSgenome.Dmelanogaster.UCSC.dm6chrs.fa")
