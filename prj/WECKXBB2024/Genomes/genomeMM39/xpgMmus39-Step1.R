# 1. Load the BSgenome object for mouse (assembly mm39):
#
library(BSgenome.Mmusculus.UCSC.mm39)

write("\n\nLOG: The downloaded BSgenome object BSgenome.Mmusculus.UCSC.mm39 is as follows:\n", stdout())
BSgenome.Mmusculus.UCSC.mm39


# 2. We create an object for a mouse chromosome-only reference genome:
#
Chromosomes <- paste0("chr", c(1:19, "X", "Y", "M"))
ChrSeq <- lapply(Chromosomes, function(x) BSgenome.Mmusculus.UCSC.mm39[[x]])
names(ChrSeq) <- Chromosomes
ChrSeqSet <- DNAStringSet(ChrSeq)
write("\n\nLOG: Extracting only the full-length chromosomes gives:\n", stdout())
ChrSeqSet


# 3. Save the sequences to a genome file on disk:
#
write("\n\nLOG: The FASTA-formatted genome file is saved as BSgenome.Mmusculus.UCSC.mm39chrs.fa.", stdout())
writeXStringSet(ChrSeqSet, "BSgenome.Mmusculus.UCSC.mm39chrs.fa")
