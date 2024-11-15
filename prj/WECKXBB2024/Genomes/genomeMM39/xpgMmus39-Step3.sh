#!/bin/bash
#

# Using the NCBI EDirect utilities (https://www.nlm.nih.gov/dataguide/edirect/install.html),
#  we download the following rRNA gene accessions from NCBI Nucleotide into file
#  "Mmus_rRNA.fa" (used for screening out rRNA tags in the ChIPseq and RNAseq data, with
#  tagdust):
#
# V00849.1 Mouse gene for 45S ribosomal RNA
# X00525.1 Mouse 28S ribosomal RNA
# X00686.1 Mouse gene for 18S rRNA
# J01871.1 Mouse 5.8S ribosomal RNA
# K02235.1 Mouse 5S rRNA gene M4
#
esearch -db nuccore -query "V00849.1[ACCN] or X00525.1[ACCN] or X00686.1[ACCN] or J01871.1[ACCN] or K02235.1[ACCN]" \
	| efetch -format fasta > Mmus_rRNA.fa
