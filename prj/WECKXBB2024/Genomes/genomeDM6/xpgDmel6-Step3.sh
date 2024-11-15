#!/bin/bash
#

# This script gets the dm6 annotation file from NCBI and adjusts the sequence names to UCSC style:
#
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz
wget -O - https://hgdownload.cse.ucsc.edu/goldenPath/dm6/database/chromAlias.txt.gz | gunzip > dm6_chromAlias.txt


if ! command -v chromToUcsc &> /dev/null
then
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc
	chmod a+x chromToUcsc
	./chromToUcsc -a dm6_chromAlias.txt -i GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz -o dm6.gff
else
	chromToUcsc -a dm6_chromAlias.txt -i GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz -o dm6.gff
fi

# We also create a table of gene name versus GeneID entries:
#
grep -E "	gene	" dm6.gff | cut -f9 | sed -e "s/,/;/" | cut -d";" -f1,2 \
	| sed -e "s/^ID=gene-//" | sed -e "s/;Dbxref=FLYBASE:/	/" > dm6GeneID

n=`cat dm6GeneID | wc -l`
echo -e "\nHead of dm6GeneID ($n entries):\n"
head -5 dm6GeneID


# Using the NCBI EDirect utilities (https://www.nlm.nih.gov/dataguide/edirect/install.html),
#  we download the following rRNA gene accessions from NCBI Nucleotide into file
#  "Dmel_rRNA.fa" (used for screening out rRNA tags in the ChIPseq and RNAseq data, with
#  tagdust):
#
# NR_133562.1 Drosophila gene for 29S ribosomal RNA
# NR_133599.1 Drosophila gene for 18S ribosomal RNA
# NR_001870.1 Drosophila gene for 5S ribosomal RNA
# NR_133552.1 Drosophila gene for 2S ribosomal RNA
#
esearch -db nuccore -query "NR_133562 [ACCN] or NR_133559 [ACCN] or NR_001870 [ACCN] or NR_133552 [ACCN]" \
	| efetch -format fasta > Dmel_rRNA.fa
