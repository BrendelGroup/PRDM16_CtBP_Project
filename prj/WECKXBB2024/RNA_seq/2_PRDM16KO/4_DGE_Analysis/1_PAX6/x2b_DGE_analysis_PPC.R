###Calling the required libraries

suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))


###Experiment design and differential gene expression analysis

#Assigning a variable for the raw counts matrix
cat("Reading the raw counts table.\n\n\n")
counts_data <- read.table(file = './v1a_PRDM16_PPC_raw_counts.tsv', row.names = 1, sep = '\t', header = TRUE)

#Assigning a variable for the samples to conditions matrix
cat("Reading the samples to conditions table.\n\n\n")
colData <- read.table(file = './v1b_PRDM16_PPC_samples_to_conditions.tsv', row.names = 1, sep = '\t', header = TRUE)

#Check the column names
cat("The column names are:\n")
colnames(counts_data)
cat("\n\n\n")

#This should return TRUE if the column names of the raw counts matrix is the same row names of the samples to conditions (it has to be this way)
cat("The column names in the raw counts matrix are the row names in the samples to conditions matrix:")
all(colnames(counts_data) %in% rownames(colData)) 
cat("\n\n")

#This should return TRUE if the column names of the raw counts matrix is in the same order of the row names of the samples to conditions (it has to be this way)
cat("The column names in the raw counts matrix are in the same order as the row names in the samples to conditions matrix:")
all(colnames(counts_data) == rownames(colData))
cat("\n\n\n")

#Setting up the matrix for DESeq2
cat("Setting up the matrix for DESeq2.\n")
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ condition)
cat("\n\n\n")

#Checking it out
cat("Done. Let's check it out:\n\n")
dds
cat("\n\n\n")

#Keep the genes that have 10 or more reads mapped to it
cat("Keeping the entries with 10 or more reads.\n\n\n")
keep <- rowSums(counts(dds)) >= 10

#Keep the genes that have 10 or more reads mapped to it in dds
cat("Applying that to the DESeq2 matrix.\n\n\n")
dds <- dds[keep, ]

#To make sure that the control condition is the reference (otherwise, if the Knockdown condition alphabetically comes before CTRL, the knockdown condition would be considered the reference and the volcano plot would look inverted)
cat("Making sure the control condition is the reference.\n\n\n")
dds$condition <- relevel(dds$condition, ref = "WT")

#Check the list of conditions
cat("Checking out the list of conditions:\n")
dds$condition
cat("\n\n\n")

#Perform the differential gene expression analysis
cat("Performing the differntial gene expression analysis:\n\n")
dds <- DESeq(dds)
cat("\n\n\n")

#Check it out
cat("Done. Let's check it out:\n\n")
dds
cat("\n\n\n")

#Get the differential gene expression matrix that will be used for some of the upcomming plots
cat("Moving the results to a data frame [res].\n\n\n")
res <- results(dds)

#Get a summary on the matrix
cat("The summary res:\n")
summary(res)
cat("\n\n\n")

#Export the DGE matrix
cat("Exporting the results table.\n\n\n")
write.table(res,file="./v2a_PRDM16_PPC_DGE_Matrix.tsv")

#Remove NA values as they affect upcoming plots badly
cat("Removing the NA values from the results table.\n\n\n")
res2 <- na.omit(res)

#Export the DGE matrix without the NA values
cat("Exporting the results table without the NA values [AKA res2].\n\n\n")
write.table(res2,file="./v2b_PRDM16_PPC_DGE_Matrix_NA_omitted.tsv")

#Keep the upregulated or downregulated genes with a typical log2FoldChange threshold
cat("Filtering out the genes that do not cross the fold change value [AKA res3].\n\n\n")
res3 <- res2[res2$log2FoldChange > 0.6 | res2$log2FoldChange < -0.6 ,]

#After the log2FoldChange threshold, keep the genes with a p value of 0.05 or below
cat("Filtering out the genes that cross the false adjusted P value threshold [AKA res4].\n\n\n")
res4 <- res3[res3$pvalue < 0.05 ,]

#Keep the significantly upregulated genes
cat("Filtering out the significantly upregulated genes [AKA res5].\n\n\n")
res5 <- res4[res4$log2FoldChange > 0.6 ,]

#Keep the significantly downregulated genes
cat("Filtering out the significantly downregulated genes [AKA res6].\n\n\n")
res6 <- res4[res4$log2FoldChange < -0.6 ,]

#Export the DGE matrix without the NA values and with the genes that cross the log2FoldChange and FDR thresholds (in other words, the significantly differentially expressed genes [bidirectional, up, and down])
cat("Exporting the filtered results table.\n\n\n")
write.table(res4,file="./v2c_PRDM16_PPC_DGE_Matrix_NA_omitted_filtered.tsv")
write.table(res5,file="./v2d_PRDM16_PPC_DGE_Matrix_significantly_upReg.tsv")
write.table(res6,file="./v2e_PRDM16_PPC_DGE_Matrix_significantly_downReg.tsv")

#Check how many rows and columns are there in each matrix
cat("The dimensions of res :")
dim(res)
cat("The dimensions of res2:")
dim(res2)
cat("The dimensions of res3:")
dim(res3)
cat("The dimensions of res4:")
dim(res4)
cat("The dimensions of res5:")
dim(res5)
cat("The dimensions of res6:")
dim(res6)
cat("\n\n\n")
