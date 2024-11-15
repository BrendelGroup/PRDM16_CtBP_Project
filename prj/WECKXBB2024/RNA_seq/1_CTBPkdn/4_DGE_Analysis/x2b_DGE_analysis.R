###Calling the required libraries

suppressPackageStartupMessages(library(rio))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggrepel))


###Experiment design and differential gene expression analysis

#Assigning a variable for the raw counts matrix
cat("Reading the raw counts table.\n\n\n")
counts_data <- read.table(file = './v1a_CTBPkdn_raw_counts.tsv', row.names = 1, sep = '\t', header = TRUE)

#Assigning a variable for the samples to conditions matrix
cat("Reading the samples to conditions table.\n\n\n")
colData <- read.table(file = './v1b_CTBPkdn_samples_to_conditions.tsv', row.names = 1, sep = '\t', header = TRUE)

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

#Keep the genes that have 10 or more reads mapped to it (This is useful for the PCA plot)
cat("Keeping the entries with 10 or more reads.\n\n\n")
keep <- rowSums(counts(dds)) >= 10

#Keep the genes that have 10 or more reads mapped to it in dds
cat("Applying that to the DESeq2 matrix.\n\n\n")
dds <- dds[keep, ]

#To make sure that the control condition is the reference (otherwise, if the Knockdown condition alphabetically comes before CTRL, the knockdown condition would be considered the reference and the volcano plot would look inverted)
cat("Making sure the control condition is the reference.\n\n\n")
dds$condition <- relevel(dds$condition, ref = "CTRL")

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
write.table(res,file="./v2a_CTBPkdn_DGE_Matrix.tsv")

#Remove NA values
cat("Removing the NA values from the results table.\n\n\n")
res2 <- na.omit(res)

#Export the DGE matrix without the NA values
cat("Exporting the results table without the NA values [AKA res2].\n\n\n")
write.table(res2,file="./v2b_CTBPkdn_DGE_Matrix_NA_omitted.tsv")

#Keep the upregulated or downregulated genes with a typical log2FoldChange threshold
cat("Filtering out the genes that do not cross the fold change value [AKA res3].\n\n\n")
res3 <- res2[res2$log2FoldChange > 1.0 | res2$log2FoldChange < -1.0 ,]

#After the log2FoldChange threshold, keep the genes with a FDR value of 0.01 or below (We could use p value threshold instead by replacing the padj with pvalue)
cat("Filtering out the genes that cross the false adjusted P value threshold [AKA res4].\n\n\n")
res4 <- res3[res3$padj < 0.01 ,]

#Keep the significantly upregulated genes
cat("Filtering out the significantly upregulated genes [AKA res5].\n\n\n")
res5 <- res4[res4$log2FoldChange > 1.0 ,]

#Keep the significantly downregulated genes
cat("Filtering out the significantly downregulated genes [AKA res6].\n\n\n")
res6 <- res4[res4$log2FoldChange < -1.0 ,]

#Export the DGE matrix without the NA values and with the genes that cross the log2FoldChange and FDR thresholds (in other words, the significantly differentially expressed genes [bidirectional, up, and down])
cat("Exporting the filtered results table.\n\n\n")
write.table(res4,file="./v2c_CTBPkdn_DGE_Matrix_NA_omitted_filtered.tsv")
write.table(res5,file="./v2d_CTBPkdn_DGE_Matrix_significantly_upReg.tsv")
write.table(res6,file="./v2e_CTBPkdn_DGE_Matrix_significantly_downReg.tsv")

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

#create a normalized matrix and calculate the z-scores for an RNA-Seq heat map
cat("Nomalizing the counts and calculating their z-scores for the heat map.\n\n\n")
normalized_counts <- counts(dds, normalized = TRUE)
rlog_counts <- assay(rlog(dds))
normalized_counts2 <- normalized_counts[rownames(normalized_counts) %in% rownames(res4),]
normalized_counts3 <- t(scale(t(normalized_counts2)))

#Generate the supplementary figure 2 heat map
cat("Generating the heat map for the CtBP1/2 genes.\n\n\n")
CTBP_Genes_to_extract <- c("Ctbp1","Ctbp2")
CTBP_genes <- normalized_counts3[rownames(normalized_counts3) %in% CTBP_Genes_to_extract, ]
jpeg(filename = "./v3b_CTBPkdn_heat_map_CtBP_genes.jpeg", width = 700, height = 200)
Heatmap(as.matrix(CTBP_genes), show_row_names = TRUE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, row_dend_reorder = TRUE, column_dend_reorder = TRUE,  row_title = "Genes", column_title = "Samples", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Z-score"), row_title_gp = gpar(fontsize = 18), column_title_gp = gpar(fontsize = 18), row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16))

###Generating a volcanoplot using ggplot
#Giving a pseudo FDR value for the genes that have an FDR value of 0 (otherwise they will be plotted to infinity on the y-axis as the log2 of a zero value is not defined)
cat("Making pseudo adjusted P values for genes with an adjusted P value of zero.\n\n\n")
resppa <- res
min_nonzero_padj <- min(resppa$padj[resppa$padj > 0], na.rm = TRUE)
pseudo_padj <- min_nonzero_padj * 0.1
resppa$padj[resppa$padj == 0] <- pseudo_padj
resppa2 <- na.omit(resppa)

#Create a new matrix that we will edit for ggplot and the volcano plot
vpDFPath <- resppa2

#Give the column of gene names a name instead of "X"
vpDFPath <- cbind(rownames(vpDFPath), data.frame(vpDFPath, row.names=NULL))
colnames(vpDFPath)[1] <- "gene_id"

#Create a new column and give values of "NO", "UP", or "DOWN" based on how the genes cross the thresholds
vpDFPath$diffexpressed <- "NO"
vpDFPath$diffexpressed[vpDFPath$log2FoldChange > 1.0 & vpDFPath$padj < 0.01] <- "UP"
vpDFPath$diffexpressed[vpDFPath$log2FoldChange < -1.0 & vpDFPath$padj < 0.01] <- "DOWN"


#Generate a volcano plot without labels
cat("Generating a volcano plot without labels.\n\n\n")
jpeg(filename = "./v4a_CTBPkdn_Volcano_Plot_no_labels.jpeg", res = 200, width = 2500, height = 2000)
ggplot(data = vpDFPath, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) + geom_vline(xintercept = c(-1.0,1.0), col = "gray", linetype = 'dashed') + geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + geom_point(size = 1.5) + theme_classic(base_size = 20) + theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'), axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'), plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c("blue", "grey", "red"), labels = c("Downregulated", "Not significant", "Upregulated")) + coord_cartesian(ylim = c(0, 220), xlim = c(-14, 14)) + labs(color = 'Significance', x = expression("log"[2]*"(FC)"), y = expression("-log"[10]*"(p-adj)")) + scale_x_continuous(breaks = seq(-12, 12, 2))


#Generate a volcano plot with Ctbp1/2 labels
cat("Generating a volcano plot with Ctbp1/2 labels.\n\n\n")
vpDFPath$delabel <- ifelse(vpDFPath$gene_id %in% c("Ctbp1"), vpDFPath$gene_id, NA)
jpeg(filename = "./v4b_CTBPkdn_Volcano_Plot_with_labels_Arrow.jpeg", res = 200, width = 2500, height = 2000)
ggplot(data = vpDFPath, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) + geom_vline(xintercept = c(-1.0,1.0), col = "gray", linetype = 'dashed') + geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + geom_point(size = 1.5) + theme_classic(base_size = 20) + theme(axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'), axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'), plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = c("blue", "grey", "red"), labels = c("Downregulated", "Not significant", "Upregulated")) + coord_cartesian(ylim = c(0, 220), xlim = c(-14, 14)) + labs(color = 'Significance', x = expression("log"[2]*"(FC)"), y = expression("-log"[10]*"(p-adj)")) + scale_x_continuous(breaks = seq(-12, 12, 2)) + geom_text_repel(max.overlaps = Inf , size = 2.2, vjust = -1) + annotate("text", x = -6, y = 175, label = sum(vpDFPath$diffexpressed == "DOWN"), size = 5, color = "blue") + annotate("text", x = -6, y = 165, label = "Downregulated genes", size = 4, color = "blue") + annotate("text", x = 6, y = 175, label = sum(vpDFPath$diffexpressed == "UP"), size = 5, color = "red") + annotate("text", x = 6, y = 165, label = "Upregulated genes", size = 4, color = "red") + geom_segment(data = subset(vpDFPath, gene_id == "Ctbp2"), aes(x = log2FoldChange, y = -log10(padj), xend = log2FoldChange - 3, yend = -log10(padj) + 20),arrow = arrow(length = unit(0.2, "cm")), color = "black")
