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
counts_data <- read.table(file = '../RNA_seq/1_CtbpKD/4_DGE_Analysis/v1a_CtbpKD_raw_counts.tsv', row.names = 1, sep = '\t', header = TRUE)

#Assigning a variable for the samples to conditions matrix
cat("Reading the samples to conditions table.\n\n\n")
colData <- read.table(file = '../RNA_seq/1_CtbpKD/4_DGE_Analysis/v1b_CtbpKD_samples_to_conditions.tsv', row.names = 1, sep = '\t', header = TRUE)

#Check the column names if needed
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

#Remove NA values as they affect upcoming plots badly
cat("Removing the NA values from the results table.\n\n\n")
res2 <- na.omit(res)

#Keep the upregulated or downregulated genes with a typical log2FoldChange threshold
cat("Filtering out the genes that do not cross the fold change value [AKA res3].\n\n\n")
res3 <- res2[res2$log2FoldChange > 1.0 | res2$log2FoldChange < -1.0 ,]

#After the log2FoldChange threshold, keep the genes with a FDR value of 0.01 or below (We could use p value threshold instead by replacing the padj with pvalue)
cat("Filtering out the genes that cross the false adjusted P value threshold [AKA res4].\n\n\n")
res4 <- res3[res3$padj < 0.01 ,]

### creating dot plots for genes of interest
# Load genes from a text file
cmg_genes <- readLines("./a2a_CTBP2_cmg_direct_targets_sorted.txt")

# Subset the res4 matrix for the selected genes
subset_cmg_genes <- res4[rownames(res4) %in% cmg_genes, ]

# Create a new column for coloring based on log2FoldChange direction
subset_cmg_genes$color <- ifelse(subset_cmg_genes$log2FoldChange > 0, "red", "blue")

# Sort the genes by log2FoldChange (upregulated genes at the top, downregulated at the bottom)
subset_cmg_genes <- subset_cmg_genes[order(subset_cmg_genes$log2FoldChange, decreasing = TRUE), ]

subset_cmg_genes$regulation <- ifelse(subset_cmg_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

# Dot plot with sorted gene names on the y-axis and dot size/color based on log2 fold change
jpeg(filename = "./3_Heat_Maps_and_Plots/Fig_2E.jpeg", res = 820, width = 3500, height = 7800)
ggplot(subset_cmg_genes, aes(x = 1, y = reorder(rownames(subset_cmg_genes), log2FoldChange))) +
    geom_point(aes(size = abs(log2FoldChange), color = regulation)) +  # Map color to regulation status
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +  # Set colors for the legend
    scale_size_continuous(name = "Log2 Fold Change") +
    labs(x = NULL, y = "Cell migration genes") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Hide x-axis text since it's just a placeholder
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 16),  # Bold and increase font size of y-axis labels
          legend.title = element_text(face = "bold"),  # Make legend title bold
          legend.text = element_text(face = "bold")) +  # Make legend text bold
    guides(color = guide_legend(title = "Gene Regulation"))  # Title for the color legend


# Load genes from a text file
ngg_genes <- readLines("./a2b_CTBP2_ngg_direct_targets_sorted.txt")

# Subset the res4 matrix for the selected genes
subset_ngg_genes <- res4[rownames(res4) %in% ngg_genes, ]

# Create a new column for coloring based on log2FoldChange direction
subset_ngg_genes$color <- ifelse(subset_ngg_genes$log2FoldChange > 0, "red", "blue")

# Sort the genes by log2FoldChange (upregulated genes at the top, downregulated at the bottom)
subset_ngg_genes <- subset_ngg_genes[order(subset_ngg_genes$log2FoldChange, decreasing = TRUE), ]

subset_ngg_genes$regulation <- ifelse(subset_ngg_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

# Dot plot with sorted gene names on the y-axis and dot size/color based on log2 fold change
jpeg(filename = "./3_Heat_Maps_and_Plots/Fig_S2C.jpeg", res = 800, width = 3000, height = 8000)
ggplot(subset_ngg_genes, aes(x = 1, y = reorder(rownames(subset_ngg_genes), log2FoldChange))) +
    geom_point(aes(size = abs(log2FoldChange), color = regulation)) +  # Map color to regulation status
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +  # Set colors for the legend
    scale_size_continuous(name = "Log2 Fold Change") +
    labs(x = NULL, y = "Neurogenesis genes") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Hide x-axis text since it's just a placeholder
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 14),  # Bold and increase font size of y-axis labels
          legend.title = element_text(face = "bold"),  # Make legend title bold
          legend.text = element_text(face = "bold")) +  # Make legend text bold
    guides(color = guide_legend(title = "Gene Regulation"))  # Title for the color legend

#create a normalized matrix and calculate the z-scores for an RNA-Seq heat map
cat("Nomalizing the counts and calculating their z-scores for the heat map.\n\n\n")
normalized_counts <- counts(dds, normalized = TRUE)
rlog_counts <- assay(rlog(dds))
normalized_counts2 <- normalized_counts[rownames(normalized_counts) %in% rownames(res4),]
normalized_counts3 <- t(scale(t(normalized_counts2)))

#Generate the heat map for the direct CtBP1/2 activation targets
cat("Generating the heat map for the CtBP1/2 activation targets.\n\n\n")
CtBP_actvation_targets_to_read <- readLines("./1_Intersections/CtBP2_activation_target_genes.txt")
CtBP_actvation_targets <- normalized_counts3[rownames(normalized_counts3) %in% CtBP_actvation_targets_to_read, ]
jpeg(filename = "./3_Heat_Maps_and_Plots/Fig_S4A.jpeg", res = 250, width = 3640, height = 4200)
Heatmap(as.matrix(CtBP_actvation_targets), show_row_names = FALSE, show_column_names = TRUE, cluster_rows = TRUE, cluster_columns = FALSE, row_dend_reorder = TRUE, column_dend_reorder = TRUE,  row_title = "Genes", column_title = "Samples", show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "Z-score"))



###Experiment design and differential gene expression analysis for PAX6

#Assigning a variable for the raw counts matrix
cat("Reading the raw counts table.\n\n\n")
counts_data <- read.table(file = '../RNA_seq/2_Prdm16KOinvivo/4_DGE_Analysis/1_PAX6/v1a_Prdm16KO_PPC_raw_counts.tsv', row.names = 1, sep = '\t', header = TRUE)

#Assigning a variable for the samples to conditions matrix
cat("Reading the samples to conditions table.\n\n\n")
colData <- read.table(file = '../RNA_seq/2_Prdm16KOinvivo/4_DGE_Analysis/1_PAX6/v1b_Prdm16KO_PPC_samples_to_conditions.tsv', row.names = 1, sep = '\t', header = TRUE)

#Check the column names if needed
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

#To make sure that the control condition is the reference (otherwise, if the Knockdown condition alphabetically comes before CTRL, the knockout condition would be considered the reference and the volcano plot would look inverted)
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

#Remove NA values as they affect upcoming plots badly
cat("Removing the NA values from the results table.\n\n\n")
res2 <- na.omit(res)

#Keep the upregulated or downregulated genes with a typical log2FoldChange threshold
cat("Filtering out the genes that do not cross the fold change value [AKA res3].\n\n\n")
res3 <- res2[res2$log2FoldChange > 0.6 | res2$log2FoldChange < -0.6 ,]

#After the log2FoldChange threshold, keep the genes with a FDR value of 0.01 or below (We could use p value threshold instead by replacing the padj with pvalue)
cat("Filtering out the genes that cross the false adjusted P value threshold [AKA res4].\n\n\n")
res4 <- res3[res3$pvalue < 0.05 ,]

### creating dot plots for genes of interest
# Load genes from a text file
cmg_genes <- readLines("./a2c_PRDM16_cmg_direct_targets_sorted.txt")

# Subset the res4 matrix for the selected genes
subset_cmg_genes <- res4[rownames(res4) %in% cmg_genes, ]

# Create a new column for coloring based on log2FoldChange direction
subset_cmg_genes$color <- ifelse(subset_cmg_genes$log2FoldChange > 0, "red", "blue")

# Sort the genes by log2FoldChange (upregulated genes at the top, downregulated at the bottom)
subset_cmg_genes <- subset_cmg_genes[order(subset_cmg_genes$log2FoldChange, decreasing = TRUE), ]

subset_cmg_genes$regulation <- ifelse(subset_cmg_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

# Dot plot with sorted gene names on the y-axis and dot size/color based on log2 fold change
jpeg(filename = "./3_Heat_Maps_and_Plots/Fig_S3B_migration.jpeg", res = 820, width = 3500, height = 7000)
ggplot(subset_cmg_genes, aes(x = 1, y = reorder(rownames(subset_cmg_genes), log2FoldChange))) +
    geom_point(aes(size = abs(log2FoldChange), color = regulation)) +  # Map color to regulation status
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +  # Set colors for the legend
    scale_size_continuous(name = "Log2 Fold Change") +
    labs(x = NULL, y = "Cell migration genes") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Hide x-axis text since it's just a placeholder
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 16),  # Bold and increase font size of y-axis labels
          legend.title = element_text(face = "bold"),  # Make legend title bold
          legend.text = element_text(face = "bold")) +  # Make legend text bold
    guides(color = guide_legend(title = "Gene Regulation"))  # Title for the color legend


# Load genes from a text file
ngg_genes <- readLines("./a2d_PRDM16_ngg_direct_targets_sorted.txt")

# Subset the res4 matrix for the selected genes
subset_ngg_genes <- res4[rownames(res4) %in% ngg_genes, ]

# Create a new column for coloring based on log2FoldChange direction
subset_ngg_genes$color <- ifelse(subset_ngg_genes$log2FoldChange > 0, "red", "blue")

# Sort the genes by log2FoldChange (upregulated genes at the top, downregulated at the bottom)
subset_ngg_genes <- subset_ngg_genes[order(subset_ngg_genes$log2FoldChange, decreasing = TRUE), ]

subset_ngg_genes$regulation <- ifelse(subset_ngg_genes$log2FoldChange > 0, "Upregulated", "Downregulated")

# Dot plot with sorted gene names on the y-axis and dot size/color based on log2 fold change
jpeg(filename = "./3_Heat_Maps_and_Plots/Fig_S3B_neurogenesis.jpeg", res = 800, width = 3000, height = 7500)
ggplot(subset_ngg_genes, aes(x = 1, y = reorder(rownames(subset_ngg_genes), log2FoldChange))) +
    geom_point(aes(size = abs(log2FoldChange), color = regulation)) +  # Map color to regulation status
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +  # Set colors for the legend
    scale_size_continuous(name = "Log2 Fold Change") +
    labs(x = NULL, y = "Neurogenesis genes") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Hide x-axis text since it's just a placeholder
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face = "bold", size = 14),  # Bold and increase font size of y-axis labels
          legend.title = element_text(face = "bold"),  # Make legend title bold
          legend.text = element_text(face = "bold")) +  # Make legend text bold
    guides(color = guide_legend(title = "Gene Regulation"))  # Title for the color legend
