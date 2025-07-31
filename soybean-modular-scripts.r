## ===== 00_functions.R =====
# Custom functions for soybean RNA-seq analysis
# Load this file first in each script

# Function to create directory if it doesn't exist
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("Created directory:", path, "\n")
  }
}

# Function to save checkpoint with timestamp
save_checkpoint <- function(file_name, objects = NULL) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  versioned_name <- sub("\\.RData$", paste0("_", timestamp, ".RData"), file_name)
  
  if (is.null(objects)) {
    save.image(file_name)
    cat("Saved workspace to:", file_name, "\n")
  } else {
    save(list = objects, file = file_name)
    cat("Saved objects to:", file_name, "\n")
  }
  
  # Also save versioned copy
  file.copy(file_name, versioned_name)
  cat("Versioned copy saved as:", versioned_name, "\n")
}

# Function to load checkpoint
load_checkpoint <- function(file_name) {
  if (file.exists(file_name)) {
    load(file_name, envir = .GlobalEnv)
    cat("Loaded checkpoint:", file_name, "\n")
  } else {
    stop("Checkpoint file not found:", file_name)
  }
}

## ===== 01_setup_and_import.R =====
#-------------------------------------------------------#
#     01. Setup and Data Import                         #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#


# Load required libraries
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(rtracklayer)
  library(magrittr)
  library(ggplot2)
})

# Set working directory (modify as needed)
# setwd("/path/to/your/project")

# Create directory structure
create_dir("data")
create_dir("results")
create_dir("results/figures")
create_dir("results/tables")
create_dir("results/checkpoints")

# Read sample information
cat("\nReading sample information...\n")
targets <- readTargets("Targets_Final.txt")
targets <- targets[1:60,]

# Add group information
targets$group <- paste(targets$Line, targets$Timepoint, sep = "_")
targets$GpF <- factor(targets$group, levels = c(
  "PI532462A_TP0","PI532462A_TP1","PI532462A_TP2","PI532462A_TP3","PI532462A_TP4",
  "LD112170_TP0","LD112170_TP1","LD112170_TP2","LD112170_TP3","LD112170_TP4",
  "PI612713B_TP0","PI612713B_TP1","PI612713B_TP2","PI612713B_TP3","PI612713B_TP4",
  "PI547745_TP0","PI547745_TP1","PI547745_TP2","PI547745_TP3","PI547745_TP4"
))
targets$col <- as.numeric(targets$GpF)
targets$Batch <- as.factor(targets$Batch)

# Summary statistics
cat("\nDataset summary:\n")
cat("Total samples:", nrow(targets), "\n")
cat("Total reads:", sum(targets$NumReads)/1e9, "billion\n")
cat("Leaf types:", unique(targets$Leaf_type), "\n")
cat("Timepoints:", unique(targets$Timepoint), "\n")
cat("Lines:", unique(targets$Line), "\n")

# Load Salmon count data
cat("\nLoading Salmon count data...\n")
temp <- load("SalmonSummarizedOutput.RData")
cat("Loaded objects:", temp, "\n")

# Fix sample names
rownames(meta_info) <- gsub("salmo", "", rownames(meta_info))

# Match sample order
temp_order <- match(targets$Sample, rownames(meta_info))
meta_info <- meta_info[temp_order,]
tx.all$abundance <- tx.all$abundance[,temp_order]
tx.all$counts <- tx.all$counts[,temp_order]
tx.all$length <- tx.all$length[,temp_order]

# Update column names
rownames(meta_info) <- targets$Label
colnames(tx.all$abundance) <- targets$Label
colnames(tx.all$counts) <- targets$Label
colnames(tx.all$length) <- targets$Label

# Load annotation for transcript to gene mapping
cat("\nLoading gene annotation...\n")
annotation <- read.delim("Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")
all_tx_info <- annotation[,c("transcriptName", "locusName")]

# Check transcript coverage
cat("Transcripts in annotation:", nrow(all_tx_info), "\n")
cat("Transcripts in data:", nrow(tx.all$counts), "\n")
cat("Matching transcripts:", sum(rownames(tx.all$counts) %in% all_tx_info$transcriptName), "\n")

# Summarize to gene level
cat("\nSummarizing to gene level...\n")
library(tximport)
gene_all <- tximport::summarizeToGene(tx.all, all_tx_info, countsFromAbundance = "lengthScaledTPM")

# Create DGEList
d_salmon <- DGEList(gene_all$counts, samples = targets)
d_salmon$genes <- data.frame(GeneID = rownames(d_salmon$counts))

cat("\nGenes in dataset:", nrow(d_salmon$counts), "\n")
cat("Samples in dataset:", ncol(d_salmon$counts), "\n")

# Save checkpoint
save_checkpoint("results/checkpoints/01_data_imported.RData")

cat("\n=== Data import complete ===\n")

## ===== 02_quality_control.R =====
#-------------------------------------------------------#
#     02. Quality Control Analysis                      #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint from previous step
cat("Loading data checkpoint...\n")
load_checkpoint("results/checkpoints/01_data_imported.RData")

# Load additional libraries
library(gplots)

# Create QC output directory
create_dir("results/figures/QC")

## Read Fate Analysis
cat("\nCalculating read fate statistics...\n")
ReadFateSalmon <- data.frame(mapped = meta_info$num_mapped) / targets$NumReads * 100
ReadFateSalmon$Label <- factor(targets$Label, levels = unique(targets$Label), ordered = TRUE)
ReadFateSalmon$Group <- factor(targets$GpF, levels = unique(targets$GpF), ordered = TRUE)

# Summary statistics
cat("\nMapping rate summary:\n")
print(summary(ReadFateSalmon$mapped))

# Plot read fates
tiff("results/figures/QC/Salmon_ReadFates.tiff", res = 300, unit = "in", height = 6, width = 8)
barplot(ReadFateSalmon$mapped,
        ylim = c(0,100),
        ylab = "Percent of mapped reads",
        xlab = "Samples",
        las = 2, col = targets$Batch,
        main = "Salmon Read Mapping Rates", 
        cex.axis = 0.8,
        names.arg = targets$Label)
dev.off()

## Batch Effect Visualization
cat("\nVisualizing batch effects...\n")

# Prepare batch data
col_num <- which(colnames(targets) %in% c("Label", "Batch", "Assigned_to_gene", "Total_accounted", "Leaf_type"))
batch <- targets[,col_num]
batch.sort <- batch[order(batch$Batch, decreasing = FALSE),]
batch.sort$Total_accounted <- as.numeric(batch.sort$Total_accounted/1000000)

# Plot batch effects
tiff("results/figures/QC/Batch_Effect_Total_reads.tiff", res = 300, unit = "in", height = 6, width = 8)
barplot(batch.sort$Total,
        ylim = c(0,50),
        ylab = "Millions of reads",
        xlab = "",
        las = 2, col = batch.sort$Batch,
        main = "Batch effect - Read counts", 
        cex.axis = 0.8,
        names.arg = batch.sort$Label)
legend("topright", legend = c("Batch1", "Batch2"), fill = c("black", "red"))
dev.off()

## Raw Count Distributions
cat("\nPlotting count distributions...\n")

# Log transform for visualization
raw_log <- log2(d_salmon$counts + 0.1)

tiff("results/figures/QC/Log_Counts_Distribution.tiff", res = 300, unit = "in", height = 6, width = 6)
plotDensities(raw_log, group = d_salmon$samples$GpF, 
              col = 1:4, main = "log2(raw counts + 0.1)", legend = FALSE)
dev.off()

## Library Sizes
tiff("results/figures/QC/Library_sizes.tiff", res = 300, unit = "in", height = 6, width = 6)
barplot(d_salmon$samples$lib.size / 1e6,
        ylim = c(0,50),
        ylab = "Total number of reads (million)",
        xlab = "Samples",
        las = 2, col = d_salmon$samples$col,
        main = "Library sizes", 
        cex.axis = 0.8,
        names.arg = targets$Label)
dev.off()

## Initial Clustering
cat("\nPerforming initial clustering...\n")

tiff("results/figures/QC/MDS_plot_raw.tiff", res = 300, unit = "in", height = 6, width = 6)
plotMDS(d_salmon, main = "MDS plot - Raw data",
        col = d_salmon$samples$col,
        top = 5000, prior.count = 2)
dev.off()

# Interactive MDS plot
library(Glimma)
glMDSPlot(d_salmon, top = 5000, labels = d_salmon$samples$Label,
          groups = d_salmon$samples[,c("Leaf_type","Timepoint","Rep","group", "Batch")],
          folder = "results/figures/QC/glimma",
          html = "MDS-Plot_raw",
          launch = FALSE, prior.count = 2)

## Housekeeping Gene Analysis
cat("\nAnalyzing housekeeping genes...\n")

library(readxl)
hk.genes <- read_excel("Housekeeping_genes_Machado_et_al_2019.xlsx", sheet = 1)
hk.genes$Gene <- as.character(hk.genes$Gene)

# Check how many housekeeping genes are in our dataset
hk_present <- sum(rownames(d_salmon$counts) %in% hk.genes$Gene)
cat("Housekeeping genes found:", hk_present, "out of", nrow(hk.genes), "\n")

# Save QC results
qc_results <- list(
  read_fate = ReadFateSalmon,
  mapping_summary = summary(ReadFateSalmon$mapped),
  library_sizes = d_salmon$samples$lib.size,
  batch_info = batch.sort,
  hk_genes_found = hk_present
)

save(qc_results, file = "results/tables/qc_results.RData")

# Save checkpoint
save_checkpoint("results/checkpoints/02_qc_complete.RData")

cat("\n=== Quality control complete ===\n")



#-------------------------------------------------------#
#     03. Batch Effect Correction                       #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading QC checkpoint...\n")
load_checkpoint("results/checkpoints/02_qc_complete.RData")

# Load additional libraries
library(sva)
library(Biobase)
library(pvca)

# Create output directory
create_dir("results/figures/batch_correction")

## Gene Filtering (needed before batch correction)
cat("\nFiltering genes...\n")

# Calculate CPM values
cpm_values <- cpm(d_salmon$counts)

# Filter genes with at least 0.5 CPM in at least 3 samples
above0.5cpm <- rowSums(cpm_values >= 0.5)
keep_genes <- above0.5cpm >= 3

cat("Genes before filtering:", nrow(d_salmon), "\n")
cat("Genes after filtering:", sum(keep_genes), "\n")
cat("Proportion kept:", mean(keep_genes), "\n") #0.7867202 

# Apply filtering
d_filt <- d_salmon[keep_genes, , keep.lib.sizes = FALSE]

# Check proportion of reads kept
cat("Proportion of reads kept:", mean(d_filt$samples$lib.size / d_salmon$samples$lib.size), "\n")

## Combat-Seq Batch Correction
cat("\nPerforming ComBat-seq batch correction...\n")

# Extract count matrix
count_matrix <- d_filt$counts

# Get batch information
batch <- d_filt$samples$Batch

# Apply ComBat-seq
adjusted_count_matrix <- ComBat_seq(count_matrix, batch = batch, group = NULL)

# Create new DGEList with adjusted counts
d_filt_combatseq <- DGEList(adjusted_count_matrix, samples = targets)

## TMM Normalization
cat("\nPerforming TMM normalization...\n")

d_filt_combatseq <- calcNormFactors(d_filt_combatseq)

# Plot normalization factors
tiff("results/figures/batch_correction/TMM_normalization_factors.tiff", 
     res = 300, unit = "in", height = 6, width = 6)
barplot(d_filt_combatseq$samples$norm.factors, col = targets$col,
        ylim = c(0,1.5),
        las = 2,
        ylab = "TMM norm values",
        main = "TMM normalization factors",
        names.arg = targets$Label)
abline(h = 1)
dev.off()

## PVCA Analysis - Before and After
cat("\nPerforming PVCA analysis...\n")

# Calculate logCPM for both datasets
logCPM_before <- cpm(d_filt, log = TRUE, prior.count = 2)
logCPM_combatseq <- cpm(d_filt_combatseq, log = TRUE, prior.count = 2)

# BEFORE batch correction
pData_before <- d_filt$samples[,c("Batch", "Leaf_type", "Timepoint")]
myExpressionSet_before <- ExpressionSet(
  assayData = as.matrix(logCPM_before),
  phenoData = AnnotatedDataFrame(pData_before)
)

pvcaObj_before <- pvcaBatchAssess(myExpressionSet_before, 
                                  c("Batch", "Leaf_type", "Timepoint"), 0.6)

# AFTER batch correction
pData <- d_filt_combatseq$samples[,c("Batch", "Leaf_type", "Timepoint")]
myExpressionSet <- ExpressionSet(
  assayData = as.matrix(logCPM_combatseq),
  phenoData = AnnotatedDataFrame(pData)
)

pvcaObj <- pvcaBatchAssess(myExpressionSet, 
                           c("Batch", "Leaf_type", "Timepoint"), 0.6)

# Create comparison plot with value annotations
tiff("results/figures/batch_correction/PVCA_comparison.tiff", 
     res = 300, units = "in", height = 6, width = 12)
par(mfrow = c(1, 2), mar = c(8, 4, 4, 2))

# Function for coloring
effect_colors <- function(labels) {
  colors <- rep("gray", length(labels))
  colors[grep("Batch", labels)] <- "red"
  colors[labels == "Timepoint"] <- "darkblue"
  colors[labels == "Leaf_type"] <- "blue"
  colors[grep(":", labels)] <- "orange"
  colors[labels == "resid"] <- "darkgray"
  return(colors)
}

# Before correction
bp1 <- barplot(pvcaObj_before$dat, 
               names.arg = pvcaObj_before$label,
               main = "Before ComBat-seq",
               ylab = "Proportion of Variance",
               ylim = c(0, max(c(pvcaObj_before$dat, pvcaObj$dat)) * 1.2),
               las = 2,
               cex.names = 0.8,
               col = effect_colors(pvcaObj_before$label))

# Add value annotations
text(x = bp1, y = pvcaObj_before$dat + 0.02, 
     labels = sprintf("%.1f%%", pvcaObj_before$dat * 100),
     cex = 0.7, pos = 3)

# After correction
bp2 <- barplot(pvcaObj$dat, 
               names.arg = pvcaObj$label,
               main = "After ComBat-seq",
               ylab = "Proportion of Variance",
               ylim = c(0, max(c(pvcaObj_before$dat, pvcaObj$dat)) * 1.2),
               las = 2,
               cex.names = 0.8,
               col = effect_colors(pvcaObj$label))

# Add value annotations
text(x = bp2, y = pvcaObj$dat + 0.02, 
     labels = sprintf("%.1f%%", pvcaObj$dat * 100),
     cex = 0.7, pos = 3)

dev.off()

## MDS Plot After Correction
tiff("results/figures/batch_correction/MDS_plot_combatseq.tiff", 
     res = 300, unit = "in", height = 6, width = 6)
plotMDS(d_filt_combatseq, main = "MDS plot after ComBat-seq",
        col = d_filt_combatseq$samples$col,
        top = 5000, prior.count = 2)
dev.off()

# Save batch correction results
batch_results <- list(
  genes_filtered = sum(keep_genes),
  proportion_kept = mean(keep_genes),
  pvca_before = pvcaObj_before,
  pvca_after = pvcaObj,
  batch_variance_before = sum(pvcaObj_before$dat[grep("Batch", pvcaObj_before$label)]),
  batch_variance_after = sum(pvcaObj$dat[grep("Batch", pvcaObj$label)])
)

cat("\n=== Batch correction complete ===\n")
cat("Batch variance reduced from", 
    round(batch_results$batch_variance_before * 100, 1), 
    "% to", 
    round(batch_results$batch_variance_after * 100, 1), "%\n")

## Correlation Analysis Between Batches - Before and After Correction
cat("\n========== CORRELATION ANALYSIS ==========\n")

# Identify batch samples
unique_batches <- unique(d_filt$samples$Batch)
batch1_samples <- which(d_filt$samples$Batch == unique_batches[1])
batch2_samples <- which(d_filt$samples$Batch == unique_batches[2])

cat("Batch 1 samples:", length(batch1_samples), "\n")
cat("Batch 2 samples:", length(batch2_samples), "\n")

# BEFORE correction
logCounts_before <- log2(d_filt$counts + 1)

# Within-batch correlations
cor_within_batch1_before <- cor(logCounts_before[, batch1_samples])
mean_within_batch1_before <- mean(cor_within_batch1_before[upper.tri(cor_within_batch1_before)])

cor_within_batch2_before <- cor(logCounts_before[, batch2_samples])
mean_within_batch2_before <- mean(cor_within_batch2_before[upper.tri(cor_within_batch2_before)])

# Between-batch correlations
cor_between_batches_before <- cor(logCounts_before[, batch1_samples], 
                                  logCounts_before[, batch2_samples])
mean_between_batches_before <- mean(cor_between_batches_before)

cat("\nBEFORE ComBat-seq:\n")
cat(sprintf("Mean within-batch correlation: %.4f\n", 
            mean(c(mean_within_batch1_before, mean_within_batch2_before))))
cat(sprintf("Mean between-batch correlation: %.4f\n", mean_between_batches_before))
cat(sprintf("Difference (within - between): %.4f\n", 
            mean(c(mean_within_batch1_before, mean_within_batch2_before)) - mean_between_batches_before))

# AFTER correction
cor_within_batch1_after <- cor(logCPM_combatseq[, batch1_samples])
mean_within_batch1_after <- mean(cor_within_batch1_after[upper.tri(cor_within_batch1_after)])

cor_within_batch2_after <- cor(logCPM_combatseq[, batch2_samples])
mean_within_batch2_after <- mean(cor_within_batch2_after[upper.tri(cor_within_batch2_after)])

cor_between_batches_after <- cor(logCPM_combatseq[, batch1_samples], 
                                 logCPM_combatseq[, batch2_samples])
mean_between_batches_after <- mean(cor_between_batches_after)

cat("\nAFTER ComBat-seq:\n")
cat(sprintf("Mean within-batch correlation: %.4f\n", 
            mean(c(mean_within_batch1_after, mean_within_batch2_after))))
cat(sprintf("Mean between-batch correlation: %.4f\n", mean_between_batches_after))
cat(sprintf("Difference (within - between): %.4f\n", 
            mean(c(mean_within_batch1_after, mean_within_batch2_after)) - mean_between_batches_after))

# Create comparison barplot with annotations and legend outside
tiff("results/figures/batch_correction/Correlation_comparison.tiff", 
     res = 300, units = "in", height = 6, width = 10)

# Adjust margins to make room for legend on the right
par(mar = c(5, 4, 4, 8), xpd = TRUE)

comparison_data <- matrix(
  c(mean(c(mean_within_batch1_before, mean_within_batch2_before)),
    mean_between_batches_before,
    mean(c(mean_within_batch1_after, mean_within_batch2_after)),
    mean_between_batches_after),
  nrow = 2,
  dimnames = list(c("Within Batch", "Between Batch"),
                  c("Before Correction", "After Correction"))
)

bp <- barplot(comparison_data,
              beside = TRUE,
              col = c("lightblue", "pink"),
              legend = FALSE,  # Turn off automatic legend
              ylab = "Mean Correlation",
              main = "Sample Correlations: Before vs After Batch Correction",
              ylim = c(0, 1.1))

# Add value annotations on bars
for (i in 1:ncol(comparison_data)) {
  for (j in 1:nrow(comparison_data)) {
    text(x = bp[j, i], 
         y = comparison_data[j, i] + 0.02,
         labels = sprintf("%.3f", comparison_data[j, i]),
         cex = 0.8, pos = 3)
  }
}

# Add legend outside the plot area
legend("topright", 
       inset = c(-0.15, 0),
       legend = rownames(comparison_data),
       fill = c("lightblue", "pink"),
       bty = "n")

dev.off()

## RLE (Relative Log Expression) Plots
cat("\n========== RLE ANALYSIS ==========\n")

# Function to calculate RLE
calculate_RLE <- function(data) {
  median_expr <- apply(data, 1, median)
  rle_data <- sweep(data, 1, median_expr, "-")
  return(rle_data)
}

# Before correction
rle_before <- calculate_RLE(logCounts_before)

# After correction  
rle_after <- calculate_RLE(logCPM_combatseq)

# Create RLE plots
tiff("results/figures/batch_correction/RLE_comparison.tiff", 
     res = 300, units = "in", height = 8, width = 14)
par(mfrow = c(1, 2))

# Before
boxplot(rle_before, 
        outline = FALSE,
        col = as.numeric(as.factor(d_filt$samples$Batch)),
        main = "RLE Plot - Before Batch Correction",
        ylab = "Relative Log Expression",
        xlab = "Samples",
        las = 2,
        names = targets$Label)
abline(h = 0, lty = 2, col = "darkgray", lwd = 2)
legend("topright", 
       legend = paste("Batch", unique_batches),
       fill = 1:length(unique_batches),
       bty = "n")

# After
boxplot(rle_after, 
        outline = FALSE,
        col = as.numeric(as.factor(d_filt_combatseq$samples$Batch)),
        main = "RLE Plot - After Batch Correction",
        ylab = "Relative Log Expression",
        xlab = "Samples",
        las = 2,
        names = targets$Label)
abline(h = 0, lty = 2, col = "darkgray", lwd = 2)

dev.off()

# Calculate RLE statistics
rle_stats_before <- data.frame(
  Sample = colnames(rle_before),
  Batch = d_filt$samples$Batch,
  Median = apply(rle_before, 2, median),
  IQR = apply(rle_before, 2, IQR)
)

rle_stats_after <- data.frame(
  Sample = colnames(rle_after),
  Batch = d_filt_combatseq$samples$Batch,
  Median = apply(rle_after, 2, median),
  IQR = apply(rle_after, 2, IQR)
)

cat("\nRLE Statistics Summary:\n")
cat("BEFORE - Median SD:", sd(rle_stats_before$Median), "\n")
cat("AFTER - Median SD:", sd(rle_stats_after$Median), "\n")
cat("Improvement:", round((1 - sd(rle_stats_after$Median)/sd(rle_stats_before$Median)) * 100, 1), "%\n")

save(batch_results, file = "results/tables/batch_correction_results.RData")

# Save key objects for next steps
save(d_filt, d_filt_combatseq, logCPM_combatseq, 
     file = "results/checkpoints/03_batch_corrected.RData")


## ===== 04_preprocessing.R =====
#-------------------------------------------------------#
#     04. Data Preprocessing                            #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading batch-corrected data...\n")
load_checkpoint("results/checkpoints/03_batch_corrected.RData")

# Note: Filtering and normalization were already done in script 03
# This script focuses on additional preprocessing and data preparation

## Rename for consistency
d_filt <- d_filt_combatseq

## Calculate various expression matrices
cat("\nCalculating expression matrices...\n")

# Raw counts
raw_counts <- d_filt$counts

# CPM values
cpm_values <- cpm(d_filt)

# LogCPM values (already calculated but ensure it's available)
if (!exists("logCPM_combatseq")) {
  logCPM_combatseq <- cpm(d_filt, log = TRUE, prior.count = 2)
}

# Rename for consistency
logCPM <- logCPM_combatseq

## Expression summaries by group
cat("\nCalculating group-wise expression summaries...\n")

# Calculate mean expression for each group
group_means <- sapply(unique(targets$group), function(grp) {
  samples <- which(targets$group == grp)
  rowMeans(logCPM[, samples])
})

# Calculate averages by timepoint across all lines
timepoint_means <- sapply(unique(targets$Timepoint), function(tp) {
  samples <- which(targets$Timepoint == tp)
  rowMeans(logCPM[, samples])
})

# Calculate averages by leaf type
leaftype_means <- sapply(unique(targets$Leaf_type), function(lt) {
  samples <- which(targets$Leaf_type == lt)
  rowMeans(logCPM[, samples])
})

## Identify highly variable genes
cat("\nIdentifying highly variable genes...\n")

# Calculate coefficient of variation for each gene
gene_cv <- apply(logCPM, 1, function(x) sd(x) / mean(x))
highly_variable <- names(sort(gene_cv, decreasing = TRUE)[1:1000])

## Write out key expression matrices
cat("\nSaving expression matrices...\n")

write.csv(logCPM, file = "results/tables/logCPM_all_samples.csv")
write.csv(group_means, file = "results/tables/logCPM_group_means.csv")
write.csv(cpm_values, file = "results/tables/CPM_values.csv")

# Save preprocessing summary
preprocessing_summary <- list(
  total_genes = nrow(d_filt),
  total_samples = ncol(d_filt),
  library_sizes = d_filt$samples$lib.size,
  norm_factors = d_filt$samples$norm.factors,
  highly_variable_genes = highly_variable,
  expression_range = range(logCPM)
)

save(preprocessing_summary, file = "results/tables/preprocessing_summary.RData")

# Save checkpoint with all processed data
save(d_filt, logCPM, logCPM_combatseq, cpm_values, raw_counts,
     group_means, timepoint_means, leaftype_means,
     file = "results/checkpoints/04_preprocessed_data.RData")

cat("\n=== Data preprocessing complete ===\n")
cat("Total genes for analysis:", nrow(d_filt), "\n") #38045
cat("Total samples:", ncol(d_filt), "\n") #60 



## Housekeeping Gene Analysis with Normalized Data
cat("\nAnalyzing housekeeping genes with normalized data...\n")

# Load housekeeping genes if not already loaded
if (!exists("hk.genes")) {
  library(readxl)
  hk.genes <- read_excel("Housekeeping_genes_Machado_et_al_2019.xlsx", sheet = 1)
  hk.genes$Gene <- as.character(hk.genes$Gene)
}

# Get housekeeping genes present in our filtered, normalized data
hk_genes_present <- hk.genes$Gene[hk.genes$Gene %in% rownames(logCPM)]
cat("Housekeeping genes in filtered data:", length(hk_genes_present), "\n")

if (length(hk_genes_present) > 0) {
  # Extract expression data
  hk_genes_logCPM <- logCPM[hk_genes_present, ]
  
  # Order samples
  i_names <- order(targets$Timepoint, targets$Leaf_type)
  
  # Create heatmap
  library(gplots)
  color_scale <- colorpanel(100, "blue", "white", "red")
  
  tiff("results/figures/Housekeeping_genes_normalized_heatmap.tiff", 
       res = 300, unit = "in", height = 8, width = 10)
  
  heatmap.2(hk_genes_logCPM[, i_names], 
            col = color_scale, 
            density.info = "none", 
            trace = "none", 
            labRow = ifelse(length(hk_genes_present) > 50, "", rownames(hk_genes_logCPM)),
            colsep = c(12, 24, 36, 48),
            margins = c(7, 2), 
            key.xlab = "logCPM (normalized)",
            Colv = FALSE, 
            dendrogram = "row",
            main = paste("Housekeeping Genes (Normalized) -", length(hk_genes_present), "genes"))
  
  dev.off()
  
  cat("Housekeeping gene heatmap created\n")
}


## ===== 05_statistical_design.R =====
#-------------------------------------------------------#
#     05. Statistical Design Setup                      #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading preprocessed data...\n")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")

# Create output directory
create_dir("results/tables/design")

## Design Matrix
cat("\nCreating design matrix...\n")

# Create design matrix with one coefficient for each group
design <- model.matrix(~ 0 + d_filt$samples$group)

# Clean up column names
colnames(design) <- levels(d_filt$samples$group)
rownames(design) <- d_filt$samples$Label

cat("Design matrix dimensions:", dim(design), "\n")
cat("Groups in design:", colnames(design), "\n")

# Save design matrix
write.csv(design, file = "results/tables/design/design_matrix.csv")

## Contrast Matrix
cat("\nCreating contrast matrix...\n")

# Define all contrasts
cont_matrix <- makeContrasts(
  # PI532462A (Broad) comparisons
  B_462A_1vs0 = PI532462A_TP1 - PI532462A_TP0, 
  B_462A_2vs0 = PI532462A_TP2 - PI532462A_TP0,
  B_462A_3vs0 = PI532462A_TP3 - PI532462A_TP0,
  B_462A_4vs0 = PI532462A_TP4 - PI532462A_TP0,
  
  # PI612713B (Narrow) comparisons
  N_713B_1vs0 = PI612713B_TP1 - PI612713B_TP0, 
  N_713B_2vs0 = PI612713B_TP2 - PI612713B_TP0,
  N_713B_3vs0 = PI612713B_TP3 - PI612713B_TP0,
  N_713B_4vs0 = PI612713B_TP4 - PI612713B_TP0,
  
  # PI547745 (Narrow) comparisons
  N_745_1vs0 = PI547745_TP1 - PI547745_TP0, 
  N_745_2vs0 = PI547745_TP2 - PI547745_TP0,
  N_745_3vs0 = PI547745_TP3 - PI547745_TP0,
  N_745_4vs0 = PI547745_TP4 - PI547745_TP0,
  
  # LD112170 (Broad) comparisons
  B_LD11_1vs0 = LD112170_TP1 - LD112170_TP0, 
  B_LD11_2vs0 = LD112170_TP2 - LD112170_TP0,
  B_LD11_3vs0 = LD112170_TP3 - LD112170_TP0,
  B_LD11_4vs0 = LD112170_TP4 - LD112170_TP0,
  
  # Between line comparisons at TP0
  B_LD11.vs.N_713B_TP0 = PI612713B_TP0 - LD112170_TP0,
  B_LD11.vs.N_745_TP0 = PI547745_TP0 - LD112170_TP0,
  B_462A.vs.N_713B_TP0 = PI612713B_TP0 - PI532462A_TP0,
  B_462A.vs.N_745_TP0 = PI547745_TP0 - PI532462A_TP0,
  B_LD11.vs.B_462A_TP0 = PI532462A_TP0 - LD112170_TP0,
  N_713B.vs.N_745_TP0 = PI547745_TP0 - PI612713B_TP0,
  
  # Additional timepoint comparisons
  B_462A_2vs1 = PI532462A_TP2 - PI532462A_TP1,
  B_462A_3vs1 = PI532462A_TP3 - PI532462A_TP1,
  B_462A_4vs1 = PI532462A_TP4 - PI532462A_TP1,
  B_462A_3vs2 = PI532462A_TP3 - PI532462A_TP2,
  B_462A_4vs2 = PI532462A_TP4 - PI532462A_TP2,
  B_462A_4vs3 = PI532462A_TP4 - PI532462A_TP3,
  
  N_713B_2vs1 = PI612713B_TP2 - PI612713B_TP1, 
  N_713B_3vs1 = PI612713B_TP3 - PI612713B_TP1,
  N_713B_4vs1 = PI612713B_TP4 - PI612713B_TP1,
  N_713B_3vs2 = PI612713B_TP3 - PI612713B_TP2, 
  N_713B_4vs2 = PI612713B_TP4 - PI612713B_TP2,
  N_713B_4vs3 = PI612713B_TP4 - PI612713B_TP3,
  
  N_745_2vs1 = PI547745_TP2 - PI547745_TP1, 
  N_745_3vs1 = PI547745_TP3 - PI547745_TP1,
  N_745_4vs1 = PI547745_TP4 - PI547745_TP1,
  N_745_3vs2 = PI547745_TP3 - PI547745_TP2,
  N_745_4vs2 = PI547745_TP4 - PI547745_TP2,
  N_745_4vs3 = PI547745_TP4 - PI547745_TP3,
  
  B_LD11_2vs1 = LD112170_TP2 - LD112170_TP1, 
  B_LD11_3vs1 = LD112170_TP3 - LD112170_TP1,
  B_LD11_4vs1 = LD112170_TP4 - LD112170_TP1,
  B_LD11_3vs2 = LD112170_TP3 - LD112170_TP2,
  B_LD11_4vs2 = LD112170_TP4 - LD112170_TP2,
  B_LD11_4vs3 = LD112170_TP4 - LD112170_TP3,
  
  levels = design
)

cat("Number of contrasts:", ncol(cont_matrix), "\n")

# Save contrast matrix
write.csv(cont_matrix, file = "results/tables/design/contrast_matrix.csv")

## Estimate Dispersions
cat("\nEstimating dispersions...\n")

d_filt <- estimateDisp(d_filt, design, robust = TRUE)

# Check dispersion values
cat("Common dispersion:", d_filt$common.dispersion, "\n") #0.04131481
cat("Biological coefficient of variation (BCV):", sqrt(d_filt$common.dispersion), "\n") #0.2032604 

## Fit Statistical Model
cat("\nFitting statistical model...\n")

fit_edgeR <- glmQLFit(d_filt, design, robust = TRUE)

# Plot quasi-likelihood dispersions
pdf("results/figures/QL_dispersion_plot.pdf")
plotQLDisp(fit_edgeR)
dev.off()

# Save statistical models
save(design, cont_matrix, d_filt, fit_edgeR,
     file = "results/checkpoints/05_statistical_models.RData")

cat("\n=== Statistical design complete ===\n")
cat("Ready for differential expression testing\n")

#-------------------------------------------------------#
#     06. Differential Expression Analysis              #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading statistical models...\n")
load_checkpoint("results/checkpoints/05_statistical_models.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")  # Ensure logCPM is available

# Create output directory
create_dir("results/tables/DE")

## Run Differential Expression Tests
cat("\nRunning differential expression tests...\n")
cat("Total contrasts to test:", ncol(cont_matrix), "\n\n")

# Initialize list to store results
de_results <- list()
de_results_treat <- list()

# Define contrast names
contrast_names <- colnames(cont_matrix)

# Run tests for each contrast
for (i in 1:ncol(cont_matrix)) {
  contrast_name <- contrast_names[i]
  cat("Testing contrast", i, "of", ncol(cont_matrix), ":", contrast_name, "\n")
  
  # Standard test
  de_results[[contrast_name]] <- glmQLFTest(fit_edgeR, contrast = cont_matrix[, i])
  
  # Test with fold-change threshold
  de_results_treat[[contrast_name]] <- glmTreat(fit_edgeR, 
                                                contrast = cont_matrix[, i], 
                                                lfc = log2(1.2))
}

## Extract and Format Results
cat("\nExtracting detailed results...\n")

# Function to extract and format results
extract_results <- function(de_obj, contrast_name) {
  results <- topTags(de_obj, n = Inf, sort.by = "none")$table
  results$FC <- 2^abs(results$logFC) * sign(results$logFC)
  results$contrast <- contrast_name
  # Ensure GeneID column exists
  if (!"GeneID" %in% colnames(results)) {
    results$GeneID <- rownames(results)
  }
  return(results)
}

# Extract all results
all_results <- list()
all_results_treat <- list()

for (contrast_name in contrast_names) {
  all_results[[contrast_name]] <- extract_results(de_results[[contrast_name]], 
                                                  contrast_name)
  all_results_treat[[contrast_name]] <- extract_results(de_results_treat[[contrast_name]], 
                                                        contrast_name)
}

## Create Decision Matrix
cat("\nCreating decision matrices...\n")

# Standard tests
edgeR_coded <- do.call(cbind, lapply(de_results, function(x) {
  decideTests(x)[, 1]
}))

# Tests with FC threshold
edgeR_tr_coded <- do.call(cbind, lapply(de_results_treat, function(x) {
  decideTests(x)[, 1]
}))

# IMPORTANT: Set column names before converting to TestResults
colnames(edgeR_coded) <- names(de_results)
colnames(edgeR_tr_coded) <- names(de_results_treat)

# Convert to TestResults object
edgeR_coded <- new("TestResults", edgeR_coded)
edgeR_tr_coded <- new("TestResults", edgeR_tr_coded)

# Set attributes
attr(edgeR_coded, "labels") <- c("Down", "NotSig", "Up")
attr(edgeR_tr_coded, "labels") <- c("Down", "NotSig", "Up")

## Summary Statistics
cat("\nDifferential expression summary (FC > 1.2):\n")
de_summary <- summary(edgeR_tr_coded)
print(de_summary)

# Save summary
write.csv(de_summary, file = "results/tables/DE/DE_summary_FC1.2.csv")

## Create Comprehensive DE Results Matrix
cat("\nCreating comprehensive DE results matrix...\n")

# Get all genes
all_genes <- rownames(d_filt)

# Initialize matrices for different statistics
logFC_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
pvalue_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
fdr_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
logCPM_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))

# Set row and column names
rownames(logFC_matrix) <- all_genes
rownames(pvalue_matrix) <- all_genes
rownames(fdr_matrix) <- all_genes
rownames(logCPM_matrix) <- all_genes

colnames(logFC_matrix) <- paste0(contrast_names, "_logFC")
colnames(pvalue_matrix) <- paste0(contrast_names, "_PValue")
colnames(fdr_matrix) <- paste0(contrast_names, "_FDR")
colnames(logCPM_matrix) <- paste0(contrast_names, "_logCPM")

# Fill matrices with results
for (i in 1:length(contrast_names)) {
  contrast_name <- contrast_names[i]
  results <- all_results_treat[[contrast_name]]
  
  # Match genes - handle both cases where results might have GeneID column or use rownames
  if ("GeneID" %in% colnames(results)) {
    result_genes <- results$GeneID
  } else {
    result_genes <- rownames(results)
  }
  
  gene_idx <- match(result_genes, all_genes)
  
  # Fill in values
  logFC_matrix[gene_idx, i] <- results$logFC
  pvalue_matrix[gene_idx, i] <- results$PValue
  fdr_matrix[gene_idx, i] <- results$FDR
  logCPM_matrix[gene_idx, i] <- results$logCPM
}

# Create comprehensive data frame
comprehensive_de_results <- data.frame(
  GeneID = all_genes,
  stringsAsFactors = FALSE
)

# Add all statistics
comprehensive_de_results <- cbind(
  comprehensive_de_results,
  logFC_matrix,
  fdr_matrix,
  pvalue_matrix,
  logCPM_matrix
)

# Add summary statistics
comprehensive_de_results$Max_AbsLogFC <- apply(logFC_matrix, 1, function(x) {
  if(all(is.na(x))) return(NA)
  max(abs(x), na.rm = TRUE)
})

comprehensive_de_results$Min_FDR <- apply(fdr_matrix, 1, function(x) {
  if(all(is.na(x))) return(NA)
  min(x, na.rm = TRUE)
})

comprehensive_de_results$N_Sig_Contrasts <- rowSums(
  fdr_matrix < 0.05 & abs(logFC_matrix) > log2(1.2), 
  na.rm = TRUE
)

# Add average expression
if (exists("logCPM")) {
  comprehensive_de_results$AveExpr <- rowMeans(logCPM[all_genes, ])
}

# Add decision matrix
decision_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
rownames(decision_matrix) <- all_genes
colnames(decision_matrix) <- paste0(contrast_names, "_Decision")

for (i in 1:length(contrast_names)) {
  fdr_col <- fdr_matrix[, i]
  logfc_col <- logFC_matrix[, i]
  
  decision_matrix[, i] <- ifelse(
    is.na(fdr_col) | is.na(logfc_col), NA,
    ifelse(fdr_col < 0.05 & logfc_col > log2(1.2), "Up",
           ifelse(fdr_col < 0.05 & logfc_col < -log2(1.2), "Down", "None"))
  )
}

comprehensive_de_results <- cbind(comprehensive_de_results, decision_matrix)

## Add Gene Annotation
cat("\nAdding gene annotation...\n")

# Load annotation file
annotation_file <- "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt"
if (file.exists(annotation_file)) {
  annotation <- read.delim(annotation_file, stringsAsFactors = FALSE)
  
  # Check available columns
  cat("Available annotation columns:\n")
  print(colnames(annotation))
  
  # Use all columns except pacId (which is just an internal ID)
  exclude_cols <- c("pacId")
  
  if ("locusName" %in% colnames(annotation)) {
    # Get all columns except excluded ones
    cols_to_use <- setdiff(colnames(annotation), exclude_cols)
    gene_annotation <- annotation[, cols_to_use]
    
    # Remove duplicate genes (keep first occurrence)
    gene_annotation <- gene_annotation[!duplicated(gene_annotation$locusName), ]
    
    # Match annotation to our genes
    annotation_idx <- match(comprehensive_de_results$GeneID, gene_annotation$locusName)
    
    # Add annotation columns (excluding locusName to avoid duplication)
    ann_cols_to_add <- setdiff(colnames(gene_annotation), "locusName")
    
    # Insert annotation after GeneID
    comprehensive_de_results_annotated <- cbind(
      comprehensive_de_results[, "GeneID", drop = FALSE],
      gene_annotation[annotation_idx, ann_cols_to_add, drop = FALSE],
      comprehensive_de_results[, -1, drop = FALSE]  # All other columns except GeneID
    )
    
    # Use annotated version
    comprehensive_de_results <- comprehensive_de_results_annotated
    
    # Report annotation coverage
    n_annotated <- sum(!is.na(annotation_idx))
    cat("Genes with annotation:", n_annotated, "out of", nrow(comprehensive_de_results), 
        "(", round(n_annotated/nrow(comprehensive_de_results)*100, 1), "%)\n")
  } else {
    cat("Warning: locusName column not found in annotation file\n")
  }
} else {
  cat("Warning: Annotation file not found:", annotation_file, "\n")
  cat("Proceeding without annotation\n")
}

# Sort by number of significant contrasts and minimum FDR
comprehensive_de_results <- comprehensive_de_results[
  order(comprehensive_de_results$N_Sig_Contrasts, 
        comprehensive_de_results$Min_FDR, 
        decreasing = c(TRUE, FALSE)), 
]

## Write Results Files
cat("\nWriting result files...\n")

# Write individual contrast results
for (i in seq_along(contrast_names)) {
  contrast_name <- contrast_names[i]
  
  # Add decision codes to results - use index instead of name
  all_results_treat[[contrast_name]]$DE_code <- edgeR_tr_coded[, i]
  
  # Write to file
  filename <- paste0("results/tables/DE/", contrast_name, "_results.csv")
  write.csv(all_results_treat[[contrast_name]], file = filename, row.names = TRUE)
}

# Write comprehensive matrix - all genes
write.csv(comprehensive_de_results, 
          file = "results/tables/DE/comprehensive_DE_results_all_contrasts.csv", 
          row.names = FALSE)

# Write significant genes only
significant_genes <- comprehensive_de_results[comprehensive_de_results$N_Sig_Contrasts > 0, ]
write.csv(significant_genes, 
          file = "results/tables/DE/comprehensive_DE_results_significant_only.csv", 
          row.names = FALSE)

cat("Comprehensive DE matrix created with", nrow(comprehensive_de_results), "genes\n")
cat("Significant genes in at least one contrast:", nrow(significant_genes), "\n")

# Create simplified version with key contrasts
key_contrasts <- c("B_462A_4vs0", "N_713B_4vs0", "N_745_4vs0", "B_LD11_4vs0")
key_cols <- c("GeneID")

# Add annotation columns if they exist
if ("Best-hit-arabi-name" %in% colnames(comprehensive_de_results)) {
  key_cols <- c(key_cols, "Best-hit-arabi-name", "Best-hit-arabi-defline")
}

# Add expression and significance summary
key_cols <- c(key_cols, "AveExpr", "N_Sig_Contrasts")

# Add key contrast results
for (contrast in key_contrasts) {
  key_cols <- c(key_cols, 
                paste0(contrast, "_logFC"),
                paste0(contrast, "_FDR"),
                paste0(contrast, "_Decision"))
}

# Extract columns that exist
key_cols <- key_cols[key_cols %in% colnames(comprehensive_de_results)]
simplified_de_results <- comprehensive_de_results[, key_cols]

write.csv(simplified_de_results, 
          file = "results/tables/DE/simplified_DE_results_key_contrasts.csv", 
          row.names = FALSE)

# Create top genes summary
top_genes <- head(significant_genes, 100)
write.csv(top_genes, 
          file = "results/tables/DE/top100_DE_genes_annotated.csv", 
          row.names = FALSE)

# ## OPTIONAL: Create Original Master Gene Table (only if needed for legacy scripts)
# ## Skip this if you're running out of memory
# create_master_table <- FALSE  # Set to TRUE only if needed
# 
# if (create_master_table) {
#   cat("\nCreating master gene table (legacy format)...\n")
#   
#   # Only create for significant genes to save memory
#   sig_genes <- rownames(edgeR_tr_coded)[rowSums(abs(edgeR_tr_coded) > 0) > 0]
#   cat("Creating master table for", length(sig_genes), "significant genes only\n")
#   
#   master_gene_table <- cbind(
#     d_filt$genes[sig_genes, , drop=FALSE],
#     logFC_matrix[sig_genes, ],
#     fdr_matrix[sig_genes, ],
#     edgeR_tr_coded[sig_genes, ]
#   )
#   
#   master_gene_table$AveExpr <- rowMeans(logCPM[sig_genes, ])
#   
#   write.csv(master_gene_table, file = "results/tables/DE/master_gene_table_sig_only.csv", 
#             row.names = TRUE)
# } else {
#   cat("\nSkipping master gene table creation (set create_master_table = TRUE if needed)\n")
# }

# Save DE results
save(de_results, de_results_treat, all_results, all_results_treat,
     edgeR_coded, edgeR_tr_coded, de_summary,
     comprehensive_de_results, significant_genes,
     file = "results/checkpoints/06_de_results.RData")

cat("\n=== Differential expression analysis complete ===\n")
cat("Total genes tested:", nrow(edgeR_tr_coded), "\n")
cat("Files written to: results/tables/DE/\n")
cat("\nKey output files:\n")
cat("- comprehensive_DE_results_all_contrasts.csv: Full results with annotation\n")
cat("- comprehensive_DE_results_significant_only.csv: Filtered to significant genes\n")
cat("- simplified_DE_results_key_contrasts.csv: Key contrasts only\n")
cat("- top100_DE_genes_annotated.csv: Top differentially expressed genes\n")


#-------------------------------------------------------#
#     07. Differential Expression Visualization         #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoints
cat("Loading DE results...\n")
load_checkpoint("results/checkpoints/06_de_results.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")
load_checkpoint("results/checkpoints/05_statistical_models.RData")  # Contains d_filt and cont_matrix

# Ensure targets exists
if (!exists("targets")) {
  targets <- d_filt$samples
}

# Load additional libraries
library(Glimma)
library(reshape2)
library (limma)
library (ggplot2)
library(UpSetR)  # For handling >5 set comparisons

# Create output directories
create_dir("results/figures/DE/volcano_plots")
create_dir("results/figures/DE/MA_plots")
create_dir("results/figures/DE/venn_diagrams")
create_dir("results/figures/DE/barplots")
create_dir("results/figures/DE/glimma")

## 1. DEG Bar Charts
cat("\nCreating DEG bar charts...\n")

# Function to create DEG barplot
create_deg_barplot <- function(line_name, contrasts, filename) {
  # Extract relevant data
  deg_data <- t(summary(edgeR_tr_coded)[c(1,3), contrasts])
  
  # Melt for ggplot
  deg_melted <- melt(deg_data, varnames = c("Comparison", "Direction"),
                     value.name = "gene.no")
  
  # Make down values negative
  deg_melted$gene.no <- ifelse(deg_melted$Direction == "Down", 
                               -abs(deg_melted$gene.no), 
                               deg_melted$gene.no)
  
  # Create plot
  p <- ggplot(deg_melted, aes(x = Comparison, y = gene.no, fill = Direction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
    labs(y = "Number of genes", x = "Comparison", title = line_name) +
    geom_text(aes(label = abs(gene.no), 
                  y = ifelse(gene.no > 0, gene.no + 100, gene.no - 100)),
              color = "black", size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# Create barplots for each line
create_deg_barplot("PI532462A (Broad)", 1:4, 
                   "results/figures/DE/barplots/DEG_barplot_462A.tiff")
create_deg_barplot("PI612713B (Narrow)", 5:8, 
                   "results/figures/DE/barplots/DEG_barplot_713B.tiff")
create_deg_barplot("PI547745 (Narrow)", 9:12, 
                   "results/figures/DE/barplots/DEG_barplot_745.tiff")
create_deg_barplot("LD112170 (Broad)", 13:16, 
                   "results/figures/DE/barplots/DEG_barplot_LD11.tiff")

## 2. MA/MD Plots
cat("\nCreating MA plots...\n")

# Function to create MA plot
create_ma_plot <- function(de_obj, de_coded, contrast_idx, contrast_name) {
  filename <- paste0("results/figures/DE/MA_plots/MA_", contrast_name, ".tiff")
  
  tiff(filename, res = 300, units = "in", height = 6, width = 6)
  plotMD(de_obj, status = de_coded[, contrast_idx], 
         legend = "topright", hl.cex = 0.7, 
         main = gsub("_", " ", contrast_name))
  dev.off()
}

# Create MA plots for key contrasts
key_contrasts <- c("B_462A_1vs0", "B_462A_4vs0", "N_713B_1vs0", "N_713B_4vs0",
                   "N_745_1vs0", "N_745_4vs0", "B_LD11_1vs0", "B_LD11_4vs0")

for (i in seq_along(key_contrasts)) {
  contrast_name <- key_contrasts[i]
  contrast_idx <- which(colnames(cont_matrix) == contrast_name)
  if (length(contrast_idx) > 0) {
    create_ma_plot(de_results_treat[[contrast_name]], edgeR_tr_coded, 
                   contrast_idx, contrast_name)
  }
}

## 3. Volcano Plots
cat("\nCreating volcano plots...\n")

# Function to create volcano plot
create_volcano_plot <- function(results, contrast_name) {
  # Prepare data
  volcano_data <- data.frame(
    logFC = results$logFC,
    negLogP = -log10(results$PValue),
    FDR = results$FDR,
    gene = rownames(results)
  )
  
  # Define significance
  volcano_data$sig <- "Not Sig"
  volcano_data$sig[volcano_data$FDR < 0.05 & volcano_data$logFC > log2(1.2)] <- "Up"
  volcano_data$sig[volcano_data$FDR < 0.05 & volcano_data$logFC < -log2(1.2)] <- "Down"
  
  # Create plot
  p <- ggplot(volcano_data, aes(x = logFC, y = negLogP, color = sig)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "gray")) +
    geom_vline(xintercept = c(-log2(1.2), log2(1.2)), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(x = "log2 Fold Change", y = "-log10(p-value)", 
         title = gsub("_", " ", contrast_name)) +
    theme_minimal()
  
  filename <- paste0("results/figures/DE/volcano_plots/volcano_", contrast_name, ".pdf")
  ggsave(filename, plot = p, width = 8, height = 6)
}

# Create volcano plots for key contrasts
for (contrast_name in key_contrasts) {
  create_volcano_plot(all_results_treat[[contrast_name]], contrast_name)
}

## 4. Venn Diagrams
cat("\nCreating Venn diagrams...\n")

# 4-way Venn for each line
create_line_venn <- function(line_name, contrast_indices, filename) {
  tiff(filename, res = 300, units = "in", height = 8, width = 8)
  vennDiagram(edgeR_tr_coded[, contrast_indices], 
              include = "both",
              main = paste(line_name, "- DEGs across timepoints"))
  dev.off()
}

create_line_venn("PI532462A", 1:4, 
                 "results/figures/DE/venn_diagrams/venn_462A_DEGs.tiff")
create_line_venn("PI612713B", 5:8, 
                 "results/figures/DE/venn_diagrams/venn_713B_DEGs.tiff")
create_line_venn("PI547745", 9:12, 
                 "results/figures/DE/venn_diagrams/venn_745_DEGs.tiff")
create_line_venn("LD112170", 13:16, 
                 "results/figures/DE/venn_diagrams/venn_LD11_DEGs.tiff")

# Venn for leaf type comparisons at TP0
# Split into two 3-way Venns since vennDiagram can't handle >5 sets
cat("Creating Venn diagrams for TP0 comparisons (split due to >5 sets)...\n")

tiff("results/figures/DE/venn_diagrams/venn_TP0_comparisons_part1.tiff", 
     res = 300, units = "in", height = 8, width = 8)
vennDiagram(edgeR_tr_coded[, 17:19], 
            include = "both",
            main = "Leaf type comparisons at TP0 (Part 1)")
dev.off()

tiff("results/figures/DE/venn_diagrams/venn_TP0_comparisons_part2.tiff", 
     res = 300, units = "in", height = 8, width = 8)
vennDiagram(edgeR_tr_coded[, 20:22], 
            include = "both",
            main = "Leaf type comparisons at TP0 (Part 2)")
dev.off()

# Alternative: Create UpSet plot for all 6 comparisons
if (requireNamespace("UpSetR", quietly = TRUE)) {
  # Convert to binary matrix for UpSetR
  upset_data <- as.data.frame(edgeR_tr_coded[, 17:22] != 0)
  colnames(upset_data) <- colnames(cont_matrix)[17:22]
  
  # Remove rows with no DEGs
  upset_data <- upset_data[rowSums(upset_data) > 0, ]
  
  pdf("results/figures/DE/venn_diagrams/upset_TP0_comparisons.pdf", width = 12, height = 8)
  print(upset(upset_data, 
              sets = colnames(upset_data),
              order.by = "freq",
              main.bar.color = "darkblue",
              sets.bar.color = "darkred",
              text.scale = 1.5,
              set_size.show = TRUE))
  dev.off()
}

## 5. Interactive Glimma Plots
cat("\nCreating interactive Glimma plots...\n")

# Function to create Glimma MD plot
create_glimma_plot <- function(de_obj, contrast_idx, contrast_name) {
  glMDPlot(de_obj, 
           counts = logCPM, 
           status = edgeR_tr_coded[, contrast_idx], 
           groups = d_filt$samples$group, 
           samples = d_filt$samples$Label, 
           sample.cols = d_filt$samples$col,
           folder = "results/figures/DE/glimma", 
           html = paste0("MDplot_", contrast_name), 
           main = gsub("_", " ", contrast_name),
           launch = FALSE)
}

# Create Glimma plots for key contrasts
for (contrast_name in key_contrasts) {
  contrast_idx <- which(colnames(cont_matrix) == contrast_name)
  if (length(contrast_idx) > 0) {
    create_glimma_plot(de_results_treat[[contrast_name]], contrast_idx, contrast_name)
  }
}

## 6. Summary Heatmap of Top DEGs
cat("\nCreating summary heatmap...\n")

# Get top DEGs from each contrast
top_degs <- unique(unlist(lapply(all_results_treat, function(x) {
  x <- x[order(x$FDR), ]
  # Make sure we don't exceed the number of rows
  n_top <- min(50, nrow(x))
  if (n_top > 0) {
    rownames(x)[1:n_top]
  } else {
    NULL
  }
})))

# Remove any NA values
top_degs <- top_degs[!is.na(top_degs)]

if (length(top_degs) > 0) {
  # Create heatmap
  library(gplots)
  color_scale <- colorpanel(100, "blue", "white", "red")
  
  # Scale expression data
  heatmap_data <- logCPM[top_degs, ]
  heatmap_scaled <- t(scale(t(heatmap_data)))
  
  # Get sample colors
  sample_colors <- as.character(d_filt$samples$col)
  
  pdf("results/figures/DE/top_DEGs_heatmap.pdf", width = 10, height = 12)
  heatmap.2(heatmap_scaled, 
            col = color_scale, 
            trace = "none", 
            density.info = "none",
            margins = c(10, 8),
            ColSideColors = sample_colors,
            labRow = FALSE,
            main = paste("Top DEGs per contrast (n =", length(top_degs), ")"))
  dev.off()
  
  # Also create a smaller heatmap with gene labels for top 100 genes
  if (length(top_degs) > 100) {
    top_100_degs <- top_degs[1:100]
    heatmap_data_100 <- logCPM[top_100_degs, ]
    heatmap_scaled_100 <- t(scale(t(heatmap_data_100)))
    
    pdf("results/figures/DE/top100_DEGs_heatmap_labeled.pdf", width = 12, height = 14)
    heatmap.2(heatmap_scaled_100, 
              col = color_scale, 
              trace = "none", 
              density.info = "none",
              margins = c(10, 10),
              ColSideColors = sample_colors,
              cexRow = 0.5,
              main = "Top 100 DEGs with labels")
    dev.off()
  }
} else {
  cat("Warning: No top DEGs found for heatmap creation\n")
}

# Save visualization checkpoint
save_checkpoint("results/checkpoints/07_de_visualizations_complete.RData")

cat("\n=== DE visualization complete ===\n")
cat("Figures saved to: results/figures/DE/\n")
cat("\nKey outputs:\n")
cat("- DEG barplots for each line\n")
cat("- MA plots for key contrasts\n")
cat("- Volcano plots for key contrasts\n")
cat("- Venn diagrams for timepoint comparisons\n")
cat("- UpSet plot for TP0 comparisons (if >5 sets)\n")
cat("- Interactive Glimma plots\n")
cat("- Heatmaps of top DEGs\n")



## ===== run_all.R =====
#-------------------------------------------------------#
#     Master Script - Run Complete Pipeline             #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# This script runs the complete analysis pipeline
# Each step can also be run independently

# Set up environment
cat("Starting Soybean RNA-Seq Analysis Pipeline\n")
cat("=======================================\n\n")

# Record start time
start_time <- Sys.time()

# Source utility functions
source("00_functions.R")

# Run each analysis step
tryCatch({
  
  # Step 1: Setup and Import
  cat("\n>>> STEP 1: Data Import <<<\n")
  source("01_setup_and_import.R")
  
  # Step 2: Quality Control
  cat("\n>>> STEP 2: Quality Control <<<\n")
  source("02_quality_control.R")
  
  # Step 3: Batch Correction
  cat("\n>>> STEP 3: Batch Correction <<<\n")
  source("03_batch_correction.R")
  
  # Step 4: Data Preprocessing
  cat("\n>>> STEP 4: Data Preprocessing <<<\n")
  source("04_preprocessing.R")
  
  # Step 5: Statistical Design
  cat("\n>>> STEP 5: Statistical Design <<<\n")
  source("05_statistical_design.R")
  
  # Step 6: Differential Expression
  cat("\n>>> STEP 6: Differential Expression Analysis <<<\n")
  source("06_differential_expression.R")
  
  # Step 7: DE Visualization
  cat("\n>>> STEP 7: DE Visualization <<<\n")
  source("07_de_visualization.R")
  
  # Additional analyses can be added here:
  # source("08_wgcna_prep.R")
  # source("09_wgcna_analysis.R")
  # source("10_wgcna_visualization.R")
  # source("11_gene_sets.R")
  # source("12_integration.R")
  # source("13_JAG1_analysis.R")
  # source("14_reporting.R")
  
  cat("\n\n=== PIPELINE COMPLETED SUCCESSFULLY ===\n")
  
}, error = function(e) {
  cat("\n\n!!! ERROR IN PIPELINE !!!\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("Please check the error and run the failed step individually\n")
})

# Record end time
end_time <- Sys.time()
cat("\nTotal runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Session info
cat("\n=== Session Information ===\n")
sessionInfo()






## ===== 08_wgcna_prep.R =====
#-------------------------------------------------------#
#     08. WGCNA Preparation                             #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading DE results for WGCNA preparation...\n")
load_checkpoint("results/checkpoints/06_de_results.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")


# Install BiocManager if you haven't already
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install GO.db
BiocManager::install("GO.db")
# Install impute package
BiocManager::install("impute")

library(edgeR)

# Load WGCNA library
library(WGCNA)
enableWGCNAThreads()

# Create output directory
create_dir("results/WGCNA")
create_dir("results/figures/WGCNA")

## ANOVA for Gene Selection
cat("\nPerforming ANOVA for gene selection...\n")

# Create design matrix with intercept for ANOVA
d_filt$samples$newGroup <- relevel(d_filt$samples$group, ref = "PI532462A_TP0")
design_anova <- model.matrix(~ d_filt$samples$newGroup)

# Clean column names
colnames(design_anova)[1] <- "Intercept"
colnames(design_anova)[-1] <- levels(d_filt$samples$newGroup)[-1]

# Re-estimate dispersions for new design
d_anova <- estimateDisp(d_filt, design_anova, robust = TRUE)

# Fit model
fit_anova <- glmQLFit(d_anova, design_anova, robust = TRUE)

# Test all groups simultaneously
eR_anova <- glmQLFTest(fit_anova, coef = 2:ncol(design_anova))

# Get detailed results
eR_anova_detailed <- topTags(eR_anova, n = Inf, sort.by = "none")$table

# Examine FDR distribution
pdf("results/figures/WGCNA/FDR_distribution_ANOVA.pdf")
hist(eR_anova_detailed$FDR, breaks = 100, 
     main = "FDR Distribution - ANOVA",
     xlab = "FDR", col = "lightblue")
abline(v = c(0.01, 0.05, 0.1), col = c("red", "orange", "green"), lty = 2)
dev.off()

# Gene selection for WGCNA
fdr_threshold <- 0.01
genes_for_wgcna <- rownames(eR_anova_detailed)[eR_anova_detailed$FDR < fdr_threshold]

cat("\nGene selection summary:\n")
cat("Total genes tested:", nrow(eR_anova_detailed), "\n")
cat("Genes passing FDR <", fdr_threshold, ":", length(genes_for_wgcna), "\n")
cat("Proportion selected:", round(length(genes_for_wgcna)/nrow(eR_anova_detailed), 3), "\n")

## Prepare Expression Data for WGCNA
cat("\nPreparing expression data for WGCNA...\n")

# Subset and transpose data (WGCNA expects samples in rows)
datExpr <- t(logCPM[genes_for_wgcna, ])

# Check for missing values
if (sum(is.na(datExpr)) > 0) {
  cat("Warning: Missing values detected. Removing genes with NAs...\n")
  good_genes <- colSums(is.na(datExpr)) == 0
  datExpr <- datExpr[, good_genes]
}

# Check for zero variance genes
gene_vars <- apply(datExpr, 2, var)
if (any(gene_vars == 0)) {
  cat("Warning: Zero variance genes detected. Removing...\n")
  datExpr <- datExpr[, gene_vars > 0]
}

cat("Final dimensions for WGCNA:\n")
cat("Samples:", nrow(datExpr), "\n")
cat("Genes:", ncol(datExpr), "\n")

## Soft Threshold Selection
cat("\nCalculating soft threshold power...\n")

# Define powers to test
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Calculate soft threshold
sft <- pickSoftThreshold(datExpr, 
                         powerVector = powers,
                         networkType = "signed hybrid",
                         corFnc = "bicor",
                         verbose = 5)

# Plot soft threshold results
pdf("results/figures/WGCNA/soft_threshold_selection.pdf", width = 9, height = 5)
par(mfrow = c(1,2))

# Scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "red")
abline(h = 0.90, col = "blue")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], 
     labels = powers, cex = 0.9, col = "red")

dev.off()

# Select power
if (nrow(datExpr) > 20) {
  recommended_power <- 9  # For signed hybrid with >20 samples
} else {
  recommended_power <- sft$powerEstimate
}

cat("\nRecommended soft threshold power:", recommended_power, "\n") #9

## Create Sample Dendrogram
cat("\nCreating sample dendrogram...\n")

pdf("results/figures/WGCNA/sample_dendrogram.pdf", width = 12, height = 6)
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

# Save WGCNA preparation results
wgcna_prep <- list(
  genes_selected = genes_for_wgcna,
  n_genes = length(genes_for_wgcna),
  soft_threshold = recommended_power,
  sft_results = sft,
  anova_results = eR_anova_detailed
)

save(datExpr, wgcna_prep, genes_for_wgcna, recommended_power, eR_anova_detailed,
     file = "results/checkpoints/08_wgcna_prep.RData")

cat("\n=== WGCNA preparation complete ===\n")
cat("Genes selected for WGCNA:", ncol(datExpr), "\n") #36350 
cat("Recommended power:", recommended_power, "\n")

## ===== 09_wgcna_analysis.R =====
#-------------------------------------------------------#
#     09. WGCNA Network Analysis                        #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoint
cat("Loading WGCNA preparation data...\n")
load_checkpoint("results/checkpoints/08_wgcna_prep.RData")

# Set up WGCNA
library(WGCNA)
enableWGCNAThreads()

# Create output directories
create_dir("results/WGCNA/networks")
create_dir("results/figures/WGCNA/modules")

## One-Step Network Construction
cat("\nConstructing gene network and identifying modules...\n")
cat("This may take several minutes depending on the number of genes...\n")

# Run blockwiseModules
net <- blockwiseModules(
  datExpr,
  power = recommended_power,
  maxBlockSize = 35000,  # Should be larger than number of genes
  TOMType = "signed",
  networkType = "signed hybrid",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.2,  # Start with 0.2, can adjust later
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "results/WGCNA/networks/TOM",
  loadTOMs = TRUE,
  corType = "bicor",
  maxPOutliers = 0.1,
  deepSplit = 2,
  verbose = 3
)

cat("\nModule detection complete!\n")
cat("Number of modules found:", max(net$colors), "\n") #36 

# Module size summary
module_sizes <- table(net$colors)
cat("\nModule sizes:\n")
print(module_sizes)

## Examine Module Merging
cat("\nExamining module merging...\n")

# Get module eigengenes
MEs <- net$MEs
MEDiss <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot module eigengene dendrogram
pdf("results/figures/WGCNA/modules/module_eigengene_clustering.pdf", width = 10, height = 6)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h = 0.2, col = "red", lty = 2)
abline(h = 0.15, col = "blue", lty = 2)
abline(h = 0.1, col = "green", lty = 2)
dev.off()

## Try Different Merge Heights
cat("\nTesting different module merge heights...\n")

merge_heights <- c(0.1, 0.15, 0.2)
net_list <- list()

for (merge_h in merge_heights) {
  cat("Testing merge height:", merge_h, "\n")
  
  net_temp <- blockwiseModules(
    datExpr,
    power = recommended_power,
    maxBlockSize = 35000,
    TOMType = "signed",
    networkType = "signed hybrid",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = merge_h,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    loadTOM = TRUE,
    saveTOMs = FALSE,
    saveTOMFileBase = "results/WGCNA/networks/TOM",
    corType = "bicor",
    maxPOutliers = 0.1,
    deepSplit = 2,
    verbose = 0
  )
  
  net_list[[as.character(merge_h)]] <- net_temp
  cat("  Modules found:", max(net_temp$colors), "\n")
}

## Module Colors
cat("\nAssigning module colors...\n")

# Get the colors for only the genes in the dendrogram
# First, identify which genes are in the dendrogram
dendroGenes <- net$blockGenes[[1]]  # Genes in the first (and likely only) block

# Get colors for just these genes
moduleColors <- labels2colors(net$colors[dendroGenes])

# Verify the length matches
cat("Number of genes in dendrogram:", length(net$dendrograms[[1]]$order), "\n")
cat("Number of module colors:", length(moduleColors), "\n")

# Now plot
pdf("results/figures/WGCNA/modules/gene_dendrogram_modules.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], 
                    moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Also update the mergedColors variable for consistency
mergedColors <- moduleColors

# Check the distribution of modules
table(mergedColors)




## Module-Trait Relationships
# Load the targets data
targets <- read.csv("Targets_Final.csv", row.names = 1)

# Verify it loaded correctly
cat("Dimensions of targets:", dim(targets), "\n")
cat("Column names:", colnames(targets), "\n")

# Now create the trait data
moduleTraitData <- data.frame(
  Narrow = as.numeric(targets$Leaf_type == "Narrow"),
  Broad = as.numeric(targets$Leaf_type == "Broad"),
  TP0 = as.numeric(targets$Timepoint == "TP0"),
  TP1 = as.numeric(targets$Timepoint == "TP1"),
  TP2 = as.numeric(targets$Timepoint == "TP2"),
  TP3 = as.numeric(targets$Timepoint == "TP3"),
  TP4 = as.numeric(targets$Timepoint == "TP4")
)

# Make sure the rows match your expression data
rownames(moduleTraitData) <- rownames(targets)

# Verify the trait data
cat("\nTrait data dimensions:", dim(moduleTraitData), "\n")
cat("First few rows:\n")
head(moduleTraitData)


# Add continuous traits if available
if (exists("logCPM")) {
  # JAG1 expression
  if ("Glyma.20G116200" %in% rownames(logCPM)) {
    moduleTraitData$JAG1_expr <- logCPM["Glyma.20G116200", ]
  }
}

# Calculate correlations
moduleTraitCor <- cor(net$MEs, moduleTraitData, use = "p", method = "pearson")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Create module-trait heatmap
pdf("results/figures/WGCNA/modules/module_trait_relationships.pdf", width = 10, height = 8)
par(mar = c(6, 8.5, 3, 3))

# Display correlations and p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(moduleTraitData),
               yLabels = names(net$MEs),
               ySymbols = names(net$MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()


# Create interaction traits (Leaf type × Timepoint)
moduleTraitData_interactions <- data.frame(
  # Original traits
  moduleTraitData,
  
  # Narrow leaf at each timepoint
  Narrow_TP0 = as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP0"),
  Narrow_TP1 = as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP1"),
  Narrow_TP2 = as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP2"),
  Narrow_TP3 = as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP3"),
  Narrow_TP4 = as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP4"),
  
  # Broad leaf at each timepoint
  Broad_TP0 = as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP0"),
  Broad_TP1 = as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP1"),
  Broad_TP2 = as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP2"),
  Broad_TP3 = as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP3"),
  Broad_TP4 = as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP4")
)

# Recalculate correlations with interactions
moduleTraitCor_full <- cor(net$MEs, moduleTraitData_interactions, use = "p")
moduleTraitPvalue_full <- corPvalueStudent(moduleTraitCor_full, nrow(datExpr))

# Fix for the first PDF
pdf("results/figures/WGCNA/modules/module_trait_relationships_detailed.pdf", width = 16, height = 10)
par(mar = c(8, 8.5, 3, 3))

# Reorder columns for better visualization
col_order <- c("Narrow", "Broad", "TP0", "TP1", "TP2", "TP3", "TP4", 
               "Narrow_TP0", "Narrow_TP1", "Narrow_TP2", "Narrow_TP3", "Narrow_TP4",
               "Broad_TP0", "Broad_TP1", "Broad_TP2", "Broad_TP3", "Broad_TP4")

# Add JAG1 if it exists
if("JAG1_expr" %in% colnames(moduleTraitCor_full)) {
  col_order <- c(col_order, "JAG1_expr")
}

# Reorder the correlation matrix
moduleTraitCor_ordered <- moduleTraitCor_full[, col_order]
moduleTraitPvalue_ordered <- moduleTraitPvalue_full[, col_order]

# Create text matrix with the reordered data
textMatrix_ordered <- paste(signif(moduleTraitCor_ordered, 2), "\n(",
                            signif(moduleTraitPvalue_ordered, 1), ")", sep = "")
dim(textMatrix_ordered) <- dim(moduleTraitCor_ordered)

labeledHeatmap(Matrix = moduleTraitCor_ordered,
               xLabels = colnames(moduleTraitCor_ordered),
               yLabels = names(net$MEs),
               ySymbols = names(net$MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_ordered,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = "Module-trait relationships (detailed with interactions)",
               cex.lab.x = 0.8,
               cex.lab.y = 0.7)

# Add vertical lines to separate trait groups
abline(v = c(2.5, 7.5, 12.5, 17.5), lty = 2, col = "gray")

dev.off()

# Create a second version with timepoints grouped (Narrow and Broad together)
pdf("results/figures/WGCNA/modules/module_trait_interactions_grouped.pdf", width = 12, height = 10)
par(mar = c(10, 8.5, 3, 3))

# Group by timepoint instead of leaf type
grouped_cols <- c("Narrow_TP0", "Broad_TP0", 
                  "Narrow_TP1", "Broad_TP1",
                  "Narrow_TP2", "Broad_TP2",
                  "Narrow_TP3", "Broad_TP3",
                  "Narrow_TP4", "Broad_TP4")

moduleTraitCor_grouped <- moduleTraitCor_full[, grouped_cols]
moduleTraitPvalue_grouped <- moduleTraitPvalue_full[, grouped_cols]

textMatrix_grouped <- paste(signif(moduleTraitCor_grouped, 2), "\n(",
                            signif(moduleTraitPvalue_grouped, 1), ")", sep = "")
dim(textMatrix_grouped) <- dim(moduleTraitCor_grouped)

# Create custom labels that are shorter
custom_labels <- c("N-TP0", "B-TP0", "N-TP1", "B-TP1", 
                   "N-TP2", "B-TP2", "N-TP3", "B-TP3", 
                   "N-TP4", "B-TP4")

labeledHeatmap(Matrix = moduleTraitCor_grouped,
               xLabels = custom_labels,
               yLabels = names(net$MEs),
               ySymbols = names(net$MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_grouped,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = "Module-trait relationships: Leaf type effects at each timepoint",
               cex.lab.x = 1.0,
               cex.lab.y = 0.7)

# Add vertical lines to separate timepoints
abline(v = seq(2.5, 8.5, by = 2), lty = 2, col = "black", lwd = 1.5)

# Add text labels for timepoints
mtext("TP0", side = 1, at = 1.5, line = 7, cex = 0.8, font = 2)
mtext("TP1", side = 1, at = 3.5, line = 7, cex = 0.8, font = 2)
mtext("TP2", side = 1, at = 5.5, line = 7, cex = 0.8, font = 2)
mtext("TP3", side = 1, at = 7.5, line = 7, cex = 0.8, font = 2)
mtext("TP4", side = 1, at = 9.5, line = 7, cex = 0.8, font = 2)

dev.off()

# Create a difference heatmap to show leaf type effects more clearly
pdf("results/figures/WGCNA/modules/module_trait_leaf_differences.pdf", width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))

# Calculate differences (Narrow - Broad) at each timepoint
diff_matrix <- matrix(0, nrow = nrow(moduleTraitCor_full), ncol = 5)
rownames(diff_matrix) <- rownames(moduleTraitCor_full)
colnames(diff_matrix) <- c("TP0", "TP1", "TP2", "TP3", "TP4")

for(i in 0:4) {
  narrow_col <- paste0("Narrow_TP", i)
  broad_col <- paste0("Broad_TP", i)
  diff_matrix[, paste0("TP", i)] <- moduleTraitCor_full[, narrow_col] - moduleTraitCor_full[, broad_col]
}

# Create text showing the difference values
textMatrix_diff <- signif(diff_matrix, 2)
dim(textMatrix_diff) <- dim(diff_matrix)

labeledHeatmap(Matrix = diff_matrix,
               xLabels = colnames(diff_matrix),
               yLabels = rownames(diff_matrix),
               ySymbols = rownames(diff_matrix),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix_diff,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = "Leaf type effect (Narrow - Broad) across timepoints",
               cex.lab.x = 1.2,
               cex.lab.y = 0.7)

dev.off()

# Print summary
cat("\n=== Modules with consistent leaf type effects across all timepoints ===\n")
consistent_modules <- rownames(diff_matrix)[apply(abs(diff_matrix) > 0.3, 1, all)]
if(length(consistent_modules) > 0) {
  for(mod in consistent_modules) {
    cat(sprintf("%s: ", mod))
    cat(sprintf("TP0=%.2f, TP1=%.2f, TP2=%.2f, TP3=%.2f, TP4=%.2f\n", 
                diff_matrix[mod, "TP0"], diff_matrix[mod, "TP1"], 
                diff_matrix[mod, "TP2"], diff_matrix[mod, "TP3"], 
                diff_matrix[mod, "TP4"]))
  }
} else {
  cat("No modules show consistent leaf type effects across all timepoints\n")
}

# Print summary of interesting patterns
cat("\n=== Modules with strong leaf-type × timepoint interactions ===\n")
for(tp in c("TP0", "TP1", "TP2", "TP3", "TP4")) {
  narrow_col <- paste0("Narrow_", tp)
  broad_col <- paste0("Broad_", tp)
  
  # Find modules with opposite effects in narrow vs broad at this timepoint
  diff <- moduleTraitCor_full[, narrow_col] - moduleTraitCor_full[, broad_col]
  strong_diff <- names(diff)[abs(diff) > 0.5]
  
  if(length(strong_diff) > 0) {
    cat("\nAt", tp, "- modules with strong leaf type differences:\n")
    for(mod in strong_diff) {
      cat(sprintf("  %s: Narrow=%.2f, Broad=%.2f (diff=%.2f)\n", 
                  mod, 
                  moduleTraitCor_full[mod, narrow_col],
                  moduleTraitCor_full[mod, broad_col],
                  diff[mod]))
    }
  }
}




## Calculate Module Membership (kME)
cat("\nCalculating module membership values...\n")

# First, let's check if kMEs were calculated correctly
cat("Checking kME calculation...\n")
cat("kME dimensions:", dim(kMEs), "\n")
cat("First few gene names:", head(rownames(kMEs)), "\n")

# Calculate kME for all genes if not done yet
if (!exists("kMEs") || is.null(kMEs)) {
  kMEs <- signedKME(datExpr, net$MEs, corFnc = "bicor")
}

# Hub genes list
hub_genes <- list()

# Get unique module numbers (excluding 0 which is grey/unassigned)
module_numbers <- unique(net$colors)
module_numbers <- module_numbers[module_numbers != 0]

cat("\nProcessing", length(module_numbers), "modules...\n")

# Track successful modules
successful_modules <- 0

for (mod_num in module_numbers) {
  # Get module color
  mod_color <- labels2colors(mod_num)
  
  # Get genes in this module
  module_genes <- colnames(datExpr)[net$colors == mod_num]
  
  if (length(module_genes) == 0) {
    cat("Warning: Module", mod_num, "(", mod_color, ") has no genes\n")
    next
  }
  
  # The kME column name uses the module number
  kME_col <- paste0("kME", mod_num)
  
  # Check if kME column exists
  if (!kME_col %in% colnames(kMEs)) {
    cat("Warning: No kME column", kME_col, "for module", mod_num, "(", mod_color, ")\n")
    next
  }
  
  # Get kME values for module genes
  module_kMEs <- kMEs[module_genes, kME_col, drop = FALSE]
  
  # Create a data frame to preserve names
  kme_df <- data.frame(
    Gene = module_genes,
    kME = module_kMEs[, 1],
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA kME values
  kme_df <- kme_df[!is.na(kme_df$kME), ]
  
  if (nrow(kme_df) == 0) {
    cat("Warning: No valid kME values for module", mod_num, "\n")
    next
  }
  
  # Sort by kME value
  kme_df <- kme_df[order(kme_df$kME, decreasing = TRUE), ]
  
  # Get top hub genes (up to 10)
  n_hubs <- min(10, nrow(kme_df))
  top_hubs_df <- kme_df[1:n_hubs, ]
  
  # Store hub genes
  hub_genes[[paste0("Module_", mod_color)]] <- data.frame(
    Gene = top_hubs_df$Gene,
    kME = top_hubs_df$kME,
    Module_number = mod_num,
    Module_color = mod_color,
    stringsAsFactors = FALSE
  )
  
  successful_modules <- successful_modules + 1
}

cat("\nSuccessfully processed", successful_modules, "modules\n")

# Save hub genes to file if we have any
if (length(hub_genes) > 0) {
  hub_genes_df <- do.call(rbind, hub_genes)
  rownames(hub_genes_df) <- NULL
  
  # Create directory if it doesn't exist
  if (!dir.exists("results/tables/WGCNA")) {
    dir.create("results/tables/WGCNA", recursive = TRUE)
  }
  
  write.csv(hub_genes_df, "results/tables/WGCNA/module_hub_genes.csv", row.names = FALSE)
  
  # Display summary
  cat("\nHub genes identified for", length(hub_genes), "modules\n")
  cat("Total hub genes:", nrow(hub_genes_df), "\n")
  
  # Show example hub genes
  cat("\nExample hub genes from first 3 modules:\n")
  print(head(hub_genes_df, 15))
  
  # Show top hub genes for modules of interest
  if (exists("moduleTraitCor")) {
    # Find modules most correlated with Narrow leaves
    narrow_cors <- moduleTraitCor[, "Narrow"]
    
    # Get the top 3 modules by absolute correlation
    me_names <- names(sort(abs(narrow_cors), decreasing = TRUE)[1:min(3, length(narrow_cors))])
    
    cat("\n\nTop hub genes in modules associated with narrow leaves:\n")
    for (me_name in me_names) {
      # Extract module number from ME name
      mod_num <- as.numeric(gsub("ME", "", me_name))
      mod_color <- labels2colors(mod_num)
      mod_df_name <- paste0("Module_", mod_color)
      
      if (mod_df_name %in% names(hub_genes)) {
        cat("\n", me_name, "(", mod_color, ") - correlation with Narrow: r =", 
            round(narrow_cors[me_name], 3), "\n")
        mod_hub_genes <- hub_genes[[mod_df_name]]
        print(head(mod_hub_genes[, c("Gene", "kME")], 5))
      }
    }
  }
} else {
  cat("\nNo hub genes identified. Check the kME calculation and module assignments.\n")
}

# Create a summary of module sizes
cat("\n\nModule size summary:\n")
module_sizes <- table(labels2colors(net$colors))
module_sizes <- module_sizes[names(module_sizes) != "grey"]
module_sizes <- sort(module_sizes, decreasing = TRUE)
print(head(module_sizes, 10))


# Show which modules have kME values
cat("\n\nModules with kME columns available:\n")
available_kme <- gsub("kME", "", colnames(kMEs))
available_colors <- labels2colors(as.numeric(available_kme))
cat(paste(available_colors, collapse = ", "), "\n")

# Save hub genes
hub_genes_df <- do.call(rbind, hub_genes)
write.csv(hub_genes_df, file = "results/WGCNA/hub_genes_per_module.csv", row.names = FALSE)


# Explore Module Membership patterns
cat("\n=== MODULE MEMBERSHIP ANALYSIS ===\n")

# 1. Check kME distribution for each module
cat("\nkME value distributions by module:\n")
for (mod_num in 1:5) {  # Show first 5 modules as example
  if (mod_num %in% module_numbers) {
    mod_color <- labels2colors(mod_num)
    kME_col <- paste0("kME", mod_num)
    module_genes <- colnames(datExpr)[net$colors == mod_num]
    
    if (kME_col %in% colnames(kMEs) && length(module_genes) > 0) {
      kme_values <- kMEs[module_genes, kME_col]
      cat(sprintf("\nModule %d (%s) - %d genes:\n", mod_num, mod_color, length(module_genes)))
      cat(sprintf("  Mean kME: %.3f\n", mean(kme_values, na.rm = TRUE)))
      cat(sprintf("  Min kME: %.3f, Max kME: %.3f\n", 
                  min(kme_values, na.rm = TRUE), 
                  max(kme_values, na.rm = TRUE)))
    }
  }
}

# 2. Visualize kME distributions
pdf("results/figures/WGCNA/module_membership_distributions.pdf", width = 12, height = 8)
par(mfrow = c(2, 3))

# Plot kME distributions for top 6 modules by size
top_modules <- names(sort(table(net$colors[net$colors != 0]), decreasing = TRUE)[1:6])

for (mod_num in as.numeric(top_modules)) {
  mod_color <- labels2colors(mod_num)
  kME_col <- paste0("kME", mod_num)
  module_genes <- colnames(datExpr)[net$colors == mod_num]
  
  if (kME_col %in% colnames(kMEs) && length(module_genes) > 0) {
    kme_values <- kMEs[module_genes, kME_col]
    
    hist(kme_values, 
         breaks = 30,
         col = mod_color,
         border = "black",
         main = paste("Module", mod_num, "(", mod_color, ")"),
         xlab = "Module Membership (kME)",
         ylab = "Number of genes")
    
    # Add vertical line at 0.7 (common threshold for hub genes)
    abline(v = 0.7, col = "red", lty = 2, lwd = 2)
    
    # Add statistics
    text(0.3, max(hist(kme_values, plot = FALSE)$counts) * 0.9,
         paste("n =", length(module_genes), "\n",
               "Hub genes (kME>0.7):", sum(kme_values > 0.7, na.rm = TRUE)),
         pos = 4)
  }
}
dev.off()

# 3. Create a summary table of module membership
module_kme_summary <- data.frame()

for (mod_num in module_numbers) {
  mod_color <- labels2colors(mod_num)
  kME_col <- paste0("kME", mod_num)
  module_genes <- colnames(datExpr)[net$colors == mod_num]
  
  if (kME_col %in% colnames(kMEs) && length(module_genes) > 0) {
    kme_values <- kMEs[module_genes, kME_col]
    
    summary_row <- data.frame(
      Module_Number = mod_num,
      Module_Color = mod_color,
      Total_Genes = length(module_genes),
      Mean_kME = mean(kme_values, na.rm = TRUE),
      Median_kME = median(kme_values, na.rm = TRUE),
      Min_kME = min(kme_values, na.rm = TRUE),
      Max_kME = max(kme_values, na.rm = TRUE),
      Hub_Genes_07 = sum(kme_values > 0.7, na.rm = TRUE),
      Hub_Genes_08 = sum(kme_values > 0.8, na.rm = TRUE),
      Hub_Genes_09 = sum(kme_values > 0.9, na.rm = TRUE)
    )
    
    module_kme_summary <- rbind(module_kme_summary, summary_row)
  }
}

# Sort by mean kME (shows module cohesiveness)
module_kme_summary <- module_kme_summary[order(module_kme_summary$Mean_kME, decreasing = TRUE), ]

write.csv(module_kme_summary, "results/tables/WGCNA/module_membership_summary.csv", row.names = FALSE)

cat("\n\nModule Membership Summary (top 10 most cohesive modules):\n")
print(head(module_kme_summary, 10))

# 4. Identify genes with high membership in multiple modules
cat("\n\nChecking for genes with high membership in multiple modules...\n")

# Find genes with kME > 0.6 in multiple modules
high_kme_threshold <- 0.6
multi_module_genes <- c()

for (gene in rownames(kMEs)) {
  high_modules <- c()
  
  for (mod_num in module_numbers) {
    kME_col <- paste0("kME", mod_num)
    if (kME_col %in% colnames(kMEs)) {
      if (!is.na(kMEs[gene, kME_col]) && kMEs[gene, kME_col] > high_kme_threshold) {
        high_modules <- c(high_modules, mod_num)
      }
    }
  }
  
  if (length(high_modules) > 1) {
    multi_module_genes <- c(multi_module_genes, gene)
  }
}

cat("Genes with high membership (kME >", high_kme_threshold, ") in multiple modules:", 
    length(multi_module_genes), "\n")

# 5. Compare assigned module vs highest kME
cat("\n\nChecking module assignment accuracy...\n")
misassigned <- 0
for (gene in colnames(datExpr)[1:1000]) {  # Check first 1000 genes as example
  assigned_module <- net$colors[gene]
  if (assigned_module == 0) next
  
  # Find module with highest kME
  gene_kmes <- kMEs[gene, paste0("kME", module_numbers)]
  if (any(!is.na(gene_kmes))) {
    best_module <- module_numbers[which.max(gene_kmes)]
    
    if (best_module != assigned_module) {
      misassigned <- misassigned + 1
    }
  }
}

cat("Genes checked: 1000\n")
cat("Genes in different module than highest kME:", misassigned, 
    "(", round(misassigned/1000*100, 1), "%)\n")


## Module Preservation (if comparing conditions)
cat("\nChecking module preservation between leaf types...\n")

# Separate data by leaf type
broad_samples <- rownames(datExpr)[targets$Leaf_type == "Broad"]
narrow_samples <- rownames(datExpr)[targets$Leaf_type == "Narrow"]

multiExpr <- list(
  Broad = list(data = datExpr[broad_samples, ]),
  Narrow = list(data = datExpr[narrow_samples, ])
)

# Set up multiColor
multiColor <- list(Broad = net$colors, Narrow = net$colors)

# Note: Module preservation analysis is computationally intensive
# Uncomment below to run (may take 30+ minutes)

mp <- modulePreservation(multiExpr,
                        multiColor,
                        referenceNetworks = 1,
                        nPermutations = 100,
                        randomSeed = 12345,
                        quickCor = 0,
                        verbose = 3)

# Save WGCNA results
wgcna_results <- list(
  module_colors = mergedColors,
  module_eigengenes = net$MEs,
  module_trait_cor = moduleTraitCor,
  module_trait_pval = moduleTraitPvalue,
  kME_values = kMEs,
  hub_genes = hub_genes_df,
  module_sizes = module_sizes
)

save(net, wgcna_results, kMEs, mergedColors, moduleTraitCor, moduleTraitPvalue,
     file = "results/checkpoints/09_wgcna_analysis.RData")

cat("\n=== WGCNA analysis complete ===\n")
cat("Total modules:", max(net$colors), "\n")
cat("Module-trait correlations calculated\n")
cat("Hub genes identified\n")

## ===== 10_wgcna_visualization.R =====
#-------------------------------------------------------#
#     10. WGCNA Visualization                           #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoints
cat("Loading WGCNA results...\n")
load_checkpoint("results/checkpoints/09_wgcna_analysis.RData")
load_checkpoint("results/checkpoints/08_wgcna_prep.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")

# Load visualization functions
source("src/JDplot.mat.R")  # If available
source("src/JDmoduleHeatmaps.R")  # If available

# Create output directories
create_dir("results/figures/WGCNA/heatmaps")
create_dir("results/figures/WGCNA/barplots")
create_dir("results/figures/WGCNA/networks")

## Module Heatmaps
cat("\nCreating module heatmaps...\n")

# Function to create module heatmap - FIXED VERSION
create_module_heatmap <- function(module_num, expr_data, module_colors, sample_info) {
  # Get module genes
  module_genes <- names(module_colors)[module_colors == module_num]
  
  if (length(module_genes) < 2) return(NULL)
  
  # Subset and scale expression data
  module_data <- expr_data[module_genes, ]
  module_scaled <- t(scale(t(module_data)))
  
  # Create heatmap
  filename <- paste0("results/figures/WGCNA/heatmaps/module_", module_num, "_heatmap.pdf")
  pdf(filename, width = 10, height = 8)
  
  # Order samples by timepoint
  sample_order <- order(sample_info$Timepoint, sample_info$Leaf_type)
  
  # Define colors for leaf types
  leaf_colors <- c("Broad" = "orange", "Narrow" = "purple")
  col_colors <- leaf_colors[sample_info$Leaf_type[sample_order]]
  
  # Define colors for timepoints
  time_colors <- rainbow(5)[as.numeric(factor(sample_info$Timepoint))]
  time_colors <- time_colors[sample_order]
  
  # Create color matrix for multiple annotation rows
  col_annotations <- rbind(
    LeafType = col_colors,
    Timepoint = time_colors
  )
  
  # Create heatmap with proper labels
  labCol <- sample_info$Sample[sample_order]
  if ("Label" %in% colnames(sample_info)) {
    labCol <- sample_info$Label[sample_order]
  }
  
  # Create heatmap
  heatmap.2(module_scaled[, sample_order],
            col = colorRampPalette(c("blue", "white", "red"))(100),
            trace = "none",
            density.info = "none",
            labRow = ifelse(length(module_genes) > 50, "", rownames(module_scaled)),
            labCol = labCol,
            ColSideColors = col_colors,  # Now it's a character vector
            margins = c(10, 8),
            main = paste("Module", module_num, "-", length(module_genes), "genes"),
            cexCol = 0.8)
  
  # Add legend
  legend("topright", 
         legend = c("Broad", "Narrow"), 
         fill = c("orange", "purple"),
         title = "Leaf Type",
         bty = "n")
  
  dev.off()
}


# Load required libraries
library(gplots)
library(ggplot2)  # For the barplots later
# Load targets data
targets <- read.csv("Targets_Final.csv", row.names = 1)
# Create heatmaps for all modules
for (mod in unique(net$colors)) {
  if (mod == 0) next  # Skip grey module
  cat("Creating heatmap for module", mod, "\n")
  create_module_heatmap(mod, logCPM, net$colors, targets)
}











## Module Eigengene Barplots
cat("\nCreating module eigengene barplots...\n")



# Function to create ME barplot - DEBUGGED VERSION
create_ME_barplot <- function(ME_data, module_num, sample_info) {
  ME_col <- paste0("ME", module_num)
  if (!ME_col %in% colnames(ME_data)) {
    cat("  Warning: Column", ME_col, "not found in ME_data\n")
    return(NULL)
  }
  
  # Debug: Check dimensions
  cat("  ME_data rows:", nrow(ME_data), ", sample_info rows:", nrow(sample_info), "\n")
  
  # Make sure we have matching samples
  if (nrow(ME_data) != nrow(sample_info)) {
    cat("  Warning: Row count mismatch\n")
    return(NULL)
  }
  
  # Prepare data - ensure all required columns exist
  plot_data <- data.frame(
    ME = ME_data[, ME_col],
    Sample = rownames(ME_data),
    Timepoint = sample_info$Timepoint,
    LeafType = sample_info$Leaf_type,
    stringsAsFactors = FALSE
  )
  
  # Check if plot_data is valid
  if (nrow(plot_data) == 0 || any(is.na(plot_data$LeafType))) {
    cat("  Warning: Invalid plot data\n")
    return(NULL)
  }
  
  # Simple bar plot without faceting first
  p1 <- ggplot(plot_data, aes(x = Sample, y = ME, fill = LeafType)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          legend.position = "top") +
    labs(title = paste("Module", module_num, "Eigengene Expression"),
         y = "Module Eigengene",
         x = "") +
    scale_fill_manual(values = c("Broad" = "#E69F00", "Narrow" = "#56B4E9"))
  
  filename1 <- paste0("results/figures/WGCNA/barplots/ME_barplot_simple_module_", module_num, ".pdf")
  ggsave(filename1, plot = p1, width = 12, height = 6)
  
  # Boxplot by timepoint
  p2 <- ggplot(plot_data, aes(x = Timepoint, y = ME, fill = LeafType)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               size = 2, alpha = 0.8) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(title = paste("Module", module_num, "Eigengene Expression by Timepoint"),
         y = "Module Eigengene",
         x = "Timepoint") +
    scale_fill_manual(values = c("Broad" = "#E69F00", "Narrow" = "#56B4E9"))
  
  filename2 <- paste0("results/figures/WGCNA/barplots/ME_boxplot_module_", module_num, ".pdf")
  ggsave(filename2, plot = p2, width = 8, height = 6)
  
  # Line plot showing temporal pattern
  mean_data <- aggregate(ME ~ Timepoint + LeafType, data = plot_data, FUN = mean)
  
  p3 <- ggplot(mean_data, aes(x = Timepoint, y = ME, color = LeafType, group = LeafType)) +
    geom_line(size = 1.5) +
    geom_point(size = 3) +
    theme_minimal() +
    theme(legend.position = "top") +
    labs(title = paste("Module", module_num, "Average Eigengene Expression Over Time"),
         y = "Mean Module Eigengene",
         x = "Timepoint") +
    scale_color_manual(values = c("Broad" = "#E69F00", "Narrow" = "#56B4E9"))
  
  filename3 <- paste0("results/figures/WGCNA/barplots/ME_lineplot_module_", module_num, ".pdf")
  ggsave(filename3, plot = p3, width = 8, height = 6)
  
  cat("  Successfully created plots for module", module_num, "\n")
}

# Test with one module first
test_module <- as.numeric(sub("ME", "", rownames(moduleTraitCor)[significant_modules[1]]))
cat("\nTesting with module", test_module, "\n")
create_ME_barplot(net$MEs, test_module, targets)

# If successful, run for all
if (file.exists(paste0("results/figures/WGCNA/barplots/ME_barplot_simple_module_", test_module, ".pdf"))) {
  cat("\nTest successful! Creating plots for all significant modules...\n")
  
  for (i in significant_modules) {
    module_num <- as.numeric(sub("ME", "", rownames(moduleTraitCor)[i]))
    if (module_num == 0) next
    
    tryCatch({
      create_ME_barplot(net$MEs, module_num, targets)
    }, error = function(e) {
      cat("  Error with module", module_num, ":", e$message, "\n")
    })
  }
}





## Module Eigengene Relationships
cat("\nVisualizing module eigengene relationships...\n")

# Calculate ME adjacency
ME_adj <- adjacency(net$MEs, power = 1, type = "signed")

# Create ME network plot
pdf("results/figures/WGCNA/ME_network_plot.pdf", width = 10, height = 10)
par(cex = 1.0)
plotEigengeneNetworks(net$MEs, 
                      "Module Eigengene Network",
                      marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(3, 4, 1, 2),
                      plotDendrograms = TRUE,
                      xLabelsAngle = 90)
dev.off()

## Gene Significance vs Module Membership Plots
cat("\nCreating GS vs MM plots for significant modules...\n")

# Load DE results for gene significance
load_checkpoint("results/checkpoints/06_de_results.RData")

# Function to create GS vs MM plot
create_GS_MM_plot <- function(module_num, GS_values, kME_values, module_colors) {
  # Get module genes
  module_genes <- names(module_colors)[module_colors == module_num]
  
  # Get kME for this module
  kME_col <- paste0("kME", module_num)
  
  # Prepare data
  plot_data <- data.frame(
    Gene = module_genes,
    MM = kME_values[module_genes, kME_col],
    GS = GS_values[module_genes]
  )
  
  # Remove NAs
  plot_data <- na.omit(plot_data)
  
  # Calculate correlation
  cor_value <- cor(plot_data$MM, plot_data$GS, use = "p")
  
  # Create plot
  p <- ggplot(plot_data, aes(x = MM, y = GS)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(title = paste("Module", module_num, "- Gene Significance vs Module Membership"),
         subtitle = paste("Correlation:", round(cor_value, 3)),
         x = paste("Module Membership in", module_num),
         y = "Gene Significance")
  
  filename <- paste0("results/figures/WGCNA/GS_vs_MM_module_", module_num, ".pdf")
  ggsave(filename, plot = p, width = 8, height = 6)
}

# Calculate gene significance (using p-values from a key contrast)
key_contrast <- "B_462A_4vs0"  # Modify as needed
if (key_contrast %in% names(all_results_treat)) {
  GS_pvalue <- -log10(all_results_treat[[key_contrast]]$PValue)
  names(GS_pvalue) <- rownames(all_results_treat[[key_contrast]])
  
  # Create plots for significant modules
  for (i in significant_modules) {
    module_num <- as.numeric(sub("ME", "", rownames(moduleTraitCor)[i]))
    if (module_num == 0) next
    
    cat("Creating GS vs MM plot for module", module_num, "\n")
    create_GS_MM_plot(module_num, GS_pvalue, kMEs, net$colors)
  }
}













## Module Network Visualization
cat("\nCreating module network visualizations...\n")

# First, check if TOM file exists and load it properly
TOM_file <- "results/WGCNA/networks/TOM-block.1.RData"
TOM_loaded <- FALSE

if (file.exists(TOM_file)) {
  cat("Loading TOM file...\n")
  # Load and check what objects are in the file
  load_env <- new.env()
  load(TOM_file, envir = load_env)
  cat("Objects in TOM file:", ls(load_env), "\n")
  
  # The TOM matrix might have different names
  if ("TOM" %in% ls(load_env)) {
    TOM <- load_env$TOM
    TOM_loaded <- TRUE
  } else if ("dissTOM" %in% ls(load_env)) {
    TOM <- 1 - load_env$dissTOM  # Convert dissimilarity to similarity
    TOM_loaded = TRUE
  } else {
    # Try to find the TOM object
    tom_objects <- grep("TOM", ls(load_env), value = TRUE)
    if (length(tom_objects) > 0) {
      TOM <- get(tom_objects[1], envir = load_env)
      TOM_loaded <- TRUE
      cat("Loaded TOM as:", tom_objects[1], "\n")
    }
  }
  
  if (TOM_loaded) {
    cat("TOM dimensions:", dim(TOM), "\n")
    cat("TOM class:", class(TOM), "\n")
  }
} else {
  cat("TOM file not found at:", TOM_file, "\n")
}

# Create network visualizations if TOM is loaded
if (TOM_loaded) {
  # Select top modules to visualize
  module_nums <- unique(net$colors)
  module_nums <- module_nums[module_nums != 0]  # Remove grey
  
  # Get top 5 modules by size
  module_sizes <- table(net$colors)
  top_modules <- as.numeric(names(sort(module_sizes[names(module_sizes) != "0"], 
                                       decreasing = TRUE)[1:min(5, length(module_sizes)-1)]))
  
  for (mod in top_modules) {
    cat("Creating network visualization for module", mod, "\n")
    
    # Get module genes
    module_genes <- colnames(datExpr)[net$colors == mod]
    
    # Limit to top hub genes if module is large
    if (length(module_genes) > 50) {
      kME_col <- paste0("kME", mod)
      if (kME_col %in% colnames(kMEs)) {
        module_kMEs <- kMEs[module_genes, kME_col]
        module_kMEs <- module_kMEs[!is.na(module_kMEs)]
        if (length(module_kMEs) > 0) {
          top_genes <- names(sort(module_kMEs, decreasing = TRUE)[1:min(50, length(module_kMEs))])
          module_genes <- top_genes
        }
      }
    }
    
    # Get indices of module genes in the TOM matrix
    # The TOM matrix should have the same order as datExpr
    all_genes <- colnames(datExpr)
    
    # For blockwise modules, only genes in blocks are in TOM
    if (exists("net$blockGenes")) {
      block_genes <- all_genes[net$blockGenes[[1]]]
      module_indices <- which(block_genes %in% module_genes)
      
      if (length(module_indices) > 1) {
        # Subset TOM
        module_TOM <- as.matrix(TOM[module_indices, module_indices])
        
        # Create network plot
        pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_network.pdf"), 
            width = 10, height = 10)
        
        # Set gene names
        gene_names <- block_genes[module_indices]
        rownames(module_TOM) <- gene_names
        colnames(module_TOM) <- gene_names
        
        # Create heatmap with better visualization
        col_scheme <- colorRampPalette(c("white", "yellow", "red"))(100)
        
        # If too many genes, don't show labels
        show_labels <- length(module_indices) <= 30
        
        heatmap(module_TOM,
                Rowv = TRUE, Colv = TRUE,  # Allow clustering
                scale = "none",
                col = col_scheme,
                main = paste("Module", mod, "TOM -", length(module_indices), "genes"),
                labRow = if(show_labels) gene_names else FALSE,
                labCol = if(show_labels) gene_names else FALSE,
                margins = c(8, 8))
        
        dev.off()
        cat("  Created network plot for module", mod, "\n")
      } else {
        cat("  Too few genes in module", mod, "for network plot\n")
      }
    }
  }
} else {
  cat("\nSkipping network visualization - TOM not available\n")
  cat("To create network plots, ensure TOM is saved during network construction\n")
}




# Alternative: Create correlation-based network for top modules
cat("\nCreating correlation-based networks for top modules...\n")

# First, let's check module sizes
cat("Top module sizes:\n")
for (mod in top_modules[1:min(5, length(top_modules))]) {
  module_size <- sum(net$colors == mod)
  cat("  Module", mod, ":", module_size, "genes\n")
}

# Now create correlation plots with better size handling
for (mod in top_modules[1:min(3, length(top_modules))]) {
  # Get module genes
  module_genes <- colnames(datExpr)[net$colors == mod]
  original_size <- length(module_genes)
  
  cat("\nProcessing module", mod, "with", original_size, "genes\n")
  
  # Adjust size limits and select top hub genes if needed
  if (length(module_genes) > 50) {
    # Limit to top 50 hub genes
    kME_col <- paste0("kME", mod)
    if (kME_col %in% colnames(kMEs)) {
      # Get kME values ensuring we maintain gene names
      module_kMEs <- kMEs[module_genes, kME_col, drop = TRUE]
      
      # Remove NAs and keep track of gene names
      valid_indices <- !is.na(module_kMEs)
      valid_genes <- module_genes[valid_indices]
      valid_kMEs <- module_kMEs[valid_indices]
      
      if (length(valid_kMEs) > 0) {
        # Sort and select top genes
        sort_order <- order(valid_kMEs, decreasing = TRUE)
        n_select <- min(50, length(valid_kMEs))
        module_genes <- valid_genes[sort_order[1:n_select]]
        cat("  Selected top", n_select, "hub genes based on kME\n")
      } else {
        # If all kME values are NA, just take first 50
        module_genes <- module_genes[1:min(50, length(module_genes))]
        cat("  No valid kME values, selected first 50 genes\n")
      }
    } else {
      # If no kME column exists, just take first 50
      module_genes <- module_genes[1:min(50, length(module_genes))]
      cat("  No kME column found, selected first 50 genes\n")
    }
  }
  
  # Check final size
  if (length(module_genes) >= 5) {
    cat("  Creating plots with", length(module_genes), "genes\n")
    
    # Calculate correlations
    module_expr <- t(datExpr[, module_genes])
    module_cor <- cor(t(module_expr), use = "p")
    
    # Create heatmap
    pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_correlation.pdf"),
        width = 12, height = 12)
    
    # Determine if we should show labels
    show_labels <- length(module_genes) <= 30
    
    # Handle potential issues with gene names
    if (is.null(rownames(module_cor))) {
      rownames(module_cor) <- module_genes
      colnames(module_cor) <- module_genes
    }
    
    heatmap.2(module_cor,
              trace = "none",
              density.info = "none",
              col = colorRampPalette(c("blue", "white", "red"))(100),
              main = paste("Module", mod, "Gene Correlations\n",
                           length(module_genes), "genes shown out of", original_size, "total"),
              labRow = if(show_labels) rownames(module_cor) else FALSE,
              labCol = if(show_labels) colnames(module_cor) else FALSE,
              margins = if(show_labels) c(10, 10) else c(5, 5),
              key.title = "Correlation",
              symkey = TRUE,
              symbreaks = TRUE)
    
    dev.off()
    cat("  Created correlation heatmap\n")
    
    # Also create a simplified visualization for large modules
    if (length(module_genes) > 30) {
      # Create a summary plot showing distribution of correlations
      pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_correlation_summary.pdf"),
          width = 10, height = 8)
      
      # Get upper triangle of correlation matrix (excluding diagonal)
      cor_values <- module_cor[upper.tri(module_cor)]
      
      par(mfrow = c(2, 1))
      
      # Histogram of correlations
      hist(cor_values, 
           breaks = 50, 
           main = paste("Module", mod, "- Distribution of Gene-Gene Correlations"),
           xlab = "Correlation coefficient",
           col = "lightblue",
           xlim = c(-1, 1))
      abline(v = mean(cor_values), col = "red", lwd = 2, lty = 2)
      abline(v = median(cor_values), col = "blue", lwd = 2, lty = 2)
      legend("topright", 
             legend = c(paste("Mean:", round(mean(cor_values), 3)),
                        paste("Median:", round(median(cor_values), 3))),
             col = c("red", "blue"), lty = 2, lwd = 2)
      
      # Density plot
      plot(density(cor_values), 
           main = paste("Module", mod, "- Density of Gene-Gene Correlations"),
           xlab = "Correlation coefficient",
           xlim = c(-1, 1))
      polygon(density(cor_values), col = rgb(0, 0, 1, 0.3))
      
      dev.off()
      cat("  Created correlation summary plot\n")
    }
    
  } else {
    cat("  Module", mod, "has too few genes (", length(module_genes), ") for correlation plot\n")
  }
}

cat("\nModule network visualization complete!\n")




# Extended correlation-based network analysis for more modules
cat("\n\n=== Extended Module Network Analysis ===\n")

# Get more modules - let's do top 10 by size
module_sizes <- table(net$colors)
module_sizes <- module_sizes[names(module_sizes) != "0"]  # Remove grey
top_modules_extended <- as.numeric(names(sort(module_sizes, decreasing = TRUE)[1:min(10, length(module_sizes))]))

cat("Processing top 10 modules by size:\n")
for (i in 1:length(top_modules_extended)) {
  cat(sprintf("  %d. Module %d: %d genes (%s)\n", 
              i, 
              top_modules_extended[i], 
              module_sizes[as.character(top_modules_extended[i])],
              labels2colors(top_modules_extended[i])))
}

# Process all top 10 modules
for (mod in top_modules_extended) {
  # Get module genes
  module_genes <- colnames(datExpr)[net$colors == mod]
  original_size <- length(module_genes)
  mod_color <- labels2colors(mod)
  
  cat(sprintf("\nProcessing module %d (%s) with %d genes\n", mod, mod_color, original_size))
  
  # Select genes based on module size
  if (length(module_genes) > 50) {
    # For large modules, select top hub genes
    kME_col <- paste0("kME", mod)
    if (kME_col %in% colnames(kMEs)) {
      module_kMEs <- kMEs[module_genes, kME_col, drop = TRUE]
      valid_indices <- !is.na(module_kMEs)
      valid_genes <- module_genes[valid_indices]
      valid_kMEs <- module_kMEs[valid_indices]
      
      if (length(valid_kMEs) > 0) {
        sort_order <- order(valid_kMEs, decreasing = TRUE)
        n_select <- min(50, length(valid_kMEs))
        module_genes <- valid_genes[sort_order[1:n_select]]
        cat(sprintf("  Selected top %d hub genes (kME range: %.3f to %.3f)\n", 
                    n_select, 
                    max(valid_kMEs), 
                    valid_kMEs[sort_order[n_select]]))
      }
    }
  }
  
  if (length(module_genes) >= 5) {
    # Calculate correlations
    module_expr <- t(datExpr[, module_genes])
    module_cor <- cor(t(module_expr), use = "p")
    
    # Set row and column names
    rownames(module_cor) <- module_genes
    colnames(module_cor) <- module_genes
    
    # 1. Create correlation heatmap
    pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_", mod_color, "_correlation.pdf"),
        width = 12, height = 12)
    
    show_labels <- length(module_genes) <= 30
    
    heatmap.2(module_cor,
              trace = "none",
              density.info = "none",
              col = colorRampPalette(c("navy", "white", "darkred"))(100),
              main = sprintf("Module %d (%s) Gene Correlations\n%d genes shown out of %d total", 
                             mod, mod_color, length(module_genes), original_size),
              labRow = if(show_labels) rownames(module_cor) else FALSE,
              labCol = if(show_labels) colnames(module_cor) else FALSE,
              margins = if(show_labels) c(10, 10) else c(5, 5),
              key.title = "Correlation",
              symkey = TRUE,
              symbreaks = TRUE)
    
    dev.off()
    cat("  ✓ Created correlation heatmap\n")
    
    # 2. Create network graph for smaller modules
    if (length(module_genes) <= 40) {
      library(igraph)
      
      # Try different correlation thresholds
      thresholds <- c(0.8, 0.7, 0.6)
      threshold_used <- NULL
      
      for (cor_threshold in thresholds) {
        adj_matrix <- module_cor
        adj_matrix[abs(adj_matrix) < cor_threshold] <- 0
        diag(adj_matrix) <- 0
        
        # Check if we have enough edges
        n_edges <- sum(adj_matrix != 0) / 2
        if (n_edges >= 10 && n_edges <= 200) {
          threshold_used <- cor_threshold
          break
        }
      }
      
      if (!is.null(threshold_used)) {
        g <- graph_from_adjacency_matrix(adj_matrix, 
                                         mode = "undirected", 
                                         weighted = TRUE,
                                         diag = FALSE)
        
        if (ecount(g) > 0) {
          pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_", mod_color, "_network_graph.pdf"),
              width = 12, height = 12)
          
          # Set layout
          set.seed(123)  # For reproducible layout
          if (vcount(g) <= 20) {
            layout <- layout_with_kk(g)
          } else {
            layout <- layout_with_fr(g)
          }
          
          # Set node properties
          V(g)$label <- V(g)$name
          
          # Color nodes by kME if available
          if (paste0("kME", mod) %in% colnames(kMEs)) {
            node_kME <- kMEs[V(g)$name, paste0("kME", mod)]
            node_kME[is.na(node_kME)] <- 0
            
            # Create color gradient
            color_palette <- colorRampPalette(c("lightgray", "yellow", "orange", "red"))(100)
            node_colors <- color_palette[cut(node_kME, breaks = 100, labels = FALSE)]
            V(g)$color <- node_colors
            
            # Size by kME
            V(g)$size <- 5 + 25 * ((node_kME - min(node_kME)) / (max(node_kME) - min(node_kME)))
          } else {
            V(g)$color <- "lightblue"
            V(g)$size <- 15
          }
          
          # Set edge properties
          E(g)$width <- abs(E(g)$weight) * 5
          edge_colors <- ifelse(E(g)$weight > 0, 
                                rgb(1, 0, 0, abs(E(g)$weight)), 
                                rgb(0, 0, 1, abs(E(g)$weight)))
          E(g)$color <- edge_colors
          
          # Plot
          plot(g,
               layout = layout,
               vertex.label.cex = 0.7,
               vertex.label.color = "black",
               vertex.frame.color = "black",
               edge.curved = 0.1,
               main = sprintf("Module %d (%s) Network\n%d genes, %d edges (|cor| > %.1f)",
                              mod, mod_color, vcount(g), ecount(g), threshold_used))
          
          # Add legend for node colors
          if (paste0("kME", mod) %in% colnames(kMEs)) {
            legend("bottomright",
                   legend = c("Low kME", "High kME"),
                   fill = c("lightgray", "red"),
                   title = "Module\nMembership",
                   bty = "n")
          }
          
          dev.off()
          cat(sprintf("  ✓ Created network graph (%d nodes, %d edges)\n", 
                      vcount(g), ecount(g)))
        }
      } else {
        cat("  - Skipped network graph (no suitable threshold found)\n")
      }
    }
    
    # 3. Create summary statistics
    pdf(paste0("results/figures/WGCNA/networks/module_", mod, "_", mod_color, "_summary.pdf"),
        width = 12, height = 10)
    
    par(mfrow = c(2, 2))
    
    # Get correlation values
    cor_values <- module_cor[upper.tri(module_cor)]
    
    # a. Histogram
    hist(cor_values, 
         breaks = 50, 
         main = sprintf("Module %d (%s) - Correlation Distribution", mod, mod_color),
         xlab = "Correlation coefficient",
         col = mod_color,
         xlim = c(-1, 1))
    abline(v = mean(cor_values), col = "red", lwd = 2, lty = 2)
    abline(v = 0, col = "black", lwd = 1)
    
    # b. Boxplot comparison with other modules
    if (length(top_modules_extended) > 1) {
      all_cors <- list()
      for (other_mod in top_modules_extended[1:min(5, length(top_modules_extended))]) {
        other_genes <- colnames(datExpr)[net$colors == other_mod]
        if (length(other_genes) > 50) {
          other_genes <- other_genes[1:50]
        }
        if (length(other_genes) >= 5) {
          other_cor <- cor(t(t(datExpr[, other_genes])), use = "p")
          all_cors[[paste0("M", other_mod)]] <- other_cor[upper.tri(other_cor)]
        }
      }
      
      boxplot(all_cors,
              main = "Correlation Comparison Across Modules",
              ylab = "Correlation coefficient",
              col = rainbow(length(all_cors)),
              las = 2)
      abline(h = 0, lty = 2)
    }
    
    # c. Module eigengene plot
    me_col <- paste0("ME", mod)
    if (me_col %in% colnames(net$MEs)) {
      me_values <- net$MEs[, me_col]
      plot(me_values, 
           type = "h",
           main = sprintf("Module %d (%s) Eigengene Values", mod, mod_color),
           xlab = "Sample",
           ylab = "Eigengene value",
           col = mod_color,
           lwd = 2)
      abline(h = 0, lty = 2)
    }
    
    # d. Connectivity distribution
    connectivity <- rowSums(abs(module_cor)) - 1  # Subtract 1 for self-correlation
    hist(connectivity,
         breaks = 20,
         main = sprintf("Module %d (%s) - Gene Connectivity", mod, mod_color),
         xlab = "Sum of absolute correlations",
         col = mod_color)
    
    # Add text summary
    mtext(sprintf("Mean cor: %.3f | SD: %.3f | Density: %.1f%%", 
                  mean(cor_values), 
                  sd(cor_values),
                  100 * sum(abs(cor_values) > 0.5) / length(cor_values)),
          side = 3, line = -2, outer = TRUE)
    
    dev.off()
    cat("  ✓ Created summary plots\n")
    
  } else {
    cat(sprintf("  - Module %d has too few genes (%d) for analysis\n", mod, length(module_genes)))
  }
}

cat("\n\nModule network analysis complete!\n")
cat(sprintf("Created visualizations for %d modules\n", length(top_modules_extended)))

# Summary report
cat("\n=== Summary of Created Files ===\n")
cat("Location: results/figures/WGCNA/networks/\n")
cat("\nFor each module:\n")
cat("  - module_X_COLOR_correlation.pdf: Gene-gene correlation heatmap\n")
cat("  - module_X_COLOR_network_graph.pdf: Network graph (if ≤40 genes)\n")
cat("  - module_X_COLOR_summary.pdf: Statistical summary plots\n")


## Summary Report
cat("\nCreating WGCNA summary report...\n")

# Module summary statistics
module_summary <- data.frame(
  Module = 0:max(net$colors),
  Size = as.vector(table(net$colors)),
  Primary_Trait = NA,
  Max_Correlation = NA,
  P_value = NA
)

# Find primary trait association for each module
for (i in 1:nrow(module_summary)) {
  if (module_summary$Module[i] == 0) next
  
  me_row <- paste0("ME", module_summary$Module[i])
  if (me_row %in% rownames(moduleTraitCor)) {
    max_idx <- which.max(abs(moduleTraitCor[me_row, ]))
    module_summary$Primary_Trait[i] <- colnames(moduleTraitCor)[max_idx]
    module_summary$Max_Correlation[i] <- moduleTraitCor[me_row, max_idx]
    module_summary$P_value[i] <- moduleTraitPvalue[me_row, max_idx]
  }
}

# Save summary
write.csv(module_summary, file = "results/WGCNA/module_summary_statistics.csv", row.names = FALSE)

# Save visualization checkpoint
save_checkpoint("results/checkpoints/10_wgcna_visualization_complete.RData")

cat("\n=== WGCNA visualization complete ===\n")
cat("Heatmaps created for all modules\n")
cat("Significant module-trait relationships visualized\n")
cat("Results saved to: results/figures/WGCNA/\n")

## ===== 11_gene_sets.R =====
#--------------------------------------------------------#
#     11. Gene Set Analysis                             #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoints
cat("Loading data for gene set analysis...\n")
load_checkpoint("results/checkpoints/06_de_results.RData")
load_checkpoint("results/checkpoints/09_wgcna_analysis.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")

# Load additional libraries
library(gplots)
library(RColorBrewer)

# Create output directories
create_dir("results/gene_sets")
create_dir("results/figures/gene_sets")

## Define Gene Sets of Interest
cat("\nDefining gene sets...\n")

# 1. Core leaf shape genes (significant in multiple comparisons)
broad_core <- Reduce(intersect, list(
  rownames(all_results_treat$B_462A_4vs0)[all_results_treat$B_462A_4vs0$FDR < 0.05],
  rownames(all_results_treat$B_LD11_4vs0)[all_results_treat$B_LD11_4vs0$FDR < 0.05]
))

narrow_core <- Reduce(intersect, list(
  rownames(all_results_treat$N_713B_4vs0)[all_results_treat$N_713B_4vs0$FDR < 0.05],
  rownames(all_results_treat$N_745_4vs0)[all_results_treat$N_745_4vs0$FDR < 0.05]
))

# 2. Timepoint-specific genes
tp_specific_genes <- list()
for (tp in c("1vs0", "2vs0", "3vs0", "4vs0")) {
  contrasts <- grep(tp, names(all_results_treat), value = TRUE)
  
  # Get genes significant in at least 2 lines
  sig_genes <- lapply(contrasts, function(x) {
    rownames(all_results_treat[[x]])[all_results_treat[[x]]$FDR < 0.05]
  })
  
  tp_specific_genes[[tp]] <- Reduce(union, sig_genes)
}

# 3. Early response genes (TP1 specific)
early_response <- setdiff(tp_specific_genes[["1vs0"]], 
                          unlist(tp_specific_genes[c("2vs0", "3vs0", "4vs0")]))

# 4. Late response genes (TP4 specific)
late_response <- setdiff(tp_specific_genes[["4vs0"]], 
                         unlist(tp_specific_genes[c("1vs0", "2vs0", "3vs0")]))

# 5. Module-specific gene sets
module_genes <- split(names(net$colors), net$colors)

cat("\nGene set summary:\n")
cat("Broad core genes:", length(broad_core), "\n")
cat("Narrow core genes:", length(narrow_core), "\n")
cat("Early response genes:", length(early_response), "\n")
cat("Late response genes:", length(late_response), "\n")
cat("WGCNA modules:", length(module_genes), "\n")











## Temporal Expression Patterns
cat("\nAnalyzing temporal expression patterns...\n")

# Fix the sample name mismatch
# Option 1: Change targets rownames to match logCPM
targets_fixed <- targets
rownames(targets_fixed) <- targets$Label

# Verify the fix
cat("Fixed sample names match:", all(rownames(targets_fixed) %in% colnames(logCPM)), "\n")

# Now use the fixed targets for plotting
plot_temporal_pattern_simple <- function(genes, title, filename, targets_data = targets_fixed) {
  if (length(genes) == 0) {
    cat("No genes to plot for:", title, "\n")
    return(NULL)
  }
  
  # Filter genes that exist in logCPM
  genes <- intersect(genes, rownames(logCPM))
  if (length(genes) == 0) {
    cat("No valid genes found in expression data for:", title, "\n")
    return(NULL)
  }
  
  cat("Processing", title, "with", length(genes), "genes\n")
  
  # Limit genes if too many
  if (length(genes) > 500) {
    cat("  Selecting 500 random genes from", length(genes), "total\n")
    set.seed(123)
    genes <- sample(genes, 500)
  }
  
  # Get expression data
  expr_data <- logCPM[genes, , drop = FALSE]
  
  # Create plot
  pdf(filename, width = 10, height = 8)
  
  # Calculate mean expression for each condition
  plot_data <- data.frame()
  
  for (leaf in c("Broad", "Narrow")) {
    for (tp in paste0("TP", 0:4)) {
      # Find samples for this condition
      condition_samples <- rownames(targets_data)[
        targets_data$Leaf_type == leaf & targets_data$Timepoint == tp
      ]
      
      if (length(condition_samples) > 0) {
        # Get expression for these samples
        condition_expr <- expr_data[, condition_samples, drop = FALSE]
        
        # Calculate mean across samples for each gene
        gene_means <- rowMeans(condition_expr, na.rm = TRUE)
        
        # Calculate overall mean and SE
        plot_data <- rbind(plot_data, 
                           data.frame(
                             Timepoint = as.numeric(sub("TP", "", tp)),
                             Expression = mean(gene_means, na.rm = TRUE),
                             SE = sd(gene_means, na.rm = TRUE) / sqrt(length(gene_means)),
                             LeafType = leaf,
                             n_samples = length(condition_samples),
                             n_genes = length(gene_means)
                           ))
      }
    }
  }
  
  # Create line plot
  library(ggplot2)
  p <- ggplot(plot_data, aes(x = Timepoint, y = Expression, color = LeafType)) +
    geom_line(size = 2) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = Expression - SE, ymax = Expression + SE), 
                  width = 0.2, size = 1) +
    scale_color_manual(values = c("Broad" = "#E69F00", "Narrow" = "#56B4E9")) +
    labs(title = paste(title, "-", length(genes), "genes"),
         x = "Timepoint",
         y = "Mean Expression (logCPM)") +
    theme_bw() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold")) +
    scale_x_continuous(breaks = 0:4)
  
  print(p)
  dev.off()
  cat("  Created line plot:", filename, "\n")
  
  # Create heatmap version
  create_temporal_heatmap(genes, title, expr_data, targets_data, filename)
}

# Helper function for heatmap
create_temporal_heatmap <- function(genes, title, expr_data, targets_data, base_filename) {
  heatmap_file <- sub(".pdf", "_heatmap.pdf", base_filename)
  pdf(heatmap_file, width = 10, height = 12)
  
  # Calculate mean expression matrix
  mean_matrix <- matrix(NA, nrow = length(genes), ncol = 10)
  col_names <- c()
  col_idx <- 1
  
  for (leaf in c("Broad", "Narrow")) {
    for (tp in paste0("TP", 0:4)) {
      condition_samples <- rownames(targets_data)[
        targets_data$Leaf_type == leaf & targets_data$Timepoint == tp
      ]
      
      if (length(condition_samples) > 0) {
        mean_matrix[, col_idx] <- rowMeans(expr_data[, condition_samples, drop = FALSE], 
                                           na.rm = TRUE)
        col_names <- c(col_names, paste(substr(leaf, 1, 1), sub("TP", "", tp), sep = ""))
        col_idx <- col_idx + 1
      }
    }
  }
  
  # Clean up matrix
  mean_matrix <- mean_matrix[, 1:(col_idx-1), drop = FALSE]
  colnames(mean_matrix) <- col_names
  
  # Scale by row
  mean_matrix_scaled <- t(scale(t(mean_matrix)))
  
  # Remove rows with NaN (genes with no variance)
  valid_rows <- !apply(mean_matrix_scaled, 1, function(x) any(is.nan(x)))
  mean_matrix_scaled <- mean_matrix_scaled[valid_rows, , drop = FALSE]
  
  if (nrow(mean_matrix_scaled) > 100) {
    # For large gene sets, cluster and show top variable genes
    gene_vars <- apply(mean_matrix_scaled, 1, var, na.rm = TRUE)
    top_var_genes <- order(gene_vars, decreasing = TRUE)[1:100]
    mean_matrix_scaled <- mean_matrix_scaled[top_var_genes, ]
    subtitle <- paste("Top 100 most variable genes of", sum(valid_rows))
  } else {
    subtitle <- paste(nrow(mean_matrix_scaled), "genes")
  }
  
  # Define column colors
  col_colors <- c(rep("#E69F00", 5), rep("#56B4E9", 5))  # Broad=orange, Narrow=blue
  
  heatmap.2(mean_matrix_scaled,
            trace = "none",
            density.info = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = paste(title, "\n", subtitle),
            labRow = FALSE,
            ColSideColors = col_colors,
            margins = c(8, 5),
            cexCol = 1.2,
            dendrogram = "row",
            Colv = FALSE,
            key.title = "Z-score",
            key.xlab = "Expression",
            keysize = 1.2)
  
  # Add legend
  legend("topright", 
         legend = c("Broad", "Narrow"),
         fill = c("#E69F00", "#56B4E9"),
         bty = "n",
         cex = 0.8)
  
  dev.off()
  cat("  Created heatmap:", heatmap_file, "\n")
}

# Now create the plots
plot_temporal_pattern_simple(broad_core, "Broad Core Genes", 
                             "results/figures/gene_sets/temporal_broad_core.pdf")

plot_temporal_pattern_simple(narrow_core, "Narrow Core Genes", 
                             "results/figures/gene_sets/temporal_narrow_core.pdf")

plot_temporal_pattern_simple(early_response, "Early Response Genes", 
                             "results/figures/gene_sets/temporal_early_response.pdf")

plot_temporal_pattern_simple(late_response, "Late Response Genes", 
                             "results/figures/gene_sets/temporal_late_response.pdf")

# Also create plots for some interesting WGCNA modules
if (exists("module_genes")) {
  # Plot top 3 modules by size
  module_sizes <- sapply(module_genes, length)
  top_modules <- names(sort(module_sizes, decreasing = TRUE))[2:4]  # Skip grey (0)
  
  for (mod in top_modules) {
    if (mod != "0") {
      mod_color <- labels2colors(as.numeric(mod))
      plot_temporal_pattern_simple(module_genes[[mod]], 
                                   paste("Module", mod, "(", mod_color, ")"), 
                                   paste0("results/figures/gene_sets/temporal_module_", mod, ".pdf"))
    }
  }
}







## Heatmaps of Key Gene Sets
cat("\nCreating gene set heatmaps...\n")

# Function to create gene set heatmap
create_geneset_heatmap <- function(genes, title, filename, max_genes = 100) {
  if (length(genes) == 0) return(NULL)
  
  # Limit number of genes if too many
  if (length(genes) > max_genes) {
    # Select most variable genes
    gene_vars <- apply(logCPM[genes, ], 1, var)
    genes <- names(sort(gene_vars, decreasing = TRUE)[1:max_genes])
  }
  
  # Get and scale expression data
  expr_data <- logCPM[genes, ]
  expr_scaled <- t(scale(t(expr_data)))
  
  # Order samples
  sample_order <- order(targets$Timepoint, targets$Leaf_type, targets$Line)
  
  # Create annotation colors
  ann_colors <- list(
    Timepoint = rainbow(5),
    LeafType = c("Broad" = "darkgreen", "Narrow" = "purple"),
    Line = rainbow(4)
  )
  names(ann_colors$Timepoint) <- paste0("TP", 0:4)
  names(ann_colors$Line) <- unique(targets$Line)
  
  # Sample annotations
  sample_ann <- data.frame(
    Timepoint = targets$Timepoint[sample_order],
    LeafType = targets$Leaf_type[sample_order],
    Line = targets$Line[sample_order]
  )
  
  # Create heatmap
  pdf(filename, width = 12, height = 10)
  
  # Create color palette
  colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  
  # Plot heatmap
  heatmap.2(expr_scaled[, sample_order],
            col = colors,
            trace = "none",
            density.info = "none",
            dendrogram = "both",
            labRow = ifelse(length(genes) > 50, "", rownames(expr_scaled)),
            labCol = targets$Label[sample_order],
            ColSideColors = as.character(as.numeric(factor(sample_ann$LeafType))),
            margins = c(10, 8),
            main = title,
            key.title = "Z-score",
            key.xlab = "Expression",
            cexCol = 0.8)
  
  dev.off()
}

# Create heatmaps
create_geneset_heatmap(broad_core, "Broad Leaf Core Genes", 
                       "results/figures/gene_sets/heatmap_broad_core.pdf")
create_geneset_heatmap(narrow_core, "Narrow Leaf Core Genes", 
                       "results/figures/gene_sets/heatmap_narrow_core.pdf")

# Combine broad and narrow core genes
leaf_shape_genes <- union(broad_core, narrow_core)
create_geneset_heatmap(leaf_shape_genes, "All Leaf Shape Core Genes", 
                       "results/figures/gene_sets/heatmap_leaf_shape_core.pdf")









## Gene Set Overlap Analysis
cat("\nAnalyzing gene set overlaps...\n")

# Create upset plot or venn diagram
library(VennDiagram)

# Overlap between leaf types
pdf("results/figures/gene_sets/core_genes_overlap.pdf", width = 8, height = 8)
venn.plot <- draw.pairwise.venn(
  area1 = length(broad_core),
  area2 = length(narrow_core),
  cross.area = length(intersect(broad_core, narrow_core)),
  category = c("Broad Core", "Narrow Core"),
  fill = c("darkgreen", "purple"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5
)
dev.off()



## Module Enrichment in Gene Sets
cat("\nTesting module enrichment in gene sets...\n")

# Function to test enrichment
test_module_enrichment <- function(gene_set, module_genes) {
  enrichment_results <- data.frame(
    Module = integer(),
    GeneSet_Size = integer(),
    Module_Size = integer(),
    Overlap = integer(),
    P_value = numeric(),
    Odds_Ratio = numeric()
  )
  
  all_genes <- unique(unlist(module_genes))
  
  for (mod in names(module_genes)) {
    if (mod == "0") next  # Skip grey module
    
    # Create contingency table
    in_set_in_mod <- length(intersect(gene_set, module_genes[[mod]]))
    in_set_not_mod <- length(setdiff(gene_set, module_genes[[mod]]))
    not_set_in_mod <- length(setdiff(module_genes[[mod]], gene_set))
    not_set_not_mod <- length(all_genes) - in_set_in_mod - in_set_not_mod - not_set_in_mod
    
    # Fisher's exact test
    cont_table <- matrix(c(in_set_in_mod, in_set_not_mod,
                           not_set_in_mod, not_set_not_mod), 
                         nrow = 2)
    
    test_result <- fisher.test(cont_table)
    
    enrichment_results <- rbind(enrichment_results, data.frame(
      Module = as.integer(mod),
      GeneSet_Size = length(gene_set),
      Module_Size = length(module_genes[[mod]]),
      Overlap = in_set_in_mod,
      P_value = test_result$p.value,
      Odds_Ratio = test_result$estimate
    ))
  }
  
  enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
  return(enrichment_results)
}

# Test enrichments
broad_enrichment <- test_module_enrichment(broad_core, module_genes)
narrow_enrichment <- test_module_enrichment(narrow_core, module_genes)

# Add gene set labels
broad_enrichment$GeneSet <- "Broad_Core"
narrow_enrichment$GeneSet <- "Narrow_Core"

# Combine and save
all_enrichments <- rbind(broad_enrichment, narrow_enrichment)
write.csv(all_enrichments, 
          file = "results/gene_sets/module_geneset_enrichments.csv", 
          row.names = FALSE)

## Save Gene Sets
cat("\nSaving gene sets...\n")

gene_sets <- list(
  broad_core = broad_core,
  narrow_core = narrow_core,
  leaf_shape_core = leaf_shape_genes,
  early_response = early_response,
  late_response = late_response,
  timepoint_specific = tp_specific_genes
)

# Save as RData
save(gene_sets, file = "results/gene_sets/gene_sets.RData")

# Save as text files
for (set_name in names(gene_sets)) {
  if (is.list(gene_sets[[set_name]])) {
    # For nested lists
    for (sub_name in names(gene_sets[[set_name]])) {
      write.table(gene_sets[[set_name]][[sub_name]], 
                  file = paste0("results/gene_sets/", set_name, "_", sub_name, ".txt"),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
  } else {
    write.table(gene_sets[[set_name]], 
                file = paste0("results/gene_sets/", set_name, ".txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

# Save checkpoint
save_checkpoint("results/checkpoints/11_gene_sets_complete.RData")

cat("\n=== Gene set analysis complete ===\n")
cat("Gene sets identified and saved\n")
cat("Module enrichments calculated\n")
cat("Results saved to: results/gene_sets/\n")

## ===== 12_integration.R =====
#-------------------------------------------------------#
#     12. Results Integration                           #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load all analysis results
cat("Loading all analysis results...\n")
load_checkpoint("results/checkpoints/06_de_results.RData")
load_checkpoint("results/checkpoints/09_wgcna_analysis.RData")
load_checkpoint("results/checkpoints/11_gene_sets_complete.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")

# Create output directories
create_dir("results/integrated")
create_dir("results/figures/integrated")

## Integrate DE and WGCNA Results
cat("\nIntegrating differential expression and WGCNA results...\n")

# Create master gene table
master_gene_info <- data.frame(
  GeneID = rownames(logCPM),
  stringsAsFactors = FALSE
)

# Add average expression
master_gene_info$AvgExpression <- rowMeans(logCPM)

# Add coefficient of variation
master_gene_info$CV <- apply(logCPM, 1, function(x) sd(x)/mean(x))

# Add WGCNA module assignment
if (exists("genes_for_wgcna")) {
  master_gene_info$WGCNA_Module <- NA
  master_gene_info$WGCNA_Module[match(names(net$colors), master_gene_info$GeneID)] <- net$colors
  
  # Add module color
  master_gene_info$Module_Color <- NA
  master_gene_info$Module_Color[!is.na(master_gene_info$WGCNA_Module)] <- 
    labels2colors(master_gene_info$WGCNA_Module[!is.na(master_gene_info$WGCNA_Module)])
}

# Add DE summary statistics
de_summary_stats <- data.frame(
  GeneID = rownames(edgeR_tr_coded),
  N_DE_Contrasts = rowSums(edgeR_tr_coded != 0),
  N_Up = rowSums(edgeR_tr_coded > 0),
  N_Down = rowSums(edgeR_tr_coded < 0)
)

master_gene_info <- merge(master_gene_info, de_summary_stats, by = "GeneID", all.x = TRUE)

# Add gene set membership
master_gene_info$Is_Broad_Core <- master_gene_info$GeneID %in% gene_sets$broad_core
master_gene_info$Is_Narrow_Core <- master_gene_info$GeneID %in% gene_sets$narrow_core
master_gene_info$Is_Early_Response <- master_gene_info$GeneID %in% gene_sets$early_response
master_gene_info$Is_Late_Response <- master_gene_info$GeneID %in% gene_sets$late_response

# Add maximum fold changes
max_fc_data <- data.frame(
  GeneID = rownames(logCPM),
  Max_LogFC_Broad = NA,
  Max_LogFC_Narrow = NA,
  Max_LogFC_Overall = NA
)

# Calculate max fold changes
broad_contrasts <- grep("B_462A|B_LD11", names(all_results_treat), value = TRUE)
narrow_contrasts <- grep("N_713B|N_745", names(all_results_treat), value = TRUE)

for (gene in max_fc_data$GeneID) {
  # Broad max FC
  broad_fcs <- sapply(broad_contrasts, function(x) {
    if (gene %in% rownames(all_results_treat[[x]])) {
      return(all_results_treat[[x]][gene, "logFC"])
    } else {
      return(0)
    }
  })
  max_fc_data[max_fc_data$GeneID == gene, "Max_LogFC_Broad"] <- max(abs(broad_fcs))
  
  # Narrow max FC
  narrow_fcs <- sapply(narrow_contrasts, function(x) {
    if (gene %in% rownames(all_results_treat[[x]])) {
      return(all_results_treat[[x]][gene, "logFC"])
    } else {
      return(0)
    }
  })
  max_fc_data[max_fc_data$GeneID == gene, "Max_LogFC_Narrow"] <- max(abs(narrow_fcs))
  
  # Overall max FC
  max_fc_data[max_fc_data$GeneID == gene, "Max_LogFC_Overall"] <- max(abs(c(broad_fcs, narrow_fcs)))
}

master_gene_info <- merge(master_gene_info, max_fc_data, by = "GeneID")

## Module-DE Integration
cat("\nIntegrating module information with DE results...\n")

# For each module, summarize DE patterns
module_de_summary <- data.frame()

for (mod in unique(na.omit(master_gene_info$WGCNA_Module))) {
  if (mod == 0) next
  
  module_genes <- master_gene_info$GeneID[master_gene_info$WGCNA_Module == mod & !is.na(master_gene_info$WGCNA_Module)]
  
  # Calculate DE enrichment for each contrast
  for (contrast_name in colnames(edgeR_tr_coded)) {
    de_genes <- rownames(edgeR_tr_coded)[edgeR_tr_coded[, contrast_name] != 0]
    
    # Fisher's exact test
    in_mod_de <- sum(module_genes %in% de_genes)
    in_mod_not_de <- sum(!module_genes %in% de_genes)
    not_mod_de <- sum(!rownames(edgeR_tr_coded) %in% module_genes & rownames(edgeR_tr_coded) %in% de_genes)
    not_mod_not_de <- sum(!rownames(edgeR_tr_coded) %in% module_genes & !rownames(edgeR_tr_coded) %in% de_genes)
    
    test_result <- fisher.test(matrix(c(in_mod_de, in_mod_not_de, not_mod_de, not_mod_not_de), nrow = 2))
    
    module_de_summary <- rbind(module_de_summary, data.frame(
      Module = mod,
      Contrast = contrast_name,
      Module_Size = length(module_genes),
      DE_in_Module = in_mod_de,
      Percent_DE = round(in_mod_de / length(module_genes) * 100, 2),
      P_value = test_result$p.value,
      Odds_Ratio = test_result$estimate
    ))
  }
}

module_de_summary$FDR <- p.adjust(module_de_summary$P_value, method = "BH")

# Save module-DE summary
write.csv(module_de_summary, 
          file = "results/integrated/module_DE_enrichment.csv", 
          row.names = FALSE)

## Create Integration Visualizations
cat("\nCreating integration visualizations...\n")

# 1. Module-trait-DE heatmap
significant_modules <- unique(module_de_summary$Module[module_de_summary$FDR < 0.05])

if (length(significant_modules) > 0) {
  # Prepare data for heatmap
  module_trait_de <- data.frame(
    Module = paste0("ME", 0:max(net$colors))
  )
  
  # Add trait correlations
  for (trait in colnames(moduleTraitCor)) {
    module_trait_de[, paste0(trait, "_cor")] <- moduleTraitCor[, trait]
    module_trait_de[, paste0(trait, "_pval")] <- moduleTraitPvalue[, trait]
  }
  
  # Add DE enrichments (simplified - just count significant contrasts)
  de_counts <- aggregate(FDR ~ Module, 
                         data = module_de_summary[module_de_summary$FDR < 0.05, ],
                         FUN = length)
  
  module_trait_de$N_DE_Enriched <- 0
  module_trait_de$N_DE_Enriched[match(de_counts$Module, 0:max(net$colors))] <- de_counts$FDR
  
  # Create combined heatmap
  pdf("results/figures/integrated/module_trait_DE_summary.pdf", width = 12, height = 10)
  # Implementation depends on specific visualization needs
  dev.off()
}

# 2. Gene network for top modules
cat("\nCreating gene networks for top modules...\n")

# Identify top modules based on trait correlations
top_modules <- which(apply(abs(moduleTraitCor), 1, max) > 0.7)

for (i in top_modules[1:min(3, length(top_modules))]) {
  module_num <- as.numeric(sub("ME", "", rownames(moduleTraitCor)[i]))
  if (module_num == 0) next
  
  # Get module hub genes
  module_genes <- master_gene_info$GeneID[master_gene_info$WGCNA_Module == module_num & !is.na(master_gene_info$WGCNA_Module)]
  
  # Get top 20 hub genes
  if (paste0("kME", module_num) %in% colnames(kMEs)) {
    hub_scores <- kMEs[module_genes, paste0("kME", module_num)]
    top_hubs <- names(sort(hub_scores, decreasing = TRUE)[1:min(20, length(hub_scores))])
    
    # Create hub gene summary
    hub_summary <- master_gene_info[master_gene_info$GeneID %in% top_hubs, 
                                    c("GeneID", "AvgExpression", "N_DE_Contrasts", 
                                      "Is_Broad_Core", "Is_Narrow_Core")]
    
    write.csv(hub_summary, 
              file = paste0("results/integrated/module_", module_num, "_hub_genes.csv"),
              row.names = FALSE)
  }
}

## Identify Key Regulatory Genes
cat("\nIdentifying key regulatory genes...\n")

# Criteria for key genes:
# 1. High module membership (hub genes)
# 2. Significant in multiple DE contrasts
# 3. Member of core gene sets
# 4. High expression variance

key_genes <- master_gene_info[
  master_gene_info$N_DE_Contrasts >= 5 &  # DE in at least 5 contrasts
    (master_gene_info$Is_Broad_Core | master_gene_info$Is_Narrow_Core) &  # Core gene
    master_gene_info$CV > quantile(master_gene_info$CV, 0.75, na.rm = TRUE) &  # High variance
    !is.na(master_gene_info$WGCNA_Module) & master_gene_info$WGCNA_Module != 0,  # In a module
]

# Add hub gene status
key_genes$Is_Hub <- FALSE
for (i in 1:nrow(key_genes)) {
  mod <- key_genes$WGCNA_Module[i]
  if (paste0("kME", mod) %in% colnames(kMEs)) {
    kme_value <- kMEs[key_genes$GeneID[i], paste0("kME", mod)]
    key_genes$Is_Hub[i] <- kme_value > 0.8
  }
}

# Sort by importance
key_genes <- key_genes[order(key_genes$N_DE_Contrasts, decreasing = TRUE), ]

# Save key genes
write.csv(key_genes, file = "results/integrated/key_regulatory_genes.csv", row.names = FALSE)

## Create Summary Report
cat("\nCreating integration summary report...\n")

summary_stats <- list(
  Total_Genes = nrow(master_gene_info),
  Genes_in_WGCNA = sum(!is.na(master_gene_info$WGCNA_Module)),
  N_Modules = length(unique(na.omit(master_gene_info$WGCNA_Module))) - 1,  # Exclude grey
  DE_Genes_Any_Contrast = sum(master_gene_info$N_DE_Contrasts > 0, na.rm = TRUE),
  Broad_Core_Genes = sum(master_gene_info$Is_Broad_Core),
  Narrow_Core_Genes = sum(master_gene_info$Is_Narrow_Core),
  Key_Regulatory_Genes = nrow(key_genes),
  Modules_Trait_Associated = sum(apply(moduleTraitPvalue, 1, min) < 0.05),
  Modules_DE_Enriched = length(unique(module_de_summary$Module[module_de_summary$FDR < 0.05]))
)

# Save summary
saveRDS(summary_stats, file = "results/integrated/analysis_summary_stats.rds")

# Print summary
cat("\n=== INTEGRATION SUMMARY ===\n")
for (stat_name in names(summary_stats)) {
  cat(stat_name, ":", summary_stats[[stat_name]], "\n")
}

## Save Master Gene Table
cat("\nSaving master gene table...\n")

# Order columns logically
col_order <- c("GeneID", "AvgExpression", "CV", "WGCNA_Module", "Module_Color",
               "N_DE_Contrasts", "N_Up", "N_Down", "Max_LogFC_Overall",
               "Max_LogFC_Broad", "Max_LogFC_Narrow",
               "Is_Broad_Core", "Is_Narrow_Core", "Is_Early_Response", "Is_Late_Response")

master_gene_info <- master_gene_info[, col_order]

# Save as CSV
write.csv(master_gene_info, 
          file = "results/integrated/master_gene_table.csv", 
          row.names = FALSE)

# Save checkpoint
save(master_gene_info, module_de_summary, key_genes, summary_stats,
     file = "results/checkpoints/12_integration_complete.RData")

cat("\n=== Results integration complete ===\n")
cat("Master gene table created with", nrow(master_gene_info), "genes\n")
cat("Key regulatory genes identified:", nrow(key_genes), "\n")
cat("Results saved to: results/integrated/\n")

## ===== 13_JAG1_analysis.R =====
#-------------------------------------------------------#
#     13. JAG1-Specific Analysis                        #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load checkpoints
cat("Loading data for JAG1 analysis...\n")
load_checkpoint("results/checkpoints/12_integration_complete.RData")
load_checkpoint("results/checkpoints/04_preprocessed_data.RData")
load_checkpoint("results/checkpoints/06_de_results.RData")

# Create output directories
create_dir("results/JAG1_analysis")
create_dir("results/figures/JAG1_analysis")

# Define JAG1 gene ID
JAG1_ID <- "Glyma.20G116200"
JAG2_ID <- "Glyma.10G273800"  # If you want to analyze JAG2 as well

# Check if JAG1 is in dataset
if (!JAG1_ID %in% rownames(logCPM)) {
  stop("JAG1 gene not found in dataset!")
}

## JAG1 Expression Profile
cat("\nAnalyzing JAG1 expression profile...\n")

# Extract JAG1 expression
JAG1_expr <- logCPM[JAG1_ID, ]

# Create expression plot by group
jag1_data <- data.frame(
  Expression = JAG1_expr,
  Sample = names(JAG1_expr),
  Line = targets$Line,
  Timepoint = targets$Timepoint,
  LeafType = targets$Leaf_type,
  Group = targets$group
)

# Plot JAG1 expression
library(ggplot2)
p1 <- ggplot(jag1_data, aes(x = Timepoint, y = Expression, color = Line, group = Line)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ LeafType) +
  theme_minimal() +
  labs(title = "JAG1 Expression Across Development",
       y = "Expression (logCPM)",
       x = "Timepoint") +
  theme(legend.position = "bottom")

ggsave("results/figures/JAG1_analysis/JAG1_expression_profile.pdf", 
       plot = p1, width = 10, height = 6)

# Statistical analysis of JAG1 expression
jag1_stats <- aggregate(Expression ~ LeafType + Timepoint, data = jag1_data, 
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))

write.csv(jag1_stats, file = "results/JAG1_analysis/JAG1_expression_stats.csv")

## JAG1 Co-expression Network
cat("\nBuilding JAG1 co-expression network...\n")

# Focus on TP0 (meristem) where JAG1 is expected to be active
TP0_samples <- which(targets$Timepoint == "TP0")

# Calculate correlations with all genes at TP0
cor_with_JAG1_TP0 <- cor(JAG1_expr[TP0_samples], 
                         t(logCPM[, TP0_samples]), 
                         method = "spearman")

# Get highly correlated genes
cor_threshold <- 0.7
JAG1_network_genes <- names(which(abs(cor_with_JAG1_TP0) > cor_threshold))
JAG1_network_genes <- setdiff(JAG1_network_genes, JAG1_ID)  # Remove JAG1 itself

cat("Genes correlated with JAG1 at TP0 (|r| >", cor_threshold, "):", 
    length(JAG1_network_genes), "\n")

# Separate positive and negative correlations
JAG1_pos_cor <- names(which(cor_with_JAG1_TP0 > cor_threshold))
JAG1_neg_cor <- names(which(cor_with_JAG1_TP0 < -cor_threshold))

cat("Positively correlated:", length(JAG1_pos_cor), "\n")
cat("Negatively correlated:", length(JAG1_neg_cor), "\n")

# Create network data
jag1_network_data <- data.frame(
  Gene = JAG1_network_genes,
  Correlation = cor_with_JAG1_TP0[JAG1_network_genes],
  Direction = ifelse(cor_with_JAG1_TP0[JAG1_network_genes] > 0, "Positive", "Negative")
)

# Add gene information
jag1_network_data <- merge(jag1_network_data, 
                           master_gene_info[, c("GeneID", "AvgExpression", "N_DE_Contrasts", "WGCNA_Module")],
                           by.x = "Gene", by.y = "GeneID", all.x = TRUE)

# Save network genes
write.csv(jag1_network_data, 
          file = "results/JAG1_analysis/JAG1_network_genes_TP0.csv", 
          row.names = FALSE)

## JAG1 Network Dynamics
cat("\nAnalyzing JAG1 network dynamics across timepoints...\n")

# Calculate correlations at each timepoint
timepoint_cors <- list()

for (tp in unique(targets$Timepoint)) {
  tp_samples <- which(targets$Timepoint == tp)
  
  # Calculate correlations
  cors <- cor(JAG1_expr[tp_samples], 
              t(logCPM[JAG1_network_genes, tp_samples]), 
              method = "spearman")
  
  timepoint_cors[[tp]] <- as.vector(cors)
}

# Create heatmap of correlation dynamics
cor_dynamics <- do.call(cbind, timepoint_cors)
rownames(cor_dynamics) <- JAG1_network_genes
colnames(cor_dynamics) <- paste0("TP", 0:4)

# Plot correlation dynamics
pdf("results/figures/JAG1_analysis/JAG1_network_correlation_dynamics.pdf", 
    width = 8, height = 10)

heatmap.2(cor_dynamics,
          col = colorRampPalette(c("blue", "white", "red"))(100),
          trace = "none",
          density.info = "none",
          margins = c(8, 10),
          main = "JAG1 Network Correlation Dynamics",
          xlab = "Timepoint",
          key.title = "Correlation",
          cexCol = 1.2)

dev.off()



## Differential Network Analysis
cat("\nPerforming differential network analysis...\n")

# Compare JAG1 network between leaf types
broad_samples <- which(targets$Leaf_type == "Broad")
narrow_samples <- which(targets$Leaf_type == "Narrow")

# Calculate leaf-type specific correlations
broad_cors <- cor(JAG1_expr[broad_samples], 
                  t(logCPM[, broad_samples]), 
                  method = "spearman")

narrow_cors <- cor(JAG1_expr[narrow_samples], 
                   t(logCPM[, narrow_samples]), 
                   method = "spearman")

# Find differential relationships
cor_diff <- broad_cors - narrow_cors
names(cor_diff) <- rownames(logCPM)

# Get genes with different relationships
diff_threshold <- 0.5
diff_network_genes <- names(which(abs(cor_diff) > diff_threshold))

cat("Genes with differential correlation to JAG1 (|diff| >", diff_threshold, "):", 
    length(diff_network_genes), "\n")

# Create differential network summary
diff_network_data <- data.frame(
  Gene = diff_network_genes,
  Broad_Correlation = broad_cors[diff_network_genes],
  Narrow_Correlation = narrow_cors[diff_network_genes],
  Correlation_Difference = cor_diff[diff_network_genes],
  Switch_Type = ifelse(
    sign(broad_cors[diff_network_genes]) != sign(narrow_cors[diff_network_genes]),
    "Sign_Switch",
    "Magnitude_Change"
  )
)

# Add DE information
for (gene in diff_network_data$Gene) {
  # Check if DE between leaf types
  de_contrasts <- grep("vs.*TP0", names(all_results_treat), value = TRUE)
  
  de_status <- sapply(de_contrasts, function(x) {
    if (gene %in% rownames(all_results_treat[[x]])) {
      return(all_results_treat[[x]][gene, "FDR"] < 0.05)
    }
    return(FALSE)
  })
  
  diff_network_data[diff_network_data$Gene == gene, "DE_LeafType"] <- any(de_status)
}

write.csv(diff_network_data, 
          file = "results/JAG1_analysis/JAG1_differential_network.csv", 
          row.names = FALSE)

## JAG1 Target Prediction
cat("\nPredicting JAG1 target genes...\n")

# Criteria for JAG1 targets:
# 1. High correlation at TP0
# 2. Expression changes after TP0
# 3. Different expression pattern in leaf types

jag1_targets <- data.frame(
  Gene = character(),
  TP0_Correlation = numeric(),
  Expression_Change_TP1 = numeric(),
  LeafType_Difference = numeric(),
  Score = numeric(),
  stringsAsFactors = FALSE
)

for (gene in JAG1_network_genes) {
  # TP0 correlation
  tp0_cor <- cor_with_JAG1_TP0[gene]
  
  # Expression change from TP0 to TP1
  tp0_expr <- mean(logCPM[gene, targets$Timepoint == "TP0"])
  tp1_expr <- mean(logCPM[gene, targets$Timepoint == "TP1"])
  expr_change <- tp1_expr - tp0_expr
  
  # Leaf type difference
  broad_expr <- mean(logCPM[gene, targets$Leaf_type == "Broad"])
  narrow_expr <- mean(logCPM[gene, targets$Leaf_type == "Narrow"])
  leaf_diff <- broad_expr - narrow_expr
  
  # Calculate score
  score <- abs(tp0_cor) * abs(expr_change) * abs(leaf_diff)
  
  jag1_targets <- rbind(jag1_targets, data.frame(
    Gene = gene,
    TP0_Correlation = tp0_cor,
    Expression_Change_TP1 = expr_change,
    LeafType_Difference = leaf_diff,
    Score = score
  ))
}

# Sort by score
jag1_targets <- jag1_targets[order(jag1_targets$Score, decreasing = TRUE), ]

# Add functional annotation if available
jag1_targets <- merge(jag1_targets, 
                      master_gene_info[, c("GeneID", "WGCNA_Module", "N_DE_Contrasts")],
                      by.x = "Gene", by.y = "GeneID", all.x = TRUE)

# Save top targets
write.csv(head(jag1_targets, 100), 
          file = "results/JAG1_analysis/JAG1_predicted_targets_top100.csv", 
          row.names = FALSE)

## JAG1 Module Analysis
cat("\nAnalyzing JAG1's WGCNA module...\n")

# Find JAG1's module
jag1_module_info <- master_gene_info[master_gene_info$GeneID == JAG1_ID, ]

if (!is.na(jag1_module_info$WGCNA_Module) && jag1_module_info$WGCNA_Module != 0) {
  jag1_module <- jag1_module_info$WGCNA_Module
  cat("JAG1 is in module:", jag1_module, "\n")
  
  # Get all genes in JAG1's module
  module_genes <- master_gene_info$GeneID[
    master_gene_info$WGCNA_Module == jag1_module & !is.na(master_gene_info$WGCNA_Module)
  ]
  
  cat("Total genes in JAG1's module:", length(module_genes), "\n")
  
  # Check enrichment of JAG1 network in its module
  network_in_module <- sum(JAG1_network_genes %in% module_genes)
  
  cat("JAG1 network genes in same module:", network_in_module, 
      "(", round(network_in_module/length(JAG1_network_genes)*100, 1), "%)\n")
  
  # Module eigengene correlation with JAG1
  if (paste0("ME", jag1_module) %in% colnames(net$MEs)) {
    me_jag1_cor <- cor(JAG1_expr, net$MEs[, paste0("ME", jag1_module)])
    cat("JAG1 correlation with module eigengene:", round(me_jag1_cor, 3), "\n")
  }
}

## Expression Trajectory Analysis
cat("\nAnalyzing JAG1-dependent expression trajectories...\n")

# Define JAG1-high and JAG1-low samples
jag1_median <- median(JAG1_expr)
jag1_high_samples <- which(JAG1_expr > jag1_median)
jag1_low_samples <- which(JAG1_expr <= jag1_median)

# For top JAG1 targets, compare trajectories
top_targets <- head(jag1_targets$Gene, 20)

# Create trajectory plots
pdf("results/figures/JAG1_analysis/JAG1_target_trajectories.pdf", 
    width = 12, height = 10)
par(mfrow = c(4, 5))

for (target in top_targets) {
  # Get expression data
  target_expr <- logCPM[target, ]
  
  # Calculate means by timepoint and JAG1 status
  trajectory_data <- aggregate(target_expr, 
                               by = list(Timepoint = targets$Timepoint,
                                         JAG1_Status = ifelse(1:length(target_expr) %in% jag1_high_samples, 
                                                              "JAG1_High", "JAG1_Low")),
                               FUN = mean)
  
  # Reshape for plotting
  trajectory_wide <- reshape(trajectory_data, 
                             idvar = "Timepoint", 
                             timevar = "JAG1_Status", 
                             direction = "wide")
  
  # Plot
  tp_numeric <- as.numeric(sub("TP", "", trajectory_wide$Timepoint))
  plot(tp_numeric, trajectory_wide$x.JAG1_High, 
       type = "b", col = "red", pch = 19,
       ylim = range(c(trajectory_wide$x.JAG1_High, trajectory_wide$x.JAG1_Low)),
       xlab = "Timepoint", ylab = "Expression",
       main = target, xaxt = "n")
  
  lines(tp_numeric, trajectory_wide$x.JAG1_Low, 
        type = "b", col = "blue", pch = 19)
  
  axis(1, at = tp_numeric, labels = trajectory_wide$Timepoint)
  legend("topright", legend = c("JAG1 High", "JAG1 Low"), 
         col = c("red", "blue"), lty = 1, pch = 19, cex = 0.8)
}

dev.off()

## JAG1 Regulatory Model
cat("\nBuilding JAG1 regulatory model...\n")

# Summarize JAG1 regulatory patterns
jag1_model <- list(
  JAG1_Expression = list(
    Peak_Timepoint = names(which.max(tapply(JAG1_expr, targets$Timepoint, mean))),
    Broad_vs_Narrow_FC = mean(JAG1_expr[targets$Leaf_type == "Broad"]) - 
      mean(JAG1_expr[targets$Leaf_type == "Narrow"]),
    CV = sd(JAG1_expr) / mean(JAG1_expr)
  ),
  
  Network_Summary = list(
    Total_Correlated_Genes = length(JAG1_network_genes),
    Positive_Targets = length(JAG1_pos_cor),
    Negative_Targets = length(JAG1_neg_cor),
    Differential_Network_Genes = nrow(diff_network_data)
  ),
  
  Top_Targets = head(jag1_targets[, c("Gene", "Score")], 10),
  
  Module_Info = if (!is.na(jag1_module_info$WGCNA_Module)) {
    list(
      Module = jag1_module_info$WGCNA_Module,
      Module_Size = sum(master_gene_info$WGCNA_Module == jag1_module_info$WGCNA_Module, na.rm = TRUE),
      Module_Trait_Association = moduleTraitCor[paste0("ME", jag1_module_info$WGCNA_Module), ]
    )
  } else {
    "JAG1 not in a WGCNA module"
  }
)

# Save model
saveRDS(jag1_model, file = "results/JAG1_analysis/JAG1_regulatory_model.rds")

## Create JAG1 Summary Report
cat("\nCreating JAG1 analysis summary report...\n")

# Create summary figure
pdf("results/figures/JAG1_analysis/JAG1_summary_figure.pdf", width = 14, height = 10)
layout(matrix(c(1,1,2,2,3,3,4,5,5,6,6,6), nrow = 2, byrow = TRUE))

# 1. JAG1 expression profile
barplot(tapply(JAG1_expr, targets$group, mean),
        las = 2, col = targets$col[!duplicated(targets$group)],
        main = "JAG1 Expression by Group",
        ylab = "Expression (logCPM)")

# 2. Network size by timepoint
network_sizes <- sapply(names(timepoint_cors), function(tp) {
  sum(abs(timepoint_cors[[tp]]) > 0.5)
})
barplot(network_sizes, names.arg = names(timepoint_cors),
        main = "JAG1 Network Size by Timepoint",
        ylab = "Number of Correlated Genes",
        col = "skyblue")

# 3. Top targets heatmap
top_target_expr <- logCPM[head(jag1_targets$Gene, 10), ]
heatmap(as.matrix(top_target_expr), 
        scale = "row",
        main = "Top 10 JAG1 Targets",
        margins = c(10, 8))

# 4. Correlation distribution
hist(cor_with_JAG1_TP0, breaks = 50,
     main = "JAG1 Correlation Distribution (TP0)",
     xlab = "Correlation", col = "lightgreen")
abline(v = c(-cor_threshold, cor_threshold), col = "red", lty = 2)

# 5. Differential network
plot(diff_network_data$Broad_Correlation, 
     diff_network_data$Narrow_Correlation,
     pch = 19, col = ifelse(diff_network_data$Switch_Type == "Sign_Switch", "red", "blue"),
     main = "Differential JAG1 Network",
     xlab = "Broad Leaf Correlation",
     ylab = "Narrow Leaf Correlation")
abline(0, 1, lty = 2)
abline(h = 0, v = 0, lty = 3, col = "gray")

# 6. Module-trait relationships if available
if (!is.na(jag1_module_info$WGCNA_Module) && jag1_module_info$WGCNA_Module != 0) {
  barplot(moduleTraitCor[paste0("ME", jag1_module_info$WGCNA_Module), ],
          las = 2,
          main = paste("JAG1 Module", jag1_module_info$WGCNA_Module, "Trait Correlations"),
          ylim = c(-1, 1))
  abline(h = 0)
}

dev.off()

# Save all JAG1 analysis results
jag1_results <- list(
  expression_profile = jag1_data,
  network_genes = jag1_network_data,
  differential_network = diff_network_data,
  predicted_targets = jag1_targets,
  regulatory_model = jag1_model
)

save(jag1_results, file = "results/checkpoints/13_JAG1_analysis_complete.RData")

cat("\n=== JAG1 analysis complete ===\n")
cat("JAG1 network genes identified:", length(JAG1_network_genes), "\n")
cat("Differential network genes:", nrow(diff_network_data), "\n")
cat("Top predicted targets saved\n")
cat("Results saved to: results/JAG1_analysis/\n")

## ===== 14_reporting.R =====
#-------------------------------------------------------#
#     14. Final Reporting and Export                    #
#     Soybean RNA-Seq Analysis                          #
#-------------------------------------------------------#

# Load all results
cat("Loading all analysis results for reporting...\n")
load_checkpoint("results/checkpoints/12_integration_complete.RData")
load_checkpoint("results/checkpoints/13_JAG1_analysis_complete.RData")
load_checkpoint("results/checkpoints/11_gene_sets_complete.RData")

# Create final output directories
create_dir("results/final_report")
create_dir("results/final_report/figures")
create_dir("results/final_report/tables")
create_dir("results/final_report/supplementary")

## Generate Publication Figures
cat("\nGenerating publication-ready figures...\n")

# Figure 1: Experimental Overview
pdf("results/final_report/figures/Figure1_experimental_overview.pdf", 
    width = 12, height = 10)
layout(matrix(c(1,1,2,3,4,4), nrow = 3, byrow = TRUE))

# A. Sample distribution
par(mar = c(5, 4, 4, 2))
plot_data <- table(targets$Leaf_type, targets$Timepoint)
barplot(plot_data, beside = TRUE, 
        col = c("darkgreen", "purple"),
        main = "A. Experimental Design",
        xlab = "Timepoint", ylab = "Number of Samples",
        legend.text = rownames(plot_data))

# B. PCA plot
load_checkpoint("results/checkpoints/03_batch_corrected.RData")
mds_data <- plotMDS(d_filt_combatseq, plot = FALSE)
plot(mds_data$x, mds_data$y,
     col = targets$col, pch = 19, cex = 1.5,
     xlab = "PC1", ylab = "PC2",
     main = "B. Sample Clustering (Post Batch Correction)")
legend("topright", legend = unique(targets$group), 
       col = unique(targets$col), pch = 19, cex = 0.8, ncol = 2)

# C. Batch effect correction
load("results/tables/batch_correction_results.RData")
barplot(c(batch_results$batch_variance_before, batch_results$batch_variance_after),
        names.arg = c("Before", "After"),
        col = c("red", "green"),
        main = "C. Batch Effect Reduction",
        ylab = "Proportion of Variance",
        ylim = c(0, 0.4))

# D. Gene filtering summary
pie_data <- c(
  Kept = batch_results$genes_filtered,
  Filtered = nrow(d_salmon) - batch_results$genes_filtered
)
pie(pie_data, 
    labels = paste(names(pie_data), "\n", pie_data, " genes"),
    col = c("lightblue", "lightgray"),
    main = "D. Gene Filtering")

dev.off()

# Figure 2: Differential Expression Overview
cat("Creating Figure 2: Differential Expression Overview...\n")

pdf("results/final_report/figures/Figure2_differential_expression.pdf", 
    width = 14, height = 10)
layout(matrix(c(1,2,3,4,5,5), nrow = 2, byrow = TRUE))

# A-D: Volcano plots for each line at TP4 vs TP0
contrasts_to_plot <- c("B_462A_4vs0", "N_713B_4vs0", "N_745_4vs0", "B_LD11_4vs0")
line_names <- c("PI532462A", "PI612713B", "PI547745", "LD112170")

for (i in 1:4) {
  results <- all_results_treat[[contrasts_to_plot[i]]]
  
  plot(results$logFC, -log10(results$PValue),
       pch = 19, col = rgb(0,0,0,0.3),
       xlab = "log2 Fold Change",
       ylab = "-log10(p-value)",
       main = paste(LETTERS[i], ". ", line_names[i], " (TP4 vs TP0)", sep = ""),
       xlim = c(-10, 10))
  
  # Highlight significant genes
  sig_genes <- results$FDR < 0.05 & abs(results$logFC) > log2(1.2)
  points(results$logFC[sig_genes], -log10(results$PValue[sig_genes]),
         pch = 19, col = ifelse(results$logFC[sig_genes] > 0, "red", "blue"))
  
  abline(v = c(-log2(1.2), log2(1.2)), lty = 2, col = "gray")
  abline(h = -log10(0.05), lty = 2, col = "gray")
  
  # Add counts
  n_up <- sum(results$FDR < 0.05 & results$logFC > log2(1.2))
  n_down <- sum(results$FDR < 0.05 & results$logFC < -log2(1.2))
  
  legend("topright", 
         legend = c(paste("Up:", n_up), paste("Down:", n_down)),
         text.col = c("red", "blue"),
         bty = "n")
}

# E: Summary heatmap of DEGs across all contrasts
de_summary_matrix <- as.matrix(summary(edgeR_tr_coded))
de_summary_matrix <- de_summary_matrix[c(1,3), ]  # Keep only Up and Down

# Create heatmap
image(1:ncol(de_summary_matrix), 1:2, t(de_summary_matrix),
      col = colorRampPalette(c("white", "darkred"))(100),
      xlab = "Contrast", ylab = "",
      main = "E. Differential Expression Summary",
      axes = FALSE)

axis(1, at = 1:ncol(de_summary_matrix), 
     labels = colnames(de_summary_matrix), las = 2, cex.axis = 0.7)
axis(2, at = 1:2, labels = c("Down", "Up"), las = 1)

# Add numbers
for (i in 1:ncol(de_summary_matrix)) {
  for (j in 1:2) {
    text(i, j, de_summary_matrix[j, i], cex = 0.6)
  }
}

dev.off()

# Figure 3: WGCNA Results
cat("Creating Figure 3: WGCNA Results...\n")

if (exists("net") && exists("moduleTraitCor")) {
  pdf("results/final_report/figures/Figure3_WGCNA_results.pdf", 
      width = 12, height = 10)
  
  # Module-trait relationships heatmap
  # (Simplified version - full version in WGCNA visualization script)
  
  # Select significant modules
  sig_modules <- which(apply(moduleTraitPvalue, 1, min) < 0.01)
  
  if (length(sig_modules) > 0) {
    heatmap(moduleTraitCor[sig_modules, ],
            scale = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            margins = c(10, 8),
            main = "Module-Trait Relationships")
  }
  
  dev.off()
}

# Figure 4: JAG1 Analysis Summary
cat("Creating Figure 4: JAG1 Analysis Summary...\n")

# (This would be a more polished version of the JAG1 summary figure from script 13)
# Copy the best elements from the JAG1 analysis

## Generate Summary Tables
cat("\nGenerating summary tables...\n")

# Table 1: Sample Information
write.csv(targets[, c("Sample", "Label", "Line", "Leaf_type", "Timepoint", "Batch")],
          file = "results/final_report/tables/Table1_sample_information.csv",
          row.names = FALSE)

# Table 2: Top DEGs
top_degs_list <- list()
for (contrast in c("B_462A_4vs0", "N_713B_4vs0", "N_745_4vs0", "B_LD11_4vs0")) {
  top_genes <- head(all_results_treat[[contrast]][order(all_results_treat[[contrast]]$FDR), ], 20)
  top_genes$Contrast <- contrast
  top_degs_list[[contrast]] <- top_genes[, c("GeneID", "logFC", "FDR", "Contrast")]
}

top_degs_combined <- do.call(rbind, top_degs_list)
write.csv(top_degs_combined,
          file = "results/final_report/tables/Table2_top_DEGs.csv",
          row.names = FALSE)

# Table 3: Key Regulatory Genes
if (exists("key_genes")) {
  write.csv(head(key_genes, 50),
            file = "results/final_report/tables/Table3_key_regulatory_genes.csv",
            row.names = FALSE)
}

# Table 4: JAG1 Network Summary
if (exists("jag1_results")) {
  write.csv(head(jag1_results$predicted_targets, 20),
            file = "results/final_report/tables/Table4_JAG1_top_targets.csv",
            row.names = FALSE)
}

## Generate Supplementary Data
cat("\nGenerating supplementary data files...\n")

# Supplementary Data 1: Complete DE results
save(all_results_treat, edgeR_tr_coded,
     file = "results/final_report/supplementary/SuppData1_complete_DE_results.RData")

# Supplementary Data 2: WGCNA results
if (exists("net")) {
  save(net, moduleTraitCor, moduleTraitPvalue, kMEs,
       file = "results/final_report/supplementary/SuppData2_WGCNA_results.RData")
}

# Supplementary Data 3: Gene sets
save(gene_sets,
     file = "results/final_report/supplementary/SuppData3_gene_sets.RData")

# Supplementary Data 4: Normalized expression matrix
write.csv(logCPM,
          file = "results/final_report/supplementary/SuppData4_normalized_expression.csv.gz")

## Generate Methods Text
cat("\nGenerating methods text...\n")

methods_text <- "
RNA-Seq Data Analysis Methods

Data Processing and Quality Control:
RNA-seq data from 60 soybean leaf samples (2 broad-leaved lines, 2 narrow-leaved lines, 
5 timepoints, 3 replicates) were processed using the following pipeline:

1. Read Alignment and Quantification:
   - Salmon (v1.9.0) was used for quasi-mapping and transcript quantification
   - Transcripts were summarized to gene level using tximport (countsFromAbundance='lengthScaledTPM')

2. Quality Control:
   - Samples were assessed for mapping rates and library sizes
   - Initial clustering revealed batch effects between sequencing runs

3. Batch Correction:
   - ComBat-seq was applied to correct for batch effects
   - PVCA analysis confirmed reduction in batch-associated variance from 34.5% to 13.6%

4. Data Normalization:
   - Genes with <0.5 CPM in fewer than 3 samples were filtered
   - TMM normalization was applied using edgeR
   - Final dataset contained 38,045 genes

Differential Expression Analysis:
- edgeR (v3.36.0) was used with a quasi-likelihood framework
- Design matrix included all group combinations (4 lines × 5 timepoints)
- Contrasts tested timepoint comparisons within lines and line comparisons at TP0
- Genes with FDR < 0.05 and |log2FC| > log2(1.2) were considered differentially expressed

WGCNA Analysis:
- Genes significant in ANOVA across all groups (FDR < 0.01, n=24,900) were selected
- Soft threshold power of 9 was used for signed hybrid networks
- Dynamic tree cutting with deepSplit=2 and minModuleSize=30
- Module merging at height 0.2 resulted in 28 modules
- Module-trait correlations were calculated using Pearson correlation

JAG1 Network Analysis:
- Co-expression network constructed using Spearman correlation at TP0
- Genes with |r| > 0.7 were considered part of the JAG1 network
- Differential network analysis compared correlations between leaf types
- Target genes were prioritized based on correlation, expression change, and leaf-type specificity

Statistical Software:
All analyses were performed in R (v4.5.1) using Bioconductor packages.
"

writeLines(methods_text, "results/final_report/methods_text.txt")

## Generate Analysis Summary Report
cat("\nGenerating final analysis summary...\n")

summary_report <- paste0(
  "Soybean Leaf Shape RNA-Seq Analysis Summary\n",
  "==========================================\n\n",
  
  "Dataset Overview:\n",
  "- Total samples: ", nrow(targets), "\n",
  "- Genotypes: ", paste(unique(targets$Line), collapse = ", "), "\n",
  "- Timepoints: ", paste(unique(targets$Timepoint), collapse = ", "), "\n",
  "- Total genes analyzed: ", nrow(master_gene_info), "\n\n",
  
  "Key Findings:\n",
  "1. Differential Expression:\n",
  "   - Genes DE in any contrast: ", sum(master_gene_info$N_DE_Contrasts > 0, na.rm = TRUE), "\n",
  "   - Broad leaf core genes: ", sum(master_gene_info$Is_Broad_Core), "\n",
  "   - Narrow leaf core genes: ", sum(master_gene_info$Is_Narrow_Core), "\n\n",
  
  "2. WGCNA Analysis:\n",
  "   - Modules identified: ", ifelse(exists("net"), max(net$colors), "N/A"), "\n",
  "   - Modules with trait associations: ", 
  ifelse(exists("moduleTraitPvalue"), sum(apply(moduleTraitPvalue, 1, min) < 0.05), "N/A"), "\n\n",
  
  "3. JAG1 Analysis:\n",
  "   - JAG1 network genes at TP0: ", 
  ifelse(exists("jag1_results"), length(jag1_results$network_genes$Gene), "N/A"), "\n",
  "   - Predicted JAG1 targets: ", 
  ifelse(exists("jag1_results"), nrow(jag1_results$predicted_targets), "N/A"), "\n\n",
  
  "4. Key Regulatory Genes:\n",
  "   - Total identified: ", ifelse(exists("key_genes"), nrow(key_genes), "N/A"), "\n",
  "   - Hub genes: ", ifelse(exists("key_genes"), sum(key_genes$Is_Hub), "N/A"), "\n\n",
  
  "Output Files:\n",
  "- Figures: results/final_report/figures/\n",
  "- Tables: results/final_report/tables/\n",
  "- Supplementary: results/final_report/supplementary/\n\n",
  
  "Analysis completed: ", Sys.Date(), "\n"
)

writeLines(summary_report, "results/final_report/analysis_summary.txt")

## Create a README for the results
readme_text <- "
Soybean Leaf Shape RNA-Seq Analysis Results
===========================================

This directory contains the final results from the comprehensive RNA-seq analysis 
of soybean leaf shape determination.

Directory Structure:
-------------------
- figures/: Publication-ready figures
  - Figure1_experimental_overview.pdf: Study design and QC
  - Figure2_differential_expression.pdf: DE analysis summary
  - Figure3_WGCNA_results.pdf: Gene co-expression modules
  - Figure4_JAG1_analysis.pdf: JAG1 regulatory network

- tables/: Main result tables
  - Table1_sample_information.csv: Sample metadata
  - Table2_top_DEGs.csv: Top differentially expressed genes
  - Table3_key_regulatory_genes.csv: Key regulatory genes
  - Table4_JAG1_top_targets.csv: Predicted JAG1 targets

- supplementary/: Complete datasets
  - SuppData1_complete_DE_results.RData: All DE results
  - SuppData2_WGCNA_results.RData: WGCNA analysis objects
  - SuppData3_gene_sets.RData: Defined gene sets
  - SuppData4_normalized_expression.csv.gz: Expression matrix

Key Files:
----------
- analysis_summary.txt: Overview of results
- methods_text.txt: Methods description for publication

For questions about this analysis, contact: [Your contact info]
"

writeLines(readme_text, "results/final_report/README.txt")

## Final session info
cat("\nSaving session information...\n")
sink("results/final_report/session_info.txt")
sessionInfo()
sink()

cat("\n=== Final reporting complete ===\n")
cat("All results exported to: results/final_report/\n")
cat("Ready for publication and sharing!\n")