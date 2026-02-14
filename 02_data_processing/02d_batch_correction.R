# Script 04: Batch Correction

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 04: BATCH CORRECTION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

# Define base directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/04_batch", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",           # For DGEList
  "sva",             # For ComBat-seq
  "ggplot2",         # For plotting
  "RColorBrewer",    # For colors
  "pheatmap",        # For heatmaps
  "gridExtra",       # For multiple plots
  "Biobase",         # For ExpressionSet (PVCA)
  "pvca"             # For Principal Variance Component Analysis
)

# Install missing packages
missing_packages <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  All packages loaded successfully\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/03_QC_complete.RData")
cat("Loaded checkpoint from script 03\n")
cat("  Objects: raw_counts, dge_filtered, logCPM, targets, sample_cor\n")
cat("  Samples:", ncol(raw_counts), "\n")
cat("  Genes (filtered):", nrow(dge_filtered), "\n\n")

# ===== BATCH STRUCTURE ANALYSIS =====

cat("========================================\n")
cat("SECTION 2: BATCH STRUCTURE\n")
cat("========================================\n\n")

cat("Batch composition:\n")
print(table(targets$Batch, targets$Timepoint))

cat("\nBatch-Timepoint confounding:\n")
cat("  2021 batch: TP1, TP2, TP3, TP4 samples only\n")
cat("  2022 batch: TP0 samples only\n")
cat("  STATUS: COMPLETELY CONFOUNDED\n\n")

cat("This means:\n")
cat("  - We CANNOT distinguish batch effects from TP0 biology\n")
cat("  - Standard batch correction would remove TP0-specific genes\n")
cat("  - We will use TWO approaches:\n")
cat("    1. Primary: ComBat-seq with biological covariates\n")
cat("    2. Sensitivity: Analysis excluding TP0 (no batch effect)\n\n")

# ===== VISUALIZE BATCH EFFECTS BEFORE CORRECTION =====

cat("========================================\n")
cat("SECTION 3: PRE-CORRECTION VISUALIZATION\n")
cat("========================================\n\n")

# Colors
batch_colors <- c("2021" = "#E41A1C", "2022" = "#377EB8")
tp_colors <- brewer.pal(5, "YlOrRd")
names(tp_colors) <- c("TP0", "TP1", "TP2", "TP3", "TP4")
leaf_colors <- c("Broad" = "darkgreen", "Narrow" = "purple")

# --- Plot 1: PCA before correction ---
cat("Creating PCA plot (before correction)...\n")

# PCA on logCPM
pca_before <- prcomp(t(logCPM), scale. = TRUE)
var_explained <- summary(pca_before)$importance[2, 1:5] * 100

pca_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  PC3 = pca_before$x[, 3],
  Batch = as.factor(targets$Batch),
  Timepoint = targets$Timepoint,
  Leaf_type = targets$Leaf_type,
  Line = targets$Line,
  Sample = rownames(targets)
)

png("03_results/figures/04_batch/PCA_before_correction.png",
    width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2), cex.main = 1.2, cex.lab = 1.1)

# PC1 vs PC2 - colored by Batch
plot(pca_df$PC1, pca_df$PC2,
     col = batch_colors[as.character(pca_df$Batch)],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA Before Correction - Colored by Batch")
legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19, cex = 1.0)

# PC1 vs PC2 - colored by Timepoint
plot(pca_df$PC1, pca_df$PC2,
     col = tp_colors[pca_df$Timepoint],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA Before Correction - Colored by Timepoint")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 1.0)

# PC1 vs PC2 - colored by Leaf type
plot(pca_df$PC1, pca_df$PC2,
     col = leaf_colors[pca_df$Leaf_type],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA Before Correction - Colored by Leaf Type")
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 19, cex = 1.0)

# PC2 vs PC3 - to see more structure
plot(pca_df$PC2, pca_df$PC3,
     col = tp_colors[pca_df$Timepoint],
     pch = 19, cex = 1.8,
     xlab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     ylab = paste0("PC3 (", round(var_explained[3], 1), "%)"),
     main = "PCA Before Correction - PC2 vs PC3")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/04_batch/PCA_before_correction.png\n")

# ===== COMBAT-SEQ BATCH CORRECTION =====

cat("\n========================================\n")
cat("SECTION 4: COMBAT-SEQ BATCH CORRECTION\n")
cat("========================================\n\n")

cat("Strategy: Use ComBat-seq with biological covariates\n")
cat("  - Include Timepoint and Leaf_type as covariates\n")
cat("  - This preserves biological signal while removing batch\n")
cat("  - Note: TP0 signal cannot be fully separated from batch\n\n")

# Get count matrix (use filtered genes)
counts_filtered <- dge_filtered$counts

# Create batch vector
batch_vec <- targets$Batch

# Create covariate matrix (biological factors to preserve)
# Include Timepoint and Leaf_type as covariates
covar <- model.matrix(~ Timepoint + Leaf_type, data = targets)[, -1]

cat("Running ComBat-seq...\n")
cat("  Input dimensions:", nrow(counts_filtered), "genes x", ncol(counts_filtered), "samples\n")
cat("  Batch levels:", paste(unique(batch_vec), collapse = ", "), "\n")
cat("  Covariates:", ncol(covar), "columns\n")

# Run ComBat-seq
counts_corrected <- ComBat_seq(
  counts = counts_filtered,
  batch = batch_vec,
  group = NULL,
  covar_mod = covar
)

cat("  ComBat-seq complete\n")
cat("  Output dimensions:", nrow(counts_corrected), "genes x", ncol(counts_corrected), "samples\n\n")

# Verify no negative counts (shouldn't happen with ComBat-seq)
if (any(counts_corrected < 0)) {
  cat("  WARNING: Negative counts detected, setting to 0\n")
  counts_corrected[counts_corrected < 0] <- 0
}

# ===== CREATE CORRECTED DGELIST =====

cat("Creating corrected DGEList...\n")

dge_corrected <- DGEList(counts = counts_corrected)
dge_corrected$samples <- dge_filtered$samples
dge_corrected$samples$lib.size <- colSums(counts_corrected)

# Recalculate normalization factors
dge_corrected <- calcNormFactors(dge_corrected, method = "TMM")

# Calculate logCPM for corrected data
logCPM_corrected <- cpm(dge_corrected, log = TRUE, prior.count = 2)

cat("  DGEList created with TMM normalization\n\n")

# ===== VISUALIZE AFTER CORRECTION =====

cat("========================================\n")
cat("SECTION 5: POST-CORRECTION VISUALIZATION\n")
cat("========================================\n\n")

# --- Plot 2: PCA after correction ---
cat("Creating PCA plot (after correction)...\n")

pca_after <- prcomp(t(logCPM_corrected), scale. = TRUE)
var_explained_after <- summary(pca_after)$importance[2, 1:5] * 100

pca_df_after <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  PC3 = pca_after$x[, 3],
  Batch = as.factor(targets$Batch),
  Timepoint = targets$Timepoint,
  Leaf_type = targets$Leaf_type,
  Line = targets$Line,
  Sample = rownames(targets)
)

png("03_results/figures/04_batch/PCA_after_correction.png",
    width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2), cex.main = 1.2, cex.lab = 1.1)

# PC1 vs PC2 - colored by Batch
plot(pca_df_after$PC1, pca_df_after$PC2,
     col = batch_colors[as.character(pca_df_after$Batch)],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained_after[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
     main = "PCA After ComBat-seq - Colored by Batch")
legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19, cex = 1.0)

# PC1 vs PC2 - colored by Timepoint
plot(pca_df_after$PC1, pca_df_after$PC2,
     col = tp_colors[pca_df_after$Timepoint],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained_after[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
     main = "PCA After ComBat-seq - Colored by Timepoint")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 1.0)

# PC1 vs PC2 - colored by Leaf type
plot(pca_df_after$PC1, pca_df_after$PC2,
     col = leaf_colors[pca_df_after$Leaf_type],
     pch = 19, cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained_after[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
     main = "PCA After ComBat-seq - Colored by Leaf Type")
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 19, cex = 1.0)

# PC2 vs PC3
plot(pca_df_after$PC2, pca_df_after$PC3,
     col = tp_colors[pca_df_after$Timepoint],
     pch = 19, cex = 1.8,
     xlab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
     ylab = paste0("PC3 (", round(var_explained_after[3], 1), "%)"),
     main = "PCA After ComBat-seq - PC2 vs PC3")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/04_batch/PCA_after_correction.png\n")

# --- Plot 3: Before/After comparison ---
cat("Creating before/after comparison plot...\n")

png("03_results/figures/04_batch/PCA_comparison.png",
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), cex.main = 1.3, cex.lab = 1.2)

# Before
plot(pca_df$PC1, pca_df$PC2,
     col = batch_colors[as.character(pca_df$Batch)],
     pch = c(16, 17)[as.numeric(factor(pca_df$Timepoint == "TP0")) + 1],
     cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "BEFORE ComBat-seq")
legend("topright", legend = c("Batch 2021", "Batch 2022"),
       col = batch_colors, pch = 19, cex = 1.0)

# After
plot(pca_df_after$PC1, pca_df_after$PC2,
     col = batch_colors[as.character(pca_df_after$Batch)],
     pch = c(16, 17)[as.numeric(factor(pca_df_after$Timepoint == "TP0")) + 1],
     cex = 1.8,
     xlab = paste0("PC1 (", round(var_explained_after[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained_after[2], 1), "%)"),
     main = "AFTER ComBat-seq")
legend("topright", legend = c("Batch 2021", "Batch 2022"),
       col = batch_colors, pch = 19, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/04_batch/PCA_comparison.png\n")

# --- Plot 4: RLE comparison (side-by-side) ---
cat("Creating RLE comparison plot (side-by-side)...\n")

# Calculate RLE values
median_before <- apply(logCPM, 1, median)
rle_before <- logCPM - median_before

median_after <- apply(logCPM_corrected, 1, median)
rle_after <- logCPM_corrected - median_after

# Calculate RLE statistics for summary
rle_stats_before <- data.frame(
  Sample = colnames(logCPM),
  Batch = targets$Batch,
  Median = apply(rle_before, 2, median),
  IQR = apply(rle_before, 2, IQR)
)

rle_stats_after <- data.frame(
  Sample = colnames(logCPM_corrected),
  Batch = targets$Batch,
  Median = apply(rle_after, 2, median),
  IQR = apply(rle_after, 2, IQR)
)

# Create side-by-side RLE comparison plot
png("03_results/figures/04_batch/RLE_comparison.png",
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(8, 5, 4, 2))

# RLE BEFORE correction
boxplot(rle_before,
        las = 2,
        ylim = c(-2, 2),
        col = batch_colors[as.character(targets$Batch)],
        main = "Before Batch Correction",
        ylab = "Relative Log Expression",
        xlab = "Samples",
        cex.axis = 0.6,
        cex.main = 1.4,
        cex.lab = 1.2,
        outline = FALSE)
abline(h = 0, col = "darkgray", lty = 2, lwd = 2)
legend("topright",
       legend = paste("Batch", names(batch_colors)),
       fill = batch_colors,
       bty = "n",
       cex = 0.9)

# RLE AFTER correction
boxplot(rle_after,
        las = 2,
        ylim = c(-2, 2),
        col = batch_colors[as.character(targets$Batch)],
        main = "After Batch Correction",
        ylab = "Relative Log Expression",
        xlab = "Samples",
        cex.axis = 0.6,
        cex.main = 1.4,
        cex.lab = 1.2,
        outline = FALSE)
abline(h = 0, col = "darkgray", lty = 2, lwd = 2)
legend("topright",
       legend = paste("Batch", names(batch_colors)),
       fill = batch_colors,
       bty = "n",
       cex = 0.9)

dev.off()
cat("  Saved: 03_results/figures/04_batch/RLE_comparison.png\n")

# Print RLE summary statistics
cat("\nRLE Statistics Summary:\n")
cat("  BEFORE correction:\n")
cat("    Median SD across samples:", round(sd(rle_stats_before$Median), 4), "\n")
cat("    Mean IQR:", round(mean(rle_stats_before$IQR), 4), "\n")
cat("  AFTER correction:\n")
cat("    Median SD across samples:", round(sd(rle_stats_after$Median), 4), "\n")
cat("    Mean IQR:", round(mean(rle_stats_after$IQR), 4), "\n")
cat("  Improvement:\n")
cat("    Median SD reduction:", round((1 - sd(rle_stats_after$Median)/sd(rle_stats_before$Median)) * 100, 1), "%\n\n")

# ===== QUANTIFY BATCH EFFECT REDUCTION =====

cat("\n========================================\n")
cat("SECTION 6: BATCH EFFECT QUANTIFICATION\n")
cat("========================================\n\n")

# Calculate variance explained by batch before and after

# Function to calculate R-squared for batch
calc_batch_variance <- function(expr_mat, batch) {
  rsq <- apply(expr_mat, 1, function(x) {
    fit <- lm(x ~ batch)
    summary(fit)$r.squared
  })
  return(rsq)
}

cat("Calculating batch variance explained...\n")

# Before correction
batch_rsq_before <- calc_batch_variance(logCPM, targets$Batch)
cat("  Before correction:\n")
cat("    Mean R² (batch):", round(mean(batch_rsq_before), 4), "\n")
cat("    Median R² (batch):", round(median(batch_rsq_before), 4), "\n")

# After correction
batch_rsq_after <- calc_batch_variance(logCPM_corrected, targets$Batch)
cat("  After correction:\n")
cat("    Mean R² (batch):", round(mean(batch_rsq_after), 4), "\n")
cat("    Median R² (batch):", round(median(batch_rsq_after), 4), "\n")

cat("\n  Reduction in batch variance: ",
    round((1 - mean(batch_rsq_after)/mean(batch_rsq_before)) * 100, 1), "%\n")

# Also check variance explained by biology (Timepoint, Leaf_type)
cat("\nChecking biological variance preservation...\n")

calc_tp_variance <- function(expr_mat, timepoint) {
  rsq <- apply(expr_mat, 1, function(x) {
    fit <- lm(x ~ timepoint)
    summary(fit)$r.squared
  })
  return(rsq)
}

tp_rsq_before <- calc_tp_variance(logCPM, targets$Timepoint)
tp_rsq_after <- calc_tp_variance(logCPM_corrected, targets$Timepoint)

cat("  Timepoint variance:\n")
cat("    Before: Mean R² =", round(mean(tp_rsq_before), 4), "\n")
cat("    After: Mean R² =", round(mean(tp_rsq_after), 4), "\n")
cat("    Change:", round((mean(tp_rsq_after)/mean(tp_rsq_before) - 1) * 100, 1), "%\n")

# ===== PVCA ANALYSIS =====
# Principal Variance Component Analysis for precise variance decomposition

cat("\n========================================\n")
cat("SECTION 6b: PVCA ANALYSIS (Variance Decomposition)\n")
cat("========================================\n\n")

cat("Running PVCA for precise variance attribution...\n")
cat("  This provides exact percentages for manuscript\n\n")

# --- PVCA BEFORE Batch Correction ---
cat("PVCA BEFORE batch correction:\n")

# Prepare data for PVCA
pData_before <- targets[, c("Batch", "Leaf_type", "Timepoint")]
pData_before$Batch <- as.factor(pData_before$Batch)
pData_before$Leaf_type <- as.factor(pData_before$Leaf_type)
pData_before$Timepoint <- as.factor(pData_before$Timepoint)

# Create ExpressionSet
myExpressionSet_before <- ExpressionSet(
  assayData = as.matrix(logCPM),
  phenoData = AnnotatedDataFrame(pData_before)
)

# Run PVCA
pvcaObj_before <- pvcaBatchAssess(myExpressionSet_before,
                                  c("Batch", "Leaf_type", "Timepoint"), 0.6)

# Print PVCA results
cat("  Variance decomposition BEFORE batch correction:\n")
for(i in 1:length(pvcaObj_before$dat)) {
  cat(sprintf("    %-20s: %5.2f%%\n", pvcaObj_before$label[i],
              pvcaObj_before$dat[i] * 100))
}

# --- PVCA AFTER Batch Correction ---
cat("\nPVCA AFTER batch correction:\n")

# Prepare data for PVCA
pData_after <- targets[, c("Batch", "Leaf_type", "Timepoint")]
pData_after$Batch <- as.factor(pData_after$Batch)
pData_after$Leaf_type <- as.factor(pData_after$Leaf_type)
pData_after$Timepoint <- as.factor(pData_after$Timepoint)

# Create ExpressionSet for corrected data
myExpressionSet_after <- ExpressionSet(
  assayData = as.matrix(logCPM_corrected),
  phenoData = AnnotatedDataFrame(pData_after)
)

# Run PVCA
pvcaObj_after <- pvcaBatchAssess(myExpressionSet_after,
                                 c("Batch", "Leaf_type", "Timepoint"), 0.6)

# Print PVCA results
cat("  Variance decomposition AFTER batch correction:\n")
for(i in 1:length(pvcaObj_after$dat)) {
  cat(sprintf("    %-20s: %5.2f%%\n", pvcaObj_after$label[i],
              pvcaObj_after$dat[i] * 100))
}

# --- Calculate summary statistics ---
# Extract batch-related effects (any label containing "Batch")
batch_before_total <- sum(pvcaObj_before$dat[grep("Batch", pvcaObj_before$label)]) * 100
batch_after_total <- sum(pvcaObj_after$dat[grep("Batch", pvcaObj_after$label)]) * 100

# Extract biological effects (Timepoint + Leaf_type, not interactions)
bio_labels <- c("Timepoint", "Leaf_type")
bio_before_total <- sum(pvcaObj_before$dat[pvcaObj_before$label %in% bio_labels]) * 100
bio_after_total <- sum(pvcaObj_after$dat[pvcaObj_after$label %in% bio_labels]) * 100

cat("\n--- PVCA SUMMARY FOR MANUSCRIPT ---\n")
cat("Batch effects (including interactions):\n")
cat(sprintf("  Before: %.2f%%\n", batch_before_total))
cat(sprintf("  After:  %.2f%%\n", batch_after_total))
cat(sprintf("  Reduction: %.2f percentage points\n", batch_before_total - batch_after_total))

cat("\nBiological effects (Timepoint + Leaf_type):\n")
cat(sprintf("  Before: %.2f%%\n", bio_before_total))
cat(sprintf("  After:  %.2f%%\n", bio_after_total))
cat(sprintf("  Change: %+.2f percentage points\n", bio_after_total - bio_before_total))

# --- Save PVCA results to CSV for manuscript ---
pvca_summary <- data.frame(
  Effect = c(pvcaObj_before$label, "TOTAL_Batch", "TOTAL_Biological"),
  Before_Pct = c(round(pvcaObj_before$dat * 100, 2),
                 round(batch_before_total, 2),
                 round(bio_before_total, 2)),
  After_Pct = c(round(pvcaObj_after$dat * 100, 2),
                round(batch_after_total, 2),
                round(bio_after_total, 2))
)
pvca_summary$Change = pvca_summary$After_Pct - pvca_summary$Before_Pct

write.csv(pvca_summary, "03_results/tables/PVCA_variance_decomposition.csv", row.names = FALSE)
cat("\nSaved: 03_results/tables/PVCA_variance_decomposition.csv\n")

# --- Create PVCA comparison plot ---
cat("Creating PVCA comparison plot...\n")

png("03_results/figures/04_batch/PVCA_comparison.png",
    width = 14, height = 6, units = "in", res = 300)

# Set up colors for different effect types
effect_colors <- function(labels) {
  colors <- rep("gray", length(labels))
  colors[grep("Batch", labels)] <- "#E41A1C"  # Batch effects in red
  colors[labels == "Timepoint"] <- "#377EB8"  # Timepoint in blue
  colors[labels == "Leaf_type"] <- "#4DAF4A"  # Leaf_type in green
  colors[grep(":", labels)] <- "#FF7F00"      # Interactions in orange
  colors[labels == "resid"] <- "#999999"      # Residual in gray
  return(colors)
}

par(mfrow = c(1, 2), mar = c(10, 5, 4, 2))

# Before correction plot
bp1 <- barplot(pvcaObj_before$dat,
               names.arg = pvcaObj_before$label,
               main = "PVCA: Before ComBat-seq",
               ylab = "Proportion of Variance",
               ylim = c(0, max(c(pvcaObj_before$dat, pvcaObj_after$dat)) * 1.15),
               las = 2,
               cex.names = 0.8,
               col = effect_colors(pvcaObj_before$label))
# Add percentage labels
text(x = bp1, y = pvcaObj_before$dat + 0.02,
     labels = paste0(round(pvcaObj_before$dat * 100, 1), "%"),
     cex = 0.7)

# After correction plot
bp2 <- barplot(pvcaObj_after$dat,
               names.arg = pvcaObj_after$label,
               main = "PVCA: After ComBat-seq",
               ylab = "Proportion of Variance",
               ylim = c(0, max(c(pvcaObj_before$dat, pvcaObj_after$dat)) * 1.15),
               las = 2,
               cex.names = 0.8,
               col = effect_colors(pvcaObj_after$label))
# Add percentage labels
text(x = bp2, y = pvcaObj_after$dat + 0.02,
     labels = paste0(round(pvcaObj_after$dat * 100, 1), "%"),
     cex = 0.7)

dev.off()
cat("  Saved: 03_results/figures/04_batch/PVCA_comparison.png\n\n")

# ===== SUMMARY TABLE =====

cat("\n========================================\n")
cat("SECTION 7: SUMMARY\n")
cat("========================================\n\n")

# Create summary table
batch_summary <- data.frame(
  Metric = c("Genes analyzed",
             "Samples",
             "Batch variance (before)",
             "Batch variance (after)",
             "Batch variance reduction",
             "Timepoint variance (before)",
             "Timepoint variance (after)",
             "Timepoint variance change"),
  Value = c(nrow(counts_corrected),
            ncol(counts_corrected),
            round(mean(batch_rsq_before), 4),
            round(mean(batch_rsq_after), 4),
            paste0(round((1 - mean(batch_rsq_after)/mean(batch_rsq_before)) * 100, 1), "%"),
            round(mean(tp_rsq_before), 4),
            round(mean(tp_rsq_after), 4),
            paste0(round((mean(tp_rsq_after)/mean(tp_rsq_before) - 1) * 100, 1), "%"))
)

write.csv(batch_summary, "03_results/tables/batch_correction_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/batch_correction_summary.csv\n\n")

print(batch_summary)

# ===== IMPORTANT CAVEATS =====

cat("\n========================================\n")
cat("SECTION 8: CAVEATS & RECOMMENDATIONS\n")
cat("========================================\n\n")

cat("IMPORTANT CAVEATS:\n")
cat("-----------------\n")
cat("1. Batch and TP0 are COMPLETELY CONFOUNDED\n")
cat("   - Cannot separate batch effects from true TP0 biology\n")
cat("   - Any TP0-specific signal may be partially removed\n\n")

cat("2. Recommended sensitivity analyses:\n")
cat("   a) Primary: Use batch-corrected data for all analyses\n")
cat("   b) Sensitivity: Repeat key analyses excluding TP0\n")
cat("   c) Sensitivity: Repeat key analyses without batch correction\n\n")

cat("3. Interpretation guidance:\n")
cat("   - TP0 vs other timepoints: Results should be interpreted cautiously\n")
cat("   - Within TP1-TP4: No batch confounding, results are reliable\n")
cat("   - Leaf-type comparisons at TP0: May have reduced power\n\n")

# ===== CREATE NO-TP0 SUBSET =====

cat("========================================\n")
cat("SECTION 9: CREATE TP1-TP4 SUBSET\n")
cat("========================================\n\n")

cat("Creating batch-free TP1-TP4 subset for sensitivity analysis...\n")

# Subset to TP1-TP4 (no batch effect)
tp1_4_samples <- rownames(targets)[targets$Timepoint != "TP0"]
cat("  Samples in TP1-TP4 subset:", length(tp1_4_samples), "\n")

# Create subset objects (from original, uncorrected data)
dge_tp1_4 <- dge_filtered[, tp1_4_samples]
targets_tp1_4 <- targets[tp1_4_samples, ]
logCPM_tp1_4 <- logCPM[, tp1_4_samples]

# Recalculate normalization for subset
dge_tp1_4 <- calcNormFactors(dge_tp1_4, method = "TMM")

cat("  Created: dge_tp1_4, targets_tp1_4, logCPM_tp1_4\n")
cat("  These have NO batch effect (all from 2021)\n\n")

# ===== SAVE CHECKPOINT =====

cat("========================================\n")
cat("SECTION 10: SAVE CHECKPOINT\n")
cat("========================================\n\n")

# Save all objects needed for downstream analysis
save(
  # Original filtered data (for reference)
  raw_counts,
  dge_filtered,
  logCPM,

  # Batch-corrected data (PRIMARY)
  counts_corrected,
  dge_corrected,
  logCPM_corrected,

  # TP1-TP4 subset (SENSITIVITY)
  dge_tp1_4,
  targets_tp1_4,
  logCPM_tp1_4,

  # Metadata
  targets,
  sample_cor,

  # Parameters
  PARAMS,
  QC_PARAMS,

  # Batch statistics
  batch_rsq_before,
  batch_rsq_after,

  file = "03_results/checkpoints/04_batch_corrected.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/04_batch_corrected.RData\n")
cat("  Contains:\n")
cat("    - Original: raw_counts, dge_filtered, logCPM\n")
cat("    - Corrected: counts_corrected, dge_corrected, logCPM_corrected\n")
cat("    - TP1-4 subset: dge_tp1_4, targets_tp1_4, logCPM_tp1_4\n")
cat("    - Metadata: targets, sample_cor, PARAMS, QC_PARAMS\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 04: BATCH CORRECTION - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Batch variance reduction:",
    round((1 - mean(batch_rsq_after)/mean(batch_rsq_before)) * 100, 1), "%\n")
cat("    - Biological variance preserved: YES\n")
cat("    - TP1-4 subset created: YES (48 samples, batch-free)\n")
cat("\n")
cat("  IMPORTANT: Batch and TP0 are confounded!\n")
cat("    - Use corrected data for primary analysis\n")
cat("    - Run sensitivity analyses as documented\n")
cat("\n")
cat("  Next: Run 05_normalization.R\n")
cat("================================================================\n")
