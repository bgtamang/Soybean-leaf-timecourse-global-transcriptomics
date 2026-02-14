# Script 07: PCA and Clustering

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 07: PCA AND CLUSTERING\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/07_PCA", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "ggrepel",
  "RColorBrewer",
  "pheatmap",
  "dendextend",
  "factoextra",
  "dplyr",
  "tidyr"
)

# Install missing
missing <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (pkg in missing) {
    tryCatch({
      if (pkg %in% c("edgeR", "limma")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }, error = function(e) cat("Could not install:", pkg, "\n"))
  }
}

invisible(lapply(required_packages, function(p) {
  tryCatch(library(p, character.only = TRUE), error = function(e) NULL)
}))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded checkpoint from script 06\n")
cat("  Using PRIMARY dataset (all 60 samples)\n\n")

# Use the primary dataset
expr_mat <- v_primary$E
targets <- targets_primary

cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n\n")

# ===== DEFINE COLORS =====

# Consistent color scheme
batch_colors <- c("2021" = "#E41A1C", "2022" = "#377EB8")
leaf_colors <- c("Broad" = "#2E8B57", "Narrow" = "#9932CC")
line_colors <- c("PI532462A" = "#1B9E77", "LD112170" = "#D95F02",
                 "PI612713B" = "#7570B3", "PI547745" = "#E7298A")
tp_colors <- c("TP0" = "#FFFFB2", "TP1" = "#FECC5C", "TP2" = "#FD8D3C",
               "TP3" = "#F03B20", "TP4" = "#BD0026")

# ===== SECTION 2: PCA ANALYSIS =====

cat("========================================\n")
cat("SECTION 2: PCA ANALYSIS\n")
cat("========================================\n\n")

# Perform PCA
cat("Running PCA...\n")
pca_result <- prcomp(t(expr_mat), scale. = TRUE, center = TRUE)

# Variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100
cum_var <- summary(pca_result)$importance[3, ] * 100

cat("Variance explained by top PCs:\n")
for (i in 1:5) {
  cat("  PC", i, ":", round(var_explained[i], 2), "% (cumulative:",
      round(cum_var[i], 2), "%)\n")
}

# Create PCA data frame
pca_df <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  PC4 = pca_result$x[, 4],
  PC5 = pca_result$x[, 5],
  Line = targets$Line,
  Leaf_type = targets$Leaf_type,
  Timepoint = targets$Timepoint,
  Batch = as.factor(targets$Batch),
  Group = targets$Group
)

# ===== SECTION 3: PCA VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 3: PCA VISUALIZATIONS\n")
cat("========================================\n\n")

# --- Plot 1: PC1 vs PC2 multi-panel ---
cat("Creating main PCA plots...\n")

png("03_results/figures/07_PCA/PCA_PC1_PC2_multipanel.png",
    width = 14, height = 12, units = "in", res = 300)

par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), cex.main = 1.3, cex.lab = 1.2)

# By Leaf type
plot(pca_df$PC1, pca_df$PC2,
     col = leaf_colors[pca_df$Leaf_type],
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Leaf Type")
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 19, cex = 1.0)

# By Timepoint
plot(pca_df$PC1, pca_df$PC2,
     col = tp_colors[pca_df$Timepoint],
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Timepoint")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.9)

# By Line
plot(pca_df$PC1, pca_df$PC2,
     col = line_colors[pca_df$Line],
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Line")
legend("topright", legend = names(line_colors), col = line_colors, pch = 19, cex = 0.9)

# By Batch
plot(pca_df$PC1, pca_df$PC2,
     col = batch_colors[as.character(pca_df$Batch)],
     pch = 19, cex = 2,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA - Colored by Batch")
legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/07_PCA/PCA_PC1_PC2_multipanel.png\n")

# --- Plot 2: PC1-PC2 with shapes for timepoint, color for leaf type ---
cat("Creating combined PCA plot...\n")

p_combined <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Leaf_type, shape = Timepoint)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = leaf_colors) +
  scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
  theme_bw(base_size = 14) +
  labs(title = "PCA - Leaf Type and Developmental Stage",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right")

ggsave("03_results/figures/07_PCA/PCA_combined.png", p_combined,
       width = 10, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/07_PCA/PCA_combined.png\n")

# --- Plot 3: Faceted by Line ---
cat("Creating faceted PCA plot...\n")

p_facet <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Timepoint)) +
  geom_point(size = 3) +
  scale_color_manual(values = tp_colors) +
  facet_wrap(~Line, ncol = 2) +
  theme_bw(base_size = 14) +
  labs(title = "PCA by Line - Developmental Trajectory",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/07_PCA/PCA_by_line.png", p_facet,
       width = 12, height = 10, dpi = 300)
cat("  Saved: 03_results/figures/07_PCA/PCA_by_line.png\n")

# --- Plot 4: PC2 vs PC3 (often shows developmental trajectory) ---
cat("Creating PC2 vs PC3 plot...\n")

png("03_results/figures/07_PCA/PCA_PC2_PC3.png",
    width = 10, height = 8, units = "in", res = 300)

par(cex.main = 1.3, cex.lab = 1.2)
plot(pca_df$PC2, pca_df$PC3,
     col = tp_colors[pca_df$Timepoint],
     pch = c(15, 16, 17, 18)[as.numeric(factor(pca_df$Line))],
     cex = 2,
     xlab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     ylab = paste0("PC3 (", round(var_explained[3], 1), "%)"),
     main = "PCA - PC2 vs PC3 (Developmental Trajectory)")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.9,
       title = "Timepoint")
legend("bottomright", legend = names(line_colors), pch = c(15, 16, 17, 18), cex = 0.9,
       title = "Line")

dev.off()
cat("  Saved: 03_results/figures/07_PCA/PCA_PC2_PC3.png\n")

# --- Plot 5: Scree plot ---
cat("Creating scree plot...\n")

png("03_results/figures/07_PCA/PCA_scree.png",
    width = 10, height = 6, units = "in", res = 300)

par(mfrow = c(1, 2), cex.main = 1.2, cex.lab = 1.1)

# Variance explained
barplot(var_explained[1:10],
        names.arg = paste0("PC", 1:10),
        col = "steelblue",
        main = "Variance Explained by PC",
        ylab = "% Variance",
        las = 2)

# Cumulative variance
plot(1:10, cum_var[1:10], type = "b",
     pch = 19, col = "darkred",
     xlab = "Principal Component",
     ylab = "Cumulative % Variance",
     main = "Cumulative Variance Explained")
abline(h = 80, lty = 2, col = "gray")
text(8, 82, "80% threshold", cex = 0.9)

dev.off()
cat("  Saved: 03_results/figures/07_PCA/PCA_scree.png\n")

# ===== SECTION 4: HIERARCHICAL CLUSTERING =====

cat("\n========================================\n")
cat("SECTION 4: HIERARCHICAL CLUSTERING\n")
cat("========================================\n\n")

# Calculate sample distances
cat("Calculating sample distances...\n")
sample_cor <- cor(expr_mat, method = "pearson")
sample_dist <- as.dist(1 - sample_cor)

# Hierarchical clustering
hc <- hclust(sample_dist, method = "average")

# --- Plot 6: Dendrogram ---
cat("Creating dendrogram...\n")

png("03_results/figures/07_PCA/dendrogram_full.png",
    width = 16, height = 8, units = "in", res = 300)

# Create dendrogram
dend <- as.dendrogram(hc)

# Color labels by leaf type
label_colors <- leaf_colors[targets$Leaf_type[order.dendrogram(dend)]]
labels_colors(dend) <- label_colors

par(mar = c(12, 5, 4, 2))
plot(dend, main = "Sample Clustering (Average Linkage)",
     ylab = "Distance (1 - correlation)",
     cex.main = 1.3)
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 15, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/07_PCA/dendrogram_full.png\n")

# --- Plot 7: Heatmap with clustering ---
cat("Creating correlation heatmap...\n")

# Annotation
annotation_col <- data.frame(
  Leaf_type = targets$Leaf_type,
  Timepoint = targets$Timepoint,
  Line = targets$Line,
  Batch = as.factor(targets$Batch),
  row.names = colnames(expr_mat)
)

annotation_colors <- list(
  Leaf_type = leaf_colors,
  Timepoint = tp_colors,
  Line = line_colors,
  Batch = batch_colors
)

png("03_results/figures/07_PCA/correlation_heatmap.png",
    width = 14, height = 12, units = "in", res = 300)

pheatmap(sample_cor,
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         breaks = seq(0.85, 1, length.out = 101),
         main = "Sample Correlation Heatmap",
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         clustering_method = "average")

dev.off()
cat("  Saved: 03_results/figures/07_PCA/correlation_heatmap.png\n")

# ===== SECTION 5: PCA LOADINGS =====

cat("\n========================================\n")
cat("SECTION 5: PCA LOADINGS\n")
cat("========================================\n\n")

# Get top genes contributing to each PC
cat("Identifying top genes for each PC...\n")

loadings <- pca_result$rotation
n_top <- 50

# Top genes for PC1
pc1_top <- names(sort(abs(loadings[, 1]), decreasing = TRUE))[1:n_top]
pc2_top <- names(sort(abs(loadings[, 2]), decreasing = TRUE))[1:n_top]
pc3_top <- names(sort(abs(loadings[, 3]), decreasing = TRUE))[1:n_top]

# Save loadings
loadings_df <- data.frame(
  Gene = rownames(loadings),
  PC1 = loadings[, 1],
  PC2 = loadings[, 2],
  PC3 = loadings[, 3],
  PC4 = loadings[, 4],
  PC5 = loadings[, 5]
)
loadings_df <- loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ]

write.csv(loadings_df, "03_results/tables/PCA_loadings.csv", row.names = FALSE)
cat("  Saved: 03_results/tables/PCA_loadings.csv\n")

# Check if JAG1/JAG2 are top contributors
cat("\nJAG genes in PCA loadings:\n")
if (PARAMS$JAG1 %in% rownames(loadings)) {
  jag1_rank_pc1 <- which(names(sort(abs(loadings[, 1]), decreasing = TRUE)) == PARAMS$JAG1)
  jag1_rank_pc2 <- which(names(sort(abs(loadings[, 2]), decreasing = TRUE)) == PARAMS$JAG1)
  cat("  JAG1 rank - PC1:", jag1_rank_pc1, ", PC2:", jag1_rank_pc2, "\n")
}

# ===== SECTION 6: SUMMARY =====

cat("\n========================================\n")
cat("SECTION 6: SUMMARY\n")
cat("========================================\n\n")

# Create summary table
pca_summary <- data.frame(
  Component = paste0("PC", 1:10),
  Variance_Pct = round(var_explained[1:10], 2),
  Cumulative_Pct = round(cum_var[1:10], 2)
)

write.csv(pca_summary, "03_results/tables/PCA_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/PCA_summary.csv\n\n")

cat("PCA Summary:\n")
print(pca_summary[1:5, ])

cat("\nKey observations:\n")
cat("  - PC1 explains", round(var_explained[1], 1), "% of variance\n")
cat("  - First 5 PCs explain", round(cum_var[5], 1), "% of variance\n")
cat("  - Samples cluster primarily by: (check plots)\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 07: PCA AND CLUSTERING - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Outputs:\n")
cat("    - PCA plots: 07_PCA/*.png\n")
cat("    - Loadings: PCA_loadings.csv\n")
cat("    - Summary: PCA_summary.csv\n")
cat("\n")
cat("  Next: Run 08_variance_analysis.R\n")
cat("================================================================\n")
