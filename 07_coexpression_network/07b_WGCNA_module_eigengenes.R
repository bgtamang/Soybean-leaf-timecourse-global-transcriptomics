# Script 19: Module Eigengene Analysis

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 19: MODULE EIGENGENE ANALYSIS\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/19_WGCNA_eigengenes", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "WGCNA",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "reshape2"
)

invisible(lapply(required_packages, library, character.only = TRUE))
enableWGCNAThreads()
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINTS =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/18_WGCNA_prep.RData")
cat("Loaded WGCNA preparation data\n")
cat("  Genes:", ncol(datExpr), "\n")
cat("  Samples:", nrow(datExpr), "\n")
cat("  Modules:", max(net$colors), "\n\n")

# ===== MODULE EIGENGENES =====

cat("========================================\n")
cat("SECTION 2: MODULE EIGENGENES\n")
cat("========================================\n\n")

# Get module eigengenes from network
MEs <- net$MEs

cat("Module eigengenes calculated\n")
cat("  Number of MEs:", ncol(MEs), "\n")
cat("  Samples:", nrow(MEs), "\n\n")

# Order MEs by module number
ME_order <- order(as.numeric(gsub("ME", "", colnames(MEs))))
MEs <- MEs[, ME_order]

# Add sample information
ME_data <- data.frame(
  Sample = rownames(MEs),
  stringsAsFactors = FALSE
)
ME_data <- cbind(ME_data, MEs)

# Parse sample names to get Leaf_type and Timepoint
ME_data$Leaf_type <- targets$Leaf_type
ME_data$Timepoint <- targets$Timepoint
ME_data$Group <- paste(targets$Leaf_type, targets$Timepoint, sep = "_")

cat("Sample metadata added to ME data\n\n")

# ===== MODULE SIZE ANALYSIS =====

cat("========================================\n")
cat("SECTION 3: MODULE SIZE ANALYSIS\n")
cat("========================================\n\n")

# Calculate module sizes
module_sizes <- table(net$colors)
module_info <- data.frame(
  Module = as.numeric(names(module_sizes)),
  Size = as.vector(module_sizes),
  Color = labels2colors(as.numeric(names(module_sizes)))
)
module_info <- module_info[order(module_info$Module), ]

# Calculate percentage of total genes
module_info$Percentage <- round(module_info$Size / sum(module_info$Size) * 100, 2)

cat("Module size summary:\n")
print(head(module_info, 15))

write.csv(module_info, "03_results/tables/WGCNA/module_sizes.csv", row.names = FALSE)
cat("\nSaved: module_sizes.csv\n")

# Plot module sizes
png("03_results/figures/19_WGCNA_eigengenes/module_sizes.png",
    width = 1000, height = 600, res = 120)

module_info_plot <- module_info[module_info$Module != 0, ]  # Exclude grey
module_info_plot <- module_info_plot[order(-module_info_plot$Size), ]

par(mar = c(7, 4, 4, 2))
barplot(
  module_info_plot$Size,
  names.arg = module_info_plot$Color,
  col = module_info_plot$Color,
  las = 2,
  main = "WGCNA Module Sizes",
  ylab = "Number of Genes",
  border = "white"
)
dev.off()
cat("Saved: module_sizes.png\n")

# ===== MODULE EIGENGENE CLUSTERING =====

cat("\n========================================\n")
cat("SECTION 4: MODULE EIGENGENE RELATIONSHIPS\n")
cat("========================================\n\n")

# Calculate ME correlations (excluding grey module ME0)
MEs_no_grey <- MEs[, !grepl("ME0", colnames(MEs))]

# Hierarchical clustering of MEs
ME_diss <- 1 - cor(MEs_no_grey, use = "pairwise.complete.obs")
ME_tree <- hclust(as.dist(ME_diss), method = "average")

# Plot ME dendrogram
png("03_results/figures/19_WGCNA_eigengenes/ME_clustering.png",
    width = 1000, height = 600, res = 120)
par(mar = c(5, 4, 4, 2))
plot(ME_tree, main = "Module Eigengene Clustering",
     xlab = "", sub = "")
abline(h = 0.25, col = "red", lty = 2)
abline(h = 0.20, col = "blue", lty = 2)
legend("topright", legend = c("h = 0.25", "h = 0.20"),
       col = c("red", "blue"), lty = 2)
dev.off()
cat("Saved: ME_clustering.png\n")

# ===== ME HEATMAP ACROSS SAMPLES =====

cat("\n========================================\n")
cat("SECTION 5: ME EXPRESSION PATTERNS\n")
cat("========================================\n\n")

# Create sample annotation
sample_anno <- data.frame(
  Leaf_type = targets$Leaf_type,
  Timepoint = targets$Timepoint,
  row.names = rownames(MEs)
)

# Color annotations
anno_colors <- list(
  Leaf_type = c("Broad" = "#E69F00", "Narrow" = "#56B4E9"),
  Timepoint = c("TP0" = "#1b9e77", "TP1" = "#d95f02", "TP2" = "#7570b3",
                "TP3" = "#e7298a", "TP4" = "#66a61e")
)

# Order samples by Leaf_type and Timepoint
sample_order <- order(sample_anno$Leaf_type, sample_anno$Timepoint)
MEs_ordered <- MEs_no_grey[sample_order, ]
sample_anno_ordered <- sample_anno[sample_order, ]

# Create heatmap
png("03_results/figures/19_WGCNA_eigengenes/ME_heatmap_all_samples.png",
    width = 1200, height = 800, res = 120)

pheatmap(
  t(MEs_ordered),
  annotation_col = sample_anno_ordered,
  annotation_colors = anno_colors,
  cluster_cols = FALSE,  # Keep sample order
  cluster_rows = TRUE,
  show_colnames = FALSE,
  main = "Module Eigengene Expression Across Samples",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize = 10
)

dev.off()
cat("Saved: ME_heatmap_all_samples.png\n")

# ===== ME PATTERNS BY CONDITION =====

cat("\n========================================\n")
cat("SECTION 6: ME PATTERNS BY CONDITION\n")
cat("========================================\n\n")

# Calculate mean ME for each condition
ME_means <- aggregate(MEs_no_grey,
                      by = list(Group = ME_data$Group),
                      FUN = mean)
rownames(ME_means) <- ME_means$Group
ME_means <- ME_means[, -1]

# Parse group info
ME_means$Leaf_type <- sapply(strsplit(rownames(ME_means), "_"), `[`, 1)
ME_means$Timepoint <- sapply(strsplit(rownames(ME_means), "_"), `[`, 2)

# Create ordered levels
timepoint_order <- c("TP0", "TP1", "TP2", "TP3", "TP4")
ME_means$Timepoint <- factor(ME_means$Timepoint, levels = timepoint_order)

# Plot mean ME heatmap
ME_means_matrix <- ME_means[, !colnames(ME_means) %in% c("Leaf_type", "Timepoint")]

# Annotation for condition heatmap
condition_anno <- data.frame(
  Leaf_type = ME_means$Leaf_type,
  Timepoint = as.character(ME_means$Timepoint),
  row.names = rownames(ME_means)
)

png("03_results/figures/19_WGCNA_eigengenes/ME_heatmap_conditions.png",
    width = 1000, height = 800, res = 120)

pheatmap(
  t(as.matrix(ME_means_matrix)),
  annotation_col = condition_anno,
  annotation_colors = anno_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  main = "Mean Module Eigengene by Condition",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize = 10
)

dev.off()
cat("Saved: ME_heatmap_conditions.png\n")

# ===== TEMPORAL PATTERNS =====

cat("\n========================================\n")
cat("SECTION 7: TEMPORAL PATTERNS\n")
cat("========================================\n\n")

# Identify top modules by variance across conditions
ME_var <- apply(ME_means_matrix, 2, var)
top_variable_modules <- names(sort(ME_var, decreasing = TRUE)[1:min(12, ncol(ME_means_matrix))])

cat("Top variable modules:\n")
for (i in 1:length(top_variable_modules)) {
  mod_num <- gsub("ME", "", top_variable_modules[i])
  mod_color <- labels2colors(as.numeric(mod_num))
  cat(sprintf("  %s (%s): variance = %.4f\n",
              top_variable_modules[i], mod_color, ME_var[top_variable_modules[i]]))
}

# Create temporal pattern plots for top modules
png("03_results/figures/19_WGCNA_eigengenes/ME_temporal_patterns.png",
    width = 1400, height = 1000, res = 120)

par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))

for (me_name in top_variable_modules) {
  mod_num <- as.numeric(gsub("ME", "", me_name))
  mod_color <- labels2colors(mod_num)

  # Prepare data
  plot_data <- data.frame(
    ME = ME_data[[me_name]],
    Timepoint = ME_data$Timepoint,
    Leaf_type = ME_data$Leaf_type
  )

  # Calculate means and SE
  means <- aggregate(ME ~ Timepoint + Leaf_type, data = plot_data, FUN = mean)

  # Plot
  broad_means <- means[means$Leaf_type == "Broad", ]
  narrow_means <- means[means$Leaf_type == "Narrow", ]

  # Order by timepoint
  broad_means <- broad_means[order(broad_means$Timepoint), ]
  narrow_means <- narrow_means[order(narrow_means$Timepoint), ]

  ylim_range <- range(c(broad_means$ME, narrow_means$ME))
  ylim_range <- ylim_range + c(-0.1, 0.1) * diff(ylim_range)

  plot(1:5, broad_means$ME, type = "b", pch = 16, col = "#E69F00",
       ylim = ylim_range, xlab = "Timepoint", ylab = "Module Eigengene",
       main = paste0(me_name, " (", mod_color, ")"),
       xaxt = "n", lwd = 2)
  axis(1, at = 1:5, labels = paste0("TP", 0:4))
  lines(1:5, narrow_means$ME, type = "b", pch = 17, col = "#56B4E9", lwd = 2)
  abline(h = 0, lty = 2, col = "gray")

  if (me_name == top_variable_modules[1]) {
    legend("topright", legend = c("Broad", "Narrow"),
           col = c("#E69F00", "#56B4E9"), pch = c(16, 17),
           lwd = 2, cex = 0.8)
  }
}

dev.off()
cat("\nSaved: ME_temporal_patterns.png\n")

# ===== ME NETWORK VISUALIZATION =====

cat("\n========================================\n")
cat("SECTION 8: ME NETWORK\n")
cat("========================================\n\n")

# Plot eigengene network
png("03_results/figures/19_WGCNA_eigengenes/ME_network.png",
    width = 1000, height = 1000, res = 120)

par(cex = 0.8)
plotEigengeneNetworks(
  MEs_no_grey,
  "Module Eigengene Network",
  marDendro = c(0, 4, 1, 2),
  marHeatmap = c(3, 4, 1, 2),
  cex.lab = 0.7,
  plotDendrograms = TRUE,
  xLabelsAngle = 90
)

dev.off()
cat("Saved: ME_network.png\n")

# ===== IDENTIFY MODULES OF INTEREST =====

cat("\n========================================\n")
cat("SECTION 9: MODULES OF INTEREST\n")
cat("========================================\n\n")

# Identify modules with strong leaf type effects
# Compare Broad vs Narrow at each timepoint

leaf_effect <- data.frame(
  Module = gsub("ME", "", colnames(ME_means_matrix)),
  stringsAsFactors = FALSE
)

# Calculate difference at each timepoint
for (tp in paste0("TP", 0:4)) {
  broad_row <- paste0("Broad_", tp)
  narrow_row <- paste0("Narrow_", tp)

  if (broad_row %in% rownames(ME_means_matrix) && narrow_row %in% rownames(ME_means_matrix)) {
    leaf_effect[[paste0("Diff_", tp)]] <- as.numeric(ME_means_matrix[narrow_row, ] - ME_means_matrix[broad_row, ])
  }
}

# Calculate overall leaf effect (mean absolute difference)
diff_cols <- grep("^Diff_", colnames(leaf_effect), value = TRUE)
leaf_effect$Mean_Abs_Diff <- rowMeans(abs(leaf_effect[, diff_cols]))

# Also calculate TP0 difference (when JAG1 is most active)
leaf_effect$TP0_Diff <- leaf_effect$Diff_TP0

# Sort by TP0 difference
leaf_effect <- leaf_effect[order(abs(leaf_effect$TP0_Diff), decreasing = TRUE), ]

cat("Modules with largest Narrow vs Broad difference at TP0:\n")
print(head(leaf_effect[, c("Module", "TP0_Diff", "Mean_Abs_Diff")], 10))

write.csv(leaf_effect, "03_results/tables/WGCNA/module_leaf_effects.csv", row.names = FALSE)
cat("\nSaved: module_leaf_effects.csv\n")

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 10: SAVE CHECKPOINT\n")
cat("========================================\n\n")

save(
  MEs,
  ME_data,
  ME_means,
  module_info,
  leaf_effect,
  ME_tree,
  # Keep previous objects
  datExpr,
  net,
  wgcna_prep,
  recommended_power,
  module_colors,
  targets,
  logCPM,
  file = "03_results/checkpoints/19_WGCNA_eigengenes.RData"
)

cat("Checkpoint saved: 19_WGCNA_eigengenes.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 19 COMPLETE: MODULE EIGENGENE ANALYSIS\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat("  - Module eigengenes calculated for", ncol(MEs), "modules\n")
cat("  - Largest module:", max(module_info$Size[module_info$Module != 0]), "genes\n")
cat("  - Most variable module:", top_variable_modules[1], "\n")
cat("  - Strongest TP0 leaf effect: Module",
    leaf_effect$Module[1], "(diff =", round(leaf_effect$TP0_Diff[1], 3), ")\n\n")

cat("OUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/19_WGCNA_eigengenes.RData\n")
cat("  - Tables: 03_results/tables/WGCNA/\n")
cat("    - module_sizes.csv\n")
cat("    - module_leaf_effects.csv\n")
cat("  - Figures: 03_results/figures/19_WGCNA_eigengenes/\n")
cat("    - module_sizes.png\n")
cat("    - ME_clustering.png\n")
cat("    - ME_heatmap_all_samples.png\n")
cat("    - ME_heatmap_conditions.png\n")
cat("    - ME_temporal_patterns.png\n")
cat("    - ME_network.png\n\n")

cat("NEXT STEP: Run Script 20 for module-trait correlation analysis\n\n")
