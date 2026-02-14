# Script 05: Normalization

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 05: NORMALIZATION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("03_results/figures/05_normalization", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",           # For DGEList and normalization
  "limma",           # For voom transformation
  "ggplot2",         # For plotting
  "RColorBrewer",    # For colors
  "pheatmap",        # For heatmaps
  "dplyr"            # For data manipulation
)

missing <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing)
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/04_batch_corrected.RData")
cat("Loaded checkpoint from script 04\n")
cat("  Objects loaded:\n")
cat("    - dge_corrected: Batch-corrected counts\n")
cat("    - logCPM_corrected: Log-transformed CPM (corrected)\n")
cat("    - dge_tp1_4: TP1-TP4 subset (batch-free)\n")
cat("    - targets: Sample metadata\n\n")

cat("Data dimensions:\n")
cat("  All samples: ", nrow(dge_corrected), " genes x ", ncol(dge_corrected), " samples\n", sep = "")
cat("  TP1-4 subset: ", nrow(dge_tp1_4), " genes x ", ncol(dge_tp1_4), " samples\n\n", sep = "")

# ===== SECTION 2: VERIFY NORMALIZATION =====

cat("========================================\n")
cat("SECTION 2: VERIFY NORMALIZATION\n")
cat("========================================\n\n")

# Check TMM normalization factors
cat("TMM Normalization Factors (batch-corrected data):\n")
norm_factors <- dge_corrected$samples$norm.factors
cat("  Mean:", round(mean(norm_factors), 4), "\n")
cat("  SD:", round(sd(norm_factors), 4), "\n")
cat("  Range:", round(min(norm_factors), 4), "to", round(max(norm_factors), 4), "\n\n")

# Library sizes after correction
cat("Library Sizes (batch-corrected):\n")
lib_sizes <- dge_corrected$samples$lib.size
cat("  Mean:", format(round(mean(lib_sizes)), big.mark = ","), "\n")
cat("  Median:", format(round(median(lib_sizes)), big.mark = ","), "\n")
cat("  Range:", format(min(lib_sizes), big.mark = ","), "to",
    format(max(lib_sizes), big.mark = ","), "\n\n")

# ===== SECTION 3: CREATE DESIGN MATRICES =====

cat("========================================\n")
cat("SECTION 3: CREATE DESIGN MATRICES\n")
cat("========================================\n\n")

# Design matrix for full analysis (all 60 samples)
# Using Group (Line_Timepoint) for maximum flexibility
cat("Creating design matrix for full dataset...\n")

# Ensure Group is a factor with proper levels
targets$Group <- factor(targets$Group)
cat("  Groups:", length(levels(targets$Group)), "levels\n")

# Design: ~0 + Group (means model for flexibility in contrasts)
design_full <- model.matrix(~0 + Group, data = targets)
colnames(design_full) <- gsub("Group", "", colnames(design_full))

cat("  Design matrix dimensions:", nrow(design_full), "x", ncol(design_full), "\n")
cat("  Columns:", paste(head(colnames(design_full), 5), collapse = ", "), "...\n\n")

# Design matrix for TP1-4 subset (batch-free)
cat("Creating design matrix for TP1-4 subset...\n")
targets_tp1_4$Group <- factor(targets_tp1_4$Group)
design_tp1_4 <- model.matrix(~0 + Group, data = targets_tp1_4)
colnames(design_tp1_4) <- gsub("Group", "", colnames(design_tp1_4))

cat("  Design matrix dimensions:", nrow(design_tp1_4), "x", ncol(design_tp1_4), "\n\n")

# ===== SECTION 4: VOOM TRANSFORMATION =====

cat("========================================\n")
cat("SECTION 4: VOOM TRANSFORMATION\n")
cat("========================================\n\n")

cat("Performing voom transformation (full dataset)...\n")

# Voom with quality weights for full dataset
png("03_results/figures/05_normalization/voom_mean_variance_full.png",
    width = 10, height = 8, units = "in", res = 300)
par(cex.main = 1.3, cex.lab = 1.2)
v_full <- voom(dge_corrected, design_full, plot = TRUE)
title("Voom Mean-Variance Trend (All 60 Samples)")
dev.off()
cat("  Saved: 03_results/figures/05_normalization/voom_mean_variance_full.png\n")

cat("  Voom object created for full dataset\n")
cat("    - E: Expression matrix (log2-CPM)\n")
cat("    - weights: Precision weights\n")
cat("    - design: Design matrix\n\n")

# Voom for TP1-4 subset
cat("Performing voom transformation (TP1-4 subset)...\n")

png("03_results/figures/05_normalization/voom_mean_variance_tp1_4.png",
    width = 10, height = 8, units = "in", res = 300)
par(cex.main = 1.3, cex.lab = 1.2)
v_tp1_4 <- voom(dge_tp1_4, design_tp1_4, plot = TRUE)
title("Voom Mean-Variance Trend (TP1-4 Only, Batch-Free)")
dev.off()
cat("  Saved: 03_results/figures/05_normalization/voom_mean_variance_tp1_4.png\n\n")

# ===== SECTION 5: EXPRESSION DISTRIBUTIONS =====

cat("========================================\n")
cat("SECTION 5: EXPRESSION DISTRIBUTIONS\n")
cat("========================================\n\n")

# Colors
batch_colors <- c("2021" = "#E41A1C", "2022" = "#377EB8")
leaf_colors <- c("Broad" = "darkgreen", "Narrow" = "purple")
tp_colors <- brewer.pal(5, "YlOrRd")
names(tp_colors) <- c("TP0", "TP1", "TP2", "TP3", "TP4")

# --- Plot 1: Density distributions ---
cat("Creating density distribution plots...\n")

png("03_results/figures/05_normalization/expression_density.png",
    width = 12, height = 8, units = "in", res = 300)

par(mfrow = c(1, 2), cex.main = 1.2, cex.lab = 1.1)

# Before voom (logCPM)
plotDensities(logCPM_corrected,
              col = tp_colors[targets$Timepoint],
              main = "Expression Density (logCPM)",
              legend = FALSE)
legend("topright", legend = names(tp_colors), col = tp_colors, lty = 1, cex = 0.8)

# After voom
plotDensities(v_full$E,
              col = tp_colors[targets$Timepoint],
              main = "Expression Density (voom)",
              legend = FALSE)
legend("topright", legend = names(tp_colors), col = tp_colors, lty = 1, cex = 0.8)

dev.off()
cat("  Saved: 03_results/figures/05_normalization/expression_density.png\n")

# --- Plot 2: Boxplots by group ---
cat("Creating expression boxplots...\n")

# Prepare data for ggplot
expr_long <- data.frame(
  Sample = rep(colnames(v_full$E), each = nrow(v_full$E)),
  Expression = as.vector(v_full$E),
  Gene = rep(rownames(v_full$E), ncol(v_full$E))
)

# Sample-level summary
sample_expr <- data.frame(
  Sample = colnames(v_full$E),
  Median = apply(v_full$E, 2, median),
  Mean = apply(v_full$E, 2, mean),
  SD = apply(v_full$E, 2, sd),
  Timepoint = targets$Timepoint,
  Leaf_type = targets$Leaf_type,
  Line = targets$Line,
  Batch = targets$Batch
)

p1 <- ggplot(sample_expr, aes(x = Timepoint, y = Median, fill = Leaf_type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = leaf_colors) +
  facet_wrap(~Batch, scales = "free_x") +
  theme_bw(base_size = 14) +
  labs(title = "Median Expression by Timepoint and Batch",
       subtitle = "Voom-normalized log2-CPM",
       y = "Median Expression",
       x = "Timepoint") +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/05_normalization/expression_by_group.png", p1,
       width = 12, height = 6, dpi = 300)
cat("  Saved: 03_results/figures/05_normalization/expression_by_group.png\n")

# ===== SECTION 6: JAG1 EXPRESSION =====

cat("\n========================================\n")
cat("SECTION 6: JAG1 EXPRESSION\n")
cat("========================================\n\n")

# Check JAG1 and JAG2 expression
jag1_id <- PARAMS$JAG1
jag2_id <- PARAMS$JAG2

cat("JAG1 (", jag1_id, "):\n", sep = "")
if (jag1_id %in% rownames(v_full$E)) {
  jag1_expr <- v_full$E[jag1_id, ]
  cat("  Mean expression:", round(mean(jag1_expr), 3), "log2-CPM\n")
  cat("  Range:", round(min(jag1_expr), 3), "to", round(max(jag1_expr), 3), "\n")

  # Expression by leaf type
  cat("  By leaf type:\n")
  cat("    Broad:", round(mean(jag1_expr[targets$Leaf_type == "Broad"]), 3), "\n")
  cat("    Narrow:", round(mean(jag1_expr[targets$Leaf_type == "Narrow"]), 3), "\n")
} else {
  cat("  WARNING: JAG1 not found in expression matrix\n")
}

cat("\nJAG2 (", jag2_id, "):\n", sep = "")
if (jag2_id %in% rownames(v_full$E)) {
  jag2_expr <- v_full$E[jag2_id, ]
  cat("  Mean expression:", round(mean(jag2_expr), 3), "log2-CPM\n")
  cat("  Range:", round(min(jag2_expr), 3), "to", round(max(jag2_expr), 3), "\n")
} else {
  cat("  WARNING: JAG2 not found in expression matrix\n")
}

# --- Plot JAG1/JAG2 expression ---
cat("\nCreating JAG expression plots...\n")

# Prepare JAG data
jag_data <- data.frame(
  Sample = colnames(v_full$E),
  Timepoint = targets$Timepoint,
  Leaf_type = targets$Leaf_type,
  Line = targets$Line,
  Batch = targets$Batch,
  JAG1 = if(jag1_id %in% rownames(v_full$E)) v_full$E[jag1_id, ] else NA,
  JAG2 = if(jag2_id %in% rownames(v_full$E)) v_full$E[jag2_id, ] else NA
)

# JAG1 expression plot
if (!all(is.na(jag_data$JAG1))) {
  p_jag1 <- ggplot(jag_data, aes(x = Timepoint, y = JAG1, color = Line, group = Line)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    facet_wrap(~Leaf_type) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = 14) +
    labs(title = "GmJAG1 Expression Across Development",
         subtitle = paste("Gene:", jag1_id),
         y = "Expression (log2-CPM)",
         x = "Timepoint") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  ggsave("03_results/figures/05_normalization/JAG1_expression.png", p_jag1,
         width = 10, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/05_normalization/JAG1_expression.png\n")
}

# JAG1 vs JAG2 comparison
if (!all(is.na(jag_data$JAG1)) && !all(is.na(jag_data$JAG2))) {
  jag_long <- tidyr::pivot_longer(jag_data,
                                   cols = c(JAG1, JAG2),
                                   names_to = "Gene",
                                   values_to = "Expression")

  p_jag_both <- ggplot(jag_long, aes(x = Timepoint, y = Expression,
                                      color = Gene, linetype = Leaf_type)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    scale_color_manual(values = c("JAG1" = "blue", "JAG2" = "red")) +
    theme_bw(base_size = 14) +
    labs(title = "GmJAG1 vs GmJAG2 Expression",
         y = "Expression (log2-CPM)",
         x = "Timepoint") +
    theme(plot.title = element_text(size = 16, face = "bold"))

  ggsave("03_results/figures/05_normalization/JAG1_vs_JAG2.png", p_jag_both,
         width = 10, height = 6, dpi = 300)
  cat("  Saved: 03_results/figures/05_normalization/JAG1_vs_JAG2.png\n")
}

# ===== SECTION 7: SAMPLE CLUSTERING (POST-NORMALIZATION) =====

cat("\n========================================\n")
cat("SECTION 7: SAMPLE CLUSTERING\n")
cat("========================================\n\n")

cat("Creating MDS plot (voom-normalized)...\n")

png("03_results/figures/05_normalization/MDS_voom.png",
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(2, 2), cex.main = 1.2, cex.lab = 1.1)

# MDS from voom object
plotMDS(v_full, col = leaf_colors[targets$Leaf_type],
        pch = 19, cex = 1.5, main = "MDS (Voom) - Leaf Type")
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 19)

plotMDS(v_full, col = tp_colors[targets$Timepoint],
        pch = 19, cex = 1.5, main = "MDS (Voom) - Timepoint")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.8)

plotMDS(v_full, col = batch_colors[as.character(targets$Batch)],
        pch = 19, cex = 1.5, main = "MDS (Voom) - Batch")
legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19)

# Dim 3 vs 4
plotMDS(v_full, dim.plot = c(3, 4), col = tp_colors[targets$Timepoint],
        pch = 19, cex = 1.5, main = "MDS (Voom) - Dim 3 vs 4")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 0.8)

dev.off()
cat("  Saved: 03_results/figures/05_normalization/MDS_voom.png\n")

# ===== SECTION 8: NORMALIZATION SUMMARY =====

cat("\n========================================\n")
cat("SECTION 8: NORMALIZATION SUMMARY\n")
cat("========================================\n\n")

# Create summary table
norm_summary <- data.frame(
  Metric = c(
    "Total genes",
    "Total samples",
    "Samples (full)",
    "Samples (TP1-4)",
    "Mean library size",
    "Mean norm factor",
    "Expression range (log2-CPM)",
    "JAG1 mean expression",
    "JAG2 mean expression",
    "Groups (design)",
    "Voom weights range"
  ),
  Value = c(
    nrow(v_full$E),
    ncol(v_full$E),
    "60 (all)",
    "48 (batch-free)",
    format(round(mean(lib_sizes)), big.mark = ","),
    round(mean(norm_factors), 4),
    paste(round(min(v_full$E), 2), "to", round(max(v_full$E), 2)),
    if(jag1_id %in% rownames(v_full$E)) round(mean(v_full$E[jag1_id, ]), 3) else "NA",
    if(jag2_id %in% rownames(v_full$E)) round(mean(v_full$E[jag2_id, ]), 3) else "NA",
    ncol(design_full),
    paste(round(min(v_full$weights), 3), "to", round(max(v_full$weights), 3))
  )
)

write.csv(norm_summary, "03_results/tables/normalization_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/normalization_summary.csv\n\n")

print(norm_summary)

# ===== SECTION 9: SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 9: SAVE CHECKPOINT\n")
cat("========================================\n\n")

# Save all objects for downstream analysis
save(
  # Voom objects (PRIMARY for DE analysis)
  v_full,                # Voom object for all 60 samples
  v_tp1_4,               # Voom object for TP1-4 (batch-free)

  # Design matrices
  design_full,           # Design matrix for all samples
  design_tp1_4,          # Design matrix for TP1-4

  # DGEList objects
  dge_corrected,         # Batch-corrected DGEList
  dge_tp1_4,             # TP1-4 DGEList

  # Expression matrices
  logCPM_corrected,      # Log-CPM (batch-corrected)

  # Metadata
  targets,               # Full metadata
  targets_tp1_4,         # TP1-4 metadata

  # Parameters
  PARAMS,
  QC_PARAMS,

  file = "03_results/checkpoints/05_normalized.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/05_normalized.RData\n")
cat("  Contains:\n")
cat("    - Voom objects: v_full, v_tp1_4\n")
cat("    - Design matrices: design_full, design_tp1_4\n")
cat("    - DGEList: dge_corrected, dge_tp1_4\n")
cat("    - Expression: logCPM_corrected\n")
cat("    - Metadata: targets, targets_tp1_4\n")
cat("    - Parameters: PARAMS, QC_PARAMS\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 05: NORMALIZATION - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Voom transformation complete\n")
cat("    - Full dataset: ", ncol(v_full$E), " samples, ", nrow(v_full$E), " genes\n", sep = "")
cat("    - TP1-4 subset: ", ncol(v_tp1_4$E), " samples (batch-free)\n", sep = "")
cat("    - Design matrix groups:", ncol(design_full), "\n")
cat("    - JAG1 detected: ", jag1_id %in% rownames(v_full$E), "\n", sep = "")
cat("\n")
cat("  Data ready for:\n")
cat("    - Differential expression (limma)\n")
cat("    - Mixed-effects modeling\n")
cat("    - Co-expression analysis (WGCNA)\n")
cat("\n")
cat("  Next: Run 06_sample_validation.R\n")
cat("================================================================\n")
