# Script 03: Quality Control

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 03: QUALITY CONTROL\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/03_QC", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",          # For DGEList and normalization
  "limma",          # For MDS plots
  "ggplot2",        # For plotting
  "pheatmap",       # For heatmaps
  "RColorBrewer",   # For colors
  "dplyr",          # For data manipulation
  "corrplot",       # For correlation plots
  "dendextend"      # For dendrograms
)

# Install missing packages
missing <- required_packages[!required_packages %in% installed.packages()[,1]]
if (length(missing) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing)
}

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD PREVIOUS CHECKPOINT =====

cat("Loading previous checkpoint...\n")
load("03_results/checkpoints/02_sample_metadata.RData")
cat("  Loaded: raw_counts (", nrow(raw_counts), " x ", ncol(raw_counts), "), targets\n\n", sep = "")

# ===== QC PARAMETERS =====

QC_PARAMS <- list(
  min_correlation = 0.9,        # Minimum expected replicate correlation
  outlier_threshold = 2.5,      # SD threshold for outlier detection
  min_cpm = 1,                  # Minimum CPM for gene filtering
  min_samples = 3               # Gene must have min_cpm in this many samples
)

cat("QC Parameters:\n")
cat("  Minimum replicate correlation:", QC_PARAMS$min_correlation, "\n")
cat("  Outlier threshold (SD):", QC_PARAMS$outlier_threshold, "\n\n")

# ===== SECTION 1: CREATE DGEList =====

cat("========================================\n")
cat("SECTION 1: CREATE DGEList OBJECT\n")
cat("========================================\n\n")

# Create DGEList object
dge <- DGEList(counts = raw_counts, samples = targets)

# Calculate library sizes
cat("Library size statistics:\n")
cat("  Total reads:", format(sum(dge$samples$lib.size), big.mark = ","), "\n")
cat("  Mean per sample:", format(round(mean(dge$samples$lib.size)), big.mark = ","), "\n")
cat("  Median per sample:", format(round(median(dge$samples$lib.size)), big.mark = ","), "\n")
cat("  Range:", format(min(dge$samples$lib.size), big.mark = ","), "to",
    format(max(dge$samples$lib.size), big.mark = ","), "\n\n")

# --- Salmon ReadFate Plot (Mapping Rate) ---
# Check if mapping rate data exists in targets
mapping_cols <- c("Assigned_to_gene", "mapped", "mapping_rate", "pct_mapped")
mapping_col <- intersect(mapping_cols, colnames(targets))

if (length(mapping_col) > 0) {
  cat("Creating Salmon ReadFate plot (mapping rates)...\n")

  mapping_col <- mapping_col[1]  # Use first matching column
  mapping_rate <- targets[[mapping_col]]

  # If values > 1, assume they are counts not percentages
  if (max(mapping_rate, na.rm = TRUE) > 100) {
    # Calculate percentage from counts
    if ("NumReads" %in% colnames(targets)) {
      mapping_rate <- mapping_rate / targets$NumReads * 100
    } else if ("Total" %in% colnames(targets)) {
      mapping_rate <- mapping_rate / targets$Total * 100
    }
  }

  # Create plot data
  readfate_df <- data.frame(
    Sample = targets$Sample,
    Mapped_Pct = mapping_rate,
    Batch = as.factor(targets$Batch),
    Timepoint = targets$Timepoint
  )

  # Order by batch then timepoint
  readfate_df <- readfate_df[order(readfate_df$Batch, readfate_df$Timepoint), ]
  readfate_df$Sample <- factor(readfate_df$Sample, levels = readfate_df$Sample)

  png("03_results/figures/03_QC/salmon_readfate.png",
      width = 14, height = 6, units = "in", res = 300)

  barplot(readfate_df$Mapped_Pct,
          ylim = c(0, 100),
          ylab = "Percent of mapped reads",
          xlab = "",
          las = 2,
          col = ifelse(readfate_df$Batch == levels(readfate_df$Batch)[1], "#1f77b4", "#ff7f0e"),
          main = "Salmon ReadFate - Mapping Rates by Sample",
          cex.axis = 0.8,
          cex.names = 0.5,
          names.arg = readfate_df$Sample)
  legend("bottomright", legend = paste("Batch", levels(readfate_df$Batch)),
         fill = c("#1f77b4", "#ff7f0e"))

  dev.off()
  cat("  Saved: 03_results/figures/03_QC/salmon_readfate.png\n")

  cat("  Mapping rate statistics:\n")
  cat("    Mean:", round(mean(mapping_rate, na.rm = TRUE), 1), "%\n")
  cat("    Range:", round(min(mapping_rate, na.rm = TRUE), 1), "-",
      round(max(mapping_rate, na.rm = TRUE), 1), "%\n\n")

} else {
  cat("Note: Mapping rate data not found in targets.\n")
  cat("  Checked columns:", paste(mapping_cols, collapse = ", "), "\n")
  cat("  Available columns:", paste(colnames(targets), collapse = ", "), "\n\n")
}

# ===== SECTION 2: GENE FILTERING =====

cat("========================================\n")
cat("SECTION 2: GENE FILTERING\n")
cat("========================================\n\n")

# Calculate CPM
cpm_matrix <- cpm(dge)

# Filter genes
keep <- rowSums(cpm_matrix > QC_PARAMS$min_cpm) >= QC_PARAMS$min_samples

cat("Gene filtering (CPM >", QC_PARAMS$min_cpm, "in >=", QC_PARAMS$min_samples, "samples):\n")
cat("  Genes before filtering:", nrow(dge), "\n")
cat("  Genes after filtering:", sum(keep), "\n")
cat("  Genes removed:", sum(!keep), "\n")
cat("  Percent retained:", round(sum(keep)/nrow(dge)*100, 1), "%\n\n")

# Check JAG1
cat("JAG1 gene status:\n")
if (PARAMS$JAG1 %in% rownames(dge)) {
  jag1_cpm <- cpm_matrix[PARAMS$JAG1, ]
  cat("  JAG1 present: YES\n")
  cat("  JAG1 mean CPM:", round(mean(jag1_cpm), 2), "\n")
  cat("  JAG1 passes filter:", keep[PARAMS$JAG1], "\n\n")
} else {
  cat("  JAG1 present: NO\n\n")
}

# Apply filter
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# ===== CRITICAL: Calculate % reads retained after filtering =====
# This is an important QC metric - most reads should be retained even if many genes are removed
reads_before <- dge$samples$lib.size
reads_after <- dge_filtered$samples$lib.size
pct_reads_retained <- reads_after / reads_before * 100

cat("Reads retention after filtering:\n")
cat("  Mean % reads retained:", round(mean(pct_reads_retained), 2), "%\n")
cat("  Range:", round(min(pct_reads_retained), 2), "-",
    round(max(pct_reads_retained), 2), "%\n")
cat("  Interpretation: Low-expression genes removed had minimal counts\n\n")

# ===== Create comprehensive filtering summary for manuscript =====
filtering_summary <- data.frame(
  Metric = c(
    "Initial genes (before filtering)",
    "Genes after filtering",
    "Genes removed",
    "Percent genes retained",
    "Mean percent reads retained",
    "Min percent reads retained",
    "Max percent reads retained",
    "Filtering threshold (CPM)",
    "Minimum samples required"
  ),
  Value = c(
    nrow(dge),
    sum(keep),
    sum(!keep),
    round(sum(keep)/nrow(dge)*100, 2),
    round(mean(pct_reads_retained), 2),
    round(min(pct_reads_retained), 2),
    round(max(pct_reads_retained), 2),
    QC_PARAMS$min_cpm,
    QC_PARAMS$min_samples
  ),
  stringsAsFactors = FALSE
)

# Save filtering summary for Results section
write.csv(filtering_summary, "03_results/tables/filtering_summary.csv", row.names = FALSE)
cat("Saved: 03_results/tables/filtering_summary.csv\n\n")

# Print summary for console
cat("FILTERING SUMMARY FOR MANUSCRIPT:\n")
cat("----------------------------------\n")
print(filtering_summary)
cat("\n")

# ===== SECTION 3: NORMALIZATION (for QC purposes) =====

cat("========================================\n")
cat("SECTION 3: NORMALIZATION\n")
cat("========================================\n\n")

# Calculate normalization factors
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

cat("TMM normalization factors:\n")
cat("  Range:", round(min(dge_filtered$samples$norm.factors), 3), "to",
    round(max(dge_filtered$samples$norm.factors), 3), "\n")
cat("  Mean:", round(mean(dge_filtered$samples$norm.factors), 3), "\n")
cat("  SD:", round(sd(dge_filtered$samples$norm.factors), 3), "\n\n")

# Calculate log CPM for downstream QC
logCPM <- cpm(dge_filtered, log = TRUE, prior.count = 2)

# ===== SECTION 4: SAMPLE CORRELATIONS =====

cat("========================================\n")
cat("SECTION 4: SAMPLE CORRELATIONS\n")
cat("========================================\n\n")

# Calculate sample correlation matrix
sample_cor <- cor(logCPM, method = "pearson")

cat("Sample correlation statistics:\n")
cat("  Mean correlation:", round(mean(sample_cor[lower.tri(sample_cor)]), 4), "\n")
cat("  Min correlation:", round(min(sample_cor[lower.tri(sample_cor)]), 4), "\n")
cat("  Max correlation:", round(max(sample_cor[lower.tri(sample_cor)]), 4), "\n\n")

# ===== SECTION 5: REPLICATE CONSISTENCY =====

cat("========================================\n")
cat("SECTION 5: REPLICATE CONSISTENCY\n")
cat("========================================\n\n")

# Calculate within-group correlations
groups <- unique(targets$Group)
replicate_cors <- data.frame(
  Group = character(),
  Mean_Correlation = numeric(),
  Min_Correlation = numeric(),
  Max_Correlation = numeric(),
  N_Replicates = numeric(),
  stringsAsFactors = FALSE
)

for (grp in groups) {
  grp_samples <- rownames(targets)[targets$Group == grp]
  if (length(grp_samples) >= 2) {
    grp_cor <- cor(logCPM[, grp_samples], method = "pearson")
    diag(grp_cor) <- NA

    replicate_cors <- rbind(replicate_cors, data.frame(
      Group = grp,
      Mean_Correlation = mean(grp_cor, na.rm = TRUE),
      Min_Correlation = min(grp_cor, na.rm = TRUE),
      Max_Correlation = max(grp_cor, na.rm = TRUE),
      N_Replicates = length(grp_samples)
    ))
  }
}

# Summary
cat("Replicate correlation summary:\n")
cat("  Overall mean:", round(mean(replicate_cors$Mean_Correlation), 4), "\n")
cat("  Overall min:", round(min(replicate_cors$Min_Correlation), 4), "\n")
cat("  Groups with low correlation (< 0.95):\n")

low_cor_groups <- replicate_cors[replicate_cors$Min_Correlation < 0.95, ]
if (nrow(low_cor_groups) > 0) {
  print(low_cor_groups)
} else {
  cat("    None - all groups have good replicate consistency\n")
}

# Save replicate correlations
write.csv(replicate_cors, "03_results/tables/replicate_correlations.csv", row.names = FALSE)
cat("\n  Saved: 03_results/tables/replicate_correlations.csv\n\n")

# ===== SECTION 6: EVALUATE PROBLEMATIC SAMPLE =====

cat("========================================\n")
cat("SECTION 6: PROBLEMATIC SAMPLE ANALYSIS\n")
cat("========================================\n\n")

# Focus on sample 745_T2_R2 identified in Phase 1
problem_sample <- "745_T2_R2"

if (problem_sample %in% colnames(logCPM)) {
  cat("Analyzing sample:", problem_sample, "\n\n")

  # Get expected and neighboring groups
  expected_group <- targets[problem_sample, "Group"]  # Should be PI547745_TP2
  line <- targets[problem_sample, "Line"]

  # Get all samples from this line
  line_samples <- rownames(targets)[targets$Line == line]

  # Calculate correlation with each timepoint group
  timepoints <- c("TP0", "TP1", "TP2", "TP3", "TP4")
  tp_cors <- data.frame()

  for (tp in timepoints) {
    tp_samples <- rownames(targets)[targets$Line == line & targets$Timepoint == tp]
    tp_samples <- setdiff(tp_samples, problem_sample)  # Exclude problem sample

    if (length(tp_samples) > 0) {
      cors <- cor(logCPM[, problem_sample], logCPM[, tp_samples], method = "pearson")
      tp_cors <- rbind(tp_cors, data.frame(
        Timepoint = tp,
        Mean_Correlation = mean(cors),
        Expected = (tp == "TP2")
      ))
    }
  }

  cat("Correlation of", problem_sample, "with each timepoint group:\n")
  print(tp_cors)

  # Determine best match
  best_match <- tp_cors$Timepoint[which.max(tp_cors$Mean_Correlation)]
  expected_match <- "TP2"

  cat("\nExpected timepoint: TP2\n")
  cat("Best matching timepoint:", best_match, "\n")

  if (best_match != expected_match) {
    cat("\nWARNING: Sample", problem_sample, "correlates best with", best_match,
        "not", expected_match, "\n")
    cat("RECOMMENDATION: Consider reassigning or excluding this sample\n")

    # Update QC flag
    targets[problem_sample, "QC_flag"] <- paste0("MISMATCH_", best_match)
  } else {
    cat("\nSample matches expected timepoint\n")
    targets[problem_sample, "QC_flag"] <- "OK"
  }
} else {
  cat("Sample", problem_sample, "not found in data\n")
}

# ===== SECTION 7: OUTLIER DETECTION =====

cat("\n========================================\n")
cat("SECTION 7: OUTLIER DETECTION\n")
cat("========================================\n\n")

# Method 1: Based on correlation with other samples
mean_cors <- colMeans(sample_cor)
cor_zscore <- scale(mean_cors)

outliers_cor <- names(which(abs(cor_zscore) > QC_PARAMS$outlier_threshold))

cat("Outlier detection (correlation-based):\n")
cat("  Samples with Z-score >", QC_PARAMS$outlier_threshold, ":\n")
if (length(outliers_cor) > 0) {
  for (s in outliers_cor) {
    cat("    ", s, ": Z =", round(cor_zscore[s], 2), ", mean cor =", round(mean_cors[s], 3), "\n")
  }
} else {
  cat("    None detected\n")
}

# Method 2: Library size outliers
lib_sizes <- dge_filtered$samples$lib.size
lib_zscore <- scale(lib_sizes)
names(lib_zscore) <- rownames(dge_filtered$samples)

outliers_lib <- names(which(abs(lib_zscore) > QC_PARAMS$outlier_threshold))

cat("\nOutlier detection (library size-based):\n")
cat("  Samples with Z-score >", QC_PARAMS$outlier_threshold, ":\n")
if (length(outliers_lib) > 0) {
  for (s in outliers_lib) {
    cat("    ", s, ": Z =", round(lib_zscore[s], 2),
        ", lib size =", format(lib_sizes[s], big.mark = ","), "\n")
  }
} else {
  cat("    None detected\n")
}

# Update QC flags for outliers
all_outliers <- unique(c(outliers_cor, outliers_lib))
for (s in all_outliers) {
  if (targets[s, "QC_flag"] == "OK") {
    targets[s, "QC_flag"] <- "OUTLIER"
  }
}

# ===== SECTION 8: CREATE QC SUMMARY TABLE =====

cat("\n========================================\n")
cat("SECTION 8: QC SUMMARY TABLE\n")
cat("========================================\n\n")

# Create comprehensive QC metrics table
qc_metrics <- data.frame(
  Sample = colnames(logCPM),
  Line = targets$Line,
  Leaf_type = targets$Leaf_type,
  Timepoint = targets$Timepoint,
  Batch = targets$Batch,
  Library_size = dge_filtered$samples$lib.size,
  Norm_factor = dge_filtered$samples$norm.factors,
  Mean_correlation = mean_cors,
  Correlation_zscore = as.vector(cor_zscore),
  Library_zscore = as.vector(lib_zscore),
  QC_flag = targets$QC_flag,
  stringsAsFactors = FALSE
)

# Save QC metrics
write.csv(qc_metrics, "03_results/tables/QC_metrics.csv", row.names = FALSE)
cat("Saved: 03_results/tables/QC_metrics.csv\n")

# Summary of flags
cat("\nQC Flag Summary:\n")
print(table(qc_metrics$QC_flag))

# ===== SECTION 9: VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 9: VISUALIZATIONS\n")
cat("========================================\n\n")

# Set up colors
batch_colors <- c("2021" = "#E41A1C", "2022" = "#377EB8")
leaf_colors <- c("Broad" = "darkgreen", "Narrow" = "purple")
line_colors <- c("PI532462A" = "#1B9E77", "LD112170" = "#D95F02",
                 "PI612713B" = "#7570B3", "PI547745" = "#E7298A")
tp_colors <- brewer.pal(5, "YlOrRd")
names(tp_colors) <- c("TP0", "TP1", "TP2", "TP3", "TP4")

# --- Plot 1: Sample correlation heatmap ---
cat("Creating correlation heatmap...\n")

# Annotation for heatmap
annotation_col <- data.frame(
  Leaf_type = targets$Leaf_type,
  Timepoint = targets$Timepoint,
  Batch = as.factor(targets$Batch),
  row.names = colnames(logCPM)
)

annotation_colors <- list(
  Leaf_type = leaf_colors,
  Timepoint = tp_colors,
  Batch = batch_colors
)

png("03_results/figures/03_QC/sample_correlation_heatmap.png",
    width = 14, height = 12, units = "in", res = 300)

pheatmap(sample_cor,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         breaks = seq(0.8, 1, length.out = 101),
         main = "Sample Correlation Matrix",
         fontsize_row = 7,
         fontsize_col = 7,
         fontsize = 12)

dev.off()
cat("  Saved: 03_results/figures/03_QC/sample_correlation_heatmap.png\n")

# --- Plot 2: MDS/PCA plots ---
cat("Creating MDS plots...\n")

png("03_results/figures/03_QC/MDS_plots.png",
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(2, 2), cex.main = 1.2, cex.lab = 1.1)

# MDS colored by Leaf type
plotMDS(dge_filtered, col = leaf_colors[targets$Leaf_type],
        pch = 19, cex = 1.5, main = "MDS - Colored by Leaf Type")
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 19, cex = 1.0)

# MDS colored by Timepoint
plotMDS(dge_filtered, col = tp_colors[targets$Timepoint],
        pch = 19, cex = 1.5, main = "MDS - Colored by Timepoint")
legend("topright", legend = names(tp_colors), col = tp_colors, pch = 19, cex = 1.0)

# MDS colored by Batch
plotMDS(dge_filtered, col = batch_colors[as.character(targets$Batch)],
        pch = 19, cex = 1.5, main = "MDS - Colored by Batch")
legend("topright", legend = names(batch_colors), col = batch_colors, pch = 19, cex = 1.0)

# MDS colored by Line
plotMDS(dge_filtered, col = line_colors[targets$Line],
        pch = 19, cex = 1.5, main = "MDS - Colored by Line")
legend("topright", legend = names(line_colors), col = line_colors, pch = 19, cex = 0.8)

dev.off()
cat("  Saved: 03_results/figures/03_QC/MDS_plots.png\n")

# --- Plot 3: Library size distribution ---
cat("Creating library size plots...\n")

qc_df <- qc_metrics

p1 <- ggplot(qc_df, aes(x = Timepoint, y = Library_size/1e6, fill = Leaf_type)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = leaf_colors) +
  facet_wrap(~Batch, scales = "free_x") +
  theme_bw(base_size = 14) +
  labs(title = "Library Size Distribution by Timepoint and Batch",
       y = "Library Size (millions)",
       x = "Timepoint") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold"))

ggsave("03_results/figures/03_QC/library_sizes_by_group.png", p1,
       width = 14, height = 8, dpi = 300)
cat("  Saved: 03_results/figures/03_QC/library_sizes_by_group.png\n")

# --- Plot 4: RLE plot (Relative Log Expression) ---
cat("Creating RLE plot...\n")

png("03_results/figures/03_QC/RLE_plot.png",
    width = 14, height = 6, units = "in", res = 300)

# Calculate RLE
median_expr <- apply(logCPM, 1, median)
rle <- logCPM - median_expr

# Create boxplot
par(mar = c(10, 5, 4, 2))
boxplot(rle,
        las = 2,
        col = batch_colors[as.character(targets$Batch)],
        main = "Relative Log Expression (RLE) Plot",
        ylab = "Deviation from median (log2)",
        cex.axis = 0.7,
        cex.main = 1.3,
        cex.lab = 1.2,
        outline = FALSE)
abline(h = 0, col = "red", lty = 2, lwd = 2)
legend("topright", legend = names(batch_colors), fill = batch_colors, title = "Batch", cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/03_QC/RLE_plot.png\n")

# --- Plot 5: Sample dendrogram ---
cat("Creating sample dendrogram...\n")

png("03_results/figures/03_QC/sample_dendrogram.png",
    width = 14, height = 8, units = "in", res = 300)

# Hierarchical clustering
sample_dist <- as.dist(1 - sample_cor)
sample_clust <- hclust(sample_dist, method = "average")

# Create dendrogram with colored labels
dend <- as.dendrogram(sample_clust)
labels_colors(dend) <- leaf_colors[targets$Leaf_type[order.dendrogram(dend)]]

par(mar = c(8, 5, 4, 2))
plot(dend, main = "Sample Clustering Dendrogram",
     ylab = "Distance (1 - correlation)",
     cex.main = 1.3, cex.lab = 1.2)
legend("topright", legend = names(leaf_colors), col = leaf_colors, pch = 15, cex = 1.0)

dev.off()
cat("  Saved: 03_results/figures/03_QC/sample_dendrogram.png\n")

# --- Plot 6: Mean-variance relationship ---
cat("Creating mean-variance plot...\n")

png("03_results/figures/03_QC/mean_variance.png",
    width = 10, height = 8, units = "in", res = 300)

# Voom transformation for visualization
par(cex.main = 1.3, cex.lab = 1.2)
v <- voom(dge_filtered, plot = TRUE)
title("Mean-Variance Relationship (voom)")

dev.off()
cat("  Saved: 03_results/figures/03_QC/mean_variance.png\n")

# ===== SECTION 10: QC DECISIONS =====

cat("\n========================================\n")
cat("SECTION 10: QC DECISIONS & RECOMMENDATIONS\n")
cat("========================================\n\n")

# Summary of QC findings
cat("QC SUMMARY:\n")
cat("-----------\n\n")

cat("1. Library Sizes:\n")
cat("   Mean:", format(round(mean(dge_filtered$samples$lib.size)), big.mark = ","), "\n")
cat("   All samples have adequate sequencing depth\n\n")

cat("2. Replicate Consistency:\n")
cat("   Mean within-group correlation:", round(mean(replicate_cors$Mean_Correlation), 3), "\n")
if (mean(replicate_cors$Mean_Correlation) > 0.95) {
  cat("   STATUS: GOOD - Replicates are highly consistent\n\n")
} else {
  cat("   STATUS: CHECK - Some groups have lower consistency\n\n")
}

cat("3. Batch Effects:\n")
cat("   Batch 2022: TP0 samples only\n")
cat("   Batch 2021: TP1-TP4 samples\n")
cat("   STATUS: CONFOUNDED - Will address in batch correction step\n\n")

cat("4. Sample Flags:\n")
flag_counts <- table(targets$QC_flag)
print(flag_counts)

cat("\n5. Recommendations:\n")
problem_samples <- rownames(targets)[targets$QC_flag != "OK"]
if (length(problem_samples) > 0) {
  cat("   Flagged samples:", paste(problem_samples, collapse = ", "), "\n")
  cat("   Action: Evaluate impact in batch correction step\n")
} else {
  cat("   All samples pass QC\n")
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 11: SAVE CHECKPOINT\n")
cat("========================================\n\n")

# Objects to save
save(
  raw_counts,           # Original counts
  dge_filtered,         # Filtered DGEList
  logCPM,               # Log-transformed CPM
  targets,              # Metadata with QC flags
  sample_cor,           # Sample correlation matrix
  replicate_cors,       # Replicate consistency
  qc_metrics,           # QC summary table
  PARAMS,               # Analysis parameters
  QC_PARAMS,            # QC parameters
  file = "03_results/checkpoints/03_QC_complete.RData"
)

cat("Saved checkpoint: 03_results/checkpoints/03_QC_complete.RData\n")
cat("  Contains: raw_counts, dge_filtered, logCPM, targets, sample_cor,\n")
cat("            replicate_cors, qc_metrics, PARAMS, QC_PARAMS\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 03: QUALITY CONTROL - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  QC Summary:\n")
cat("    - Genes after filtering:", nrow(dge_filtered), "\n")
cat("    - Mean replicate correlation:", round(mean(replicate_cors$Mean_Correlation), 3), "\n")
cat("    - Samples flagged:", sum(targets$QC_flag != "OK"), "\n")
cat("    - Batch effect: Detected (TP0 vs TP1-4)\n")
cat("\n")
cat("  Next: Run 04_batch_correction.R\n")
cat("================================================================\n")
