# Script 18: WGCNA Network Construction

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 18: WGCNA NETWORK CONSTRUCTION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
#base_dir <-"C:/Users/bgtam/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/"

setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/WGCNA", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/18_WGCNA_network", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/WGCNA", recursive = TRUE, showWarnings = FALSE)  # For TOM files

# Track overall timing
script_start_time <- Sys.time()

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "WGCNA",
  "edgeR",
  "limma",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer"
)

# Check and install WGCNA if needed
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  cat("Installing WGCNA and dependencies...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("impute")
  BiocManager::install("preprocessCore")
  install.packages("WGCNA")
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Enable multi-threading
enableWGCNAThreads()
cat("  Packages loaded\n")
cat("  WGCNA threads enabled\n\n")

# ===== LOAD CHECKPOINTS =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded validated expression data\n")

load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n\n")

# Extract expression matrix (log2 CPM)
logCPM <- v_primary$E
targets <- targets_primary

cat("Expression matrix dimensions:", nrow(logCPM), "genes x", ncol(logCPM), "samples\n\n")

# ===== GENE SELECTION FOR WGCNA =====

cat("========================================\n")
cat("SECTION 2: GENE SELECTION\n")
cat("========================================\n\n")

# Strategy: Use ANOVA to select genes with significant variation across conditions
# This reduces noise and focuses on biologically relevant genes

cat("Performing ANOVA for gene selection...\n")
cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n")
anova_start_time <- Sys.time()

# Create design matrix for ANOVA (test all group differences)
group <- factor(paste(targets$Leaf_type, targets$Timepoint, sep = "_"))
design_anova <- model.matrix(~0 + group)
colnames(design_anova) <- levels(group)

# Fit model
y_anova <- new("EList", list(E = logCPM))
fit_anova <- lmFit(y_anova, design_anova)

# F-test for any difference between groups
contrasts_anova <- makeContrasts(
  # Test if any group differs from overall mean
  # Using contrasts between all groups
  Broad_TP0 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Broad_TP1 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Broad_TP2 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Broad_TP3 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Broad_TP4 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Narrow_TP0 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Narrow_TP1 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Narrow_TP2 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Narrow_TP3 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  Narrow_TP4 - (Broad_TP0 + Broad_TP1 + Broad_TP2 + Broad_TP3 + Broad_TP4 + Narrow_TP0 + Narrow_TP1 + Narrow_TP2 + Narrow_TP3 + Narrow_TP4)/10,
  levels = design_anova
)

fit_anova2 <- contrasts.fit(fit_anova, contrasts_anova)
fit_anova2 <- eBayes(fit_anova2)

anova_end_time <- Sys.time()
cat("  ANOVA completed in", round(difftime(anova_end_time, anova_start_time, units = "secs"), 1), "seconds\n")

# Get F-statistic results
f_results <- topTable(fit_anova2, number = Inf, sort.by = "none")

# Select genes based on FDR threshold
fdr_threshold <- 0.01  # Strict threshold for WGCNA

genes_for_wgcna <- rownames(f_results)[f_results$adj.P.Val < fdr_threshold]

cat("\nGene selection summary:\n")
cat("  Total genes tested:", nrow(f_results), "\n")
cat("  Genes passing FDR <", fdr_threshold, ":", length(genes_for_wgcna), "\n")
cat("  Proportion selected:", round(length(genes_for_wgcna)/nrow(f_results) * 100, 1), "%\n")

# Plot FDR distribution
png("03_results/figures/18_WGCNA_network/FDR_distribution.png",
    width = 800, height = 600, res = 120)
hist(f_results$adj.P.Val, breaks = 100,
     main = "FDR Distribution - ANOVA for Gene Selection",
     xlab = "Adjusted P-value (FDR)",
     col = "lightblue", border = "white")
abline(v = c(0.01, 0.05, 0.1), col = c("red", "orange", "green"), lty = 2, lwd = 2)
legend("topright", legend = c("FDR = 0.01", "FDR = 0.05", "FDR = 0.10"),
       col = c("red", "orange", "green"), lty = 2, lwd = 2)
dev.off()
cat("  Saved: FDR_distribution.png\n")

# Ensure JAG1 targets are included
# target_table has Confidence_Tier column - filter for actual targets (not "Not_Target")
JAG1_targets_genes <- unique(target_table$GeneID[target_table$Confidence_Tier != "Not_Target"])
cat("JAG1 target genes to ensure inclusion:", length(JAG1_targets_genes), "\n")
extra_genes <- setdiff(JAG1_targets_genes, genes_for_wgcna)
if (length(extra_genes) > 0) {
  cat("\nAdding", length(extra_genes), "JAG1 target genes not in ANOVA selection\n")
  genes_for_wgcna <- unique(c(genes_for_wgcna, extra_genes))
}

cat("\nFinal gene count for WGCNA:", length(genes_for_wgcna), "\n")

# ===== PREPARE EXPRESSION MATRIX =====

cat("\n========================================\n")
cat("SECTION 3: PREPARE EXPRESSION DATA\n")
cat("========================================\n")
cat("  Time elapsed:", round(difftime(Sys.time(), script_start_time, units = "mins"), 1), "minutes\n\n")

# Subset expression data for selected genes
# WGCNA expects samples in rows, genes in columns
datExpr <- t(logCPM[genes_for_wgcna, ])

cat("Expression matrix for WGCNA:\n")
cat("  Samples (rows):", nrow(datExpr), "\n")
cat("  Genes (columns):", ncol(datExpr), "\n")

# Check for missing values
na_count <- sum(is.na(datExpr))
if (na_count > 0) {
  cat("\nWarning:", na_count, "missing values detected\n")
  # Remove genes with missing values
  good_genes <- apply(datExpr, 2, function(x) sum(is.na(x)) == 0)
  datExpr <- datExpr[, good_genes]
  cat("  Removed genes with NAs. Remaining:", ncol(datExpr), "genes\n")
}

# Check for zero variance genes
gene_vars <- apply(datExpr, 2, var)
zero_var <- sum(gene_vars == 0)
if (zero_var > 0) {
  cat("\nWarning:", zero_var, "zero-variance genes detected\n")
  datExpr <- datExpr[, gene_vars > 0]
  cat("  Removed zero-variance genes. Remaining:", ncol(datExpr), "genes\n")
}

# Check sample clustering for outliers
cat("\nChecking for sample outliers...\n")
sampleTree <- hclust(dist(datExpr), method = "average")

png("03_results/figures/18_WGCNA_network/sample_dendrogram.png",
    width = 1200, height = 600, res = 120)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers",
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Add a cut line at height that would remove outliers (if any)
abline(h = 100, col = "red", lty = 2)
dev.off()
cat("  Saved: sample_dendrogram.png\n")

# ===== SOFT THRESHOLD SELECTION =====

cat("\n========================================\n")
cat("SECTION 4: SOFT THRESHOLD SELECTION\n")
cat("========================================\n")
cat("  Time elapsed:", round(difftime(Sys.time(), script_start_time, units = "mins"), 1), "minutes\n\n")

cat("Calculating soft threshold power...\n")
cat("(This may take a few minutes...)\n")
cat("  Started at:", format(Sys.time(), "%H:%M:%S"), "\n\n")
sft_start_time <- Sys.time()

# Define powers to test
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))

# Calculate soft threshold
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed hybrid",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1),
  verbose = 3
)

sft_end_time <- Sys.time()
cat("\n  Soft threshold calculation completed at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  Duration:", round(difftime(sft_end_time, sft_start_time, units = "mins"), 1), "minutes\n\n")

# Determine recommended power
# Look for power where R^2 > 0.8 and mean connectivity is reasonable
sft_df <- sft$fitIndices
r2_threshold <- 0.80

# Find powers meeting R^2 criterion
valid_powers <- sft_df$Power[sft_df$SFT.R.sq > r2_threshold]

if (length(valid_powers) > 0) {
  # Choose the lowest power meeting criterion
  recommended_power <- min(valid_powers)
} else {
  # If no power meets criterion, use default of 9 (common for signed networks)
  recommended_power <- 9
  cat("Note: No power achieves R^2 >", r2_threshold, "\n")
  cat("Using default power of 9 for signed hybrid network\n")
}

cat("\nSoft threshold analysis results:\n")
cat("  Recommended power:", recommended_power, "\n")
cat("  Scale-free R^2 at this power:",
    round(sft_df$SFT.R.sq[sft_df$Power == recommended_power], 3), "\n")
cat("  Mean connectivity:",
    round(sft_df$mean.k.[sft_df$Power == recommended_power], 1), "\n")

# Plot soft threshold results
png("03_results/figures/18_WGCNA_network/soft_threshold_selection.png",
    width = 1000, height = 500, res = 120)
par(mfrow = c(1, 2))

# Scale-free topology fit
plot(sft_df$Power, -sign(sft_df$slope) * sft_df$SFT.R.sq,
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     main = "Scale Independence",
     type = "n")
text(sft_df$Power, -sign(sft_df$slope) * sft_df$SFT.R.sq,
     labels = powers, col = "red", cex = 0.9)
abline(h = 0.80, col = "red", lty = 2)
abline(v = recommended_power, col = "blue", lty = 2)

# Mean connectivity
plot(sft_df$Power, sft_df$mean.k.,
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean Connectivity",
     type = "n")
text(sft_df$Power, sft_df$mean.k.,
     labels = powers, col = "red", cex = 0.9)
abline(v = recommended_power, col = "blue", lty = 2)

dev.off()
cat("  Saved: soft_threshold_selection.png\n")

# ===== NETWORK CONSTRUCTION =====

cat("\n========================================\n")
cat("SECTION 5: NETWORK CONSTRUCTION\n")
cat("========================================\n")
cat("  Time elapsed:", round(difftime(Sys.time(), script_start_time, units = "mins"), 1), "minutes\n\n")

cat("Constructing gene co-expression network...\n")
cat("Using signed hybrid network with bicor correlation\n")
cat("----------------------------------------------\n")
cat("  PROGRESS TRACKING:\n")
cat("  - Step 1: Calculate correlations (may take 5-10 min)\n")
cat("  - Step 2: Build TOM matrix (may take 5-15 min)\n")
cat("  - Step 3: Hierarchical clustering\n")
cat("  - Step 4: Module detection and merging\n")
cat("----------------------------------------------\n")
cat("  Network construction started at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  Processing", ncol(datExpr), "genes across", nrow(datExpr), "samples\n")
cat("  Estimated time: 10-30 minutes\n\n")
net_start_time <- Sys.time()

# Set maximum block size based on available memory
# For ~20,000 genes, one block should work
max_block <- min(ncol(datExpr), 30000)

# Construct network and detect modules
net <- blockwiseModules(
  datExpr,
  power = recommended_power,
  maxBlockSize = max_block,
  TOMType = "signed",
  networkType = "signed hybrid",
  minModuleSize = 30,        # Minimum genes per module
  reassignThreshold = 0,
  mergeCutHeight = 0.20,     # Merge similar modules
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "03_results/WGCNA/TOM",
  corType = "bicor",
  maxPOutliers = 0.1,
  deepSplit = 2,             # Module detection sensitivity
  verbose = 3
)

net_end_time <- Sys.time()
net_duration <- difftime(net_end_time, net_start_time, units = "mins")

cat("\n\n========================================\n")
cat("  NETWORK CONSTRUCTION COMPLETE!\n")
cat("========================================\n")
cat("  Finished at:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  Duration:", round(net_duration, 1), "minutes\n")
cat("========================================\n\n")

# Module summary
module_colors <- labels2colors(net$colors)
n_modules <- max(net$colors)
module_sizes <- table(net$colors)

cat("\nModule summary:\n")
cat("  Total modules detected:", n_modules, "(excluding grey/unassigned)\n")
cat("  Genes in modules:", sum(net$colors != 0), "\n")
cat("  Unassigned genes (grey):", sum(net$colors == 0), "\n")

# Print module sizes
cat("\nModule sizes (top 10):\n")
module_size_df <- data.frame(
  Module = names(sort(module_sizes, decreasing = TRUE)),
  Size = as.vector(sort(module_sizes, decreasing = TRUE))
)
module_size_df$Color <- labels2colors(as.numeric(module_size_df$Module))
print(head(module_size_df, 10))

# Plot gene dendrogram with module colors
png("03_results/figures/18_WGCNA_network/gene_dendrogram_modules.png",
    width = 1200, height = 800, res = 120)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Dendrogram and Module Colors"
)
dev.off()
cat("\nSaved: gene_dendrogram_modules.png\n")

# ===== SAVE GENE SELECTION SUMMARY =====

gene_selection_summary <- data.frame(
  Metric = c(
    "Total genes in expression matrix",
    "Genes passing ANOVA FDR threshold",
    "FDR threshold used",
    "JAG1 targets added",
    "Final genes for WGCNA",
    "Genes after QC filtering",
    "Soft threshold power",
    "Scale-free R^2",
    "Mean connectivity",
    "Modules detected",
    "Genes in modules",
    "Unassigned genes"
  ),
  Value = c(
    nrow(logCPM),
    length(rownames(f_results)[f_results$adj.P.Val < fdr_threshold]),
    fdr_threshold,
    length(extra_genes),
    length(genes_for_wgcna),
    ncol(datExpr),
    recommended_power,
    round(sft_df$SFT.R.sq[sft_df$Power == recommended_power], 3),
    round(sft_df$mean.k.[sft_df$Power == recommended_power], 1),
    n_modules,
    sum(net$colors != 0),
    sum(net$colors == 0)
  )
)

write.csv(gene_selection_summary,
          "03_results/tables/WGCNA/gene_selection_summary.csv",
          row.names = FALSE)
cat("Saved: gene_selection_summary.csv\n")

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 6: SAVE CHECKPOINT\n")
cat("========================================\n")
cat("  Time elapsed:", round(difftime(Sys.time(), script_start_time, units = "mins"), 1), "minutes\n\n")

wgcna_prep <- list(
  genes_selected = colnames(datExpr),
  n_genes = ncol(datExpr),
  soft_threshold = recommended_power,
  sft_results = sft,
  network_type = "signed hybrid",
  correlation_type = "bicor"
)

save(
  datExpr,
  net,
  wgcna_prep,
  recommended_power,
  module_colors,
  f_results,
  targets,
  logCPM,
  file = "03_results/checkpoints/18_WGCNA_prep.RData"
)

cat("Checkpoint saved: 18_WGCNA_prep.RData\n")

# ===== SUMMARY =====

script_end_time <- Sys.time()
total_duration <- difftime(script_end_time, script_start_time, units = "mins")

cat("\n================================================================\n")
cat("  SCRIPT 18 COMPLETE: WGCNA NETWORK CONSTRUCTION\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  TOTAL RUNTIME:", round(total_duration, 1), "minutes\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat("  - Genes selected for WGCNA:", ncol(datExpr), "\n")
cat("  - Soft threshold power:", recommended_power, "\n")
cat("  - Modules detected:", n_modules, "\n")
cat("  - Largest module:", max(module_sizes), "genes\n")
cat("  - Unassigned genes:", sum(net$colors == 0), "\n\n")

cat("OUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/18_WGCNA_prep.RData\n")
cat("  - Tables: 03_results/tables/WGCNA/gene_selection_summary.csv\n")
cat("  - Figures: 03_results/figures/18_WGCNA_network/\n")
cat("    - FDR_distribution.png\n")
cat("    - sample_dendrogram.png\n")
cat("    - soft_threshold_selection.png\n")
cat("    - gene_dendrogram_modules.png\n\n")

cat("NEXT STEP: Run Script 19 for module eigengene analysis\n\n")
