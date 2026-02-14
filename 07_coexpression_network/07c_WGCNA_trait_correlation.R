# Script 20: Module-Trait Correlation Analysis

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 20: MODULE-TRAIT CORRELATION ANALYSIS\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/figures/20_WGCNA_traits", recursive = TRUE, showWarnings = FALSE)

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

load("03_results/checkpoints/19_WGCNA_eigengenes.RData")
cat("Loaded WGCNA eigengene data\n")
cat("  Modules:", ncol(MEs), "\n")
cat("  Samples:", nrow(MEs), "\n\n")

# ===== LOAD PHENOTYPE DATA =====

cat("========================================\n")
cat("SECTION 2: LOAD LEAF PHENOTYPE DATA\n")
cat("========================================\n\n")

phenotype_file <- "01_data/soybean_leaf_data_long_format.csv"

if (file.exists(phenotype_file)) {
  pheno_data <- read.csv(phenotype_file, stringsAsFactors = FALSE)
  cat("Loaded phenotype data:", nrow(pheno_data), "measurements\n")
  cat("  Columns:", paste(colnames(pheno_data), collapse = ", "), "\n")
  cat("  Leaf types:", paste(unique(pheno_data$Leaf), collapse = ", "), "\n")
  cat("  Leaf positions:", paste(unique(pheno_data$Leaf_Type), collapse = ", "), "\n")
  cat("  Measurements:", paste(unique(pheno_data$Measurement), collapse = ", "), "\n")
  cat("  Day range:", min(pheno_data$Day), "-", max(pheno_data$Day), "\n")
  cat("  Replicates:", length(unique(pheno_data$Rep)), "\n\n")

  # Convert to wide format for easier processing
  pheno_wide <- reshape(pheno_data,
                        idvar = c("Leaf", "Rep", "Day", "Leaf_Type"),
                        timevar = "Measurement",
                        direction = "wide")
  colnames(pheno_wide) <- gsub("value.", "", colnames(pheno_wide))

  # Focus on V1 and V2 leaves (where Broad vs Narrow differences occur)
  # VC leaves (cotyledons) don't show the phenotypic difference
  pheno_v1v2 <- pheno_wide[pheno_wide$Leaf_Type %in% c("V1", "V2"), ]

  cat("V1/V2 leaf measurements:", nrow(pheno_v1v2), "observations\n\n")

  # ===== CALCULATE PHENOTYPE SUMMARY STATISTICS =====

  cat("========================================\n")
  cat("SECTION 2B: PHENOTYPE SUMMARY\n")
  cat("========================================\n\n")

  # Calculate mature leaf characteristics (using measurements from later days)
  # V1 leaves - use days when V1 is fully developed (later measurements)
  v1_mature <- pheno_v1v2[pheno_v1v2$Leaf_Type == "V1" & !is.na(pheno_v1v2$Ratio), ]

  # Get the last available measurement for each rep (mature leaf)
  v1_mature_final <- do.call(rbind, lapply(split(v1_mature,
                                                  paste(v1_mature$Leaf, v1_mature$Rep)),
                                            function(x) x[which.max(x$Day), ]))

  # Summary by leaf type
  v1_summary <- aggregate(cbind(Length, Width, Ratio) ~ Leaf,
                          data = v1_mature_final,
                          FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                              sd = sd(x, na.rm = TRUE)))

  # Flatten the summary
  v1_summary_flat <- data.frame(
    Leaf = v1_summary$Leaf,
    V1_Length_mean = v1_summary$Length[, "mean"],
    V1_Length_sd = v1_summary$Length[, "sd"],
    V1_Width_mean = v1_summary$Width[, "mean"],
    V1_Width_sd = v1_summary$Width[, "sd"],
    V1_Ratio_mean = v1_summary$Ratio[, "mean"],
    V1_Ratio_sd = v1_summary$Ratio[, "sd"]
  )

  cat("V1 Leaf Summary (mature):\n")
  print(v1_summary_flat)
  cat("\n")

  # Calculate V2 summary if available
  v2_mature <- pheno_v1v2[pheno_v1v2$Leaf_Type == "V2" & !is.na(pheno_v1v2$Ratio), ]

  if (nrow(v2_mature) > 0) {
    v2_mature_final <- do.call(rbind, lapply(split(v2_mature,
                                                    paste(v2_mature$Leaf, v2_mature$Rep)),
                                              function(x) x[which.max(x$Day), ]))

    v2_summary <- aggregate(cbind(Length, Width, Ratio) ~ Leaf,
                            data = v2_mature_final,
                            FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                sd = sd(x, na.rm = TRUE)))

    v2_summary_flat <- data.frame(
      Leaf = v2_summary$Leaf,
      V2_Length_mean = v2_summary$Length[, "mean"],
      V2_Length_sd = v2_summary$Length[, "sd"],
      V2_Width_mean = v2_summary$Width[, "mean"],
      V2_Width_sd = v2_summary$Width[, "sd"],
      V2_Ratio_mean = v2_summary$Ratio[, "mean"],
      V2_Ratio_sd = v2_summary$Ratio[, "sd"]
    )

    cat("V2 Leaf Summary (mature):\n")
    print(v2_summary_flat)
    cat("\n")
  }

  # Combine into phenotype summary
  phenotype_summary <- v1_summary_flat
  if (exists("v2_summary_flat")) {
    phenotype_summary <- merge(phenotype_summary, v2_summary_flat, by = "Leaf")
  }

  # Calculate difference between Broad and Narrow
  broad_idx <- which(phenotype_summary$Leaf == "Broad")
  narrow_idx <- which(phenotype_summary$Leaf == "Narrow")

  cat("KEY PHENOTYPE DIFFERENCES:\n")
  cat("  V1 Ratio: Broad =", round(phenotype_summary$V1_Ratio_mean[broad_idx], 2),
      "vs Narrow =", round(phenotype_summary$V1_Ratio_mean[narrow_idx], 2), "\n")
  if (exists("v2_summary_flat")) {
    cat("  V2 Ratio: Broad =", round(phenotype_summary$V2_Ratio_mean[broad_idx], 2),
        "vs Narrow =", round(phenotype_summary$V2_Ratio_mean[narrow_idx], 2), "\n")
  }
  cat("\n")

  # Save phenotype summary
  write.csv(phenotype_summary, "03_results/tables/WGCNA/phenotype_summary.csv", row.names = FALSE)
  cat("Saved: phenotype_summary.csv\n\n")

  phenotype_available <- TRUE

} else {
  cat("Phenotype data file not found:", phenotype_file, "\n")
  cat("Continuing with categorical traits only...\n\n")
  phenotype_available <- FALSE
}

# ===== CREATE TRAIT MATRIX =====

cat("========================================\n")
cat("SECTION 3: CREATE TRAIT MATRIX\n")
cat("========================================\n\n")

# Create binary trait variables
traitData <- data.frame(
  # Leaf type (Narrow = 1, Broad = 0)
  Narrow = as.numeric(targets$Leaf_type == "Narrow"),
  Broad = as.numeric(targets$Leaf_type == "Broad"),

  # Timepoints (binary indicators)
  TP0 = as.numeric(targets$Timepoint == "TP0"),
  TP1 = as.numeric(targets$Timepoint == "TP1"),
  TP2 = as.numeric(targets$Timepoint == "TP2"),
  TP3 = as.numeric(targets$Timepoint == "TP3"),
  TP4 = as.numeric(targets$Timepoint == "TP4"),

  row.names = rownames(targets)
)

# Add interaction terms (Narrow at specific timepoints)
traitData$Narrow_TP0 <- as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP0")
traitData$Narrow_TP1 <- as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP1")
traitData$Narrow_TP2 <- as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP2")
traitData$Narrow_TP3 <- as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP3")
traitData$Narrow_TP4 <- as.numeric(targets$Leaf_type == "Narrow" & targets$Timepoint == "TP4")

# Broad at specific timepoints
traitData$Broad_TP0 <- as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP0")
traitData$Broad_TP1 <- as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP1")
traitData$Broad_TP2 <- as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP2")
traitData$Broad_TP3 <- as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP3")
traitData$Broad_TP4 <- as.numeric(targets$Leaf_type == "Broad" & targets$Timepoint == "TP4")

# ===== ADD LINE-SPECIFIC L:W RATIO (ACTUAL PHENOTYPE) =====
# These are the actual measured leaf length:width ratios for each line
# Narrow lines: PI612713B = 2.36, PI547745 = 2.16
# Broad lines: LD112170 = 1.59, PI532462A = 1.55

cat("Adding line-specific L:W ratios (actual phenotype measurements)...\n")

line_lw_ratios <- c(
  "PI612713B" = 2.36,  # Narrow
  "PI547745" = 2.16,   # Narrow
  "LD112170" = 1.59,   # Broad
  "PI532462A" = 1.55   # Broad
)

traitData$Line_LW_Ratio <- line_lw_ratios[as.character(targets$Line)]

cat("  Line-specific L:W ratios assigned:\n")
cat("    PI612713B (Narrow): 2.36\n")
cat("    PI547745 (Narrow): 2.16\n")
cat("    LD112170 (Broad): 1.59\n")
cat("    PI532462A (Broad): 1.55\n\n")

# ===== ADD PHENOTYPE MEASUREMENTS AS CONTINUOUS TRAITS =====

if (phenotype_available) {
  cat("Adding continuous phenotype traits...\n\n")

  # Assign V1 leaf ratio to samples based on leaf type
  # This gives each sample a continuous value representing its expected mature leaf shape
  traitData$V1_Ratio <- ifelse(targets$Leaf_type == "Broad",
                                phenotype_summary$V1_Ratio_mean[broad_idx],
                                phenotype_summary$V1_Ratio_mean[narrow_idx])

  traitData$V1_Length <- ifelse(targets$Leaf_type == "Broad",
                                 phenotype_summary$V1_Length_mean[broad_idx],
                                 phenotype_summary$V1_Length_mean[narrow_idx])

  traitData$V1_Width <- ifelse(targets$Leaf_type == "Broad",
                                phenotype_summary$V1_Width_mean[broad_idx],
                                phenotype_summary$V1_Width_mean[narrow_idx])

  # Add V2 traits if available
  if (exists("v2_summary_flat")) {
    traitData$V2_Ratio <- ifelse(targets$Leaf_type == "Broad",
                                  phenotype_summary$V2_Ratio_mean[broad_idx],
                                  phenotype_summary$V2_Ratio_mean[narrow_idx])

    traitData$V2_Length <- ifelse(targets$Leaf_type == "Broad",
                                   phenotype_summary$V2_Length_mean[broad_idx],
                                   phenotype_summary$V2_Length_mean[narrow_idx])

    traitData$V2_Width <- ifelse(targets$Leaf_type == "Broad",
                                  phenotype_summary$V2_Width_mean[broad_idx],
                                  phenotype_summary$V2_Width_mean[narrow_idx])
  }

  # Create developmental stage-specific phenotype traits
  # These capture the interaction between timepoint and phenotype
  # At early timepoints (TP0, TP1), leaves haven't developed yet
  # At later timepoints (TP2+), V1 is developed
  # V2 develops even later

  # V1 phenotype is relevant at TP1+ (when V1 starts developing)
  traitData$V1_Ratio_TP1plus <- traitData$V1_Ratio *
    as.numeric(targets$Timepoint %in% c("TP1", "TP2", "TP3", "TP4"))

  # V1 phenotype at mature stages (TP2+)
  traitData$V1_Ratio_TP2plus <- traitData$V1_Ratio *
    as.numeric(targets$Timepoint %in% c("TP2", "TP3", "TP4"))

  # V2 phenotype is relevant at TP2+ (when V2 starts developing)
  if (exists("v2_summary_flat")) {
    traitData$V2_Ratio_TP2plus <- traitData$V2_Ratio *
      as.numeric(targets$Timepoint %in% c("TP2", "TP3", "TP4"))

    traitData$V2_Ratio_TP3plus <- traitData$V2_Ratio *
      as.numeric(targets$Timepoint %in% c("TP3", "TP4"))
  }

  cat("Continuous phenotype traits added:\n")
  phenotype_cols <- grep("V[12]_", colnames(traitData), value = TRUE)
  cat("  ", paste(phenotype_cols, collapse = ", "), "\n\n")
}

# Ensure row names match
traitData <- traitData[rownames(MEs), ]

cat("Trait matrix created with", ncol(traitData), "traits:\n")
cat("  Categorical:", paste(colnames(traitData)[1:17], collapse = ", "), "\n")
if (phenotype_available) {
  cat("  Phenotype:", paste(grep("V[12]_", colnames(traitData), value = TRUE), collapse = ", "), "\n")
}
cat("\n")

# ===== CALCULATE CORRELATIONS =====

cat("========================================\n")
cat("SECTION 4: MODULE-TRAIT CORRELATIONS\n")
cat("========================================\n\n")

# Calculate correlations
nSamples <- nrow(MEs)
moduleTraitCor <- cor(MEs, traitData, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

cat("Calculated correlations for", nrow(moduleTraitCor), "modules x",
    ncol(moduleTraitCor), "traits\n\n")

# ===== MAIN HEATMAP =====

cat("========================================\n")
cat("SECTION 5: TRAIT CORRELATION HEATMAPS\n")
cat("========================================\n\n")

# Create text matrix with correlation and p-value
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Main effects only heatmap
png("03_results/figures/20_WGCNA_traits/module_trait_main_effects.png",
    width = 800, height = 1000, res = 120)

par(mar = c(6, 10, 3, 2))

# Subset to main effects
main_effects <- c("Narrow", "Broad", "TP0", "TP1", "TP2", "TP3", "TP4")
cor_main <- moduleTraitCor[, main_effects]
pval_main <- moduleTraitPvalue[, main_effects]

# Text matrix for main effects
textMatrix_main <- paste(signif(cor_main, 2), "\n(",
                         signif(pval_main, 1), ")", sep = "")
dim(textMatrix_main) <- dim(cor_main)

labeledHeatmap(
  Matrix = cor_main,
  xLabels = colnames(cor_main),
  yLabels = rownames(cor_main),
  ySymbols = rownames(cor_main),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_main,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = "Module-Trait Correlations\n(Main Effects)"
)

dev.off()
cat("Saved: module_trait_main_effects.png\n")

# Full heatmap with interactions
png("03_results/figures/20_WGCNA_traits/module_trait_full.png",
    width = 1200, height = 1000, res = 120)

par(mar = c(8, 10, 3, 2))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(moduleTraitCor),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Module-Trait Correlations\n(All Effects and Interactions)"
)

dev.off()
cat("Saved: module_trait_full.png\n")

# ===== NARROW vs BROAD COMPARISON =====

cat("\n========================================\n")
cat("SECTION 6: NARROW vs BROAD MODULE DIFFERENCES\n")
cat("========================================\n\n")

# Calculate difference in correlation (Narrow - Broad) at each timepoint
diff_matrix <- matrix(0, nrow = nrow(moduleTraitCor), ncol = 5)
rownames(diff_matrix) <- rownames(moduleTraitCor)
colnames(diff_matrix) <- paste0("TP", 0:4)

for (i in 0:4) {
  narrow_col <- paste0("Narrow_TP", i)
  broad_col <- paste0("Broad_TP", i)
  diff_matrix[, paste0("TP", i)] <- moduleTraitCor[, narrow_col] - moduleTraitCor[, broad_col]
}

# Plot difference heatmap
png("03_results/figures/20_WGCNA_traits/module_trait_leaf_differences.png",
    width = 600, height = 1000, res = 120)

par(mar = c(6, 10, 3, 2))

labeledHeatmap(
  Matrix = diff_matrix,
  xLabels = colnames(diff_matrix),
  yLabels = rownames(diff_matrix),
  ySymbols = rownames(diff_matrix),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-Trait Correlation Differences\n(Narrow - Broad at Each Timepoint)"
)

dev.off()
cat("Saved: module_trait_leaf_differences.png\n")

# ===== IDENTIFY SIGNIFICANT MODULES =====

cat("\n========================================\n")
cat("SECTION 7: SIGNIFICANT MODULES\n")
cat("========================================\n\n")

# Find modules significantly associated with Narrow leaf type
sig_threshold <- 0.05
bonferroni_threshold <- 0.05 / nrow(moduleTraitCor)

# For Narrow trait
narrow_cor <- moduleTraitCor[, "Narrow"]
narrow_pval <- moduleTraitPvalue[, "Narrow"]

sig_narrow <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = narrow_cor,
  P_value = narrow_pval,
  FDR = p.adjust(narrow_pval, method = "BH"),
  Significant = narrow_pval < sig_threshold,
  Significant_Bonf = narrow_pval < bonferroni_threshold
)
sig_narrow <- sig_narrow[order(sig_narrow$P_value), ]

cat("Modules significantly correlated with Narrow leaf type (p < 0.05):\n")
print(sig_narrow[sig_narrow$Significant, c("Module", "Correlation", "P_value")])

# For TP0 (meristem stage)
tp0_cor <- moduleTraitCor[, "TP0"]
tp0_pval <- moduleTraitPvalue[, "TP0"]

sig_tp0 <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = tp0_cor,
  P_value = tp0_pval,
  FDR = p.adjust(tp0_pval, method = "BH"),
  Significant = tp0_pval < sig_threshold
)
sig_tp0 <- sig_tp0[order(sig_tp0$P_value), ]

cat("\nModules significantly correlated with TP0 (p < 0.05):\n")
print(sig_tp0[sig_tp0$Significant, c("Module", "Correlation", "P_value")])

# For Narrow_TP0 (JAG1 most active)
narrow_tp0_cor <- moduleTraitCor[, "Narrow_TP0"]
narrow_tp0_pval <- moduleTraitPvalue[, "Narrow_TP0"]

sig_narrow_tp0 <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = narrow_tp0_cor,
  P_value = narrow_tp0_pval,
  FDR = p.adjust(narrow_tp0_pval, method = "BH"),
  Significant = narrow_tp0_pval < sig_threshold
)
sig_narrow_tp0 <- sig_narrow_tp0[order(sig_narrow_tp0$P_value), ]

cat("\nModules significantly correlated with Narrow_TP0 (p < 0.05):\n")
print(sig_narrow_tp0[sig_narrow_tp0$Significant, c("Module", "Correlation", "P_value")])

# ===== CREATE SUMMARY TABLE =====

cat("\n========================================\n")
cat("SECTION 8: SUMMARY TABLES\n")
cat("========================================\n\n")

# Combine all correlations and p-values
trait_summary <- data.frame(
  Module = rownames(moduleTraitCor),

  # Main effects
  Cor_Narrow = moduleTraitCor[, "Narrow"],
  P_Narrow = moduleTraitPvalue[, "Narrow"],
  Cor_Broad = moduleTraitCor[, "Broad"],
  P_Broad = moduleTraitPvalue[, "Broad"],
  Cor_TP0 = moduleTraitCor[, "TP0"],
  P_TP0 = moduleTraitPvalue[, "TP0"],

  # Line-specific L:W ratio (actual phenotype)
  Cor_Line_LW = moduleTraitCor[, "Line_LW_Ratio"],
  P_Line_LW = moduleTraitPvalue[, "Line_LW_Ratio"],

  # Key interaction
  Cor_Narrow_TP0 = moduleTraitCor[, "Narrow_TP0"],
  P_Narrow_TP0 = moduleTraitPvalue[, "Narrow_TP0"],
  Cor_Broad_TP0 = moduleTraitCor[, "Broad_TP0"],
  P_Broad_TP0 = moduleTraitPvalue[, "Broad_TP0"],

  # Difference at TP0
  Diff_TP0 = diff_matrix[, "TP0"],

  stringsAsFactors = FALSE
)

# Add module size
trait_summary$Module_Size <- module_info$Size[match(
  gsub("ME", "", trait_summary$Module),
  as.character(module_info$Module)
)]

# Add significance flags
trait_summary$Sig_Narrow <- trait_summary$P_Narrow < 0.05
trait_summary$Sig_TP0 <- trait_summary$P_TP0 < 0.05
trait_summary$Sig_Narrow_TP0 <- trait_summary$P_Narrow_TP0 < 0.05
trait_summary$Sig_Line_LW <- trait_summary$P_Line_LW < 0.05

# Sort by correlation with Narrow
trait_summary <- trait_summary[order(abs(trait_summary$Cor_Narrow), decreasing = TRUE), ]

write.csv(trait_summary, "03_results/tables/WGCNA/module_trait_correlations.csv",
          row.names = FALSE)
cat("Saved: module_trait_correlations.csv\n")

# ===== COMPARISON: BINARY vs LINE-SPECIFIC L:W RATIO =====
cat("\n========================================\n")
cat("SECTION 8B: BINARY vs LINE-SPECIFIC PHENOTYPE COMPARISON\n")
cat("========================================\n\n")

cat("Comparing module correlations:\n")
cat("  Binary (Narrow): 2 unique values (0 or 1)\n")
cat("  Line L:W Ratio: 4 unique values (1.55, 1.59, 2.16, 2.36)\n\n")

# Correlation between the two methods
cor_comparison <- cor(trait_summary$Cor_Narrow, trait_summary$Cor_Line_LW)
cat("Correlation between Binary and Line-specific: r =", round(cor_comparison, 4), "\n\n")

# Create comparison plot
png("03_results/figures/20_WGCNA_traits/binary_vs_line_lw_comparison.png",
    width = 1200, height = 500, res = 120)

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Plot 1: Binary vs Line-specific correlation scatter
plot(trait_summary$Cor_Narrow, trait_summary$Cor_Line_LW,
     pch = 16, col = ifelse(trait_summary$Sig_Narrow, "red", "gray60"),
     xlab = "Correlation with Binary (Narrow=1, Broad=0)",
     ylab = "Correlation with Line-specific L:W Ratio",
     main = "Module-Trait Correlations:\nBinary vs Actual Phenotype",
     cex = 1.2)
abline(0, 1, lty = 2, col = "blue", lwd = 2)
abline(h = 0, v = 0, lty = 3, col = "gray")

# Add module labels for significant ones
sig_idx <- which(trait_summary$Sig_Narrow | trait_summary$Sig_Line_LW)
if (length(sig_idx) > 0) {
  text(trait_summary$Cor_Narrow[sig_idx],
       trait_summary$Cor_Line_LW[sig_idx],
       labels = trait_summary$Module[sig_idx],
       pos = 3, cex = 0.7, col = "darkblue")
}

legend("bottomright",
       legend = c(paste("r =", round(cor_comparison, 3)), "1:1 line"),
       lty = c(NA, 2), col = c("black", "blue"), pch = c(NA, NA), bty = "n")

# Plot 2: Ranked comparison
plot_order <- order(trait_summary$Cor_Narrow, decreasing = TRUE)
barplot_data <- rbind(
  trait_summary$Cor_Narrow[plot_order],
  trait_summary$Cor_Line_LW[plot_order]
)
colnames(barplot_data) <- trait_summary$Module[plot_order]

barplot(barplot_data, beside = TRUE,
        col = c("steelblue", "coral"),
        las = 2, cex.names = 0.6,
        ylab = "Correlation",
        main = "Module Correlations Ranked by Binary")
legend("topright", legend = c("Binary", "Line L:W"),
       fill = c("steelblue", "coral"), cex = 0.8)
abline(h = 0, lty = 2)

# Plot 3: Show the difference
diff_cor <- trait_summary$Cor_Line_LW - trait_summary$Cor_Narrow
names(diff_cor) <- trait_summary$Module
diff_cor <- sort(diff_cor)

barplot(diff_cor, las = 2, cex.names = 0.6,
        col = ifelse(diff_cor > 0, "forestgreen", "firebrick"),
        ylab = "Difference (Line L:W - Binary)",
        main = "Correlation Difference:\nLine-specific vs Binary")
abline(h = 0, lty = 2)

dev.off()
cat("Saved: binary_vs_line_lw_comparison.png\n")

# Print top modules by each method
cat("\nTop 6 modules by BINARY (Narrow) correlation:\n")
top_binary <- head(trait_summary[order(abs(trait_summary$Cor_Narrow), decreasing = TRUE),
                                  c("Module", "Cor_Narrow", "P_Narrow")], 6)
print(top_binary)

cat("\nTop 6 modules by LINE-SPECIFIC L:W correlation:\n")
top_line <- head(trait_summary[order(abs(trait_summary$Cor_Line_LW), decreasing = TRUE),
                                c("Module", "Cor_Line_LW", "P_Line_LW")], 6)
print(top_line)

# ===== VISUALIZE TOP MODULES =====

cat("\n========================================\n")
cat("SECTION 9: TOP MODULE VISUALIZATIONS\n")
cat("========================================\n\n")

# Select top modules by various criteria
top_narrow <- head(trait_summary$Module[order(abs(trait_summary$Cor_Narrow), decreasing = TRUE)], 6)
top_tp0 <- head(trait_summary$Module[order(abs(trait_summary$Cor_TP0), decreasing = TRUE)], 6)
top_diff <- head(trait_summary$Module[order(abs(trait_summary$Diff_TP0), decreasing = TRUE)], 6)

top_modules <- unique(c(top_narrow, top_tp0, top_diff))
top_modules <- head(top_modules, 12)

cat("Top modules selected for visualization:\n")
cat("  ", paste(top_modules, collapse = ", "), "\n\n")

# Create detailed plots for top modules
png("03_results/figures/20_WGCNA_traits/top_modules_temporal.png",
    width = 1400, height = 1000, res = 120)

par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))

for (me_name in top_modules) {
  mod_num <- as.numeric(gsub("ME", "", me_name))
  mod_color <- labels2colors(mod_num)
  mod_size <- trait_summary$Module_Size[trait_summary$Module == me_name]

  # Prepare data
  plot_data <- data.frame(
    ME = MEs[, me_name],
    Timepoint = factor(targets$Timepoint, levels = paste0("TP", 0:4)),
    Leaf_type = targets$Leaf_type
  )

  # Calculate means
  means <- aggregate(ME ~ Timepoint + Leaf_type, data = plot_data, FUN = mean)

  broad_means <- means[means$Leaf_type == "Broad", ]
  narrow_means <- means[means$Leaf_type == "Narrow", ]

  # Order
  broad_means <- broad_means[order(broad_means$Timepoint), ]
  narrow_means <- narrow_means[order(narrow_means$Timepoint), ]

  ylim_range <- range(c(broad_means$ME, narrow_means$ME))
  ylim_range <- ylim_range + c(-0.1, 0.1) * diff(ylim_range)

  # Get correlation values
  cor_narrow <- round(trait_summary$Cor_Narrow[trait_summary$Module == me_name], 2)
  p_narrow <- round(trait_summary$P_Narrow[trait_summary$Module == me_name], 4)

  plot(1:5, broad_means$ME, type = "b", pch = 16, col = "#E69F00",
       ylim = ylim_range, xlab = "Timepoint", ylab = "Module Eigengene",
       main = paste0(me_name, " (", mod_color, ", n=", mod_size, ")\nr=", cor_narrow, ", p=", p_narrow),
       xaxt = "n", lwd = 2, cex.main = 0.9)
  axis(1, at = 1:5, labels = paste0("TP", 0:4))
  lines(1:5, narrow_means$ME, type = "b", pch = 17, col = "#56B4E9", lwd = 2)
  abline(h = 0, lty = 2, col = "gray")

  if (me_name == top_modules[1]) {
    legend("topright", legend = c("Broad", "Narrow"),
           col = c("#E69F00", "#56B4E9"), pch = c(16, 17), lwd = 2, cex = 0.7)
  }
}

dev.off()
cat("Saved: top_modules_temporal.png\n")

# ===== BOXPLOTS BY CONDITION =====

png("03_results/figures/20_WGCNA_traits/top_modules_boxplots.png",
    width = 1200, height = 1000, res = 120)

par(mfrow = c(3, 4), mar = c(5, 4, 3, 1))

for (me_name in top_modules) {
  mod_num <- as.numeric(gsub("ME", "", me_name))
  mod_color <- labels2colors(mod_num)

  plot_data <- data.frame(
    ME = MEs[, me_name],
    Condition = paste(targets$Leaf_type, targets$Timepoint, sep = "_")
  )

  # Order conditions
  cond_order <- c(paste("Broad", paste0("TP", 0:4), sep = "_"),
                  paste("Narrow", paste0("TP", 0:4), sep = "_"))
  plot_data$Condition <- factor(plot_data$Condition, levels = cond_order)

  # Colors
  box_cols <- rep(c("#E69F00", "#56B4E9"), each = 5)

  boxplot(ME ~ Condition, data = plot_data,
          col = box_cols, las = 2, cex.axis = 0.6,
          main = paste0(me_name, " (", mod_color, ")"),
          ylab = "Module Eigengene")
  abline(h = 0, lty = 2, col = "gray")
  abline(v = 5.5, lty = 2, col = "black")
}

dev.off()
cat("Saved: top_modules_boxplots.png\n")

# ===== CORRELATION PLOT =====

cat("\n========================================\n")
cat("SECTION 10: CORRELATION SUMMARY PLOT\n")
cat("========================================\n\n")

png("03_results/figures/20_WGCNA_traits/correlation_summary.png",
    width = 1000, height = 800, res = 120)

par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

# Plot 1: Narrow correlation vs TP0 correlation
plot(trait_summary$Cor_Narrow, trait_summary$Cor_TP0,
     pch = 16, col = ifelse(trait_summary$Sig_Narrow, "red", "gray"),
     xlab = "Correlation with Narrow",
     ylab = "Correlation with TP0",
     main = "Module-Trait Correlation Patterns")
abline(h = 0, v = 0, lty = 2, col = "gray")
text(trait_summary$Cor_Narrow[trait_summary$Sig_Narrow],
     trait_summary$Cor_TP0[trait_summary$Sig_Narrow],
     labels = trait_summary$Module[trait_summary$Sig_Narrow],
     pos = 3, cex = 0.6)
legend("topleft", legend = c("Significant (Narrow)", "Not significant"),
       col = c("red", "gray"), pch = 16)

# Plot 2: Narrow vs Narrow_TP0
plot(trait_summary$Cor_Narrow, trait_summary$Cor_Narrow_TP0,
     pch = 16, col = ifelse(trait_summary$Sig_Narrow_TP0, "blue", "gray"),
     xlab = "Correlation with Narrow",
     ylab = "Correlation with Narrow_TP0",
     main = "Narrow Effect: Overall vs TP0-Specific")
abline(h = 0, v = 0, lty = 2, col = "gray")
abline(0, 1, lty = 3, col = "red")
text(trait_summary$Cor_Narrow[trait_summary$Sig_Narrow_TP0],
     trait_summary$Cor_Narrow_TP0[trait_summary$Sig_Narrow_TP0],
     labels = trait_summary$Module[trait_summary$Sig_Narrow_TP0],
     pos = 3, cex = 0.6)
legend("topleft", legend = c("Significant (Narrow_TP0)", "Not significant"),
       col = c("blue", "gray"), pch = 16)

dev.off()
cat("Saved: correlation_summary.png\n")

# ===== SAVE CHECKPOINT =====

# ===== PHENOTYPE-SPECIFIC ANALYSIS =====

if (phenotype_available) {
  cat("\n========================================\n")
  cat("SECTION 11: PHENOTYPE TRAIT CORRELATIONS\n")
  cat("========================================\n\n")

  # Create heatmap for phenotype traits only
  pheno_traits <- grep("V[12]_", colnames(moduleTraitCor), value = TRUE)

  if (length(pheno_traits) > 0) {
    png("03_results/figures/20_WGCNA_traits/module_phenotype_correlations.png",
        width = 900, height = 1000, res = 120)

    par(mar = c(8, 10, 3, 2))

    pheno_cor <- moduleTraitCor[, pheno_traits, drop = FALSE]
    pheno_pval <- moduleTraitPvalue[, pheno_traits, drop = FALSE]

    textMatrix_pheno <- paste(signif(pheno_cor, 2), "\n(",
                               signif(pheno_pval, 1), ")", sep = "")
    dim(textMatrix_pheno) <- dim(pheno_cor)

    labeledHeatmap(
      Matrix = pheno_cor,
      xLabels = colnames(pheno_cor),
      yLabels = rownames(pheno_cor),
      ySymbols = rownames(pheno_cor),
      colorLabels = FALSE,
      colors = blueWhiteRed(50),
      textMatrix = textMatrix_pheno,
      setStdMargins = FALSE,
      cex.text = 0.5,
      zlim = c(-1, 1),
      main = "Module-Phenotype Correlations\n(V1/V2 Leaf Measurements)"
    )

    dev.off()
    cat("Saved: module_phenotype_correlations.png\n")

    # Identify modules most correlated with V1_Ratio
    if ("V1_Ratio" %in% colnames(moduleTraitCor)) {
      v1_ratio_cor <- moduleTraitCor[, "V1_Ratio"]
      v1_ratio_pval <- moduleTraitPvalue[, "V1_Ratio"]

      sig_v1_ratio <- data.frame(
        Module = rownames(moduleTraitCor),
        Correlation = v1_ratio_cor,
        P_value = v1_ratio_pval,
        FDR = p.adjust(v1_ratio_pval, method = "BH"),
        Significant = v1_ratio_pval < 0.05
      )
      sig_v1_ratio <- sig_v1_ratio[order(abs(sig_v1_ratio$Correlation), decreasing = TRUE), ]

      cat("\nTop modules correlated with V1 Leaf Ratio:\n")
      print(head(sig_v1_ratio[, c("Module", "Correlation", "P_value")], 10))
      cat("\n")

      # Save phenotype correlations
      write.csv(sig_v1_ratio, "03_results/tables/WGCNA/module_V1_ratio_correlations.csv",
                row.names = FALSE)
    }

    # Plot V1_Ratio vs Narrow correlation
    if ("V1_Ratio" %in% colnames(moduleTraitCor)) {
      png("03_results/figures/20_WGCNA_traits/narrow_vs_v1ratio_correlation.png",
          width = 800, height = 700, res = 120)

      plot(moduleTraitCor[, "Narrow"], moduleTraitCor[, "V1_Ratio"],
           pch = 16, col = ifelse(moduleTraitPvalue[, "V1_Ratio"] < 0.05, "red", "gray"),
           xlab = "Correlation with Narrow (categorical)",
           ylab = "Correlation with V1 Ratio (continuous)",
           main = "Module Correlations: Categorical vs Continuous Phenotype")
      abline(h = 0, v = 0, lty = 2, col = "gray")
      abline(0, 1, lty = 3, col = "blue")

      # Add R2 annotation
      cor_test <- cor.test(moduleTraitCor[, "Narrow"], moduleTraitCor[, "V1_Ratio"])
      legend("bottomright",
             legend = paste("r =", round(cor_test$estimate, 3),
                          "\np =", format(cor_test$p.value, digits = 2)),
             bty = "n")

      # Label significant modules
      sig_idx <- which(moduleTraitPvalue[, "V1_Ratio"] < 0.05)
      if (length(sig_idx) > 0) {
        text(moduleTraitCor[sig_idx, "Narrow"],
             moduleTraitCor[sig_idx, "V1_Ratio"],
             labels = rownames(moduleTraitCor)[sig_idx],
             pos = 3, cex = 0.6)
      }

      dev.off()
      cat("Saved: narrow_vs_v1ratio_correlation.png\n")
    }
  }
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 12: SAVE CHECKPOINT\n")
cat("========================================\n\n")

# Prepare objects to save
objects_to_save <- c(
  "traitData",
  "moduleTraitCor",
  "moduleTraitPvalue",
  "diff_matrix",
  "trait_summary",
  "sig_narrow",
  "sig_tp0",
  "sig_narrow_tp0",
  # Keep previous objects
  "MEs",
  "ME_data",
  "module_info",
  "datExpr",
  "net",
  "wgcna_prep",
  "recommended_power",
  "module_colors",
  "targets",
  "logCPM"
)

# Add phenotype objects if available
if (phenotype_available) {
  objects_to_save <- c(objects_to_save,
                       "phenotype_summary",
                       "pheno_data",
                       "phenotype_available")
  if (exists("sig_v1_ratio")) {
    objects_to_save <- c(objects_to_save, "sig_v1_ratio")
  }
}

save(
  list = objects_to_save,
  file = "03_results/checkpoints/20_WGCNA_traits.RData"
)

cat("Checkpoint saved: 20_WGCNA_traits.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 20 COMPLETE: MODULE-TRAIT CORRELATION ANALYSIS\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat("  - Modules significantly correlated with Narrow:",
    sum(trait_summary$Sig_Narrow), "\n")
cat("  - Modules significantly correlated with TP0:",
    sum(trait_summary$Sig_TP0), "\n")
cat("  - Modules significantly correlated with Narrow_TP0:",
    sum(trait_summary$Sig_Narrow_TP0), "\n")

if (phenotype_available) {
  cat("\nPHENOTYPE DATA INTEGRATED:\n")
  cat("  - V1 Ratio (Broad):", round(phenotype_summary$V1_Ratio_mean[broad_idx], 3), "\n")
  cat("  - V1 Ratio (Narrow):", round(phenotype_summary$V1_Ratio_mean[narrow_idx], 3), "\n")
  if (exists("sig_v1_ratio")) {
    n_sig_v1 <- sum(sig_v1_ratio$Significant)
    cat("  - Modules correlated with V1_Ratio (p < 0.05):", n_sig_v1, "\n")
  }
}

# Print top correlations
cat("\nTop 5 modules correlated with Narrow leaf type:\n")
top5 <- head(trait_summary[order(abs(trait_summary$Cor_Narrow), decreasing = TRUE), ], 5)
for (i in 1:5) {
  cat(sprintf("  %s: r = %.3f, p = %.4f\n",
              top5$Module[i], top5$Cor_Narrow[i], top5$P_Narrow[i]))
}

cat("\nOUTPUT FILES:\n")
cat("  - Checkpoint: 03_results/checkpoints/20_WGCNA_traits.RData\n")
cat("  - Tables: 03_results/tables/WGCNA/\n")
cat("    - module_trait_correlations.csv\n")
if (phenotype_available) {
  cat("    - phenotype_summary.csv\n")
  cat("    - module_V1_ratio_correlations.csv\n")
}
cat("  - Figures: 03_results/figures/20_WGCNA_traits/\n")
cat("    - module_trait_main_effects.png\n")
cat("    - module_trait_full.png\n")
cat("    - module_trait_leaf_differences.png\n")
cat("    - top_modules_temporal.png\n")
cat("    - top_modules_boxplots.png\n")
cat("    - correlation_summary.png\n")
if (phenotype_available) {
  cat("    - module_phenotype_correlations.png\n")
  cat("    - narrow_vs_v1ratio_correlation.png\n")
}
cat("\n")

cat("NEXT STEP: Run Script 21 for JAG1 module analysis and hub genes\n\n")
