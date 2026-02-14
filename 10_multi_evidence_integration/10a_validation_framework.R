# Script 26: Target Validation Framework

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 26: TARGET VALIDATION FRAMEWORK\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/26_validation", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD DATA =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

# Load integrated functional data
load("03_results/checkpoints/24_functional_integrated.RData")
cat("Loaded functional integration data\n")

# Load JAG1 targets
load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n")

# Try to load WGCNA
wgcna_available <- FALSE
if (file.exists("03_results/checkpoints/21_WGCNA_JAG1.RData")) {
  load("03_results/checkpoints/21_WGCNA_JAG1.RData")
  wgcna_available <- TRUE
  cat("Loaded WGCNA JAG1 data\n")
}

# Try to load co-expression
coexp_available <- FALSE
if (file.exists("03_results/checkpoints/15_JAG1_network.RData")) {
  load("03_results/checkpoints/15_JAG1_network.RData")
  coexp_available <- TRUE
  cat("Loaded JAG1 co-expression data\n")
}

jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("\nJAG1 targets:", nrow(jag1_targets), "\n")

# ===== LITERATURE REFERENCE GENES =====

cat("\n========================================\n")
cat("SECTION 2: LITERATURE REFERENCE\n")
cat("========================================\n\n")

# Known JAG1/JAGGED-related genes from literature
# Based on Arabidopsis JAGGED and soybean studies

literature_genes <- data.frame(
  Gene_Family = c(
    # Cell cycle related (JAGGED regulates cell proliferation)
    "Cyclin", "Cyclin", "CDK",
    # Cell wall (leaf expansion)
    "Expansin", "XTH", "Pectin",
    # Transcription factors (JAG1 regulatory network)
    "TCP", "TCP", "CUC", "AS1",
    # Hormone signaling
    "Auxin_TF", "GA_receptor",
    # Boundary/meristem
    "NAC", "KNOX"
  ),
  Expected_Pattern = c(
    # Cell cycle - expected up in narrow (de-repressed)
    "Up_in_Narrow", "Up_in_Narrow", "Up_in_Narrow",
    # Cell wall
    "Up_in_Narrow", "Up_in_Narrow", "Variable",
    # TFs
    "Up_in_Narrow", "Up_in_Narrow", "Variable", "Variable",
    # Hormones
    "Up_in_Narrow", "Variable",
    # Boundary
    "Variable", "Variable"
  ),
  Description = c(
    "Cell cycle progression", "Cell cycle progression", "Cell division",
    "Cell wall loosening", "Cell wall modification", "Cell wall structure",
    "Leaf development TF", "Leaf development TF", "Boundary formation", "Leaf polarity",
    "Auxin signaling", "GA perception",
    "Boundary/stress TF", "Meristem maintenance"
  ),
  stringsAsFactors = FALSE
)

cat("Reference gene families from literature:\n")
print(literature_genes)

write.csv(literature_genes, "03_results/tables/validation/literature_reference.csv",
          row.names = FALSE)

# ===== VALIDATION SCORING =====

cat("\n========================================\n")
cat("SECTION 3: VALIDATION SCORING\n")
cat("========================================\n\n")

# Create validation score for each target
validation_df <- data.frame(
  GeneID = jag1_targets$GeneID,
  Confidence_Tier = jag1_targets$Confidence_Tier,
  stringsAsFactors = FALSE
)

# Score 1: DE Consistency (from tier)
validation_df$DE_Score <- ifelse(validation_df$Confidence_Tier == "Gold", 3,
                                  ifelse(validation_df$Confidence_Tier == "Silver", 2, 1))

# Score 2: logFC magnitude
if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  validation_df$logFC <- jag1_targets$Mean_logFC_Pairwise
  validation_df$logFC_Score <- pmin(abs(validation_df$logFC), 3)  # Cap at 3
} else {
  validation_df$logFC_Score <- 0
}

# Score 3: Co-expression with JAG1
if (coexp_available && exists("tp0_network")) {
  # Check if genes are in JAG1 co-expression network
  if ("GeneID" %in% colnames(tp0_network)) {
    validation_df$In_Coexp_Network <- validation_df$GeneID %in% tp0_network$GeneID
    validation_df$Coexp_Score <- as.numeric(validation_df$In_Coexp_Network)
  } else {
    validation_df$Coexp_Score <- 0
  }
} else {
  validation_df$Coexp_Score <- 0
}

# Score 4: WGCNA module
if (wgcna_available && exists("target_modules")) {
  validation_df$WGCNA_Module <- target_modules[match(validation_df$GeneID, names(target_modules))]
  validation_df$In_Module <- !is.na(validation_df$WGCNA_Module) & validation_df$WGCNA_Module != 0
  validation_df$Module_Score <- as.numeric(validation_df$In_Module)
} else {
  validation_df$Module_Score <- 0
}

# Score 5: Functional relevance
if ("Category" %in% colnames(master_table)) {
  validation_df$Category <- master_table$Category[match(validation_df$GeneID, master_table$GeneID)]
  # Higher score for development/hormone categories
  relevant_categories <- c("Transcription", "Development", "Hormone_signaling", "Cell_wall")
  validation_df$Function_Score <- ifelse(validation_df$Category %in% relevant_categories, 1.5,
                                          ifelse(validation_df$Category != "Uncategorized", 0.5, 0))
} else {
  validation_df$Function_Score <- 0
}

# Total validation score
validation_df$Validation_Score <- validation_df$DE_Score +
  validation_df$logFC_Score +
  validation_df$Coexp_Score +
  validation_df$Module_Score +
  validation_df$Function_Score

# Rank by validation score
validation_df <- validation_df[order(-validation_df$Validation_Score), ]

cat("Validation score summary:\n")
cat("  Score range:", round(min(validation_df$Validation_Score), 2), "-",
    round(max(validation_df$Validation_Score), 2), "\n")
cat("  Mean score:", round(mean(validation_df$Validation_Score), 2), "\n\n")

# ===== VALIDATION CATEGORIES =====

cat("========================================\n")
cat("SECTION 4: VALIDATION CATEGORIES\n")
cat("========================================\n\n")

# Categorize by validation confidence
score_q75 <- quantile(validation_df$Validation_Score, 0.75)
score_q50 <- quantile(validation_df$Validation_Score, 0.50)

validation_df$Validation_Level <- ifelse(validation_df$Validation_Score >= score_q75, "High",
                                          ifelse(validation_df$Validation_Score >= score_q50, "Medium", "Low"))

cat("Validation levels:\n")
print(table(validation_df$Validation_Level))

# Cross-tabulate with confidence tier
tier_validation <- table(validation_df$Confidence_Tier, validation_df$Validation_Level)
cat("\nTier vs Validation Level:\n")
print(tier_validation)

# ===== TOP VALIDATED TARGETS =====

cat("\n========================================\n")
cat("SECTION 5: TOP VALIDATED TARGETS\n")
cat("========================================\n\n")

# Get top 50 validated targets
top_validated <- head(validation_df, 50)

cat("Top 10 validated targets:\n")
print_cols <- intersect(c("GeneID", "Confidence_Tier", "Validation_Score", "Category"),
                        colnames(validation_df))
print(head(validation_df[, print_cols], 10))

write.csv(top_validated, "03_results/tables/validation/top50_validated.csv",
          row.names = FALSE)
cat("\nSaved: top50_validated.csv\n")

# ===== VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 6: VISUALIZATIONS\n")
cat("========================================\n\n")

# Plot 1: Validation score distribution by tier
png("03_results/figures/26_validation/validation_by_tier.png",
    width = 800, height = 600, res = 120)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Boxplot
boxplot(Validation_Score ~ Confidence_Tier, data = validation_df,
        col = c("gray", "#FFD700", "#C0C0C0"),
        main = "Validation Score by Tier",
        xlab = "Confidence Tier",
        ylab = "Validation Score")

# Histogram
hist(validation_df$Validation_Score, breaks = 20,
     main = "Validation Score Distribution",
     xlab = "Validation Score",
     col = "steelblue", border = "white")
abline(v = c(score_q50, score_q75), col = c("orange", "red"), lty = 2, lwd = 2)
legend("topright", legend = c("Median", "75th %ile"),
       col = c("orange", "red"), lty = 2)

dev.off()
cat("Saved: validation_by_tier.png\n")

# Plot 2: Evidence type contribution
png("03_results/figures/26_validation/evidence_contribution.png",
    width = 800, height = 600, res = 120)

evidence_means <- c(
  DE = mean(validation_df$DE_Score),
  logFC = mean(validation_df$logFC_Score),
  Coexp = mean(validation_df$Coexp_Score),
  Module = mean(validation_df$Module_Score),
  Function = mean(validation_df$Function_Score)
)

barplot(evidence_means,
        main = "Mean Evidence Contribution to Validation Score",
        ylab = "Mean Score",
        col = c("darkblue", "blue", "lightblue", "green", "orange"),
        border = "white")

dev.off()
cat("Saved: evidence_contribution.png\n")

# Plot 3: High-confidence target summary
png("03_results/figures/26_validation/high_confidence_summary.png",
    width = 1000, height = 600, res = 120)

high_conf <- validation_df[validation_df$Validation_Level == "High", ]

par(mfrow = c(1, 2))

# Tier distribution in high confidence
tier_counts <- table(high_conf$Confidence_Tier)
pie(tier_counts,
    main = paste("High Validation Targets (n=", nrow(high_conf), ")\nby Confidence Tier"),
    col = c("#CD7F32", "#FFD700", "#C0C0C0"))

# Category distribution
if ("Category" %in% colnames(high_conf)) {
  cat_counts <- sort(table(high_conf$Category), decreasing = TRUE)
  cat_counts <- head(cat_counts, 8)

  par(mar = c(8, 4, 4, 2))
  barplot(cat_counts, las = 2,
          main = "High Validation Targets\nby Functional Category",
          col = "darkgreen",
          border = "white")
}

dev.off()
cat("Saved: high_confidence_summary.png\n")

# ===== SAVE RESULTS =====

cat("\n========================================\n")
cat("SECTION 7: SAVE RESULTS\n")
cat("========================================\n\n")

# Save full validation table
write.csv(validation_df, "03_results/tables/validation/all_targets_validated.csv",
          row.names = FALSE)
cat("Saved: all_targets_validated.csv\n")

# Save high-confidence list
write.csv(high_conf, "03_results/tables/validation/high_confidence_targets.csv",
          row.names = FALSE)
cat("Saved: high_confidence_targets.csv\n")

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 8: SAVE CHECKPOINT\n")
cat("========================================\n\n")

validation_results <- list(
  validation_df = validation_df,
  high_confidence = high_conf,
  score_thresholds = c(median = score_q50, q75 = score_q75),
  literature_genes = literature_genes
)

save(
  validation_results,
  validation_df,
  file = "03_results/checkpoints/26_validation.RData"
)

cat("Checkpoint saved: 26_validation.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 26 COMPLETE: TARGET VALIDATION FRAMEWORK\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("KEY RESULTS:\n")
cat("  - Validated", nrow(validation_df), "JAG1 targets\n")
cat("  - High-confidence targets:", nrow(high_conf), "\n")
cat("  - Score range:", round(min(validation_df$Validation_Score), 2), "-",
    round(max(validation_df$Validation_Score), 2), "\n")
cat("  - Top validated gene:", validation_df$GeneID[1], "\n\n")

cat("NEXT STEP: Run Script 27 for cross-validation summary\n\n")
