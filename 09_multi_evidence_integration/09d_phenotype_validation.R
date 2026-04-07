# Script 32: Phenotype-Enhanced Validation

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 32: PHENOTYPE-ENHANCED VALIDATION\n")
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

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")
library(dplyr)
library(tidyr)
cat("  Packages loaded\n\n")

# ===== LOAD DATA =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

# Load JAG1 targets
if (file.exists("03_results/checkpoints/14_JAG1_targets.RData")) {
  load("03_results/checkpoints/14_JAG1_targets.RData")
  jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
  cat("Loaded JAG1 targets:", nrow(jag1_targets), "\n")
} else {
  stop("JAG1 targets not found. Run Script 14 first.")
}

# Load previous validation (if available)
if (file.exists("03_results/checkpoints/27_cross_validation.RData")) {
  load("03_results/checkpoints/27_cross_validation.RData")
  cat("Loaded cross-validation results\n")
  has_validation <- TRUE
} else {
  cat("Cross-validation not found. Will create new validation.\n")
  has_validation <- FALSE
}

# Load phenotype integration
if (file.exists("03_results/checkpoints/31_phenotype_integration.RData")) {
  load("03_results/checkpoints/31_phenotype_integration.RData")
  cat("Loaded phenotype integration results\n")
  has_phenotype <- TRUE

  # Extract gene-phenotype correlations
  if (exists("phenotype_results") && !is.null(phenotype_results$gene_pheno_cor)) {
    gene_pheno_cor <- phenotype_results$gene_pheno_cor
    cat("  Gene-phenotype correlations:", nrow(gene_pheno_cor), "genes\n")
  } else if (exists("gene_pheno_cor")) {
    cat("  Gene-phenotype correlations:", nrow(gene_pheno_cor), "genes\n")
  } else {
    has_phenotype <- FALSE
    cat("  Warning: gene_pheno_cor not found in phenotype results\n")
  }
} else {
  cat("Phenotype integration not found. Run Script 31 first.\n")
  has_phenotype <- FALSE
}

cat("\n")

# ===== DEFINE VALIDATION WEIGHTS =====

cat("========================================\n")
cat("SECTION 2: DEFINE VALIDATION WEIGHTS\n")
cat("========================================\n\n")

# Updated weights including phenotype
# Total should sum to 1.0

weights <- list(
  DE_Tier = 0.30,        # Differential expression tier (was 0.40)
  Coexpression = 0.15,   # Co-expression with JAG1 (was 0.20)
  Module = 0.15,         # WGCNA module membership (was 0.20)
  Functional = 0.15,     # GO/pathway evidence (was 0.20)
  Phenotype = 0.25       # ** NEW: Phenotype correlation **
)

cat("Validation weights:\n")
for (w in names(weights)) {
  cat("  ", w, ":", weights[[w]] * 100, "%\n")
}
cat("  Total:", sum(unlist(weights)) * 100, "%\n\n")

# ===== CREATE ENHANCED VALIDATION TABLE =====

cat("========================================\n")
cat("SECTION 3: CREATE VALIDATION TABLE\n")
cat("========================================\n\n")

# Start with JAG1 targets
validation_enhanced <- data.frame(
  GeneID = jag1_targets$GeneID,
  Confidence_Tier = jag1_targets$Confidence_Tier,
  stringsAsFactors = FALSE
)

# Add DE evidence score (from tier)
tier_scores <- c("Gold" = 1.0, "Silver" = 0.67, "Bronze" = 0.33)
validation_enhanced$DE_Score <- tier_scores[validation_enhanced$Confidence_Tier]

cat("DE scores assigned based on tier\n")

# Add mean logFC if available
if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  validation_enhanced$Mean_logFC <- jag1_targets$Mean_logFC_Pairwise
}

# ===== ADD PHENOTYPE EVIDENCE =====

cat("\n========================================\n")
cat("SECTION 4: ADD PHENOTYPE EVIDENCE\n")
cat("========================================\n\n")

if (has_phenotype && exists("gene_pheno_cor")) {

  # Merge phenotype correlations
  validation_enhanced <- merge(
    validation_enhanced,
    gene_pheno_cor[, c("GeneID", "Cor_LW_Ratio", "FDR_LW_Ratio", "Pheno_Association")],
    by = "GeneID",
    all.x = TRUE
  )

  # Calculate phenotype score
  # Higher positive correlation with L:W ratio = better (narrow leaf trait)
  # Because JAG1 targets should be UP in narrow (high L:W) leaves

  validation_enhanced$Pheno_Score <- 0

  # Significant positive correlation with L:W ratio
  sig_pos <- !is.na(validation_enhanced$FDR_LW_Ratio) &
             validation_enhanced$FDR_LW_Ratio < 0.05 &
             validation_enhanced$Cor_LW_Ratio > 0

  # Score based on correlation strength
  validation_enhanced$Pheno_Score[sig_pos] <-
    pmin(validation_enhanced$Cor_LW_Ratio[sig_pos] / 0.5, 1)  # Max score at r=0.5

  # Bonus for strong correlation (r > 0.4)
  strong_cor <- !is.na(validation_enhanced$Cor_LW_Ratio) &
                validation_enhanced$Cor_LW_Ratio > 0.4

  validation_enhanced$Pheno_Score[strong_cor] <- 1.0

  cat("Phenotype scores calculated:\n")
  cat("  Genes with phenotype data:", sum(!is.na(validation_enhanced$Cor_LW_Ratio)), "\n")
  cat("  Significant positive correlation:", sum(sig_pos), "\n")
  cat("  Strong correlation (r > 0.4):", sum(strong_cor, na.rm = TRUE), "\n")

} else {
  cat("No phenotype data available. Phenotype score set to NA.\n")
  validation_enhanced$Cor_LW_Ratio <- NA
  validation_enhanced$FDR_LW_Ratio <- NA
  validation_enhanced$Pheno_Association <- NA
  validation_enhanced$Pheno_Score <- NA
}

# ===== ADD OTHER EVIDENCE SCORES =====

cat("\n========================================\n")
cat("SECTION 5: ADD OTHER EVIDENCE\n")
cat("========================================\n\n")

# Placeholder scores if previous validation not available
# These would be filled from previous scripts

if (!("Coexpr_Score" %in% colnames(validation_enhanced))) {
  validation_enhanced$Coexpr_Score <- 0.5  # Default moderate
  cat("Co-expression score: using default (0.5)\n")
}

if (!("Module_Score" %in% colnames(validation_enhanced))) {
  validation_enhanced$Module_Score <- 0.5  # Default moderate
  cat("Module score: using default (0.5)\n")
}

if (!("Functional_Score" %in% colnames(validation_enhanced))) {
  validation_enhanced$Functional_Score <- 0.5  # Default moderate
  cat("Functional score: using default (0.5)\n")
}

# ===== CALCULATE COMPOSITE SCORE =====

cat("\n========================================\n")
cat("SECTION 6: CALCULATE COMPOSITE SCORE\n")
cat("========================================\n\n")

# Calculate weighted composite score
validation_enhanced$Composite_Score <- 0

# Add each component
validation_enhanced$Composite_Score <- validation_enhanced$Composite_Score +
  weights$DE_Tier * validation_enhanced$DE_Score

validation_enhanced$Composite_Score <- validation_enhanced$Composite_Score +
  weights$Coexpression * validation_enhanced$Coexpr_Score

validation_enhanced$Composite_Score <- validation_enhanced$Composite_Score +
  weights$Module * validation_enhanced$Module_Score

validation_enhanced$Composite_Score <- validation_enhanced$Composite_Score +
  weights$Functional * validation_enhanced$Functional_Score

# Add phenotype if available
if (has_phenotype) {
  # For genes with phenotype data
  has_pheno <- !is.na(validation_enhanced$Pheno_Score)
  validation_enhanced$Composite_Score[has_pheno] <-
    validation_enhanced$Composite_Score[has_pheno] +
    weights$Phenotype * validation_enhanced$Pheno_Score[has_pheno]

  # For genes without phenotype data, redistribute weight
  no_pheno <- is.na(validation_enhanced$Pheno_Score)
  if (sum(no_pheno) > 0) {
    # Normalize remaining scores
    remaining_weight <- 1 - weights$Phenotype
    validation_enhanced$Composite_Score[no_pheno] <-
      validation_enhanced$Composite_Score[no_pheno] / remaining_weight
  }
}

# Ensure score is 0-1
validation_enhanced$Composite_Score <- pmin(pmax(validation_enhanced$Composite_Score, 0), 1)

cat("Composite scores calculated\n")
cat("  Score range:", round(min(validation_enhanced$Composite_Score), 3), "-",
    round(max(validation_enhanced$Composite_Score), 3), "\n")
cat("  Mean score:", round(mean(validation_enhanced$Composite_Score), 3), "\n")

# ===== CLASSIFY VALIDATION LEVELS =====

cat("\n========================================\n")
cat("SECTION 7: CLASSIFY VALIDATION LEVELS\n")
cat("========================================\n\n")

# Create validation tiers based on composite score
validation_enhanced$Validation_Level <- cut(
  validation_enhanced$Composite_Score,
  breaks = c(0, 0.4, 0.6, 0.8, 1.0),
  labels = c("Low", "Moderate", "High", "Very_High"),
  include.lowest = TRUE
)

# Create phenotype-enhanced tier
validation_enhanced$Pheno_Enhanced_Tier <- validation_enhanced$Confidence_Tier

# Upgrade tier if strong phenotype support
if (has_phenotype) {
  # Genes with strong phenotype correlation get upgraded
  strong_pheno <- !is.na(validation_enhanced$Pheno_Score) &
                  validation_enhanced$Pheno_Score >= 0.8

  # Silver -> Gold if strong phenotype
  upgrade_silver <- strong_pheno & validation_enhanced$Confidence_Tier == "Silver"
  validation_enhanced$Pheno_Enhanced_Tier[upgrade_silver] <- "Gold*"

  # Bronze -> Silver if strong phenotype
  upgrade_bronze <- strong_pheno & validation_enhanced$Confidence_Tier == "Bronze"
  validation_enhanced$Pheno_Enhanced_Tier[upgrade_bronze] <- "Silver*"

  cat("Tier upgrades based on phenotype:\n")
  cat("  Silver -> Gold*:", sum(upgrade_silver), "\n")
  cat("  Bronze -> Silver*:", sum(upgrade_bronze), "\n")
}

# Summary by validation level
cat("\nValidation level distribution:\n")
print(table(validation_enhanced$Validation_Level))

cat("\nEnhanced tier distribution:\n")
print(table(validation_enhanced$Pheno_Enhanced_Tier))

# ===== IDENTIFY TOP CANDIDATES =====

cat("\n========================================\n")
cat("SECTION 8: TOP CANDIDATES\n")
cat("========================================\n\n")

# Sort by composite score
validation_enhanced <- validation_enhanced[order(-validation_enhanced$Composite_Score), ]

# Top 20 candidates
top_20 <- head(validation_enhanced, 20)

cat("TOP 20 JAG1 TARGET CANDIDATES:\n")
cat("(Ranked by composite validation score)\n\n")

print(top_20[, c("GeneID", "Confidence_Tier", "Pheno_Enhanced_Tier",
                  "Composite_Score", "Cor_LW_Ratio")])

# ===== PHENOTYPE-PRIORITIZED LIST =====

cat("\n========================================\n")
cat("SECTION 9: PHENOTYPE-PRIORITIZED TARGETS\n")
cat("========================================\n\n")

if (has_phenotype) {
  # Genes with BOTH high DE tier AND phenotype correlation
  pheno_prioritized <- validation_enhanced %>%
    filter(
      Confidence_Tier %in% c("Gold", "Silver"),
      !is.na(Cor_LW_Ratio),
      Cor_LW_Ratio > 0.3,
      FDR_LW_Ratio < 0.05
    ) %>%
    arrange(desc(Composite_Score))

  cat("Phenotype-prioritized targets:\n")
  cat("  (Gold/Silver tier + significant positive L:W correlation)\n")
  cat("  Total:", nrow(pheno_prioritized), "genes\n\n")

  if (nrow(pheno_prioritized) > 0) {
    print(head(pheno_prioritized[, c("GeneID", "Confidence_Tier",
                                      "Mean_logFC", "Cor_LW_Ratio", "Composite_Score")], 15))

    write.csv(pheno_prioritized,
              "03_results/tables/validation/phenotype_prioritized_targets.csv",
              row.names = FALSE)
  }
}

# ===== SAVE RESULTS =====

cat("\n========================================\n")
cat("SECTION 10: SAVE RESULTS\n")
cat("========================================\n\n")

# Save full validation table
write.csv(validation_enhanced,
          "03_results/tables/validation/phenotype_enhanced_validation.csv",
          row.names = FALSE)
cat("Saved: phenotype_enhanced_validation.csv\n")

# Save high-confidence list
high_confidence <- validation_enhanced %>%
  filter(Validation_Level %in% c("High", "Very_High"))

write.csv(high_confidence,
          "03_results/tables/validation/high_confidence_targets.csv",
          row.names = FALSE)
cat("Saved: high_confidence_targets.csv (", nrow(high_confidence), " genes)\n")

# ===== VISUALIZATION =====

cat("\n========================================\n")
cat("SECTION 11: VISUALIZATION\n")
cat("========================================\n\n")

dir.create("03_results/figures/32_validation", recursive = TRUE, showWarnings = FALSE)

# Plot: Composite score distribution by tier
png("03_results/figures/32_validation/composite_score_by_tier.png",
    width = 10, height = 6, units = "in", res = 150)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Boxplot by tier
boxplot(Composite_Score ~ Confidence_Tier, data = validation_enhanced,
        col = c("#CD7F32", "#FFD700", "#C0C0C0"),
        main = "Composite Score by DE Tier",
        ylab = "Composite Validation Score",
        xlab = "Confidence Tier")

# Phenotype correlation vs DE score
if (has_phenotype) {
  plot(validation_enhanced$DE_Score, validation_enhanced$Cor_LW_Ratio,
       pch = 16, col = rgb(0, 0, 1, 0.5),
       xlab = "DE Score (from tier)",
       ylab = "Correlation with L:W Ratio",
       main = "DE Evidence vs Phenotype Evidence")
  abline(h = 0, lty = 2, col = "gray")
  abline(h = 0.3, lty = 2, col = "red")

  # Highlight top candidates
  top_idx <- validation_enhanced$Composite_Score >= 0.8
  points(validation_enhanced$DE_Score[top_idx],
         validation_enhanced$Cor_LW_Ratio[top_idx],
         pch = 1, col = "red", cex = 1.5, lwd = 2)
  legend("bottomright", legend = c("All targets", "Top candidates"),
         pch = c(16, 1), col = c(rgb(0, 0, 1, 0.5), "red"), cex = 0.8)
}

dev.off()
cat("Saved: composite_score_by_tier.png\n")

# Plot: Validation evidence heatmap for top genes
if (nrow(top_20) >= 10) {
  png("03_results/figures/32_validation/top_targets_evidence.png",
      width = 8, height = 10, units = "in", res = 150)

  # Create evidence matrix
  evidence_cols <- c("DE_Score", "Coexpr_Score", "Module_Score",
                     "Functional_Score", "Pheno_Score")
  evidence_cols <- evidence_cols[evidence_cols %in% colnames(top_20)]

  evidence_matrix <- as.matrix(top_20[1:min(20, nrow(top_20)), evidence_cols])
  rownames(evidence_matrix) <- top_20$GeneID[1:min(20, nrow(top_20))]

  # Replace NA with 0 for visualization
  evidence_matrix[is.na(evidence_matrix)] <- 0

  # Simple heatmap
  par(mar = c(8, 10, 3, 2))
  image(t(evidence_matrix[nrow(evidence_matrix):1, ]),
        col = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
        axes = FALSE,
        main = "Evidence Scores for Top 20 Targets")

  axis(1, at = seq(0, 1, length.out = ncol(evidence_matrix)),
       labels = evidence_cols, las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, 1, length.out = nrow(evidence_matrix)),
       labels = rev(rownames(evidence_matrix)), las = 1, cex.axis = 0.7)

  dev.off()
  cat("Saved: top_targets_evidence.png\n")
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 12: SAVE CHECKPOINT\n")
cat("========================================\n\n")

save(
  validation_enhanced,
  weights,
  file = "03_results/checkpoints/32_phenotype_validation.RData"
)

cat("Checkpoint saved: 32_phenotype_validation.RData\n")

# ===== SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 32 COMPLETE: PHENOTYPE-ENHANCED VALIDATION\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("VALIDATION SUMMARY:\n")
cat("  Total JAG1 targets validated:", nrow(validation_enhanced), "\n")
cat("  Phenotype data available:", has_phenotype, "\n\n")

cat("VALIDATION WEIGHTS:\n")
for (w in names(weights)) {
  cat("  ", w, ":", weights[[w]] * 100, "%\n")
}

cat("\nVALIDATION LEVEL DISTRIBUTION:\n")
print(table(validation_enhanced$Validation_Level))

if (has_phenotype) {
  cat("\nPHENOTYPE EVIDENCE:\n")
  cat("  Targets with phenotype correlation:",
      sum(!is.na(validation_enhanced$Cor_LW_Ratio)), "\n")
  cat("  Significant positive correlation:",
      sum(validation_enhanced$Cor_LW_Ratio > 0.3 &
          validation_enhanced$FDR_LW_Ratio < 0.05, na.rm = TRUE), "\n")
}

cat("\nOUTPUT FILES:\n")
cat("  - 03_results/tables/validation/phenotype_enhanced_validation.csv\n")
cat("  - 03_results/tables/validation/high_confidence_targets.csv\n")
if (has_phenotype) {
  cat("  - 03_results/tables/validation/phenotype_prioritized_targets.csv\n")
}
cat("  - 03_results/figures/32_validation/*.png\n\n")

cat("INTERPRETATION:\n")
cat("  - Composite Score 0.8-1.0: Very High confidence targets\n")
cat("  - Composite Score 0.6-0.8: High confidence targets\n")
cat("  - Phenotype correlation based on V1 leaf L:W ratio\n")
cat("  - Positive correlation with L:W ratio = associated with narrow phenotype\n")
cat("  - Genes with high L:W ratio correlation are candidate JAG1 targets\n\n")
