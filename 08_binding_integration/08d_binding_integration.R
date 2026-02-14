# Script 25e: Published Binding Data Integration
# Updated: 2026-01-25
# Now uses ACTUAL PEAK COORDINATES and BINDING REGION information

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 25e: BINDING DATA INTEGRATION\n")
cat("  Using ACTUAL PEAK COORDINATES from ChIP-seq and DAP-seq\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/binding_integration", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/25e_binding", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c("dplyr", "tidyr", "ggplot2", "VennDiagram", "RColorBrewer", "grid")

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    cat("  Loaded:", pkg, "\n")
  } else {
    cat("  Missing:", pkg, "\n")
  }
}
cat("\n")

# ===== LOAD RNA-SEQ DATA =====

cat("========================================\n")
cat("SECTION 1: LOAD RNA-SEQ RESULTS\n")
cat("========================================\n\n")

# Load JAG1 targets from your analysis
if (file.exists("03_results/checkpoints/14_JAG1_targets.RData")) {
  load("03_results/checkpoints/14_JAG1_targets.RData")
  jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
  cat("Loaded JAG1 targets from RNA-seq:", nrow(jag1_targets), "genes\n")

  cat("\nRNA-seq target tiers:\n")
  print(table(jag1_targets$Confidence_Tier))

} else {
  stop("JAG1 targets not found. Run Script 14 first.")
}

# Get all DE genes (for background calculation)
all_de_genes <- target_table$GeneID
cat("\nTotal genes tested:", length(all_de_genes), "\n")

# ===== LOAD BINDING DATA WITH PEAK COORDINATES =====

cat("\n========================================\n")
cat("SECTION 2: LOAD BINDING DATA (PEAK-BASED)\n")
cat("========================================\n\n")

# --- Huang et al. 2021 ChIP-seq (with binding region info) ---
huang_summary_file <- "01_data/published/Huang_2021_ChIPseq_gene_summary.csv"
huang_peaks_file <- "01_data/published/peaks/Huang_2021_ChIPseq_peaks_full.csv"

if (file.exists(huang_summary_file)) {
  huang_summary <- read.csv(huang_summary_file, stringsAsFactors = FALSE)
  cat("Loaded Huang ChIP-seq gene summary:", nrow(huang_summary), "genes\n")

  cat("\nChIP-seq binding region distribution:\n")
  cat("  Promoter binding:", sum(huang_summary$Has_Promoter_Binding), "genes\n")
  cat("  Genic (gene body) binding:", sum(huang_summary$Has_Genic_Binding), "genes\n")
  cat("  Downstream binding:", sum(huang_summary$Has_Downstream_Binding), "genes\n")

  cat("\nPrimary binding region:\n")
  print(table(huang_summary$Primary_Binding_Region))

  has_huang <- TRUE
} else {
  cat("ERROR: Huang gene summary not found.\n")
  cat("Run 25a_extract_binding_data.R first.\n")
  has_huang <- FALSE
}

# --- Wang et al. 2024 DAP-seq ---
wang_summary_file <- "01_data/published/Wang_2024_DAPseq_gene_summary.csv"

if (file.exists(wang_summary_file)) {
  wang_summary <- read.csv(wang_summary_file, stringsAsFactors = FALSE)
  cat("\nLoaded Wang DAP-seq gene summary:", nrow(wang_summary), "genes\n")
  cat("  With DAP-seq peak:", sum(wang_summary$Has_DAPseq_Peak), "\n")
  cat("  With motif:", sum(wang_summary$Has_Motif), "\n")
  has_wang <- TRUE
} else {
  cat("ERROR: Wang gene summary not found.\n")
  cat("Run 25a_extract_binding_data.R first.\n")
  has_wang <- FALSE
}

if (!has_wang || !has_huang) {
  stop("Missing binding data files. Run 25a_extract_binding_data.R first.")
}

# ===== CREATE INTEGRATED BINDING EVIDENCE =====

cat("\n========================================\n")
cat("SECTION 3: CREATE INTEGRATED EVIDENCE TABLE\n")
cat("========================================\n\n")

# Start with JAG1 targets
integrated <- data.frame(
  GeneID = jag1_targets$GeneID,
  DE_Tier = jag1_targets$Confidence_Tier,
  stringsAsFactors = FALSE
)

# Add logFC if available
if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  integrated$Mean_logFC <- jag1_targets$Mean_logFC_Pairwise
}

# --- Merge Huang ChIP-seq data ---
integrated <- merge(integrated, huang_summary[, c("GeneID", "Has_Promoter_Binding",
                                                   "Has_Genic_Binding", "Has_Downstream_Binding",
                                                   "Primary_Binding_Region", "N_Peaks",
                                                   "Closest_Peak_Distance")],
                    by = "GeneID", all.x = TRUE)

# Rename columns for clarity
colnames(integrated)[colnames(integrated) == "Has_Promoter_Binding"] <- "ChIP_Promoter"
colnames(integrated)[colnames(integrated) == "Has_Genic_Binding"] <- "ChIP_Genic"
colnames(integrated)[colnames(integrated) == "Has_Downstream_Binding"] <- "ChIP_Downstream"
colnames(integrated)[colnames(integrated) == "Primary_Binding_Region"] <- "ChIP_Primary_Region"
colnames(integrated)[colnames(integrated) == "N_Peaks"] <- "ChIP_N_Peaks"
colnames(integrated)[colnames(integrated) == "Closest_Peak_Distance"] <- "ChIP_Closest_Distance"

# Fill NAs with FALSE/0
integrated$ChIP_Promoter[is.na(integrated$ChIP_Promoter)] <- FALSE
integrated$ChIP_Genic[is.na(integrated$ChIP_Genic)] <- FALSE
integrated$ChIP_Downstream[is.na(integrated$ChIP_Downstream)] <- FALSE
integrated$ChIP_N_Peaks[is.na(integrated$ChIP_N_Peaks)] <- 0
integrated$ChIP_Primary_Region[is.na(integrated$ChIP_Primary_Region)] <- "No_ChIP_binding"

# Add overall ChIP-seq flag
integrated$Has_ChIPseq <- integrated$ChIP_Promoter | integrated$ChIP_Genic | integrated$ChIP_Downstream

# --- Merge Wang DAP-seq data ---
integrated <- merge(integrated, wang_summary[, c("GeneID", "Has_DAPseq_Peak", "Has_Motif")],
                    by = "GeneID", all.x = TRUE)

colnames(integrated)[colnames(integrated) == "Has_DAPseq_Peak"] <- "DAP_Peak"
colnames(integrated)[colnames(integrated) == "Has_Motif"] <- "DAP_Motif"

integrated$DAP_Peak[is.na(integrated$DAP_Peak)] <- FALSE
integrated$DAP_Motif[is.na(integrated$DAP_Motif)] <- FALSE

# Add overall DAP-seq flag
integrated$Has_DAPseq <- integrated$DAP_Peak | integrated$DAP_Motif

cat("Integrated data created:", nrow(integrated), "genes\n")

# ===== BINDING EVIDENCE CLASSIFICATION =====

cat("\n========================================\n")
cat("SECTION 4: EVIDENCE-BASED CLASSIFICATION\n")
cat("========================================\n\n")

# New classification based on actual binding data
#
# TIER 1 (GOLD): Both ChIP-seq AND DAP-seq evidence
#   - Strongest evidence for direct binding
#   - In vivo (ChIP) and in vitro (DAP) concordance
#
# TIER 2a: ChIP-seq only (in vivo binding)
#   - May include indirect targets (co-factor mediated)
#   - Further subdivided by binding location
#
# TIER 2b: DAP-seq only (in vitro binding)
#   - Direct sequence binding, but may not occur in vivo
#
# TIER 3: DE only (no binding evidence)
#   - Could be downstream/indirect targets

integrated$Binding_Class <- "DE_Only"

# Tier 1: Both methods
tier1_idx <- integrated$Has_ChIPseq & integrated$Has_DAPseq
integrated$Binding_Class[tier1_idx] <- "ChIP_AND_DAP"

# Tier 2a: ChIP-seq only
tier2a_idx <- integrated$Has_ChIPseq & !integrated$Has_DAPseq
integrated$Binding_Class[tier2a_idx] <- "ChIP_Only"

# Tier 2b: DAP-seq only
tier2b_idx <- !integrated$Has_ChIPseq & integrated$Has_DAPseq
integrated$Binding_Class[tier2b_idx] <- "DAP_Only"

cat("Binding evidence classification:\n")
print(table(integrated$Binding_Class))

# ===== DETAILED ChIP-seq BINDING LOCATION ANALYSIS =====

cat("\n========================================\n")
cat("SECTION 5: ChIP-seq BINDING LOCATION ANALYSIS\n")
cat("========================================\n\n")

# For genes with ChIP-seq binding, analyze WHERE binding occurs
chipseq_genes <- integrated[integrated$Has_ChIPseq, ]

cat("ChIP-seq binding location for JAG1 targets:\n")
cat("  Total targets with ChIP-seq binding:", nrow(chipseq_genes), "\n")
cat("  - Promoter binding:", sum(chipseq_genes$ChIP_Promoter), "\n")
cat("  - Genic (gene body) binding:", sum(chipseq_genes$ChIP_Genic), "\n")
cat("  - Downstream binding:", sum(chipseq_genes$ChIP_Downstream), "\n")

# Key insight: genes with ONLY genic binding
only_genic_binding <- chipseq_genes[chipseq_genes$ChIP_Genic &
                                     !chipseq_genes$ChIP_Promoter, ]
cat("\n  >>> Targets with ONLY gene body binding (no promoter):", nrow(only_genic_binding), "\n")
cat("  >>> These would be MISSED by promoter-only motif searches!\n")

# Binding location breakdown for targets
cat("\nPrimary binding region for targets with ChIP-seq:\n")
print(table(chipseq_genes$ChIP_Primary_Region))

# ===== CREATE CONFIDENCE TIERS =====

cat("\n========================================\n")
cat("SECTION 6: FINAL CONFIDENCE TIERS\n")
cat("========================================\n\n")

# Combine DE tier and binding evidence for final scoring
#
# Priority scoring:
# - DE Gold = 3, Silver = 2, Bronze = 1
# - ChIP+DAP = 3, ChIP only = 2, DAP only = 1.5, None = 0
# - Promoter binding bonus = +1

de_scores <- c("Gold" = 3, "Silver" = 2, "Bronze" = 1)
integrated$DE_Score <- de_scores[integrated$DE_Tier]
integrated$DE_Score[is.na(integrated$DE_Score)] <- 0

binding_scores <- c("ChIP_AND_DAP" = 3, "ChIP_Only" = 2, "DAP_Only" = 1.5, "DE_Only" = 0)
integrated$Binding_Score <- binding_scores[integrated$Binding_Class]

# Bonus for promoter binding (more likely canonical regulation)
integrated$Promoter_Bonus <- ifelse(integrated$ChIP_Promoter, 1, 0)

# Total priority score
integrated$Priority_Score <- integrated$DE_Score + integrated$Binding_Score + integrated$Promoter_Bonus

# Create final confidence tier
integrated$Final_Tier <- "Tier4_Exploratory"

# Tier 1: High confidence - Gold/Silver DE + both binding methods
integrated$Final_Tier[integrated$DE_Tier %in% c("Gold", "Silver") &
                       integrated$Binding_Class == "ChIP_AND_DAP"] <- "Tier1_HighConfidence"

# Tier 2: Good evidence - Gold/Silver DE + one binding method OR Bronze + both
integrated$Final_Tier[integrated$DE_Tier %in% c("Gold", "Silver") &
                       integrated$Binding_Class %in% c("ChIP_Only", "DAP_Only")] <- "Tier2_GoodEvidence"
integrated$Final_Tier[integrated$DE_Tier == "Bronze" &
                       integrated$Binding_Class == "ChIP_AND_DAP"] <- "Tier2_GoodEvidence"

# Tier 3: Moderate evidence - Bronze + one binding method
integrated$Final_Tier[integrated$DE_Tier == "Bronze" &
                       integrated$Binding_Class %in% c("ChIP_Only", "DAP_Only")] <- "Tier3_ModerateEvidence"

cat("Final confidence tier distribution:\n")
print(table(integrated$Final_Tier))

# Sort by priority score
integrated <- integrated[order(-integrated$Priority_Score, -integrated$DE_Score), ]

# ===== IDENTIFY TOP CANDIDATES =====

cat("\n========================================\n")
cat("SECTION 7: TOP DIRECT TARGET CANDIDATES\n")
cat("========================================\n\n")

# Highest confidence candidates
top_candidates <- integrated[integrated$Final_Tier == "Tier1_HighConfidence", ]

cat("TIER 1 HIGH-CONFIDENCE DIRECT TARGETS:\n")
cat("(Gold/Silver DE + ChIP-seq + DAP-seq)\n")
cat("Total:", nrow(top_candidates), "genes\n\n")

if (nrow(top_candidates) > 0) {
  display_cols <- c("GeneID", "DE_Tier", "Binding_Class", "ChIP_Primary_Region",
                    "Priority_Score")
  if ("Mean_logFC" %in% colnames(top_candidates)) {
    display_cols <- c(display_cols, "Mean_logFC")
  }
  print(head(top_candidates[, display_cols], 20))
}

# Genes with promoter binding
promoter_targets <- integrated[integrated$ChIP_Promoter & integrated$Has_DAPseq, ]
cat("\n\nTargets with PROMOTER binding + DAP-seq:", nrow(promoter_targets), "\n")

# Genes with gene body binding only
genic_only_targets <- integrated[integrated$ChIP_Genic & !integrated$ChIP_Promoter, ]
cat("Targets with GENE BODY binding only:", nrow(genic_only_targets), "\n")

# ===== OVERLAP STATISTICS =====

cat("\n========================================\n")
cat("SECTION 8: OVERLAP STATISTICS\n")
cat("========================================\n\n")

rna_seq_targets <- integrated$GeneID
huang_genes <- huang_summary$GeneID
wang_genes <- wang_summary$GeneID[wang_summary$Has_DAPseq_Peak]

overlap_huang <- intersect(rna_seq_targets, huang_genes)
overlap_wang <- intersect(rna_seq_targets, wang_genes)
triple_overlap <- intersect(overlap_huang, overlap_wang)

cat("Overlap with published binding data:\n")
cat("  Your DE targets:", length(rna_seq_targets), "\n")
cat("  Huang ChIP-seq genes:", length(huang_genes), "\n")
cat("  Wang DAP-seq genes (with peaks):", length(wang_genes), "\n\n")

cat("  DE targets with ChIP-seq:", length(overlap_huang),
    sprintf("(%.1f%%)", 100 * length(overlap_huang) / length(rna_seq_targets)), "\n")
cat("  DE targets with DAP-seq:", length(overlap_wang),
    sprintf("(%.1f%%)", 100 * length(overlap_wang) / length(rna_seq_targets)), "\n")
cat("  DE targets with BOTH:", length(triple_overlap),
    sprintf("(%.1f%%)", 100 * length(triple_overlap) / length(rna_seq_targets)), "\n")

# ===== ENRICHMENT ANALYSIS =====

cat("\n========================================\n")
cat("SECTION 9: ENRICHMENT ANALYSIS\n")
cat("========================================\n\n")

# Are DE genes enriched in binding datasets?
total_genes <- length(all_de_genes)
n_binding <- length(intersect(all_de_genes, huang_genes))
n_de <- length(rna_seq_targets)
n_de_and_binding <- sum(integrated$Has_ChIPseq)

a <- n_de_and_binding
b <- n_de - a
c <- n_binding - a
d <- total_genes - a - b - c

contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                      dimnames = list(c("DE_target", "Not_target"),
                                     c("ChIP_binding", "No_ChIP_binding")))

cat("Contingency table (DE vs ChIP-seq binding):\n")
print(contingency)

fisher_result <- fisher.test(contingency)

cat("\nEnrichment of DE targets in ChIP-seq binding:\n")
cat("  Odds ratio:", round(fisher_result$estimate, 2), "\n")
cat("  P-value:", format(fisher_result$p.value, scientific = TRUE), "\n")
cat("  95% CI:", round(fisher_result$conf.int[1], 2), "-",
    round(fisher_result$conf.int[2], 2), "\n")

if (fisher_result$p.value < 0.05 && fisher_result$estimate > 1) {
  cat("  >>> DE genes ARE significantly enriched in ChIP-seq binding data\n")
} else {
  cat("  >>> No significant enrichment\n")
}

# ===== SAVE RESULTS =====

cat("\n========================================\n")
cat("SECTION 10: SAVE RESULTS\n")
cat("========================================\n\n")

# Save main integrated table
write.csv(integrated,
          "03_results/tables/binding_integration/integrated_binding_analysis.csv",
          row.names = FALSE)
cat("Saved: integrated_binding_analysis.csv\n")

# Save by tier
for (tier in unique(integrated$Final_Tier)) {
  tier_genes <- integrated[integrated$Final_Tier == tier, ]
  filename <- paste0("03_results/tables/binding_integration/", tier, ".csv")
  write.csv(tier_genes, filename, row.names = FALSE)
  cat(sprintf("Saved: %s (%d genes)\n", basename(filename), nrow(tier_genes)))
}

# Save top candidates
if (nrow(top_candidates) > 0) {
  write.csv(top_candidates,
            "03_results/tables/binding_integration/TOP_high_confidence_targets.csv",
            row.names = FALSE)
  cat("Saved: TOP_high_confidence_targets.csv\n")
}

# Save genes with only gene body binding (important for mechanism!)
write.csv(genic_only_targets,
          "03_results/tables/binding_integration/targets_genic_binding_only.csv",
          row.names = FALSE)
cat("Saved: targets_genic_binding_only.csv (", nrow(genic_only_targets), "genes)\n")

# Save gene lists for external tools
writeLines(integrated$GeneID[integrated$Final_Tier == "Tier1_HighConfidence"],
           "03_results/tables/binding_integration/genes_tier1_high_confidence.txt")
writeLines(integrated$GeneID[integrated$ChIP_Promoter],
           "03_results/tables/binding_integration/genes_with_promoter_binding.txt")
writeLines(genic_only_targets$GeneID,
           "03_results/tables/binding_integration/genes_genic_binding_only.txt")

# === BACKWARD COMPATIBILITY ===
# Create legacy format output for downstream scripts (34d, 34e) that expect old column names
legacy_output <- integrated[, c("GeneID", "DE_Tier")]
legacy_output$In_Wang_DAPseq <- integrated$Has_DAPseq
legacy_output$In_Huang_ChIPseq <- integrated$Has_ChIPseq
legacy_output$Integration_Tier <- integrated$Final_Tier
legacy_output$N_Binding_Evidence <- as.integer(integrated$Has_ChIPseq) + as.integer(integrated$Has_DAPseq)
legacy_output$Priority_Score <- integrated$Priority_Score
if ("Mean_logFC" %in% colnames(integrated)) {
  legacy_output$Mean_logFC <- integrated$Mean_logFC
}

write.csv(legacy_output,
          "03_results/tables/binding_integration/integrated_targets_all.csv",
          row.names = FALSE)
cat("Saved: integrated_targets_all.csv (legacy format for backward compatibility)\n")

# ===== VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 11: VISUALIZATIONS\n")
cat("========================================\n\n")

# --- Figure 1: Binding class distribution ---
png("03_results/figures/25e_binding/binding_class_distribution.png",
    width = 10, height = 6, units = "in", res = 150)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Panel A: Binding evidence class
binding_counts <- table(integrated$Binding_Class)
binding_colors <- c("ChIP_AND_DAP" = "#2E7D32", "ChIP_Only" = "#1976D2",
                    "DAP_Only" = "#F57C00", "DE_Only" = "#757575")
barplot(binding_counts[names(binding_colors)],
        col = binding_colors,
        main = "A. Binding Evidence Classification",
        ylab = "Number of DE Targets",
        las = 2, cex.names = 0.8)

# Panel B: Final confidence tiers
tier_counts <- table(integrated$Final_Tier)
tier_colors <- c("Tier1_HighConfidence" = "#FFD700",
                 "Tier2_GoodEvidence" = "#C0C0C0",
                 "Tier3_ModerateEvidence" = "#CD7F32",
                 "Tier4_Exploratory" = "#808080")
barplot(tier_counts[names(tier_colors)],
        col = tier_colors,
        main = "B. Final Confidence Tiers",
        ylab = "Number of DE Targets",
        las = 2, cex.names = 0.7)

dev.off()
cat("Saved: binding_class_distribution.png\n")

# --- Figure 2: ChIP-seq binding location ---
png("03_results/figures/25e_binding/chipseq_binding_location.png",
    width = 10, height = 6, units = "in", res = 150)

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Panel A: Binding regions (can have multiple)
region_data <- c(
  "Promoter" = sum(integrated$ChIP_Promoter),
  "Gene Body" = sum(integrated$ChIP_Genic),
  "Downstream" = sum(integrated$ChIP_Downstream),
  "No ChIP" = sum(!integrated$Has_ChIPseq)
)
barplot(region_data, col = c("#4CAF50", "#2196F3", "#FF9800", "#9E9E9E"),
        main = "A. ChIP-seq Binding Regions\n(among DE targets)",
        ylab = "Number of Genes", las = 2)

# Panel B: Primary binding region for genes with ChIP
primary_counts <- table(integrated$ChIP_Primary_Region)
primary_colors <- c("Genic_Only" = "#2196F3", "Promoter_Only" = "#4CAF50",
                    "Promoter_and_Genic" = "#9C27B0", "Multiple" = "#FF5722",
                    "No_ChIP_binding" = "#9E9E9E")
barplot(primary_counts, col = primary_colors[names(primary_counts)],
        main = "B. Primary Binding Location",
        ylab = "Number of Genes", las = 2, cex.names = 0.7)

dev.off()
cat("Saved: chipseq_binding_location.png\n")

# --- Figure 3: Venn diagram ---
# Use ggplot2-based Venn instead to avoid VennDiagram calculation issues
png("03_results/figures/25e_binding/three_way_venn.png",
    width = 10, height = 8, units = "in", res = 150)

# Create a simple bar chart showing overlaps instead of Venn
par(mar = c(5, 8, 4, 2))

overlap_data <- c(
  "DE targets (total)" = length(rna_seq_targets),
  "With ChIP-seq binding" = sum(integrated$Has_ChIPseq),
  "With DAP-seq binding" = sum(integrated$Has_DAPseq),
  "With BOTH ChIP + DAP" = sum(integrated$Binding_Class == "ChIP_AND_DAP"),
  "ChIP-seq only" = sum(integrated$Binding_Class == "ChIP_Only"),
  "DAP-seq only" = sum(integrated$Binding_Class == "DAP_Only"),
  "DE only (no binding)" = sum(integrated$Binding_Class == "DE_Only")
)

barplot(overlap_data, horiz = TRUE, las = 1,
        col = c("gray40", "#009E73", "#56B4E9", "#D55E00", "#009E73", "#56B4E9", "gray70"),
        main = "Overlap: DE Targets with Published Binding Data",
        xlab = "Number of Genes",
        xlim = c(0, max(overlap_data) * 1.2))

# Add counts at end of bars
text(overlap_data + max(overlap_data) * 0.02,
     seq(0.7, by = 1.2, length.out = length(overlap_data)),
     labels = overlap_data, pos = 4, cex = 0.9)

dev.off()
cat("Saved: three_way_venn.png\n")

# --- Figure 4: Priority score distribution ---
png("03_results/figures/25e_binding/priority_scores.png",
    width = 8, height = 6, units = "in", res = 150)

hist(integrated$Priority_Score,
     breaks = seq(-0.5, max(integrated$Priority_Score) + 0.5, by = 0.5),
     col = "steelblue", border = "white",
     main = "Distribution of Priority Scores",
     xlab = "Priority Score (DE tier + Binding evidence + Promoter bonus)",
     ylab = "Number of Genes")
abline(v = 5.5, col = "red", lty = 2, lwd = 2)
text(5.7, par("usr")[4] * 0.9, "High\nconfidence\nthreshold", col = "red", adj = 0, cex = 0.9)

dev.off()
cat("Saved: priority_scores.png\n")

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SECTION 12: SAVE CHECKPOINT\n")
cat("========================================\n\n")

binding_integration <- list(
  integrated_results = integrated,
  top_candidates = top_candidates,
  genic_only_targets = genic_only_targets,
  overlap_stats = list(
    n_de_targets = length(rna_seq_targets),
    n_with_chipseq = sum(integrated$Has_ChIPseq),
    n_with_dapseq = sum(integrated$Has_DAPseq),
    n_with_both = sum(integrated$Binding_Class == "ChIP_AND_DAP"),
    n_promoter_binding = sum(integrated$ChIP_Promoter),
    n_genic_only = nrow(genic_only_targets)
  ),
  enrichment_test = fisher_result
)

save(binding_integration, integrated, top_candidates, genic_only_targets,
     file = "03_results/checkpoints/25e_binding_integration.RData")
cat("Checkpoint saved: 25e_binding_integration.RData\n")

# ===== FINAL SUMMARY =====

cat("\n================================================================\n")
cat("  SCRIPT 25e COMPLETE: BINDING DATA INTEGRATION\n")
cat("================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

cat("KEY FINDINGS:\n")
cat("=============\n\n")

cat("1. BINDING EVIDENCE CLASSIFICATION:\n")
print(table(integrated$Binding_Class))

cat("\n2. ChIP-seq BINDING LOCATION (among DE targets):\n")
cat("   - With promoter binding:", sum(integrated$ChIP_Promoter), "\n")
cat("   - With gene body binding:", sum(integrated$ChIP_Genic), "\n")
cat("   - With ONLY gene body binding:", nrow(genic_only_targets), "\n")
cat("   >>> These", nrow(genic_only_targets), "genes are missed by promoter motif searches!\n")

cat("\n3. FINAL CONFIDENCE TIERS:\n")
print(table(integrated$Final_Tier))

cat("\n4. TOP HIGH-CONFIDENCE TARGETS:", nrow(top_candidates), "\n")
cat("   (Gold/Silver DE + both ChIP-seq and DAP-seq evidence)\n")

cat("\n5. ENRICHMENT:\n")
cat("   DE targets are", ifelse(fisher_result$estimate > 1, "ENRICHED", "not enriched"),
    "in ChIP-seq binding\n")
cat("   Odds ratio:", round(fisher_result$estimate, 2), "\n")
cat("   P-value:", format(fisher_result$p.value, scientific = TRUE), "\n")

cat("\nOUTPUT FILES:\n")
cat("  Tables: 03_results/tables/binding_integration/\n")
cat("  Figures: 03_results/figures/25e_binding/\n")
cat("  Checkpoint: 03_results/checkpoints/25e_binding_integration.RData\n\n")

cat("IMPORTANT NOTE:\n")
cat("  This analysis uses ACTUAL PEAK COORDINATES from ChIP-seq,\n")
cat("  not arbitrary promoter windows with motif scanning.\n")
cat("  This is more accurate than the previous 3kb promoter motif approach.\n\n")
