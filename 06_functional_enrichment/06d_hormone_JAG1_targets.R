#!/usr/bin/env Rscript
# ==============================================================================
# Script 37b: Hormone Pathway Enrichment in JAG1 Targets
# ==============================================================================
# Purpose: Test whether JAG1 target genes are enriched for hormone pathway genes
#
# Question: Is JAG1 directly responsible for regulating hormone genes?
#
# Input Files:
#   - 03_results/checkpoints/37a_hormone_DE.RData: Hormone gene mappings from 37a
#   - 03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv: JAG1 targets
#   - 03_results/checkpoints/05_normalized.RData: Background genes
#
# Dependency: Run AFTER 37a (DE enrichment analysis)
#
# Output:
#   - Hormone pathway enrichment statistics for JAG1 targets
#   - Detailed gene lists with expression changes
#   - Visualization of results
#
# KEGG Reference:
#   - Pathway: gmx04075 (Plant hormone signal transduction - Glycine max)
# ==============================================================================

cat("=======================================================================\n")
cat("Script 37b: Hormone Pathway Enrichment in JAG1 Targets\n")
cat("=======================================================================\n\n")

# ------------------------------------------------------------------------------
# 0. Set Working Directory to Project Root
# ------------------------------------------------------------------------------
script_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) file.path(getwd(), "02_scripts"))
project_root <- dirname(script_dir)
setwd(project_root)
cat("Working directory:", getwd(), "\n\n")

# ------------------------------------------------------------------------------
# 1. Setup and Load Dependencies
# ------------------------------------------------------------------------------
cat("1. Loading packages and data from 37a...\n")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Create output directories
dir.create("03_results/tables/hormone_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("03_results/figures/37b_hormone_JAG1", showWarnings = FALSE, recursive = TRUE)

# Load checkpoint from 37a
if (!file.exists("03_results/checkpoints/37a_hormone_DE.RData")) {
  stop("Checkpoint from 37a not found. Please run 37a_hormone_DE_enrichment.R first.")
}

load("03_results/checkpoints/37a_hormone_DE.RData")
cat("   Loaded hormone mappings from 37a\n")
cat("   Hormone genes in background:", length(hormone_de_analysis$hormone_in_background), "\n")

# Extract variables for convenience
hormone_in_background <- hormone_de_analysis$hormone_in_background
background_genes <- hormone_de_analysis$background_genes
de_hormone_genes <- hormone_de_analysis$de_hormone_genes

# ------------------------------------------------------------------------------
# 2. Load JAG1 Targets
# ------------------------------------------------------------------------------
cat("\n2. Loading JAG1 targets...\n")

jag1_targets <- read_csv("03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv",
                         show_col_types = FALSE)
target_genes <- jag1_targets$GeneID
cat("   JAG1 targets:", length(target_genes), "genes\n")

# By tier
tier_counts <- table(jag1_targets$Confidence_Tier)
cat("   By tier:\n")
for (tier in names(tier_counts)) {
  cat("     ", tier, ":", tier_counts[tier], "\n")
}

# ------------------------------------------------------------------------------
# 3. Enrichment Analysis (Fisher's Exact Test)
# ------------------------------------------------------------------------------
cat("\n3. Performing enrichment analysis...\n")

n_targets <- length(target_genes)
n_background <- length(background_genes)

# JAG1 targets that are hormone pathway genes
targets_in_hormone <- intersect(target_genes, hormone_in_background)
n_targets_hormone <- length(targets_in_hormone)

# Build contingency table
a <- n_targets_hormone
b <- n_targets - a
c <- length(hormone_in_background) - a
d <- n_background - n_targets - c

cat("\n   Contingency table:\n")
cat("                          Hormone    Not-Hormone    Total\n")
cat(sprintf("   JAG1 target            %7d    %11d    %5d\n", a, b, a+b))
cat(sprintf("   Non-target             %7d    %11d    %5d\n", c, d, c+d))
cat(sprintf("   Total                  %7d    %11d    %5d\n", a+c, b+d, n_background))

# Expected under null hypothesis
expected <- (length(hormone_in_background) / n_background) * n_targets
fold_enrichment <- ifelse(expected > 0, a / expected, NA)

cat("\n   Expected hormone genes among targets:", round(expected, 1), "\n")
cat("   Observed hormone genes among targets:", a, "\n")
cat("   Fold enrichment:", round(fold_enrichment, 2), "\n")

# Fisher's exact test
contingency <- matrix(c(a, c, b, d), nrow = 2)
fisher_result <- fisher.test(contingency, alternative = "greater")
fisher_depletion <- fisher.test(contingency, alternative = "less")

cat("\n   Fisher's exact test (enrichment):\n")
cat("   P-value:", format(fisher_result$p.value, digits = 4), "\n")
cat("   Odds ratio:", round(fisher_result$estimate, 2), "\n")
cat("   P-value (depletion test):", format(fisher_depletion$p.value, digits = 4), "\n")

# ------------------------------------------------------------------------------
# 4. Identify Specific Hormone Genes Among JAG1 Targets
# ------------------------------------------------------------------------------
cat("\n4. Identifying hormone genes among JAG1 targets...\n")

if (length(targets_in_hormone) > 0) {

  hormone_target_details <- jag1_targets %>%
    filter(GeneID %in% targets_in_hormone) %>%
    left_join(hormone_gene_details %>% select(GlymaID, Matched_KO, Arabi_defline),
              by = c("GeneID" = "GlymaID")) %>%
    select(GeneID, Confidence_Tier, Matched_KO,
           Mean_logFC_Pairwise, NarrowvsBroad_TP0_logFC, NarrowvsBroad_TP0_FDR,
           Best_hit_arabi_defline, KEGG_Arabi_defline = Arabi_defline) %>%
    arrange(desc(abs(Mean_logFC_Pairwise)))

  write_csv(hormone_target_details,
            "03_results/tables/hormone_analysis/37b_hormone_genes_in_JAG1_targets.csv")

  cat("   Found", nrow(hormone_target_details), "hormone genes among JAG1 targets\n")

  # Summary by confidence tier
  tier_summary <- hormone_target_details %>%
    group_by(Confidence_Tier) %>%
    summarise(
      N = n(),
      Mean_logFC = round(mean(Mean_logFC_Pairwise, na.rm = TRUE), 2),
      Pct_Up = round(100 * mean(Mean_logFC_Pairwise > 0, na.rm = TRUE), 1),
      .groups = "drop"
    )

  cat("\n   By confidence tier:\n")
  print(as.data.frame(tier_summary))

} else {
  hormone_target_details <- data.frame()
  cat("   No hormone genes found among JAG1 targets\n")
}

# ------------------------------------------------------------------------------
# 5. Compare with DE Analysis (from 37a)
# ------------------------------------------------------------------------------
cat("\n5. Comparing with DE hormone genes...\n")

# Hormone genes by category
hormone_categories <- data.frame(
  Category = c("All hormone genes", "DE hormone genes", "JAG1 target hormone genes",
               "Both DE and JAG1 target"),
  Count = c(
    length(hormone_in_background),
    sum(de_hormone_genes$Significant),
    length(targets_in_hormone),
    length(intersect(targets_in_hormone, de_hormone_genes$Gene[de_hormone_genes$Significant]))
  )
)

cat("\n   Hormone gene categories:\n")
print(hormone_categories)

write_csv(hormone_categories,
          "03_results/tables/hormone_analysis/37b_hormone_category_comparison.csv")

# ------------------------------------------------------------------------------
# 6. Create Results Summary Table
# ------------------------------------------------------------------------------
cat("\n6. Creating summary tables...\n")

enrichment_summary <- data.frame(
  Pathway = "KEGG gmx04075 (Plant hormone signal transduction)",
  Hormone_genes_in_background = length(hormone_in_background),
  JAG1_targets_total = n_targets,
  JAG1_targets_in_pathway = a,
  Expected = round(expected, 2),
  Fold_enrichment = round(fold_enrichment, 2),
  P_value_enrichment = fisher_result$p.value,
  P_value_depletion = fisher_depletion$p.value,
  Odds_ratio = round(fisher_result$estimate, 2),
  Interpretation = ifelse(fisher_result$p.value < 0.05, "Significantly enriched",
                          ifelse(fisher_depletion$p.value < 0.05, "Significantly depleted",
                                 "Not significantly different from background"))
)

write_csv(enrichment_summary,
          "03_results/tables/hormone_analysis/37b_JAG1_hormone_enrichment_results.csv")

# ------------------------------------------------------------------------------
# 7. Visualization
# ------------------------------------------------------------------------------
cat("\n7. Creating visualizations...\n")

# Plot 1: Enrichment visualization
enrichment_plot_data <- data.frame(
  Category = c("Expected", "Observed"),
  Count = c(expected, a),
  Type = c("Expected\n(by chance)", "Observed\n(actual)")
)

p1 <- ggplot(enrichment_plot_data, aes(x = Type, y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = round(Count, 1)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Expected\n(by chance)" = "grey60",
                                "Observed\n(actual)" = "darkgreen")) +
  labs(
    title = "Hormone Pathway Genes Among JAG1 Targets",
    subtitle = paste0("KEGG gmx04075 | Fold enrichment: ", round(fold_enrichment, 2),
                     " | P = ", format(fisher_result$p.value, digits = 3)),
    x = "",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "none",
    axis.text = element_text(size = 12)
  ) +
  expand_limits(y = max(expected, a) * 1.15)

ggsave("03_results/figures/37b_hormone_JAG1/JAG1_hormone_enrichment_barplot.png",
       p1, width = 6, height = 6, dpi = 150)
ggsave("03_results/figures/37b_hormone_JAG1/JAG1_hormone_enrichment_barplot.pdf",
       p1, width = 6, height = 6)
cat("   Saved: JAG1_hormone_enrichment_barplot.png/pdf\n")

# Plot 2: Expression changes of hormone genes in JAG1 targets
if (nrow(hormone_target_details) >= 3) {

  plot_data <- hormone_target_details %>%
    arrange(desc(abs(Mean_logFC_Pairwise))) %>%
    head(30)

  p2 <- plot_data %>%
    ggplot(aes(x = reorder(GeneID, Mean_logFC_Pairwise),
               y = Mean_logFC_Pairwise,
               fill = Confidence_Tier)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#A0A0A0",
                                  "Bronze" = "#CD7F32")) +
    coord_flip() +
    labs(
      title = "Hormone Genes in JAG1 Targets - Expression Changes",
      subtitle = paste0("Top ", nrow(plot_data), " genes by |logFC| | Positive = upregulated in narrow leaves"),
      x = "Gene ID",
      y = "Mean log2 Fold Change (Narrow vs Broad)",
      fill = "Confidence"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 7)
    )

  plot_height <- max(6, min(12, nrow(plot_data) * 0.3))

  ggsave("03_results/figures/37b_hormone_JAG1/JAG1_hormone_gene_expression.png",
         p2, width = 10, height = plot_height, dpi = 150)
  ggsave("03_results/figures/37b_hormone_JAG1/JAG1_hormone_gene_expression.pdf",
         p2, width = 10, height = plot_height)
  cat("   Saved: JAG1_hormone_gene_expression.png/pdf\n")
}

# Plot 3: Comparison - DE vs JAG1 targets (Venn-style bar chart)
comparison_data <- data.frame(
  Category = c("DE only", "JAG1 target only", "Both DE & JAG1 target"),
  Count = c(
    sum(de_hormone_genes$Significant) - length(intersect(targets_in_hormone,
                                                          de_hormone_genes$Gene[de_hormone_genes$Significant])),
    length(targets_in_hormone) - length(intersect(targets_in_hormone,
                                                   de_hormone_genes$Gene[de_hormone_genes$Significant])),
    length(intersect(targets_in_hormone, de_hormone_genes$Gene[de_hormone_genes$Significant]))
  )
)

p3 <- ggplot(comparison_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("DE only" = "steelblue",
                                "JAG1 target only" = "darkgreen",
                                "Both DE & JAG1 target" = "purple")) +
  labs(
    title = "Hormone Genes: DE Status vs JAG1 Target Status",
    subtitle = "Comparing differentially expressed genes with JAG1 binding targets",
    x = "",
    y = "Number of Hormone Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  ) +
  expand_limits(y = max(comparison_data$Count) * 1.15)

ggsave("03_results/figures/37b_hormone_JAG1/hormone_DE_vs_JAG1_comparison.png",
       p3, width = 8, height = 6, dpi = 150)
ggsave("03_results/figures/37b_hormone_JAG1/hormone_DE_vs_JAG1_comparison.pdf",
       p3, width = 8, height = 6)
cat("   Saved: hormone_DE_vs_JAG1_comparison.png/pdf\n")

# ------------------------------------------------------------------------------
# 8. Save Checkpoint
# ------------------------------------------------------------------------------
cat("\n8. Saving checkpoint...\n")

hormone_jag1_analysis <- list(
  enrichment_summary = enrichment_summary,
  hormone_target_details = hormone_target_details,
  targets_in_hormone = targets_in_hormone,
  hormone_categories = hormone_categories,
  fisher_result = fisher_result,
  fold_enrichment = fold_enrichment
)

save(hormone_jag1_analysis, hormone_target_details,
     file = "03_results/checkpoints/37b_hormone_JAG1.RData")
cat("   Saved: 37b_hormone_JAG1.RData\n")

# ------------------------------------------------------------------------------
# 9. Print Summary Report
# ------------------------------------------------------------------------------
cat("\n")
cat("=======================================================================\n")
cat("           SCRIPT 37b SUMMARY: JAG1 TARGET HORMONE ENRICHMENT\n")
cat("=======================================================================\n\n")

cat("QUESTION: Is JAG1 directly responsible for regulating hormone genes?\n\n")

cat("RESULTS:\n")
cat("  Total JAG1 targets:", n_targets, "\n")
cat("  Hormone genes among targets:", a, "\n")
cat("  Expected by chance:", round(expected, 1), "\n")
cat("  Fold enrichment:", round(fold_enrichment, 2), "\n")
cat("  P-value (enrichment):", format(fisher_result$p.value, digits = 4), "\n")
cat("  P-value (depletion):", format(fisher_depletion$p.value, digits = 4), "\n\n")

cat("INTERPRETATION:\n")
if (fisher_result$p.value < 0.05) {
  cat("  *** JAG1 targets ARE significantly enriched for hormone genes ***\n")
  cat("  JAG1 appears to directly regulate hormone signaling.\n")
} else if (fisher_depletion$p.value < 0.05) {
  cat("  *** JAG1 targets are significantly DEPLETED for hormone genes ***\n")
  cat("  JAG1 actively avoids hormone pathway genes.\n")
} else {
  cat("  JAG1 targets are NOT significantly enriched for hormone genes.\n")
  cat("  Hormone genes are present at approximately background frequency.\n")
  cat("  --> JAG1's primary function is NOT hormone regulation.\n")
}

if (nrow(hormone_target_details) > 0) {
  cat("\n  Nonetheless,", nrow(hormone_target_details), "hormone genes ARE among JAG1 targets:\n")
  top_hormone <- hormone_target_details %>%
    arrange(desc(abs(Mean_logFC_Pairwise))) %>%
    head(5)

  for (i in 1:nrow(top_hormone)) {
    direction <- ifelse(top_hormone$Mean_logFC_Pairwise[i] > 0, "UP", "DOWN")
    cat(sprintf("    - %s (%s): %s, logFC = %.2f\n",
                top_hormone$GeneID[i],
                top_hormone$Confidence_Tier[i],
                direction,
                top_hormone$Mean_logFC_Pairwise[i]))
  }
}

cat("\nOUTPUT FILES:\n")
cat("  Tables:\n")
cat("    - 03_results/tables/hormone_analysis/37b_JAG1_hormone_enrichment_results.csv\n")
cat("    - 03_results/tables/hormone_analysis/37b_hormone_genes_in_JAG1_targets.csv\n")
cat("    - 03_results/tables/hormone_analysis/37b_hormone_category_comparison.csv\n")
cat("  Figures:\n")
cat("    - 03_results/figures/37b_hormone_JAG1/JAG1_hormone_enrichment_barplot.png/pdf\n")
cat("    - 03_results/figures/37b_hormone_JAG1/JAG1_hormone_gene_expression.png/pdf\n")
cat("    - 03_results/figures/37b_hormone_JAG1/hormone_DE_vs_JAG1_comparison.png/pdf\n")
cat("  Checkpoint:\n")
cat("    - 03_results/checkpoints/37b_hormone_JAG1.RData\n")

cat("\n=======================================================================\n")
cat("Script 37b complete. Next: Run 37c for motif analysis of hormone genes.\n")
cat("=======================================================================\n")
