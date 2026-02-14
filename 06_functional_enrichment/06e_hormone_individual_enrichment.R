#!/usr/bin/env Rscript
# ==============================================================================
# Script 37e: Individual Hormone Pathway Enrichment Analysis
# ==============================================================================
# Purpose: Test enrichment for EACH hormone pathway separately
#          (auxin, cytokinin, JA, SA, ABA, ethylene, brassinosteroid, gibberellin)
#
# Rationale: The combined pathway analysis (37a, 37b) may miss pathway-specific
#            effects if one hormone is enriched but others dilute the signal.
#
# Input Files:
#   - 03_results/checkpoints/37a_hormone_DE.RData: Hormone mappings and DE data
#   - 03_results/checkpoints/37b_hormone_JAG1.RData: JAG1 target data
#   - 01_data/gmx04075.txt: KEGG pathway file (for hormone classification)
#
# Output:
#   - Enrichment statistics for each individual hormone pathway
#   - Comparison across hormones
#   - Visualization
#
# Dependency: Run AFTER 37a and 37b
# ==============================================================================

cat("=======================================================================\n")
cat("Script 37e: Individual Hormone Pathway Enrichment Analysis\n")
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
cat("1. Loading packages and previous results...\n")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Create output directories
dir.create("03_results/tables/hormone_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("03_results/figures/37e_individual_hormones", showWarnings = FALSE, recursive = TRUE)

# Load checkpoints
if (!file.exists("03_results/checkpoints/37a_hormone_DE.RData")) {
  stop("Checkpoint from 37a not found. Please run 37a first.")
}
load("03_results/checkpoints/37a_hormone_DE.RData")
cat("   Loaded 37a checkpoint\n")

if (!file.exists("03_results/checkpoints/37b_hormone_JAG1.RData")) {
  stop("Checkpoint from 37b not found. Please run 37b first.")
}
load("03_results/checkpoints/37b_hormone_JAG1.RData")
cat("   Loaded 37b checkpoint\n")

# Extract key data
background_genes <- hormone_de_analysis$background_genes
de_genes <- hormone_de_analysis$de_genes
de_genes_up <- hormone_de_analysis$de_genes_up
de_genes_down <- hormone_de_analysis$de_genes_down

# Load JAG1 targets
jag1_targets <- read_csv("03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv",
                         show_col_types = FALSE)
target_genes <- jag1_targets$GeneID

cat("   Background genes:", length(background_genes), "\n")
cat("   DE genes:", length(de_genes), "\n")
cat("   JAG1 targets:", length(target_genes), "\n")

# ------------------------------------------------------------------------------
# 2. Define Individual Hormone Pathways Based on KEGG KO IDs
# ------------------------------------------------------------------------------
cat("\n2. Defining individual hormone pathways...\n")

# KEGG KO IDs for each hormone pathway (from gmx04075)
# Reference: https://www.kegg.jp/pathway/gmx04075

hormone_ko_mapping <- list(
  Auxin = c(
    "K13946",  # AUX1 - auxin influx carrier
    "K14484",  # AUX/IAA - auxin-responsive protein
    "K14485",  # TIR1 - transport inhibitor response 1
    "K14486",  # ARF - auxin response factor
    "K14487",  # GH3 - auxin responsive GH3
    "K14488",  # SAUR - small auxin up RNA
    "K14489"   # AHP - histidine phosphotransfer (shared with cytokinin)
  ),
  Cytokinin = c(
    "K14489",  # AHP - histidine phosphotransfer
    "K14490",  # AHP - histidine phosphotransfer protein
    "K14491",  # CRE1 - cytokinin receptor
    "K14492",  # B-ARR - type-B response regulator
    "K14493"   # A-ARR - type-A response regulator
  ),
  Gibberellin = c(
    "K14494",  # GID1 - gibberellin receptor
    "K14495"   # DELLA - DELLA protein
  ),
  ABA = c(
    "K14496",  # PYR/PYL - abscisic acid receptor
    "K14497",  # PP2C - protein phosphatase 2C
    "K14498",  # SnRK2 - SNF1-related protein kinase 2
    "K14432"   # ABF/AREB - ABA responsive element binding factor
  ),
  Ethylene = c(
    "K14509",  # ETR - ethylene receptor
    "K14510",  # CTR1 - constitutive triple response 1
    "K14511",  # MPK6 - mitogen-activated protein kinase 6
    "K14512",  # EIN2 - ethylene insensitive 2
    "K14513",  # EIN3 - ethylene insensitive 3
    "K14514"   # EBF1/2 - EIN3-binding F-box protein
  ),
  Brassinosteroid = c(
    "K13416",  # BAK1 - BRI1-associated receptor kinase 1
    "K13415",  # BRI1 - brassinosteroid insensitive 1
    "K14500",  # BSK - BR-signaling kinase
    "K14501",  # BSU1 - BRI1 suppressor 1
    "K14502",  # BIN2 - brassinosteroid insensitive 2
    "K14503",  # BZR1/BES1 - brassinazole resistant 1
    "K14504"   # TCH4 - xyloglucan endotransglucosylase
  ),
  Jasmonic_acid = c(
    "K13464",  # JAZ - jasmonate ZIM domain
    "K14505",  # COI1 - coronatine insensitive 1
    "K14506",  # JAR1 - jasmonic acid-amino synthetase
    "K14507",  # MYC2 - transcription factor MYC2
    "K14508"   # ORCA3 - AP2-domain DNA binding protein
  ),
  Salicylic_acid = c(
    "K14431",  # TGA - transcription factor TGA
    "K14432",  # NPR1 - regulatory protein NPR1 (corrected: K14508 is ORCA3/JA, not NPR1)
    "K13449"   # PR1 - pathogenesis-related protein 1
  )
)

cat("   Defined", length(hormone_ko_mapping), "hormone pathways\n")
for (h in names(hormone_ko_mapping)) {
  cat("     ", h, ":", length(hormone_ko_mapping[[h]]), "KO IDs\n")
}

# ------------------------------------------------------------------------------
# 3. Load Annotation and Map Genes to Individual Pathways
# ------------------------------------------------------------------------------
cat("\n3. Mapping genes to individual hormone pathways...\n")

# Load annotation file
annotation_file <- "01_data/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt"
annotations <- read_tsv(annotation_file,
                        col_names = c("pacId", "locusName", "transcriptName", "peptideName",
                                     "Pfam", "Panther", "ec", "KOG", "KO", "GO",
                                     "Arabi_name", "Arabi_defline",
                                     "Clamy_name", "Clamy_defline",
                                     "Rice_name", "Rice_defline"),
                        skip = 1,
                        show_col_types = FALSE)

# Get gene-level KO annotations
gene_ko <- annotations %>%
  group_by(locusName) %>%
  summarise(
    KO = paste(unique(na.omit(KO[KO != ""])), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(KO = na_if(KO, "")) %>%
  filter(!is.na(KO))

# Function to find genes with specific KO IDs
find_genes_by_ko <- function(ko_ids, gene_ko_df, background) {
  pattern <- paste(ko_ids, collapse = "|")
  genes <- gene_ko_df %>%
    filter(grepl(pattern, KO)) %>%
    pull(locusName)
  intersect(genes, background)
}

# Map genes for each hormone pathway
hormone_genes <- list()
for (hormone in names(hormone_ko_mapping)) {
  hormone_genes[[hormone]] <- find_genes_by_ko(
    hormone_ko_mapping[[hormone]],
    gene_ko,
    background_genes
  )
  cat("   ", hormone, ":", length(hormone_genes[[hormone]]), "genes\n")
}

# ------------------------------------------------------------------------------
# 4. Enrichment Analysis Function
# ------------------------------------------------------------------------------
cat("\n4. Performing enrichment analysis for each hormone pathway...\n")

perform_enrichment <- function(test_genes, pathway_genes, background_genes, pathway_name) {
  n_test <- length(test_genes)
  n_background <- length(background_genes)
  n_pathway_bg <- length(pathway_genes)

  # Overlap
  overlap <- length(intersect(test_genes, pathway_genes))

  # Expected
  expected <- (n_pathway_bg / n_background) * n_test

  # Fold enrichment
  fold <- ifelse(expected > 0, overlap / expected, NA)

  # Fisher's test (only if we have genes in the pathway)
  if (n_pathway_bg > 0 && n_test > 0) {
    a <- overlap
    b <- n_test - a
    c <- n_pathway_bg - a
    d <- n_background - n_test - c

    # Ensure no negative values
    if (all(c(a, b, c, d) >= 0)) {
      fisher_result <- fisher.test(matrix(c(a, c, b, d), nrow = 2), alternative = "greater")
      p_value <- fisher_result$p.value
      odds_ratio <- fisher_result$estimate
    } else {
      p_value <- NA
      odds_ratio <- NA
    }
  } else {
    p_value <- NA
    odds_ratio <- NA
  }

  data.frame(
    Hormone = pathway_name,
    Genes_in_background = n_pathway_bg,
    Genes_in_test_set = overlap,
    Expected = round(expected, 2),
    Fold_enrichment = round(fold, 2),
    P_value = p_value,
    Odds_ratio = round(odds_ratio, 2),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Test Enrichment in DE Genes
# ------------------------------------------------------------------------------
cat("\n5. Testing enrichment in DE genes (Narrow vs Broad)...\n")

de_enrichment <- bind_rows(lapply(names(hormone_genes), function(h) {
  perform_enrichment(de_genes, hormone_genes[[h]], background_genes, h)
}))

# Add FDR correction
de_enrichment$FDR <- p.adjust(de_enrichment$P_value, method = "BH")
de_enrichment$Significant <- de_enrichment$FDR < 0.05

# Sort by p-value
de_enrichment <- de_enrichment %>% arrange(P_value)

cat("\n   DE Gene Enrichment Results:\n")
print(as.data.frame(de_enrichment[, c("Hormone", "Genes_in_test_set", "Expected",
                                       "Fold_enrichment", "P_value", "FDR")]))

write_csv(de_enrichment,
          "03_results/tables/hormone_analysis/37e_DE_individual_hormone_enrichment.csv")

# ------------------------------------------------------------------------------
# 6. Test Enrichment in JAG1 Targets
# ------------------------------------------------------------------------------
cat("\n6. Testing enrichment in JAG1 targets...\n")

jag1_enrichment <- bind_rows(lapply(names(hormone_genes), function(h) {
  perform_enrichment(target_genes, hormone_genes[[h]], background_genes, h)
}))

# Add FDR correction
jag1_enrichment$FDR <- p.adjust(jag1_enrichment$P_value, method = "BH")
jag1_enrichment$Significant <- jag1_enrichment$FDR < 0.05

# Sort by p-value
jag1_enrichment <- jag1_enrichment %>% arrange(P_value)

cat("\n   JAG1 Target Enrichment Results:\n")
print(as.data.frame(jag1_enrichment[, c("Hormone", "Genes_in_test_set", "Expected",
                                         "Fold_enrichment", "P_value", "FDR")]))

write_csv(jag1_enrichment,
          "03_results/tables/hormone_analysis/37e_JAG1_individual_hormone_enrichment.csv")

# ------------------------------------------------------------------------------
# 7. Test Enrichment in Upregulated vs Downregulated DE Genes
# ------------------------------------------------------------------------------
cat("\n7. Testing enrichment in up/down-regulated genes separately...\n")

de_up_enrichment <- bind_rows(lapply(names(hormone_genes), function(h) {
  result <- perform_enrichment(de_genes_up, hormone_genes[[h]], background_genes, h)
  result$Direction <- "Upregulated"
  result
}))

de_down_enrichment <- bind_rows(lapply(names(hormone_genes), function(h) {
  result <- perform_enrichment(de_genes_down, hormone_genes[[h]], background_genes, h)
  result$Direction <- "Downregulated"
  result
}))

de_direction_enrichment <- bind_rows(de_up_enrichment, de_down_enrichment)
de_direction_enrichment$FDR <- p.adjust(de_direction_enrichment$P_value, method = "BH")
de_direction_enrichment$Significant <- de_direction_enrichment$FDR < 0.05

write_csv(de_direction_enrichment,
          "03_results/tables/hormone_analysis/37e_DE_direction_hormone_enrichment.csv")

# Print summary
cat("\n   Upregulated genes:\n")
up_summary <- de_up_enrichment %>% arrange(P_value) %>% head(3)
for (i in 1:nrow(up_summary)) {
  cat(sprintf("     %s: %d genes, fold=%.2f, P=%.4f\n",
              up_summary$Hormone[i], up_summary$Genes_in_test_set[i],
              up_summary$Fold_enrichment[i], up_summary$P_value[i]))
}

cat("\n   Downregulated genes:\n")
down_summary <- de_down_enrichment %>% arrange(P_value) %>% head(3)
for (i in 1:nrow(down_summary)) {
  cat(sprintf("     %s: %d genes, fold=%.2f, P=%.4f\n",
              down_summary$Hormone[i], down_summary$Genes_in_test_set[i],
              down_summary$Fold_enrichment[i], down_summary$P_value[i]))
}

# ------------------------------------------------------------------------------
# 8. Create Combined Summary Table
# ------------------------------------------------------------------------------
cat("\n8. Creating combined summary table...\n")

combined_summary <- de_enrichment %>%
  select(Hormone, DE_genes = Genes_in_test_set, DE_expected = Expected,
         DE_fold = Fold_enrichment, DE_pvalue = P_value) %>%
  left_join(
    jag1_enrichment %>%
      select(Hormone, JAG1_genes = Genes_in_test_set, JAG1_expected = Expected,
             JAG1_fold = Fold_enrichment, JAG1_pvalue = P_value),
    by = "Hormone"
  )

# Add significance markers
combined_summary <- combined_summary %>%
  mutate(
    DE_sig = ifelse(DE_pvalue < 0.05, "*", ""),
    JAG1_sig = ifelse(JAG1_pvalue < 0.05, "*", "")
  )

cat("\n   Combined Summary (DE and JAG1 targets):\n")
print(as.data.frame(combined_summary))

write_csv(combined_summary,
          "03_results/tables/hormone_analysis/37e_combined_hormone_summary.csv")

# ------------------------------------------------------------------------------
# 9. Visualization
# ------------------------------------------------------------------------------
cat("\n9. Creating visualizations...\n")

# Plot 1: Fold enrichment comparison (DE vs JAG1)
plot_data <- combined_summary %>%
  select(Hormone, DE_fold, JAG1_fold) %>%
  pivot_longer(cols = c(DE_fold, JAG1_fold),
               names_to = "Test", values_to = "Fold") %>%
  mutate(Test = ifelse(Test == "DE_fold", "DE Genes", "JAG1 Targets"))

p1 <- ggplot(plot_data, aes(x = reorder(Hormone, -Fold), y = Fold, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("DE Genes" = "steelblue", "JAG1 Targets" = "darkgreen")) +
  labs(
    title = "Individual Hormone Pathway Enrichment",
    subtitle = "Red dashed line = no enrichment (fold = 1)",
    x = "Hormone Pathway",
    y = "Fold Enrichment",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave("03_results/figures/37e_individual_hormones/individual_hormone_enrichment.png",
       p1, width = 10, height = 6, dpi = 150)
ggsave("03_results/figures/37e_individual_hormones/individual_hormone_enrichment.pdf",
       p1, width = 10, height = 6)
cat("   Saved: individual_hormone_enrichment.png/pdf\n")

# Plot 2: P-value heatmap
pvalue_data <- combined_summary %>%
  select(Hormone, DE_pvalue, JAG1_pvalue) %>%
  mutate(
    DE_log10p = -log10(DE_pvalue + 1e-10),
    JAG1_log10p = -log10(JAG1_pvalue + 1e-10)
  ) %>%
  select(Hormone, DE_log10p, JAG1_log10p) %>%
  pivot_longer(cols = c(DE_log10p, JAG1_log10p),
               names_to = "Test", values_to = "neg_log10_p") %>%
  mutate(Test = ifelse(Test == "DE_log10p", "DE Genes", "JAG1 Targets"))

p2 <- ggplot(pvalue_data, aes(x = Test, y = Hormone, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(neg_log10_p, 2)), color = "black", size = 3.5) +
  scale_fill_gradient(low = "white", high = "red",
                      name = "-log10(P)") +
  geom_hline(yintercept = seq(0.5, 8.5, 1), color = "grey80") +
  labs(
    title = "Individual Hormone Pathway Enrichment Significance",
    subtitle = "Higher values = more significant (>1.3 corresponds to P < 0.05)",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

ggsave("03_results/figures/37e_individual_hormones/hormone_pvalue_heatmap.png",
       p2, width = 7, height = 6, dpi = 150)
ggsave("03_results/figures/37e_individual_hormones/hormone_pvalue_heatmap.pdf",
       p2, width = 7, height = 6)
cat("   Saved: hormone_pvalue_heatmap.png/pdf\n")

# Plot 3: Gene counts by hormone
count_data <- combined_summary %>%
  select(Hormone, DE_genes, JAG1_genes) %>%
  pivot_longer(cols = c(DE_genes, JAG1_genes),
               names_to = "Test", values_to = "Count") %>%
  mutate(Test = ifelse(Test == "DE_genes", "DE Genes", "JAG1 Targets"))

p3 <- ggplot(count_data, aes(x = reorder(Hormone, -Count), y = Count, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("DE Genes" = "steelblue", "JAG1 Targets" = "darkgreen")) +
  labs(
    title = "Hormone Genes: Counts in DE and JAG1 Target Sets",
    x = "Hormone Pathway",
    y = "Number of Genes",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  expand_limits(y = max(count_data$Count, na.rm = TRUE) * 1.15)

ggsave("03_results/figures/37e_individual_hormones/hormone_gene_counts.png",
       p3, width = 10, height = 6, dpi = 150)
ggsave("03_results/figures/37e_individual_hormones/hormone_gene_counts.pdf",
       p3, width = 10, height = 6)
cat("   Saved: hormone_gene_counts.png/pdf\n")

# ------------------------------------------------------------------------------
# 10. List Genes by Hormone Pathway
# ------------------------------------------------------------------------------
cat("\n10. Creating gene lists by hormone pathway...\n")

# Get the hormone genes that are DE or JAG1 targets
hormone_gene_lists <- data.frame()

for (hormone in names(hormone_genes)) {
  genes <- hormone_genes[[hormone]]

  for (gene in genes) {
    is_de <- gene %in% de_genes
    is_target <- gene %in% target_genes

    if (is_de || is_target) {
      hormone_gene_lists <- rbind(hormone_gene_lists, data.frame(
        Hormone = hormone,
        GeneID = gene,
        Is_DE = is_de,
        Is_JAG1_target = is_target,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Add expression info
if (nrow(hormone_gene_lists) > 0) {
  de_results <- read_csv("03_results/tables/DE/NarrowvsBroad_TP0.csv", show_col_types = FALSE)

  hormone_gene_lists <- hormone_gene_lists %>%
    left_join(de_results %>% select(Gene, logFC, FDR, Direction),
              by = c("GeneID" = "Gene"))

  write_csv(hormone_gene_lists,
            "03_results/tables/hormone_analysis/37e_hormone_genes_by_pathway.csv")

  cat("   Saved gene lists with", nrow(hormone_gene_lists), "entries\n")
}

# ------------------------------------------------------------------------------
# 11. Save Checkpoint
# ------------------------------------------------------------------------------
cat("\n11. Saving checkpoint...\n")

individual_hormone_analysis <- list(
  hormone_ko_mapping = hormone_ko_mapping,
  hormone_genes = hormone_genes,
  de_enrichment = de_enrichment,
  jag1_enrichment = jag1_enrichment,
  de_direction_enrichment = de_direction_enrichment,
  combined_summary = combined_summary,
  hormone_gene_lists = hormone_gene_lists
)

save(individual_hormone_analysis,
     file = "03_results/checkpoints/37e_individual_hormones.RData")
cat("   Saved: 37e_individual_hormones.RData\n")

# ------------------------------------------------------------------------------
# 12. Print Summary Report
# ------------------------------------------------------------------------------
cat("\n")
cat("=======================================================================\n")
cat("      SCRIPT 37e SUMMARY: INDIVIDUAL HORMONE ENRICHMENT\n")
cat("=======================================================================\n\n")

cat("QUESTION: Are any INDIVIDUAL hormone pathways enriched?\n")
cat("          (The combined analysis may miss pathway-specific signals)\n\n")

cat("HORMONE PATHWAY SIZES (genes in background):\n")
for (h in names(hormone_genes)) {
  cat(sprintf("  %s: %d genes\n", h, length(hormone_genes[[h]])))
}

cat("\n")
cat("ENRICHMENT IN DE GENES (Narrow vs Broad):\n")
cat("=========================================\n")
for (i in 1:nrow(de_enrichment)) {
  sig <- ifelse(de_enrichment$P_value[i] < 0.05, " ***", "")
  cat(sprintf("  %-15s: %2d genes (exp: %5.1f), fold=%5.2f, P=%.4f%s\n",
              de_enrichment$Hormone[i],
              de_enrichment$Genes_in_test_set[i],
              de_enrichment$Expected[i],
              de_enrichment$Fold_enrichment[i],
              de_enrichment$P_value[i],
              sig))
}

cat("\n")
cat("ENRICHMENT IN JAG1 TARGETS:\n")
cat("===========================\n")
for (i in 1:nrow(jag1_enrichment)) {
  sig <- ifelse(jag1_enrichment$P_value[i] < 0.05, " ***", "")
  cat(sprintf("  %-15s: %2d genes (exp: %5.1f), fold=%5.2f, P=%.4f%s\n",
              jag1_enrichment$Hormone[i],
              jag1_enrichment$Genes_in_test_set[i],
              jag1_enrichment$Expected[i],
              jag1_enrichment$Fold_enrichment[i],
              jag1_enrichment$P_value[i],
              sig))
}

cat("\n")
cat("KEY FINDINGS:\n")
cat("=============\n")

# Check for significant enrichments
sig_de <- de_enrichment %>% filter(P_value < 0.05)
sig_jag1 <- jag1_enrichment %>% filter(P_value < 0.05)

if (nrow(sig_de) > 0) {
  cat("  *** SIGNIFICANT enrichment in DE genes:\n")
  for (i in 1:nrow(sig_de)) {
    cat(sprintf("      %s: fold=%.2f, P=%.4f\n",
                sig_de$Hormone[i], sig_de$Fold_enrichment[i], sig_de$P_value[i]))
  }
} else {
  cat("  No individual hormone pathway is significantly enriched in DE genes.\n")
}

if (nrow(sig_jag1) > 0) {
  cat("  *** SIGNIFICANT enrichment in JAG1 targets:\n")
  for (i in 1:nrow(sig_jag1)) {
    cat(sprintf("      %s: fold=%.2f, P=%.4f\n",
                sig_jag1$Hormone[i], sig_jag1$Fold_enrichment[i], sig_jag1$P_value[i]))
  }
} else {
  cat("  No individual hormone pathway is significantly enriched in JAG1 targets.\n")
}

# Top pathway
top_de <- de_enrichment %>% slice(1)
top_jag1 <- jag1_enrichment %>% slice(1)

cat("\n  Closest to significance:\n")
cat(sprintf("    DE genes: %s (fold=%.2f, P=%.4f)\n",
            top_de$Hormone, top_de$Fold_enrichment, top_de$P_value))
cat(sprintf("    JAG1 targets: %s (fold=%.2f, P=%.4f)\n",
            top_jag1$Hormone, top_jag1$Fold_enrichment, top_jag1$P_value))

cat("\nINTERPRETATION:\n")
if (nrow(sig_de) == 0 && nrow(sig_jag1) == 0) {
  cat("  Even at the individual hormone level, no pathway shows significant\n")
  cat("  enrichment. This confirms that hormone signaling is NOT a primary\n")
  cat("  function of JAG1 - the result is consistent across all hormone types.\n")
} else {
  cat("  Individual hormone pathway analysis reveals pathway-specific effects\n")
  cat("  that were masked in the combined analysis.\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Tables:\n")
cat("    - 37e_DE_individual_hormone_enrichment.csv\n")
cat("    - 37e_JAG1_individual_hormone_enrichment.csv\n")
cat("    - 37e_DE_direction_hormone_enrichment.csv\n")
cat("    - 37e_combined_hormone_summary.csv\n")
cat("    - 37e_hormone_genes_by_pathway.csv\n")
cat("  Figures:\n")
cat("    - individual_hormone_enrichment.png/pdf\n")
cat("    - hormone_pvalue_heatmap.png/pdf\n")
cat("    - hormone_gene_counts.png/pdf\n")

cat("\n=======================================================================\n")
cat("Script 37e complete.\n")
cat("=======================================================================\n")
