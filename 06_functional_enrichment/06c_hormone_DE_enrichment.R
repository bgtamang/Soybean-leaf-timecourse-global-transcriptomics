#!/usr/bin/env Rscript
#          leaves are enriched for plant hormone pathway genes
#
# Question: Are hormone pathways altered in leaf shape differences?
#           (regardless of whether JAG1 directly regulates them)
#
# Input Files:
#   - 01_data/gmx04075.txt: KEGG pathway file (Plant hormone signal transduction)
#   - 01_data/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt: Phytozome annotation
#   - 03_results/tables/DE/NarrowvsBroad_TP0.csv: DE results
#   - 03_results/checkpoints/05_normalized.RData: Background genes
#
# Output:
#   - Hormone pathway enrichment statistics for DE genes
#   - Comparison with random expectation
#   - Visualization of results
#
# Dependency: Run BEFORE 37b (JAG1 target analysis)
#
# KEGG Reference:
#   - Pathway: gmx04075 (Plant hormone signal transduction - Glycine max)
#   - URL: https://www.kegg.jp/pathway/gmx04075

cat("Script 37a: Hormone Pathway Enrichment in DE Genes (TP0)\n")

# 0. Set Working Directory to Project Root
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
project_root <- dirname(script_dir)
setwd(project_root)
cat("Working directory:", getwd(), "\n\n")

# 1. Setup and Load Dependencies
cat("1. Loading packages...\n")

suppressPackageStartupMessages({
  library(tidyverse)
})

# Create output directories
dir.create("03_results/tables/hormone_analysis", showWarnings = FALSE, recursive = TRUE)
dir.create("03_results/figures/37a_hormone_DE", showWarnings = FALSE, recursive = TRUE)

# 2. Parse KEGG Pathway File (gmx04075.txt)
cat("\n2. Parsing KEGG pathway file...\n")

kegg_file <- "01_data/gmx04075.txt"

if (!file.exists(kegg_file)) {
  stop("KEGG pathway file not found: ", kegg_file,
       "\n   Download from: https://rest.kegg.jp/get/gmx04075")
}

# Read the entire file
kegg_lines <- readLines(kegg_file)

# Extract GENE entries (lines starting with spaces followed by NCBI ID)
gene_lines <- kegg_lines[grepl("^\\s+\\d+\\s+", kegg_lines)]

cat("   Found", length(gene_lines), "gene entries in KEGG pathway\n")

# Parse each gene line to extract NCBI ID, Symbol, and KO
parse_kegg_gene <- function(line) {
  ncbi_id <- str_extract(line, "\\d+")
  ko_matches <- str_extract_all(line, "K\\d{5}")[[1]]
  ko_ids <- if (length(ko_matches) > 0) paste(ko_matches, collapse = ";") else NA
  symbol_match <- str_match(line, "\\d+\\s+([A-Za-z0-9_-]+)")
  symbol <- if (!is.na(symbol_match[1,2])) symbol_match[1,2] else NA
  desc_match <- str_match(line, ";\\s*(.+?)\\s*\\[KO:")
  description <- if (!is.na(desc_match[1,2])) desc_match[1,2] else NA

  data.frame(
    NCBI_GeneID = ncbi_id,
    Symbol = symbol,
    KO = ko_ids,
    Description = description,
    stringsAsFactors = FALSE
  )
}

kegg_genes <- bind_rows(lapply(gene_lines, parse_kegg_gene))
cat("   Parsed", nrow(kegg_genes), "genes\n")
cat("   Genes with KO annotation:", sum(!is.na(kegg_genes$KO)), "\n")

# Get unique KO IDs from pathway
kegg_ko_ids <- kegg_genes %>%
  filter(!is.na(KO)) %>%
  pull(KO) %>%
  str_split(";") %>%
  unlist() %>%
  unique()

cat("   Unique KO IDs in pathway:", length(kegg_ko_ids), "\n")

# 3. Load Phytozome Annotation File
cat("\n3. Loading Phytozome annotation file...\n")

annotation_file <- "01_data/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt"

if (!file.exists(annotation_file)) {
  stop("Annotation file not found: ", annotation_file)
}

annotations <- read_tsv(annotation_file,
                        col_names = c("pacId", "locusName", "transcriptName", "peptideName",
                                     "Pfam", "Panther", "ec", "KOG", "KO", "GO",
                                     "Arabi_name", "Arabi_defline",
                                     "Clamy_name", "Clamy_defline",
                                     "Rice_name", "Rice_defline"),
                        skip = 1,
                        show_col_types = FALSE)

cat("   Loaded", nrow(annotations), "transcript annotations\n")

# Get unique gene-level annotations
gene_annotations <- annotations %>%
  group_by(locusName) %>%
  summarise(
    KO = paste(unique(na.omit(KO[KO != ""])), collapse = ";"),
    GO = paste(unique(na.omit(GO[GO != ""])), collapse = " "),
    Arabi_name = first(na.omit(Arabi_name)),
    Arabi_defline = first(na.omit(Arabi_defline)),
    .groups = "drop"
  ) %>%
  mutate(KO = na_if(KO, ""))

cat("   Collapsed to", nrow(gene_annotations), "unique genes\n")
cat("   Genes with KO annotation:", sum(!is.na(gene_annotations$KO)), "\n")

# 4. Map KEGG Hormone Pathway to Glyma IDs via KO
cat("\n4. Mapping KEGG pathway genes to Glyma IDs via KO...\n")

find_genes_with_ko <- function(annotations_df, ko_ids) {
  pattern <- paste(ko_ids, collapse = "|")
  annotations_df %>%
    filter(!is.na(KO)) %>%
    filter(grepl(pattern, KO)) %>%
    pull(locusName)
}

hormone_glyma_ids <- find_genes_with_ko(gene_annotations, kegg_ko_ids)
cat("   Glyma genes with pathway KO IDs:", length(hormone_glyma_ids), "\n")

# Create detailed mapping table
hormone_gene_details <- gene_annotations %>%
  filter(locusName %in% hormone_glyma_ids) %>%
  select(GlymaID = locusName, KO, Arabi_name, Arabi_defline) %>%
  rowwise() %>%
  mutate(
    Matched_KO = {
      gene_kos <- str_split(KO, ";")[[1]]
      matched <- intersect(gene_kos, kegg_ko_ids)
      paste(matched, collapse = ";")
    }
  ) %>%
  ungroup()

# 5. Load DE Genes and Background
cat("\n5. Loading DE genes and background...\n")

# Load DE results for NarrowvsBroad_TP0
de_file <- "03_results/tables/DE/NarrowvsBroad_TP0.csv"
de_results <- read_csv(de_file, show_col_types = FALSE)

# Significant DE genes (FDR < 0.05)
de_genes_all <- de_results %>% filter(Significant == TRUE)
de_genes <- de_genes_all$Gene
de_genes_up <- de_results %>% filter(Significant == TRUE, Direction == "Up") %>% pull(Gene)
de_genes_down <- de_results %>% filter(Significant == TRUE, Direction == "Down") %>% pull(Gene)

cat("   Total DE genes (FDR < 0.05):", length(de_genes), "\n")
cat("   Upregulated in narrow:", length(de_genes_up), "\n")
cat("   Downregulated in narrow:", length(de_genes_down), "\n")

# Load background genes from checkpoint
# Try multiple checkpoint options based on what's available
if (file.exists("03_results/checkpoints/05_normalized.RData")) {
  load("03_results/checkpoints/05_normalized.RData")
  # The checkpoint contains dge_corrected, not y
  if (exists("dge_corrected")) {
    background_genes <- rownames(dge_corrected$counts)
  } else if (exists("v_full")) {
    background_genes <- rownames(v_full$E)
  } else {
    background_genes <- gene_annotations$locusName
    cat("   (Checkpoint loaded but no count matrix found, using annotation)\n")
  }
} else if (file.exists("03_results/checkpoints/03_QC_complete.RData")) {
  load("03_results/checkpoints/03_QC_complete.RData")
  if (exists("dge_filtered")) {
    background_genes <- rownames(dge_filtered$counts)
  } else {
    background_genes <- gene_annotations$locusName
  }
} else {
  background_genes <- gene_annotations$locusName
  cat("   (No checkpoint found, using annotation file as background)\n")
}
cat("   Background genes:", length(background_genes), "\n")

# Filter hormone genes to those in background
hormone_in_background <- intersect(hormone_glyma_ids, background_genes)
cat("   Hormone pathway genes in background:", length(hormone_in_background), "\n")

# 6. Enrichment Analysis (Fisher's Exact Test) - All DE Genes
cat("\n6. Performing enrichment analysis for ALL DE genes...\n")

perform_enrichment_test <- function(test_genes, background, hormone_genes, label) {
  n_test <- length(test_genes)
  n_background <- length(background)

  # Test genes that are hormone pathway genes
  test_in_hormone <- intersect(test_genes, hormone_genes)
  a <- length(test_in_hormone)
  b <- n_test - a
  c <- length(hormone_genes) - a
  d <- n_background - n_test - c

  # Expected under null
  expected <- (length(hormone_genes) / n_background) * n_test
  fold_enrichment <- ifelse(expected > 0, a / expected, NA)

  # Fisher's exact test
  contingency <- matrix(c(a, c, b, d), nrow = 2)
  fisher_enrich <- fisher.test(contingency, alternative = "greater")
  fisher_deplete <- fisher.test(contingency, alternative = "less")

  list(
    Label = label,
    Test_genes = n_test,
    Hormone_in_test = a,
    Expected = round(expected, 2),
    Fold_enrichment = round(fold_enrichment, 2),
    P_enrichment = fisher_enrich$p.value,
    P_depletion = fisher_deplete$p.value,
    Odds_ratio = round(fisher_enrich$estimate, 2),
    Genes = test_in_hormone
  )
}

# Test all DE genes
result_all <- perform_enrichment_test(de_genes, background_genes, hormone_in_background, "All_DE")
result_up <- perform_enrichment_test(de_genes_up, background_genes, hormone_in_background, "Upregulated")
result_down <- perform_enrichment_test(de_genes_down, background_genes, hormone_in_background, "Downregulated")

cat("\n   Results for ALL DE genes:\n")
cat("     Hormone genes in DE:", result_all$Hormone_in_test, "\n")
cat("     Expected by chance:", result_all$Expected, "\n")
cat("     Fold enrichment:", result_all$Fold_enrichment, "\n")
cat("     P-value (enrichment):", format(result_all$P_enrichment, digits = 4), "\n")

cat("\n   Results for UPREGULATED genes:\n")
cat("     Hormone genes:", result_up$Hormone_in_test, "\n")
cat("     Fold enrichment:", result_up$Fold_enrichment, "\n")
cat("     P-value:", format(result_up$P_enrichment, digits = 4), "\n")

cat("\n   Results for DOWNREGULATED genes:\n")
cat("     Hormone genes:", result_down$Hormone_in_test, "\n")
cat("     Fold enrichment:", result_down$Fold_enrichment, "\n")
cat("     P-value:", format(result_down$P_enrichment, digits = 4), "\n")

# 7. Identify Specific Hormone Genes in DE Set
cat("\n7. Identifying hormone genes among DE genes...\n")

de_hormone_genes <- de_results %>%
  filter(Gene %in% hormone_in_background) %>%
  left_join(hormone_gene_details %>% select(GlymaID, Matched_KO, Arabi_defline),
            by = c("Gene" = "GlymaID")) %>%
  arrange(FDR)

# Add DE status
de_hormone_genes <- de_hormone_genes %>%
  mutate(DE_status = case_when(
    Significant & Direction == "Up" ~ "Up_in_Narrow",
    Significant & Direction == "Down" ~ "Down_in_Narrow",
    TRUE ~ "Not_DE"
  ))

write_csv(de_hormone_genes,
          "03_results/tables/hormone_analysis/37a_hormone_genes_DE_status.csv")

cat("   Total hormone genes in dataset:", nrow(de_hormone_genes), "\n")
cat("   DE hormone genes:", sum(de_hormone_genes$Significant), "\n")
cat("   - Upregulated:", sum(de_hormone_genes$DE_status == "Up_in_Narrow"), "\n")
cat("   - Downregulated:", sum(de_hormone_genes$DE_status == "Down_in_Narrow"), "\n")

# 8. Create Results Summary Table
cat("\n8. Creating summary tables...\n")

enrichment_summary <- data.frame(
  Gene_Set = c("All_DE", "Upregulated", "Downregulated"),
  Genes_tested = c(result_all$Test_genes, result_up$Test_genes, result_down$Test_genes),
  Hormone_genes_found = c(result_all$Hormone_in_test, result_up$Hormone_in_test, result_down$Hormone_in_test),
  Expected = c(result_all$Expected, result_up$Expected, result_down$Expected),
  Fold_enrichment = c(result_all$Fold_enrichment, result_up$Fold_enrichment, result_down$Fold_enrichment),
  P_value_enrichment = c(result_all$P_enrichment, result_up$P_enrichment, result_down$P_enrichment),
  P_value_depletion = c(result_all$P_depletion, result_up$P_depletion, result_down$P_depletion),
  Odds_ratio = c(result_all$Odds_ratio, result_up$Odds_ratio, result_down$Odds_ratio)
)

enrichment_summary <- enrichment_summary %>%
  mutate(
    Interpretation = case_when(
      P_value_enrichment < 0.05 ~ "Significantly ENRICHED",
      P_value_depletion < 0.05 ~ "Significantly DEPLETED",
      TRUE ~ "Not significantly different"
    )
  )

write_csv(enrichment_summary,
          "03_results/tables/hormone_analysis/37a_DE_hormone_enrichment_results.csv")

# Print summary
cat("\n   Enrichment Summary:\n")
print(as.data.frame(enrichment_summary[, c("Gene_Set", "Hormone_genes_found", "Expected",
                                            "Fold_enrichment", "P_value_enrichment", "Interpretation")]))

# 9. Visualization
cat("\n9. Creating visualizations...\n")

# Plot 1: Enrichment comparison
plot_data <- enrichment_summary %>%
  pivot_longer(cols = c(Hormone_genes_found, Expected),
               names_to = "Type", values_to = "Count") %>%
  mutate(Type = ifelse(Type == "Hormone_genes_found", "Observed", "Expected"))

p1 <- ggplot(plot_data, aes(x = Gene_Set, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  geom_text(aes(label = round(Count, 1)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Expected" = "grey60", "Observed" = "steelblue")) +
  labs(
    title = "Hormone Pathway Genes in DE Genes (Narrow vs Broad, TP0)",
    subtitle = "KEGG gmx04075 - Plant hormone signal transduction",
    x = "Gene Set",
    y = "Number of Hormone Genes",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text = element_text(size = 11)
  ) +
  expand_limits(y = max(plot_data$Count) * 1.15)

ggsave("03_results/figures/37a_hormone_DE/DE_hormone_enrichment.png",
       p1, width = 8, height = 6, dpi = 150)
ggsave("03_results/figures/37a_hormone_DE/DE_hormone_enrichment.pdf",
       p1, width = 8, height = 6)
cat("   Saved: DE_hormone_enrichment.png/pdf\n")

# Plot 2: Expression of DE hormone genes
de_hormone_sig <- de_hormone_genes %>%
  filter(Significant == TRUE) %>%
  arrange(desc(abs(logFC))) %>%
  head(30)

if (nrow(de_hormone_sig) > 0) {
  p2 <- de_hormone_sig %>%
    ggplot(aes(x = reorder(Gene, logFC), y = logFC, fill = Direction)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    coord_flip() +
    labs(
      title = "DE Hormone Genes (Narrow vs Broad, TP0)",
      subtitle = paste0("Top ", nrow(de_hormone_sig), " genes by |logFC|"),
      x = "Gene ID",
      y = "log2 Fold Change (Narrow / Broad)",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 7)
    )

  plot_height <- max(6, min(12, nrow(de_hormone_sig) * 0.3))
  ggsave("03_results/figures/37a_hormone_DE/DE_hormone_genes_expression.png",
         p2, width = 10, height = plot_height, dpi = 150)
  ggsave("03_results/figures/37a_hormone_DE/DE_hormone_genes_expression.pdf",
         p2, width = 10, height = plot_height)
  cat("   Saved: DE_hormone_genes_expression.png/pdf\n")
}

# 10. Save Checkpoint
cat("\n10. Saving checkpoint...\n")

hormone_de_analysis <- list(
  enrichment_summary = enrichment_summary,
  de_hormone_genes = de_hormone_genes,
  hormone_in_background = hormone_in_background,
  hormone_glyma_ids = hormone_glyma_ids,
  kegg_ko_ids = kegg_ko_ids,
  background_genes = background_genes,
  de_genes = de_genes,
  de_genes_up = de_genes_up,
  de_genes_down = de_genes_down
)

save(hormone_de_analysis, hormone_gene_details,
     file = "03_results/checkpoints/37a_hormone_DE.RData")
cat("   Saved: 37a_hormone_DE.RData\n")

# 11. Print Summary Report
cat("\n")
cat("                SCRIPT 37a SUMMARY: HORMONE DE ENRICHMENT\n")

cat("QUESTION: Are hormone pathway genes differentially expressed between\n")
cat("          narrow and broad leaves at TP0?\n\n")

cat("DATA SOURCES:\n")
cat("  KEGG Pathway: gmx04075 (Plant hormone signal transduction)\n")
cat("  DE comparison: NarrowvsBroad_TP0\n\n")

cat("RESULTS:\n")
cat("  Background genes:", length(background_genes), "\n")
cat("  Hormone genes in background:", length(hormone_in_background), "\n")
cat("  Total DE genes:", length(de_genes), "\n")
cat("  DE hormone genes:", sum(de_hormone_genes$Significant), "\n\n")

for (i in 1:nrow(enrichment_summary)) {
  cat(sprintf("  %s:\n", enrichment_summary$Gene_Set[i]))
  cat(sprintf("    Found: %d, Expected: %.1f, Fold: %.2f, P = %.4f\n",
              enrichment_summary$Hormone_genes_found[i],
              enrichment_summary$Expected[i],
              enrichment_summary$Fold_enrichment[i],
              enrichment_summary$P_value_enrichment[i]))
  cat(sprintf("    --> %s\n\n", enrichment_summary$Interpretation[i]))
}

cat("INTERPRETATION:\n")
if (any(enrichment_summary$P_value_enrichment < 0.05)) {
  cat("  *** Hormone genes ARE significantly enriched among DE genes ***\n")
  cat("  Hormone signaling IS altered between leaf types.\n")
} else {
  cat("  Hormone genes are NOT significantly enriched among DE genes.\n")
  cat("  Hormone pathway genes are present at approximately background frequency.\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Tables:\n")
cat("    - 03_results/tables/hormone_analysis/37a_DE_hormone_enrichment_results.csv\n")
cat("    - 03_results/tables/hormone_analysis/37a_hormone_genes_DE_status.csv\n")
cat("  Figures:\n")
cat("    - 03_results/figures/37a_hormone_DE/DE_hormone_enrichment.png/pdf\n")
cat("    - 03_results/figures/37a_hormone_DE/DE_hormone_genes_expression.png/pdf\n")
cat("  Checkpoint:\n")
cat("    - 03_results/checkpoints/37a_hormone_DE.RData\n")

cat("\n=======================================================================\n")
