# 06f: Functional Integration
# Integrates GO enrichment, KEGG pathways, and WGCNA modules into a master gene table
# Inputs: 22_GO_enrichment.RData, 23_KEGG_pathways.RData, 21_WGCNA_JAG1.RData, 14_JAG1_targets.RData
# Output: 24_functional_integrated.RData

rm(list = ls())
gc()

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

dir.create("03_results/figures/24_integration", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(WGCNA)

# Load data
load("03_results/checkpoints/22_GO_enrichment.RData")
cat("Loaded GO enrichment data\n")

load("03_results/checkpoints/23_KEGG_pathways.RData")
cat("Loaded KEGG pathway data\n")

wgcna_available <- FALSE
if (file.exists("03_results/checkpoints/21_WGCNA_JAG1.RData")) {
  load("03_results/checkpoints/21_WGCNA_JAG1.RData")
  wgcna_available <- TRUE
  cat("Loaded WGCNA JAG1 data\n")
}

load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 target data\n\n")

jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]

# Create master gene table
master_table <- data.frame(
  GeneID = jag1_targets$GeneID,
  Confidence_Tier = jag1_targets$Confidence_Tier,
  stringsAsFactors = FALSE
)

if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  master_table$Mean_logFC <- jag1_targets$Mean_logFC_Pairwise
}

# Add functional category from KEGG
if (!is.null(kegg_analysis$detailed_results)) {
  all_cats <- do.call(rbind, kegg_analysis$detailed_results)
  cat_mapping <- all_cats[!duplicated(all_cats$GeneID), c("GeneID", "Category")]
  master_table <- merge(master_table, cat_mapping, by = "GeneID", all.x = TRUE)
  master_table$Category[is.na(master_table$Category)] <- "Uncategorized"
}

# Add WGCNA module
if (wgcna_available && exists("target_modules")) {
  module_df <- data.frame(
    GeneID = names(target_modules),
    WGCNA_Module = target_modules,
    stringsAsFactors = FALSE
  )
  module_df$Module_Color <- labels2colors(module_df$WGCNA_Module)
  master_table <- merge(master_table, module_df, by = "GeneID", all.x = TRUE)
}

cat("Master table created:", nrow(master_table), "genes\n")
cat("Columns:", paste(colnames(master_table), collapse = ", "), "\n\n")

# Count evidence types per gene
master_table$Evidence_Count <- 0

master_table$Has_DE <- TRUE
master_table$Evidence_Count <- master_table$Evidence_Count + 1

if ("Category" %in% colnames(master_table)) {
  master_table$Has_Function <- master_table$Category != "Uncategorized"
  master_table$Evidence_Count <- master_table$Evidence_Count + as.numeric(master_table$Has_Function)
}

if ("WGCNA_Module" %in% colnames(master_table)) {
  master_table$Has_Module <- !is.na(master_table$WGCNA_Module) & master_table$WGCNA_Module != 0
  master_table$Evidence_Count <- master_table$Evidence_Count + as.numeric(master_table$Has_Module)
}

cat("Evidence integration summary:\n")
cat("  Genes with DE evidence:", sum(master_table$Has_DE), "\n")
if ("Has_Function" %in% colnames(master_table)) {
  cat("  Genes with functional annotation:", sum(master_table$Has_Function), "\n")
}
if ("Has_Module" %in% colnames(master_table)) {
  cat("  Genes in WGCNA modules:", sum(master_table$Has_Module, na.rm = TRUE), "\n")
}

evidence_dist <- table(master_table$Evidence_Count)
cat("\nEvidence count distribution:\n")
print(evidence_dist)

# Tier vs Function crosstab
if ("Category" %in% colnames(master_table)) {
  tier_function <- table(master_table$Confidence_Tier, master_table$Category)
  cat("\nConfidence Tier vs Functional Category:\n")
  print(tier_function)

  tier_function_df <- as.data.frame.matrix(tier_function)
  tier_function_df$Tier <- rownames(tier_function_df)
  tier_function_df <- tier_function_df[, c("Tier", setdiff(colnames(tier_function_df), "Tier"))]
  write.csv(tier_function_df, "03_results/tables/functional/tier_vs_function.csv",
            row.names = FALSE)

  png("03_results/figures/24_integration/tier_function_heatmap.png",
      width = 1000, height = 600, res = 120)
  pheatmap(tier_function,
           main = "JAG1 Targets: Confidence Tier vs Functional Category",
           color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
           display_numbers = TRUE,
           number_format = "%d",
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           fontsize = 10)
  dev.off()
  cat("Saved: tier_function_heatmap.png\n")
}

# Module vs Function
if (wgcna_available && "WGCNA_Module" %in% colnames(master_table) && "Category" %in% colnames(master_table)) {
  module_genes <- master_table[!is.na(master_table$WGCNA_Module) & master_table$WGCNA_Module != 0, ]

  if (nrow(module_genes) > 0) {
    module_function <- table(module_genes$Module_Color, module_genes$Category)
    module_function <- module_function[rowSums(module_function) >= 5, ]

    if (nrow(module_function) > 0) {
      cat("\nModule vs Functional Category (top modules):\n")
      print(module_function[1:min(10, nrow(module_function)), ])

      png("03_results/figures/24_integration/module_function_heatmap.png",
          width = 1000, height = 800, res = 120)
      module_function_norm <- module_function / rowSums(module_function) * 100
      pheatmap(module_function_norm,
               main = "WGCNA Module vs Functional Category (%)",
               color = colorRampPalette(c("white", "orange", "red"))(50),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               fontsize = 9)
      dev.off()
      cat("Saved: module_function_heatmap.png\n")
    }
  }
}

# Priority gene list
master_table$Priority_Score <- 0
tier_weights <- c("Gold" = 3, "Silver" = 2, "Bronze" = 1)
master_table$Priority_Score <- master_table$Priority_Score +
  tier_weights[master_table$Confidence_Tier]
master_table$Priority_Score <- master_table$Priority_Score + master_table$Evidence_Count

if ("Mean_logFC" %in% colnames(master_table)) {
  master_table$Priority_Score <- master_table$Priority_Score +
    pmin(abs(master_table$Mean_logFC), 3)
}

master_table <- master_table[order(-master_table$Priority_Score), ]

cat("\nTop 20 priority genes:\n")
priority_cols <- intersect(c("GeneID", "Confidence_Tier", "Category", "Module_Color",
                             "Mean_logFC", "Priority_Score"), colnames(master_table))
print(head(master_table[, priority_cols], 20))

# Save results
write.csv(master_table, "03_results/tables/functional/integrated_gene_table.csv",
          row.names = FALSE)
cat("\nSaved: integrated_gene_table.csv\n")

priority_list <- head(master_table, 100)
write.csv(priority_list, "03_results/tables/functional/priority_targets_top100.csv",
          row.names = FALSE)
cat("Saved: priority_targets_top100.csv\n")

# Summary statistics
summary_stats <- data.frame(
  Metric = c(
    "Total JAG1 targets",
    "Gold tier targets",
    "Silver tier targets",
    "Bronze tier targets",
    "Targets with functional annotation",
    "Targets in WGCNA modules",
    "Multi-evidence targets (>=2)",
    "Priority score range"
  ),
  Value = c(
    nrow(master_table),
    sum(master_table$Confidence_Tier == "Gold"),
    sum(master_table$Confidence_Tier == "Silver"),
    sum(master_table$Confidence_Tier == "Bronze"),
    sum(master_table$Has_Function, na.rm = TRUE),
    sum(master_table$Has_Module, na.rm = TRUE),
    sum(master_table$Evidence_Count >= 2, na.rm = TRUE),
    paste(range(master_table$Priority_Score), collapse = " - ")
  ),
  stringsAsFactors = FALSE
)

cat("\nIntegration Summary:\n")
print(summary_stats)
write.csv(summary_stats, "03_results/tables/functional/integration_summary.csv",
          row.names = FALSE)

# Save checkpoint
functional_integration <- list(
  master_table = master_table,
  summary_stats = summary_stats,
  wgcna_available = wgcna_available
)

save(
  functional_integration,
  master_table,
  summary_stats,
  file = "03_results/checkpoints/24_functional_integrated.RData"
)

cat("Checkpoint saved: 24_functional_integrated.RData\n")
cat("\nIntegrated", nrow(master_table), "JAG1 targets\n")
cat("Top priority gene:", master_table$GeneID[1], "\n")
cat("Priority score range:", min(master_table$Priority_Score), "-",
    max(master_table$Priority_Score), "\n")
