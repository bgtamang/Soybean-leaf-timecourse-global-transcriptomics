
rm(list = ls())
gc()

cat("\n")
cat("  SCRIPT 27: CROSS-VALIDATION SUMMARY\n")


base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase2-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "VennDiagram"
)

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
  }
}
cat("  Packages loaded\n\n")
load("03_results/checkpoints/26_validation.RData")
cat("Loaded validation results\n")

load("03_results/checkpoints/14_JAG1_targets.RData")
cat("Loaded JAG1 targets\n")

# Try loading optional checkpoints
if (file.exists("03_results/checkpoints/21_WGCNA_JAG1.RData")) {
  load("03_results/checkpoints/21_WGCNA_JAG1.RData")
  cat("Loaded WGCNA results\n")
}

if (file.exists("03_results/checkpoints/24_functional_integrated.RData")) {
  load("03_results/checkpoints/24_functional_integrated.RData")
  cat("Loaded functional integration\n")
}

jag1_targets <- target_table[target_table$Confidence_Tier != "Not_Target", ]
cat("\nTotal JAG1 targets:", nrow(jag1_targets), "\n")


cat("\n========================================\n")

# Define gene sets from different methods
gene_sets <- list()

# DE-based (Gold tier)
gene_sets$DE_Gold <- jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Gold"]

# DE-based (Gold + Silver)
gene_sets$DE_GoldSilver <- jag1_targets$GeneID[jag1_targets$Confidence_Tier %in% c("Gold", "Silver")]

# High validation score
gene_sets$High_Validation <- validation_df$GeneID[validation_df$Validation_Level == "High"]

# WGCNA module-associated
if (exists("narrow_modules") && exists("target_modules")) {
  sig_modules <- narrow_modules$Module[narrow_modules$Significant]
  gene_sets$WGCNA_Associated <- names(target_modules)[target_modules %in% as.numeric(sig_modules)]
}

# Print set sizes
cat("Gene set sizes:\n")
for (set_name in names(gene_sets)) {
  cat("  ", set_name, ":", length(gene_sets[[set_name]]), "\n")
}


cat("\n========================================\n")

# Calculate pairwise overlaps
if (length(gene_sets) >= 2) {
  overlap_matrix <- matrix(0, nrow = length(gene_sets), ncol = length(gene_sets))
  rownames(overlap_matrix) <- names(gene_sets)
  colnames(overlap_matrix) <- names(gene_sets)

  for (i in 1:length(gene_sets)) {
    for (j in 1:length(gene_sets)) {
      overlap_matrix[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
    }
  }

  cat("Overlap matrix:\n")
  print(overlap_matrix)

  # Visualize
  dir.create("03_results/figures/27_cross_validation", recursive = TRUE, showWarnings = FALSE)

  png("03_results/figures/27_cross_validation/method_overlap_heatmap.png",
      width = 800, height = 700, res = 120)

  pheatmap(overlap_matrix,
           main = "JAG1 Target Overlap Across Methods",
           color = colorRampPalette(c("white", "steelblue", "darkblue"))(50),
           display_numbers = TRUE,
           number_format = "%d",
           cluster_rows = FALSE,
           cluster_cols = FALSE)

  dev.off()
  cat("\nSaved: method_overlap_heatmap.png\n")
}


cat("\n========================================\n")

# Find genes appearing in multiple gene sets
all_genes <- unique(unlist(gene_sets))
gene_membership <- sapply(all_genes, function(g) {
  sum(sapply(gene_sets, function(s) g %in% s))
})

membership_df <- data.frame(
  GeneID = names(gene_membership),
  N_Methods = as.vector(gene_membership),
  stringsAsFactors = FALSE
)

# Add tier information
membership_df$Confidence_Tier <- jag1_targets$Confidence_Tier[
  match(membership_df$GeneID, jag1_targets$GeneID)
]

# Add validation score
membership_df$Validation_Score <- validation_df$Validation_Score[
  match(membership_df$GeneID, validation_df$GeneID)
]

membership_df <- membership_df[order(-membership_df$N_Methods, -membership_df$Validation_Score), ]

cat("Multi-method support distribution:\n")
print(table(membership_df$N_Methods))

# Consensus genes (in >50% of methods)
n_methods <- length(gene_sets)
consensus_threshold <- ceiling(n_methods / 2)

consensus_genes <- membership_df$GeneID[membership_df$N_Methods >= consensus_threshold]
cat("\nConsensus genes (in >=", consensus_threshold, "methods):", length(consensus_genes), "\n")


cat("\n========================================\n")

# Create final validated target list
# Criteria: Gold tier OR (Silver tier AND high validation) OR consensus

final_validated <- unique(c(
  # All Gold tier
  jag1_targets$GeneID[jag1_targets$Confidence_Tier == "Gold"],
  # Silver tier with high validation
  validation_df$GeneID[validation_df$Confidence_Tier == "Silver" &
                         validation_df$Validation_Level == "High"],
  # Consensus genes
  consensus_genes
))

cat("Final validated targets:", length(final_validated), "\n")

# Create detailed table
final_table <- data.frame(
  GeneID = final_validated,
  stringsAsFactors = FALSE
)

# Add all available information
final_table$Confidence_Tier <- jag1_targets$Confidence_Tier[
  match(final_table$GeneID, jag1_targets$GeneID)
]

final_table$Validation_Score <- validation_df$Validation_Score[
  match(final_table$GeneID, validation_df$GeneID)
]

final_table$N_Methods <- membership_df$N_Methods[
  match(final_table$GeneID, membership_df$GeneID)
]

if (exists("master_table") && "Category" %in% colnames(master_table)) {
  final_table$Category <- master_table$Category[
    match(final_table$GeneID, master_table$GeneID)
  ]
}

# Add logFC
if ("Mean_logFC_Pairwise" %in% colnames(jag1_targets)) {
  final_table$Mean_logFC <- jag1_targets$Mean_logFC_Pairwise[
    match(final_table$GeneID, jag1_targets$GeneID)
  ]
}

# Sort by validation score
final_table <- final_table[order(-final_table$Validation_Score), ]

cat("\nTop 10 final validated targets:\n")
print(head(final_table, 10))

# Save
write.csv(final_table, "03_results/tables/validation/final_validated_targets.csv",
          row.names = FALSE)
cat("\nSaved: final_validated_targets.csv\n")


cat("\n========================================\n")

summary_stats <- data.frame(
  Category = c(
    "Total JAG1 target candidates",
    "Gold tier (DE in 4/4 comparisons)",
    "Silver tier (DE in 3/4 comparisons)",
    "Bronze tier (DE in 2/4 or pooled)",
    "High validation score",
    "Consensus (multiple methods)",
    "Final validated targets"
  ),
  Count = c(
    nrow(jag1_targets),
    sum(jag1_targets$Confidence_Tier == "Gold"),
    sum(jag1_targets$Confidence_Tier == "Silver"),
    sum(jag1_targets$Confidence_Tier == "Bronze"),
    sum(validation_df$Validation_Level == "High"),
    length(consensus_genes),
    length(final_validated)
  ),
  stringsAsFactors = FALSE
)

cat("Summary:\n")
print(summary_stats)

write.csv(summary_stats, "03_results/tables/validation/validation_summary.csv",
          row.names = FALSE)


cat("\n========================================\n")

# Summary figure
png("03_results/figures/27_cross_validation/validation_summary.png",
    width = 1000, height = 600, res = 120)

par(mfrow = c(1, 2))

# Funnel plot
funnel_data <- c(
  nrow(jag1_targets),
  sum(jag1_targets$Confidence_Tier %in% c("Gold", "Silver")),
  sum(validation_df$Validation_Level == "High"),
  length(final_validated)
)
funnel_labels <- c("All Targets", "Gold+Silver", "High Valid.", "Final")

barplot(funnel_data,
        names.arg = funnel_labels,
        col = c("lightgray", "gold", "steelblue", "darkgreen"),
        main = "Target Filtering Funnel",
        ylab = "Number of Genes")

# Final breakdown
tier_in_final <- table(final_table$Confidence_Tier)
pie(tier_in_final,
    main = paste("Final Validated (n=", length(final_validated), ")"),
    col = c("#CD7F32", "#FFD700", "#C0C0C0"))

dev.off()
cat("Saved: validation_summary.png\n")


cat("\n========================================\n")

cross_validation <- list(
  gene_sets = gene_sets,
  membership_df = membership_df,
  consensus_genes = consensus_genes,
  final_validated = final_validated,
  final_table = final_table,
  summary_stats = summary_stats
)

save(
  cross_validation,
  final_table,
  file = "03_results/checkpoints/27_cross_validation.RData"
)

cat("Checkpoint saved: 27_cross_validation.RData\n")


cat("\n================================================================\n")
cat("  Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

cat("KEY RESULTS:\n")
cat("  - Total candidates:", nrow(jag1_targets), "\n")
cat("  - Consensus genes:", length(consensus_genes), "\n")
cat("  - Final validated:", length(final_validated), "\n\n")

