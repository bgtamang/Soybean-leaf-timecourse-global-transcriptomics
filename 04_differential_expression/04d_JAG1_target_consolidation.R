# Script 14: JAG1 Target Consolidation

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 14: JAG1 TARGET CONSOLIDATION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/JAG1_targets", recursive = TRUE, showWarnings = FALSE)
dir.create("03_results/figures/14_JAG1_targets", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "edgeR",
  "limma",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "UpSetR"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINTS =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/12_DE_organized.RData")
cat("Loaded DE results\n")

load("03_results/checkpoints/06_validated.RData")
cat("Loaded expression data\n\n")

# Get expression matrix
expr_mat <- v_primary$E
targets <- targets_primary

# ===== DEFINE ANALYSIS PARAMETERS =====

cat("========================================\n")
cat("SECTION 2: ANALYSIS PARAMETERS\n")
cat("========================================\n\n")

# Key gene IDs
JAG1_ID <- "Glyma.20G116200"
JAG2_ID <- "Glyma.10G273800"

# Define the 4 pairwise comparisons at TP0
# Format: NarrowLine_vs_BroadLine_TP0
pairwise_comparisons <- c(
  "PI612713BvsPI532462A_TP0",  # Pair 1: Narrow1 vs Broad1
  "PI612713BvsLD112170_TP0",   # Pair 2: Narrow1 vs Broad2
  "PI547745vsPI532462A_TP0",   # Pair 3: Narrow2 vs Broad1
  "PI547745vsLD112170_TP0"     # Pair 4: Narrow2 vs Broad2
)

# Pooled comparison
pooled_comparison <- "NarrowvsBroad_TP0"

# Thresholds
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1.0  # 2-fold change

cat("JAG1 ID:", JAG1_ID, "\n")
cat("JAG2 ID:", JAG2_ID, "\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")
cat("logFC threshold:", LOGFC_THRESHOLD, "(", 2^LOGFC_THRESHOLD, "-fold)\n\n")

cat("Pairwise comparisons for target identification:\n")
for (pc in pairwise_comparisons) {
  exists_flag <- ifelse(pc %in% names(de_results), "OK", "MISSING")
  cat("  -", pc, "[", exists_flag, "]\n")
}
cat("\nPooled comparison:", pooled_comparison,
    "[", ifelse(pooled_comparison %in% names(de_results), "OK", "MISSING"), "]\n\n")

# ===== SECTION 3: EXTRACT UP-REGULATED GENES =====

cat("========================================\n")
cat("SECTION 3: EXTRACT UP-REGULATED GENES\n")
cat("========================================\n\n")

cat("Identifying genes UP in Narrow vs Broad (potential JAG1 targets)...\n\n")

# Function to get UP genes from a comparison
get_up_genes <- function(contrast_name, de_results_list, fdr = 0.05, logfc = 1.0) {
  if (!contrast_name %in% names(de_results_list)) {
    warning(paste("Contrast not found:", contrast_name))
    return(character(0))
  }

  results <- de_results_list[[contrast_name]]

  # Genes with FDR < threshold AND logFC > threshold (UP in Narrow)
  up_genes <- results$Gene[results$FDR < fdr & results$logFC > logfc]

  return(up_genes)
}

# Get UP genes for each pairwise comparison
up_genes_list <- list()
for (pc in pairwise_comparisons) {
  up_genes_list[[pc]] <- get_up_genes(pc, de_results, FDR_THRESHOLD, LOGFC_THRESHOLD)
  cat("  ", pc, ": ", length(up_genes_list[[pc]]), " genes UP\n", sep = "")
}

# Get UP genes for pooled comparison
up_genes_pooled <- get_up_genes(pooled_comparison, de_results, FDR_THRESHOLD, LOGFC_THRESHOLD)
cat("\n  Pooled (", pooled_comparison, "): ", length(up_genes_pooled), " genes UP\n\n", sep = "")

# ===== SECTION 4: CREATE CONFIDENCE TIERS =====

cat("========================================\n")
cat("SECTION 4: ASSIGN CONFIDENCE TIERS\n")
cat("========================================\n\n")

# Get all unique UP genes
all_up_genes <- unique(c(unlist(up_genes_list), up_genes_pooled))
cat("Total unique genes UP in any comparison:", length(all_up_genes), "\n\n")

# Create a matrix tracking which comparisons each gene is UP in
up_matrix <- matrix(0, nrow = length(all_up_genes), ncol = length(pairwise_comparisons) + 1)
rownames(up_matrix) <- all_up_genes
colnames(up_matrix) <- c(pairwise_comparisons, pooled_comparison)

# Fill matrix
for (pc in pairwise_comparisons) {
  up_matrix[up_genes_list[[pc]], pc] <- 1
}
up_matrix[up_genes_pooled, pooled_comparison] <- 1

# Count how many pairwise comparisons each gene is UP in
pairwise_count <- rowSums(up_matrix[, pairwise_comparisons, drop = FALSE])
in_pooled <- up_matrix[, pooled_comparison]

# Assign confidence tiers
# Gold: UP in all 4 pairwise comparisons
# Silver: UP in 3/4 pairwise comparisons
# Bronze: UP in pooled OR 2/4 pairwise

tiers <- rep("Not_Target", length(all_up_genes))
names(tiers) <- all_up_genes

# Assign tiers (from highest to lowest priority)
tiers[pairwise_count == 4] <- "Gold"
tiers[pairwise_count == 3 & tiers == "Not_Target"] <- "Silver"
tiers[(pairwise_count == 2 | in_pooled == 1) & tiers == "Not_Target"] <- "Bronze"

# Summary
tier_summary <- table(tiers)
cat("Confidence Tier Assignment:\n")
cat("  Gold (4/4 pairwise):  ", sum(tiers == "Gold"), " genes\n", sep = "")
cat("  Silver (3/4 pairwise):", sum(tiers == "Silver"), " genes\n", sep = "")
cat("  Bronze (2/4 or pooled):", sum(tiers == "Bronze"), " genes\n", sep = "")
cat("  -----------------------------------------\n")
cat("  Total targets:        ", sum(tiers != "Not_Target"), " genes\n\n", sep = "")

# ===== SECTION 5: BUILD TARGET TABLE =====

cat("========================================\n")
cat("SECTION 5: BUILD TARGET TABLE\n")
cat("========================================\n\n")

# Create comprehensive target table
target_table <- data.frame(
  GeneID = all_up_genes,
  Confidence_Tier = tiers,
  N_Pairwise_UP = pairwise_count,
  UP_in_Pooled = in_pooled,
  stringsAsFactors = FALSE
)

# Add logFC and FDR from each comparison
for (pc in c(pairwise_comparisons, pooled_comparison)) {
  if (pc %in% names(de_results)) {
    results <- de_results[[pc]]

    # Match genes
    idx <- match(all_up_genes, results$Gene)

    # Add columns
    target_table[[paste0(pc, "_logFC")]] <- results$logFC[idx]
    target_table[[paste0(pc, "_FDR")]] <- results$FDR[idx]
  }
}

# Calculate mean logFC across pairwise comparisons
logfc_cols <- paste0(pairwise_comparisons, "_logFC")
logfc_cols <- logfc_cols[logfc_cols %in% colnames(target_table)]

if (length(logfc_cols) > 0) {
  target_table$Mean_logFC_Pairwise <- rowMeans(target_table[, logfc_cols], na.rm = TRUE)
}

# Add expression statistics
if (all(all_up_genes %in% rownames(expr_mat))) {
  target_table$Mean_Expression <- rowMeans(expr_mat[all_up_genes, ])
} else {
  # Handle genes not in expression matrix
  expr_means <- rep(NA, length(all_up_genes))
  names(expr_means) <- all_up_genes
  found_genes <- all_up_genes[all_up_genes %in% rownames(expr_mat)]
  expr_means[found_genes] <- rowMeans(expr_mat[found_genes, , drop = FALSE])
  target_table$Mean_Expression <- expr_means
}

# Sort by tier then by mean logFC
tier_order <- c("Gold", "Silver", "Bronze", "Not_Target")
target_table$Tier_Rank <- match(target_table$Confidence_Tier, tier_order)
target_table <- target_table[order(target_table$Tier_Rank, -target_table$Mean_logFC_Pairwise), ]
target_table$Tier_Rank <- NULL

# Add rank within tier
target_table$Rank_in_Tier <- ave(seq_len(nrow(target_table)),
                                  target_table$Confidence_Tier,
                                  FUN = seq_along)

cat("Target table created with", nrow(target_table), "genes\n")
cat("Columns:", ncol(target_table), "\n\n")

# ===== SECTION 6: ADD ANNOTATION =====

cat("========================================\n")
cat("SECTION 6: ADD GENE ANNOTATION\n")
cat("========================================\n\n")

# Try to load annotation file (located in Phase1-Exploratory)
annotation_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (file.exists(annotation_file)) {
  cat("Loading gene annotation from:\n  ", annotation_file, "\n")

  annotation <- read.delim(annotation_file, stringsAsFactors = FALSE)

  # Get relevant columns
  if ("locusName" %in% colnames(annotation)) {
    # Match by gene ID
    idx <- match(target_table$GeneID, annotation$locusName)

    # Add annotation columns if they exist
    anno_cols <- c("Best.hit.arabi.name", "arabi.symbol", "arabi.defline")
    for (col in anno_cols) {
      if (col %in% colnames(annotation)) {
        target_table[[col]] <- annotation[[col]][idx]
      }
    }

    cat("Added annotation for", sum(!is.na(idx)), "genes\n\n")
  }
} else {
  cat("Annotation file not found. Skipping annotation.\n\n")
}

# ===== SECTION 7: SAVE RESULTS =====

cat("========================================\n")
cat("SECTION 7: SAVE RESULTS\n")
cat("========================================\n\n")

# Save full table
write.csv(target_table,
          file = "03_results/tables/JAG1_targets/JAG1_targets_all_tiers.csv",
          row.names = FALSE)
cat("Saved: JAG1_targets_all_tiers.csv\n")

# Save tier-specific tables
for (tier in c("Gold", "Silver", "Bronze")) {
  tier_data <- target_table[target_table$Confidence_Tier == tier, ]
  if (nrow(tier_data) > 0) {
    write.csv(tier_data,
              file = paste0("03_results/tables/JAG1_targets/JAG1_targets_", tolower(tier), ".csv"),
              row.names = FALSE)
    cat("Saved: JAG1_targets_", tolower(tier), ".csv (", nrow(tier_data), " genes)\n", sep = "")
  }
}

# ===== SECTION 8: VISUALIZATIONS =====

cat("\n========================================\n")
cat("SECTION 8: VISUALIZATIONS\n")
cat("========================================\n\n")

# 8.1 UpSet plot of pairwise comparisons
cat("Creating UpSet plot of pairwise comparisons...\n")

upset_data <- as.data.frame(up_matrix[, pairwise_comparisons])
upset_data <- upset_data[rowSums(upset_data) > 0, , drop = FALSE]

# Simplify column names for plotting
colnames(upset_data) <- gsub("vs", " vs ", gsub("_TP0", "", colnames(upset_data)))

if (nrow(upset_data) > 0 && ncol(upset_data) >= 2) {
  tryCatch({
    png("03_results/figures/14_JAG1_targets/upset_JAG1_targets_pairwise.png",
        width = 12, height = 8, units = "in", res = 300)

    print(upset(upset_data,
                sets = colnames(upset_data),
                order.by = "freq",
                main.bar.color = "#D73027",
                sets.bar.color = "#4575B4",
                text.scale = 1.5,
                set_size.show = TRUE))

    dev.off()
    cat("  Saved: upset_JAG1_targets_pairwise.png\n")
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("  Warning: UpSet plot failed -", conditionMessage(e), "\n")
  })
}

# 8.2 Tier distribution barplot
cat("Creating tier distribution plot...\n")

tier_counts <- table(factor(target_table$Confidence_Tier,
                            levels = c("Gold", "Silver", "Bronze")))

png("03_results/figures/14_JAG1_targets/barplot_tier_distribution.png",
    width = 8, height = 6, units = "in", res = 300)

par(mar = c(5, 5, 4, 2))
bp <- barplot(tier_counts,
              col = c("#FFD700", "#C0C0C0", "#CD7F32"),
              main = "JAG1 Target Candidates by Confidence Tier",
              ylab = "Number of Genes",
              xlab = "Confidence Tier",
              cex.names = 1.2,
              cex.lab = 1.2,
              ylim = c(0, max(tier_counts) * 1.2))

# Add counts on bars
text(bp, tier_counts, labels = tier_counts, pos = 3, cex = 1.2, font = 2)

# Add legend explaining tiers
legend("topright",
       legend = c("Gold: UP in 4/4 pairwise",
                  "Silver: UP in 3/4 pairwise",
                  "Bronze: UP in 2/4 or pooled"),
       fill = c("#FFD700", "#C0C0C0", "#CD7F32"),
       cex = 0.9)

dev.off()
cat("  Saved: barplot_tier_distribution.png\n")

# 8.3 Heatmap of Gold tier targets
gold_genes <- target_table$GeneID[target_table$Confidence_Tier == "Gold"]

if (length(gold_genes) > 0 && length(gold_genes) <= 100) {
  cat("Creating heatmap of Gold tier targets...\n")

  # Get TP0 samples only
  tp0_samples <- targets$Sample[targets$Timepoint == "TP0"]
  tp0_samples <- tp0_samples[tp0_samples %in% colnames(expr_mat)]

  if (length(tp0_samples) > 0) {
    # Get expression for gold genes at TP0
    gold_in_expr <- gold_genes[gold_genes %in% rownames(expr_mat)]

    if (length(gold_in_expr) > 0) {
      heat_data <- expr_mat[gold_in_expr, tp0_samples, drop = FALSE]

      # Scale rows
      heat_data_scaled <- t(scale(t(heat_data)))

      # Annotation
      tp0_targets <- targets[targets$Sample %in% tp0_samples, ]
      sample_anno <- data.frame(
        Leaf_Type = tp0_targets$Leaf_type,
        Line = tp0_targets$Line,
        row.names = tp0_targets$Sample
      )

      # Colors
      anno_colors <- list(
        Leaf_Type = c("Broad" = "#2E8B57", "Narrow" = "#D95F02"),
        Line = setNames(brewer.pal(4, "Set1"), unique(tp0_targets$Line))
      )

      png("03_results/figures/14_JAG1_targets/heatmap_gold_targets_TP0.png",
          width = 12, height = max(8, length(gold_in_expr) * 0.15 + 2),
          units = "in", res = 300)

      pheatmap(heat_data_scaled,
               color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
               annotation_col = sample_anno,
               annotation_colors = anno_colors,
               show_colnames = TRUE,
               fontsize_row = 8,
               fontsize_col = 8,
               main = paste0("Gold Tier JAG1 Targets at TP0 (n=", length(gold_in_expr), ")"),
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "ward.D2")

      dev.off()
      cat("  Saved: heatmap_gold_targets_TP0.png\n")
    }
  }
} else if (length(gold_genes) > 100) {
  cat("  Skipping Gold heatmap (", length(gold_genes), " genes - too many)\n")
}

# 8.4 Mean logFC comparison across tiers
cat("Creating logFC distribution by tier...\n")

# Prepare data for ggplot
logfc_data <- target_table[target_table$Confidence_Tier %in% c("Gold", "Silver", "Bronze"),
                            c("GeneID", "Confidence_Tier", "Mean_logFC_Pairwise")]

if (nrow(logfc_data) > 0) {
  logfc_data$Confidence_Tier <- factor(logfc_data$Confidence_Tier,
                                        levels = c("Gold", "Silver", "Bronze"))

  p <- ggplot(logfc_data, aes(x = Confidence_Tier, y = Mean_logFC_Pairwise,
                               fill = Confidence_Tier)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Gold" = "#FFD700", "Silver" = "#C0C0C0", "Bronze" = "#CD7F32")) +
    theme_bw(base_size = 14) +
    labs(title = "Fold Change Distribution by Confidence Tier",
         subtitle = "Mean logFC across pairwise Narrow vs Broad comparisons",
         x = "Confidence Tier",
         y = "Mean log2 Fold Change") +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

  ggsave("03_results/figures/14_JAG1_targets/boxplot_logFC_by_tier.png",
         p, width = 8, height = 6, dpi = 300)
  cat("  Saved: boxplot_logFC_by_tier.png\n")
}

# ===== SECTION 9: JAG1/JAG2 STATUS CHECK =====

cat("\n========================================\n")
cat("SECTION 9: JAG1/JAG2 STATUS CHECK\n")
cat("========================================\n\n")

# Check if JAG1 itself appears as a target (it shouldn't if comparing Narrow vs Broad)
if (JAG1_ID %in% target_table$GeneID) {
  jag1_tier <- target_table$Confidence_Tier[target_table$GeneID == JAG1_ID]
  cat("NOTE: JAG1 itself appears in target list (Tier:", jag1_tier, ")\n")
  cat("  This may indicate expression differences due to genotype.\n")
} else {
  cat("JAG1 (", JAG1_ID, ") is NOT in the target list.\n", sep = "")
  cat("  (Expected - JAG1 is the regulator, not a target)\n")
}

if (JAG2_ID %in% target_table$GeneID) {
  jag2_tier <- target_table$Confidence_Tier[target_table$GeneID == JAG2_ID]
  cat("\nJAG2 (", JAG2_ID, ") IS in target list (Tier: ", jag2_tier, ")\n", sep = "")
  cat("  This could indicate cross-regulation between JAG1 and JAG2.\n")
} else {
  cat("\nJAG2 (", JAG2_ID, ") is NOT in the target list.\n", sep = "")
}

# ===== SAVE CHECKPOINT =====

cat("\n========================================\n")
cat("SAVING CHECKPOINT\n")
cat("========================================\n\n")

# Save key objects
save(target_table,
     up_genes_list,
     up_genes_pooled,
     up_matrix,
     pairwise_comparisons,
     pooled_comparison,
     JAG1_ID, JAG2_ID,
     FDR_THRESHOLD, LOGFC_THRESHOLD,
     file = "03_results/checkpoints/14_JAG1_targets.RData")

cat("Saved checkpoint: 14_JAG1_targets.RData\n")

# ===== SUMMARY =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 14: JAG1 TARGET CONSOLIDATION - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  SUMMARY OF JAG1 TARGET CANDIDATES:\n")
cat("  -----------------------------------------\n")
cat("    Gold tier (4/4 comparisons):  ", sum(target_table$Confidence_Tier == "Gold"), "\n", sep = "")
cat("    Silver tier (3/4 comparisons):", sum(target_table$Confidence_Tier == "Silver"), "\n", sep = "")
cat("    Bronze tier (2/4 or pooled):  ", sum(target_table$Confidence_Tier == "Bronze"), "\n", sep = "")
cat("    -----------------------------------------\n")
cat("    Total candidates:             ", sum(target_table$Confidence_Tier != "Not_Target"), "\n\n", sep = "")
cat("  Files saved:\n")
cat("    - 03_results/tables/JAG1_targets/JAG1_targets_*.csv\n")
cat("    - 03_results/figures/14_JAG1_targets/*.png\n")
cat("    - 03_results/checkpoints/14_JAG1_targets.RData\n")
cat("\n")
cat("  Next: Run Script 15 for JAG1 co-expression analysis\n")
cat("================================================================\n")
