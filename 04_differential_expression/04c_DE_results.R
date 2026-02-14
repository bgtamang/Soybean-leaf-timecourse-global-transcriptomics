# Script 12: Differential Expression Results

# ===== CLEAR ENVIRONMENT =====
rm(list = ls())
gc()

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 12: DE RESULTS EXTRACTION\n")
cat("  GmJAG1 Soybean RNA-Seq Analysis\n")
cat("================================================================\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

# ===== SETUP =====

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))
cat("Working directory:", getwd(), "\n\n")

# Create output directories
dir.create("03_results/tables/DE", recursive = TRUE, showWarnings = FALSE)

# ===== LOAD REQUIRED PACKAGES =====

cat("Loading required packages...\n")

required_packages <- c(
  "dplyr",
  "tidyr"
)

invisible(lapply(required_packages, library, character.only = TRUE))
cat("  Packages loaded\n\n")

# ===== LOAD CHECKPOINT =====

cat("========================================\n")
cat("SECTION 1: LOAD DATA\n")
cat("========================================\n\n")

load("03_results/checkpoints/11_DE_results.RData")
cat("Loaded checkpoint from script 11\n")
cat("  DE results for", length(de_results), "contrasts\n\n")

# ===== SECTION 2: EXPORT INDIVIDUAL CONTRAST RESULTS =====

cat("========================================\n")
cat("SECTION 2: EXPORT INDIVIDUAL RESULTS\n")
cat("========================================\n\n")

cat("Saving results for each contrast...\n")

# Significance thresholds
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1

for (contrast_name in names(de_results)) {
  results <- de_results[[contrast_name]]

  # Add significance column
  results$Significant <- results$FDR < FDR_THRESHOLD & abs(results$logFC) > LOGFC_THRESHOLD
  results$Direction <- ifelse(results$logFC > 0, "Up", "Down")
  results$Direction[!results$Significant] <- "NS"

  # Reorder columns
  results <- results[, c("Gene", "logFC", "logCPM", "F", "PValue", "FDR",
                         "Significant", "Direction")]

  # Sort by FDR
  results <- results[order(results$FDR), ]

  # Save full results
  filename <- paste0("03_results/tables/DE/", contrast_name, ".csv")
  write.csv(results, filename, row.names = FALSE)

  # Update de_results with enhanced version
  de_results[[contrast_name]] <- results
}

cat("  Saved", length(de_results), "individual result files\n\n")

# ===== SECTION 2B: CREATE COMPREHENSIVE RESULTS MATRIX =====

cat("========================================\n")
cat("SECTION 2B: COMPREHENSIVE RESULTS MATRIX\n")
cat("========================================\n\n")

cat("Creating wide-format matrix with all contrasts...\n")

# Get all genes and contrast names
all_genes <- rownames(de_results[[1]])
contrast_names <- names(de_results)

# Initialize matrices for different statistics
logFC_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
pvalue_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
fdr_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))
logCPM_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(contrast_names))

# Set row and column names
rownames(logFC_matrix) <- all_genes
rownames(pvalue_matrix) <- all_genes
rownames(fdr_matrix) <- all_genes
rownames(logCPM_matrix) <- all_genes

colnames(logFC_matrix) <- paste0(contrast_names, "_logFC")
colnames(pvalue_matrix) <- paste0(contrast_names, "_PValue")
colnames(fdr_matrix) <- paste0(contrast_names, "_FDR")
colnames(logCPM_matrix) <- paste0(contrast_names, "_logCPM")

# Fill matrices with results
for (i in seq_along(contrast_names)) {
  contrast_name <- contrast_names[i]
  results <- de_results[[contrast_name]]

  # Match genes
  gene_idx <- match(results$Gene, all_genes)

  logFC_matrix[gene_idx, i] <- results$logFC
  pvalue_matrix[gene_idx, i] <- results$PValue
  fdr_matrix[gene_idx, i] <- results$FDR
  logCPM_matrix[gene_idx, i] <- results$logCPM
}

# Create comprehensive data frame
comprehensive_de_results <- data.frame(
  GeneID = all_genes,
  stringsAsFactors = FALSE
)

# Add summary statistics first
comprehensive_de_results$Max_AbsLogFC <- apply(logFC_matrix, 1, function(x) {
  if(all(is.na(x))) return(NA)
  max(abs(x), na.rm = TRUE)
})

comprehensive_de_results$Min_FDR <- apply(fdr_matrix, 1, function(x) {
  if(all(is.na(x))) return(NA)
  min(x, na.rm = TRUE)
})

comprehensive_de_results$N_Sig_Contrasts <- rowSums(
  fdr_matrix < FDR_THRESHOLD & abs(logFC_matrix) > LOGFC_THRESHOLD,
  na.rm = TRUE
)

comprehensive_de_results$AveExpr <- rowMeans(logCPM_matrix, na.rm = TRUE)

cat("  Matrix created:", nrow(comprehensive_de_results), "genes x", length(contrast_names), "contrasts\n")

# ===== SECTION 2C: ADD GENE ANNOTATION =====

cat("\n========================================\n")
cat("SECTION 2C: ADD GENE ANNOTATION\n")
cat("========================================\n\n")

cat("Loading gene annotation...\n")

# Path to annotation file (located in Phase1-Exploratory)
annotation_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if (file.exists(annotation_file)) {
  annotation <- read.delim(annotation_file, stringsAsFactors = FALSE)

  cat("  Loaded annotation with", nrow(annotation), "entries\n")
  cat("  Columns available:", paste(colnames(annotation), collapse = ", "), "\n\n")

  # Get unique gene-level annotation (use locusName)
  if ("locusName" %in% colnames(annotation)) {
    # Select useful columns
    ann_cols <- c("locusName")
    if ("Best.hit.arabi.name" %in% colnames(annotation)) ann_cols <- c(ann_cols, "Best.hit.arabi.name")
    if ("Best-hit-arabi-name" %in% colnames(annotation)) ann_cols <- c(ann_cols, "Best-hit-arabi-name")
    if ("Best.hit.arabi.defline" %in% colnames(annotation)) ann_cols <- c(ann_cols, "Best.hit.arabi.defline")
    if ("Best-hit-arabi-defline" %in% colnames(annotation)) ann_cols <- c(ann_cols, "Best-hit-arabi-defline")
    if ("arabi.symbol" %in% colnames(annotation)) ann_cols <- c(ann_cols, "arabi.symbol")
    if ("arabi-symbol" %in% colnames(annotation)) ann_cols <- c(ann_cols, "arabi-symbol")

    # Get existing columns only
    ann_cols <- ann_cols[ann_cols %in% colnames(annotation)]

    gene_annotation <- annotation[, ann_cols, drop = FALSE]

    # Remove duplicate genes (keep first occurrence)
    gene_annotation <- gene_annotation[!duplicated(gene_annotation$locusName), ]

    # Match annotation to our genes
    annotation_idx <- match(comprehensive_de_results$GeneID, gene_annotation$locusName)

    # Add annotation columns
    for (col in setdiff(ann_cols, "locusName")) {
      comprehensive_de_results[[col]] <- gene_annotation[[col]][annotation_idx]
    }

    # Report coverage
    n_annotated <- sum(!is.na(annotation_idx))
    cat("  Genes with annotation:", n_annotated, "out of", nrow(comprehensive_de_results),
        "(", round(n_annotated/nrow(comprehensive_de_results)*100, 1), "%)\n")
  }
} else {
  cat("  Warning: Annotation file not found at:", annotation_file, "\n")
  cat("  Proceeding without annotation\n")
}

# Add all statistics matrices
comprehensive_de_results <- cbind(
  comprehensive_de_results,
  logFC_matrix,
  fdr_matrix
)

# Create decision matrix
decision_matrix <- matrix("NS", nrow = length(all_genes), ncol = length(contrast_names))
rownames(decision_matrix) <- all_genes
colnames(decision_matrix) <- paste0(contrast_names, "_Direction")

for (i in seq_along(contrast_names)) {
  fdr_col <- fdr_matrix[, i]
  logfc_col <- logFC_matrix[, i]

  decision_matrix[, i] <- ifelse(
    is.na(fdr_col) | is.na(logfc_col), NA,
    ifelse(fdr_col < FDR_THRESHOLD & logfc_col > LOGFC_THRESHOLD, "Up",
           ifelse(fdr_col < FDR_THRESHOLD & logfc_col < -LOGFC_THRESHOLD, "Down", "NS"))
  )
}

comprehensive_de_results <- cbind(comprehensive_de_results, decision_matrix)

# Sort by number of significant contrasts and minimum FDR
comprehensive_de_results <- comprehensive_de_results[
  order(-comprehensive_de_results$N_Sig_Contrasts,
        comprehensive_de_results$Min_FDR),
]

# Save comprehensive results - all genes
write.csv(comprehensive_de_results,
          "03_results/tables/DE/comprehensive_DE_results_all_contrasts.csv",
          row.names = FALSE)
cat("\nSaved: 03_results/tables/DE/comprehensive_DE_results_all_contrasts.csv\n")

# Save significant genes only
significant_genes_df <- comprehensive_de_results[comprehensive_de_results$N_Sig_Contrasts > 0, ]
write.csv(significant_genes_df,
          "03_results/tables/DE/comprehensive_DE_results_significant_only.csv",
          row.names = FALSE)
cat("Saved: 03_results/tables/DE/comprehensive_DE_results_significant_only.csv\n")
cat("  Significant genes:", nrow(significant_genes_df), "\n")

# Save top 100 DE genes
top_de_genes <- head(comprehensive_de_results[comprehensive_de_results$N_Sig_Contrasts > 0, ], 100)
write.csv(top_de_genes,
          "03_results/tables/DE/top100_DE_genes_annotated.csv",
          row.names = FALSE)
cat("Saved: 03_results/tables/DE/top100_DE_genes_annotated.csv\n\n")

# ===== SECTION 3: CREATE COMBINED RESULTS TABLE =====

cat("========================================\n")
cat("SECTION 3: COMBINED RESULTS\n")
cat("========================================\n\n")

cat("Creating combined DE results table...\n")

# Combine all significant results
all_sig_results <- data.frame()

for (contrast_name in names(de_results)) {
  results <- de_results[[contrast_name]]
  sig_results <- results[results$Significant, ]

  if (nrow(sig_results) > 0) {
    sig_results$Contrast <- contrast_name

    # Add contrast type
    type_info <- contrast_desc[contrast_desc$Contrast == contrast_name, "Type"]
    sig_results$Type <- if (length(type_info) > 0) type_info else "Unknown"

    all_sig_results <- rbind(all_sig_results, sig_results)
  }
}

cat("  Total significant results:", nrow(all_sig_results), "\n")
cat("  Unique genes:", length(unique(all_sig_results$Gene)), "\n\n")

# Save combined results
write.csv(all_sig_results, "03_results/tables/DE_all_significant.csv", row.names = FALSE)
cat("Saved: 03_results/tables/DE_all_significant.csv\n\n")

# ===== SECTION 4: MULTI-CONTRAST ANALYSIS =====

cat("========================================\n")
cat("SECTION 4: MULTI-CONTRAST ANALYSIS\n")
cat("========================================\n\n")

# Count how many contrasts each gene is DE in
gene_contrast_count <- all_sig_results %>%
  group_by(Gene) %>%
  summarize(
    N_contrasts = n(),
    Contrasts = paste(Contrast, collapse = "; "),
    Directions = paste(Direction, collapse = "; "),
    Mean_logFC = mean(logFC),
    Min_FDR = min(FDR),
    .groups = "drop"
  ) %>%
  arrange(desc(N_contrasts), Min_FDR)

cat("Genes by number of DE contrasts:\n")
print(table(gene_contrast_count$N_contrasts))

# Save multi-contrast results
write.csv(gene_contrast_count, "03_results/tables/DE_genes_multicontrast.csv", row.names = FALSE)
cat("\nSaved: 03_results/tables/DE_genes_multicontrast.csv\n")

# Top genes (DE in most contrasts)
cat("\n\nTop 20 genes (DE in most contrasts):\n")
print(head(gene_contrast_count[, c("Gene", "N_contrasts", "Mean_logFC", "Min_FDR")], 20))

# ===== SECTION 5: KEY CONTRAST ANALYSIS (JAG1 TARGETS) =====

cat("\n========================================\n")
cat("SECTION 5: KEY CONTRAST ANALYSIS\n")
cat("========================================\n\n")

# Analyze Narrow vs Broad at TP0 (genes up in narrow are potential JAG1 targets)
cat("Analyzing Narrow vs Broad comparisons at TP0...\n")
cat("  (Genes UP in Narrow = potential JAG1 targets - derepressed)\n\n")

# Get genes UP in narrow for each comparison
narrow_up_genes <- list()

for (contrast_name in key_contrasts) {
  if (contrast_name %in% names(de_results)) {
    results <- de_results[[contrast_name]]
    up_genes <- results$Gene[results$Significant & results$Direction == "Up"]
    narrow_up_genes[[contrast_name]] <- up_genes
    cat("  ", contrast_name, ":", length(up_genes), "genes UP in Narrow\n")
  }
}

# Find genes consistently UP in narrow across all comparisons
if (length(narrow_up_genes) > 0) {
  all_narrow_up <- Reduce(intersect, narrow_up_genes)
  cat("\n  Genes UP in Narrow in ALL", length(narrow_up_genes), "comparisons:",
      length(all_narrow_up), "\n")

  # Find genes UP in majority (>= 3/4)
  gene_up_count <- table(unlist(narrow_up_genes))
  majority_up <- names(gene_up_count[gene_up_count >= 3])
  cat("  Genes UP in Narrow in >= 3/4 comparisons:", length(majority_up), "\n")

  # Create JAG1 candidate table
  jag1_candidates <- data.frame(
    Gene = names(gene_up_count),
    N_comparisons_up = as.numeric(gene_up_count)
  ) %>%
    arrange(desc(N_comparisons_up))

  # Add expression info from pooled comparison if available
  if ("NarrowvsBroad_TP0" %in% names(de_results)) {
    pooled_results <- de_results[["NarrowvsBroad_TP0"]]
    jag1_candidates <- merge(jag1_candidates, pooled_results[, c("Gene", "logFC", "FDR")],
                             by = "Gene", all.x = TRUE)
    colnames(jag1_candidates)[colnames(jag1_candidates) == "logFC"] <- "Pooled_logFC"
    colnames(jag1_candidates)[colnames(jag1_candidates) == "FDR"] <- "Pooled_FDR"
    jag1_candidates <- jag1_candidates %>% arrange(desc(N_comparisons_up), Pooled_FDR)
  }

  # Save JAG1 candidates
  write.csv(jag1_candidates, "03_results/tables/JAG1_candidates_preliminary.csv", row.names = FALSE)
  cat("\nSaved: 03_results/tables/JAG1_candidates_preliminary.csv\n")

  cat("\nTop 20 preliminary JAG1 target candidates:\n")
  print(head(jag1_candidates, 20))
}

# ===== SECTION 6: TEMPORAL DE SUMMARY =====

cat("\n========================================\n")
cat("SECTION 6: TEMPORAL DE SUMMARY\n")
cat("========================================\n\n")

# Analyze temporal patterns
temporal_contrasts <- de_summary[de_summary$Type == "Temporal_vsTP0", ]

if (nrow(temporal_contrasts) > 0) {
  temporal_contrasts <- temporal_contrasts %>%
    mutate(
      Line = gsub("_TP[0-9]vsTP0", "", Contrast),
      Timepoint = gsub(".*_(TP[0-9])vsTP0", "\\1", Contrast)
    )

  cat("DE genes per line and timepoint:\n\n")

  temporal_wide <- temporal_contrasts %>%
    select(Line, Timepoint, Total_DE, Up, Down) %>%
    arrange(Line, Timepoint)

  print(as.data.frame(temporal_wide))

  # Save temporal summary
  write.csv(temporal_wide, "03_results/tables/DE_temporal_summary.csv", row.names = FALSE)
  cat("\nSaved: 03_results/tables/DE_temporal_summary.csv\n")
}

# ===== SECTION 7: SAVE OUTPUTS =====

cat("\n========================================\n")
cat("SECTION 7: SAVE OUTPUTS\n")
cat("========================================\n\n")

# Save organized checkpoint
save(de_results, de_results_treat, de_summary, contrast_desc,
     edgeR_coded, edgeR_tr_coded,  # TestResults objects from script 11
     comprehensive_de_results, significant_genes_df,
     all_sig_results, gene_contrast_count,
     logFC_matrix, fdr_matrix,  # Matrices for downstream use
     key_contrasts, broad_lines, narrow_lines,
     FDR_THRESHOLD, LOGFC_THRESHOLD, LFC_THRESHOLD,
     file = "03_results/checkpoints/12_DE_organized.RData")
cat("Saved: 03_results/checkpoints/12_DE_organized.RData\n")
cat("  Includes: comprehensive_de_results, TestResults objects, matrices\n")

# ===== SESSION INFO =====

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

print(sessionInfo())

# ===== COMPLETION =====

cat("\n")
cat("================================================================\n")
cat("  SCRIPT 12: DE RESULTS - COMPLETE\n")
cat("================================================================\n")
cat("  Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")
cat("  Summary:\n")
cat("    - Individual contrast files:", length(de_results), "\n")
cat("    - Total significant DE events:", nrow(all_sig_results), "\n")
cat("    - Unique DE genes:", length(unique(all_sig_results$Gene)), "\n")
cat("    - Genes DE in >= 5 contrasts:", sum(gene_contrast_count$N_contrasts >= 5), "\n")
cat("\n")
cat("  Output files:\n")
cat("    - 03_results/tables/DE/*.csv (individual contrasts)\n")
cat("    - 03_results/tables/DE/comprehensive_DE_results_all_contrasts.csv\n")
cat("    - 03_results/tables/DE/comprehensive_DE_results_significant_only.csv\n")
cat("    - 03_results/tables/DE/top100_DE_genes_annotated.csv\n")
cat("    - 03_results/tables/DE_all_significant.csv\n")
cat("    - 03_results/tables/DE_genes_multicontrast.csv\n")
cat("    - 03_results/tables/JAG1_candidates_preliminary.csv\n")
cat("    - 03_results/tables/DE_temporal_summary.csv\n")
cat("\n")
cat("  Next: Run 13_DE_visualization.R\n")
cat("================================================================\n")
