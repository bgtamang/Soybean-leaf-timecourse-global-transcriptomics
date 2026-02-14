# Regenerate Dashboard Data from Gene-Level Expression (Phase3)
# This script updates the dashboard data to use gene-level (not transcript-level) expression
# Run this script from the GmJAG1_Dashboard directory

library(dplyr)

# Set paths
base_dir <- dirname(getwd())  # Phase3-Refined-Analysis
checkpoint_dir <- file.path(base_dir, "03_results", "checkpoints")
data_dir <- "data"

cat("==============================================\n")
cat("Regenerating Dashboard Data (Gene-Level)\n")
cat("==============================================\n\n")

# ============================================================
# 1. Load normalized gene-level expression matrix
# ============================================================
cat("Loading normalized expression data...\n")

norm_file <- file.path(checkpoint_dir, "05_normalized.RData")
if (!file.exists(norm_file)) {
  stop("Normalized data file not found: ", norm_file)
}

load(norm_file)
cat("  Loaded: 05_normalized.RData\n")

# The normalized data should contain voom objects or logCPM matrix
# Check what objects are available
cat("  Objects in checkpoint: ", paste(ls(), collapse = ", "), "\n")

# Extract expression matrix - try various naming conventions
if (exists("v_full")) {
  expr_mat <- v_full$E
  cat("  Using v_full$E (voom-normalized log2-CPM, all samples)\n")
} else if (exists("v_primary")) {
  expr_mat <- v_primary$E
  cat("  Using v_primary$E (voom-normalized log2-CPM)\n")
} else if (exists("v")) {
  expr_mat <- v$E
  cat("  Using v$E (voom-normalized log2-CPM)\n")
} else if (exists("logCPM_corrected")) {
  expr_mat <- logCPM_corrected
  cat("  Using logCPM_corrected (batch-corrected log2-CPM)\n")
} else if (exists("logcpm")) {
  expr_mat <- logcpm
  cat("  Using logcpm matrix\n")
} else {
  # List available objects for debugging
  stop("Could not find expression matrix. Available objects: ", paste(ls(), collapse = ", "))
}

cat("  Expression matrix dimensions: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples\n")
cat("  Sample gene IDs: ", paste(head(rownames(expr_mat), 3), collapse = ", "), "...\n")

# Verify gene-level (no version suffix like .1, .2)
sample_ids <- head(rownames(expr_mat), 10)
has_version <- any(grepl("\\.[0-9]+$", sample_ids))
if (has_version) {
  cat("\n  WARNING: Gene IDs appear to have version suffixes!\n")
  cat("  Example: ", sample_ids[grepl("\\.[0-9]+$", sample_ids)][1], "\n")
  cat("  Stripping version suffixes...\n")

  # Strip version suffix and aggregate by summing (for counts) or averaging (for normalized)
  gene_ids <- sub("\\.[0-9]+$", "", rownames(expr_mat))

  # Average expression for duplicate gene IDs (in case of multiple transcripts)
  expr_df <- as.data.frame(expr_mat)
  expr_df$GeneID <- gene_ids
  expr_mat <- expr_df %>%
    group_by(GeneID) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()
  rownames(expr_mat) <- expr_mat$GeneID
  expr_mat$GeneID <- NULL
  expr_mat <- as.matrix(expr_mat)

  cat("  After aggregation: ", nrow(expr_mat), " genes\n")
} else {
  cat("  Gene IDs are already at gene level (no version suffix)\n")
}

# ============================================================
# 2. Save updated expression matrix
# ============================================================
cat("\nSaving gene-level expression matrix...\n")

saveRDS(expr_mat, file.path(data_dir, "expression_matrix.rds"))
cat("  Saved: data/expression_matrix.rds\n")
cat("  Dimensions: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples\n")

# ============================================================
# 3. Update gene_list.rds
# ============================================================
cat("\nUpdating gene list...\n")

gene_list <- rownames(expr_mat)
saveRDS(gene_list, file.path(data_dir, "gene_list.rds"))
cat("  Saved: data/gene_list.rds (", length(gene_list), " genes)\n")

# ============================================================
# 4. Verify other data files have matching gene IDs
# ============================================================
cat("\nVerifying gene ID consistency...\n")

# Check JAG1 targets
if (file.exists(file.path(data_dir, "jag1_targets.rds"))) {
  jag1 <- readRDS(file.path(data_dir, "jag1_targets.rds"))
  jag1_in_expr <- sum(jag1$GeneID %in% gene_list)
  cat("  JAG1 targets: ", jag1_in_expr, "/", nrow(jag1), " found in expression matrix\n")

  # Check if JAG1 targets have version suffixes
  if (any(grepl("\\.[0-9]+$", jag1$GeneID))) {
    cat("  NOTE: JAG1 targets have version suffixes, stripping...\n")
    jag1$GeneID <- sub("\\.[0-9]+$", "", jag1$GeneID)
    saveRDS(jag1, file.path(data_dir, "jag1_targets.rds"))
    cat("  Updated: data/jag1_targets.rds\n")
  }
}

# Check DE results
if (file.exists(file.path(data_dir, "de_results.rds"))) {
  de <- readRDS(file.path(data_dir, "de_results.rds"))
  if ("GeneID" %in% names(de)) {
    de_in_expr <- sum(de$GeneID %in% gene_list)
    cat("  DE results: ", de_in_expr, "/", length(unique(de$GeneID)), " genes found in expression matrix\n")

    if (any(grepl("\\.[0-9]+$", de$GeneID))) {
      cat("  NOTE: DE results have version suffixes, stripping...\n")
      de$GeneID <- sub("\\.[0-9]+$", "", de$GeneID)
      saveRDS(de, file.path(data_dir, "de_results.rds"))
      cat("  Updated: data/de_results.rds\n")
    }
  }
}

# Check hub genes
if (file.exists(file.path(data_dir, "hub_genes.rds"))) {
  hubs <- readRDS(file.path(data_dir, "hub_genes.rds"))
  if ("GeneID" %in% names(hubs) || "Gene" %in% names(hubs)) {
    gene_col <- if ("GeneID" %in% names(hubs)) "GeneID" else "Gene"
    if (any(grepl("\\.[0-9]+$", hubs[[gene_col]]))) {
      cat("  NOTE: Hub genes have version suffixes, stripping...\n")
      hubs[[gene_col]] <- sub("\\.[0-9]+$", "", hubs[[gene_col]])
      saveRDS(hubs, file.path(data_dir, "hub_genes.rds"))
      cat("  Updated: data/hub_genes.rds\n")
    }
  }
}

# Check gene-phenotype correlations
if (file.exists(file.path(data_dir, "gene_phenotype_cors.rds"))) {
  cors <- readRDS(file.path(data_dir, "gene_phenotype_cors.rds"))
  if ("GeneID" %in% names(cors)) {
    if (any(grepl("\\.[0-9]+$", cors$GeneID))) {
      cat("  NOTE: Gene-phenotype correlations have version suffixes, stripping...\n")
      cors$GeneID <- sub("\\.[0-9]+$", "", cors$GeneID)
      saveRDS(cors, file.path(data_dir, "gene_phenotype_cors.rds"))
      cat("  Updated: data/gene_phenotype_cors.rds\n")
    }
  }
}

# ============================================================
# 5. Copy housekeeping genes to data folder
# ============================================================
cat("\nCopying housekeeping genes file...\n")

hk_source <- file.path(base_dir, "03_results", "tables", "housekeeping_stability.csv")
hk_dest <- file.path(data_dir, "housekeeping_stability.csv")

if (file.exists(hk_source)) {
  file.copy(hk_source, hk_dest, overwrite = TRUE)
  cat("  Copied: housekeeping_stability.csv to data folder\n")

  # Verify housekeeping genes match expression matrix
  hk <- read.csv(hk_dest, stringsAsFactors = FALSE)
  hk_in_expr <- sum(hk$Gene %in% gene_list)
  cat("  Housekeeping genes: ", hk_in_expr, "/", nrow(hk), " found in expression matrix\n")
} else {
  cat("  WARNING: housekeeping_stability.csv not found at: ", hk_source, "\n")
}

# ============================================================
# 6. Summary
# ============================================================
cat("\n==============================================\n")
cat("Dashboard Data Regeneration Complete!\n")
cat("==============================================\n")
cat("\nFinal expression matrix:\n")
cat("  - Genes: ", nrow(expr_mat), "\n")
cat("  - Samples: ", ncol(expr_mat), "\n")
cat("  - Gene ID format: ", rownames(expr_mat)[1], " (no version suffix)\n")
cat("\nYou can now restart the Shiny app.\n")
