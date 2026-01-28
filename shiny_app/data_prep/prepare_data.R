#!/usr/bin/env Rscript
# Preprocess data for Shiny dashboard

library(dplyr)
library(readr)
library(tidyr)

# Paths
base_path <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase2-Refined-Analysis"
tables_path <- file.path(base_path, "03_results/tables")
data_path <- file.path(base_path, "01_data")
output_path <- file.path(base_path, "GmJAG1_Dashboard/data")

cat("=== Preparing Shiny Dashboard Data ===\n\n")

# -----------------------------------------------------------------------------
# 1. Experimental Design (Sample Metadata)
# -----------------------------------------------------------------------------
cat("1. Processing experimental design...\n")
exp_design <- read_csv(file.path(tables_path, "experimental_design.csv"), show_col_types = FALSE)
saveRDS(exp_design, file.path(output_path, "experimental_design.rds"), compress = "xz")
cat("   Saved:", nrow(exp_design), "samples\n")

# -----------------------------------------------------------------------------
# 2. Expression Matrix (from Salmon output)
# -----------------------------------------------------------------------------
cat("\n2. Processing expression data...\n")
load(file.path(data_path, "SalmonSummarizedOutput.RData"))

# tx.all contains counts and abundance (TPM)
if (exists("tx.all")) {
  # TPM matrix for expression viewer
  tpm_matrix <- tx.all$abundance

  # Clean column names - strip "salmo" prefix to match experimental design
  colnames(tpm_matrix) <- gsub("^salmo", "", colnames(tpm_matrix))

  # Keep top 20,000 most variable genes for performance
  gene_var <- apply(tpm_matrix, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(20000, length(gene_var))]
  tpm_subset <- tpm_matrix[top_genes, ]

  saveRDS(tpm_subset, file.path(output_path, "expression_matrix.rds"), compress = "xz")
  cat("   Saved:", nrow(tpm_subset), "genes x", ncol(tpm_subset), "samples\n")
  cat("   Sample names (first 3):", paste(head(colnames(tpm_subset), 3), collapse = ", "), "\n")

  # Full gene list for lookup
  gene_list <- data.frame(GeneID = rownames(tpm_matrix))
  saveRDS(gene_list, file.path(output_path, "gene_list.rds"), compress = "xz")
  cat("   Full gene list:", nrow(gene_list), "genes\n")
} else {
  cat("   WARNING: tx.all not found in RData file\n")
}

# -----------------------------------------------------------------------------
# 3. JAG1 Targets
# -----------------------------------------------------------------------------
cat("\n3. Processing JAG1 targets...\n")
jag1_targets <- read_csv(file.path(tables_path, "JAG1_targets/JAG1_targets_FINAL.csv"), show_col_types = FALSE)

# Select key columns for dashboard
jag1_slim <- jag1_targets %>%
  select(
    Rank, GeneID, Confidence_Tier, Final_Confidence, Confidence_Score,
    N_Pairwise_UP, Mean_logFC_Pairwise, Pattern,
    Best_hit_arabi_name, Best_hit_arabi_defline,
    JAG1_Correlation_TP0, Mean_Expression,
    starts_with("FC_TP"),
    NarrowvsBroad_TP0_logFC, NarrowvsBroad_TP0_FDR,
    UP_in_Pooled, Correlated_with_JAG1
  )

saveRDS(jag1_slim, file.path(output_path, "jag1_targets.rds"), compress = "xz")
cat("   Saved:", nrow(jag1_slim), "targets with", ncol(jag1_slim), "columns\n")

# Tier summary
tier_summary <- jag1_targets %>%
  group_by(Confidence_Tier) %>%
  summarise(Count = n(), .groups = "drop")
saveRDS(tier_summary, file.path(output_path, "tier_summary.rds"))
cat("   Tier summary:", paste(tier_summary$Confidence_Tier, tier_summary$Count, sep = ":", collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 4. DE Results Summary
# -----------------------------------------------------------------------------
cat("\n4. Processing DE results...\n")
de_file <- file.path(tables_path, "DE/comprehensive_DE_results_significant_only.csv")

if (file.exists(de_file)) {
  de_full <- read_csv(de_file, show_col_types = FALSE)

  # Columns are named like: PI532462A_TP1vsTP0_logFC and PI532462A_TP1vsTP0_FDR
  logfc_cols <- grep("_logFC$", names(de_full), value = TRUE)
  fdr_cols <- grep("_FDR$", names(de_full), value = TRUE)

  if (length(logfc_cols) > 0) {
    # Get annotation columns (may be named differently)
    anno_cols <- intersect(c("Best.hit.arabi.name", "Best.hit.arabi.defline",
                             "Best_hit_arabi_name", "Best_hit_arabi_defline"), names(de_full))

    # Pivot logFC columns
    de_logfc <- de_full %>%
      select(GeneID, all_of(anno_cols), all_of(logfc_cols)) %>%
      pivot_longer(
        cols = all_of(logfc_cols),
        names_to = "Contrast",
        values_to = "logFC"
      ) %>%
      mutate(Contrast = gsub("_logFC$", "", Contrast))

    # Pivot FDR columns
    de_fdr <- de_full %>%
      select(GeneID, all_of(fdr_cols)) %>%
      pivot_longer(
        cols = all_of(fdr_cols),
        names_to = "Contrast",
        values_to = "FDR"
      ) %>%
      mutate(Contrast = gsub("_FDR$", "", Contrast))

    # Join and filter
    de_long <- de_logfc %>%
      left_join(de_fdr, by = c("GeneID", "Contrast")) %>%
      filter(!is.na(logFC), !is.na(FDR))

    # Rename annotation columns for consistency
    if ("Best.hit.arabi.name" %in% names(de_long)) {
      de_long <- de_long %>% rename(Best_hit_arabi_name = Best.hit.arabi.name)
    }
    if ("Best.hit.arabi.defline" %in% names(de_long)) {
      de_long <- de_long %>% rename(Best_hit_arabi_defline = Best.hit.arabi.defline)
    }

    saveRDS(de_long, file.path(output_path, "de_results.rds"), compress = "xz")
    cat("   Saved:", nrow(de_long), "gene-contrast pairs\n")

    # Contrast list
    contrasts <- unique(de_long$Contrast)
    saveRDS(contrasts, file.path(output_path, "contrast_list.rds"))
    cat("   Contrasts:", length(contrasts), "\n")
  } else {
    cat("   WARNING: No logFC columns found in DE results\n")
  }
} else {
  cat("   WARNING: DE results file not found\n")
}

# -----------------------------------------------------------------------------
# 5. WGCNA Module Data
# -----------------------------------------------------------------------------
cat("\n5. Processing WGCNA data...\n")

# Module sizes first (to create mapping)
module_sizes_raw <- read_csv(file.path(tables_path, "WGCNA/module_sizes.csv"), show_col_types = FALSE)

# Create module number to color mapping
module_color_map <- module_sizes_raw %>%
  mutate(Module_Num = Module) %>%
  select(Module_Num, Color, Size)

# Module-trait correlations - convert MExx to color names
module_traits_raw <- read_csv(file.path(tables_path, "WGCNA/module_trait_correlations.csv"), show_col_types = FALSE)

# Extract module number from MExx format
module_traits <- module_traits_raw %>%
  mutate(Module_Num = as.numeric(gsub("ME", "", Module))) %>%
  left_join(module_color_map %>% select(Module_Num, Color), by = "Module_Num") %>%
  mutate(Module = ifelse(!is.na(Color), Color, Module)) %>%
  select(-Module_Num, -Color)

saveRDS(module_traits, file.path(output_path, "module_traits.rds"), compress = "xz")
cat("   Module-trait:", nrow(module_traits), "modules\n")

# Module sizes - use Color as Module name for consistency
module_sizes <- module_sizes_raw %>%
  mutate(Module = Color) %>%
  select(Module, Size, Color, Percentage)

saveRDS(module_sizes, file.path(output_path, "module_sizes.rds"))
cat("   Module sizes:", nrow(module_sizes), "modules\n")

# Hub genes (if exists)
hub_file <- file.path(tables_path, "WGCNA/hub_genes.csv")
if (file.exists(hub_file)) {
  hub_genes_raw <- read_csv(hub_file, show_col_types = FALSE)

  # Convert module to color if needed
  if ("Module" %in% names(hub_genes_raw) && is.numeric(hub_genes_raw$Module)) {
    hub_genes_raw <- hub_genes_raw %>%
      left_join(module_color_map %>% select(Module_Num, Color), by = c("Module" = "Module_Num")) %>%
      mutate(Module = ifelse(!is.na(Color), Color, as.character(Module))) %>%
      select(-Color)
  }

  # Keep top 20 per module
  hub_top <- hub_genes_raw %>%
    group_by(Module) %>%
    slice_head(n = 20) %>%
    ungroup()
  saveRDS(hub_top, file.path(output_path, "hub_genes.rds"), compress = "xz")
  cat("   Hub genes: top 20 per module\n")
}

# JAG1 module enrichment - convert to color names
jag1_enrich_raw <- read_csv(file.path(tables_path, "WGCNA/JAG1_module_enrichment.csv"), show_col_types = FALSE)

# Use Module_Color if available, otherwise map from Module number
if ("Module_Color" %in% names(jag1_enrich_raw)) {
  jag1_enrich <- jag1_enrich_raw %>%
    mutate(Module = Module_Color)
} else if (is.numeric(jag1_enrich_raw$Module)) {
  jag1_enrich <- jag1_enrich_raw %>%
    left_join(module_color_map %>% select(Module_Num, Color), by = c("Module" = "Module_Num")) %>%
    mutate(Module = ifelse(!is.na(Color), Color, as.character(Module))) %>%
    select(-Color)
} else {
  jag1_enrich <- jag1_enrich_raw
}

saveRDS(jag1_enrich, file.path(output_path, "jag1_module_enrichment.rds"))
cat("   JAG1 enrichment:", nrow(jag1_enrich), "modules\n")

# -----------------------------------------------------------------------------
# 6. GO Enrichment
# -----------------------------------------------------------------------------
cat("\n6. Processing GO enrichment...\n")

# Use combined module GO file for more comprehensive data
go_combined_file <- file.path(tables_path, "functional/modules/GO_all_modules_combined.csv")
go_targets_file <- file.path(tables_path, "functional/GO_all_targets.csv")

if (file.exists(go_combined_file)) {
  go_results <- read_csv(go_combined_file, show_col_types = FALSE)
  # Add Description column using GO_ID if not present
  if (!"Description" %in% names(go_results) && "GO_ID" %in% names(go_results)) {
    go_results$Description <- go_results$GO_ID
  }
  # Add Term column as alias
  if (!"Term" %in% names(go_results)) {
    go_results$Term <- go_results$GO_ID
  }
  saveRDS(go_results, file.path(output_path, "go_enrichment.rds"), compress = "xz")
  cat("   GO terms (combined modules):", nrow(go_results), "\n")
} else if (file.exists(go_targets_file)) {
  go_results <- read_csv(go_targets_file, show_col_types = FALSE)
  if (!"Description" %in% names(go_results) && "GO_ID" %in% names(go_results)) {
    go_results$Description <- go_results$GO_ID
  }
  saveRDS(go_results, file.path(output_path, "go_enrichment.rds"), compress = "xz")
  cat("   GO terms:", nrow(go_results), "\n")
}

# GO by tier
for (tier in c("gold", "silver", "bronze")) {
  tier_file <- file.path(tables_path, paste0("functional/GO_", tier, "_tier.csv"))
  if (file.exists(tier_file)) {
    tier_go <- read_csv(tier_file, show_col_types = FALSE)
    if (!"Description" %in% names(tier_go) && "GO_ID" %in% names(tier_go)) {
      tier_go$Description <- tier_go$GO_ID
    }
    saveRDS(tier_go, file.path(output_path, paste0("go_", tier, ".rds")))
  }
}

# -----------------------------------------------------------------------------
# 7. Phenotype Data
# -----------------------------------------------------------------------------
cat("\n7. Processing phenotype data...\n")

# Phenotype traits
pheno_file <- file.path(tables_path, "phenotype/phenotype_traits.csv")
if (file.exists(pheno_file)) {
  pheno_traits <- read_csv(pheno_file, show_col_types = FALSE)
  saveRDS(pheno_traits, file.path(output_path, "phenotype_traits.rds"))
  cat("   Phenotype traits loaded\n")
}

# Gene-phenotype correlations (top 5000 for performance)
gene_pheno_file <- file.path(tables_path, "phenotype/gene_phenotype_correlations.csv")
if (file.exists(gene_pheno_file)) {
  gene_pheno <- read_csv(gene_pheno_file, show_col_types = FALSE)
  # Keep top correlated genes (use Cor_LW_Ratio column)
  gene_pheno_top <- gene_pheno %>%
    arrange(desc(abs(Cor_LW_Ratio))) %>%
    head(5000)
  saveRDS(gene_pheno_top, file.path(output_path, "gene_phenotype_cors.rds"), compress = "xz")
  cat("   Gene-phenotype: top 5000 correlations\n")
}

# -----------------------------------------------------------------------------
# 8. QC Metrics
# -----------------------------------------------------------------------------
cat("\n8. Processing QC data...\n")
qc_file <- file.path(tables_path, "QC_metrics.csv")
if (file.exists(qc_file)) {
  qc_metrics <- read_csv(qc_file, show_col_types = FALSE)
  saveRDS(qc_metrics, file.path(output_path, "qc_metrics.rds"))
  cat("   QC metrics:", nrow(qc_metrics), "samples\n")
}

# PCA data from checkpoint
pca_checkpoint <- file.path(base_path, "03_results/checkpoints/07_PCA_complete.RData")
if (file.exists(pca_checkpoint)) {
  load(pca_checkpoint)
  if (exists("pca_results")) {
    pca_data <- data.frame(
      Sample = rownames(pca_results$x),
      PC1 = pca_results$x[, 1],
      PC2 = pca_results$x[, 2],
      PC3 = pca_results$x[, 3]
    )
    pca_var <- summary(pca_results)$importance[2, 1:3] * 100

    saveRDS(list(coords = pca_data, variance = pca_var), file.path(output_path, "pca_data.rds"))
    cat("   PCA data saved\n")
  }
}

# -----------------------------------------------------------------------------
# 9. Binding Integration
# -----------------------------------------------------------------------------
cat("\n9. Processing binding data...\n")
binding_file <- file.path(tables_path, "binding_integration/integrated_targets_all.csv")
if (file.exists(binding_file)) {
  binding_data <- read_csv(binding_file, show_col_types = FALSE)
  saveRDS(binding_data, file.path(output_path, "binding_integration.rds"), compress = "xz")
  cat("   Binding integration:", nrow(binding_data), "genes\n")
}

# -----------------------------------------------------------------------------
# 10. Temporal Classification
# -----------------------------------------------------------------------------
cat("\n10. Processing temporal data...\n")
temporal_file <- file.path(tables_path, "temporal/temporal_classification_all.csv")
if (file.exists(temporal_file)) {
  temporal_data <- read_csv(temporal_file, show_col_types = FALSE)

  # Rename gene_id to GeneID if needed
  if ("gene_id" %in% names(temporal_data) && !"GeneID" %in% names(temporal_data)) {
    temporal_data <- temporal_data %>% rename(GeneID = gene_id)
  }

  # Merge with JAG1 targets to get FC columns
  fc_cols <- c("FC_TP0", "FC_TP1", "FC_TP2", "FC_TP3", "FC_TP4", "Best_hit_arabi_name")
  fc_cols_present <- fc_cols[fc_cols %in% names(jag1_slim)]

  if (length(fc_cols_present) > 0) {
    temporal_data <- temporal_data %>%
      left_join(
        jag1_slim %>% select(GeneID, all_of(fc_cols_present)),
        by = "GeneID"
      )
  }

  # Check for NA values after join
  fc_non_na <- sum(!is.na(temporal_data$FC_TP0))
  cat("   Genes with FC values:", fc_non_na, "of", nrow(temporal_data), "\n")

  # Filter out genes without FC values for cleaner data
  if (fc_non_na < nrow(temporal_data)) {
    temporal_data <- temporal_data %>% filter(!is.na(FC_TP0))
    cat("   After filtering NA:", nrow(temporal_data), "genes\n")
  }

  saveRDS(temporal_data, file.path(output_path, "temporal_classification.rds"), compress = "xz")
  cat("   Temporal classification:", nrow(temporal_data), "genes\n")
  cat("   Columns:", paste(names(temporal_data), collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# 11. PCA Data (compute from expression matrix)
# -----------------------------------------------------------------------------
cat("\n11. Computing PCA from expression data...\n")
if (exists("tpm_subset")) {
  # Print sample names for debugging
  cat("   Expression matrix sample names (first 3):", paste(head(colnames(tpm_subset), 3), collapse = ", "), "\n")
  cat("   Experimental design Sample (first 3):", paste(head(exp_design$Sample, 3), collapse = ", "), "\n")

  # Remove genes with zero variance
  gene_var <- apply(tpm_subset, 1, var, na.rm = TRUE)

  # Run PCA on top 5000 most variable genes for speed
  top_var_genes <- names(sort(gene_var[gene_var > 0], decreasing = TRUE))[1:min(5000, sum(gene_var > 0))]
  expr_pca <- t(tpm_subset[top_var_genes, ])

  pca_result <- prcomp(expr_pca, scale. = TRUE, center = TRUE)

  # Create PCA data for dashboard - fix sample names
  pca_samples_raw <- rownames(pca_result$x)

  # Expression matrix uses "salmoXXX" format, experimental design uses "XXX"
  # Strip "salmo" prefix to match experimental design Sample column
  pca_samples_clean <- gsub("^salmo", "", pca_samples_raw)

  pca_data <- data.frame(
    Sample = pca_samples_clean,
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3]
  )

  # Variance explained
  var_explained <- summary(pca_result)$importance[2, 1:3] * 100

  # Check if sample names now match experimental design
  match_sample <- sum(pca_samples_clean %in% exp_design$Sample)
  cat("   Sample name matches after cleaning:", match_sample, "of", length(pca_samples_clean), "\n")

  saveRDS(list(coords = pca_data, variance = var_explained), file.path(output_path, "pca_data.rds"))
  cat("   PCA computed:", nrow(pca_data), "samples, variance explained:",
      paste(round(var_explained, 1), "%", collapse = ", "), "\n")
} else {
  cat("   WARNING: Expression matrix not available for PCA\n")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n=== Data Preparation Complete ===\n")
cat("Output directory:", output_path, "\n")

# List created files
files <- list.files(output_path, pattern = "\\.rds$", full.names = TRUE)
file_sizes <- file.info(files)$size / 1024  # KB
cat("\nCreated files:\n")
for (i in seq_along(files)) {
  cat(sprintf("  %s: %.1f KB\n", basename(files[i]), file_sizes[i]))
}
cat(sprintf("\nTotal size: %.1f MB\n", sum(file_sizes) / 1024))
