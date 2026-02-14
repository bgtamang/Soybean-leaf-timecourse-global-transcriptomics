#!/usr/bin/env Rscript
# Preprocess data for Shiny dashboard (Phase3)

library(dplyr)
library(readr)
library(tidyr)

# Paths — Phase3
base_path <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase3-Refined-Analysis"
tables_path <- file.path(base_path, "03_results/tables")
data_path <- file.path(base_path, "01_data")
output_path <- file.path(base_path, "GmJAG1_Dashboard/data")

cat("=== Preparing Shiny Dashboard Data (Phase3) ===\n\n")

# -----------------------------------------------------------------------------
# 1. Experimental Design (Sample Metadata)
# -----------------------------------------------------------------------------
cat("1. Processing experimental design...\n")
exp_design <- read_csv(file.path(tables_path, "experimental_design.csv"), show_col_types = FALSE)
saveRDS(exp_design, file.path(output_path, "experimental_design.rds"), compress = "xz")
cat("   Saved:", nrow(exp_design), "samples\n")

# -----------------------------------------------------------------------------
# 2. Expression Matrix (voom-transformed, batch-corrected — matches publication)
# -----------------------------------------------------------------------------
cat("\n2. Processing expression data (voom-transformed, batch-corrected)...\n")
checkpoint_path <- file.path(base_path, "03_results/checkpoints/06_validated.RData")
load(checkpoint_path)

if (exists("v_primary")) {
  expr_matrix <- v_primary$E

  cat("   Loaded v_primary$E:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  cat("   Sample names (first 3):", paste(head(colnames(expr_matrix), 3), collapse = ", "), "\n")
  cat("   Value range:", round(min(expr_matrix), 2), "to", round(max(expr_matrix), 2), "\n")

  saveRDS(expr_matrix, file.path(output_path, "expression_matrix.rds"), compress = "xz")
  cat("   Saved: expression_matrix.rds\n")

  # Full gene list for lookup
  gene_list <- data.frame(GeneID = rownames(expr_matrix))
  saveRDS(gene_list, file.path(output_path, "gene_list.rds"), compress = "xz")
  cat("   Full gene list:", nrow(gene_list), "genes\n")
} else {
  cat("   WARNING: v_primary not found in checkpoint\n")
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
go_targets_file <- file.path(tables_path, "functional/GO_all_targets_combined.csv")

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

# GO by tier — Phase3 has separate BP/CC/MF files per tier
for (tier in c("gold", "silver", "bronze")) {
  # Try combined first, then fall back to BP
  tier_file <- file.path(tables_path, paste0("functional/GO_", tier, "_tier.csv"))
  tier_file_bp <- file.path(tables_path, paste0("functional/GO_", tier, "_tier_BP.csv"))

  if (file.exists(tier_file)) {
    tier_go <- read_csv(tier_file, show_col_types = FALSE)
  } else if (file.exists(tier_file_bp)) {
    # Combine BP + CC + MF for this tier
    # Force all columns to character first to avoid type conflicts
    coerce_cols <- function(df) {
      if (is.null(df)) return(NULL)
      df %>% mutate(across(where(is.numeric), as.character))
    }
    bp <- tryCatch(read_csv(tier_file_bp, show_col_types = FALSE, col_types = cols(.default = "c")), error = function(e) NULL)
    cc <- tryCatch(read_csv(file.path(tables_path, paste0("functional/GO_", tier, "_tier_CC.csv")), show_col_types = FALSE, col_types = cols(.default = "c")), error = function(e) NULL)
    mf <- tryCatch(read_csv(file.path(tables_path, paste0("functional/GO_", tier, "_tier_MF.csv")), show_col_types = FALSE, col_types = cols(.default = "c")), error = function(e) NULL)
    tier_go <- bind_rows(bp, cc, mf)
    # Convert numeric columns back
    num_cols <- c("Count", "GeneCount", "p.adjust", "pvalue", "qvalue", "FDR", "Fold_Enrichment", "RichFactor", "BgRatio_num", "GeneRatio_num")
    for (nc in intersect(num_cols, names(tier_go))) {
      tier_go[[nc]] <- as.numeric(tier_go[[nc]])
    }
  } else {
    next
  }

  if (!is.null(tier_go) && nrow(tier_go) > 0) {
    if (!"Description" %in% names(tier_go) && "GO_ID" %in% names(tier_go)) {
      tier_go$Description <- tier_go$GO_ID
    }
    saveRDS(tier_go, file.path(output_path, paste0("go_", tier, ".rds")))
    cat("   GO", tier, "tier:", nrow(tier_go), "terms\n")
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
  # Keep top correlated genes (use Cor_LW_Ratio column if present, else first cor column)
  cor_col <- intersect(c("Cor_LW_Ratio", "Cor_V1_Ratio"), names(gene_pheno))[1]
  if (!is.na(cor_col)) {
    gene_pheno_top <- gene_pheno %>%
      arrange(desc(abs(.data[[cor_col]]))) %>%
      head(5000)
  } else {
    gene_pheno_top <- head(gene_pheno, 5000)
  }
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
# 11. PCA Data (from voom-transformed expression — matches publication Fig 2A)
# -----------------------------------------------------------------------------
cat("\n11. Computing PCA from voom-transformed expression...\n")
if (exists("expr_matrix")) {
  # Use ALL genes (same as publication Fig2A_PCA.R)
  expr_pca <- t(expr_matrix)

  pca_result <- prcomp(expr_pca, scale. = TRUE, center = TRUE)

  pca_data <- data.frame(
    Sample = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3]
  )

  var_explained <- summary(pca_result)$importance[2, 1:3] * 100

  # Verify sample name matching
  match_sample <- sum(pca_data$Sample %in% exp_design$Sample)
  cat("   Sample name matches:", match_sample, "of", nrow(pca_data), "\n")

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
