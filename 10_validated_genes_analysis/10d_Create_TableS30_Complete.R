# =============================================================================
# Create Complete Table S30: High-Confidence Validated Targets
# Includes binding, phenotype correlation, WGCNA, and annotation
# =============================================================================

rm(list = ls())
gc()

library(tidyverse)

base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ"
setwd(file.path(base_dir, "Phase3-Refined-Analysis"))

fig6_dir <- "Manuscript/Figure/Figure6"
supp_dir <- "03_results/tables/supplementary_v6"

cat("\n=== Creating Complete Table S30 ===\n\n")

# =============================================================================
# 1. LOAD BASE DATA (validated genes with binding info)
# =============================================================================

validated_genes <- read_csv(file.path(fig6_dir, "Fig6D_Final_Validated_Genes.csv"),
                            show_col_types = FALSE)
cat("Loaded", nrow(validated_genes), "validated genes\n")

# =============================================================================
# 2. ADD PHENOTYPE CORRELATION
# =============================================================================

pheno_cor <- read_csv(file.path(fig6_dir, "Fig6D_All_Target_Phenotype_Correlations.csv"),
                      show_col_types = FALSE)
cat("Loaded phenotype correlations for", nrow(pheno_cor), "genes\n")

# Merge phenotype correlation
validated_genes <- validated_genes %>%
  left_join(pheno_cor %>% select(GeneID, Cor_LW_Ratio, Pval_LW_Ratio, Correlation_Strength),
            by = "GeneID")

# =============================================================================
# 3. ADD WGCNA MODULE INFO
# =============================================================================

hub_genes <- read_csv("03_results/tables/WGCNA/hub_genes.csv", show_col_types = FALSE)
cat("Loaded WGCNA module info for", nrow(hub_genes), "genes\n")

# Merge WGCNA info
validated_genes <- validated_genes %>%
  left_join(hub_genes %>% select(GeneID, Module, Module_Color, kME),
            by = "GeneID")

# =============================================================================
# 4. ADD ARABIDOPSIS ANNOTATION
# =============================================================================

annot_file <- file.path(base_dir, "Phase1-Exploratory", "Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")
if (file.exists(annot_file)) {
  annot <- read_delim(annot_file, delim = "\t", show_col_types = FALSE)

  # Find the relevant columns
  gene_col <- intersect(c("locusName", "transcript", "gene_id", "Gene"), names(annot))[1]

  annot_subset <- annot %>%
    select(any_of(c(gene_col, "Best-hit-arabi-name", "Best-hit-arabi-defline", "Pfam", "GO")))

  if (!is.null(gene_col)) {
    names(annot_subset)[1] <- "GeneID"
  }

  # CRITICAL: Remove duplicate genes (multiple transcripts per gene)
  annot_subset <- annot_subset %>%
    distinct(GeneID, .keep_all = TRUE)
  cat("Annotation deduplicated to", nrow(annot_subset), "unique genes\n")

  # Standardize column names to match requested output
  if ("Best-hit-arabi-name" %in% names(annot_subset)) {
    annot_subset <- annot_subset %>%
      rename(Arabi_Name = `Best-hit-arabi-name`)
  }
  if ("Best-hit-arabi-defline" %in% names(annot_subset)) {
    annot_subset <- annot_subset %>%
      rename(Arabi_Defline = `Best-hit-arabi-defline`)
  }

  # Merge annotation
  validated_genes <- validated_genes %>%
    left_join(annot_subset, by = "GeneID")

  cat("Added Arabidopsis annotation\n")
} else {
  cat("WARNING: Annotation file not found\n")
}

# =============================================================================
# 5. ORGANIZE FINAL TABLE
# =============================================================================

# Select specific columns in requested order
table_s30 <- validated_genes %>%
  select(
    GeneID,
    DE_Tier,
    Mean_logFC,
    Binding_Class,
    Has_ChIPseq,
    Has_DAPseq,
    ChIP_Primary_Region,
    ChIP_N_Peaks,
    ChIP_Closest_Distance,
    any_of(c("ChIP_Promoter")),
    any_of(c("ChIP_Genic")),
    Cor_LW_Ratio,
    Pval_LW_Ratio,
    Correlation_Strength,
    Module,
    Module_Color,
    kME,
    Priority_Score,
    Final_Tier,
    any_of(c("ChIP_Downstream")),
    any_of(c("DAP_Peak")),
    any_of(c("DAP_Motif")),
    any_of(c("DE_Score")),
    any_of(c("Binding_Score")),
    any_of(c("Promoter_Bonus")),
    Arabi_Name,
    Pfam,
    GO,
    any_of(c("Arabi_Defline"))
  ) %>%
  arrange(desc(Priority_Score), desc(abs(Cor_LW_Ratio)))

# =============================================================================
# 6. VERIFY AND SAVE TABLE
# =============================================================================

# Deduplicate (in case joins produced extra rows)
table_s30 <- table_s30 %>% distinct(GeneID, .keep_all = TRUE)
n_genes <- nrow(table_s30)
cat("\nVERIFIED: Table has", n_genes, "high-confidence genes\n")

# Save to supplementary folder
out_name <- paste0("TableS30_HighConfidence_", n_genes, "_Targets.csv")
write_csv(table_s30, file.path(supp_dir, out_name))
cat("Saved:", out_name, "\n")

# Also save to Figure 6 folder
write_csv(table_s30, file.path(fig6_dir, out_name))
cat("Saved copy to Figure6 folder\n")

# =============================================================================
# 7. SUMMARY STATISTICS
# =============================================================================

cat("\n=== Table S30 Summary ===\n")
cat("Total genes:", nrow(table_s30), "\n")

cat("\nBinding breakdown:\n")
print(table(table_s30$Binding_Class))

cat("\nDE Tier breakdown:\n")
print(table(table_s30$DE_Tier))

cat("\nWGCNA Module breakdown:\n")
module_summary <- table_s30 %>%
  group_by(Module_Color) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(desc(N))
print(as.data.frame(module_summary))

cat("\nPhenotype correlation range:\n")
cat("  Min r:", round(min(table_s30$Cor_LW_Ratio, na.rm = TRUE), 3), "\n")
cat("  Max r:", round(max(table_s30$Cor_LW_Ratio, na.rm = TRUE), 3), "\n")
cat("  Mean |r|:", round(mean(abs(table_s30$Cor_LW_Ratio), na.rm = TRUE), 3), "\n")

cat("\n=== COMPLETE ===\n")
