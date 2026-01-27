
# Set base directory and working directory
base_dir <- "C:/Users/bgtamang/OneDrive - University of Illinois - Urbana/Desktop/Soybean-RNASEQ/Phase2-Refined-Analysis"
setwd(base_dir)

library(tidyverse)
library(edgeR)
library(limma)
library(pheatmap)
library(RColorBrewer)

# 1. Load Required Data

# Load normalized expression data (contains v_full, targets, design_full, etc.)
load(file.path(base_dir, "03_results/checkpoints/05_normalized.RData"))

# Load DE results
load(file.path(base_dir, "03_results/checkpoints/12_DE_organized.RData"))

# Load JAG1 targets
jag1_targets <- read_csv(file.path(base_dir, "03_results/tables/JAG1_targets/JAG1_targets_FINAL.csv"),
                         show_col_types = FALSE)

# 2. Define KRP and TPR Gene Families

# Soybean KRP homologs (ICK/KIP-related proteins - cell cycle inhibitors)
# SOURCE: Guo et al. 2023, Frontiers in Plant Science
# "Characterization of the soybean KRP gene family reveals a key role for
#  GmKRP2a in root development"
# DOI: 10.3389/fpls.2023.1096467
# PMCID: PMC9911667
#
# In Arabidopsis, JAG directly represses KRP2/KRP4 - if conserved,
# these genes should be UPREGULATED in narrow (loss of JAG1 repression)
krp_genes <- tribble(
  ~gene_id,            ~family,    ~arabidopsis_ortholog, ~description,
  "Glyma.01G185400",   "GmKRP1a",  "AT2G23430 (KRP1)",   "Cyclin-dependent kinase inhibitor",
  "Glyma.11G056700",   "GmKRP1b",  "AT2G23430 (KRP1)",   "Cyclin-dependent kinase inhibitor",
  "Glyma.05G085000",   "GmKRP2a",  "AT3G50630 (KRP2)",   "Cell cycle inhibitor - JAG target in At",
  "Glyma.17G175700",   "GmKRP2b",  "AT3G50630 (KRP2)",   "Cell cycle inhibitor - JAG target in At",
  "Glyma.08G354300",   "GmKRP3",   "AT5G48820 (KRP3)",   "Cyclin-dependent kinase inhibitor",
  "Glyma.20G198800",   "GmKRP4",   "AT2G32710 (KRP4)",   "Cell cycle inhibitor - JAG target in At",
  "Glyma.18G170800",   "GmKRP5",   "AT3G24810 (KRP5)",   "Cyclin-dependent kinase inhibitor",
  "Glyma.02G133700",   "GmKRP6",   "AT3G61630 (KRP6)",   "Cyclin-dependent kinase inhibitor",
  "Glyma.07G211300",   "GmKRP7",   "AT1G49620 (KRP7)",   "Cyclin-dependent kinase inhibitor"
)

# Soybean TPL/TPR homologs (TOPLESS co-repressors)
# These interact with EAR-containing TFs like JAG1
#
# SOURCE: Yu et al. 2025, Scientific Reports
# "Global identification and characterization of soybean TPR genes with
#  expression analysis under photoperiod variations"
# DOI: 10.1038/s41598-025-10368-5
#
# Gene IDs from Table 1 of the paper (VERIFIED):
# - Identified using LisH, CTLH, and WD40 domains
# - Collinearity with Arabidopsis: AtTPL/AtTPR1 → GmTPR1,2,4,10; AtTPR3 → GmTPR6,8
tpr_genes <- tribble(
  ~gene_id,            ~family,     ~arabidopsis_ortholog, ~description,
  "Glyma.08G214600",   "GmTPR1",    "AT1G15750 (TPL)",    "TOPLESS protein 1",
  "Glyma.07G028000",   "GmTPR2",    "AT1G15750 (TPL)",    "TOPLESS protein 2",
  "Glyma.13G367300",   "GmTPR3",    "AT1G80490 (TPR1)",   "TOPLESS-related protein 3",
  "Glyma.15G006100",   "GmTPR4",    "AT1G80490 (TPR1)",   "TOPLESS-related protein 4",
  "Glyma.03G233400",   "GmTPR5",    "AT3G16830 (TPR2)",   "TOPLESS-related protein 5",
  "Glyma.19G230500",   "GmTPR6",    "AT5G27030 (TPR3)",   "TOPLESS-related protein 6",
  "Glyma.20G238400",   "GmTPR7",    "AT3G15880 (TPR4)",   "TOPLESS-related protein 7",
  "Glyma.10G150000",   "GmTPR8",    "AT5G27030 (TPR3)",   "TOPLESS-related protein 8",
  "Glyma.17G112500",   "GmTPR9",    "TPL family",         "TOPLESS-related protein 9",
  "Glyma.13G158800",   "GmTPR10",   "AT1G80490 (TPR1)",   "TOPLESS-related protein 10",
  "Glyma.04G065000",   "GmTPR11",   "Soybean-specific",   "TOPLESS-related protein 11",
  "Glyma.06G066200",   "GmTPR12",   "Soybean-specific",   "TOPLESS-related protein 12"
)

cat("Defined", nrow(krp_genes), "KRP genes and", nrow(tpr_genes), "TPR genes\n\n")

# 2b. Alternative Method: Search Annotation File
# This section provides an independent method to identify KRP and TPR genes
# using the Phytozome annotation file, as a cross-validation approach


# Path to Phytozome annotation file (used in Script 22 for GO analysis)
annotation_file <- file.path(dirname(base_dir), "Phase1-Exploratory/Gmax_880_Wm82.a6.v1.P14.annotation_info.txt")

if(file.exists(annotation_file)) {
  cat("Loading annotation file for cross-validation...\n")

  # Read annotation file
  annotation <- read_delim(annotation_file, delim = "\t",
                           col_types = cols(.default = "c"),
                           show_col_types = FALSE)

  cat("Annotation file contains", nrow(annotation), "entries\n")

  # 2b.1 Search for KRP genes by Pfam domain (PF02234 = CDK inhibitor)
  cat("\nSearching for KRP genes by Pfam domain PF02234...\n")

  krp_from_annotation <- annotation %>%
    filter(grepl("PF02234", Pfam, ignore.case = TRUE) |
           grepl("CYCLIN-DEPENDENT.*KINASE.*INHIBITOR|CDK.*INHIBITOR|KIP.*RELATED|ICK",
                 `Best-hit-arabi-defline`, ignore.case = TRUE)) %>%
    select(locusName, `Best-hit-arabi-name`, `Best-hit-arabi-defline`, Pfam) %>%
    distinct(locusName, .keep_all = TRUE)

  cat("Found", nrow(krp_from_annotation), "potential KRP genes from annotation\n")

  # Cross-validate with literature list
  literature_krp_ids <- krp_genes$gene_id
  annotation_krp_ids <- krp_from_annotation$locusName

  krp_in_both <- intersect(literature_krp_ids, annotation_krp_ids)
  krp_lit_only <- setdiff(literature_krp_ids, annotation_krp_ids)
  krp_anno_only <- setdiff(annotation_krp_ids, literature_krp_ids)

  cat("\nKRP Cross-validation:\n")
  cat("- In both literature and annotation:", length(krp_in_both), "\n")
  cat("- Literature only:", length(krp_lit_only),
      ifelse(length(krp_lit_only) > 0, paste0(" (", paste(krp_lit_only, collapse=", "), ")"), ""), "\n")
  cat("- Annotation only:", length(krp_anno_only),
      ifelse(length(krp_anno_only) > 0, paste0(" (", paste(krp_anno_only, collapse=", "), ")"), ""), "\n")

  # 2b.2 Search for TPR genes by keyword (TOPLESS)
  cat("\nSearching for TPR genes by 'TOPLESS' keyword...\n")

  tpr_from_annotation <- annotation %>%
    filter(grepl("TOPLESS", `Best-hit-arabi-defline`, ignore.case = TRUE) |
           grepl("TPL|TPR[0-9]", `Best-hit-arabi-name`, ignore.case = TRUE)) %>%
    select(locusName, `Best-hit-arabi-name`, `Best-hit-arabi-defline`, Pfam) %>%
    distinct(locusName, .keep_all = TRUE)

  cat("Found", nrow(tpr_from_annotation), "potential TPR genes from annotation\n")

  # Cross-validate with literature list
  literature_tpr_ids <- tpr_genes$gene_id
  annotation_tpr_ids <- tpr_from_annotation$locusName

  tpr_in_both <- intersect(literature_tpr_ids, annotation_tpr_ids)
  tpr_lit_only <- setdiff(literature_tpr_ids, annotation_tpr_ids)
  tpr_anno_only <- setdiff(annotation_tpr_ids, literature_tpr_ids)

  cat("\nTPR Cross-validation:\n")
  cat("- In both literature and annotation:", length(tpr_in_both), "\n")
  cat("- Literature only:", length(tpr_lit_only),
      ifelse(length(tpr_lit_only) > 0, paste0(" (", paste(tpr_lit_only, collapse=", "), ")"), ""), "\n")
  cat("- Annotation only:", length(tpr_anno_only),
      ifelse(length(tpr_anno_only) > 0, paste0(" (", paste(tpr_anno_only, collapse=", "), ")"), ""), "\n")

  # Note: The annotation file may not find all TPR genes because Yu et al. 2025
  # used domain-based identification (LisH + CTLH + WD40) which may not be
  # captured by keyword searches. The 12 genes from the paper are authoritative.

  cat("\nNOTE: TPR gene identification from Yu et al. 2025 used structural domain\n")
  cat("analysis (LisH, CTLH, WD40 domains) which is more comprehensive than\n")
  cat("annotation keyword search. Literature-based IDs are authoritative.\n")

  # Store annotation search results for report
  annotation_search_results <- list(
    krp = list(
      from_annotation = krp_from_annotation,
      cross_validation = list(
        in_both = krp_in_both,
        lit_only = krp_lit_only,
        anno_only = krp_anno_only
      )
    ),
    tpr = list(
      from_annotation = tpr_from_annotation,
      cross_validation = list(
        in_both = tpr_in_both,
        lit_only = tpr_lit_only,
        anno_only = tpr_anno_only
      )
    )
  )

} else {
  cat("WARNING: Annotation file not found at:", annotation_file, "\n")
  cat("Proceeding with literature-based gene lists only.\n")
  annotation_search_results <- NULL
}

cat("\n")

# 3. Check Which Genes Are in Dataset

all_genes <- rownames(v_full$E)

# Check KRP presence
krp_genes <- krp_genes %>%
  mutate(in_dataset = gene_id %in% all_genes)

cat("KRP genes in dataset:", sum(krp_genes$in_dataset), "/", nrow(krp_genes), "\n")
cat("Missing KRP genes:", paste(krp_genes$gene_id[!krp_genes$in_dataset], collapse = ", "), "\n\n")

# Check TPR presence
tpr_genes <- tpr_genes %>%
  mutate(in_dataset = gene_id %in% all_genes)

cat("TPR genes in dataset:", sum(tpr_genes$in_dataset), "/", nrow(tpr_genes), "\n")
cat("Missing TPR genes:", paste(tpr_genes$gene_id[!tpr_genes$in_dataset], collapse = ", "), "\n\n")

# 4. Extract Expression Data

# Get expression for present genes
krp_present <- krp_genes %>% filter(in_dataset)
tpr_present <- tpr_genes %>% filter(in_dataset)

# Extract expression matrices
if(nrow(krp_present) > 0) {
  krp_expr <- v_full$E[krp_present$gene_id, , drop = FALSE]
  cat("KRP expression matrix:", nrow(krp_expr), "genes x", ncol(krp_expr), "samples\n")
}

if(nrow(tpr_present) > 0) {
  tpr_expr <- v_full$E[tpr_present$gene_id, , drop = FALSE]
  cat("TPR expression matrix:", nrow(tpr_expr), "genes x", ncol(tpr_expr), "samples\n")
}

# 5. Check KRP Genes in JAG1 Targets

cat("\n=== KRP Genes in JAG1 Target List ===\n")

krp_in_targets <- krp_present %>%
  mutate(is_JAG1_target = gene_id %in% jag1_targets$GeneID)

krp_target_info <- krp_in_targets %>%
  left_join(jag1_targets %>% select(GeneID, Confidence_Tier, Mean_logFC_Pairwise),
            by = c("gene_id" = "GeneID"))

print(krp_target_info %>% select(gene_id, family, arabidopsis_ortholog,
                                  is_JAG1_target, Confidence_Tier, Mean_logFC_Pairwise))

cat("\nKRP genes that are JAG1 targets:", sum(krp_target_info$is_JAG1_target, na.rm = TRUE), "\n")

# 6. Check TPR Genes in JAG1 Targets

cat("\n=== TPR Genes in JAG1 Target List ===\n")

tpr_in_targets <- tpr_present %>%
  mutate(is_JAG1_target = gene_id %in% jag1_targets$GeneID)

tpr_target_info <- tpr_in_targets %>%
  left_join(jag1_targets %>% select(GeneID, Confidence_Tier, Mean_logFC_Pairwise),
            by = c("gene_id" = "GeneID"))

print(tpr_target_info %>% select(gene_id, family, arabidopsis_ortholog,
                                  is_JAG1_target, Confidence_Tier, Mean_logFC_Pairwise))

cat("\nTPR genes that are JAG1 targets:", sum(tpr_target_info$is_JAG1_target, na.rm = TRUE), "\n")

# 7. Expression Analysis by Line

cat("\n=== Expression by Line ===\n")

# Add sample metadata
sample_info <- targets %>%
  mutate(sample_id = rownames(targets)) %>%
  select(sample_id, Leaf_type, Timepoint, Line)

# Calculate mean expression by leaf type for KRP genes
if(nrow(krp_present) > 0) {
  krp_expr_long <- krp_expr %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "expression") %>%
    left_join(sample_info, by = "sample_id") %>%
    left_join(krp_present %>% select(gene_id, family), by = "gene_id")

  krp_summary <- krp_expr_long %>%
    group_by(gene_id, family, Leaf_type) %>%
    summarize(mean_expr = mean(expression, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Leaf_type, values_from = mean_expr) %>%
    mutate(narrow_vs_broad = Narrow - Broad)

  cat("\nKRP Expression Summary (Narrow vs Broad):\n")
  print(krp_summary)

  cat("\nInterpretation:\n")
  cat("- Positive narrow_vs_broad = Higher in Narrow (expected if JAG1 represses KRPs)\n")
  cat("- Negative narrow_vs_broad = Lower in Narrow (unexpected)\n")
}

# Same for TPR genes
if(nrow(tpr_present) > 0) {
  tpr_expr_long <- tpr_expr %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "expression") %>%
    left_join(sample_info, by = "sample_id") %>%
    left_join(tpr_present %>% select(gene_id, family), by = "gene_id")

  tpr_summary <- tpr_expr_long %>%
    group_by(gene_id, family, Leaf_type) %>%
    summarize(mean_expr = mean(expression, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Leaf_type, values_from = mean_expr) %>%
    mutate(narrow_vs_broad = Narrow - Broad)

  cat("\n\nTPR Expression Summary (Narrow vs Broad):\n")
  print(tpr_summary)
}

# 8. Temporal Expression Patterns

cat("\n=== Temporal Expression Patterns ===\n")

# KRP temporal patterns
if(nrow(krp_present) > 0) {
  krp_temporal <- krp_expr_long %>%
    group_by(gene_id, family, Leaf_type, Timepoint) %>%
    summarize(mean_expr = mean(expression, na.rm = TRUE),
              .groups = "drop")

  # Plot
  p_krp <- ggplot(krp_temporal, aes(x = Timepoint, y = mean_expr,
                                     color = Leaf_type, group = Leaf_type)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~gene_id + family, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("Broad" = "#2166AC", "Narrow" = "#B2182B")) +
    theme_bw() +
    labs(title = "KRP Gene Expression Across Timepoints",
         subtitle = "Expected: Higher in Narrow if JAG1 represses KRPs",
         x = "Timepoint", y = "Log2 Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(base_dir, "03_results/figures/34_KRP_TPR/KRP_temporal_expression.pdf"),
         p_krp, width = 14, height = 10)
  cat("Saved KRP temporal expression plot\n")
}

# TPR temporal patterns
if(nrow(tpr_present) > 0) {
  tpr_temporal <- tpr_expr_long %>%
    group_by(gene_id, family, Leaf_type, Timepoint) %>%
    summarize(mean_expr = mean(expression, na.rm = TRUE),
              .groups = "drop")

  p_tpr <- ggplot(tpr_temporal, aes(x = Timepoint, y = mean_expr,
                                     color = Leaf_type, group = Leaf_type)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    facet_wrap(~gene_id + family, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("Broad" = "#2166AC", "Narrow" = "#B2182B")) +
    theme_bw() +
    labs(title = "TPR/TPL Gene Expression Across Timepoints",
         x = "Timepoint", y = "Log2 Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(base_dir, "03_results/figures/34_KRP_TPR/TPR_temporal_expression.pdf"),
         p_tpr, width = 14, height = 10)
  cat("Saved TPR temporal expression plot\n")
}

# 9. Heatmap of KRP and TPR Genes

cat("\n=== Creating Expression Heatmaps ===\n")

# Create output directory
dir.create(file.path(base_dir, "03_results/figures/34_KRP_TPR"), showWarnings = FALSE, recursive = TRUE)

# Combined heatmap
combined_genes <- c(krp_present$gene_id, tpr_present$gene_id)
combined_info <- bind_rows(
  krp_present %>% mutate(gene_family = "KRP"),
  tpr_present %>% mutate(gene_family = "TPR")
)

if(length(combined_genes) > 0) {
  combined_expr <- v_full$E[combined_genes, , drop = FALSE]

  # Scale by row
  combined_scaled <- t(scale(t(combined_expr)))

  # Annotation for rows
  row_anno <- combined_info %>%
    select(gene_id, gene_family, family) %>%
    column_to_rownames("gene_id")

  # Annotation for columns
  col_anno <- targets %>%
    select(Leaf_type, Timepoint) %>%
    as.data.frame()
  rownames(col_anno) <- colnames(combined_scaled)

  # Colors
  anno_colors <- list(
    Leaf_type = c("Broad" = "#2166AC", "Narrow" = "#B2182B"),
    Timepoint = c("TP0" = "#FEE5D9", "TP1" = "#FCAE91", "TP2" = "#FB6A4A",
                  "TP3" = "#DE2D26", "TP4" = "#A50F15"),
    gene_family = c("KRP" = "#1B9E77", "TPR" = "#D95F02")
  )

  pdf(file.path(base_dir, "03_results/figures/34_KRP_TPR/KRP_TPR_heatmap.pdf"), width = 12, height = 8)
  pheatmap(combined_scaled,
           annotation_row = row_anno,
           annotation_col = col_anno,
           annotation_colors = anno_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = FALSE,
           main = "KRP and TPR Gene Expression\n(Row-scaled)")
  dev.off()
  cat("Saved combined heatmap\n")
}

# 10. Statistical Testing

cat("\n=== Statistical Testing: Narrow vs Broad ===\n")

# Function to test differential expression for specific genes
test_DE <- function(gene_ids, expr_matrix, design_info) {
  results <- data.frame()

  for(gene in gene_ids) {
    if(gene %in% rownames(expr_matrix)) {
      expr <- expr_matrix[gene, ]
      broad_expr <- expr[design_info$Leaf_type == "Broad"]
      narrow_expr <- expr[design_info$Leaf_type == "Narrow"]

      # t-test
      test <- t.test(narrow_expr, broad_expr)

      results <- rbind(results, data.frame(
        gene_id = gene,
        mean_broad = mean(broad_expr),
        mean_narrow = mean(narrow_expr),
        logFC = mean(narrow_expr) - mean(broad_expr),
        pvalue = test$p.value
      ))
    }
  }

  if(nrow(results) > 0) {
    results$FDR <- p.adjust(results$pvalue, method = "BH")
    results$significant <- results$FDR < 0.05 & abs(results$logFC) > 1
  }

  return(results)
}

# Test KRP genes
krp_DE <- test_DE(krp_present$gene_id, v_full$E, targets)
krp_DE <- krp_DE %>%
  left_join(krp_present %>% select(gene_id, family, arabidopsis_ortholog), by = "gene_id")

cat("\nKRP Differential Expression Results:\n")
print(krp_DE %>% arrange(pvalue))

# Test TPR genes
tpr_DE <- test_DE(tpr_present$gene_id, v_full$E, targets)
tpr_DE <- tpr_DE %>%
  left_join(tpr_present %>% select(gene_id, family, arabidopsis_ortholog), by = "gene_id")

cat("\nTPR Differential Expression Results:\n")
print(tpr_DE %>% arrange(pvalue))

# 11. Summary and Interpretation

cat("\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n=== SUMMARY AND BIOLOGICAL INTERPRETATION ===\n")
cat("=" %>% rep(70) %>% paste(collapse = ""))
cat("\n")

# KRP summary
krp_up_in_narrow <- sum(krp_DE$logFC > 0 & krp_DE$significant, na.rm = TRUE)
krp_down_in_narrow <- sum(krp_DE$logFC < 0 & krp_DE$significant, na.rm = TRUE)

cat("\nKRP GENES (Cell Cycle Inhibitors):\n")
cat("- Expected pattern: UPREGULATED in Narrow (loss of JAG1 repression)\n")
cat("- Significantly UP in Narrow:", krp_up_in_narrow, "\n")
cat("- Significantly DOWN in Narrow:", krp_down_in_narrow, "\n")
cat("- In JAG1 target list:", sum(krp_target_info$is_JAG1_target, na.rm = TRUE), "\n")

if(krp_up_in_narrow > 0) {
  cat("\n*** SUPPORTS Arabidopsis JAG-KRP mechanism! ***\n")
} else if(krp_down_in_narrow > 0) {
  cat("\n*** OPPOSITE of expected - novel mechanism? ***\n")
} else {
  cat("\n*** No significant KRP changes - mechanism may differ from Arabidopsis ***\n")
}

# TPR summary
tpr_up_in_narrow <- sum(tpr_DE$logFC > 0 & tpr_DE$significant, na.rm = TRUE)
tpr_down_in_narrow <- sum(tpr_DE$logFC < 0 & tpr_DE$significant, na.rm = TRUE)

cat("\nTPR/TPL GENES (Co-repressors):\n")
cat("- These interact with JAG1's EAR motif\n")
cat("- Significantly UP in Narrow:", tpr_up_in_narrow, "\n")
cat("- Significantly DOWN in Narrow:", tpr_down_in_narrow, "\n")
cat("- In JAG1 target list:", sum(tpr_target_info$is_JAG1_target, na.rm = TRUE), "\n")

# 12. Save Results

dir.create(file.path(base_dir, "03_results/tables/mechanistic"), showWarnings = FALSE, recursive = TRUE)

# Save KRP results
write_csv(krp_DE, file.path(base_dir, "03_results/tables/mechanistic/KRP_expression_analysis.csv"))
write_csv(krp_target_info, file.path(base_dir, "03_results/tables/mechanistic/KRP_JAG1_target_status.csv"))

# Save TPR results
write_csv(tpr_DE, file.path(base_dir, "03_results/tables/mechanistic/TPR_expression_analysis.csv"))
write_csv(tpr_target_info, file.path(base_dir, "03_results/tables/mechanistic/TPR_JAG1_target_status.csv"))

# Save combined summary
combined_summary <- bind_rows(
  krp_DE %>% mutate(gene_family = "KRP"),
  tpr_DE %>% mutate(gene_family = "TPR")
)
write_csv(combined_summary, file.path(base_dir, "03_results/tables/mechanistic/KRP_TPR_combined_results.csv"))

# Save annotation search results if available
if(!is.null(annotation_search_results)) {
  # KRP annotation search
  write_csv(annotation_search_results$krp$from_annotation,
            file.path(base_dir, "03_results/tables/mechanistic/KRP_annotation_search.csv"))

  # TPR annotation search
  write_csv(annotation_search_results$tpr$from_annotation,
            file.path(base_dir, "03_results/tables/mechanistic/TPR_annotation_search.csv"))

  # Cross-validation summary
  cross_val_summary <- data.frame(
    gene_family = c("KRP", "KRP", "KRP", "TPR", "TPR", "TPR"),
    category = c("In both sources", "Literature only", "Annotation only",
                 "In both sources", "Literature only", "Annotation only"),
    count = c(
      length(annotation_search_results$krp$cross_validation$in_both),
      length(annotation_search_results$krp$cross_validation$lit_only),
      length(annotation_search_results$krp$cross_validation$anno_only),
      length(annotation_search_results$tpr$cross_validation$in_both),
      length(annotation_search_results$tpr$cross_validation$lit_only),
      length(annotation_search_results$tpr$cross_validation$anno_only)
    ),
    genes = c(
      paste(annotation_search_results$krp$cross_validation$in_both, collapse = "; "),
      paste(annotation_search_results$krp$cross_validation$lit_only, collapse = "; "),
      paste(annotation_search_results$krp$cross_validation$anno_only, collapse = "; "),
      paste(annotation_search_results$tpr$cross_validation$in_both, collapse = "; "),
      paste(annotation_search_results$tpr$cross_validation$lit_only, collapse = "; "),
      paste(annotation_search_results$tpr$cross_validation$anno_only, collapse = "; ")
    )
  )
  write_csv(cross_val_summary, file.path(base_dir, "03_results/tables/mechanistic/gene_identification_cross_validation.csv"))
}

cat("\nResults saved to 03_results/tables/mechanistic/\n")

# 13. Save Checkpoint

save(krp_genes, tpr_genes, krp_present, tpr_present,
     krp_DE, tpr_DE, krp_target_info, tpr_target_info,
     krp_summary, tpr_summary, annotation_search_results,
     file = file.path(base_dir, "03_results/checkpoints/34_KRP_TPR_analysis.RData"))

# 14. Gene Identification Method Documentation

cat("\n=== GENE IDENTIFICATION METHODS DOCUMENTATION ===\n")
cat("\nKRP GENES:\n")
cat("Primary source: Guo et al. 2023, Frontiers in Plant Science\n")
cat("DOI: 10.3389/fpls.2023.1096467\n")
cat("Method: Comprehensive identification using phylogenetics and Pfam domain PF02234\n")
cat("Result: 9 GmKRP genes (GmKRP1a/b, GmKRP2a/b, GmKRP3-7)\n")

cat("\nTPR GENES:\n")
cat("Primary source: Yu et al. 2025, Scientific Reports\n")
cat("DOI: 10.1038/s41598-025-10368-5\n")
cat("Method: Domain-based search (LisH + CTLH + WD40) and phylogenetic analysis\n")
cat("Result: 12 GmTPR genes (GmTPR1-12)\n")

cat("\nCROSS-VALIDATION METHOD:\n")
cat("Secondary source: Phytozome annotation file (Gmax_880_Wm82.a6.v1)\n")
cat("KRP search: Pfam domain PF02234 + keyword matching\n")
cat("TPR search: TOPLESS keyword in description\n")
cat("Purpose: Independent validation of gene family membership\n")

cat("\nCheckpoint saved: 34_KRP_TPR_analysis.RData\n")
cat("\n=== Script 34 Complete ===\n")
